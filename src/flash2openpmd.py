#!/usr/bin/env python
"""
This file is part of Flash2OpenPMD software, which
converts the FLASH files to openPMD
Copyright (C) 2022, Flash2OpenPMD contributors
Author: Borsos Andrei Paul, Jian Fuh Ong

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import yt
import os
import numpy as np
import sys
import argparse
from scipy import constants as sc
import openpmd_api  as io
import glob
from scipy import interpolate


def list_available_fields(run_directory, filename):
    """
    The "gas"-typed field names yt exposes for this FLASH file (on-disk and
    derived), suitable for the --fields/Field parameter.
    """

    ds = yt.load(os.path.join(run_directory, filename))
    return sorted({name for ftype, name in ds.derived_field_list if ftype == "gas"})


def format_bytes(n):
    """Human-readable byte count, e.g. 8.01 MB."""

    n = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(n) < 1024:
            return f"{n:.2f} {unit}"
        n /= 1024
    return f"{n:.2f} PB"


def resolve_device(requested="auto"):
    """
    Resolve a requested device ("auto"/"cuda"/"cpu") to an actual device,
    falling back to CPU when torch or CUDA is unavailable.

    torch.cuda.is_available() can return True while every real kernel
    launch fails (e.g. a pip PyTorch build whose compute-capability targets
    don't include an older GPU), so this runs one trivial op on the device
    to confirm CUDA actually works before recommending it.
    """

    if requested == "cpu":
        return "cpu"

    try:
        import torch
    except ImportError:
        if requested == "cuda":
            print("PyTorch not installed; falling back to CPU.", file=sys.stderr)
        return "cpu"

    if torch.cuda.is_available():
        try:
            (torch.zeros(1, device="cuda") + 1).cpu()
            return "cuda"
        except RuntimeError as exc:
            print(f"CUDA reported available but a test kernel failed ({exc}); "
                  "falling back to CPU.", file=sys.stderr)
            return "cpu"

    if requested == "cuda":
        print("CUDA not available; falling back to CPU.", file=sys.stderr)
    return "cpu"


def _interp_radial_torch(rho, r, r_grid_xy, chunk_z=128):
    """
    GPU-vectorized equivalent of the per-slice np.interp loop in
    Convert.cylindricalRotateAxial: linear radial interpolation of every
    z-slice at once via torch.searchsorted, instead of looping per slice.
    """

    import torch

    device = torch.device("cuda")
    r_t = torch.from_numpy(np.asarray(r)).to(device, dtype=torch.float32)
    rho_t = torch.from_numpy(np.asarray(rho[:, :, 0])).to(device, dtype=torch.float32)
    q = torch.from_numpy(r_grid_xy.ravel()).to(device, dtype=torch.float32)

    idx = torch.searchsorted(r_t, q).clamp(1, r_t.shape[0] - 1)
    r_lo = r_t[idx - 1]
    r_hi = r_t[idx]
    w = ((q - r_lo) / (r_hi - r_lo)).clamp(0.0, 1.0)

    Nz = rho_t.shape[1]
    Ny, Nx = r_grid_xy.shape
    out = np.empty((Nz, Ny, Nx), dtype=np.float32)

    while True:
        try:
            for start in range(0, Nz, chunk_z):
                end = min(start + chunk_z, Nz)
                lo = rho_t[idx - 1, start:end]
                hi = rho_t[idx, start:end]
                vals = lo * (1.0 - w).unsqueeze(1) + hi * w.unsqueeze(1)
                out[start:end] = vals.T.reshape(end - start, Ny, Nx).cpu().numpy()
            break
        except torch.cuda.OutOfMemoryError:
            torch.cuda.empty_cache()
            if chunk_z <= 8:
                raise
            chunk_z = max(chunk_z // 2, 8)

    return out


class Convert(object):
    def __init__(self,run_directory,filename):
        """
        Parameters
        ----------
        run_directory : string 
            The path to the directory where the output of FLASH are. 

        filename : string
            The output filename of FLASH. (e.g. lasslab_hdf5_plt_cnt_0043)
        """

        if run_directory is None:
            raise ValueError('The run_directory parameter can not be None!')
            
        self.run_directory = run_directory
        self.path_to_sim = os.path.join(run_directory,filename)

    def read4Flash(self,fields,level,geometry=None,device="cpu"):
        """
        Read the FLASH output

        Parameters
        ----------
        x_min : float
            The minimum boundary in x direction in the unit of centimeter (unit of FLASH code)

        y_min : float
            The minimum boundary in y direction in the unit of centimeter (unit of FLASH code)

        z_min : float
            The minimum boundary in y direction in the unit of centimeter (unit of FLASH code)
            z_min=0 for 2D simulation

        field : string
        	The field data (unit of FLASH code), i.e. fields="El_number_density"

        level : int
            The level of refinement using yt-project. Currently fixed to level=4

        geometry : string or None
            Override the geometry detected by yt ("cartesian" or "cylindrical").
            None (default) uses the detected geometry.

        device : string
            "cpu" or "cuda". Only affects the cylindrical rotation step.

        Returns
        -------
        A 2d array containing the required density, shape and maximum boundary in each direction.

        """

        ds = yt.load(self.path_to_sim)

        ND = ds.dimensionality
        geom = geometry or ds.geometry

        print(f"Data geometry: {ds.geometry}" + (f" (overridden to {geom})" if geometry else ""))
        print(f"Data dimensionality: {ND}D")

        if ND == 2:
            scale = np.array([2**level,2**level,1])
        elif ND == 3:
            scale = np.array([2**level,2**level,2**level])

        ds.force_periodicity()
        all_data = ds.smoothed_covering_grid(level=level,
                                    left_edge=ds.domain_left_edge,
                                    dims=ds.domain_dimensions * scale,
                                    )

        if geom == 'cylindrical':
            density = self.cylindricalRotateAxial(all_data,fields=f"{fields}",device=device)
        else:
            density = all_data[('gas',f"{fields}")]

        return density, ND, geom
    
    def cylindricalRotateAxial(self,array,fields,device="cpu"):
        """
        Rotate the 2D cylindrical data to 3D Cartesian data

        Parameters
        ----------
        array : 2D array
            The 2D cylindrical data

        device : string
            "cpu" runs the original np.interp loop. "cuda" vectorizes the
            same linear interpolation across all z-slices on the GPU, with
            automatic fallback to "cpu" if torch/CUDA is unavailable or runs
            out of memory.

        Returns
        -------
        A 3D array containing the rotated Cartesian data.

        """

        rho = array[(f'{fields}')]
        r = array[("index","r")].v[:,0].ravel()
        z = array[("index","z")].v[0,:].ravel()

        Nx = len(2*r)
        Ny = len(2*r)
        X = np.linspace(-r.max(), r.max(), Nx)
        Y = np.linspace(-r.max(), r.max(), Ny)
        Z_out = z.copy()

        xx, yy = np.meshgrid(X, Y, indexing='xy')   # shape (Ny, Nx)
        r_grid_xy = np.sqrt(xx**2 + yy**2)          # radial distance at each (x,y)

        Nz = len(Z_out)

        if not np.all(np.diff(r) > 0):
            raise ValueError("r must be strictly increasing")

        if device == "cuda":
            try:
                print("Converting cylindrical to 3D Cartesian coordinates on GPU...")
                return _interp_radial_torch(rho, r, r_grid_xy)
            except (ImportError, RuntimeError) as exc:
                print(f"GPU interpolation failed ({exc}); falling back to CPU.", file=sys.stderr)

        print("Converting cylindrical to 3D Cartesian coordinates...")

        vol = np.zeros((Nz, Ny, Nx), dtype=rho.dtype)  # (z, y, x)
        rflat = r_grid_xy.ravel()
        for iz in range(Nz):
            slice_r = rho[:, iz, 0]
            vals_flat = np.interp(rflat, r, slice_r, left=slice_r[0], right=slice_r[-1])
            vol[iz, :, :] = vals_flat.reshape(Ny, Nx)

        return vol

    def get_data(self,fields,level,dtype=None,device="cpu",geometry=None):

        dens, ND, geom = self.read4Flash(fields,level,geometry=geometry,device=device)

        if ND == 2:
            if geom == 'cylindrical':
                density = np.array(dens.T, order='C', dtype=dtype)
            else:
                density = np.array(dens[:,:,0], order='C', dtype=dtype)
        elif ND == 3:
            density = np.array(dens.T, order='C', dtype=dtype)

        return density

    def write2openpmd(self,density_input,species,output_name,author="Your Name <Your@email>"):

        os.makedirs(self.run_directory, exist_ok=True)
        output_path = os.path.join(self.run_directory, f"{species}_{output_name}_0.h5")

        series_out = io.Series(output_path, io.Access.create)

        k = series_out.iterations[0]

        series_out.author = author
        # record - again,important to specify as scalar
        n_e_out = k.meshes[f"{species}_density"]
        n_e_mrc = n_e_out[io.Mesh_Record_Component.SCALAR]

        dataset = io.Dataset(
            density_input.dtype,
            density_input.shape)
        n_e_mrc.reset_dataset(dataset)

        n_e_out.set_attribute('dataOrder','C')
        n_e_out.set_attribute('axisLabels',['z','y','x'])
        n_e_mrc.store_chunk(density_input)
        # After registering a data chunk such as x_data and y_data,
        # it MUST NOT be modified or deleted until the flush() step is performed!
        series_out.flush()

        del series_out

        return output_path


def _make_arg_parser():
    parser = argparse.ArgumentParser(
        prog="flash2openpmd",
        description="Convert FLASH plot files to openPMD.",
    )
    parser.add_argument("run_directory", help="Directory containing FLASH output")
    parser.add_argument("filename", help="FLASH plot file, e.g. slab_hdf5_plt_cnt_0047")
    parser.add_argument("--fields", default="El_number_density",
                         help="yt field name (default: %(default)s)")
    parser.add_argument("--level", type=int, default=4,
                         help="AMR refinement level (default: %(default)s)")
    parser.add_argument("--dtype", choices=["float32", "float64"], default="float32",
                         help="Output array dtype (default: %(default)s)")
    parser.add_argument("--species", default="e",
                         help="Species name in output mesh/filename (default: %(default)s)")
    parser.add_argument("--output-name", default="flash2openpmd",
                         help="Output basename -> <species>_<output-name>_0.h5 (default: %(default)s)")
    parser.add_argument("--geometry", choices=["auto", "cartesian", "cylindrical"], default="auto",
                         help="Override detected geometry (default: %(default)s)")
    parser.add_argument("--device", choices=["auto", "cuda", "cpu"], default="auto",
                         help="Interpolation device (default: %(default)s)")
    parser.add_argument("--normalize", action=argparse.BooleanOptionalAction, default=True,
                         help="Divide by global max (default: enabled)")
    parser.add_argument("--threshold", type=float, default=1e-4,
                         help="Zero out values <= threshold after normalization; 0 disables (default: %(default)s)")
    parser.add_argument("--author", default="Your Name <Your@email>",
                         help="openPMD author attribute")
    parser.add_argument("--preview", action=argparse.BooleanOptionalAction, default=True,
                         help="Show plot before writing (default: enabled)")
    parser.add_argument("--yes", "-y", action="store_true",
                         help="Write without confirmation prompt")
    parser.add_argument("--log-scale", action=argparse.BooleanOptionalAction, default=True,
                         help="Preview with LogNorm color scale (default: enabled)")
    parser.add_argument("--cmap", default="jet",
                         help="Preview colormap (default: %(default)s)")
    return parser


def main(argv=None):
    args = _make_arg_parser().parse_args(argv)

    device = resolve_device(args.device)
    print(f"Using device: {device}")

    input_path = os.path.join(args.run_directory, args.filename)
    print(f"Input file size (before refine): {format_bytes(os.path.getsize(input_path))}")

    ts = Convert(args.run_directory, args.filename)
    geometry = None if args.geometry == "auto" else args.geometry
    density = ts.get_data(fields=args.fields, level=args.level, dtype=args.dtype,
                           device=device, geometry=geometry)

    print(f"Refined data size (after refine, level={args.level}): {format_bytes(density.nbytes)}")

    if args.normalize:
        density = density / density.max()
    if args.threshold > 0:
        density[density <= args.threshold] = 0

    output_path = os.path.join(args.run_directory, f"{args.species}_{args.output_name}_0.h5")

    if args.preview:
        try:
            import matplotlib.pyplot as plt

            sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib", "python"))
            from plotting import make_preview_figure

            fig = plt.figure(figsize=(8, 6))
            make_preview_figure(density, log_scale=args.log_scale, cmap=args.cmap,
                                 title=args.filename, fig=fig)
            plt.show(block=True)
        except Exception as exc:
            print(f"Preview unavailable ({exc}); continuing without it.", file=sys.stderr)

    if os.path.exists(output_path):
        print(f"Warning: {output_path} already exists and will be overwritten.")

    if not args.yes:
        if not sys.stdin.isatty():
            print("Refusing to write without confirmation on non-interactive stdin; pass -y/--yes.",
                  file=sys.stderr)
            sys.exit(1)
        answer = input(f"Write {output_path}? [y/N] ").strip().lower()
        if not answer.startswith("y"):
            print("Aborted.")
            sys.exit(1)

    written_path = ts.write2openpmd(density, species=args.species, output_name=args.output_name,
                                     author=args.author)
    print(f"Wrote {written_path}")


if __name__ == "__main__":
    main()
