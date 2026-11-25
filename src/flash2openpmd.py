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


CHARGE_STATE_FIELD = "charge_state"

def round_charge_state(zbar):
    """
    Integer charge state from the mean ionisation Zbar, rounded to the
    nearest integer with halves rounded up (0.4 -> 0, 0.5 -> 1, 2.5 -> 3).
    """

    zbar = np.asarray(zbar)
    return np.floor(zbar + 0.5).clip(min=0).astype(zbar.dtype)

# k_B / (1e3 e): multiplies a temperature in kelvin to give it in keV.
KELVIN_TO_KEV = sc.k / (1e3 * sc.e)


def add_charge_state_field(ds):
    """
    Register the mean ionisation (charge state) Zbar as the derived field
    ("gas", CHARGE_STATE_FIELD), if the FLASH file has what it needs.

    A FLASH "zbar" variable is used directly when present. Otherwise Zbar is
    derived from the mass-scalar variables ye = Zbar/Abar and sumy = 1/Abar
    as Zbar = ye / sumy (0 where sumy <= 0). The field is the continuous
    Zbar; read4Flash turns it into an integer after refinement (see
    round_charge_state).

    Returns True if the field was registered, False if the file has neither
    "zbar" nor both "ye" and "sumy".
    """

    on_disk = set(ds.field_list)

    if ("flash", "zbar") in on_disk:
        def _charge_state(field, data):
            return data.ds.arr(data["flash", "zbar"].d, "")
    elif ("flash", "ye") in on_disk and ("flash", "sumy") in on_disk:
        def _charge_state(field, data):
            ye = data["flash", "ye"].d
            sumy = data["flash", "sumy"].d
            return data.ds.arr(np.divide(ye, sumy, out=np.zeros_like(ye), where=sumy > 0), "")
    else:
        return False

    ds.add_field(("gas", CHARGE_STATE_FIELD), function=_charge_state,
                 sampling_type="cell", units="", force_override=True)
    return True


def _field_key(ds, field):
    """
    The yt field tuple for `field`: the "gas" field of that name if yt has
    one, else the raw on-disk FLASH variable (e.g. "tele", "ye").
    """

    if ("gas", field) in ds.derived_field_list:
        return ("gas", field)
    if ("flash", field) in ds.field_list:
        return ("flash", field)
    raise ValueError(f"Unknown field {field!r} in {ds.basename}")


def temperature_to_kev_factor(ds, key):
    """
    The factor converting field `key` (a yt field tuple) from its native
    units to keV if it is a temperature (e.g. FLASH "tele", "tion", "trad",
    "temp"), else None.
    """

    units = ds.field_info[key].units
    if not units:
        return None
    quantity = ds.quan(1.0, units)
    if quantity.units.dimensions != yt.units.dimensions.temperature:
        return None
    return float(quantity.to("K").d) * KELVIN_TO_KEV


def list_available_fields(run_directory, filename):
    """
    The "gas"-typed field names yt exposes for this FLASH file (on-disk and
    derived), the raw FLASH variables, and CHARGE_STATE_FIELD when it can
    be derived; suitable for the --fields/Field parameter.
    """

    ds = yt.load(os.path.join(run_directory, filename))
    names = {name for ftype, name in ds.derived_field_list if ftype == "gas"}
    names.update(name for ftype, name in ds.field_list if ftype == "flash")
    if add_charge_state_field(ds):
        names.add(CHARGE_STATE_FIELD)
    return sorted(names)


def parse_field_specs(specs):
    """
    Parse field specifications into an ordered list of (field, mesh_name)
    pairs, mesh_name being None when the field is not renamed.

    `specs` is a string or an iterable of strings; each string may hold
    several comma-separated items, and each item is "field" or
    "field=mesh_name". Duplicate fields are dropped (first one wins).
    """

    if isinstance(specs, str):
        specs = [specs]

    parsed = []
    seen = set()
    for spec in specs:
        for item in spec.split(","):
            item = item.strip()
            if not item:
                continue
            field, sep, name = item.partition("=")
            field, name = field.strip(), name.strip()
            if not field or (sep and not name):
                raise ValueError(f"Invalid field specification {item!r}; use 'field' or 'field=name'")
            if field in seen:
                continue
            seen.add(field)
            parsed.append((field, name or None))

    if not parsed:
        raise ValueError("No fields given")
    return parsed


def resolve_mesh_names(specs, species):
    """
    Map each field of `specs` (from parse_field_specs) to its openPMD mesh
    name: the rename if given, else the field name. A single field without
    a rename keeps the historical name "<species>_density".
    """

    if len(specs) == 1 and specs[0][1] is None:
        return {specs[0][0]: f"{species}_density"}

    mesh_names = {field: (name or field) for field, name in specs}
    names = list(mesh_names.values())
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        raise ValueError(f"Duplicate mesh name(s): {', '.join(duplicates)}")
    return mesh_names


def is_density_field(field):
    """Whether `field` is normalized/thresholded by apply_normalization."""

    return "density" in field.lower()


def apply_normalization(data, normalize=True, threshold=0.0):
    """
    Normalize and threshold the density fields of `data` ({field: array}),
    in place, returning `data`.

    All density fields (see is_density_field) share one divisor, so their
    ratios are preserved: `normalize` if it is a number, else the maximum
    of the first density field when it is True; False disables it. Values
    <= threshold are then zeroed in the density fields (threshold <= 0
    disables it). Other fields, e.g. CHARGE_STATE_FIELD, are untouched.
    """

    density_fields = [f for f in data if is_density_field(f)]
    if not density_fields:
        return data

    if normalize is not False and normalize is not None:
        divisor = data[density_fields[0]].max() if normalize is True else normalize
        if divisor == 0:
            raise ValueError(f"Cannot normalize: divisor is 0 (max of {density_fields[0]})")
        for field in density_fields:
            data[field] = data[field] / divisor

    if threshold > 0:
        for field in density_fields:
            data[field][data[field] <= threshold] = 0

    return data


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
        fields : string or list of strings
        	The field data (unit of FLASH code), i.e. fields="El_number_density".
        	A list reads several fields from the same refined grid.
        	CHARGE_STATE_FIELD ("charge_state") is available when the
        	file allows it (see add_charge_state_field) and is made an
        	integer after refinement (see round_charge_state). Temperature fields are
        	converted from K to keV.

        level : int
            The level of refinement using yt-project.

        geometry : string or None
            Override the geometry detected by yt ("cartesian" or "cylindrical").
            None (default) uses the detected geometry.

        device : string
            "cpu" or "cuda". Only affects the cylindrical rotation step.

        Returns
        -------
        The raw field array (a {field: array} dict when `fields` is a list,
        temperatures in keV), the dimensionality and the geometry.

        """

        single = isinstance(fields, str)
        field_list = [fields] if single else list(fields)

        ds = yt.load(self.path_to_sim)

        has_charge_state = add_charge_state_field(ds)
        if CHARGE_STATE_FIELD in field_list and not has_charge_state:
            raise ValueError(f"{CHARGE_STATE_FIELD} is unavailable: the FLASH file has "
                             "neither 'zbar' nor both 'ye' and 'sumy'")
        keys = {field: _field_key(ds, field) for field in field_list}
        kev_factors = {field: temperature_to_kev_factor(ds, key) for field, key in keys.items()}

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

        data = {}
        for field in field_list:
            print(f"Reading field: {field}")
            if geom == 'cylindrical':
                data[field] = self.cylindricalRotateAxial(all_data,fields=keys[field],device=device)
            else:
                data[field] = np.asarray(all_data[keys[field]])
            if kev_factors[field] is not None:
                print(f"Converting {field} from K to keV")
                data[field] = data[field] * kev_factors[field]
            if field == CHARGE_STATE_FIELD:
                print(f"Rounding {field} to the nearest integer")
                data[field] = round_charge_state(data[field])

        return (data[fields] if single else data), ND, geom
    
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

        if isinstance(fields, str):
            fields = ('gas', fields)
        rho = array[fields]
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

    @staticmethod
    def _orient(raw, ND, geom, dtype=None):
        """Reorder a raw field array into the C-ordered (z,)y,x layout written to openPMD."""

        if ND == 2:
            if geom == 'cylindrical':
                return np.array(np.rot90(raw), order='C', dtype=dtype)
            return np.array(raw[:,:,0].T, order='C', dtype=dtype)
        return np.array(raw.T, order='C', dtype=dtype)

    def get_fields(self,fields,level,dtype=None,device="cpu",geometry=None):
        """
        Read several fields from the same refined grid.

        Returns an ordered {field: array} dict, each array laid out as
        get_data returns it.
        """

        raw, ND, geom = self.read4Flash(list(fields),level,geometry=geometry,device=device)
        return {field: self._orient(arr, ND, geom, dtype) for field, arr in raw.items()}

    def get_data(self,fields,level,dtype=None,device="cpu",geometry=None):

        return self.get_fields([fields],level,dtype=dtype,device=device,geometry=geometry)[fields]

    def write2openpmd(self,data,species,output_name,author="Your Name <Your@email>",mesh_names=None):
        """
        Write to <run_directory>/<species>_<output_name>_0.h5.

        `data` is either a single array, written as the mesh
        "<species>_density", or a {field: array} dict, each entry written
        as its own scalar mesh in the same iteration. `mesh_names` maps
        fields to mesh names (default: the field names).
        """

        if isinstance(data, np.ndarray):
            meshes = {f"{species}_density": data}
        else:
            mesh_names = mesh_names or {}
            meshes = {mesh_names.get(field, field): arr for field, arr in data.items()}
            if len(meshes) != len(data):
                raise ValueError("Duplicate mesh names in mesh_names")

        os.makedirs(self.run_directory, exist_ok=True)
        output_path = os.path.join(self.run_directory, f"{species}_{output_name}_0.h5")

        series_out = io.Series(output_path, io.Access.create)

        k = series_out.iterations[0]

        series_out.author = author

        for mesh_name, array in meshes.items():
            # record - again,important to specify as scalar
            mesh = k.meshes[mesh_name]
            mrc = mesh[io.Mesh_Record_Component.SCALAR]

            mrc.reset_dataset(io.Dataset(array.dtype, array.shape))

            axis_labels = ['z', 'y', 'x'] if array.ndim == 3 else ['y', 'x']

            mesh.set_attribute('dataOrder','C')
            mesh.set_attribute('axisLabels', axis_labels)
            mrc.store_chunk(array)
        # After registering a data chunk,
        # it MUST NOT be modified or deleted until the flush() step is performed!
        series_out.flush()

        del series_out

        return output_path


class _NormalizeAction(argparse.Action):
    """--normalize [VALUE] enables normalization (dividing by VALUE, or by the
    global max if VALUE is omitted); --no-normalize disables it."""

    def __init__(self, option_strings, dest, default=True, help=None):
        opts = list(option_strings)
        for opt in option_strings:
            if opt.startswith("--"):
                opts.append("--no-" + opt[2:])
        super().__init__(opts, dest, nargs="?", type=float, default=default, help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        if option_string.startswith("--no-"):
            setattr(namespace, self.dest, False)
        else:
            setattr(namespace, self.dest, True if values is None else values)


def _make_arg_parser():
    parser = argparse.ArgumentParser(
        prog="flash2openpmd",
        description="Convert FLASH plot files to openPMD.",
    )
    parser.add_argument("run_directory", help="Directory containing FLASH output")
    parser.add_argument("filename", help="FLASH plot file, e.g. slab_hdf5_plt_cnt_0047")
    parser.add_argument("--fields", nargs="+", default=["El_number_density"], metavar="FIELD",
                         help="yt field name(s), space- or comma-separated, all written to one file "
                              "on the same grid. FIELD=NAME sets the mesh name (default: the field "
                              "name; a single unrenamed field is written as <species>_density). "
                              "Default: El_number_density")
    parser.add_argument("--charge-state", action="store_true",
                         help=f"Also write the mean ionisation, rounded to the nearest integer, as "
                              f"field '{CHARGE_STATE_FIELD}' (FLASH zbar, or ye/sumy); skipped with a warning if unavailable")
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
    parser.add_argument("--normalize", action=_NormalizeAction, default=True,
                         help="Divide all density fields by the max of the first density field, "
                              "or by VALUE if given; other fields (e.g. charge state) are left "
                              "unchanged, except temperatures, which are written in keV "
                              "(default: enabled, uses max)")
    parser.add_argument("--threshold", type=float, default=1e-4,
                         help="Zero out density-field values <= threshold after normalization; "
                              "0 disables (default: %(default)s)")
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


def charge_state_available(run_directory, filename):
    """Whether CHARGE_STATE_FIELD can be derived for this FLASH file."""

    return add_charge_state_field(yt.load(os.path.join(run_directory, filename)))


def main(argv=None):
    parser = _make_arg_parser()
    args = parser.parse_args(argv)

    try:
        specs = parse_field_specs(args.fields)
    except ValueError as exc:
        parser.error(str(exc))

    if args.charge_state and CHARGE_STATE_FIELD not in [f for f, _ in specs]:
        if charge_state_available(args.run_directory, args.filename):
            specs.append((CHARGE_STATE_FIELD, None))
        else:
            print(f"Warning: charge state unavailable (no 'zbar', or 'ye' and 'sumy', in "
                  f"{args.filename}); skipping {CHARGE_STATE_FIELD}.", file=sys.stderr)

    try:
        mesh_names = resolve_mesh_names(specs, args.species)
    except ValueError as exc:
        parser.error(str(exc))
    fields = [f for f, _ in specs]

    device = resolve_device(args.device)
    print(f"Using device: {device}")

    input_path = os.path.join(args.run_directory, args.filename)
    print(f"Input file size (before refine): {format_bytes(os.path.getsize(input_path))}")

    ts = Convert(args.run_directory, args.filename)
    geometry = None if args.geometry == "auto" else args.geometry
    data = ts.get_fields(fields, level=args.level, dtype=args.dtype,
                         device=device, geometry=geometry)

    total_bytes = sum(arr.nbytes for arr in data.values())
    print(f"Refined data size (after refine, level={args.level}, {len(data)} field(s)): "
          f"{format_bytes(total_bytes)}")

    apply_normalization(data, normalize=args.normalize, threshold=args.threshold)

    output_path = os.path.join(args.run_directory, f"{args.species}_{args.output_name}_0.h5")

    if args.preview:
        try:
            import matplotlib.pyplot as plt

            sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib", "python"))
            from plotting import make_preview_figure

            for field, arr in data.items():
                fig = plt.figure(figsize=(8, 6))
                make_preview_figure(arr, log_scale=args.log_scale, cmap=args.cmap,
                                     title=f"{args.filename}: {field}", fig=fig)
            plt.show(block=True)
        except Exception as exc:
            print(f"Preview unavailable ({exc}); continuing without it.", file=sys.stderr)

    if os.path.exists(output_path):
        print(f"Warning: {output_path} already exists and will be overwritten.")

    print("Meshes: " + ", ".join(f"{field} -> {mesh_names[field]}" for field in fields))

    if not args.yes:
        if not sys.stdin.isatty():
            print("Refusing to write without confirmation on non-interactive stdin; pass -y/--yes.",
                  file=sys.stderr)
            sys.exit(1)
        answer = input(f"Write {output_path}? [y/N] ").strip().lower()
        if not answer.startswith("y"):
            print("Aborted.")
            sys.exit(1)

    written_path = ts.write2openpmd(data, species=args.species, output_name=args.output_name,
                                     author=args.author, mesh_names=mesh_names)
    print(f"Wrote {written_path}")


if __name__ == "__main__":
    main()
