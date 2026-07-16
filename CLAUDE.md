# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Converts FLASH simulation output (HDF5 plot files) to the openPMD format, for use by tools like PIConGPU. Built for ELI-NP, Romania.

## Environment

There is no requirements.txt/pyproject.toml — dependencies are defined only in `misc/environment.yml` (a conda env named `flash2openpmd`: numpy, scipy, matplotlib, yt, h5py, openpmd-api, scienceplots, pytorch).

```
conda env create -f misc/environment.yml
```

`python src/launch.py` does this automatically: it checks whether the `flash2openpmd` conda env exists, offers to create it if missing, and then re-execs into that env's Python to launch the GUI. It also self-heals a known HDF5/h5py version mismatch (conda's hdf5 build vs. the pip-installed openpmd-api wheel) by detecting the "HDF5 library version mismatched" error and offering to run `pip install h5py --upgrade --no-dependencies --force-reinstall`.

**GPU note**: `misc/environment.yml` pins `cuda-version=12`, not the conda-forge default. CUDA-13 PyTorch builds drop `sm_50/60/70` kernel support; on older GPUs (e.g. Pascal-generation, compute capability 6.1) this makes `torch.cuda.is_available()` return `True` while every real kernel launch fails. Don't remove this pin without checking `torch.cuda.get_arch_list()` against the target GPU's compute capability.

There are no tests, linter, or CI configured in this repo.

## Running

```
python src/flash2openpmd.py <run_directory> <filename> [options]   # CLI
python src/launch.py                                                 # GUI (env-aware launcher)
python src/gui.py                                                    # GUI (if env already active)
```

CLI shows a preview plot before writing by default (`--no-preview` to skip, `-y`/`--yes` to skip the write confirmation). `--help` lists all options.

## Architecture

- `src/flash2openpmd.py` — core library (`Convert` class) + CLI entry point (`main()`). No GUI dependencies; this module is imported by both `gui.py` and used standalone.
- `src/gui.py` — Tkinter GUI wrapping the same `Convert` class. Imports `flash2openpmd` and `lib/python/plotting.py` directly (adds `lib/python` to `sys.path`).
- `src/launch.py` — conda-env-aware launcher; not part of the conversion logic itself.
- `lib/python/plotting.py` — preview-figure rendering shared by both the CLI (`--preview`) and the GUI, so both draw the density array identically. Applies `scienceplots` styling at import time (`science`, `ieee`, `no-latex`), then resets `figure.dpi` back to 100 — the `ieee` style sets `figure.dpi=600` for print export, which breaks on-screen/embedded (Tk canvas) rendering if left as-is.
- `misc/environment.yml` — the only dependency manifest.
- `Data/` — sample FLASH input file, used for manual testing (there's no automated test suite).
- `src/output/` — example output images and a sample converted openPMD file, referenced by README.md.

### Conversion pipeline (`Convert` in `flash2openpmd.py`)

1. `read4Flash` — loads the FLASH file via `yt.load()`, builds a `smoothed_covering_grid` at the requested AMR refinement `level` (resolution scales as `2**level` per axis).
2. Geometry branch:
   - **cartesian**: data used as-is.
   - **cylindrical**: `cylindricalRotateAxial` rotates the 2D (r, z) slice through 360° into a 3D Cartesian volume via radial linear interpolation. Two implementations kept in sync:
     - CPU: a per-z-slice `np.interp` loop (the original implementation — kept byte-identical, don't touch without re-verifying numerics).
     - GPU: `_interp_radial_torch`, a vectorized `torch.searchsorted`-based equivalent that chunks over z (halving the chunk size on `torch.cuda.OutOfMemoryError` and retrying).
3. `get_data` normalizes array ordering/dtype and returns the density array.
4. `write2openpmd` writes it out via `openpmd_api` as a single scalar mesh (`<species>_density`) with `axisLabels=['z','y','x']`, to `<run_directory>/<species>_<output_name>_0.h5`.

`resolve_device(requested)` is the device-selection entry point used by both CLI and GUI: never trust `torch.cuda.is_available()` alone — it does a real kernel probe (`torch.zeros(1, device="cuda") + 1`) before recommending `"cuda"`, since `is_available()` can be `True` on a GPU whose compute capability isn't covered by the installed PyTorch build (see GPU note above). Falls back to `"cpu"` on any failure.

### Slicing (`plotting.py`)

`make_preview_figure` reduces a 2D or 3D density array to something plottable by slicing along zero or more named axes (`slice_indices: {label: index}`), leaving 1 or 2 remaining dimensions (1D line plot or 2D image respectively). Axis labels are `("z", "y", "x")` for 3D arrays and `("y", "x")` for 2D — this must match `write2openpmd`'s `axisLabels`. `default_slice_indices` reproduces the historical default view (mid-y slice for 3D, full 2D plot for 2D). The GUI's slice panel enforces the "leave 1 or 2 dims" invariant via checkbox enable/disable logic in `gui.py`'s `_update_slice_checkbox_states`.
