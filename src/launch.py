#!/usr/bin/env python
"""
Conda-environment-aware launcher for the FLASH to openPMD GUI.

Checks whether the `flash2openpmd` conda environment (defined in
misc/environment.yml) exists, offers to create it if not, then runs the
GUI under that environment's Python.
"""

import json
import os
import subprocess
import sys

ENV_NAME = "flash2openpmd"
HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(HERE)
ENV_FILE = os.path.join(REPO_ROOT, "misc", "environment.yml")
GUI_PATH = os.path.join(HERE, "gui.py")


def _conda_envs():
    """Return {env_name: env_path} for existing conda environments, or {} if conda isn't available."""
    try:
        result = subprocess.run(["conda", "env", "list", "--json"],
                                 capture_output=True, text=True, check=True)
    except (FileNotFoundError, subprocess.CalledProcessError):
        return {}
    return {os.path.basename(p): p for p in json.loads(result.stdout)["envs"]}


def _env_python(env_path):
    for candidate in (os.path.join(env_path, "bin", "python"), os.path.join(env_path, "python.exe")):
        if os.path.exists(candidate):
            return candidate
    return None


def _ask(prompt):
    return input(prompt).strip().lower().startswith("y")


def _fix_hdf5_mismatch(python_bin):
    """
    yt's conda hdf5/h5py build can end up linked against a different
    libhdf5 than the pip-installed openpmd-api wheel, raising an "HDF5
    library version mismatched" error. Detect and offer the known fix.
    """

    check = subprocess.run([python_bin, "-c", "import openpmd_api"], capture_output=True, text=True)
    if check.returncode == 0:
        return
    if "HDF5" not in check.stderr or "mismatch" not in check.stderr.lower():
        return

    print("Detected an HDF5 library version mismatch between conda's h5py and openpmd-api.")
    if not _ask("Fix by running 'pip install h5py --upgrade --no-dependencies --force-reinstall'? [y/N] "):
        print("Skipping fix; the GUI may fail to import openpmd_api.")
        return

    subprocess.run([python_bin, "-m", "pip", "install", "h5py",
                     "--upgrade", "--no-dependencies", "--force-reinstall"], check=True)


def ensure_env():
    envs = _conda_envs()

    if ENV_NAME not in envs:
        print(f"Conda environment '{ENV_NAME}' not found.")
        if not _ask(f"Create it now from {ENV_FILE}? [y/N] "):
            print(f"Cannot continue without the environment. Create it manually with:\n"
                  f"  conda env create -f {ENV_FILE}")
            return None

        subprocess.run(["conda", "env", "create", "-f", ENV_FILE], check=True)
        envs = _conda_envs()
        if ENV_NAME not in envs:
            print("Environment creation did not produce the expected environment name.")
            return None

        python_bin = _env_python(envs[ENV_NAME])
        _fix_hdf5_mismatch(python_bin)
        return python_bin

    python_bin = _env_python(envs[ENV_NAME])
    _fix_hdf5_mismatch(python_bin)
    return python_bin


def main():
    python_bin = ensure_env()
    if python_bin is None:
        sys.exit(1)

    if os.path.abspath(python_bin) != os.path.abspath(sys.executable):
        os.execv(python_bin, [python_bin, GUI_PATH])

    sys.path.insert(0, HERE)
    from gui import main as gui_main
    gui_main()


if __name__ == "__main__":
    main()
