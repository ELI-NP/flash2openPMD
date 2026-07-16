#!/usr/bin/env python
"""
Tkinter GUI for the FLASH to openPMD converter: pick a run directory and
file, set conversion parameters, generate a preview plot, and only write
the openPMD file after the preview has been confirmed.
"""

import glob
import os
import sys
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "lib", "python"))

from flash2openpmd import Convert, resolve_device, list_available_fields, format_bytes
from plotting import make_preview_figure, axis_labels_for, default_slice_indices

_ALL_AXES = ("z", "y", "x")


class Flash2OpenPMDGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("FLASH to openPMD Converter")
        self.root.geometry("1200x800")

        self.density = None
        self.density_params = None
        self.file_sizes = (0, 0)

        self._setup_gui()

    # ------------------------------------------------------------------
    # Layout
    # ------------------------------------------------------------------

    def _setup_gui(self):
        main = ttk.Frame(self.root, padding="10")
        main.grid(row=0, column=0, sticky="nsew")
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main.columnconfigure(1, weight=1)
        main.rowconfigure(0, weight=1)

        left = ttk.Frame(main)
        left.grid(row=0, column=0, sticky="ns", padx=(0, 10))

        right = ttk.Frame(main)
        right.grid(row=0, column=1, sticky="nsew")
        right.columnconfigure(0, weight=1)
        right.rowconfigure(0, weight=1)

        self._setup_left_panel(left)
        self._setup_right_panel(right)

        bottom = ttk.Frame(main)
        bottom.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(10, 0))
        bottom.columnconfigure(0, weight=1)

        self.progress = ttk.Progressbar(bottom, mode="indeterminate")
        self.progress.grid(row=0, column=0, sticky="ew", pady=(0, 5))

        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(bottom, textvariable=self.status_var).grid(row=1, column=0, sticky="w")

    def _setup_left_panel(self, parent):
        dir_frame = ttk.LabelFrame(parent, text="Input", padding="10")
        dir_frame.grid(row=0, column=0, sticky="ew", pady=(0, 10))
        dir_frame.columnconfigure(0, weight=1)

        ttk.Label(dir_frame, text="Run directory:").grid(row=0, column=0, sticky="w")
        self.run_dir_var = tk.StringVar()
        ttk.Entry(dir_frame, textvariable=self.run_dir_var, width=30).grid(row=1, column=0, sticky="ew")
        ttk.Button(dir_frame, text="Browse...", command=self._browse_dir).grid(row=1, column=1, padx=(5, 0))

        ttk.Label(dir_frame, text="Filename:").grid(row=2, column=0, sticky="w", pady=(5, 0))
        self.filename_var = tk.StringVar()
        self.filename_combo = ttk.Combobox(dir_frame, textvariable=self.filename_var, width=28)
        self.filename_combo.grid(row=3, column=0, sticky="ew")
        self.filename_combo.bind("<<ComboboxSelected>>", lambda e: self._on_filename_changed())
        ttk.Button(dir_frame, text="Refresh", command=self._refresh_filenames).grid(row=3, column=1, padx=(5, 0))

        params_frame = ttk.LabelFrame(parent, text="Parameters", padding="10")
        params_frame.grid(row=1, column=0, sticky="ew", pady=(0, 10))

        self.fields_var = tk.StringVar(value="El_number_density")
        self.level_var = tk.IntVar(value=4)
        self.dtype_var = tk.StringVar(value="float32")
        self.species_var = tk.StringVar(value="e")
        self.output_name_var = tk.StringVar(value="flash2openpmd")
        self.author_var = tk.StringVar(value="Your Name <Your@email>")
        self.device_var = tk.StringVar(value="auto")
        self.geometry_var = tk.StringVar(value="auto")

        ttk.Label(params_frame, text="Field:").grid(row=0, column=0, sticky="w", pady=2)
        self.fields_combo = ttk.Combobox(params_frame, textvariable=self.fields_var, width=18)
        self.fields_combo.grid(row=0, column=1, sticky="ew", pady=2)
        self._add_labeled_spinbox(params_frame, 1, "Level:", self.level_var, 0, 8)
        self._add_labeled_combo(params_frame, 2, "Dtype:", self.dtype_var, ["float32", "float64"])
        self._add_labeled_entry(params_frame, 3, "Species:", self.species_var)
        self._add_labeled_entry(params_frame, 4, "Output name:", self.output_name_var)
        self._add_labeled_entry(params_frame, 5, "Author:", self.author_var)
        self._add_labeled_combo(params_frame, 6, "Device:", self.device_var, ["auto", "cuda", "cpu"])
        self._add_labeled_combo(params_frame, 7, "Geometry:", self.geometry_var,
                                 ["auto", "cartesian", "cylindrical"])

        post_frame = ttk.LabelFrame(parent, text="Post-processing", padding="10")
        post_frame.grid(row=2, column=0, sticky="ew", pady=(0, 10))
        post_frame.columnconfigure(1, weight=1)

        self.normalize_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(post_frame, text="Normalize to max", variable=self.normalize_var,
                         command=self._invalidate_preview).grid(row=0, column=0, columnspan=2, sticky="w")

        self.threshold_var = tk.StringVar(value="1e-4")
        self._add_labeled_entry(post_frame, 1, "Threshold:", self.threshold_var)

        self.log_scale_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(post_frame, text="Log scale", variable=self.log_scale_var,
                         command=self._on_display_option_changed).grid(row=2, column=0, columnspan=2, sticky="w")

        self.cmap_var = tk.StringVar(value="jet")
        self._add_labeled_combo(post_frame, 3, "Colormap:", self.cmap_var,
                                 ["jet", "viridis", "plasma", "inferno"], on_change=self._on_display_option_changed)

        for var in (self.fields_var, self.species_var, self.output_name_var, self.author_var, self.threshold_var):
            var.trace_add("write", lambda *_: self._invalidate_preview())
        for var in (self.level_var, self.dtype_var, self.device_var, self.geometry_var):
            var.trace_add("write", lambda *_: self._invalidate_preview())

        self._setup_slice_panel(parent)

        button_frame = ttk.Frame(parent)
        button_frame.grid(row=4, column=0, sticky="ew")
        button_frame.columnconfigure(0, weight=1)

        self.preview_btn = ttk.Button(button_frame, text="Generate Preview", command=self._on_generate_preview)
        self.preview_btn.grid(row=0, column=0, sticky="ew", pady=(0, 5))

        self.write_btn = ttk.Button(button_frame, text="Write to openPMD", command=self._on_write, state="disabled")
        self.write_btn.grid(row=1, column=0, sticky="ew")

    def _add_labeled_entry(self, parent, row, label, var):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=2)
        ttk.Entry(parent, textvariable=var, width=20).grid(row=row, column=1, sticky="ew", pady=2)

    def _add_labeled_spinbox(self, parent, row, label, var, lo, hi):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=2)
        ttk.Spinbox(parent, textvariable=var, from_=lo, to=hi, width=18).grid(row=row, column=1, sticky="ew", pady=2)

    def _add_labeled_combo(self, parent, row, label, var, values, on_change=None):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=2)
        combo = ttk.Combobox(parent, textvariable=var, values=values, state="readonly", width=18)
        combo.grid(row=row, column=1, sticky="ew", pady=2)
        if on_change:
            combo.bind("<<ComboboxSelected>>", lambda e: on_change())

    def _setup_slice_panel(self, parent):
        slice_frame = ttk.LabelFrame(parent, text="Slice (2D/3D -> plottable)", padding="10")
        slice_frame.grid(row=3, column=0, sticky="ew", pady=(0, 10))
        slice_frame.columnconfigure(1, weight=1)

        self.slice_vars = {}
        self.slice_widgets = {}

        for row, label in enumerate(_ALL_AXES):
            enabled_var = tk.BooleanVar(value=False)
            position_var = tk.IntVar(value=0)
            self.slice_vars[label] = {"enabled": enabled_var, "position": position_var}

            check = ttk.Checkbutton(slice_frame, text=f"Slice {label}", variable=enabled_var,
                                     command=lambda l=label: self._on_slice_toggle(l))
            check.grid(row=row, column=0, sticky="w")

            scale = ttk.Scale(slice_frame, from_=0, to=1, orient="horizontal",
                               command=lambda v, l=label: self._on_slice_position_changed(l, v))
            scale.grid(row=row, column=1, sticky="ew", padx=(5, 5))

            value_label = ttk.Label(slice_frame, text="-", width=4)
            value_label.grid(row=row, column=2, sticky="w")

            self.slice_widgets[label] = {"check": check, "scale": scale, "value_label": value_label}

        self._set_slice_controls_enabled(False)

    def _setup_right_panel(self, parent):
        self.figure = Figure(figsize=(8, 6))
        ax = self.figure.add_subplot(111)
        ax.text(0.5, 0.5, "No preview yet", ha="center", va="center", transform=ax.transAxes)

        self.canvas = FigureCanvasTkAgg(self.figure, parent)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")

        toolbar_frame = ttk.Frame(parent)
        toolbar_frame.grid(row=1, column=0, sticky="ew")
        toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        toolbar.update()

    # ------------------------------------------------------------------
    # File selection
    # ------------------------------------------------------------------

    def _browse_dir(self):
        directory = filedialog.askdirectory(title="Select FLASH run directory")
        if directory:
            self.run_dir_var.set(directory)
            self._refresh_filenames()

    def _refresh_filenames(self):
        directory = self.run_dir_var.get()
        if not directory or not os.path.isdir(directory):
            return
        matches = sorted(glob.glob(os.path.join(directory, "*_hdf5_plt_cnt_*")))
        names = [os.path.basename(m) for m in matches]
        self.filename_combo["values"] = names
        if names and not self.filename_var.get():
            self.filename_var.set(names[0])
        if names:
            self._on_filename_changed()

    def _on_filename_changed(self):
        run_directory = self.run_dir_var.get()
        filename = self.filename_var.get()
        if not run_directory or not filename:
            return
        self.status_var.set(f"Loading field list for {filename}...")
        threading.Thread(target=self._fields_worker, args=(run_directory, filename), daemon=True).start()

    def _fields_worker(self, run_directory, filename):
        try:
            fields = list_available_fields(run_directory, filename)
            self.root.after(0, self._fields_loaded, fields)
        except Exception as exc:
            self.root.after(0, self.status_var.set, f"Field list unavailable: {exc}")

    def _fields_loaded(self, fields):
        self.fields_combo["values"] = fields
        if fields and self.fields_var.get() not in fields:
            self.fields_var.set(fields[0])
        self.status_var.set(f"Loaded {len(fields)} fields")

    # ------------------------------------------------------------------
    # Parameter collection / preview invalidation
    # ------------------------------------------------------------------

    def _collect_params(self):
        return dict(
            run_directory=self.run_dir_var.get(),
            filename=self.filename_var.get(),
            fields=self.fields_var.get(),
            level=int(self.level_var.get()),
            dtype=self.dtype_var.get(),
            species=self.species_var.get(),
            output_name=self.output_name_var.get(),
            author=self.author_var.get(),
            device=self.device_var.get(),
            geometry=self.geometry_var.get(),
            normalize=self.normalize_var.get(),
            threshold=float(self.threshold_var.get() or 0),
            log_scale=self.log_scale_var.get(),
            cmap=self.cmap_var.get(),
        )

    def _invalidate_preview(self):
        self.density = None
        self.write_btn["state"] = "disabled"
        self._set_slice_controls_enabled(False)

    def _on_display_option_changed(self):
        self._redraw_preview()

    def _redraw_preview(self):
        if self.density is None:
            return
        make_preview_figure(self.density, log_scale=self.log_scale_var.get(),
                             cmap=self.cmap_var.get(), title=self.density_params.get("filename", ""),
                             fig=self.figure, slice_indices=self._current_slice_indices())
        self.canvas.draw()

    # ------------------------------------------------------------------
    # Slice controls
    # ------------------------------------------------------------------

    def _set_slice_controls_enabled(self, enabled):
        for label in _ALL_AXES:
            state = "normal" if enabled else "disabled"
            self.slice_widgets[label]["check"]["state"] = state
            self.slice_widgets[label]["scale"]["state"] = "disabled"

    def _configure_slice_controls(self, density):
        labels = axis_labels_for(density.ndim)
        defaults = default_slice_indices(density)

        for label in _ALL_AXES:
            widgets = self.slice_widgets[label]
            if label not in labels:
                self.slice_vars[label]["enabled"].set(False)
                widgets["check"]["state"] = "disabled"
                widgets["scale"]["state"] = "disabled"
                widgets["value_label"]["text"] = "-"
                continue

            dim = labels.index(label)
            size = density.shape[dim]
            widgets["scale"].configure(from_=0, to=max(size - 1, 0))

            default_pos = defaults.get(label, size // 2)
            self.slice_vars[label]["position"].set(default_pos)
            self.slice_vars[label]["enabled"].set(label in defaults)
            widgets["value_label"]["text"] = str(default_pos)

        self._update_slice_checkbox_states()

    def _update_slice_checkbox_states(self):
        if self.density is None:
            return
        labels = axis_labels_for(self.density.ndim)
        max_checked = self.density.ndim - 1
        checked_count = sum(1 for l in labels if self.slice_vars[l]["enabled"].get())

        for label in _ALL_AXES:
            widgets = self.slice_widgets[label]
            if label not in labels:
                continue
            is_checked = self.slice_vars[label]["enabled"].get()
            if checked_count >= max_checked and not is_checked:
                widgets["check"]["state"] = "disabled"
            else:
                widgets["check"]["state"] = "normal"
            widgets["scale"]["state"] = "normal" if is_checked else "disabled"

    def _current_slice_indices(self):
        if self.density is None:
            return {}
        labels = axis_labels_for(self.density.ndim)
        return {label: self.slice_vars[label]["position"].get()
                for label in labels if self.slice_vars[label]["enabled"].get()}

    def _on_slice_toggle(self, label):
        self._update_slice_checkbox_states()
        self._redraw_preview()

    def _on_slice_position_changed(self, label, value):
        idx = int(round(float(value)))
        self.slice_vars[label]["position"].set(idx)
        self.slice_widgets[label]["value_label"]["text"] = str(idx)
        if self.slice_vars[label]["enabled"].get():
            self._redraw_preview()

    # ------------------------------------------------------------------
    # Busy/idle state
    # ------------------------------------------------------------------

    def _set_busy(self, message):
        self.preview_btn["state"] = "disabled"
        self.write_btn["state"] = "disabled"
        self.progress.start(10)
        self.status_var.set(message)

    def _set_idle(self, message):
        self.progress.stop()
        self.preview_btn["state"] = "normal"
        self.status_var.set(message)

    # ------------------------------------------------------------------
    # Preview
    # ------------------------------------------------------------------

    def _on_generate_preview(self):
        try:
            params = self._collect_params()
        except ValueError as exc:
            messagebox.showerror("Error", f"Invalid parameter: {exc}")
            return

        if not params["run_directory"] or not params["filename"]:
            messagebox.showwarning("Warning", "Please select a run directory and filename.")
            return

        self._set_busy("Loading FLASH data (yt)...")
        threading.Thread(target=self._preview_worker, args=(params,), daemon=True).start()

    def _preview_worker(self, params):
        try:
            device = resolve_device(params["device"])
            ts = Convert(params["run_directory"], params["filename"])
            geometry = None if params["geometry"] == "auto" else params["geometry"]

            input_path = os.path.join(params["run_directory"], params["filename"])
            input_size = os.path.getsize(input_path)
            print(f"Input file size (before refine): {format_bytes(input_size)}")

            self.root.after(0, self.status_var.set, f"Interpolating on {device}...")
            density = ts.get_data(fields=params["fields"], level=params["level"], dtype=params["dtype"],
                                   device=device, geometry=geometry)

            print(f"Refined data size (after refine, level={params['level']}): {format_bytes(density.nbytes)}")
            self.file_sizes = (input_size, density.nbytes)

            if params["normalize"]:
                density = density / density.max()
            if params["threshold"] > 0:
                density[density <= params["threshold"]] = 0

            self.root.after(0, self._preview_done, density, params)
        except Exception as exc:
            self.root.after(0, self._on_error, str(exc))

    def _preview_done(self, density, params):
        self.density = density
        self.density_params = params

        self._set_slice_controls_enabled(True)
        self._configure_slice_controls(density)
        self._redraw_preview()

        before_size, after_size = self.file_sizes
        self._set_idle(f"Preview ready - shape={density.shape}, max={density.max():.3e}  |  "
                        f"file size before refine: {format_bytes(before_size)}, "
                        f"after refine: {format_bytes(after_size)}")
        self.write_btn["state"] = "normal"

    # ------------------------------------------------------------------
    # Write
    # ------------------------------------------------------------------

    def _on_write(self):
        if self.density is None:
            messagebox.showwarning("Warning", "Generate a preview first.")
            return

        params = self.density_params
        output_path = os.path.join(params["run_directory"], f"{params['species']}_{params['output_name']}_0.h5")

        if os.path.exists(output_path):
            prompt = f"{output_path} already exists and will be overwritten. Continue?"
        else:
            prompt = f"Write {output_path}?"
        if not messagebox.askyesno("Confirm", prompt):
            return

        self._set_busy(f"Writing {output_path}...")
        threading.Thread(target=self._write_worker, args=(params,), daemon=True).start()

    def _write_worker(self, params):
        try:
            ts = Convert(params["run_directory"], params["filename"])
            output_path = ts.write2openpmd(self.density, species=params["species"],
                                            output_name=params["output_name"], author=params["author"])
            self.root.after(0, self._write_done, output_path)
        except Exception as exc:
            self.root.after(0, self._on_error, str(exc))

    def _write_done(self, output_path):
        self._set_idle(f"Wrote {output_path}")
        self.write_btn["state"] = "normal"
        messagebox.showinfo("Success", f"Wrote {output_path}")

    def _on_error(self, message):
        self._set_idle("Error")
        self.write_btn["state"] = "disabled" if self.density is None else "normal"
        messagebox.showerror("Error", message)


def main():
    root = tk.Tk()
    Flash2OpenPMDGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
