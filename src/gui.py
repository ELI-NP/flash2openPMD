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

from flash2openpmd import (Convert, resolve_device, list_available_fields, format_bytes,
                           resolve_mesh_names, apply_normalization, is_density_field,
                           CHARGE_STATE_FIELD)
from plotting import make_preview_figure, axis_labels_for, default_slice_indices

_ALL_AXES = ("z", "y", "x")
_DEFAULT_FIELD = "El_number_density"


class Flash2OpenPMDGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("FLASH to openPMD Converter")
        self.root.geometry("1200x900")

        self.data = None           # {field: array} from the last preview
        self.mesh_names = None     # {field: mesh name} for the last preview
        self.density = None        # the array currently displayed (one of self.data)
        self.density_params = None
        self.file_sizes = (0, 0)

        self.available_fields = [_DEFAULT_FIELD]
        self.field_rows = []       # one {field_var, name_var, widgets...} per field row

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
        left.rowconfigure(0, weight=1)

        right = ttk.Frame(main)
        right.grid(row=0, column=1, sticky="nsew")
        right.columnconfigure(0, weight=1)
        right.rowconfigure(0, weight=1)

        self._setup_left_panel(self._make_scrollable(left), left)
        self._setup_right_panel(right)

        bottom = ttk.Frame(main)
        bottom.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(10, 0))
        bottom.columnconfigure(0, weight=1)

        self.progress = ttk.Progressbar(bottom, mode="indeterminate")
        self.progress.grid(row=0, column=0, sticky="ew", pady=(0, 5))

        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(bottom, textvariable=self.status_var).grid(row=1, column=0, sticky="w")

    def _make_scrollable(self, parent):
        """
        A frame inside a vertically scrolling canvas at row 0 of `parent`, so
        the settings stay reachable however many field rows are added.
        """

        background = ttk.Style().lookup("TFrame", "background") or None
        canvas = tk.Canvas(parent, highlightthickness=0, borderwidth=0, background=background)
        scrollbar = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.grid(row=0, column=0, sticky="ns")
        scrollbar.grid(row=0, column=1, sticky="ns")

        inner = ttk.Frame(canvas)
        canvas.create_window((0, 0), window=inner, anchor="nw")
        inner.bind("<Configure>", lambda e: canvas.configure(
            scrollregion=canvas.bbox("all"), width=inner.winfo_reqwidth()))

        def scroll(event):
            if event.num == 4 or getattr(event, "delta", 0) > 0:
                canvas.yview_scroll(-1, "units")
            else:
                canvas.yview_scroll(1, "units")

        def bind_wheel(_):
            for seq in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
                canvas.bind_all(seq, scroll)

        def unbind_wheel(_):
            for seq in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
                canvas.unbind_all(seq)

        inner.bind("<Enter>", bind_wheel)
        inner.bind("<Leave>", unbind_wheel)
        return inner

    def _setup_left_panel(self, parent, button_parent):
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

        self._setup_fields_panel(parent)

        params_frame = ttk.LabelFrame(parent, text="Parameters", padding="10")
        params_frame.grid(row=2, column=0, sticky="ew", pady=(0, 10))

        self.level_var = tk.IntVar(value=4)
        self.dtype_var = tk.StringVar(value="float32")
        self.species_var = tk.StringVar(value="e")
        self.output_name_var = tk.StringVar(value="flash2openpmd")
        self.author_var = tk.StringVar(value="Your Name <Your@email>")
        self.device_var = tk.StringVar(value="auto")
        self.geometry_var = tk.StringVar(value="auto")

        self._add_labeled_spinbox(params_frame, 0, "Level:", self.level_var, 0, 8)
        self._add_labeled_combo(params_frame, 1, "Dtype:", self.dtype_var, ["float32", "float64"])
        self._add_labeled_entry(params_frame, 2, "Species:", self.species_var)
        self._add_labeled_entry(params_frame, 3, "Output name:", self.output_name_var)
        self._add_labeled_entry(params_frame, 4, "Author:", self.author_var)
        self._add_labeled_combo(params_frame, 5, "Device:", self.device_var, ["auto", "cuda", "cpu"])
        self._add_labeled_combo(params_frame, 6, "Geometry:", self.geometry_var,
                                 ["auto", "cartesian", "cylindrical"])

        post_frame = ttk.LabelFrame(parent, text="Post-processing", padding="10")
        post_frame.grid(row=3, column=0, sticky="ew", pady=(0, 10))
        post_frame.columnconfigure(1, weight=1)

        self.normalize_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(post_frame, text="Normalize (density fields only)", variable=self.normalize_var,
                         command=self._invalidate_preview).grid(row=0, column=0, columnspan=2, sticky="w")

        self.normalize_value_var = tk.StringVar(value="max")
        self._add_labeled_entry(post_frame, 1, "Normalize to:", self.normalize_value_var)

        self.threshold_var = tk.StringVar(value="1e-4")
        self._add_labeled_entry(post_frame, 2, "Threshold:", self.threshold_var)

        self.density_hint_var = tk.StringVar(value="")
        ttk.Label(post_frame, textvariable=self.density_hint_var, foreground="gray",
                  wraplength=240, justify="left").grid(row=3, column=0, columnspan=2, sticky="w", pady=(0, 4))

        self.log_scale_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(post_frame, text="Log scale", variable=self.log_scale_var,
                         command=self._on_display_option_changed).grid(row=4, column=0, columnspan=2, sticky="w")

        self.cmap_var = tk.StringVar(value="jet")
        self._add_labeled_combo(post_frame, 5, "Colormap:", self.cmap_var,
                                 ["jet", "viridis", "plasma", "inferno"], on_change=self._on_display_option_changed)

        self.show_field_var = tk.StringVar(value="")
        ttk.Label(post_frame, text="Show field:").grid(row=6, column=0, sticky="w", pady=2)
        self.show_field_combo = ttk.Combobox(post_frame, textvariable=self.show_field_var,
                                              state="disabled", width=18)
        self.show_field_combo.grid(row=6, column=1, sticky="ew", pady=2)
        self.show_field_combo.bind("<<ComboboxSelected>>", lambda e: self._on_show_field_changed())

        self._update_density_hint()

        for var in (self.species_var, self.output_name_var, self.author_var,
                    self.normalize_value_var, self.threshold_var):
            var.trace_add("write", lambda *_: self._invalidate_preview())
        for var in (self.level_var, self.dtype_var, self.device_var, self.geometry_var):
            var.trace_add("write", lambda *_: self._invalidate_preview())

        self._setup_slice_panel(parent)

        button_frame = ttk.Frame(button_parent)
        button_frame.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(10, 0))
        button_frame.columnconfigure(0, weight=1)

        self.preview_btn = ttk.Button(button_frame, text="Generate Preview", command=self._on_generate_preview)
        self.preview_btn.grid(row=0, column=0, sticky="ew", pady=(0, 5))

        self.write_btn = ttk.Button(button_frame, text="Write to openPMD", command=self._on_write, state="disabled")
        self.write_btn.grid(row=1, column=0, sticky="ew")

    def _setup_fields_panel(self, parent):
        fields_frame = ttk.LabelFrame(parent, text="Fields", padding="10")
        fields_frame.grid(row=1, column=0, sticky="ew", pady=(0, 10))
        fields_frame.columnconfigure(0, weight=1)

        header = ttk.Frame(fields_frame)
        header.grid(row=0, column=0, sticky="ew")
        ttk.Label(header, text="Field", width=20).grid(row=0, column=0, sticky="w")
        ttk.Label(header, text="Mesh name (optional)").grid(row=0, column=1, sticky="w", padx=(5, 0))

        self.field_rows_frame = ttk.Frame(fields_frame)
        self.field_rows_frame.grid(row=1, column=0, sticky="ew")
        self.field_rows_frame.columnconfigure(0, weight=1)

        ttk.Button(fields_frame, text="+ Add field", command=self._add_field_row).grid(
            row=2, column=0, sticky="w", pady=(4, 0))

        self.charge_state_var = tk.BooleanVar(value=False)
        self.charge_state_check = ttk.Checkbutton(
            fields_frame, text=f"Charge state ({CHARGE_STATE_FIELD})", variable=self.charge_state_var,
            command=self._on_fields_changed)
        self.charge_state_check.grid(row=3, column=0, sticky="w", pady=(4, 0))

        self._add_field_row(_DEFAULT_FIELD, notify=False)

    def _add_field_row(self, field=None, notify=True):
        """Append a field row: a field dropdown, an optional mesh name and a remove button."""

        if field is None:
            used = set(self._selected_fields())
            field = next((f for f in self.available_fields if f not in used), self.available_fields[0])

        field_var = tk.StringVar(value=field)
        name_var = tk.StringVar(value="")
        combo = ttk.Combobox(self.field_rows_frame, textvariable=field_var,
                              values=self.available_fields, state="readonly", width=20)
        combo.bind("<<ComboboxSelected>>", lambda e: self._on_fields_changed())
        name_entry = ttk.Entry(self.field_rows_frame, textvariable=name_var, width=12)
        name_var.trace_add("write", lambda *_: self._invalidate_preview())

        row = {"field_var": field_var, "name_var": name_var, "combo": combo, "name_entry": name_entry}
        row["remove_btn"] = ttk.Button(self.field_rows_frame, text="X", width=2,
                                        command=lambda: self._remove_field_row(row))
        self.field_rows.append(row)
        self._regrid_field_rows()
        if notify:
            self._on_fields_changed()

    def _remove_field_row(self, row):
        if len(self.field_rows) <= 1:
            return
        self.field_rows.remove(row)
        for key in ("combo", "name_entry", "remove_btn"):
            row[key].destroy()
        self._regrid_field_rows()
        self._on_fields_changed()

    def _regrid_field_rows(self):
        for i, row in enumerate(self.field_rows):
            row["combo"].grid(row=i, column=0, sticky="ew", pady=1)
            row["name_entry"].grid(row=i, column=1, sticky="ew", padx=(5, 0), pady=1)
            row["remove_btn"].grid(row=i, column=2, padx=(5, 0), pady=1)
            row["remove_btn"]["state"] = "normal" if len(self.field_rows) > 1 else "disabled"

    def _on_fields_changed(self):
        self._update_density_hint()
        self._invalidate_preview()

    def _update_density_hint(self):
        """Show which of the chosen fields Normalize/Threshold act on."""

        fields = self._selected_fields()
        if self.charge_state_var.get() and CHARGE_STATE_FIELD not in fields:
            fields.append(CHARGE_STATE_FIELD)
        density = [f for f in dict.fromkeys(fields) if is_density_field(f)]
        if density:
            text = "Normalize and threshold apply only to: " + ", ".join(density)
        else:
            text = "No density field chosen: normalize and threshold have no effect."
        self.density_hint_var.set(text + "\nOther fields (charge state, temperatures in keV) are not normalized.")

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
        slice_frame.grid(row=4, column=0, sticky="ew", pady=(0, 10))
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
        if fields:
            self.available_fields = fields
        fallback = _DEFAULT_FIELD if _DEFAULT_FIELD in self.available_fields else self.available_fields[0]
        for row in self.field_rows:
            row["combo"]["values"] = self.available_fields
            if row["field_var"].get() not in self.available_fields:
                row["field_var"].set(fallback)

        has_charge_state = CHARGE_STATE_FIELD in fields
        self.charge_state_check["state"] = "normal" if has_charge_state else "disabled"
        if not has_charge_state:
            self.charge_state_var.set(False)

        self._on_fields_changed()
        self.status_var.set(f"Loaded {len(fields)} fields" +
                            ("" if has_charge_state else " (charge state unavailable)"))

    def _selected_fields(self):
        return [row["field_var"].get() for row in self.field_rows if row["field_var"].get()]

    # ------------------------------------------------------------------
    # Parameter collection / preview invalidation
    # ------------------------------------------------------------------

    def _collect_params(self):
        return dict(
            run_directory=self.run_dir_var.get(),
            filename=self.filename_var.get(),
            field_specs=self._field_specs(),
            level=int(self.level_var.get()),
            dtype=self.dtype_var.get(),
            species=self.species_var.get(),
            output_name=self.output_name_var.get(),
            author=self.author_var.get(),
            device=self.device_var.get(),
            geometry=self.geometry_var.get(),
            normalize=self.normalize_var.get(),
            normalize_value=self._parse_normalize_value(),
            threshold=float(self.threshold_var.get() or 0),
            log_scale=self.log_scale_var.get(),
            cmap=self.cmap_var.get(),
        )

    def _field_specs(self):
        """(field, mesh name or None) for each field row, plus charge state if ticked."""

        specs = []
        for row in self.field_rows:
            field = row["field_var"].get()
            if not field:
                continue
            if field in [f for f, _ in specs]:
                raise ValueError(f"field {field!r} is chosen more than once")
            specs.append((field, row["name_var"].get().strip() or None))
        if self.charge_state_var.get() and CHARGE_STATE_FIELD not in [f for f, _ in specs]:
            specs.append((CHARGE_STATE_FIELD, None))
        if not specs:
            raise ValueError("select at least one field")
        return specs

    def _parse_normalize_value(self):
        raw = self.normalize_value_var.get().strip()
        if not raw or raw.lower() == "max":
            return None
        value = float(raw)
        if value == 0:
            raise ValueError("normalize value must be non-zero")
        return value

    def _invalidate_preview(self):
        self.data = None
        self.density = None
        self.write_btn["state"] = "disabled"
        self.show_field_combo["state"] = "disabled"
        self._set_slice_controls_enabled(False)

    def _on_show_field_changed(self):
        if self.data is None:
            return
        self.density = self.data[self.show_field_var.get()]
        self._redraw_preview()

    def _on_display_option_changed(self):
        self._redraw_preview()

    def _redraw_preview(self):
        if self.density is None:
            return
        title = f"{self.density_params.get('filename', '')}: {self.show_field_var.get()}"
        make_preview_figure(self.density, log_scale=self.log_scale_var.get(),
                             cmap=self.cmap_var.get(), title=title,
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
            mesh_names = resolve_mesh_names(params["field_specs"], params["species"])
            fields = [f for f, _ in params["field_specs"]]
            data = ts.get_fields(fields, level=params["level"], dtype=params["dtype"],
                                 device=device, geometry=geometry)

            total_bytes = sum(arr.nbytes for arr in data.values())
            print(f"Refined data size (after refine, level={params['level']}, {len(data)} field(s)): "
                  f"{format_bytes(total_bytes)}")
            self.file_sizes = (input_size, total_bytes)

            if params["normalize"]:
                norm_value = params["normalize_value"]
                normalize = True if norm_value is None else norm_value
            else:
                normalize = False
            apply_normalization(data, normalize=normalize, threshold=params["threshold"])

            self.root.after(0, self._preview_done, data, mesh_names, params)
        except Exception as exc:
            self.root.after(0, self._on_error, str(exc))

    def _preview_done(self, data, mesh_names, params):
        self.data = data
        self.mesh_names = mesh_names
        self.density_params = params

        fields = list(data)
        shown = self.show_field_var.get()
        if shown not in data:
            shown = fields[0]
        self.show_field_var.set(shown)
        self.show_field_combo["values"] = fields
        self.show_field_combo["state"] = "readonly"
        self.density = data[shown]

        self._set_slice_controls_enabled(True)
        self._configure_slice_controls(self.density)
        self._redraw_preview()

        before_size, after_size = self.file_sizes
        maxima = ", ".join(f"{f}={arr.max():.3e}" for f, arr in data.items())
        self._set_idle(f"Preview ready - shape={self.density.shape}, max: {maxima}  |  "
                        f"file size before refine: {format_bytes(before_size)}, "
                        f"after refine: {format_bytes(after_size)}")
        self.write_btn["state"] = "normal"

    # ------------------------------------------------------------------
    # Write
    # ------------------------------------------------------------------

    def _on_write(self):
        if self.data is None:
            messagebox.showwarning("Warning", "Generate a preview first.")
            return

        params = self.density_params
        output_path = os.path.join(params["run_directory"], f"{params['species']}_{params['output_name']}_0.h5")

        meshes = "\n".join(f"  {f} -> {self.mesh_names[f]}" for f in self.data)
        if os.path.exists(output_path):
            prompt = f"{output_path} already exists and will be overwritten.\n\nMeshes:\n{meshes}\n\nContinue?"
        else:
            prompt = f"Write {output_path}?\n\nMeshes:\n{meshes}"
        if not messagebox.askyesno("Confirm", prompt):
            return

        self._set_busy(f"Writing {output_path}...")
        threading.Thread(target=self._write_worker, args=(params,), daemon=True).start()

    def _write_worker(self, params):
        try:
            ts = Convert(params["run_directory"], params["filename"])
            output_path = ts.write2openpmd(self.data, species=params["species"],
                                            output_name=params["output_name"], author=params["author"],
                                            mesh_names=self.mesh_names)
            self.root.after(0, self._write_done, output_path)
        except Exception as exc:
            self.root.after(0, self._on_error, str(exc))

    def _write_done(self, output_path):
        self._set_idle(f"Wrote {output_path}")
        self.write_btn["state"] = "normal"
        messagebox.showinfo("Success", f"Wrote {output_path}")

    def _on_error(self, message):
        self._set_idle("Error")
        self.write_btn["state"] = "disabled" if self.data is None else "normal"
        messagebox.showerror("Error", message)


def main():
    root = tk.Tk()
    Flash2OpenPMDGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
