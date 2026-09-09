import numpy as np
import os
import threading
import time
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib import gridspec
import tkinter as tk
from tkinter import ttk, font, filedialog, messagebox

import SSS_functions as sss
import mosaic_to_geotiff as ma

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
SPEED_CONVERSION_FACTOR = 0.51 * 0.10000000149011612    #Knots to m/s
DEFAULT_16BIT_MAX = 32768                               #(2^16) / 2 max raw value
GRID_PADDING_METERS = 100.0                             #Metric padding around UTM grid

#blend mode key used by generate_mosaic_matrix for aerial view
BLEND_MODE_OPTIONS = {
    "Sum": "sum",
    "Average": "average",
    "Overlay (latest on top)": "overlay",
}

# =============================================================================
# ETA & PERCENTAGE TRACKING
# =============================================================================
class ProgressModal(tk.Toplevel):
    """Non-blocking popup dialog providing determinate progress updates and live ETA."""
    def __init__(self, parent, total_steps=100, title_text="Processing...", mode="determinate"):
        super().__init__(parent)
        self.title("Please Wait")
        w, h = 360, 140
        self.resizable(False, False)
        self.transient(parent)
        self.grab_set()

        # Center modal over parent window
        parent.update_idletasks()
        px = parent.winfo_rootx()
        py = parent.winfo_rooty()
        pw = parent.winfo_width()
        ph = parent.winfo_height()
        x = px + max(0, (pw - w) // 2)
        y = py + max(0, (ph - h) // 2)
        self.geometry(f"{w}x{h}+{x}+{y}")

        self.total_steps = max(1, total_steps)
        self.start_time = time.time()
        self.mode = mode

        self.status_lbl = ttk.Label(self, text=title_text, font=("Helvetica", 9, "bold"), anchor="center")
        self.status_lbl.pack(pady=(12, 4), padx=10, fill="x")

        self.progress = ttk.Progressbar(self, mode=self.mode, maximum=100)
        self.progress.pack(fill="x", padx=20, pady=4)

        if self.mode == "determinate":
            self.eta_lbl = ttk.Label(
                self, text="Estimated time remaining: Calculating...",
                font=("Helvetica", 8), foreground="gray", anchor="center"
            )
            self.eta_lbl.pack(pady=(2, 10), padx=10, fill="x")
        else:
            self.progress.start(12)
            self.eta_lbl = None

    def update_progress(self, current_step, status_text=None):
        """Thread-safe update method to update progress bar and ETA label."""
        if self.mode != "determinate":
            if status_text:
                self.after(0, lambda: self.status_lbl.config(text=status_text) if self.winfo_exists() else None)
            return

        def _update():
            if not self.winfo_exists():
                return
            step = min(current_step, self.total_steps)
            pct = (step / self.total_steps) * 100
            self.progress['value'] = pct

            elapsed = time.time() - self.start_time
            if step > 0:
                rate = elapsed / step
                remaining_sec = rate * (self.total_steps - step)
                mins, secs = divmod(int(remaining_sec), 60)
                eta_str = f"{mins:02d}:{secs:02d}" if remaining_sec > 0 else "00:00"
            else:
                eta_str = "Calculating..."

            if status_text:
                self.status_lbl.config(text=status_text)
            else:
                self.status_lbl.config(text=f"Processing... {step}/{self.total_steps} ({pct:.0f}%)")

            if self.eta_lbl:
                self.eta_lbl.config(text=f"Estimated time remaining: {eta_str}")

        self.after(0, _update)

    def close(self):
        def _close():
            if self.winfo_exists():
                if self.mode == "indeterminate":
                    self.progress.stop()
                self.grab_release()
                self.destroy()
        self.after(0, _close)


# =============================================================================
# HORIZONTAL HISTOGRAM + DUAL-HANDLE RANGE SLIDER
# =============================================================================
class HistogramRangeSlider(ttk.Frame):
    """
    A single horizontal histogram of a data distribution with two draggable
    vertical handles marking a min/max display range. Whichever handle is
    closer to a click/drag gets moved; handles cannot cross each other.

    set_data(values, default_min, default_max) rebuilds the histogram and
    (re)initializes the handle positions -- call this whenever the underlying
    data changes (new file, new mode, etc).

    on_change(min_value, max_value) fires (in data units, not %) every time
    the range changes, including once right after set_data().
    """
    def __init__(self, parent, on_change=None, height=90, bins=120,
                 bar_color="#c8791a", handle_color="#2060c0", cmap=None, **kwargs):
        super().__init__(parent, **kwargs)
        self.on_change = on_change
        self.bins = bins
        self.bar_color = bar_color
        self.handle_color = handle_color
        self.cmap = cmap
        self._height = height

        # data_min/data_max are the DISPLAY domain used for pixel mapping.
        # hist_data_min/hist_data_max are the fixed bounds the histogram bars
        # were actually binned over. The display domain starts out equal to
        # the histogram's own bounds, but set_range() can temporarily widen
        # it (e.g. "fit to entire dataset" landing outside this line's range)
        # without touching the bars themselves.
        self.data_min, self.data_max = 0.0, 1.0
        self._hist_data_min, self._hist_data_max = 0.0, 1.0
        self.min_value, self.max_value = 0.0, 1.0
        self._counts = None
        self._bin_edges = None
        self._drag_target = None  # 'min' or 'max' while a drag is in progress

        self.canvas = tk.Canvas(self, height=height, bg="white",
                                 highlightthickness=1, highlightbackground="#999999",
                                 cursor="sb_h_double_arrow")
        self.canvas.pack(fill="x", expand=True)

        self.range_lbl = ttk.Label(self, text="")
        self.range_lbl.pack(anchor="w", pady=(4, 0))

        self.canvas.bind("<Configure>", lambda e: self._redraw())
        self.canvas.bind("<ButtonPress-1>", self._on_press)
        self.canvas.bind("<B1-Motion>", self._on_drag)
        self.canvas.bind("<ButtonRelease-1>", self._on_release)

    # -- public API -----------------------------------------------------
    def set_data(self, data, default_min=None, default_max=None):
        """Rebuild the histogram bars from `data` and reset the handles."""
        data = np.asarray(data, dtype=float)
        data = data[np.isfinite(data)]

        if data.size == 0:
            self._hist_data_min, self._hist_data_max = 0.0, 1.0
            self._counts, self._bin_edges = None, None
        else:
            self._hist_data_min = float(np.min(data))
            self._hist_data_max = float(np.max(data))
            if self._hist_data_max <= self._hist_data_min:
                self._hist_data_max = self._hist_data_min + 1.0
            self._counts, self._bin_edges = np.histogram(
                data, bins=self.bins, range=(self._hist_data_min, self._hist_data_max)
            )

        self.data_min, self.data_max = self._hist_data_min, self._hist_data_max

        lo = float(default_min) if default_min is not None else self.data_min
        hi = float(default_max) if default_max is not None else self.data_max
        self.min_value = min(max(lo, self.data_min), self.data_max)
        self.max_value = min(max(hi, self.data_min), self.data_max)

        self._redraw()
        self._notify()

    def set_range(self, min_value, max_value):
        """
        Reposition the handles to an explicit (min_value, max_value) WITHOUT
        touching the histogram bars -- used by the "fit to best"/"fit to
        data" buttons. If the requested range falls outside the current
        histogram's own bounds, the display domain widens just enough to
        show it (the bars compress toward one side rather than being
        redrawn), and snaps back to the histogram's own bounds next time
        set_range() is called with a narrower request.
        """
        min_value, max_value = float(min_value), float(max_value)
        if min_value > max_value:
            min_value, max_value = max_value, min_value

        self.data_min = min(self._hist_data_min, min_value)
        self.data_max = max(self._hist_data_max, max_value)
        self.min_value = min_value
        self.max_value = max_value

        self._redraw()
        self._notify()

    def get_range(self):
        return self.min_value, self.max_value

    # -- internals --------------------------------------------------------
    def _value_to_px(self, value, width):
        if self.data_max <= self.data_min:
            return 0
        frac = (value - self.data_min) / (self.data_max - self.data_min)
        return int(np.clip(frac, 0.0, 1.0) * width)

    def _px_to_value(self, px, width):
        width = max(width, 1)
        frac = float(np.clip(px / width, 0.0, 1.0))
        return self.data_min + frac * (self.data_max - self.data_min)

    def _redraw(self):
        self.canvas.delete("all")
        width = self.canvas.winfo_width()
        if width <= 1:
            width = 260
        height = self._height

        if self._counts is not None and self._bin_edges is not None and self._counts.max() > 0:
            max_count = self._counts.max()
            bar_colors = self._bar_colors()

            for i, c in enumerate(self._counts):
                if c <= 0:
                    continue
                x0 = self._value_to_px(self._bin_edges[i], width)
                x1 = max(self._value_to_px(self._bin_edges[i + 1], width), x0 + 1)
                bar_h = (c / max_count) * (height - 4)
                y1 = height - 2
                y0 = y1 - bar_h
                self.canvas.create_rectangle(x0, y0, x1, y1, fill=bar_colors[i], outline="")

        min_px = self._value_to_px(self.min_value, width)
        max_px = self._value_to_px(self.max_value, width)

        for px in (min_px, max_px):
            self.canvas.create_line(px, 0, px, height, fill=self.handle_color, width=3)

        self.range_lbl.config(text=f"Min: {self.min_value:.0f}    Max: {self.max_value:.0f}")

    def _bar_colors(self):
        """
        One fill color per histogram bin. With a colormap set, each bin's color
        is the SAME color that value would render as in the actual plot -- i.e.
        the colormap normalized over [min_value, max_value] and clipped at the
        ends, exactly like matplotlib's vmin/vmax. Without a colormap, every bin
        just uses the flat bar_color fallback.
        """
        if self.cmap is None:
            return [self.bar_color] * len(self._counts)

        bin_centers = (self._bin_edges[:-1] + self._bin_edges[1:]) / 2.0
        span = max(self.max_value - self.min_value, 1e-9)
        norm = np.clip((bin_centers - self.min_value) / span, 0.0, 1.0)
        rgba = self.cmap(norm)
        return [f"#{int(r * 255):02x}{int(g * 255):02x}{int(b * 255):02x}" for r, g, b, _ in rgba]

    def _on_press(self, event):
        width = self.canvas.winfo_width()
        min_px = self._value_to_px(self.min_value, width)
        max_px = self._value_to_px(self.max_value, width)
        self._drag_target = 'min' if abs(event.x - min_px) <= abs(event.x - max_px) else 'max'
        self._on_drag(event)

    def _on_drag(self, event):
        if self._drag_target is None:
            return
        width = self.canvas.winfo_width()
        value = self._px_to_value(event.x, width)

        if self._drag_target == 'min':
            self.min_value = min(value, self.max_value)
        else:
            self.max_value = max(value, self.min_value)

        self._redraw()
        self._notify()

    def _on_release(self, event):
        self._drag_target = None

    def _notify(self):
        if self.on_change:
            self.on_change(self.min_value, self.max_value)


# =============================================================================
# DATA CALCULATION & PROCESSING ROUTINES
# =============================================================================
def calculate_utm_grids(filepath, resampling, method_nav, pings_smooth, pings_per_file,
                        Xposition, Yposition, ship_gyro, sensor_speed, sensor_heading,
                        slant_port, port_chan, stbd_chan, progress_callback=None, step_offset=0):
    """
    Computes UTM coordinate grids (utmX, utmY) and sonar data arrays for all files.
    Accumulates smoothed trackline arrays across all files to prevent overwriting.
    """
    max_no_pings = int(np.max(pings_per_file))
    abc = np.arange(resampling - 1, -1, -1)

    utmX = np.zeros((max_no_pings, 2 * resampling, len(filepath)))
    utmY = np.zeros((max_no_pings, 2 * resampling, len(filepath)))
    sonar_data = np.zeros((2 * resampling, max_no_pings, len(filepath)))
    Xsmooth = np.zeros((max_no_pings, len(filepath)))
    Ysmooth = np.zeros((max_no_pings, len(filepath)))

    # Trackline accumulators to prevent overwriting earlier files
    all_Xsmooth_utm = []
    all_Ysmooth_utm = []

    for f_idx in range(len(filepath)):
        file_name = os.path.basename(filepath[f_idx])
        if progress_callback:
            progress_callback(step_offset + f_idx, f"Calculating UTM grids ({f_idx + 1}/{len(filepath)}): {file_name}")

        num_pings = int(pings_per_file[f_idx + 1, 0])
        max_range = np.max(slant_port[:, f_idx]) if slant_port.ndim > 1 else np.max(slant_port)

        if method_nav:
            # Preserved exact 7-parameter signature for sss.smooth_data
            (Xsmooth[:num_pings, f_idx], Ysmooth[:num_pings, f_idx],
             factor_stbd_y, factor_stbd_x, factor_port_y,
             factor_port_x) = sss.smooth_data(Xposition, Yposition, max_no_pings, pings_per_file, pings_smooth, f_idx, ship_gyro)

            # Preserved exact 3-parameter signature for sss.Lat_Long
            Xsmooth_utm, Ysmooth_utm = sss.Lat_Long(Xsmooth[:num_pings, f_idx], Ysmooth[:num_pings, f_idx], False)

            stbd_utmY = np.reshape(Ysmooth_utm, (num_pings, 1)) + np.multiply.outer(factor_stbd_y, abc) * max_range / resampling
            stbd_utmX = np.reshape(Xsmooth_utm, (num_pings, 1)) + np.multiply.outer(factor_stbd_x, abc) * max_range / resampling
            port_utmY = np.reshape(Ysmooth_utm, (num_pings, 1)) + np.multiply.outer(factor_port_y, resampling - 1 - abc) * max_range / resampling
            port_utmX = np.reshape(Xsmooth_utm, (num_pings, 1)) + np.multiply.outer(factor_port_x, resampling - 1 - abc) * max_range / resampling

            utmX[:num_pings, :, f_idx] = np.c_[stbd_utmX, port_utmX]
            utmY[:num_pings, :, f_idx] = np.c_[stbd_utmY, port_utmY]
        else:
            x_init, y_init = sss.Lat_Long(Xposition[num_pings - 1, f_idx], Yposition[num_pings - 1, f_idx], False)
            Xsmooth_utm, Ysmooth_utm = sss.Lat_Long(Xposition[:num_pings, f_idx], Yposition[:num_pings, f_idx], False)
            factor_y, factor_x = sss.Angle(ship_gyro[:num_pings, f_idx])

            distance_y = sensor_speed[:num_pings, f_idx] * SPEED_CONVERSION_FACTOR * factor_y
            distance_x = sensor_speed[:num_pings, f_idx] * SPEED_CONVERSION_FACTOR * factor_x

            lat = y_init + np.flip(np.cumsum(np.flip(distance_y)))
            long = x_init + np.flip(np.cumsum(np.flip(distance_x)))

            heading_index = np.where(sensor_heading[:, f_idx] > 270)
            sensor_heading[heading_index, f_idx] -= 360
            factor_stbd_y, factor_stbd_x = sss.Angle(sensor_heading[:num_pings, f_idx] + 90)
            sensor_heading[heading_index, f_idx] += 360

            heading_index = np.where(sensor_heading[:, f_idx] < 90)
            sensor_heading[heading_index, f_idx] += 360
            factor_port_y, factor_port_x = sss.Angle(sensor_heading[:num_pings, f_idx] - 90)

            stbd_utmY = lat[:, None] + np.multiply.outer(factor_stbd_y, abc) * max_range / resampling
            stbd_utmX = long[:, None] + np.multiply.outer(factor_stbd_x, abc) * max_range / resampling
            port_utmY = lat[:, None] + np.multiply.outer(factor_port_y, resampling - 1 - abc) * max_range / resampling
            port_utmX = long[:, None] + np.multiply.outer(factor_port_x, resampling - 1 - abc) * max_range / resampling

            utmX[:num_pings, :, f_idx] = np.c_[stbd_utmX, port_utmX]
            utmY[:num_pings, :, f_idx] = np.c_[stbd_utmY, port_utmY]

        # Accumulate smoothed UTM tracklines for this file
        all_Xsmooth_utm.append(Xsmooth_utm)
        all_Ysmooth_utm.append(Ysmooth_utm)

        ping_sonar_data = np.r_[
            np.flip(stbd_chan[:, :num_pings, f_idx], axis=0),
            np.flip(port_chan[:, :num_pings, f_idx], axis=0)
        ]
        sonar_data[:, :num_pings, f_idx] = ping_sonar_data

    # Concatenate all files into comprehensive arrays spanning the entire survey extent
    Xsmooth_utm_all = np.concatenate(all_Xsmooth_utm)
    Ysmooth_utm_all = np.concatenate(all_Ysmooth_utm)

    return utmX, utmY, sonar_data, Xsmooth_utm_all, Ysmooth_utm_all


def generate_mosaic_matrix(filepath, utmX, utmY, sonar_data, pings_per_file, resampling, mosaic_resolution,
                           blend_mode="sum", progress_callback=None, step_offset=0):
    """
    Interpolates spatial grid points into a continuous high-resolution raster grid.

    blend_mode controls how pixels covered by more than one line are combined:
      - "sum"     : amplitudes are added together (original behavior).
      - "average" : amplitudes are averaged over the number of lines that actually
                    cover each pixel (empty background pixels are excluded).
      - "overlay" : the most recently processed line wins wherever it has coverage;
                    pixels only covered by earlier lines are left untouched.
    """
    utm_non_zeros_x = np.argwhere(utmX != 0)
    X = np.arange(
        np.min(utmX[utm_non_zeros_x[:, 0], utm_non_zeros_x[:, 1], utm_non_zeros_x[:, 2]]) - GRID_PADDING_METERS,
        np.max(utmX) + GRID_PADDING_METERS,
        mosaic_resolution
    )

    utm_non_zeros_y = np.argwhere(utmY != 0)
    Y = np.arange(
        np.min(utmY[utm_non_zeros_y[:, 0], utm_non_zeros_y[:, 1], utm_non_zeros_y[:, 2]]) - GRID_PADDING_METERS,
        np.max(utmY) + GRID_PADDING_METERS,
        mosaic_resolution
    )

    high_res_grid = np.zeros((len(Y), len(X)))
    coverage_count = np.zeros((len(Y), len(X))) if blend_mode == "average" else None

    data_max = sonar_data.max() if blend_mode == "sum" else None

    for f_idx in range(len(filepath)):
        file_name = os.path.basename(filepath[f_idx])
        if progress_callback:
            progress_callback(step_offset + f_idx, f"Gridding mosaic ({f_idx + 1}/{len(filepath)}): {file_name}")

        num_pings = int(pings_per_file[f_idx + 1, 0])
        file_grid = ma.grid_data(X, Y, utmX, utmY, sonar_data, num_pings, resampling, mosaic_resolution, f_idx, filepath)
        file_mask = file_grid > 0

        if blend_mode == "overlay":
            high_res_grid[file_mask] = file_grid[file_mask]
        elif blend_mode == "average":
            high_res_grid += file_grid
            coverage_count += file_mask
        else:  # "sum" (default, preserves original behavior)
            high_res_grid += file_grid
            high_res_grid.clip(0, data_max, out=high_res_grid)

    if blend_mode == "average":
        high_res_grid = np.divide(
            high_res_grid, coverage_count,
            out=np.zeros_like(high_res_grid),
            where=coverage_count > 0
        )

    return high_res_grid


def compute_gain_filter(port_data, stbd_data, settings, resampling, pings_per_file, progress_callback=None):
    """
    Applies Automatic Gain Control (AGC) or Median Filter algorithms.

    settings layout (see get_current_display_settings):
    [0] display_min, [1] display_max, [2] first_ping, [3] last_ping,
    [4] is_corrected, [5] gain_type, [6] window_size, [7] intensity, [8] file_number
    """
    f_idx = settings[8]
    num_pings = int(pings_per_file[f_idx + 1, 0])

    if settings[5] == "AGC":
        out_data = []
        for i in range(num_pings):
            if progress_callback and i % 10 == 0:
                progress_callback(i, f"Applying AGC Filter ({i}/{num_pings} pings)")
            out_data.append(sss.Automatic_Gain_Control(port_data[:, i, f_idx], stbd_data[:, i, f_idx], settings[6], settings[7]))
        data = np.asarray(out_data).T
        return data[:int(resampling), :], data[int(resampling):, :]

    elif settings[5] == "MedFilt":
        out_data = []
        for i in range(num_pings):
            if progress_callback and i % 10 == 0:
                progress_callback(i, f"Applying Median Filter ({i}/{num_pings} pings)")
            out_data.append(sss.Median_filter(port_data[:, i, f_idx], stbd_data[:, i, f_idx], settings[6]))
        data = np.asarray(out_data).T
        return data[:int(resampling), :], data[int(resampling):, :]

    return None, None


# =============================================================================
# MAIN APPLICATION CLASS (GUI & EVENT MANAGEMENT)
# =============================================================================
class SonarProcessorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Sonar Processing")

        # Application State & Data Variables
        self.filepath = []
        self.file_number = 0
        self.resampling = 1024
        self.pings_per_file = None
        self.port_chan = None
        self.stbd_chan = None
        self.slant_port = None
        self.slant_stbd = None
        self.sonar_altitude = None

        # Grid & Mosaic Data
        self.Xsmooth_utm = None
        self.Ysmooth_utm = None
        self.high_res_grid = None
        self.scale_value = 0
        # Cached gridding inputs, kept around so blend mode can be switched
        # without re-reading the XTF files / recomputing UTM coordinates
        self.utmX = None
        self.utmY = None
        self.sonar_data = None

        # Plotting & Canvas References
        self.ax1 = None
        self.ax2 = None
        self.canvas = None
        self.port_plot_obj = None
        self.stbd_plot_obj = None

        # Bottom Tracking Overlay Lines
        self.blanking_line_port = None
        self.blanking_line_stbd = None
        self.duration_line_port = None
        self.duration_line_stbd = None
        self.threshold_line_port = None
        self.threshold_line_stbd = None

        # Mode State Flags
        self.current_mode = "raw"
        self.threshold_manually_set = False
        self.updating_programmatically = False

        # ---------------------------------------------------------------------
        # Tkinter Control Variables
        # ---------------------------------------------------------------------
        self.resampling_var = tk.IntVar(value=1024)
        self.navigation_var = tk.IntVar(value=1)
        self.course_made_good_var = tk.IntVar(value=300)

        self.blanking_var = tk.IntVar(value=0)
        self.duration_var = tk.IntVar(value=0)
        self.threshold_var = tk.IntVar(value=0)
        self.resolution_var = tk.StringVar(value="0.2")

        # Mosaic overlap blend mode ("sum" / "average" / "overlay")
        self.blend_mode_var = tk.StringVar(value="sum")
        self.blend_mode_display_var = tk.StringVar(value="Sum")

        self.first_ping_var = tk.IntVar(value=0)
        self.last_ping_var = tk.IntVar(value=600)

        self.window_size_var = tk.IntVar(value=151)
        self.intensity_var = tk.IntVar(value=20)

        # UI References
        self.file_listbox = None

        # Frame References
        self.plot_frame = None
        self.sidebar_frame = None
        self.notebook = None
        self.tab_display = None
        self.tab_processing = None
        self.actions_frame = None

        # Display Handles
        self.hist_slider = None  # HistogramRangeSlider instance, created in build_display_tab
        # The currently-displayed line's own [1st percentile, max] amplitude range
        # (recomputed each render_waterfall_plot call) -- used by "Fit to best".
        self.current_file_min = None
        self.current_file_max = None

        self.show_welcome_view()

    # =========================================================================
    # UTILITY METHODS
    # =========================================================================
    def clear_frame(self, frame):
        for widget in frame.winfo_children():
            widget.destroy()

    def upload_files(self):
        files = filedialog.askopenfilenames(
            initialdir="/",
            title="Select XTF File(s)",
            filetypes=(("XTF files", "*.xtf"), ("All files", "*.*"))
        )
        if files:
            self.filepath = list(files)
            self.refresh_file_listbox()

    def refresh_file_listbox(self):
        if self.file_listbox is None:
            return

        self.file_listbox.delete(0, tk.END)
        count = len(self.filepath)
        self.file_listbox.insert(tk.END, f"Files loaded: {count}")

        for f in self.filepath:
            self.file_listbox.insert(tk.END, f"  • {os.path.basename(f)}")

    def get_current_display_settings(self, is_corrected='no', gain_type=None, window_size=0, intensity=0):
        """
        [0] display_min, [1] display_max  -- shared port+starboard amplitude range
        [2] first_ping,  [3] last_ping
        [4] is_corrected,[5] gain_type, [6] window_size, [7] intensity, [8] file_number
        """
        if self.hist_slider is not None:
            display_min, display_max = self.hist_slider.get_range()
        else:
            display_min, display_max = 0, DEFAULT_16BIT_MAX

        return [
            display_min,
            display_max,
            self.first_ping_var.get(),
            self.last_ping_var.get(),
            is_corrected,
            gain_type,
            window_size,
            intensity,
            self.file_number
        ]

    # =========================================================================
    # VIEW 1: WELCOME SCREEN
    # =========================================================================
    def show_welcome_view(self):
        self.clear_frame(self.root)

        ttk.Label(self.root, text="Select XTF File:").grid(column=0, row=1, sticky="w", padx=5, pady=5)
        ttk.Button(self.root, text="Browse Files", command=self.upload_files).grid(column=1, row=1, sticky="w", padx=5, pady=5)

        # Scrollable Box for File Feedback
        box_frame = ttk.Frame(self.root)
        box_frame.grid(column=0, row=2, columnspan=4, sticky="ew", padx=5, pady=5)

        scrollbar = ttk.Scrollbar(box_frame, orient="vertical")
        self.file_listbox = tk.Listbox(
            box_frame,
            yscrollcommand=scrollbar.set,
            height=4,
            width=50,
            font=("Consolas", 9),
            bg="#f0f0f0",
            relief="solid",
            bd=1
        )
        scrollbar.config(command=self.file_listbox.yview)

        scrollbar.pack(side="right", fill="y")
        self.file_listbox.pack(side="left", fill="both", expand=True)

        self.refresh_file_listbox()

        ttk.Label(self.root, text="Resampling:").grid(column=0, row=3, sticky="w", padx=5, pady=5)
        ttk.Radiobutton(self.root, text='1024', variable=self.resampling_var, value=1024).grid(column=1, row=3)
        ttk.Radiobutton(self.root, text='2048', variable=self.resampling_var, value=2048).grid(column=2, row=3)
        ttk.Radiobutton(self.root, text='4096', variable=self.resampling_var, value=4096).grid(column=3, row=3)

        ttk.Label(self.root, text="Navigation:").grid(column=0, row=4, rowspan=2, sticky="w", padx=5, pady=5)
        ttk.Radiobutton(self.root, text='Use fish heading', variable=self.navigation_var, value=0).grid(column=1, row=4, rowspan=2)
        ttk.Radiobutton(self.root, text='Use course made good', variable=self.navigation_var, value=1).grid(column=2, row=4, rowspan=2)
        ttk.Label(self.root, text='(No. pings for smoothing)').grid(column=3, row=4)
        ttk.Spinbox(self.root, from_=10, to=5000, increment=10, textvariable=self.course_made_good_var, width=8).grid(column=3, row=5)

        ttk.Label(self.root, text='Mosaic resolution (m)').grid(column=0, row=6, sticky="w", padx=5, pady=5)
        ttk.Spinbox(self.root, from_=0.05, to=10.0, increment=0.05, textvariable=self.resolution_var, width=8).grid(column=1, row=6)

        ttk.Button(self.root, text="Next", command=self.process_and_show_mosaic).grid(column=0, row=7, columnspan=4, pady=15)

    # =========================================================================
    # VIEW 2: MOSAIC PROCESSING & VIEW
    # =========================================================================
    def process_and_show_mosaic(self):
        # File Selection Guard
        if not self.filepath:
            messagebox.showwarning("No Files Selected", "Please select at least one XTF file before proceeding.")
            return

        num_files = len(self.filepath)
        total_steps = 1 + (num_files * 2)  # 1 for load, num_files for UTM, num_files for Mosaic
        modal = ProgressModal(self.root, total_steps=total_steps, title_text="Processing Sonar Files...")

        def worker():
            try:
                modal.update_progress(0, "Loading XTF file data...")
                self.resampling = self.resampling_var.get()
                method_nav = self.navigation_var.get()
                pings_smooth = self.course_made_good_var.get()
                mosaic_resolution = float(self.resolution_var.get())

                (self.sonar_altitude, Xposition, Yposition, ship_gyro,
                 sensor_speed, sensor_heading, layback, self.slant_port,
                 self.slant_stbd, self.pings_per_file, self.port_chan,
                 self.stbd_chan) = sss.load_data(self.filepath, self.resampling)

                modal.update_progress(1, "Calculating UTM coordinate grids...")
                self.utmX, self.utmY, self.sonar_data, self.Xsmooth_utm, self.Ysmooth_utm = calculate_utm_grids(
                    self.filepath, self.resampling, method_nav, pings_smooth, self.pings_per_file,
                    Xposition, Yposition, ship_gyro, sensor_speed, sensor_heading,
                    self.slant_port, self.port_chan, self.stbd_chan,
                    progress_callback=modal.update_progress, step_offset=1
                )

                modal.update_progress(1 + num_files, "Rasterizing mosaic matrix...")
                self.high_res_grid = generate_mosaic_matrix(
                    self.filepath, self.utmX, self.utmY, self.sonar_data, self.pings_per_file, self.resampling, mosaic_resolution,
                    blend_mode=self.blend_mode_var.get(),
                    progress_callback=modal.update_progress, step_offset=1 + num_files
                )

                modal.update_progress(total_steps, "Complete!")
                self.file_number = 0
                self.scale_value = 0

                def on_complete():
                    modal.close()
                    self.show_mosaic_view()

                self.root.after(0, on_complete)

            except Exception as e:
                def on_error(err=e):
                    modal.close()
                    messagebox.showerror("Processing Error", f"An error occurred while generating the mosaic:\n{err}")
                self.root.after(0, on_error)

        threading.Thread(target=worker, daemon=True).start()

    def show_mosaic_view(self):
        self.clear_frame(self.root)

        blend_frame = ttk.Frame(self.root)
        blend_frame.grid(row=0, column=0, columnspan=4, sticky="w", padx=5, pady=(5, 0))
        ttk.Label(blend_frame, text="Overlap mode:").pack(side="left", padx=(0, 5))
        blend_combo = ttk.Combobox(
            blend_frame, textvariable=self.blend_mode_display_var,
            values=list(BLEND_MODE_OPTIONS.keys()), state="readonly", width=22
        )
        blend_combo.pack(side="left")
        blend_combo.bind("<<ComboboxSelected>>", self.on_blend_mode_change)

        ttk.Button(self.root, text="Waterfall view", command=self.show_waterfall_view).grid(column=0, row=1)
        ttk.Button(self.root, text="Fit color scale to data", command=lambda: self.update_mosaic_scale(1)).grid(column=1, row=1)
        ttk.Button(self.root, text="Fit color scale to best", command=lambda: self.update_mosaic_scale(2)).grid(column=2, row=1)
        ttk.Button(self.root, text="Save mosaic as geotiff", command=self.save_geotiff).grid(column=3, row=1)

        fig_canvas = self.plot_mosaic_canvas(self.scale_value)
        canvas_widget = fig_canvas.get_tk_widget()
        canvas_widget.grid(row=2, column=0, columnspan=4, sticky="nsew")

    def on_blend_mode_change(self, event=None):
        selected_label = self.blend_mode_display_var.get()
        self.blend_mode_var.set(BLEND_MODE_OPTIONS.get(selected_label, "sum"))
        self.regenerate_mosaic()

    def regenerate_mosaic(self):
        """
        Re-grids the mosaic using the currently selected blend mode, reusing the
        already-computed UTM grids/sonar data so it doesn't need to re-read the
        XTF files or recompute navigation.
        """
        if self.utmX is None or self.sonar_data is None:
            return

        blend_mode = self.blend_mode_var.get()
        mosaic_resolution = float(self.resolution_var.get())

        modal = ProgressModal(self.root, title_text="Rebuilding Mosaic...", mode="indeterminate")

        def worker():
            try:
                new_grid = generate_mosaic_matrix(
                    self.filepath, self.utmX, self.utmY, self.sonar_data,
                    self.pings_per_file, self.resampling, mosaic_resolution,
                    blend_mode=blend_mode
                )

                def on_complete():
                    self.high_res_grid = new_grid
                    modal.close()
                    self.show_mosaic_view()
                self.root.after(0, on_complete)

            except Exception as e:
                def on_error(err=e):
                    modal.close()
                    messagebox.showerror("Mosaic Error", f"Failed to rebuild mosaic:\n{err}")
                self.root.after(0, on_error)

        threading.Thread(target=worker, daemon=True).start()

    def update_mosaic_scale(self, scale_val):
        self.scale_value = scale_val
        self.show_mosaic_view()

    def plot_mosaic_canvas(self, histogram_mode):
        my_colors = sss.Color_palette()
        figure = Figure(figsize=(10, 6), dpi=100)
        ax = figure.subplots(1, 1)

        valid_mask_port = (self.port_chan > 0) & (self.port_chan != -999) & (~np.isnan(self.port_chan))
        valid_mask_stbd = (self.stbd_chan > 0) & (self.stbd_chan != -999) & (~np.isnan(self.stbd_chan))

        if histogram_mode == 0:
            minimum, maximum = 0, DEFAULT_16BIT_MAX
        elif histogram_mode == 1:
            minimum = np.percentile((self.port_chan[valid_mask_port], self.stbd_chan[valid_mask_stbd]), 0)
            maximum = np.percentile((self.port_chan[valid_mask_port], self.stbd_chan[valid_mask_stbd]), 100)
        elif histogram_mode == 2:
            minimum = np.percentile((self.port_chan[valid_mask_port], self.stbd_chan[valid_mask_stbd]), 1)
            maximum = np.percentile((self.port_chan[valid_mask_port], self.stbd_chan[valid_mask_stbd]), 100)

        color_values = np.arange(minimum, maximum)
        ax.imshow(self.high_res_grid, cmap=my_colors, aspect='auto', vmin=np.percentile(color_values, 0), vmax=np.percentile(color_values, 100))

        ax.set_xlabel('X distance [m]')
        ticks_x = np.arange(0, self.high_res_grid.shape[1], int(self.high_res_grid.shape[1] / 5))
        ticks_labels_x = np.round(np.arange(min(self.Xsmooth_utm), max(self.Xsmooth_utm), abs((max(self.Xsmooth_utm) - min(self.Xsmooth_utm)) / 5)), 2)
        if ticks_labels_x.shape[0] < ticks_x.shape[0]:
            ticks_labels_x = np.append(ticks_labels_x, np.round(max(self.Xsmooth_utm), 2))
        ax.set_xticks(ticks_x)
        ax.set_xticklabels(ticks_labels_x)

        ax.set_ylabel('Y distance [m]')
        ticks_y = np.arange(0, self.high_res_grid.shape[0], int(self.high_res_grid.shape[0] / 5))
        ticks_labels_y = np.round(np.arange(min(self.Ysmooth_utm), max(self.Ysmooth_utm), abs((max(self.Ysmooth_utm) - min(self.Ysmooth_utm)) / 5)), 2)
        if ticks_labels_y.shape[0] < ticks_y.shape[0]:
            ticks_labels_y = np.append(ticks_labels_y, np.round(max(self.Ysmooth_utm), 2))
        ax.set_yticks(ticks_y)
        ax.set_yticklabels(ticks_labels_y)

        return FigureCanvasTkAgg(figure, self.root)

    def save_geotiff(self):
        valid_mask_port = (self.port_chan > 0) & (self.port_chan != -999) & (~np.isnan(self.port_chan))
        valid_mask_stbd = (self.stbd_chan > 0) & (self.stbd_chan != -999) & (~np.isnan(self.stbd_chan))

        if histogram_mode == 0:
            minimum, maximum = 0, DEFAULT_16BIT_MAX
        elif histogram_mode == 1:
            minimum = np.percentile((self.port_chan[valid_mask_port], self.stbd_chan[valid_mask_stbd]), 0)
            maximum = np.percentile((self.port_chan[valid_mask_port], self.stbd_chan[valid_mask_stbd]), 100)
        elif histogram_mode == 2:
            minimum = np.percentile((self.port_chan[valid_mask_port], self.stbd_chan[valid_mask_stbd]), 1)
            maximum = np.percentile((self.port_chan[valid_mask_port], self.stbd_chan[valid_mask_stbd]), 100)
        

        mosaic_resolution = float(self.resolution_var.get())
        ma.save_mosaic(self.high_res_grid, self.Xsmooth_utm, self.Ysmooth_utm, mosaic_resolution, self.filepath, minimum, maximum)
        messagebox.showinfo("Export Success", f"Mosaic saved as GeoTIFF in directory:\n{os.getcwd()}")

    # =========================================================================
    # VIEW 3: WATERFALL VIEW & MODE CONTROLS
    # =========================================================================
    def show_waterfall_view(self):
        self.clear_frame(self.root)

        self.root.grid_columnconfigure(0, weight=1)
        self.root.grid_columnconfigure(1, weight=0)
        self.root.grid_rowconfigure(0, weight=1)

        self.plot_frame = ttk.Frame(self.root)
        self.plot_frame.grid(row=0, column=0, sticky="nsew", padx=5, pady=5)

        self.sidebar_frame = ttk.Frame(self.root, padding=5)
        self.sidebar_frame.grid(row=0, column=1, sticky="nsew", padx=5, pady=5)
        self.sidebar_frame.grid_rowconfigure(0, weight=1)
        self.sidebar_frame.grid_rowconfigure(1, weight=0)

        self.notebook = ttk.Notebook(self.sidebar_frame)
        self.notebook.grid(row=0, column=0, sticky="nsew")

        self.tab_display = ttk.Frame(self.notebook, padding=10)
        self.tab_processing = ttk.Frame(self.notebook, padding=10)

        self.notebook.add(self.tab_display, text="Display")
        self.notebook.add(self.tab_processing, text="Processing")

        self.build_display_tab(self.tab_display)

        self.actions_frame = ttk.LabelFrame(self.sidebar_frame, text="Actions", padding=8)
        self.actions_frame.grid(row=1, column=0, sticky="ew", pady=(10, 0))

        self.load_raw_mode()

    def build_display_tab(self, parent):
        label_font = font.Font(underline=True)
        ttk.Label(parent, text="Display Settings", font=label_font).pack(anchor="w", pady=(0, 5))

        hist_frame = ttk.LabelFrame(parent, text="Display Range", padding=8)
        hist_frame.pack(fill="x", pady=(0, 10))
        ttk.Label(hist_frame, text="Drag either handle to set the displayed min/max.").pack(anchor="w", pady=(0, 6))

        self.hist_slider = HistogramRangeSlider(hist_frame, on_change=self.update_display_range, cmap=sss.Color_palette())
        self.hist_slider.pack(fill="x")

        fit_btn_frame = ttk.Frame(hist_frame)
        fit_btn_frame.pack(fill="x", pady=(6, 0))
        ttk.Button(fit_btn_frame, text="Fit color scale to data", command=self.fit_waterfall_scale_to_data).pack(side="left", expand=True, fill="x", padx=(0, 3))
        ttk.Button(fit_btn_frame, text="Fit color scale to best", command=self.fit_waterfall_scale_to_best).pack(side="left", expand=True, fill="x", padx=(3, 0))

        ping_frame = ttk.LabelFrame(parent, text="Ping Selection", padding=8)
        ping_frame.pack(fill="x")

        ping_grid = ttk.Frame(ping_frame)
        ping_grid.pack(fill="x")

        ttk.Label(ping_grid, text="First ping:").grid(row=0, column=0, sticky="w", pady=2)
        first_spin = ttk.Spinbox(
            ping_grid, from_=0, to=100000, increment=10,
            textvariable=self.first_ping_var, width=8, command=self.update_ping_range
        )
        first_spin.grid(row=0, column=1, sticky="e")
        first_spin.bind("<KeyRelease>", self.update_ping_range)

        ttk.Label(ping_grid, text="Last ping:").grid(row=1, column=0, sticky="w", pady=2)
        last_spin = ttk.Spinbox(
            ping_grid, from_=0, to=100000, increment=10,
            textvariable=self.last_ping_var, width=8, command=self.update_ping_range
        )
        last_spin.grid(row=1, column=1, sticky="e")
        last_spin.bind("<KeyRelease>", self.update_ping_range)

    def update_display_range(self, min_val, max_val):
        if min_val >= max_val:
            return
        if self.port_plot_obj is not None:
            self.port_plot_obj.set_clim(vmin=min_val, vmax=max_val)
        if self.stbd_plot_obj is not None:
            self.stbd_plot_obj.set_clim(vmin=min_val, vmax=max_val)
        if self.canvas is not None:
            self.canvas.draw_idle()

    def fit_waterfall_scale_to_best(self):
        """Reset the handles to the CURRENTLY DISPLAYED line's own [1st percentile, max]
        range -- same range the histogram defaults to when a file is first loaded."""
        if self.hist_slider is None or self.current_file_min is None:
            return
        self.hist_slider.set_range(self.current_file_min, self.current_file_max)

    def fit_waterfall_scale_to_data(self):
        """Reset the handles to the [min, max] of every loaded line's raw data, not
        just the one currently displayed -- mirrors the mosaic view's own
        'Fit color scale to data' button."""
        if self.hist_slider is None or self.port_chan is None:
            return

        valid_mask_port = (self.port_chan > 0) & (self.port_chan != -999) & (~np.isnan(self.port_chan))
        valid_mask_stbd = (self.stbd_chan > 0) & (self.stbd_chan != -999) & (~np.isnan(self.stbd_chan))
        valid_vals = np.concatenate([self.port_chan[valid_mask_port], self.stbd_chan[valid_mask_stbd]])

        if len(valid_vals) == 0:
            return

        dataset_min = np.percentile(valid_vals, 0)
        dataset_max = np.percentile(valid_vals, 100)
        self.hist_slider.set_range(dataset_min, dataset_max)

    def update_ping_range(self, *args):
        if self.ax1 is None or self.ax2 is None or self.canvas is None:
            return

        try:
            fp = self.first_ping_var.get()
            lp = self.last_ping_var.get()
        except tk.TclError:
            return

        if fp >= lp or fp < 0:
            return

        self.ax1.set_ylim(lp, fp)
        self.ax2.set_ylim(lp, fp)

        ticks_y = np.linspace(fp, lp, 6)
        ticks_labels_y = [int(t) for t in ticks_y]

        self.ax1.set_yticks(ticks_y)
        self.ax1.set_yticklabels(ticks_labels_y)
        self.ax2.set_yticks(ticks_y)
        self.ax2.set_yticklabels(ticks_labels_y)

        if self.current_mode == "raw":
            self.update_bottom_tracking_overlays()

        self.canvas.draw_idle()

    def render_waterfall_plot(self, port_data, stbd_data, settings, computed_port_gain=None, computed_stbd_gain=None):
        self.clear_frame(self.plot_frame)

        my_colors = sss.Color_palette()

        if settings[5] in ("AGC", "MedFilt") and computed_port_gain is not None:
            port_plot, stbd_plot = computed_port_gain, computed_stbd_gain
        else:
            port_plot = port_data[:, :, settings[8]] if port_data.ndim == 3 else port_data
            stbd_plot = stbd_data[:, :, settings[8]] if stbd_data.ndim == 3 else stbd_data

        # Mask out 0, -999, and NaNs for port and starboard channels separately
        valid_p = port_plot[(port_plot > 0) & (port_plot != -999) & (~np.isnan(port_plot))]
        valid_s = stbd_plot[(stbd_plot > 0) & (stbd_plot != -999) & (~np.isnan(stbd_plot))]

        # Concatenate into a SINGLE 1D array -- this is the "whole line, both
        # channels combined" dataset that feeds the shared display histogram.
        if len(valid_p) > 0 or len(valid_s) > 0:
            valid_vals = np.concatenate([valid_p, valid_s])
            minimum = np.percentile(valid_vals, 1)
            maximum = np.percentile(valid_vals, 100)
        else:
            valid_vals = np.array([])
            minimum, maximum = 0, DEFAULT_16BIT_MAX

        file_name = os.path.splitext(os.path.basename(self.filepath[self.file_number]))[0]

        figure = Figure(figsize=(10, 6), dpi=100)
        gs = gridspec.GridSpec(1, 2, wspace=0.0)

        self.ax1 = figure.add_subplot(gs[0, 0])
        self.ax2 = figure.add_subplot(gs[0, 1])

        is_raw = (settings[4] == 'no')

        # vmin/vmax are set immediately below via hist_slider.set_data(), which
        # fires update_display_range() once it recomputes the shared histogram
        # -- no need to compute a clim here.
        self.port_plot_obj = self.ax1.imshow(port_plot.T, cmap=my_colors, aspect='auto')
        self.ax1.set_autoscalex_on(False)

        self.stbd_plot_obj = self.ax2.imshow(stbd_plot.T, cmap=my_colors, aspect='auto')
        self.ax2.set_autoscalex_on(False)

        if is_raw:
            self.duration_line_port = self.ax1.axvline(x=self.resampling, color='green', linestyle='-', linewidth=1)
            self.blanking_line_port = self.ax1.axvline(x=self.resampling, color='red', linestyle='-', linewidth=1)
            self.threshold_line_port, = self.ax1.plot([], [], color='blue', linestyle='-', linewidth=1)

            self.duration_line_stbd = self.ax2.axvline(x=0, color='green', linestyle='-', linewidth=1)
            self.blanking_line_stbd = self.ax2.axvline(x=0, color='red', linestyle='-', linewidth=1)
            self.threshold_line_stbd, = self.ax2.plot([], [], color='blue', linestyle='-', linewidth=1)

        figure.colorbar(self.port_plot_obj, ax=[self.ax1], location='left', pad=0.2)
        ticks_x = np.arange(0, self.resampling + 1, int(self.resampling / 4))
        max_hor_range = sss.Get_max_range(self.filepath, self.file_number)
        ticks_labels_x = np.arange(0, max_hor_range + 1, max_hor_range / 4)

        self.ax1.set_xticks(ticks_x)
        self.ax1.set_xticklabels(np.flip(ticks_labels_x))
        self.ax1.spines['right'].set_visible(False)
        self.ax1.set_ylabel('Along-track distance [ping number]')

        fp, lp = settings[2], settings[3]
        self.ax1.set_ylim(lp, fp)
        self.ax2.set_ylim(lp, fp)

        ticks_y = np.linspace(fp, lp, 6)
        self.ax1.set_yticks(ticks_y)
        self.ax1.set_yticklabels([int(t) for t in ticks_y])

        self.ax2.yaxis.tick_right()
        figure.colorbar(self.stbd_plot_obj, ax=[self.ax2], location='right', pad=0.2)
        self.ax2.set_xticks(ticks_x)
        self.ax2.set_xticklabels(ticks_labels_x)
        self.ax2.set_yticks(ticks_y)
        self.ax2.set_yticklabels([int(t) for t in ticks_y])
        self.ax2.spines['left'].set_visible(False)

        if is_raw:
            figure.suptitle(f'Raw data \n File name: {file_name}')
            self.ax1.set_xlabel('Port-side slant distance [m]')
            self.ax2.set_xlabel('Starboard-side slant distance [m]')
        else:
            prefix = f"{settings[5]} corrected data" if settings[5] else "Corrected data"
            figure.suptitle(f'{prefix} \n File name: {file_name}')
            self.ax1.set_xlabel('Port-side horizontal distance [m]')
            self.ax2.set_xlabel('Starboard-side horizontal distance [m]')

        self.canvas = FigureCanvasTkAgg(figure, self.plot_frame)
        canvas_widget = self.canvas.get_tk_widget()
        canvas_widget.pack(fill="both", expand=True)

        self.canvas.draw()

        # Cache this line's own [1st percentile, max] range for "Fit to best"
        self.current_file_min, self.current_file_max = minimum, maximum

        # Rebuild the shared port+starboard histogram for this dataset and
        # reset its handles to [1st percentile, max] -- this also fires
        # update_display_range() once, which applies the matching clim to
        # the plot objects created just above.
        if self.hist_slider is not None:
            self.hist_slider.set_data(valid_vals, default_min=minimum, default_max=maximum)

    # -------------------------------------------------------------------------
    # MODE 1: RAW DATA / BOTTOM TRACKING
    # -------------------------------------------------------------------------
    def load_raw_mode(self):
        self.current_mode = "raw"
        self.threshold_manually_set = False

        display_settings = self.get_current_display_settings('no')
        self.render_waterfall_plot(self.port_chan, self.stbd_chan, display_settings)

        self.clear_frame(self.tab_processing)
        bt_frame = ttk.LabelFrame(self.tab_processing, text="Bottom Tracking Controls", padding=8)
        bt_frame.pack(fill="x", pady=5)

        self.blank_lbl = ttk.Label(bt_frame, text=f"Blanking: {self.blanking_var.get()}")
        self.blank_lbl.pack(anchor="w")
        self.b_slider = ttk.Scale(
            bt_frame, from_=0, to=np.max(self.slant_port),
            variable=self.blanking_var, command=self.on_blanking_change
        )
        self.b_slider.pack(fill="x", pady=(0, 8))

        self.dur_lbl = ttk.Label(bt_frame, text=f"Duration: {self.duration_var.get()}")
        self.dur_lbl.pack(anchor="w")
        self.d_slider = ttk.Scale(
            bt_frame, from_=0, to=np.max(self.slant_port),
            variable=self.duration_var, command=self.on_duration_change
        )
        self.d_slider.pack(fill="x", pady=(0, 8))

        self.thres_lbl = ttk.Label(bt_frame, text=f"Threshold: {self.threshold_var.get()}")
        self.thres_lbl.pack(anchor="w")
        valid_p = self.port_chan[(self.port_chan > 0) & (self.port_chan != -999) & (~np.isnan(self.port_chan))]
        valid_s = self.stbd_chan[(self.stbd_chan > 0) & (self.stbd_chan != -999) & (~np.isnan(self.stbd_chan))]

        avg = np.mean(np.concatenate([valid_p, valid_s]))
        print(avg)
        self.t_slider = ttk.Scale(
            bt_frame, from_=0, to=avg,
            variable=self.threshold_var, command=self.on_threshold_change
        )
        self.t_slider.pack(fill="x")

        self.clear_frame(self.actions_frame)
        ttk.Button(self.actions_frame, text="← Back to Mosaic View", command=self.show_mosaic_view).pack(fill="x", pady=2)
        ttk.Button(self.actions_frame, text="Apply Slant Correction", command=self.apply_slant_correction).pack(fill="x", pady=2)

        # Graceful End-of-File Guard
        next_btn_state = "normal" if self.file_number + 1 < len(self.filepath) else "disabled"
        ttk.Button(self.actions_frame, text="Plot Next File", command=self.display_next_file, state=next_btn_state).pack(fill="x", pady=2)

        self.update_bottom_tracking_overlays()

    def on_blanking_change(self, val):
        if not self.threshold_manually_set:
            self.updating_programmatically = True
            blank_val = int(self.b_slider.get())
            self.t_slider.set(blank_val)
            self.threshold_var.set(blank_val)
            self.updating_programmatically = False

        self.update_bottom_tracking_overlays()

    def on_duration_change(self, val):
        self.update_bottom_tracking_overlays()

    def on_threshold_change(self, val):
        if not self.updating_programmatically:
            self.threshold_manually_set = True
        self.update_bottom_tracking_overlays()

    def update_bottom_tracking_overlays(self):
        if self.current_mode != "raw":
            return

        blank_val = int(self.b_slider.get())
        dur_val = int(self.d_slider.get())
        thres_val = int(self.t_slider.get())

        self.blank_lbl.config(text=f"Blanking: {blank_val}")
        self.dur_lbl.config(text=f"Duration: {dur_val}")
        self.thres_lbl.config(text=f"Threshold: {thres_val}")

        max_slant = np.max(self.slant_port)
        blank_px = int(blank_val * self.resampling / max_slant)
        tot_px = int((blank_val + dur_val) * self.resampling / max_slant)

        self.blanking_line_port.set_xdata([int(self.resampling - blank_px), int(self.resampling - blank_px)])
        self.blanking_line_stbd.set_xdata([blank_px, blank_px])

        self.duration_line_port.set_xdata([int(self.resampling - tot_px), int(self.resampling - tot_px)])
        self.duration_line_stbd.set_xdata([tot_px, tot_px])

        fp, lp = self.first_ping_var.get(), self.last_ping_var.get()
        if fp < lp and fp >= 0:
            p_peak, s_peak = sss.peak_displayer(
                self.port_chan, self.stbd_chan, blank_val, dur_val, thres_val,
                fp, lp, self.resampling, max_slant, self.file_number
            )
            pings_y = np.arange(fp, lp)
            self.threshold_line_stbd.set_data(s_peak, pings_y)
            self.threshold_line_port.set_data(p_peak, pings_y)

        self.canvas.draw_idle()

    def apply_slant_correction(self):
        modal = ProgressModal(self.root, title_text="Applying Slant Correction...", mode="indeterminate")

        def worker():
            try:
                chosen_blanking = int(self.b_slider.get())
                chosen_duration = int(self.d_slider.get())
                chosen_threshold = int(self.t_slider.get())

                new_port, new_stbd = sss.manual_peak_finder(
                    self.port_chan, self.stbd_chan, chosen_blanking, chosen_duration,
                    chosen_threshold, self.resampling, np.max(self.slant_port),
                    self.file_number, self.pings_per_file, self.sonar_altitude
                )

                def on_complete():
                    modal.close()
                    self.load_corrected_mode(new_port, new_stbd)

                self.root.after(0, on_complete)

            except Exception as e:
                def on_error(err=e):
                    modal.close()
                    messagebox.showerror("Slant Correction Error", f"Failed to compute slant correction:\n{err}")
                self.root.after(0, on_error)

        threading.Thread(target=worker, daemon=True).start()

    # -------------------------------------------------------------------------
    # MODE 2: SLANT CORRECTED & GAIN DATA
    # -------------------------------------------------------------------------
    def load_corrected_mode(self, port_data, stbd_data, gain_type=None):
        self.current_mode = "corrected"

        display_settings = self.get_current_display_settings(
            is_corrected='yes', gain_type=gain_type,
            window_size=self.window_size_var.get(), intensity=self.intensity_var.get()
        )

        if gain_type in ("AGC", "MedFilt"):
            num_pings = int(self.pings_per_file[self.file_number + 1, 0])
            modal = ProgressModal(self.root, total_steps=num_pings, title_text=f"Calculating {gain_type} Filter...")

            def worker():
                try:
                    p_gain, s_gain = compute_gain_filter(
                        port_data, stbd_data, display_settings, self.resampling, self.pings_per_file,
                        progress_callback=modal.update_progress
                    )

                    def on_complete():
                        modal.close()
                        self.render_waterfall_plot(port_data, stbd_data, display_settings, computed_port_gain=p_gain, computed_stbd_gain=s_gain)

                    self.root.after(0, on_complete)

                except Exception as e:
                    def on_error(err=e):
                        modal.close()
                        messagebox.showerror("Gain Calculation Error", f"Failed to apply {gain_type}:\n{err}")
                    self.root.after(0, on_error)

            threading.Thread(target=worker, daemon=True).start()
        else:
            self.render_waterfall_plot(port_data, stbd_data, display_settings)

        self.clear_frame(self.tab_processing)
        gain_frame = ttk.LabelFrame(self.tab_processing, text="Gain Settings", padding=8)
        gain_frame.pack(fill="x", pady=5)

        ttk.Label(gain_frame, text="Window size").pack(anchor="w")
        ttk.Spinbox(gain_frame, from_=3, to=999, increment=2, textvariable=self.window_size_var, width=10).pack(fill="x", pady=(0, 8))

        ttk.Label(gain_frame, text="Intensity (AGC only)").pack(anchor="w")
        ttk.Spinbox(gain_frame, from_=1, to=100, increment=1, textvariable=self.intensity_var, width=10).pack(fill="x", pady=(0, 8))

        ttk.Button(gain_frame, text="Apply AGC Gain", command=lambda: self.load_corrected_mode(port_data, stbd_data, gain_type="AGC")).pack(fill="x", pady=2)
        ttk.Button(gain_frame, text="Apply MedFilt Gain", command=lambda: self.load_corrected_mode(port_data, stbd_data, gain_type="MedFilt")).pack(fill="x", pady=2)

        if gain_type is not None:
            ttk.Button(gain_frame, text="Remove Gain", command=lambda: self.load_corrected_mode(port_data, stbd_data, gain_type=None)).pack(fill="x", pady=2)

        self.clear_frame(self.actions_frame)
        ttk.Button(self.actions_frame, text="← Back to Mosaic View", command=self.show_mosaic_view).pack(fill="x", pady=2)
        ttk.Button(self.actions_frame, text="← Back to Raw Data", command=self.load_raw_mode).pack(fill="x", pady=2)

        # Graceful End-of-File Guard
        next_btn_state = "normal" if self.file_number + 1 < len(self.filepath) else "disabled"
        ttk.Button(self.actions_frame, text="Plot Next File", command=self.display_next_file, state=next_btn_state).pack(fill="x", pady=2)

    def display_next_file(self):
        if self.file_number + 1 < len(self.filepath):
            self.file_number += 1
            self.load_raw_mode()
        else:
            messagebox.showinfo("End of Files", "You have reached the end of the loaded files list.")


# =============================================================================
# APPLICATION ENTRY POINT
# =============================================================================
if __name__ == "__main__":
    root = tk.Tk()
    app = SonarProcessorApp(root)
    root.mainloop()
