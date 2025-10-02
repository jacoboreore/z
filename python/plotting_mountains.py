# Load CSV, scale big N for plotting, and generate two clear charts.
# - Chart 1: semilogy(absZ) vs N (joined with lines)
# - Chart 2: spike view: X = (N - N0) / median_spacing, Y = -log10(absZ) (joined with lines)
# Saves PNGs for download.

import os
import math
from decimal import Decimal, getcontext
import csv

import numpy as np
import matplotlib.pyplot as plt

# --- Find valleys (local minima) of |Z| vs N -----------------------------
# A point i is a valley if Z[i-1] > Z[i] < Z[i+1].
# Optional: require some prominence so tiny wiggles are ignored.

def find_valleys(N, Z, min_prominence_ratio=1.2, window=1):
    """
    N, Z: numpy arrays (sorted by N).
    min_prominence_ratio: require both neighbors to be at least this multiple of Z[i]
                          e.g., 1.2 means neighbors >= 1.2 * Z[i]
    window: neighbors checked at +/- window indices (use 1 for immediate neighbors)
    """
    n = len(Z)
    valleys = []
    if n < 3:
        return valleys

    for i in range(window, n - window):
        left = Z[i - window]
        mid  = Z[i]
        right = Z[i + window]

        # local min shape
        if left > mid < right:
            # simple prominence: neighbors are sufficiently higher than the valley
            if (left >= min_prominence_ratio * mid) and (right >= min_prominence_ratio * mid):
                valleys.append(i)
    return np.array(valleys, dtype=int)

# tune min_prominence_ratio if you get too many/few lines
valley_idx = find_valleys(N_float, Zabs_float, min_prominence_ratio=1.15, window=1)
valley_N   = N_float[valley_idx]
valley_Z   = Zabs_float[valley_idx]

# --- Helper to annotate to x-axis tick area (small, unobtrusive) --------
def annotate_valley_lines(ax, xs, label_fmt="{:.0f}", y_pct=0.02):
    """
    Draw vertical dashed lines at xs and annotate their x values near the x-axis.
    y_pct: fraction of the y-span above the bottom to place labels.
    """
    # draw lines
    for x in xs:
        ax.axvline(x, linestyle='--', linewidth=1)

    # annotate N near the bottom of the current y-limits
    ymin, ymax = ax.get_ylim()
    y = ymin * (1.0 + 0.0) if ax.get_yscale() == 'log' else (ymin + (ymax - ymin) * y_pct)
    for x in xs:
        ax.annotate(label_fmt.format(x),
                    xy=(x, y),
                    xytext=(0, 3), textcoords="offset points",
                    ha='center', va='bottom', fontsize=8, rotation=90)

# 1) Load the CSV (t,absZ,N)
path = "../output/valley_walk_plot_20251001_135048_pid16764.csv"
rows = []
with open(path, newline="") as f:
    reader = csv.reader(f)
    header = next(reader, None)
    # Detect header
    if header and [h.strip().lower() for h in header] == ["t","absz","n"]:
        pass  # already consumed
    else:
        # No header; rewind and treat the first row as data
        f.seek(0)
        reader = csv.reader(f)
    for r in reader:
        if not r or len(r) < 3:
            continue
        rows.append((r[0].strip(), r[1].strip(), r[2].strip()))

# Basic guard
if not rows:
    raise RuntimeError("CSV appears empty or malformed (need columns: t,absZ,N).")

# 2) Convert to float for plotting + Decimal for precise scaling
getcontext().prec = 80
N_dec = [Decimal(r[2]) for r in rows]  # NOTE: r[2] is N, r[0] is t, r[1] is absZ
Z_dec = [Decimal(r[1]) for r in rows]

# Prepare float arrays (with safe clipping for |Z|)
N_float = np.array([float(n) for n in N_dec], dtype=np.float64)
eps = 1e-300  # prevent log(0)
Zabs_float = np.array([float(max(Decimal(eps), z)) for z in Z_dec], dtype=np.float64)

# (optional) sort by N in case input rows are out of order
order = np.argsort(N_float)
N_float = N_float[order]
Zabs_float = Zabs_float[order]

# 3) Chart 1 — raw semilogy: |Z| vs N (joined with lines)
plt.figure()
plt.semilogy(N_float, Zabs_float, marker='o', linestyle='-')  # <-- joined with lines
plt.xlabel("N")
plt.ylabel("|Z| (log scale)")
plt.title("Mountains & Valleys: |Z| vs N (raw)")
plt.tight_layout()
raw_png = "../output/z_plot_semilogy_raw.png"
plt.savefig(raw_png, dpi=150)
plt.show()

# 4) Chart 2 — centered & scaled X, and spike view on Y (joined with lines)
# if len(N_dec) > 1:
#     diffs = np.abs(np.diff(N_float))
#     median_spacing = np.median(diffs)
#     if median_spacing == 0 or not np.isfinite(median_spacing):
#         median_spacing = 1.0
# else:
#     median_spacing = 1.0

# N0 = float(N_dec[0])
# X_scaled = (N_float - N0) / float(median_spacing)

# # Spike view: Y = -log10(|Z|)
# Y_spike = -np.log10(Zabs_float)

# plt.figure()
# plt.plot(X_scaled, Y_spike, marker='o', linestyle='-')  # <-- joined with lines
# plt.xlabel("(N - N0) / median_spacing")
# plt.ylabel("-log10(|Z|)")
# plt.title("Spikes at zeros: centered/scaled X")
# plt.tight_layout()
# spike_png = "../output/z_plot_spikes_scaled.png"
# plt.savefig(spike_png, dpi=150)
# plt.show()

# (raw_png, spike_png)
