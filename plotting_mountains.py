# -----------------------------------------------------------------------------
# Valley Scanner © 2025 Jacob Orellana
# All rights reserved. Redistribution and use in source or binary forms
# are permitted only with explicit written permission of the author.
#
# File: plotting_mountains.py
# Purpose:
#   Plot mountains and valleys and calculated minima
# -----------------------------------------------------------------------------

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import sys

if len(sys.argv) != 2:
    raise ValueError("Specify the plotting CSV file")

try:
    import mplcursors
except ImportError as e:
    raise SystemExit("Please install mplcursors first: pip install mplcursors") from e

csv_path = sys.argv[1]
df = pd.read_csv(csv_path)

# Expected columns
if not {"t", "absZ"}.issubset(df.columns):
    raise ValueError(f"CSV must contain at least columns {{'t','absZ'}}; found {list(df.columns)}")

# Optional columns
if "minima" not in df.columns:
    df["minima"] = 0
has_t_str = "t_str" in df.columns

# Choose which numeric series to plot (t as numbers)
t_num = df["t"].astype(float)

# Smoothing (centered moving average)
n = len(df)
win = max(5, int(round(n * 0.02)))
if win % 2 == 0:
    win += 1
smooth_absZ = df["absZ"].rolling(window=win, center=True, min_periods=max(3, win//3)).mean()

# Build plot
fig, ax = plt.subplots()

# Raw series (capture handle for hover)
line_raw = ax.semilogy(t_num, df["absZ"], marker='o', linestyle='-')[0]

# Smoothed overlay (optional)
if smooth_absZ.notna().sum() > 0:
    ax.semilogy(t_num, smooth_absZ, linestyle='--')

# Vertical lines at minima
m = (df["minima"] == 1)
min_x = t_num[m].to_numpy()
min_y = df.loc[m, "absZ"].to_numpy()
for x in min_x:
    ax.axvline(x=x, linestyle=':')

# Distinct markers for minima (capture handle for hover)
line_min = None
if len(min_x) > 0:
    line_min = ax.semilogy(min_x, min_y, marker='s', linestyle='')[0]

# Axis formatting: scientific x labels + top axis with log10(t)
sf = mticker.ScalarFormatter(useMathText=True)
sf.set_powerlimits((0, 0))             # force scientific notation
ax.xaxis.set_major_formatter(sf)
ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))

def fwd(x):   # t -> log10(t)
    with np.errstate(divide='ignore'):
        return np.log10(x)
def inv(x):   # log10(t) -> t
    return 10.0 ** x

secax = ax.secondary_xaxis('top', functions=(fwd, inv))
secax.set_xlabel('log10(t)')
secax.xaxis.set_major_locator(mticker.MaxNLocator(6))

ax.set_xlabel("t")
ax.set_ylabel("|Z|")
ax.set_title("Mountains & Valleys: |Z| vs t")
ax.grid(True, which='both', axis='y')
plt.tight_layout()

# --- Hover tooltips ---
def fmt_label(idx: int) -> str:
    # Safeguard index bounds
    if idx < 0 or idx >= len(df):
        return ""
    t_label = (str(df.iloc[idx]["t_str"]) if has_t_str and pd.notna(df.iloc[idx]["t_str"])
               else f"{t_num.iloc[idx]:.6e}")
    log10_t = np.log10(float(t_num.iloc[idx])) if t_num.iloc[idx] > 0 else float('nan')
    z = df.iloc[idx]["absZ"]
    is_min = int(df.iloc[idx]["minima"]) == 1
    base = f"t = {t_label}\nlog10(t) = {log10_t:.3f}\n|Z| = {z:.3g}"
    return ("MINIMA\n" + base) if is_min else base

# Cursor for the raw line (all points)
cursor_main = mplcursors.cursor(line_raw, hover=True)
@cursor_main.connect("add")
def _on_add(sel):
    i = int(sel.target.index)
    sel.annotation.set_text(fmt_label(i))
    sel.annotation.get_bbox_patch().set(alpha=0.85)

# Cursor for minima markers (subset). Map its point index back to original df index.
if line_min is not None:
    min_indices = df.index[m].to_numpy()
    cursor_min = mplcursors.cursor(line_min, hover=True)
    @cursor_min.connect("add")
    def _on_add_min(sel):
        j = int(sel.target.index)  # index within minima subset
        i = int(min_indices[j])    # original df row index
        sel.annotation.set_text(fmt_label(i))
        sel.annotation.get_bbox_patch().set(alpha=0.9)

out_path = "output/z_plot_semilogy.png"
plt.savefig(out_path, dpi=150)
