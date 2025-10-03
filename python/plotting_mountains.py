# Load the CSV, plot |Z| vs t with minima markers & vertical lines, plus an optional smoothed curve.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Path provided by the system message
csv_path = "../output/valley_walk_plot_20251002_210601_pid33196.csv"

# Read CSV
df = pd.read_csv(csv_path)

# Basic validation of expected columns
expected_cols = {"t", "absZ"}
if not expected_cols.issubset(df.columns):
    raise ValueError(f"CSV must contain columns {expected_cols}. Found: {list(df.columns)}")

# Handle optional 'minima' column (create if missing as zeros)
if "minima" not in df.columns:
    df["minima"] = 0

# Build a smoothed series with a centered moving average
n = len(df)
# Choose a window ~2% of data length, at least 5 and odd
win = max(5, int(round(n * 0.02)))
if win % 2 == 0:
    win += 1
smooth_absZ = df["absZ"].rolling(window=win, center=True, min_periods=max(3, win//3)).mean()

# Prepare figure
plt.figure()

# Raw series
plt.semilogy(df["t"], df["absZ"], marker='o', linestyle='-')

# Smoothed series (if we have enough points for smoothing to matter)
if smooth_absZ.notna().sum() > 0:
    plt.semilogy(df["t"], smooth_absZ, linestyle='--')

# Vertical lines where minima == 1
minima_mask = (df["minima"] == 1)
minima_t_values = df.loc[minima_mask, "t"].values

for x in minima_t_values:
    plt.axvline(x=x, linestyle=':')  # vertical markers

# Optionally, also mark the points with a distinct marker
if minima_t_values.size > 0:
    # Pick the absZ values at those x to place point markers
    y_vals = df.loc[minima_mask, "absZ"].values
    plt.semilogy(minima_t_values, y_vals, marker='s', linestyle='')

plt.xlabel("N")
plt.ylabel("|Z| (log scale)")
plt.title("Mountains & Valleys: |Z| vs N (raw + smoothed)")
# plt.grid(True)
plt.tight_layout()

out_path = "../output/z_plot_semilogy_minima_and_smooth.png"
plt.savefig(out_path, dpi=150)
plt.show()

out_path
