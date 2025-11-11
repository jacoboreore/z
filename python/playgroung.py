import subprocess
import sys

def ensure_lib(package):
    try:
        __import__(package)
    except ImportError:
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", package, "--quiet"]
        )

# Ensure required libraries
ensure_lib("mpmath")
ensure_lib("matplotlib")

import mpmath as mp
import matplotlib.pyplot as plt

mp.dps = 250 # precision

# From N, we Calculate:
# m + sqrt(N - m^2) i
def N_to_complex(N, m):
    t = mp.sqrt(N - m ** 2)
    s = m + t * 1j
    return s, t

# Z
def Z(N, m):
    s, t = N_to_complex(N, m)
    z_val = mp.zeta(s)
    return t, z_val

# Dataset for plotting
dataset = []
m = mp.mpf('0.5') # <-------------- IMPORTANT: MEAN = 1/2
for n in range(5000):
    t, z = Z(n, m)
    dataset.append([t.real, abs(z)])

ts = [t for t, _ in dataset]
zs = [z for _, z in dataset]

plt.figure(figsize=(10, 5))
plt.plot(ts, zs)
plt.xlabel("t (imaginary part of s)")
plt.ylabel("|Z(s)|")
plt.title("Z Test with N -> complex transformation, mean = 1/2")
plt.legend()
plt.grid(True)
plt.show()