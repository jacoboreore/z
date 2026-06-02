# 🧭 Valley Scanner — A Numerical Framework to Explore the Riemann Zeta Landscape

**Author:** Jacob Orellana Real

**Contact:** jacoboreore@gmail.com

**Updated:** 2026/Jun/02


---

## 📘 Overview

**Valley Scanner** is a computational project designed to explore and confirm zeros of the **Riemann Zeta function** through a *homogeneous scanning method* that does **not rely on Gram points, alignment corrections, or interpolation models**.

Instead, the method performs a continuous walk across real numbers converted into the complex domain, identifying *valleys* (local minima) in $|Z(s)|$ that correspond precisely to nontrivial zeros on the critical line $\Re(s)=\frac{1}{2}$.

This project demonstrates a **direct, reproducible pathway to detect zeros**, making the exploration transparent, geometrically intuitive, and free from traditional boundary dependencies.

---

## 🎯 Motivation

Classical approaches to locating zeta zeros often rely on **Gram point alignment** or numerical phase tracking, which introduce discrete artifacts or “gaps” at scale.

The **Valley Scanner** instead follows a continuous, density-preserving process:

* It walks through $t$ values along the critical line, evaluating $|Z(1/2 + it)|$ directly.
* It detects local minima ("valleys") between successive maxima ("mountains").
* Each confirmed valley corresponds to a verified zero of $\zeta(s)$.

This homogeneous structure reproduces the correct zero density naturally, without preconditioning or fitted corrections.
In short: **the terrain itself reveals the zeros.**

---

## 🧩 What’s Published

### 🧠 Research Highlights

* Demonstration that $|Z(s)|$ terrain symmetry can reproduce zero density without Gram alignment.
* Experimental validation at multiple $t$ scales (up to $t \approx 2\times10^{15}$).
* high-t evaluations beyond $t = 10^{20}$ showing numerical stability and precision.
* Complete datasets and Python tools for reproducibility.

Repositories:

> [Datasets and exploratory tools](https://github.com/zjore/z-research/)

> [Core engine](https://github.com/jacoboreore/z)

---

## 🧪 Playground: Visualizing Mountains and Valleys

The file **`playground.py`** offers the simplest way to watch the valley scanner in action.
It plots $|Z(s)|$ against $t$ for low-height ranges, showing:

* **Mountains:** local maxima of $|Z(s)|$
* **Valleys:** confirmed zeros of $\zeta(s)$

To run it:

```bash
python3 playground.py
```

You’ll see a graphical plot where each valley corresponds to a confirmed zero.
This script is ideal for newcomers who want to visually confirm how the method detects zeros geometrically.

---

## 📂 Datasets

The repository includes:

* Datasets published in the paper (`refined_sample_*.csv`)
* Additional high-t datasets confirming reproducibility and valley consistency across independent runs.

Each dataset contains:

```
t, absZ, spacing
```

Where:

* **t** = imaginary part along the critical line
* **absZ** = $|Z(s)|$ value
* **spacing** = $\Delta t$ between consecutive zeros

Example dataset link:

> [refined_sample_1122334455.csv](https://github.com/zjore/z-research/blob/main/datasets/refined_sample_1122334455.csv)

---

## ⚙️ Technical Notes

* Calculations use **MPFR precision arithmetic** for numerical stability.
* High-t runs are distributed across multi-core EC2 instances (up to 192 CPUs per job).
* Each confirmed zero is verified by symmetric evaluation of conjugate terms.
* The process supports horizontal scaling — distributed batches can theoretically reach *cosmic-$t$* ranges (beyond $10^{22}$) with sufficient resources.

---

## 🐳 Docker Image

A prebuilt container is available on Docker Hub for direct use of the core computational engine.

Pull the image:

```bash
docker pull jacobore/jor_z:native
```

Run an interactive shell:

```bash
docker run --rm -it jacobore/jor_z:native
```

Once inside the container, you can execute the core tools through the helper script:

### 🔧 Core Execution Commands

| Command                                                   | Description                                                                      |
| --------------------------------------------------------- | -------------------------------------------------------------------------------- |
| `/app/run.sh z <t> <job_id>`                          | Compute **Z(t)** at a single point.                                              |
| `/app/run.sh z_ball <t> <job_id>`                     | Compute **Z(t)** using the alternative *ball-based* evaluation.                  |
| `/app/run.sh pre_filter <t>`                          | Perform a *pre-filter* evaluation at `t` (quick test).                           |
| `/app/run.sh valley_scanner <t> <job_id>`             | Scan a small region around `t` for candidate zeros.                              |
| `/app/run.sh valley_walk <t> <num_zeros> <job_id>`    | Walk along successive valleys starting from `t` for a specified number of zeros. |
| `/app/run.sh refine <t> <job_id>`                     | Refine a previously estimated zero around `t`.                                   |
| `/app/run.sh refine_progressive <file_path> <job_id>` | Perform a *progressive refinement* pass using an existing CSV dataset.           |
| `/app/run.sh z_os <t> <job_id>`                       | Compute **Z(t)** using the *Odlyzko–Schönhage* algorithm.                        |

---

### ⚙️ CPU and Parallelism Configuration

The container uses **OpenMP** for parallel processing.
You can fine-tune CPU concurrency and thread binding to match your hardware by setting the following environment variables when running the container:

```bash
docker run --rm -it \
  -e OMP_NUM_THREADS=32 \
  -e OMP_PROC_BIND=spread \
  -e OMP_PLACES=cores \
  jacobore/jor_z:native
```

**Variables:**

* `OMP_NUM_THREADS` → number of threads (typically equal to your CPU core count).
* `OMP_PROC_BIND` → thread placement policy (`close`, `spread`, or `master`).
* `OMP_PLACES` → where threads are placed (`cores`, `sockets`, or `threads`).

Tuning these values allows the core to **adapt CPU concurrency** for your machine’s architecture and achieve optimal performance.

---

### 🌐 Web Interface

The same functionality is available through the web interface:

👉 [**https://p56yzukrvv.us-east-1.awsapprunner.com/**](https://p56yzukrvv.us-east-1.awsapprunner.com/)

This interface allows users to launch the same computational modes from a browser-based environment with live progress tracking and dataset export capabilities.

---

🧩 *Note:*
The Docker image contains only the compiled binaries and runtime configuration.
The computation core source code is withheld until peer validation is complete, after which it will be shared for reproducibility.


## 📈 Reproducibility & Future Work

The project aims to:

1. Provide **open datasets** for peer validation.
2. Encourage **collaboration on large-scale distributed scans**.
3. Refine the theoretical interpretation of valley geometry and zero density.

---

## 💡 Quick Start

Clone the repository:

```bash
git clone https://github.com/zjore/z-research.git
cd z-research
```

Install dependencies (if required):

```bash
pip install -r requirements.txt
```

Run a quick valley scan visualization:

```bash
python3 playground.py
```

---

## 🛰️ License and Citation

This project is published under the **MIT License**.
If you reference the datasets or visualizations, please cite:

> Orellana, J. *Valley Scanner: A Continuous Method for Detecting Zeta Zeros Without Gram Alignment.* (2025).

---

### 📜 License and Usage
- Core computation binaries and algorithms: **Proprietary — All rights reserved.**
- Datasets and visualization scripts: **MIT License** (see [LICENSE](https://github.com/zjore/z-research/blob/main/LICENCE.md)).

---

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17566257.svg)](https://doi.org/10.5281/zenodo.17566257)

---

*“Walk the mountains, rest at the valleys. All is revealed with symmetry.”*
— *Valley Scanner Project, 2026*
