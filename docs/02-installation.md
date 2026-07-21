# Installation

## Requirements

| Component | Version | Notes |
|---|---|---|
| Python | ≥ 3.8 | 3.10+ recommended |
| TNTBLAST | ≥ 2.4 | Must be on your `PATH` |
| pandas | ≥ 1.0 | pandas 2.x compatible |
| numpy | any recent | |
| matplotlib | any recent | Uses Agg backend — no display required |

All Python packages are available via conda-forge or pip and are standard in Anaconda/Miniconda scientific environments.

---

## Step 1 — Install TNTBLAST

TNTBLAST must be compiled and installed separately. Source code and build instructions are at:

```
https://github.com/jgans/thermonucleotideBLAST
```

After building, confirm it is on your `PATH`:

```bash
tntblast --version
```

If this command returns a version number, TNTBLAST is correctly installed. If you get `command not found`, add the TNTBLAST binary directory to your `PATH` in your shell profile (`.bashrc`, `.bash_profile`, or equivalent).

---

## Step 2 — Install Python dependencies

**Using conda (recommended on HPC):**

```bash
conda create -n pcrstrainer python=3.10
conda activate pcrstrainer
conda install -c conda-forge pandas numpy matplotlib
```

**Using pip:**

```bash
pip install pandas numpy matplotlib
```

---

## Step 3 — Clone the repository

```bash
git clone https://github.com/BCCDC-PHL/PCR_strainer.git
cd PCR_strainer
```

To update to the latest version:

```bash
git pull origin main
```

---

## Step 4 — Verify the installation

Run PCR_strainer with no arguments to confirm it loads correctly and prints the usage message:

```bash
python pcr_strainer.py
```

You should see output beginning with `PCR_strainer v0.2.x` followed by usage instructions. If you see an `ImportError`, one of the Python dependencies is missing — install it with `pip install <package>`.

---

## HPC-specific notes

At BCCDC-PHL, PCR_strainer runs on the HPC cluster. The recommended setup:

1. Load the conda environment in your job script:
   ```bash
   conda activate pcrstrainer
   ```

2. TNTBLAST is computationally intensive. For large genome sets (tens of thousands of sequences), request at least 8 GB of memory per TNTBLAST job. PCR_strainer itself is lightweight — the memory ceiling is TNTBLAST.

3. Output HTML reports are transferred to the Windows network share via WinSCP after each run. Open them from a mapped drive letter (e.g., `Z:\results\run1_report.html`), not a UNC path (`\\server\share\...`), to avoid browser zone restrictions on local HTML files.

---

## Verifying the version

The version is stored in `pcr_strainer.py` at module level:

```python
__version__ = '0.2.9'
```

To check which version you are running:

```bash
python -c "import pcr_strainer; print(pcr_strainer.__version__)"
# or simply:
python pcr_strainer.py 2>&1 | head -3
```

Record the version used in any quality documentation. Do not rename the file — versions are tracked via git tags. To see the full change history:

```bash
git log --oneline
git tag
```
