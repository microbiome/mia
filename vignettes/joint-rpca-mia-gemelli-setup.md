# Joint-RPCA comparison: `mia` and Gemelli

This project compares Joint-RPCA results from:

- `mia` in R
- Gemelli in Python

Both implementations use the same MGX and MTX data, the same train-test split, and the same analysis settings.

## Project structure

```text
joint-rpca-mia-gemelli/
│
├── README.md
├── joint-rpca-mia-gemelli.Rproj
├── joint-rpca-mia-gemelli-comparison.qmd
├── renv.lock
├── environment.yml
│
├── data/
│   ├── mgx_raw.csv
│   ├── mtx_raw.csv
│   └── train_test_split.csv
│
├── python/
│   ├── run_gemelli.py
│   └── export_python_versions.py
│
└── results/
    ├── mia/
    └── gemelli/
```

## Required software

Install:

1. R
2. RStudio
3. Quarto
4. Windows Subsystem for Linux (WSL)
5. Miniforge or Miniconda inside WSL

## 1. Set up R

Open:

```text
joint-rpca-mia-gemelli.Rproj
```

In the RStudio Console, run:

```r
install.packages("renv")
renv::restore()
```

Install the development version of `mia`:

```r
renv::install("microbiome/mia@devel")
```

Check that `mia` loads:

```r
library(mia)
packageVersion("mia")
```

## 2. Prepare the common input data

Run the data-preparation section of the Quarto notebook.

This creates:

```text
data/mgx_raw.csv
data/mtx_raw.csv
data/train_test_split.csv
```

These files are used by both the R and Python analyses.

## 3. Run Joint-RPCA in R

Run the `mia` Joint-RPCA section in the Quarto notebook.

The main command is:

```r
mia_result <- getJointRPCA(
    mae_rclr,
    experiments = c("MGX", "MTX"),
    assay.types = c("rclr", "rclr"),
    test.set = test_ids,
    ncomponents = 3L,
    max.iterations = 10L
)
```

The R results are saved in:

```text
results/mia/
```

## 4. Set up Python and Gemelli

Open Windows Command Prompt or PowerShell and start WSL:

```bash
wsl
```

Move to the project folder:

```bash
cd /mnt/d/joint-rpca-mia-gemelli
```

Start Conda:

```bash
source ~/miniforge3/etc/profile.d/conda.sh
```

Create the environment:

```bash
conda env create -f environment.yml
```

Activate it:

```bash
conda activate mia-gemelli
```

Check Python and Gemelli:

```bash
python --version
gemelli --help
```

## 5. Run the Python analysis

Inside WSL, run:

```bash
python python/run_gemelli.py
python python/export_python_versions.py
```

The Python results are saved in:

```text
results/gemelli/
```

## 6. Create the final report

Return to Windows Command Prompt.

Move to the project folder:

```cmd
cd /d D:\joint-rpca-mia-gemelli
```

Render the report:

```cmd
quarto render joint-rpca-mia-gemelli-comparison.qmd
```

Expected output:

```text
joint-rpca-mia-gemelli-comparison.html
```

## Recommended running order

1. Open the R project.
2. Run `renv::restore()`.
3. Install the development version of `mia`.
4. Prepare the common MGX and MTX data.
5. Run the `mia` Joint-RPCA analysis.
6. Open WSL.
7. Activate the `mia-gemelli` environment.
8. Run the Gemelli Python script.
9. Export Python package versions.
10. Render the Quarto report.

## Important

Run R code in RStudio:

```r
library(mia)
```

Run terminal commands in Windows Command Prompt, PowerShell, or WSL:

```cmd
quarto render joint-rpca-mia-gemelli-comparison.qmd
```

Windows path example:

```text
D:\joint-rpca-mia-gemelli
```

WSL path example:

```text
/mnt/d/joint-rpca-mia-gemelli
```

## Expected outputs

A successful run creates:

- R Joint-RPCA results
- Gemelli Joint-RPCA results
- comparison tables
- publication-style plots
- software version information
- one final HTML report
