# Running the Pipeline

## Basic usage

```bash
python pcr_strainer.py \
    -a assays.csv \
    -g genomes.fasta \
    -o results/run_name
```

This runs all assays in `assays.csv` against `genomes.fasta` and writes output files to the `results/` directory, all prefixed with `run_name`.

---

## All arguments

### Required

| Flag | Description |
|---|---|
| `-a` | Path to the assay CSV file |
| `-g` | Path to the reference genome FASTA file |
| `-o` | Output prefix: directory path plus base name for output files. For example, `-o results/flu_2024_Q1` writes files named `results/flu_2024_Q1_assay_report.tsv`, `results/flu_2024_Q1_report.html`, etc. The directory must already exist. |

### Optional

| Flag | Default | Description |
|---|---|---|
| `-m` | 45 | Minimum melting temperature (°C) for primers and probes. TNTBLAST only reports alignments where the calculated Tm meets or exceeds this threshold. Lowering this will detect more genomes but may include non-specific alignments. |
| `-t` | 0 | Minimum prevalence threshold (%) for variants to appear in the variant report and HTML summary. Setting `-t 1` suppresses variants seen in fewer than 1% of detected genomes, reducing noise in large datasets. Must be > 0 and < 100 if set. |
| `-p` | 1 | Molar concentration of primer oligos in µM. Used in TNTBLAST's thermodynamic calculation. Change only if your assay uses non-standard concentrations. |
| `-P` | 1 | Molar concentration of probe oligos in µM. |

### Diagnostic

| Flag | Description |
|---|---|
| `--keep-tntblast-output` | Retain the raw TNTBLAST output text file for each assay (`<assay_name>_tntblast_output.txt`) in the output directory after parsing. By default these are deleted. Useful for diagnosing unexpected alignment results — see [Troubleshooting](08-troubleshooting.md). |

---

## Worked examples

### Monthly influenza surveillance run

```bash
python pcr_strainer.py \
    -a assays/influenza_assays.csv \
    -g genomes/influenza_HA_$(date +%Y%m).fasta \
    -o results/flu_HA_$(date +%Y%m) \
    -m 50 \
    -t 1
```

- `-m 50` raises the Tm threshold slightly above the default, appropriate for well-designed influenza assays and reducing false hits from low-complexity sequence regions.
- `-t 1` suppresses variants seen in fewer than 1% of detected genomes.

### SARS-CoV-2 assay check with all variants reported

```bash
python pcr_strainer.py \
    -a assays/sars2_assays.csv \
    -g genomes/sars2_sequences.fasta \
    -o results/sars2_$(date +%Y%m%d) \
    -m 45 \
    -t 0
```

`-t 0` reports every variant regardless of prevalence. Useful for a full initial evaluation but produces a larger report.

### Diagnostic run with raw TNTBLAST output retained

```bash
python pcr_strainer.py \
    -a assays/problem_assay.csv \
    -g genomes/test_genomes.fasta \
    -o results/debug_run \
    --keep-tntblast-output
```

The raw TNTBLAST output files (`<assay>_tntblast_output.txt`) will remain in `results/` after the run for manual inspection.

---

## What happens during a run

PCR_strainer processes each assay in sequence. For each assay you will see output like:

```
Running TNTBLAST for FluA_RP1 against genomes/influenza_HA_202401.fasta...
Fwd primer seq: AAGCATTCCCAATGACAAACC
Rev primer seq: TCTCCCACTTGCAAAGTTGG
Probe seq: CAGGATCACATATGGGGCCTGTCCCAG

Parsing TNTBLAST output from FluA_RP1 ...

Writing assay report...
Writing variant report...
Writing missed sequences report...
Writing PCR results...
Writing amplicon FASTA files...
  FluA_RP1: 4823 amplicons → results/run_name_FluA_RP1_amplicons.fasta

Done.
```

The HTML report is generated last, after all assays have been processed.

---

## Runtime expectations

Runtime depends primarily on genome set size and number of assays. TNTBLAST is the computationally intensive step.

| Genome set size | 1 assay | 5 assays |
|---|---|---|
| ~1,000 genomes | < 1 min | ~3 min |
| ~5,000 genomes | ~3 min | ~15 min |
| ~20,000 genomes | ~15 min | ~75 min |

These are approximate figures on a modern HPC node. Runtime scales roughly linearly with genome count and assay count.

---

## Output directory structure

After a run, your output directory will contain:

```
results/
├── run_name_assay_report.tsv
├── run_name_variant_report.tsv
├── run_name_missed_seqs_report.tsv
├── run_name_PCR_results.tsv
├── run_name_FluA_RP1_amplicons.fasta     ← one per assay
├── run_name_SARS2_E_amplicons.fasta
└── run_name_report.html
```

See [Output Reference](05-output-reference.md) for a full description of each file.

---

## Automating regular runs

For scheduled surveillance, PCR_strainer can be run as a cron job or via a workflow manager (Nextflow, Snakemake). A minimal shell script wrapper:

```bash
#!/bin/bash
# monthly_flu_check.sh
set -euo pipefail

DATE=$(date +%Y%m)
GENOME_DIR=/data/reference_genomes
ASSAY_DIR=/data/assays
OUT_DIR=/data/pcr_strainer_results

# Download latest genomes (example using NCBI datasets)
# datasets download virus genome taxon 11520 --host human ...

python pcr_strainer.py \
    -a ${ASSAY_DIR}/influenza_assays.csv \
    -g ${GENOME_DIR}/influenza_HA_${DATE}.fasta \
    -o ${OUT_DIR}/flu_HA_${DATE} \
    -m 50 \
    -t 1

# Copy report to network share
cp ${OUT_DIR}/flu_HA_${DATE}_report.html /mnt/reports/
```

Schedule this with cron or your HPC scheduler. Archive previous runs — do not overwrite; the date-stamped output prefix provides a natural audit trail.
