# Troubleshooting

Common errors and how to resolve them.

---

## Installation errors

### `tntblast: command not found`

TNTBLAST is not on your `PATH`. Either it is not installed, or the directory containing the binary is not in your shell's search path.

Check if the binary exists:
```bash
find /usr /opt /home -name "tntblast" 2>/dev/null
```

If found, add its directory to your `PATH`:
```bash
export PATH="/path/to/tntblast/directory:$PATH"
```

Add this line to your `~/.bashrc` or `~/.bash_profile` to make it persistent.

If not found, compile and install TNTBLAST from source: https://github.com/jgans/thermonucleotideBLAST

---

### `ImportError: No module named 'matplotlib'` (or pandas, numpy)

Install the missing package:
```bash
pip install matplotlib
# or
conda install -c conda-forge matplotlib
```

---

## Input errors

### `ERROR: FASTA headers in input genomes file must be unique`

Two or more sequences in your reference FASTA share the same header line. PCR_strainer uses headers as genome identifiers — duplicates are not permitted.

Find the duplicates:
```bash
grep "^>" genomes.fasta | sort | uniq -d
```

Remove or rename the duplicate sequences before rerunning.

---

### `ERROR: Line in assay file not properly formatted`

The assay CSV has a line that does not have exactly 5 or 7 comma-separated fields.

Common causes:
- A comma inside an assay or oligo name — names must not contain commas.
- A missing field — check that all required columns are present.
- Windows line endings — if the file was created on Windows, try converting: `dos2unix assays.csv`
- The file was saved in a format other than CSV — Excel's "CSV UTF-8" option is correct; "CSV (comma delimited)" may also work.

---

### `ERROR: Assay and oligo names must be unique`

Two assays have the same `assay_name`, or two oligos share the same name. All names across the entire assay CSV must be unique — including oligo names across different assays.

---

## Runtime errors

### `ValueError: All arrays must be of the same length`

This error occurs during DataFrame construction in `parse_tntblast_output`. It means that different fields in the TNTBLAST output were found a different number of times than expected, usually because TNTBLAST's output format differs slightly from what the parser expects.

**Most common cause:** TNTBLAST produced no hits for one or more assays, resulting in an empty output file. Check whether the TNTBLAST output file exists and has content:

```bash
python pcr_strainer.py ... --keep-tntblast-output
ls -la results/
cat results/<assay_name>_tntblast_output.txt | head -30
```

If the file is empty, TNTBLAST found no alignments — check your assay sequences, genome file, and Tm threshold.

**Other cause:** A TNTBLAST version with a different output format. Check:
```bash
tntblast --version
```

PCR_strainer has been tested with TNTBLAST 2.4+. Older versions may use different field names.

---

### `ValueError: Cannot set a DataFrame with multiple columns to the single column fwd_primer_site_seq`

This is a pandas 2.x compatibility issue. It occurs when `DataFrame.apply()` infers that the result should be a multi-column DataFrame rather than a Series.

Ensure you are using PCR_strainer v0.2.8 or later, which includes `result_type='reduce'` on all `.apply()` calls to fix this.

Check your version:
```bash
python pcr_strainer.py 2>&1 | head -2
```

---

### `AttributeError: Can only use .str accessor with string values, not floating`

This is also a pandas 2.x compatibility issue. It occurs when `.str.split()` is called on a column that contains NaN values (which have float type).

Ensure you are using PCR_strainer v0.2.8 or later, which uses `.apply(lambda x: x.split(...) if isinstance(x, str) else x)` instead of `.str.split()`.

---

### `TNTBLAST terminated with errors. TNTBLAST error code: N`

TNTBLAST itself exited with an error. Common causes:

- **Corrupt or malformed FASTA file.** Check with:
  ```bash
  grep -c "^>" genomes.fasta          # count sequences
  grep -v "^>" genomes.fasta | grep "[^ATGCNatgcnWSMKRYBVDHwsmkrybvdh-]" | head
  ```
  The second command looks for non-standard characters.

- **Primer sequence contains invalid characters.** TNTBLAST accepts only IUPAC nucleotide codes.

- **Out of memory.** TNTBLAST can use significant memory on large genome sets. Request more memory in your HPC job.

- **TNTBLAST binary not compatible with your system.** Try recompiling from source.

---

### `AttributeError: 'DataFrame' object has no attribute 'applymap'`

This is a pandas 2.1+ compatibility issue in a dependency. `applymap` was renamed to `map` in pandas 2.1 and removed later.

If this occurs in PCR_strainer's own code, update to the latest version. If it occurs in a dependency (e.g., a utility function outside PCR_strainer), locate the file:

```bash
grep -rn "applymap" /path/to/environment/lib/python*/site-packages/
```

and replace `df.applymap(func)` with `df.map(func)` in the relevant file.

---

## Result quality issues

### Probe site sequences look garbled — many dashes and parentheses

Example: `(CTGGGA)CAG(GCCCCGTATGT)GATC-----------CTG------`

This indicates the probe site sequence was aligned in the wrong strand orientation. Ensure you are using PCR_strainer v0.2.9 or later, which uses direct character comparison to determine probe orientation before alignment.

If you have the latest version and still see this, run with `--keep-tntblast-output` and examine the raw TNTBLAST output for the affected assay. Look for:
- The `probe range` and `amplicon range` coordinate lines
- The amplicon FASTA block (the `>` line and sequence below it)

Calculate whether `probe_range[0] - amplicon_range[0]` gives a valid positive index within the amplicon sequence. If the index is negative or exceeds the amplicon length, there is a coordinate calculation issue — please open a GitHub issue with the relevant section of the TNTBLAST output.

---

### Some or all genomes appear in the missed sequences report unexpectedly

**Check 1: Are the sequences high-N?** Open `missed_seqs_report.tsv` and check `perc_Ns`. Sequences with > 5% Ns are likely genuinely low-quality.

**Check 2: Is the Tm threshold too high?** Lower `-m` by 5°C and rerun. If previously missed sequences now appear, the calculated Tm for those genomes was marginal. Consider whether the threshold is appropriate for your assay.

**Check 3: Are the genome sequences on the correct strand?** For single-stranded RNA viruses, all sequences should be positive/coding sense. If some sequences are on the opposite strand, TNTBLAST will not find the forward primer.

**Check 4: Do the primer sequences match what is in your assay documentation?** Even a single wrong base can prevent binding. Double-check the sequences in your assay CSV against the original design document.

---

### All detected genomes have 100% mismatch at positions with degenerate bases

This is expected behaviour for the site-variant notation in `variant_report.tsv` — it is a display limitation, not an error. The mismatch counts in `PCR_results.tsv` will correctly show 0 mismatches for valid degenerate base matches. The per-position heatmap in the HTML report also handles this correctly. See [How It Works — Degenerate base handling](07-how-it-works.md#degenerate-base-handling).

---

### The HTML report does not open correctly from a network share

Open the report from a mapped drive letter (e.g., `Z:\results\run_name_report.html`) rather than a UNC path (`\\server\share\results\run_name_report.html`). Some Windows browsers classify UNC paths as the Internet zone and apply restrictions that affect local HTML rendering.

---

## Debugging with `--keep-tntblast-output`

When a result looks unexpected, retaining the raw TNTBLAST output is the most direct way to diagnose the issue:

```bash
python pcr_strainer.py -a assays.csv -g genomes.fasta -o results/debug \
    --keep-tntblast-output
```

The raw output file `results/<assay_name>_tntblast_output.txt` shows:
- Every genome that was detected and its alignment statistics
- The calculated Tm for each oligo on each genome
- The probe strand (`probe contained in forward strand (+)` or `reverse strand (-)`)
- The raw amplicon sequence for each hit

Compare one representative hit in this file against the corresponding row in `PCR_results.tsv` to verify that the parsing is producing the expected values.

---

## Getting help

- **GitHub Issues:** https://github.com/BCCDC-PHL/PCR_strainer/issues
- Include the PCR_strainer version, TNTBLAST version, the exact error message, and (if relevant) a minimal example of the assay CSV and a few lines of the TNTBLAST raw output file.
