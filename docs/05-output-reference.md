# Output Reference

PCR_strainer produces six output files per run, all sharing the prefix given to `-o`. This document describes the contents and meaning of every file and column.

---

## Summary of output files

| File | Description |
|---|---|
| `<prefix>_assay_report.tsv` | Assay-level inclusivity summary — start here |
| `<prefix>_variant_report.tsv` | Sequences of primer/probe site variants above the reporting threshold |
| `<prefix>_missed_seqs_report.tsv` | Reference genomes not detected by each assay |
| `<prefix>_PCR_results.tsv` | Full per-genome alignment results |
| `<prefix>_<assay>_amplicons.fasta` | Amplicon sequences for detected genomes (one file per assay) |
| `<prefix>_report.html` | Self-contained HTML summary report |

---

## assay_report.tsv

This is the first file to check after a run. It gives the overall inclusivity picture for each assay.

One row per assay per **total error level** (number of mismatches + gaps summed across all oligos). To get the headline inclusivity for an assay, filter to `total_errors == 0` or look at the `perc_detected` column, which is the same for all rows of the same assay.

| Column | Type | Description |
|---|---|---|
| `assay_name` | string | Assay identifier from the assay CSV |
| `total_targets` | int | Total number of sequences in the reference FASTA |
| `detected_targets` | int | Number of sequences detected by TNTBLAST at or above the Tm threshold |
| `perc_detected` | float | `detected_targets / total_targets × 100`. This is the overall inclusivity of the assay. |
| `total_errors` | int | Total mismatches + gaps across all oligos for this row |
| `target_count` | int | Number of detected sequences at this error level |
| `perc_of_detected` | float | `target_count / detected_targets × 100` |
| `perc_of_total` | float | `target_count / total_targets × 100` |

**Reading this file:** For each assay, you will see multiple rows — one for `total_errors = 0` (perfect matches), one for `total_errors = 1`, and so on. The sum of `perc_of_detected` across all rows for one assay should be approximately 100%.

---

## variant_report.tsv

Describes the specific sequences found at each oligo binding site for variants that exceeded the reporting threshold (`-t` flag, default 0%).

One row per unique oligo site sequence per assay. Only variants with at least one error (mismatch or gap) are included.

| Column | Type | Description |
|---|---|---|
| `assay_name` | string | Assay identifier |
| `oligo` | string | Which oligo: `fwd_primer`, `rev_primer`, or `probe` |
| `oligo_name` | string | Oligo identifier from the assay CSV |
| `oligo_seq` | string | Designed oligo sequence (5′→3′) |
| `total_targets` | int | Total reference sequences |
| `detected_targets` | int | Sequences detected by this assay |
| `perc_detected` | float | Overall assay inclusivity (%) |
| `oligo_site_variant` | string | Genome sequence at the binding site in site-variant notation (see below) |
| `oligo_errors` | int | Mismatches + gaps for this variant |
| `target_count` | int | Number of genomes with this exact variant sequence |
| `perc_of_detected` | float | `target_count / detected_targets × 100` |
| `perc_of_total` | float | `target_count / total_targets × 100` |

### Site-variant notation

The `oligo_site_variant` column uses a compact notation to describe how the genome sequence differs from the oligo at the binding site:

| Character | Meaning |
|---|---|
| `UPPERCASE` | This position matches the oligo exactly |
| `lowercase` | This position is a mismatch — the oligo has a different base |
| `-` | Deletion in the genome at this position (base present in oligo, absent in genome) |
| `(X)` | Insertion in the genome at this position (extra base(s) in genome not in oligo) |

**Example:** `ATGCATGCATGcATGCATGCAT` — a 20-base primer with a single mismatch at position 12 (lowercase `c`).

**Probe sequences** are always shown 5′→3′ relative to the probe as written in your assay CSV, regardless of which strand the probe binds. PCR_strainer detects probe strand orientation automatically.

**Degenerate bases** in the oligo (e.g., `R` = A or G) are handled correctly. A genome base that is a valid match for a degenerate code is shown as uppercase (match), even though the characters differ. Only genuinely mismatching bases — those not covered by the degenerate code — are shown as lowercase.

---

## missed_seqs_report.tsv

Lists reference sequences that were not detected by each assay. A sequence is missed when TNTBLAST finds no alignment meeting the Tm threshold, which can happen for two reasons: genuine assay failure (mutations at binding sites), or poor sequence quality (high N content).

| Column | Type | Description |
|---|---|---|
| `assay_name` | string | Assay identifier |
| `target` | string | FASTA header of the missed sequence |
| `target_length` | int | Length of the sequence in base pairs |
| `total_Ns` | int | Count of ambiguous N bases |
| `perc_Ns` | float | `total_Ns / target_length × 100` |

**Interpretation:** Missed sequences with > 5% Ns are likely poor-quality sequences that would not reliably amplify in wet-lab conditions either. Missed sequences with low N content warrant investigation — they may represent genuine assay inclusivity gaps.

---

## PCR_results.tsv

The most detailed output. One row per detected genome per assay. This is the source data for the HTML report and amplicon FASTA files.

### Top-level columns

| Column | Description |
|---|---|
| `assay_name` | Assay identifier |
| `target` | FASTA header of the detected genome |
| `total_errors` | Sum of mismatches + gaps across all oligos |
| `min_3prime_clamp` | Length in bp of the exact-match run at the 3′ end of the forward primer for this alignment. Longer values indicate stronger anchoring for Taq extension. |

### Per-oligo columns

These columns are repeated for each oligo, prefixed with `fwd_primer_`, `rev_primer_`, and `probe_`:

| Column | Description |
|---|---|
| `*_name` | Oligo identifier from the assay CSV |
| `*_seq` | Designed oligo sequence (5′→3′) |
| `*_site_seq` | Genome sequence at the binding site in site-variant notation (see above) |
| `*_mismatches` | Number of base mismatches |
| `*_gaps` | Number of insertion/deletion events |
| `*_errors` | `mismatches + gaps` |
| `*_tm` | Melting temperature (°C) calculated by TNTBLAST for this specific genome's binding site |

`probe_mismatches`, `probe_gaps`, `probe_errors`, and `probe_tm` are `NaN` for assays without a probe.

---

## Amplicon FASTA files

One file per assay: `<prefix>_<assay_name>_amplicons.fasta`

Contains the full amplicon sequence for every detected genome — from the start of the forward primer binding site to the end of the reverse primer binding site, inclusive, as reported by TNTBLAST.

FASTA header format:
```
><target_id> assay=<assay_name> errors=<total_errors>
```

**Uses:**
- Phylogenetic analysis of the amplicon region
- Multiple sequence alignment to visualise variant positions
- Filtering by error count to extract perfect-match amplicons

**Extracting perfect-match amplicons:**
```bash
awk '/errors=0/{p=1; print; next} /^>/{p=0} p' \
    run_name_FluA_RP1_amplicons.fasta > run_name_FluA_RP1_perfect.fasta
```

---

## HTML report

A self-contained HTML file with no external dependencies. Opens in any browser without JavaScript enabled. Safe to open on managed Windows PCs.

See [Interpreting the Report](06-interpreting-the-report.md) for a detailed guide to reading the HTML report.

**Technical notes:**
- Charts are rendered as inline SVG by matplotlib — crisp at any zoom level or print size.
- The report reads from the four TSV files described above. If the TSVs are regenerated, re-run `pcr_strainer_report.py` separately or re-run the full pipeline to update the report.
- The provenance block in the report header records the run date, PCR_strainer version, TNTBLAST version, genome file name, assay file name, and parameter settings used.
