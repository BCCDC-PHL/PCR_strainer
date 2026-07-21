# How It Works

This document describes the algorithms and technical decisions behind PCR_strainer. It is written for bioinformaticians who want to understand what the tool is actually doing, or who need to evaluate whether it is appropriate for a specific use case.

---

## The overall data flow

```
Assay CSV  ──┐
              ├──► TNTBLAST (per assay) ──► Raw text output
Genome FASTA ─┘                                    │
                                                   ▼
                                         parse_tntblast_output()
                                                   │
                                                   ▼
                                         Accumulated DataFrame
                                                   │
                          ┌────────────────────────┼─────────────────────────┐
                          ▼                        ▼                         ▼
                  assay_report.tsv       variant_report.tsv      missed_seqs_report.tsv
                          │
                          ▼
                  PCR_results.tsv ──► amplicons.fasta
                          │
                          ▼
                  pcr_strainer_report.py ──► report.html
```

---

## TNTBLAST — the alignment engine

PCR_strainer wraps [thermonucleotideBLAST (TNTBLAST)](https://github.com/jgans/thermonucleotideBLAST) developed by Jason Gans at Los Alamos National Laboratory.

TNTBLAST differs from standard sequence aligners in that it models oligonucleotide binding thermodynamically. Rather than finding the longest matching subsequence, it asks: **would this primer or probe actually bind this genomic sequence at the specified conditions?**

The key parameters passed to TNTBLAST from PCR_strainer:

- **Minimum Tm** (`-m` flag, default 45°C): TNTBLAST calculates the melting temperature for every candidate alignment using nearest-neighbour thermodynamics. Only alignments meeting or exceeding this threshold are reported. This mimics the annealing step of PCR — if the calculated Tm is below the threshold, the primer would not bind stably enough to prime extension.
- **`--best-match`**: For each genome, only the best-scoring alignment is reported, not all possible alignments. This prevents inflated hit counts from repetitive sequence regions.
- **`-m 0`**: Minimum amplicon size of 0. We want all valid primer pairs reported regardless of amplicon length; PCR_strainer handles size filtering via the Tm threshold.

### What TNTBLAST calculates

For each detected genome, TNTBLAST reports:
- The melting temperature for each oligo's binding at that specific site
- The number of mismatches and gaps for each oligo
- The genomic coordinates of the amplicon and probe binding sites
- The amplicon sequence

It does **not** perform alignment quality scoring beyond Tm and mismatch count. Positional information (is the mismatch near the 3′ end?) is calculated by PCR_strainer during parsing.

### What TNTBLAST does not model

- Modified bases (MGB, LNA, BHQ, FAM): Tm is calculated for unmodified DNA. MGB probes typically have Tm 15–20°C higher than the calculated value. See [Validation and Limitations](08-validation-and-limitations.md).
- RNA:DNA hybrid thermodynamics: for RT-PCR assays, the actual binding thermodynamics involve an RNA template, which differs from DNA:DNA. TNTBLAST uses DNA:DNA parameters throughout.
- Secondary structure in the template: folded regions of the target genome that might inhibit primer binding are not considered.

---

## Parsing TNTBLAST output

TNTBLAST produces a plain text file for each assay. PCR_strainer parses this file to extract a structured DataFrame.

The parser reads the key-value lines (`field = value`) and the FASTA block (`>header` followed by the amplicon sequence) for each hit. Fields are renamed to snake_case column names. Numeric fields (mismatches, gaps, Tm, 3′ clamp) are type-coerced from strings to integers or floats.

For assays without a probe, probe-related fields are set to `NaN` rather than left absent, so all assays have the same column structure.

---

## Extracting oligo site sequences

Once the basic alignment data is parsed, PCR_strainer reconstructs the genome sequence at each oligo binding site to enable the site-variant notation.

### Forward primer

The forward primer binding site is at the 5′ end of the amplicon sequence, so it is extracted directly:

```python
site_length = len(fwd_primer_seq) + fwd_primer_gaps
fwd_primer_site_seq = amplicon_seq[:site_length]
```

### Reverse primer

The reverse primer binds the antisense strand at the 3′ end of the amplicon. The site is extracted from the 3′ end of the amplicon sequence and reverse-complemented:

```python
site_length = len(rev_primer_seq) + rev_primer_gaps
rev_primer_site_seq = rev_comp(amplicon_seq[-site_length:])
```

After reverse complementation, the sequence reads 5′→3′ in the same direction as the primer, so the site-variant notation can be applied consistently.

### Probe

The probe can bind either the sense or antisense strand depending on the assay design. This is the most complex case.

The probe binding coordinates are given relative to the genome, and the amplicon sequence is always on the sense strand. The site is extracted by offset from the amplicon start:

```python
probe_start = probe_range[0] - amplicon_range[0]
probe_length = probe_range[1] - probe_range[0] + probe_gaps + 1
probe_site_seq = amplicon_seq[probe_start : probe_start + probe_length]
```

This gives the sense-strand sequence at the probe position. If the probe binds the antisense strand, this extracted sequence will be the reverse complement of the probe — aligning them directly would produce a garbled result.

**Strand detection:** PCR_strainer determines probe strand orientation by comparing the probe sequence to both the extracted site and its reverse complement using a direct character-by-character count before calling the alignment function:

```python
sense_matches    = sum(a == b for a, b in zip(probe_upper, probe_site_seq.upper()))
antisense_matches = sum(a == b for a, b in zip(probe_upper, rev_comp(probe_site_seq).upper()))
```

Whichever orientation produces more direct character matches is used for the alignment. For a 20–30 bp probe the margin between the correct orientation (~26–29 matches) and the wrong one (~5–12 accidental character matches) is large enough to make this determination unambiguous. The alignment function is then called once on the correctly-oriented sequence.

---

## Site-variant notation

The `write_oligo_site_variant` function generates the compact alignment notation used throughout the output. It works by finding the longest common k-mer between the oligo and site sequences, using these as anchor points, and recursively aligning the regions between anchors.

The result uses four character types:

| Character | Meaning |
|---|---|
| UPPERCASE | Match between oligo and genome at this position |
| lowercase | Mismatch — the oligo and genome differ here |
| `-` | Deletion in the genome (base present in oligo, absent in genome) |
| `(X)` | Insertion in the genome (extra base(s) not in oligo) |

The notation is always written left-to-right in the 5′→3′ direction of the oligo as supplied.

---

## Degenerate base handling

IUPAC degenerate bases in oligo sequences (R, Y, S, W, K, M, B, D, H, V, N) are handled correctly in the per-position mismatch heatmap.

TNTBLAST correctly treats degenerate bases when calculating Tm and mismatch counts — an R (A or G) at a primer position that encounters an A in the genome is counted as a match.

However, the `write_oligo_site_variant` function uses a character comparison that does not account for IUPAC codes: it compares `R` vs `A` as different characters and writes lowercase `a` (mismatch notation) even though A is a valid match for R.

PCR_strainer corrects for this in the per-position heatmap by consulting an IUPAC lookup table: if the oligo has a degenerate base at a position, and the genome base is within the set of bases that degenerate code represents, the position is not counted as a mismatch in the heatmap. All 15 standard IUPAC ambiguity codes are covered:

| Code | Bases covered |
|---|---|
| R | A, G |
| Y | C, T |
| S | G, C |
| W | A, T |
| K | G, T |
| M | A, C |
| B | C, G, T (not A) |
| D | A, G, T (not C) |
| H | A, C, T (not G) |
| V | A, C, G (not T) |
| N | A, T, G, C |

Note: the mismatch counts in `PCR_results.tsv` (`*_mismatches`) come directly from TNTBLAST and are always correct with respect to degenerate bases. Only the site-variant notation field may show degenerate positions as lowercase when they are in fact valid matches.

---

## Status thresholds

The HTML report assigns each assay a status based on the following rules, applied in order:

1. **Start with inclusivity:**
   - `perc_detected ≥ 90%` → PASS
   - `75% ≤ perc_detected < 90%` → CAUTION
   - `perc_detected < 75%` → ACTION REQUIRED

2. **Upgrade PASS → CAUTION** if any oligo has a mismatch rate ≥ 15% (i.e., 15% or more of detected genomes have at least one mismatch in that oligo).

3. **Upgrade to ACTION REQUIRED** if any variant with prevalence ≥ 15% has:
   - A mismatch within the **last 3 base pairs of the primer** (near-3′ position), or
   - **2 or more total errors** (mismatches + gaps)

The near-3′ rule captures the biologically important case where a high-prevalence variant is concentrated at the position most critical for Taq polymerase extension. Two or more errors substantially reduce binding stability regardless of position.

These thresholds are configurable via `--pass-threshold` and `--caution-threshold` in `pcr_strainer_report.py` for standalone use.

---

## The HTML report generator

`pcr_strainer_report.py` is a separate module that reads the four TSV files produced by the main pipeline and generates a self-contained HTML report. It can be run standalone after a PCR_strainer run if you need to regenerate the report without re-running TNTBLAST.

Key design decisions:

- **No JavaScript.** The report is fully static HTML and CSS. Charts are rendered server-side by matplotlib as inline SVG. This makes the report safe to open on managed Windows PCs where script execution in local HTML files may be blocked by endpoint protection.
- **No external dependencies.** No CDN calls, no fonts loaded from the internet. The file is completely self-contained.
- **CSP header.** The HTML includes a strict Content Security Policy: `default-src 'none'; style-src 'unsafe-inline'; img-src data:`. This prevents any injected content from executing, even if the report is somehow tampered with.
- **File permissions.** Output files are created with mode `0o640` (owner read/write, group read, no world access) using `os.open()` with the permission set atomically at creation time, avoiding a race window on shared HPC filesystems.
