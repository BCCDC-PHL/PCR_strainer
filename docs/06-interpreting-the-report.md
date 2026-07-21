# Interpreting the Report

This document explains how to read the HTML report and the TSV output files, and what actions to consider for each finding. It is written for scientists and laboratory directors who will use these results to make decisions about assay quality.

---

## The provenance block

The top-right corner of every report shows:

- **Run date** — when the analysis was performed
- **PCR_strainer version** — the software version used
- **TNTBLAST version** — the alignment engine version
- **Reference genomes** — file name and total sequence count
- **Assay file** — file name
- **Min Tm** and **Variant threshold** — the parameters used

**For quality records:** this block should be captured (screenshot or PDF print) along with the overall verdict whenever a report is filed in a quality system. It provides a complete record of what was run, on what data, with which software version.

---

## The overall verdict

The banner at the top of the report shows one of three statuses:

### ✓ PASS — All assays show acceptable inclusivity

All assays have overall inclusivity ≥ 90% and no individual oligo has a mismatch rate ≥ 15% among detected genomes. No immediate action is required.

**Recommended action:** File the report in the quality record. Note the run date and genome set for the audit trail.

---

### ⚠ CAUTION — One or more assays require review

One or more assays have either:
- Overall inclusivity between 75–90%, or
- An oligo with ≥ 15% of detected genomes having at least one mismatch, but without the specific high-risk variant patterns that would trigger ACTION REQUIRED.

**Recommended action:**
1. Review the highlighted assay cards in detail (see below).
2. Assess whether the mismatch is positionally significant (near the 3′ end, multiple errors) — the report provides this interpretation.
3. Determine whether the variant is in strains circulating in your region specifically, rather than globally.
4. Consider increasing monitoring frequency (e.g., monthly instead of quarterly) until the situation is resolved.
5. Document the finding and your assessment in the quality record.

---

### ⚠ ACTION REQUIRED — Significant mismatch prevalence detected

An assay is flagged ACTION REQUIRED when at least one of these conditions is met:
- Overall inclusivity < 75%, or
- A variant with prevalence ≥ 15% that has a mismatch **within the last 3 base pairs of the primer** (near-3′ end, most likely to inhibit Taq extension), or
- A variant with prevalence ≥ 15% that has **2 or more total errors** (multiple mismatches substantially reduce primer or probe binding)

**Recommended action:**
1. Review the specific variant(s) flagged and confirm the finding is genuine (not a data quality artefact).
2. Assess clinical impact: are the affected strains in circulation in your patient population? What is your current positive detection rate?
3. Escalate to the laboratory director and, if appropriate, to the assay developer or reagent supplier.
4. Consider wet-lab sensitivity testing using representative strains carrying the identified variant.
5. Initiate the laboratory's assay change management process if an update is warranted.
6. Do not modify the assay in clinical use without completing formal validation.

---

## The assay summary table

Below the overall verdict, a table lists every assay with:

| Column | What to look at |
|---|---|
| Status badge | PASS / CAUTION / ACTION REQUIRED |
| Detected / Total | Raw numbers. Low detection on a good assay suggests genome quality issues. |
| Inclusivity bar | Visual representation of the percentage |
| Oligos with errors | The worst-performing oligo and its mismatch rate |
| Peak mismatch | The position and oligo with the highest per-position mismatch rate |

Click any row to jump to that assay's detailed card.

---

## Individual assay cards

Each assay has a collapsible card. Cards for CAUTION and ACTION REQUIRED assays are open by default.

### Metric cards

The three metric cards at the top of each card show:

**Inclusivity** — the percentage of all reference genomes detected. This is the primary quality indicator. Interpret in context: a drop from 99% to 95% on an assay that has been stable is more significant than 95% at initial evaluation.

**Perfect match** — the percentage with zero errors across all oligos. The difference between Inclusivity and Perfect match tells you how many detected genomes have at least one error — they are still amplified (or hybridised to the probe), but with reduced efficiency.

**Missed genomes** — the count not detected. If this number is high, check the missed sequences report to determine whether high-N sequences account for most of the misses.

### Error distribution chart

Shows what fraction of detected genomes have 0, 1, 2, or 3+ total errors across all oligos. A healthy assay should have most genomes in the 0-error bar.

### Mismatch rate by oligo chart

Shows what percentage of detected genomes have at least one mismatch in each oligo. This helps identify which part of the assay is under pressure. Note that a genome can have a mismatch in a primer and still be detected — a single peripheral mismatch rarely abolishes amplification. Multiple mismatches, or a mismatch near the 3′ end of a primer, are the primary concerns.

### Per-position heatmap

The most diagnostically useful section. For each oligo, every position is shown as a coloured tile:

| Colour | Mismatch prevalence |
|---|---|
| Green | 0% — no mismatches |
| Yellow | < 3% |
| Orange | 3–10% |
| Pink | 10–20% |
| Red | > 20% |

The Tm range for each oligo is shown below the oligo header (where data is available), formatted as min–max°C with the mean.

**What to look for:**
- **Red or pink tiles near the 3′ end** of the forward or reverse primer are the highest-priority concern. The 3′ end of a primer is critical for Taq polymerase extension; mismatches here are far more likely to cause amplification failure than mismatches at the 5′ end.
- **Red or pink tiles in the probe** affect hybridisation efficiency and may reduce fluorescence signal, particularly for hydrolysis probes.
- **Clusters of adjacent mismatched positions** are more concerning than isolated single positions at the same prevalence.

Hover over any tile to see the exact position number and mismatch percentage.

### Sequence variants table

Lists the specific binding-site sequences found in circulating genomes, sorted by prevalence. Columns:

| Column | What it tells you |
|---|---|
| Oligo | Which primer or probe |
| Variant sequence | The genome sequence at the binding site, using site-variant notation (uppercase = match, lowercase = mismatch, `-` = deletion, `(X)` = insertion) |
| Errors | Total mismatches + gaps |
| Prevalence | % of detected genomes with this exact variant, with the absolute count |
| Tm (°C) | Melting temperature for this variant's specific genomes (single value if all the same, range if they differ) |
| Interpretation | Automated risk assessment: position, Tm implications, and urgency |

**Reading the interpretation text:** The automated interpretation uses three categories — HIGH RISK, MODERATE RISK (no label), and Monitor. These are aids to review, not clinical decisions. A trained scientist should confirm the assessment in the context of local epidemiology and wet-lab data before any action is taken.

---

## Decision guide

| Finding | Typical interpretation | Action |
|---|---|---|
| 0% mismatch at all positions | Assay fully inclusive of current strains | File report, routine monitoring |
| 1–5% mismatch at 5′ end of primer | Minor variant, low positional risk | Document, monitor trend |
| 5–15% mismatch at internal primer position | Moderate concern, watch for increase | Increase monitoring frequency, document |
| ≥ 15% mismatch at 3′ end of primer | High risk of reduced sensitivity in affected strains | Escalate to lab director, consider wet-lab investigation |
| ≥ 15% with ≥ 2 errors on one oligo | High risk regardless of position | Escalate immediately |
| Probe mismatch ≥ 15% | Reduced fluorescence signal in affected strains | Assess with Ct data; escalate if Ct shift detected |
| Inclusivity < 90% overall | Significant proportion of strains not amplified | Investigate missed sequences; correlate with clinical data |
| Inclusivity < 75% | Assay has major inclusivity gap | Action Required protocol; escalate |

---

## Using trends over time

A single report is informative; a series of reports over time is more powerful. Key things to track:

- **Overall inclusivity trend** — slow drift downward may indicate gradual accumulation of mutations across lineages.
- **Variant prevalence over time** — a variant at 2% one month and 8% the next month is more concerning than a stable 10%.
- **New variants appearing** — the first detection of a new binding-site sequence warrants heightened monitoring even at low prevalence.

Maintain a log or spreadsheet of the key metrics from each run (date, genome count, inclusivity per assay, highest mismatch rate) to enable trend analysis.

---

## Reporting to the quality system

For accredited laboratories, we recommend the following minimum documentation for each monitoring run:

1. A copy of the HTML report (or PDF printout) filed in the quality record.
2. A brief written interpretation: overall status for each assay, any specific findings, actions taken or planned.
3. A record of the genome set used (database, date downloaded, filters applied, sequence count).
4. The PCR_strainer and TNTBLAST version numbers (visible in the report provenance block).
5. The name of the scientist who reviewed and interpreted the results.
6. If the status was CAUTION or ACTION REQUIRED: the disposition — what was done, when, and by whom.
