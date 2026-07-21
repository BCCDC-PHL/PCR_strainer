# Validation and Limitations

This document describes the known limitations of PCR_strainer, guidance for validating the tool when deploying it at a new site or for a new assay, and recommendations for its use in the context of laboratory accreditation.

---

## Intended scope

PCR_strainer is a **computational monitoring tool** for assay inclusivity. It is not a diagnostic instrument, does not produce patient results, and does not replace wet-lab assay validation. Its intended use is described in [Overview](01-overview.md).

For accreditation purposes, PCR_strainer should be classified as a laboratory information tool used in quality monitoring, analogous to reagent lot tracking or external quality assurance data review — it supports the quality system but does not itself generate diagnostic results.

---

## Known limitations

### 1. Tm calculations do not account for modified bases

TNTBLAST calculates melting temperature using nearest-neighbour thermodynamics for unmodified DNA oligonucleotides. Common probe modifications substantially alter the actual Tm:

| Modification | Typical effect on Tm |
|---|---|
| Minor groove binder (MGB) | +15 to +20°C |
| Locked nucleic acid (LNA) | +4 to +8°C per LNA base |
| 2′-O-methyl RNA | +1 to +2°C per modified base |

**Practical consequence:** For MGB probes in particular, the calculated Tm reported in `PCR_results.tsv` and the HTML report may be 15–20°C below the actual experimental Tm. A genome that TNTBLAST does not detect (because the calculated Tm falls below the `-m` threshold) may in fact be detected in the wet-lab assay due to the MGB contribution.

**Mitigation:** When running PCR_strainer with MGB or LNA probes, lower the `-m` threshold to account for the expected Tm boost. For example, if your assay uses a -m 50°C threshold for unmodified probes, consider running with `-m 30°C` or lower when the probe is MGB-modified. Review the missed sequences report carefully — sequences missed by TNTBLAST but with low N content may represent genuine detections in the wet-lab assay.

Document which `-m` value you used and the rationale in your quality record.

### 2. Reference genome sequences may not represent local circulation

PCR_strainer reports inclusivity against whatever genome sequences are available in the database at the time of analysis. If strains circulating in your specific patient population are under-represented (e.g., because sequences from your region are rarely deposited), the tool may:

- **Under-report** variant prevalence if local variants are not in the database.
- **Over-report** a concern if a high-prevalence variant in the global database is not actually circulating locally.

**Mitigation:** Where possible, supplement the reference genome set with sequences generated at your laboratory. Correlate PCR_strainer findings with your laboratory's own detection rates and Ct trends. Local surveillance data takes precedence over global database estimates when assessing clinical impact.

### 3. TNTBLAST uses DNA:DNA thermodynamics throughout

For RNA targets (influenza, SARS-CoV-2, RSV, etc.), the actual assay involves an RNA:DNA hybrid at the primer extension step and, in some assay formats, at the probe hybridisation step. RNA:DNA thermodynamics differ from DNA:DNA, particularly for GC-rich sequences. TNTBLAST does not model this distinction.

**Practical consequence:** Tm values may be slightly different from actual experimental values. The effect is generally small (1–3°C) and does not meaningfully affect the identification of mismatched variants, but should be noted when interpreting borderline Tm values near the `-m` threshold.

### 4. Template secondary structure is not modelled

Folded secondary structure in the target RNA/DNA can inhibit primer binding. TNTBLAST aligns oligos to the linear sequence without considering whether that region is accessible in the folded template. Genomes that are detected computationally but fail in the wet-lab assay may have secondary structure occluding the binding site.

### 5. The amplicon size is not constrained

PCR_strainer runs TNTBLAST with a minimum amplicon size of 0 and relies on the Tm threshold alone to filter alignments. If your assay has a very short expected amplicon (< 80 bp), consider whether very long amplicons from non-specific alignments could appear in the results. The `amplicon_range` coordinates in `PCR_results.tsv` can be used to filter results by amplicon size if needed.

### 6. Site-variant notation for degenerate bases

The site-variant notation in `oligo_site_variant` fields may show lowercase characters at positions where the oligo has a degenerate base (e.g., `R`, `Y`, `M`) and the genome has a valid base covered by that code. This is a display limitation of the character comparison used to generate the notation. The mismatch count columns (`*_mismatches`) correctly report these as 0 mismatches because TNTBLAST's thermodynamic model handles degenerate codes correctly. The per-position heatmap in the HTML report also handles this correctly via an IUPAC lookup table. Only the raw `oligo_site_variant` text may be misleading at degenerate positions.

### 7. Probe strand orientation detection relies on sequence complementarity

PCR_strainer detects whether a probe binds the sense or antisense strand by comparing the probe sequence directly to the extracted site and its reverse complement, selecting the orientation with more matching characters. This approach is reliable for probes ≥ 18 bp with typical sequence composition — the margin between correct and incorrect orientation is large (typically > 15 matching characters). For very short probes (< 15 bp) or probes with highly palindromic or repetitive sequences, the detection may be less reliable. In such cases, use `--keep-tntblast-output` and inspect the raw TNTBLAST output to confirm correct orientation.

---

## Validating PCR_strainer at your site

When deploying PCR_strainer for a new assay or at a new site, the following validation steps are recommended.

### Step 1 — Verify installation

Confirm that all software versions are as expected:

```bash
python pcr_strainer.py 2>&1 | head -3
# Should show: PCR_strainer v0.2.x, TNTBLAST vX.X

tntblast --version

python -c "import pandas; print('pandas', pandas.__version__)"
python -c "import matplotlib; print('matplotlib', matplotlib.__version__)"
```

Record the output in your validation documentation.

### Step 2 — Run with a known reference set

Before using PCR_strainer for surveillance, run it against a set of sequences for which you already know the expected results. A suitable reference set would include:

- At minimum 10–20 sequences known to be detected by the assay (e.g., from positive QC material or published validation sequences).
- At minimum 2–5 sequences known to have specific mutations at primer or probe binding sites (if available).
- At minimum 1–2 sequences confirmed NOT to be detected by the assay at the standard Tm threshold (negative controls).

Compare PCR_strainer's output against your expected results and document any discrepancies.

### Step 3 — Verify degenerate base handling

If your assay primers contain IUPAC degenerate bases, construct test cases where:
- A genome has a base covered by the degenerate code → should show as uppercase (match) in the heatmap.
- A genome has a base NOT covered by the degenerate code → should show as lowercase (mismatch).

Confirm both in the HTML report heatmap and in the `fwd_primer_mismatches` / `rev_primer_mismatches` columns of `PCR_results.tsv`.

### Step 4 — Verify probe strand detection

If your assay uses a probe on the antisense strand:
- Run the assay and inspect the `probe_site_seq` values in `PCR_results.tsv`.
- The site sequences should read 5′→3′ in the probe direction and show the expected match/mismatch pattern.
- Use `--keep-tntblast-output` and compare the raw TNTBLAST output's `probe mismatches` count against the mismatch count PCR_strainer records. They should agree.

### Step 5 — Threshold verification

Run with a known genome set and verify that:
- Genomes with mismatches at known positions are captured in the variant report.
- The prevalence percentages in the assay report are consistent with your own count of matching and non-matching sequences.

---

## Monitoring in an accredited setting

### Frequency of monitoring

There is no universal standard for how often inclusivity monitoring should be performed. Consider:

- **Monthly** during periods of active pathogen evolution (e.g., during influenza season, or when new SARS-CoV-2 variants are emerging).
- **Quarterly** for stable pathogens or during low-circulation periods.
- **Triggered review** whenever a significant new variant is reported in the scientific literature or via public health alerts, or whenever unexpected assay performance is observed in routine use.

Define your monitoring schedule in your quality management system and document the rationale.

### Linking PCR_strainer findings to assay performance data

PCR_strainer results are most useful when reviewed alongside:
- **Positivity rates** from routine diagnostic work — a drop in positivity rate among samples that should be positive is a signal.
- **Ct value trends** — rising Ct values in positive samples may indicate reduced assay efficiency in affected strains.
- **External QA results** — proficiency testing panels sometimes include variant strains.

No single indicator is sufficient in isolation. A finding from PCR_strainer becomes more actionable when corroborated by one or more of these additional signals.

### Change management

If PCR_strainer identifies a concern that warrants assay modification:

1. Document the finding and the evidence base (report, genome set, version, date).
2. Obtain updated primers, probes, or reagents from the assay developer or manufacturer, or design and validate in-house modifications.
3. Perform analytical validation of the modified assay before clinical introduction (sensitivity, specificity, limit of detection, cross-reactivity as relevant).
4. Document the validation and the change management process in the quality record.
5. Notify relevant stakeholders (clinical teams, public health authorities) as appropriate.

Do not modify a clinically-used assay based solely on computational analysis, however clear the signal. Wet-lab validation is always required.

### Version control of PCR_strainer itself

PCR_strainer is version-controlled via git. When a new version is deployed:

1. Record the previous and new version numbers in the change log of your quality system.
2. Run a parallel evaluation using both versions on the same genome set to confirm that results are consistent (or document and explain any differences).
3. Update the validation documentation to reflect the new version.

The version used for each run is recorded in the HTML report provenance block and printed to the console at runtime.

---

## Summary checklist for accreditation reviewers

| Item | Location |
|---|---|
| Software version recorded | HTML report provenance block |
| TNTBLAST version recorded | HTML report provenance block |
| Reference genome set documented (database, date, filters, count) | HTML report provenance block + quality record |
| Run parameters documented (`-m`, `-t`) | HTML report provenance block |
| Results reviewed by qualified personnel | Quality record sign-off |
| Findings classified (PASS / CAUTION / ACTION REQUIRED) | HTML report overall verdict |
| Actions taken documented (if CAUTION or ACTION REQUIRED) | Quality record |
| Known limitations acknowledged and addressed in interpretation | This document + quality record |
| Monitoring frequency defined and scheduled | Quality management system |
| Validation of tool at site documented | Validation record |
