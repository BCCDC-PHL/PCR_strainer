# PCR_strainer — Overview

## What this tool does

PCR_strainer answers one question: **does the primer and probe set for a given assay still bind reliably to the pathogen genomes currently circulating in the population?**

It takes your assay design (primer and probe sequences) and a set of reference genome sequences (typically downloaded from NCBI or GISAID) and systematically checks each genome for mismatches at every oligo binding site. The results are summarised in tabular files and a self-contained HTML report that can be shared with laboratory staff without any software installed on the viewing machine.

---

## Why this matters for accredited assays

Diagnostic PCR assays are validated against pathogen sequences available at the time of design. Viruses and bacteria mutate. Over time, mutations can accumulate at primer or probe binding sites in circulating strains, reducing the sensitivity of the assay — sometimes without any obvious change in routine performance indicators until a cluster of false negatives is identified.

Regular inclusivity monitoring using PCR_strainer provides:

- **Evidence that your assay continues to perform as validated** — documented surveillance forms part of the ongoing quality record for an accredited test.
- **Early warning of emerging variants** — a shift from 0% to 5% prevalence of a probe-site mutation in circulating sequences warrants monitoring, even if clinical sensitivity is not yet visibly affected.
- **Actionable intelligence** — the report identifies exactly which position in which oligo is mutated, how prevalent that mutation is, and whether the position is near the 3′ end where it is most likely to affect amplification.

---

## What this tool is not

PCR_strainer is a **computational surveillance tool**, not a diagnostic instrument. Its outputs are based on available reference genome sequences, which may not fully represent all strains circulating in any given region or patient population.

Results from PCR_strainer:

- **Do not replace wet-lab assay validation.** Any assay modification prompted by this tool must go through your laboratory's established change management and validation process before clinical use.
- **Should be reviewed by qualified laboratory personnel** before any action is taken. The tool flags variants for review; it does not make clinical decisions.
- **Depend on the completeness of the reference genome set.** If circulating strains in your region are under-represented in public databases, variants may be missed.

---

## Intended use

PCR_strainer is intended for use by laboratory scientists and bioinformaticians in public health and clinical microbiology settings to:

1. **Monitor** the inclusivity of in-use diagnostic assays on a scheduled basis (e.g., monthly or at each pathogen database update).
2. **Investigate** unexpected assay performance issues, including potential sensitivity losses in specific outbreak strains.
3. **Evaluate** candidate assay designs before introducing a new test, by assessing how well candidate primers and probes cover the diversity of available reference sequences.
4. **Document** assay inclusivity as part of the quality record for accredited tests.

---

## Documents in this folder

| Document | Audience | Contents |
|---|---|---|
| **01-overview.md** (this file) | All | What the tool does and why |
| **02-installation.md** | Bioinformaticians | Software setup |
| **03-preparing-inputs.md** | Bioinformaticians, scientists | How to prepare assay and genome files |
| **04-running-the-pipeline.md** | Bioinformaticians | All command-line options, HPC usage |
| **05-output-reference.md** | Bioinformaticians, scientists | Detailed description of every output file and column |
| **06-interpreting-the-report.md** | Scientists, lab directors, QA | How to read results and decide what to do |
| **07-how-it-works.md** | Bioinformaticians | Algorithms, thermodynamic model, probe handling |
| **08-validation-and-limitations.md** | QA managers, accreditation | Validation approach, known limitations, QC |
| **09-troubleshooting.md** | Bioinformaticians | Common errors and solutions |

---

## Quick summary for lab directors

PCR_strainer runs automatically on your HPC at BCCDC-PHL. For each assay, it:

1. Downloads or reads a set of reference genomes (e.g., all available influenza HA sequences from NCBI).
2. Virtually "tests" each genome with your primer and probe set, recording whether the assay would amplify it and whether there are any mismatches at the binding sites.
3. Produces an HTML report showing overall inclusivity (% of genomes detected), which positions in the primers and probe have mismatches in circulating strains, and an automated flag — **PASS**, **CAUTION**, or **ACTION REQUIRED** — with a plain-language explanation.

The report is a single HTML file. No software is needed to open it — just a web browser. It can be attached to a quality record directly.
