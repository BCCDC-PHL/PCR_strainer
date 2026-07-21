# Preparing Input Files

PCR_strainer requires two input files: an **assay CSV** containing your primer and probe sequences, and a **reference genome FASTA** containing the sequences to test against.

---

## The assay CSV

### Format

A plain comma-separated text file with one assay per line and no header row.

**For assays with a probe (qPCR):**
```
assay_name,fwd_primer_name,fwd_primer_seq,rev_primer_name,rev_primer_seq,probe_name,probe_seq
```

**For assays without a probe (conventional PCR or amplicon sequencing):**
```
assay_name,fwd_primer_name,fwd_primer_seq,rev_primer_name,rev_primer_seq
```

### Example

```
FluA_RP1,FluA-RP1-F,AAGCATTCCCAATGACAAACC,FluA-RP1-R,TCTCCCACTTGCAAAGTTGG,FluA-RP1-P,CAGGATCACATATGGGGCCTGTCCCAG
SARS2_E,E-Sarbeco-F,ACAGGTACGTTAATAGTTAATAGCGT,E-Sarbeco-R,ATATTGCAGCAGTACGCACACA,E-Sarbeco-P,ACACTAGCCATCCTTACTGCGCTTCG
FluB_NS,FluB-F,TCCTCAAYTTCAGTCAAGATGTTC,FluB-R,ACTTCAATTCTGCAAAGGTGTTTTC
```

### Rules

**Sequences:**
- All sequences should be entered 5′ → 3′ as written in your assay documentation.
- Probe sequences do **not** need to be pre-reverse-complemented. PCR_strainer automatically determines which strand the probe binds and handles orientation accordingly.
- IUPAC degenerate bases are permitted: A T G C W S M K R Y B V D H N.
- Sequences are treated as DNA. For RNA-based assays (e.g., SARS-CoV-2 RT-qPCR), enter the DNA equivalent — this matches the reference genome sequences, which are also stored as DNA.

**Names:**
- All assay names and all oligo names must be unique across the entire file. Duplicate names will cause an error.
- Names can contain letters, numbers, hyphens, and underscores. Avoid commas (which would break the CSV format) and special characters.
- Assay and oligo names appear in all output files and the HTML report, so use meaningful names that match your assay documentation.

**File encoding:**
- Save as UTF-8. If you create or edit the file in Excel on Windows, save as "CSV UTF-8 (comma delimited)" — not "CSV (comma delimited)", which may omit the encoding marker. PCR_strainer handles UTF-8 BOM automatically if present.

### Adding multiple assays

Each line is one assay. You can include as many assays as needed in a single run. PCR_strainer processes them sequentially — TNTBLAST is run once per assay against the same genome set.

---

## The reference genome FASTA

### Where to download genomes

Reference genomes should come from a public database that is updated regularly:

**NCBI:**
- Influenza: [Influenza Virus Resource](https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi)
- Other pathogens: [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/)
- Download as FASTA. Use filters to select the relevant segment (e.g., HA, NA), host (human), and date range.

**GISAID:**
- SARS-CoV-2 and influenza sequences with richer metadata.
- Requires an account and agreement to the GISAID terms of service.
- **Important:** sequences downloaded from GISAID are subject to access restrictions. Do not distribute reports containing variant sequences from GISAID genomes outside your immediate team without checking the terms of service.

**GenBank sequences (NCBI) carry no redistribution restrictions.**

### Preparing the FASTA file

**FASTA headers must be unique.** PCR_strainer uses the header as the genome identifier throughout all output files. If two sequences share the same header, the tool will exit with an error.

Remove or replace special characters in headers that might cause downstream issues: spaces, slashes, and pipes can cause problems in some tools. NCBI sequences typically have clean headers; GISAID sequences may need normalising.

**All sequences should be on the same strand.** For single-stranded RNA viruses (influenza, SARS-CoV-2, RSV, etc.), all sequences in the FASTA should be the same strand sense — typically the coding/positive sense strand. Mixing strands will produce incorrect results without error messages.

**High-N sequences.** Sequences with a high proportion of ambiguous N bases may fail to be detected even if the assay is inclusive. The missed sequences report flags these separately. Consider filtering out sequences with > 5% Ns before analysis, or interpret the missed sequences report carefully.

### Recommended genome set sizes

| Pathogen | Suggested set | Update frequency |
|---|---|---|
| Influenza A (HA) | All available human sequences, ≥1000 if possible | Monthly during season |
| Influenza B | All available human sequences | Seasonally |
| SARS-CoV-2 | Representative sample by lineage, 5,000–20,000 | Monthly |
| RSV | All available with target segment | Annually or as needed |
| Other | As many as publicly available | Per major outbreak |

Larger sets give more representative prevalence estimates for variant frequencies. For rare or emerging pathogens, even 50–100 sequences are informative.

### Data governance note

Record the database, date of download, any filters applied, and the number of sequences in your quality documentation for each run. This information is captured automatically in the HTML report provenance block, but the full download record should also be maintained separately.
