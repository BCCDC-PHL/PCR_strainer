#!/usr/bin/env python3
"""
pcr_strainer_report.py
======================
Generate a self-contained, fully static HTML inclusivity report from
PCR_strainer TSV output files.

The output HTML contains no JavaScript.  Charts are rendered server-side
by matplotlib and embedded as base64-encoded PNG images.  The report opens
correctly with JavaScript disabled and will not trigger AV/endpoint-protection
alerts on managed Windows PCs.

Standalone usage (run after PCR_strainer):
    python pcr_strainer_report.py -o results/my_run [options]

Integration usage (drop into PCR_strainer's main() AFTER the write_* calls):
    from pcr_strainer_report import write_html_report
    write_html_report(
        name, out_path,
        genome_file=args['-g'],
        assay_file=args['-a'],
        min_tm=args['-m'],
        variant_threshold=args['-t'],
        tool_version=version,
    )

Dependencies
------------
    pandas, numpy     — already required by PCR_strainer
    matplotlib        — standard in scientific Python environments (Anaconda,
                        conda-forge, pip).  Uses the Agg (non-display) backend
                        so no graphical environment is needed on the HPC head node.

Windows network share / browser compatibility
---------------------------------------------
Open the report from a mapped drive letter (e.g. Z:\\...) rather than a UNC
path (\\\\server\\share\\...).  Some browsers classify UNC paths as the
Internet zone, which can affect rendering of local HTML files.

Data governance
---------------
If your reference genomes were downloaded from GISAID, the variant sequences
in this report are derived from GISAID-restricted data.  Review GISAID's terms
of service before distributing this report beyond your immediate team.
Genomes sourced from NCBI/GenBank are not subject to this restriction.

Author: BCCDC-PHL
"""

import argparse     # command-line argument parsing
import base64       # encode matplotlib PNG output for inline embedding
import io           # BytesIO buffer for matplotlib figure output
import logging      # structured logging with file handler for audit trail
import os           # file I/O, permissions (os.open), path handling
import re           # regex for input validation of sequences and names
import sys          # sys.argv, sys.exit, sys.stdout for CLI and logging
import textwrap     # dedent multi-line CLI help strings
from dataclasses import dataclass  # ReportConfig settings object
from datetime import datetime       # run date stamp in provenance block
from pathlib import Path            # path resolution and traversal checks
from typing import Dict, List, Optional         # type annotations on public functions

import matplotlib                # charting — Agg backend avoids display requirement
matplotlib.use('Agg')            # must be set before pyplot is imported
import matplotlib.pyplot as plt  # figure/axes API for chart generation
import numpy as np               # zero-initialised mismatch count arrays (np.zeros)
import pandas as pd              # TSV loading, DataFrame groupby operations


# ── Logging ───────────────────────────────────────────────────────────────────

log = logging.getLogger('pcr_strainer_report')


def _configure_logging(log_file: Optional[str] = None) -> None:
    """Configure stream + optional file logging."""
    log.setLevel(logging.DEBUG)
    fmt = logging.Formatter('%(asctime)s  %(levelname)-8s  %(message)s',
                            datefmt='%Y-%m-%d %H:%M:%S')
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.INFO)
    sh.setFormatter(fmt)
    log.addHandler(sh)
    if log_file:
        fh = logging.FileHandler(log_file, encoding='utf-8')
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        log.addHandler(fh)


# ── Configuration dataclass ───────────────────────────────────────────────────

@dataclass
class ReportConfig:
    """All settings needed to build one report."""
    genome_file:       Optional[str]   = None
    assay_file:        Optional[str]   = None
    min_tm:            Optional[float] = None
    variant_threshold: Optional[float] = None
    tool_version:      str             = 'PCR_strainer'
    pass_threshold:    float           = 90.0
    caution_threshold: float           = 75.0


# ── Constants ─────────────────────────────────────────────────────────────────

INCLUSIVITY_PASS    = 90.0
INCLUSIVITY_CAUTION = 75.0
OLIGO_CAUTION       = 5.0
OLIGO_ACTION        = 15.0

# Output file permissions: owner rw, group r, no world access.
_OUTPUT_FILE_MODE = 0o640

# TSV size limit — conservative for a shared HPC head node.
_MAX_TSV_BYTES = 100 * 1024 * 1024   # 100 MB

# Status colours (hex) used in both Python chart generation and CSS.
_STATUS_COLOUR = {
    'pass':    '#2E7D46',
    'caution': '#B07800',
    'action':  '#C0392B',
}
_STATUS_COLOUR_LIGHT = {
    'pass':    '#A8D4B6',
    'caution': '#DEB860',
    'action':  '#E98080',
}


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate a static HTML report from PCR_strainer outputs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Examples
            --------
            # Standard usage
            python pcr_strainer_report.py -o results/run1 \\
                -g genomes.fasta -a assays.csv --min-tm 50

            # With persistent audit log
            python pcr_strainer_report.py -o results/run1 --log-file run1.log
        """)
    )
    p.add_argument('-o', '--output-prefix', default=None,
                   help='Output prefix used with PCR_strainer -o flag '
                        '(e.g. results/my_run). Reads <prefix>_*.tsv and '
                        'writes <prefix>_report.html.')
    p.add_argument('-g', '--genome-file',  default=None,
                   help='Reference genome FASTA (optional; provenance only).')
    p.add_argument('-a', '--assay-file',   default=None,
                   help='Assay CSV file (optional; provenance only).')
    p.add_argument('--min-tm', type=float, default=None,
                   help='Minimum Tm threshold used (degrees C).')
    p.add_argument('--variant-threshold', type=float, default=None,
                   help='Variant prevalence threshold used (%%).')
    p.add_argument('--pass-threshold', type=float, default=INCLUSIVITY_PASS,
                   help=f'Inclusivity %% >= this is PASS (default {INCLUSIVITY_PASS})')
    p.add_argument('--caution-threshold', type=float, default=INCLUSIVITY_CAUTION,
                   help=f'Inclusivity %% >= this is CAUTION, below is ACTION '
                        f'(default {INCLUSIVITY_CAUTION})')
    p.add_argument('--tool-version', default=None,
                   help='PCR_strainer version string for provenance.')
    p.add_argument('--log-file', default=None,
                   help='Optional path for a persistent DEBUG-level log file.')
    p.add_argument('--file-mode', type=lambda x: int(x, 8),
                   default=_OUTPUT_FILE_MODE,
                   help=f'Octal permission mode for the output HTML file '
                        f'(default {oct(_OUTPUT_FILE_MODE)}).  '
                        f'Example: --file-mode 0o600 for owner-only.')
    return p.parse_args(argv)


# ── Path validation ───────────────────────────────────────────────────────────

def _validate_output_prefix(prefix: str) -> Path:
    """Resolve the output prefix; reject null bytes; warn on path traversal."""
    if '\x00' in prefix:
        raise ValueError('Output prefix contains null bytes.')
    resolved = Path(prefix).resolve()
    cwd      = Path.cwd().resolve()
    try:
        resolved.relative_to(cwd)
    except ValueError:
        log.warning(
            'Resolved output path %s is outside the current working directory %s. '
            'If intentional (e.g. an absolute scratch path) ignore this warning.',
            resolved, cwd
        )
    return resolved


# ── File loading ──────────────────────────────────────────────────────────────

def _tsv(path: str) -> pd.DataFrame:
    """Load a TSV, returning an empty DataFrame on missing or oversized file."""
    if not os.path.isfile(path):
        log.warning('Expected file not found: %s', path)
        return pd.DataFrame()
    size = os.path.getsize(path)
    if size > _MAX_TSV_BYTES:
        log.warning('%s is %.0f MB, exceeding the %.0f MB safety limit. Skipping.',
                    path, size / 1e6, _MAX_TSV_BYTES / 1e6)
        return pd.DataFrame()
    try:
        df = pd.read_csv(path, sep='\t', dtype=str)
        log.debug('Loaded %s (%d rows)', path, len(df))
        return df
    except Exception as exc:
        log.warning('Could not read %s: %s', path, exc)
        return pd.DataFrame()


def _coerce_numeric(df: pd.DataFrame, cols: List[str]) -> pd.DataFrame:
    for col in cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    return df


def load_outputs(prefix: str) -> Dict[str, pd.DataFrame]:
    """Load the four TSV files produced by PCR_strainer for the given prefix."""
    files = {
        'assay_report':  f'{prefix}_assay_report.tsv',
        'variant_report':f'{prefix}_variant_report.tsv',
        'missed_seqs':   f'{prefix}_missed_seqs_report.tsv',
        'pcr_results':   f'{prefix}_PCR_results.tsv',
    }
    dfs = {key: _tsv(path) for key, path in files.items()}
    _coerce_numeric(dfs['assay_report'],
                    ['total_targets','detected_targets','perc_detected',
                     'total_errors','target_count','perc_of_detected','perc_of_total'])
    _coerce_numeric(dfs['variant_report'],
                    ['oligo_errors','target_count','perc_of_detected','perc_of_total',
                     'total_targets','detected_targets','perc_detected'])
    _coerce_numeric(dfs['missed_seqs'],
                    ['total_Ns','target_length','perc_Ns'])
    _coerce_numeric(dfs['pcr_results'],
                    ['total_errors',
                     'fwd_primer_mismatches','fwd_primer_gaps','fwd_primer_errors',
                     'rev_primer_mismatches','rev_primer_gaps','rev_primer_errors',
                     'probe_mismatches','probe_gaps','probe_errors'])
    return dfs


# ── Input validation ──────────────────────────────────────────────────────────

_OLIGO_SEQ_RE  = re.compile(r'^[ATGCWSMKRYBVDHNatgcwsmkrybvdhn()\-]+$')
_OLIGO_NAME_RE = re.compile(r'^[\w\s\-_./:]{1,200}$')
_MAX_OLIGO_LEN  = 200
_MAX_NAME_LEN   = 200


def _sanitise_seq(s: str, field: str = 'sequence') -> str:
    """Validate a nucleotide sequence; return empty string if invalid."""
    if not isinstance(s, str):
        return ''
    s = s.strip()
    if len(s) > _MAX_OLIGO_LEN:
        log.warning('%s is unusually long (%d chars) — truncating to %d.',
                    field, len(s), _MAX_OLIGO_LEN)
        s = s[:_MAX_OLIGO_LEN]
    if s and not _OLIGO_SEQ_RE.match(s):
        log.warning('%s contains unexpected characters and will be omitted: %r',
                    field, s[:40])
        return ''
    return s


def _sanitise_name(s: str, field: str = 'name') -> str:
    """Validate a name string; strip dangerous characters if suspicious."""
    if not isinstance(s, str):
        return '(unknown)'
    s = s.strip()[:_MAX_NAME_LEN]
    if s and not _OLIGO_NAME_RE.match(s):
        log.warning('%s contains unexpected characters (will be escaped in output): %r',
                    field, s[:40])
        s = re.sub(r'[^\w\s\-_./:()]', '', s)[:_MAX_NAME_LEN]
    return s


# ── Degenerate base handling ──────────────────────────────────────────────────
#
# PCR_strainer's write_oligo_site_variant() does a character-by-character
# comparison when building the site_seq field.  This means a degenerate base
# in the oligo (e.g. R = A or G) is compared literally to the genome base
# (e.g. A), which differ as characters, so PCR_strainer writes lowercase 'a'
# — its mismatch notation — even though A is a valid match for R.
#
# TNTBLAST itself handles degenerate oligo bases correctly (it reports 0
# mismatches for R vs A), so the mismatch count columns in PCR_results.tsv
# are accurate.  The site_seq display field is the source of the discrepancy.
#
# The lookup table below is used in parse_site_seq_positions() to detect
# when a "mismatch" in site_seq is actually a valid degenerate match, and
# to avoid counting it as an error in the per-position heatmap.

_IUPAC_MATCHES: Dict[str, set] = {
    'A': {'A'},
    'T': {'T'},
    'G': {'G'},
    'C': {'C'},
    'R': {'A', 'G'},          # puRine
    'Y': {'C', 'T'},          # pYrimidine
    'S': {'G', 'C'},          # Strong
    'W': {'A', 'T'},          # Weak
    'K': {'G', 'T'},          # Keto
    'M': {'A', 'C'},          # aMino
    'B': {'C', 'G', 'T'},    # not A (B comes after A)
    'D': {'A', 'G', 'T'},    # not C (D comes after C)
    'H': {'A', 'C', 'T'},    # not G (H comes after G)
    'V': {'A', 'C', 'G'},    # not T (V comes after T)
    'N': {'A', 'T', 'G', 'C'},  # aNy
}


# ── Per-position mismatch parsing ─────────────────────────────────────────────

def parse_site_seq_positions(oligo_seq: str, site_seq: str) -> List[bool]:
    """
    Parse PCR_strainer's oligo_site_variant notation and return a boolean list
    of length len(oligo_seq) where True indicates a genuine mismatch or deletion.

    Notation in site_seq:
      UPPERCASE  = exact character match with oligo
      lowercase  = character differs from oligo (may still be a valid degenerate match)
      -          = deletion in genome (base present in oligo, absent in genome)
      (X)        = insertion in genome (does NOT consume an oligo position)

    Degenerate base handling:
      When the site_seq has a lowercase letter at a position, PCR_strainer has
      flagged it as a mismatch based on a character comparison.  However, if the
      oligo has a degenerate base at that position (e.g. R, Y, M, K ...), the
      genome base may still be a valid match.  We look up the oligo base in
      _IUPAC_MATCHES and only mark it as a true error if the genome base is not
      in the set of bases the degenerate code represents.

      Example:
        oligo    = ...R...   (R covers A and G)
        site_seq = ...a...   (PCR_strainer writes lowercase — "mismatch")
        result   = False     (A IS covered by R → not a real error)

        oligo    = ...R...
        site_seq = ...c...
        result   = True      (C is NOT covered by R → genuine mismatch)
    """
    if not isinstance(oligo_seq, str) or not isinstance(site_seq, str):
        return [False] * (len(oligo_seq) if isinstance(oligo_seq, str) else 0)

    n         = len(oligo_seq)
    errors    = [False] * n
    oligo_pos = 0
    i         = 0

    for _ in range(len(site_seq) + n):   # iteration cap against pathological input
        if i >= len(site_seq) or oligo_pos >= n:
            break
        c = site_seq[i]
        if c == '(':
            j = site_seq.find(')', i)
            i = (j + 1) if j != -1 else len(site_seq)
            # insertions do not consume an oligo position
        elif c == '-':
            errors[oligo_pos] = True   # deletion — always a real error
            oligo_pos += 1;  i += 1
        elif c.islower():
            # PCR_strainer marked this as a mismatch.  Check whether the oligo
            # has a degenerate base that legitimately covers the genome base.
            oligo_base  = oligo_seq[oligo_pos].upper()
            genome_base = c.upper()
            valid_bases = _IUPAC_MATCHES.get(oligo_base, {oligo_base})
            if genome_base not in valid_bases:
                errors[oligo_pos] = True   # genuine mismatch
            # if genome_base IS in valid_bases, it is a degenerate match —
            # leave errors[oligo_pos] as False
            oligo_pos += 1;  i += 1
        elif c.isupper():
            oligo_pos += 1;  i += 1      # confirmed match
        else:
            i += 1

    return errors


def compute_per_position_mm(
        pcr_results: pd.DataFrame
) -> Dict[str, Dict[str, dict]]:
    """Compute per-position mismatch prevalence (%) for every assay/oligo."""
    result: Dict[str, Dict[str, dict]] = {}
    if pcr_results.empty:
        return result

    for assay_name, assay_df in pcr_results.groupby('assay_name'):
        result[assay_name] = {}
        n_detected = len(assay_df)
        if n_detected == 0:
            continue

        for oligo in ('fwd_primer', 'rev_primer', 'probe'):
            seq_col  = f'{oligo}_seq'
            site_col = f'{oligo}_site_seq'
            name_col = f'{oligo}_name'

            if seq_col not in assay_df.columns:
                continue

            seqs = assay_df[seq_col].dropna().unique()
            if len(seqs) == 0 or not isinstance(seqs[0], str) or seqs[0] == '':
                continue

            oligo_seq = _sanitise_seq(seqs[0], f'{assay_name}/{oligo} seq')
            if not oligo_seq:
                continue

            raw_name  = (assay_df[name_col].dropna().iloc[0]
                         if name_col in assay_df.columns
                         and not assay_df[name_col].dropna().empty
                         else oligo)
            oligo_name = _sanitise_name(str(raw_name), f'{assay_name}/{oligo} name')

            n         = len(oligo_seq)
            mm_counts = np.zeros(n, dtype=np.int64)

            if site_col in assay_df.columns:
                for site_seq in assay_df[site_col].dropna():
                    for pos, is_err in enumerate(
                            parse_site_seq_positions(oligo_seq, str(site_seq))):
                        if pos < n and is_err:
                            mm_counts[pos] += 1

            result[assay_name][oligo] = {
                'seq':    oligo_seq,
                'name':   oligo_name,
                'mm_pct': [round(float(c) * 100.0 / n_detected, 1)
                           for c in mm_counts],
            }

    return result


# ── Status helpers ────────────────────────────────────────────────────────────

def _status(pct: float, pass_t: float, caution_t: float) -> str:
    if pct >= pass_t:    return 'pass'
    if pct >= caution_t: return 'caution'
    return 'action'


def _overall_status(assays: List[dict]) -> str:
    statuses = {a['status'] for a in assays}
    if 'action'  in statuses: return 'action'
    if 'caution' in statuses: return 'caution'
    return 'pass'


# ── Data assembly ─────────────────────────────────────────────────────────────

def build_report_data(
        dfs:    Dict[str, pd.DataFrame],
        config: ReportConfig,
) -> dict:
    """
    Assemble all report data into a plain Python dict.
    All values are Python primitives (str, int, float, list, dict, None).
    """
    ar  = dfs['assay_report']
    vr  = dfs['variant_report']
    ms  = dfs['missed_seqs']
    pcr = dfs['pcr_results']

    pass_t    = config.pass_threshold
    caution_t = config.caution_threshold

    total_genomes: Optional[int] = None
    if not ar.empty and 'total_targets' in ar.columns:
        vals = ar['total_targets'].dropna()
        if not vals.empty:
            total_genomes = int(vals.iloc[0])

    provenance = {
        'run_date':          datetime.now().strftime('%Y-%m-%d'),
        'tool_version':      _sanitise_name(config.tool_version, 'tool_version'),
        'total_genomes':     total_genomes,
        'assay_file':        os.path.basename(config.assay_file)  if config.assay_file  else None,
        'genome_file':       os.path.basename(config.genome_file) if config.genome_file else None,
        'min_tm':            float(config.min_tm)            if config.min_tm            is not None else None,
        'variant_threshold': float(config.variant_threshold) if config.variant_threshold is not None else None,
    }

    per_pos = compute_per_position_mm(pcr)

    assay_names: List[str] = []
    if not ar.empty and 'assay_name' in ar.columns:
        assay_names = ar['assay_name'].dropna().unique().tolist()
    elif not pcr.empty and 'assay_name' in pcr.columns:
        assay_names = pcr['assay_name'].dropna().unique().tolist()

    assays_data: List[dict] = []
    for assay_name in assay_names:
        assay_name_safe = _sanitise_name(str(assay_name), 'assay_name')

        ar_a = (ar[ar['assay_name'] == assay_name]
                if not ar.empty and 'assay_name' in ar.columns
                else pd.DataFrame())

        total_targets    = int(ar_a['total_targets'].iloc[0])    if not ar_a.empty and 'total_targets'    in ar_a.columns else 0
        detected_targets = int(ar_a['detected_targets'].iloc[0]) if not ar_a.empty and 'detected_targets' in ar_a.columns else 0
        perc_detected    = float(ar_a['perc_detected'].iloc[0])  if not ar_a.empty and 'perc_detected'    in ar_a.columns else 0.0

        err_dist: dict = {'labels': ['0 errors', '1 error', '2 errors', '3+ errors'],
                          'values': [0.0, 0.0, 0.0, 0.0]}
        if not ar_a.empty and 'total_errors' in ar_a.columns and 'perc_of_detected' in ar_a.columns:
            for _, row in ar_a.iterrows():
                errs = int(row['total_errors'])       if pd.notna(row.get('total_errors'))     else 0
                pct  = float(row['perc_of_detected']) if pd.notna(row.get('perc_of_detected')) else 0.0
                err_dist['values'][min(errs, 3)] += round(pct, 1)

        perfect_match_pct = err_dist['values'][0]
        missed_count      = total_targets - detected_targets

        missed_high_n_pct: Optional[float] = None
        if not ms.empty and 'assay_name' in ms.columns:
            ms_a = ms[ms['assay_name'] == assay_name]
            if not ms_a.empty and 'perc_Ns' in ms_a.columns and len(ms_a) > 0:
                missed_high_n_pct = round(
                    len(ms_a[ms_a['perc_Ns'] > 5]) * 100.0 / len(ms_a), 1)

        oligo_mm_any: Dict[str, float] = {}
        if not pcr.empty and 'assay_name' in pcr.columns:
            pcr_a = pcr[pcr['assay_name'] == assay_name]
            n_det = len(pcr_a)
            for oligo in ('fwd_primer', 'rev_primer', 'probe'):
                err_col = f'{oligo}_errors'
                if err_col in pcr_a.columns and n_det > 0:
                    n_err = int((pcr_a[err_col].fillna(0) > 0).sum())
                    oligo_mm_any[oligo] = round(n_err * 100.0 / n_det, 1)

        vr_a = (vr[vr['assay_name'] == assay_name]
                if not vr.empty and 'assay_name' in vr.columns
                else pd.DataFrame())

        oligo_labels = {'fwd_primer': 'Fwd primer',
                        'rev_primer': 'Rev primer',
                        'probe':      'Probe'}
        oligo_chart_labels: List[str]  = []
        oligo_chart_values: List[float] = []
        oligos_out: List[dict] = []

        for oligo in ('fwd_primer', 'rev_primer', 'probe'):
            pos_data = per_pos.get(assay_name, {}).get(oligo)
            if pos_data is None:
                continue

            any_mm_pct = oligo_mm_any.get(oligo, 0.0)
            oligo_chart_labels.append(oligo_labels[oligo])
            oligo_chart_values.append(any_mm_pct)

            variants_out: List[dict] = []
            if not vr_a.empty and 'oligo' in vr_a.columns:
                vr_ol = vr_a[vr_a['oligo'] == oligo].copy()
                if not vr_ol.empty and 'target_count' in vr_ol.columns:
                    vr_ol = vr_ol.sort_values('target_count', ascending=False)
                for _, vrow in vr_ol.iterrows():
                    variant_seq = _sanitise_seq(
                        str(vrow.get('oligo_site_variant', '')),
                        f'{assay_name}/{oligo} variant')
                    if not variant_seq:
                        continue
                    variants_out.append({
                        'seq':    variant_seq,
                        'errors': int(vrow['oligo_errors'])        if pd.notna(vrow.get('oligo_errors'))      else 0,
                        'count':  int(vrow['target_count'])        if pd.notna(vrow.get('target_count'))      else 0,
                        'pct':    round(float(vrow['perc_of_detected']), 1)
                                  if pd.notna(vrow.get('perc_of_detected')) else 0.0,
                    })

            oligos_out.append({
                'type':       oligo,
                'label':      oligo_labels[oligo],
                'name':       pos_data['name'],
                'seq':        pos_data['seq'],
                'any_mm_pct': any_mm_pct,
                'mm_pct':     pos_data['mm_pct'],
                'variants':   variants_out,
            })

        peak_oligo: Optional[str] = None
        peak_pct   = 0.0
        peak_pos: Optional[int]  = None
        for ol in oligos_out:
            if ol['mm_pct']:
                mx = max(ol['mm_pct'])
                if mx > peak_pct:
                    peak_pct   = mx
                    peak_oligo = ol['label']
                    peak_pos   = ol['mm_pct'].index(mx) + 1

        status = _status(perc_detected, pass_t, caution_t)
        if status == 'pass' and any(v >= OLIGO_ACTION for v in oligo_mm_any.values()):
            status = 'caution'

        action_msg: Optional[str] = None
        if status in ('caution', 'action'):
            if oligo_mm_any:
                worst  = max(oligo_mm_any, key=oligo_mm_any.get)   # type: ignore[arg-type]
                wpct   = oligo_mm_any[worst]
                wlabel = oligo_labels.get(worst, worst)
                action_msg = (
                    f'{wpct:.1f}% of detected genomes have at least one mismatch '
                    f'in the {wlabel.lower()}. Review the per-position heatmap and '
                    f'variant table to assess positional risk before continued clinical use.'
                )
            else:
                action_msg = (
                    f'Overall inclusivity is {perc_detected:.1f}%. '
                    f'Review the missed sequences report to determine whether this '
                    f'reflects genuine assay failure or poor reference genome quality.'
                )

        assays_data.append({
            'name':              assay_name_safe,
            'status':            status,
            'action_msg':        action_msg,
            'total_targets':     total_targets,
            'detected_targets':  detected_targets,
            'perc_detected':     round(perc_detected, 1),
            'perfect_match_pct': round(perfect_match_pct, 1),
            'missed_count':      missed_count,
            'missed_high_n_pct': missed_high_n_pct,
            'err_distribution':  err_dist,
            'oligo_chart':       {'labels': oligo_chart_labels,
                                  'values': oligo_chart_values},
            'oligos':            oligos_out,
            'peak_oligo':        peak_oligo,
            'peak_pct':          round(peak_pct, 1),
            'peak_pos':          peak_pos,
        })

    overall   = _overall_status(assays_data) if assays_data else 'pass'
    n_action  = sum(1 for a in assays_data if a['status'] == 'action')
    n_caution = sum(1 for a in assays_data if a['status'] == 'caution')
    n_pass    = sum(1 for a in assays_data if a['status'] == 'pass')
    n_total   = len(assays_data)
    aw        = 'assays' if n_total != 1 else 'assay'

    if overall == 'pass':
        overall_msg = f'All {n_total} {aw} show good inclusivity against the reference genome set.'
    elif overall == 'caution':
        overall_msg = (f'{n_caution} of {n_total} {aw} '
                       f'{"have" if n_caution != 1 else "has"} variants of potential concern. '
                       f'Review highlighted assays before continued use.')
    else:
        overall_msg = (f'{n_action} of {n_total} {aw} '
                       f'{"require" if n_action != 1 else "requires"} immediate review. '
                       f'Significant mismatch prevalence detected in circulating strains. '
                       f'Consult your laboratory director.')

    log.info('Report assembled: %d assays (%d pass / %d caution / %d action)',
             n_total, n_pass, n_caution, n_action)

    return {
        'provenance':      provenance,
        'overall_status':  overall,
        'overall_message': overall_msg,
        'n_pass':          n_pass,
        'n_caution':       n_caution,
        'n_action':        n_action,
        'assays':          assays_data,
        'thresholds':      {'pass': pass_t, 'caution': caution_t},
    }


# ── Chart generation (matplotlib → base64 PNG) ───────────────────────────────

def _png_b64(fig: 'plt.Figure') -> str:
    """Render a matplotlib figure to a base64-encoded PNG string, then close it."""
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=130, bbox_inches='tight')
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def _err_chart_png(err_dist: dict, status: str) -> str:
    """
    Vertical bar chart showing the percentage of detected genomes at each
    total-error level (0, 1, 2, 3+).  The zero-error bar uses the status
    colour; the remaining bars use a lighter tint.
    """
    labels = err_dist['labels']
    values = err_dist['values']
    colours = [_STATUS_COLOUR.get(status, '#888')] + \
              [_STATUS_COLOUR_LIGHT.get(status, '#ccc')] * (len(values) - 1)

    fig, ax = plt.subplots(figsize=(3.8, 2.2))
    bars = ax.bar(labels, values, color=colours, width=0.6,
                  edgecolor='white', linewidth=0.5)

    # Label each bar with its value
    for bar, val in zip(bars, values):
        if val > 0:
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.8,
                    f'{val:.1f}%',
                    ha='center', va='bottom', fontsize=7.5, color='#444')

    ax.set_ylim(0, min(100, max(values) * 1.25 + 5) if any(v > 0 for v in values) else 10)
    ax.set_ylabel('%', fontsize=8, color='#666')
    ax.tick_params(axis='x', labelsize=8, colors='#444')
    ax.tick_params(axis='y', labelsize=7, colors='#888')
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.0f}%'))
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_color('#ddd')
    ax.yaxis.grid(True, color='#eee', linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    fig.tight_layout(pad=0.4)
    return _png_b64(fig)


def _oligo_chart_png(oligo_chart: dict) -> str:
    """
    Horizontal bar chart showing the percentage of detected genomes with
    at least one mismatch in each oligo.  Each bar is coloured by severity.
    """
    labels = oligo_chart['labels']
    values = oligo_chart['values']

    if not labels:
        # No oligos found — return a tiny blank PNG
        fig, ax = plt.subplots(figsize=(3.8, 1.0))
        ax.axis('off')
        fig.patch.set_facecolor('white')
        return _png_b64(fig)

    colours = [
        _STATUS_COLOUR['action']  if v >= OLIGO_ACTION  else
        _STATUS_COLOUR['caution'] if v >= OLIGO_CAUTION else
        _STATUS_COLOUR['pass']
        for v in values
    ]

    fig, ax = plt.subplots(figsize=(3.8, 0.6 + 0.5 * len(labels)))
    bars = ax.barh(labels, values, color=colours, height=0.5,
                   edgecolor='white', linewidth=0.5)

    # Label each bar
    for bar, val in zip(bars, values):
        ax.text(bar.get_width() + 0.3,
                bar.get_y() + bar.get_height() / 2,
                f'{val:.1f}%',
                ha='left', va='center', fontsize=8, color='#444')

    x_max = max(35, max(values) * 1.3 + 5) if values else 35
    ax.set_xlim(0, x_max)
    ax.set_xlabel('%', fontsize=8, color='#666')
    ax.tick_params(axis='y', labelsize=8.5, colors='#444')
    ax.tick_params(axis='x', labelsize=7,   colors='#888')
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.0f}%'))
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_color('#ddd')
    ax.xaxis.grid(True, color='#eee', linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    fig.tight_layout(pad=0.4)
    return _png_b64(fig)


# ── HTML generation helpers ───────────────────────────────────────────────────

def _esc(s: object) -> str:
    """HTML-escape any value before inserting into an HTML attribute or element."""
    return (str(s)
            .replace('&',  '&amp;')
            .replace('<',  '&lt;')
            .replace('>',  '&gt;')
            .replace('"',  '&quot;')
            .replace("'",  '&#39;'))


def _fmt_pct(v: Optional[float]) -> str:
    return '&mdash;' if v is None else f'{v:.1f}%'


def _fmt_n(v: Optional[int]) -> str:
    return '&mdash;' if v is None else f'{v:,}'


def _mm_class(pct: float) -> str:
    if pct == 0:  return 'mm0'
    if pct < 3:   return 'mml'
    if pct < 10:  return 'mmm'
    if pct < 20:  return 'mmh'
    return 'mmc'


def _badge(status: str) -> str:
    labels = {'pass': 'PASS', 'caution': 'CAUTION', 'action': 'ACTION REQUIRED'}
    return (f'<span class="badge b-{_esc(status)}">'
            f'&bull; {_esc(labels.get(status, status))}</span>')


def _badge_sm(status: str) -> str:
    """Smaller badge for the summary table."""
    labels = {'pass': 'PASS', 'caution': 'CAUTION', 'action': 'ACTION REQ.'}
    return (f'<span class="badge b-{_esc(status)}" style="font-size:11px">'
            f'&bull; {_esc(labels.get(status, status))}</span>')


def _oligo_tag(oligo_type: str) -> str:
    cls    = {'fwd_primer': 't-fwd', 'rev_primer': 't-rev', 'probe': 't-probe'}.get(oligo_type, '')
    labels = {'fwd_primer': 'FWD',   'rev_primer': 'REV',   'probe': 'PROBE'}
    return f'<span class="otag {cls}">{labels.get(oligo_type, _esc(oligo_type))}</span>'


def _heatmap_html(seq: str, mm_pct: List[float]) -> str:
    """
    Generate a row of colour-coded base tiles.
    Uses HTML title= attributes for hover tooltips — no JavaScript needed.
    """
    tiles = ['<div class="hrow">']
    for i, (base, pct) in enumerate(zip(seq, mm_pct)):
        cls      = _mm_class(pct)
        pct_str  = f'{pct:.1f}%' if pct > 0 else ''
        is3prime = i >= len(seq) - 5
        tip      = f'Position {i+1}: {pct:.1f}% mismatch'
        if pct > 1 and is3prime:
            tip += ' — near 3\u2032 end, higher impact'
        tiles.append(
            f'<div class="bb {cls}" title="{_esc(tip)}">'
            f'<span class="bl">{_esc(base)}</span>'
            f'<span class="bp">{_esc(pct_str)}</span>'
            f'<span class="pl">{i + 1}</span>'
            f'</div>'
        )
    tiles.append('</div>')
    return ''.join(tiles)


def _variant_seq_html(seq: str) -> str:
    """
    Render PCR_strainer site_seq notation as HTML with mismatched bases
    highlighted in red.  Input is pre-validated by _sanitise_seq so only
    IUPAC chars, '-', '(', ')' are present.
    """
    out = []
    i   = 0
    while i < len(seq):
        c = seq[i]
        if c == '(':
            j   = seq.find(')', i)
            ins = seq[i+1:j] if j != -1 else seq[i+1:]
            out.append(f'<span class="mm">({_esc(ins)})</span>')
            i = j + 1 if j != -1 else len(seq)
        elif c == '-':
            out.append('<span class="mm">-</span>')
            i += 1
        elif c.islower():
            out.append(f'<span class="mm">{_esc(c.upper())}</span>')
            i += 1
        else:
            out.append(_esc(c))
            i += 1
    return ''.join(out)


def _interpret_variant(oligo_seq: str, var_seq: str,
                       pct: float, oligo_type: str) -> str:
    """Plain-text interpretation of a variant for the variants table."""
    mm_pos = []
    oligo_pos = 0
    for c in var_seq:
        if c == '(':
            continue
        if c == '-' or c.islower():
            mm_pos.append(oligo_pos)
        if c != '(':
            oligo_pos += 1

    n        = len(oligo_seq)
    is3prime = any(p >= n - 5 for p in mm_pos)
    is5prime = all(p < 5      for p in mm_pos) if mm_pos else False
    is_probe = oligo_type == 'probe'

    risk = ''
    if pct >= 15 or (pct >= 5 and is3prime):
        risk = 'HIGH RISK. '
    elif pct >= 5:
        risk = 'MODERATE RISK. '

    pos_str = ', '.join(f'pos {p+1}' for p in mm_pos) if mm_pos else ''
    detail  = f'Mismatch at {pos_str}.' if pos_str else ''
    if is3prime and not is_probe:
        detail += ' Near 3\u2032 end — can inhibit Taq extension.'
    elif is5prime and not is_probe:
        detail += ' 5\u2032 end — usually tolerated.'
    if is_probe:
        detail += ' Probe mismatch may reduce fluorescence signal.'

    threshold = ('Urgent review.'      if pct >= 20 else
                 'Review recommended.' if pct >= 10 else
                 'Monitor for increase.' if pct >= 5 else 'Monitor.')

    return f'{risk}{detail} {threshold}'.strip()


# ── HTML sections ─────────────────────────────────────────────────────────────

_CSS = """
*,*::before,*::after{box-sizing:border-box;margin:0;padding:0}
:root{
  --blue:#1F4E79;--blue-m:#2E75B6;--blue-l:#D6E4F0;--blue-p:#EBF3FB;
  --green:#1A5C2E;--green-l:#D4EDDA;
  --amber:#7D4E00;--amber-l:#FFF2CC;
  --red:#7B0000;--red-l:#FCE4D6;
  --grey-l:#F5F5F5;--grey-b:#E0E0E0;
  --text:#1A1A1A;--text-m:#555;--text-l:#888;
  --r:10px;
  font-family:system-ui,-apple-system,'Segoe UI',Arial,sans-serif;
}
body{background:#F0F4F8;color:var(--text);font-size:14px;line-height:1.5}
.page{max-width:1080px;margin:0 auto;padding:24px 20px 60px}
/* Header */
.rh{background:var(--blue);color:#fff;border-radius:var(--r);padding:24px 28px;
    margin-bottom:20px;display:flex;justify-content:space-between;
    align-items:flex-start;gap:20px}
.rh h1{font-size:22px;font-weight:600;letter-spacing:-.3px}
.rh .sub{font-size:12px;opacity:.75;margin-top:4px}
.prov{display:grid;grid-template-columns:auto auto;gap:4px 20px;font-size:12px;
      opacity:.9;text-align:right;white-space:nowrap}
.prov .lbl{opacity:.65}.prov .val{font-weight:500}
/* Verdict banners */
.vb{border-radius:var(--r);padding:14px 20px;display:flex;align-items:center;
    gap:14px;margin-bottom:20px;font-size:13px}
.vb.pass   {background:var(--green-l);border:1.5px solid #2E7D46;color:var(--green)}
.vb.caution{background:var(--amber-l);border:1.5px solid #B07800;color:var(--amber)}
.vb.action {background:var(--red-l);  border:1.5px solid #C0392B;color:var(--red)}
.vi{font-size:22px;flex-shrink:0}.vt{font-weight:600;font-size:15px}
/* Summary table */
.sm{background:#fff;border-radius:var(--r);box-shadow:0 2px 8px rgba(0,0,0,.07);
    padding:20px 24px;margin-bottom:20px}
.sec{font-size:12px;font-weight:600;color:var(--text-m);text-transform:uppercase;
     letter-spacing:.6px;margin-bottom:14px}
table.mt{width:100%;border-collapse:collapse;font-size:13px}
table.mt th{text-align:left;padding:8px 12px;background:var(--grey-l);
            font-weight:600;font-size:12px;color:var(--text-m);
            border-bottom:1px solid var(--grey-b)}
table.mt td{padding:10px 12px;border-bottom:1px solid var(--grey-b);
            vertical-align:middle}
table.mt tr:last-child td{border-bottom:none}
table.mt a{color:inherit;text-decoration:none}
table.mt a:hover{text-decoration:underline}
/* Badges */
.badge{display:inline-flex;align-items:center;gap:4px;padding:3px 10px;
       border-radius:20px;font-size:11px;font-weight:700}
.b-pass   {background:var(--green-l);color:var(--green)}
.b-caution{background:var(--amber-l);color:var(--amber)}
.b-action {background:var(--red-l);  color:var(--red)}
/* Inclusivity bar */
.ib-wrap{display:flex;align-items:center;gap:10px}
.ib{height:8px;border-radius:4px;background:var(--grey-b);flex:1;overflow:hidden}
.ib-fill{height:100%;border-radius:4px}
.f-pass{background:#2E7D46}.f-caution{background:#B07800}.f-action{background:#C0392B}
.ipct{font-size:12px;font-weight:600;min-width:38px;text-align:right}
.c-pass{color:#2E7D46}.c-caution{color:#B07800}.c-action{color:#C0392B}
.c-neu{color:var(--text)}
/* Collapsible assay cards — <details>/<summary> */
details.ac{background:#fff;border-radius:var(--r);
           box-shadow:0 2px 8px rgba(0,0,0,.07);
           margin-bottom:18px;overflow:hidden;border-top:4px solid var(--grey-b)}
details.ac.pass   {border-top-color:#2E7D46}
details.ac.caution{border-top-color:#B07800}
details.ac.action {border-top-color:#C0392B}
details.ac > summary{
  display:flex;align-items:center;justify-content:space-between;
  padding:16px 22px;cursor:pointer;user-select:none;list-style:none}
details.ac > summary::-webkit-details-marker{display:none}
details.ac > summary:hover{background:var(--grey-l)}
.ct-row{display:flex;align-items:center;gap:12px}
.ct{font-size:16px;font-weight:600}.csub{font-size:12px;color:var(--text-l);margin-top:2px}
.chev{color:var(--text-l);font-size:20px;line-height:1;transition:transform .2s}
details.ac[open] > summary .chev{transform:rotate(90deg)}
.card-body{padding:0 22px 22px}
/* Metrics row */
.mr{display:flex;gap:14px;margin:16px 0;flex-wrap:wrap}
.met{background:var(--grey-l);border-radius:8px;padding:12px 16px;flex:1;min-width:120px}
.ml{font-size:11px;color:var(--text-l);font-weight:600;text-transform:uppercase;
    letter-spacing:.5px}
.mv{font-size:22px;font-weight:600;margin-top:2px}
.ms2{font-size:11px;color:var(--text-l);margin-top:2px}
/* Charts */
.two{display:grid;grid-template-columns:1fr 1fr;gap:20px;margin-top:20px}
@media(max-width:680px){.two{grid-template-columns:1fr}}
.pt{font-size:12px;font-weight:600;color:var(--text-m);text-transform:uppercase;
    letter-spacing:.5px;margin-bottom:10px}
.chart-img{width:100%;height:auto;display:block}
/* Heatmap */
.os{margin-top:20px}.ob{margin-bottom:20px}
.oh{display:flex;align-items:center;gap:10px;margin-bottom:6px}
.otag{font-size:10px;font-weight:700;text-transform:uppercase;letter-spacing:.5px;
      padding:2px 8px;border-radius:12px}
.t-fwd  {background:#DBEAFE;color:#1E40AF}
.t-rev  {background:#EDE9FE;color:#5B21B6}
.t-probe{background:#FEF9C3;color:#92400E}
.on{font-size:12px;font-weight:600;color:var(--text-m)}
.oseq{font-family:'Courier New',monospace;font-size:11px;color:var(--text-l)}
.ends{display:flex;margin-bottom:3px;font-size:10px;color:var(--text-l);font-weight:500}
.hrow{display:flex;gap:2px;align-items:flex-end;flex-wrap:nowrap;
      overflow-x:auto;padding-bottom:28px;position:relative}
.bb{display:flex;flex-direction:column;align-items:center;width:28px;min-width:28px;
    border-radius:4px;position:relative;cursor:default}
.bb:hover{filter:brightness(.92)}
.bl{font-size:11px;font-weight:700;font-family:'Courier New',monospace;
    line-height:1;padding:5px 0 4px}
.bp{font-size:9px;line-height:1;padding:0 0 4px}
.pl{font-size:9px;color:var(--text-l);margin-top:3px;position:absolute;bottom:-22px}
.mm0{background:#EAF5ED;color:#1A5C2E}
.mml{background:#FFF3CC;color:#7D4E00}
.mmm{background:#FDDCB0;color:#7D3000}
.mmh{background:#FBBDBD;color:#7B0000}
.mmc{background:#E53E3E;color:#fff}
.leg{display:flex;align-items:center;gap:6px;margin-top:24px;font-size:11px;
     color:var(--text-l);flex-wrap:wrap}
.ls{width:14px;height:14px;border-radius:3px;flex-shrink:0}
/* Variants table */
table.vt{width:100%;border-collapse:collapse;font-size:12px;margin-top:12px}
table.vt th{text-align:left;padding:6px 10px;background:var(--grey-l);
            font-weight:600;font-size:11px;color:var(--text-m);
            border-bottom:1px solid var(--grey-b)}
table.vt td{padding:7px 10px;border-bottom:1px solid var(--grey-b);vertical-align:top}
table.vt tr:last-child td{border-bottom:none}
.vs{font-family:'Courier New',monospace;font-size:12px;letter-spacing:1px;
    word-break:break-all}
.vs .mm{color:#C0392B;font-weight:700}
.vb2{height:6px;border-radius:3px;background:var(--grey-b);margin-top:4px;overflow:hidden}
.vbf{height:100%;border-radius:3px}
/* Footer */
.rf{text-align:center;font-size:11px;color:var(--text-l);margin-top:30px;
    padding-top:16px;border-top:1px solid var(--grey-b)}
"""


def _html_head(title: str = 'PCR_strainer \u00b7 Assay Inclusivity Report') -> str:
    return (
        '<!DOCTYPE html>\n<html lang="en">\n<head>\n'
        '<meta charset="UTF-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1.0">\n'
        '<meta http-equiv="Content-Security-Policy" '
        'content="default-src \'none\'; style-src \'unsafe-inline\'; img-src data:;">\n'
        f'<title>{_esc(title)}</title>\n'
        f'<style>\n{_CSS}\n</style>\n'
        '</head>\n<body>\n<div class="page">\n'
    )


def _html_header(pv: dict) -> str:
    rows = [
        ('Run date',          pv.get('run_date')),
        ('Tool',              pv.get('tool_version')),
        ('Reference genomes', f"{pv['total_genomes']:,} sequences" if pv.get('total_genomes') else None),
        ('Assay file',        pv.get('assay_file')),
        ('Genome file',       pv.get('genome_file')),
        ('Min Tm',            f"{pv['min_tm']} \u00b0C" if pv.get('min_tm') is not None else None),
        ('Variant threshold', f"\u2265{pv['variant_threshold']}%" if pv.get('variant_threshold') is not None else None),
    ]
    prov_cells = ''.join(
        f'<span class="lbl">{_esc(label)}</span>'
        f'<span class="val">{_esc(value)}</span>'
        for label, value in rows if value is not None
    )
    return (
        '<div class="rh"><div>'
        '<h1>PCR_strainer &middot; Assay Inclusivity Report</h1>'
        '<div class="sub">Automated summary &mdash; for laboratory review only. '
        'Not for patient use.</div>'
        f'</div><div class="prov">{prov_cells}</div></div>\n'
    )


def _html_verdict(overall_status: str, overall_message: str) -> str:
    icons   = {'pass': '&#10003;', 'caution': '&#9888;', 'action': '&#9888;'}
    titles  = {
        'pass':    'All assays show acceptable inclusivity',
        'caution': 'One or more assays require review',
        'action':  'Action required \u2014 significant mismatch prevalence detected',
    }
    return (
        f'<div class="vb {_esc(overall_status)}">'
        f'<span class="vi">{icons.get(overall_status, "")}</span>'
        f'<div><div class="vt">{_esc(titles.get(overall_status, ""))}</div>'
        f'<div style="font-size:13px;margin-top:3px">{_esc(overall_message)}</div>'
        f'</div></div>\n'
    )


def _html_summary_table(assays: List[dict]) -> str:
    rows = []
    for idx, a in enumerate(assays):
        vals      = a['oligo_chart']['values']
        labels    = a['oligo_chart']['labels']
        worst_idx = vals.index(max(vals)) if vals else 0
        worst     = f"{_esc(labels[worst_idx])} ({_fmt_pct(vals[worst_idx])})" if vals else '&mdash;'
        peak      = (f"{_esc(a['peak_oligo'])} pos {_esc(str(a['peak_pos'] or '?'))} "
                     f"&middot; {_fmt_pct(a['peak_pct'])}"
                     if a['peak_oligo'] else '&mdash;')
        pct       = a['perc_detected']
        rows.append(
            f'<tr>'
            f'<td><a href="#assay-{idx}"><strong>{_esc(a["name"])}</strong></a></td>'
            f'<td>{_badge_sm(a["status"])}</td>'
            f'<td>{_fmt_n(a["detected_targets"])} / {_fmt_n(a["total_targets"])}</td>'
            f'<td><div class="ib-wrap">'
            f'<div class="ib"><div class="ib-fill f-{_esc(a["status"])}" '
            f'style="width:{min(pct, 100):.1f}%"></div></div>'
            f'<span class="ipct c-{_esc(a["status"])}">{_fmt_pct(pct)}</span>'
            f'</div></td>'
            f'<td style="font-size:12px">{worst}</td>'
            f'<td style="font-family:monospace;font-size:11px">{peak}</td>'
            f'</tr>\n'
        )
    return (
        '<div class="sm"><div class="sec">Assay summary</div>'
        '<table class="mt"><thead><tr>'
        '<th>Assay</th><th>Status</th><th>Detected</th>'
        '<th>Inclusivity</th><th>Oligos with errors</th><th>Peak mismatch</th>'
        '</tr></thead><tbody>'
        + ''.join(rows) +
        '</tbody></table></div>\n'
    )


def _html_assay_card(a: dict, thresholds: dict, idx: int) -> str:
    """Generate one collapsible <details> assay card with static charts and heatmaps."""
    status    = a['status']
    is_open   = status in ('action', 'caution') or idx == 0
    open_attr = ' open' if is_open else ''

    # ── summary row (the clickable header) ──────────────────────────────────
    summary = (
        f'<summary>'
        f'<div class="ct-row">{_badge(status)}'
        f'<div><div class="ct">{_esc(a["name"])}</div>'
        f'<div class="csub">{_fmt_n(a["detected_targets"])} detected '
        f'of {_fmt_n(a["total_targets"])} genomes</div></div></div>'
        f'<span class="chev">&#8250;</span>'
        f'</summary>\n'
    )

    parts = [f'<details class="ac {_esc(status)}" id="assay-{idx}"{open_attr}>\n',
             summary,
             '<div class="card-body">\n']

    # ── action/caution banner ────────────────────────────────────────────────
    if a['action_msg']:
        icons = {'pass': '&#10003;', 'caution': '&#9888;', 'action': '&#9888;'}
        parts.append(
            f'<div class="vb {_esc(status)}" style="margin:0 0 16px">'
            f'<span class="vi" style="font-size:18px">{icons.get(status, "")}</span>'
            f'<div style="font-size:13px">{_esc(a["action_msg"])}</div></div>\n'
        )

    # ── metric cards ────────────────────────────────────────────────────────
    missed_pct  = (a['missed_count'] / max(a['total_targets'], 1)) * 100
    missed_note = (f' &middot; {_esc(str(a["missed_high_n_pct"]))}% high-N'
                   if a['missed_high_n_pct'] is not None else '')
    pf_class    = ('c-pass'    if a['perc_detected'] >= thresholds['pass']
                   else 'c-caution' if a['perc_detected'] >= thresholds['caution']
                   else 'c-action')
    parts.append(
        f'<div class="mr">'
        f'<div class="met"><div class="ml">Inclusivity</div>'
        f'<div class="mv c-{_esc(status)}">{_fmt_pct(a["perc_detected"])}</div>'
        f'<div class="ms2">{_fmt_n(a["detected_targets"])} / {_fmt_n(a["total_targets"])} genomes</div></div>'
        f'<div class="met"><div class="ml">Perfect match</div>'
        f'<div class="mv {pf_class}">{_fmt_pct(a["perfect_match_pct"])}</div>'
        f'<div class="ms2">0 errors across all oligos</div></div>'
        f'<div class="met"><div class="ml">Missed genomes</div>'
        f'<div class="mv c-neu">{_fmt_n(a["missed_count"])}</div>'
        f'<div class="ms2">{_esc(f"{missed_pct:.1f}")}% of total{missed_note}</div></div>'
        f'</div>\n'
    )

    # ── charts (base64 PNG) ──────────────────────────────────────────────────
    err_png   = _err_chart_png(a['err_distribution'], status)
    oligo_png = _oligo_chart_png(a['oligo_chart'])
    parts.append(
        f'<div class="two">'
        f'<div><div class="pt">Error distribution</div>'
        f'<img class="chart-img" src="data:image/png;base64,{err_png}" '
        f'alt="Error distribution chart"></div>'
        f'<div><div class="pt">Mismatch rate by oligo</div>'
        f'<img class="chart-img" src="data:image/png;base64,{oligo_png}" '
        f'alt="Mismatch rate by oligo chart"></div>'
        f'</div>\n'
    )

    # ── per-position heatmaps ────────────────────────────────────────────────
    parts.append('<div class="os"><div class="pt" style="margin-top:20px">'
                 'Per-position mismatch frequency</div>\n')

    for ol in a['oligos']:
        parts.append(
            f'<div class="ob">'
            f'<div class="oh">'
            f'{_oligo_tag(ol["type"])}'
            f'<span class="on">{_esc(ol["name"])}</span>'
            f'<span class="oseq">&nbsp;&middot;&nbsp;{_esc(ol["seq"])}</span>'
            f'</div>'
            f'<div class="ends"><span>5\'</span>'
            f'<span style="flex:1"></span><span>3\'</span></div>'
            + _heatmap_html(ol['seq'], ol['mm_pct']) +
            '</div>\n'
        )

    parts.append(
        '<div class="leg">'
        '<div class="ls mm0"></div><span>0%</span>'
        '<div class="ls mml" style="margin-left:6px"></div><span>&lt;3%</span>'
        '<div class="ls mmm" style="margin-left:6px"></div><span>3\u201310%</span>'
        '<div class="ls mmh" style="margin-left:6px"></div><span>10\u201320%</span>'
        '<div class="ls mmc" style="margin-left:6px"></div><span>&gt;20%</span>'
        '<span style="margin-left:8px">mismatch prevalence at position '
        '&mdash; hover tile for detail</span>'
        '</div></div>\n'
    )

    # ── variants table ───────────────────────────────────────────────────────
    parts.append(
        '<div style="margin-top:20px">'
        '<div class="pt">Sequence variants at oligo binding sites</div>'
        '<table class="vt"><thead><tr>'
        '<th>Oligo</th><th>Variant sequence</th>'
        '<th>Errors</th><th>Prevalence</th><th>Interpretation</th>'
        '</tr></thead><tbody>\n'
    )

    has_variants = False
    for ol in a['oligos']:
        for v in ol['variants']:
            has_variants = True
            var_col = ('#C0392B' if v['pct'] >= 15 else
                       '#B07800' if v['pct'] >= 5  else '#555')
            row_bg  = ' style="background:#FEF2F0"' if v['pct'] >= 15 else ''
            fw      = '600' if v['pct'] >= 5 else '400'
            interp  = _interpret_variant(ol['seq'], v['seq'], v['pct'], ol['type'])
            bar_w   = min(v['pct'] * 3, 100)
            parts.append(
                f'<tr{row_bg}>'
                f'<td>{_oligo_tag(ol["type"])}</td>'
                f'<td><div class="vs">{_variant_seq_html(v["seq"])}</div></td>'
                f'<td>{_esc(str(v["errors"]))}</td>'
                f'<td><span style="color:{var_col};font-weight:{fw}">'
                f'{_fmt_pct(v["pct"])} ({_fmt_n(v["count"])} genomes)</span>'
                f'<div class="vb2"><div class="vbf" '
                f'style="width:{bar_w:.1f}%;background:{var_col}"></div></div></td>'
                f'<td style="font-size:11px;color:var(--text-m)">{_esc(interp)}</td>'
                f'</tr>\n'
            )

    if not has_variants:
        parts.append(
            '<tr><td colspan="5" style="text-align:center;'
            'color:var(--text-l);padding:14px">'
            'No variants above reporting threshold</td></tr>\n'
        )

    parts.append('</tbody></table></div>\n')
    parts.append('</div>\n')   # card-body
    parts.append('</details>\n')
    return ''.join(parts)


def _html_footer(pv: dict) -> str:
    return (
        f'<div class="rf">'
        f'Generated by {_esc(pv.get("tool_version", "PCR_strainer"))} '
        f'&middot; {_esc(pv.get("run_date", ""))}<br>'
        f'This report summarises computational analysis only. '
        f'Wet-lab validation is required before clinical decisions.<br>'
        f'Results are based on publicly available reference genomes and may not '
        f'represent all strains circulating in your region.'
        f'</div>\n'
        f'</div>\n'   # .page
        f'</body>\n</html>\n'
    )


# ── Top-level render ──────────────────────────────────────────────────────────

def render_html(report_data: dict) -> str:
    """
    Generate a complete, self-contained, JavaScript-free HTML report.
    Charts are rendered by matplotlib and embedded as base64 PNG images.
    """
    pv         = report_data['provenance']
    thresholds = report_data['thresholds']

    parts = [
        _html_head(),
        _html_header(pv),
        _html_verdict(report_data['overall_status'], report_data['overall_message']),
        _html_summary_table(report_data['assays']),
    ]
    for idx, assay in enumerate(report_data['assays']):
        parts.append(_html_assay_card(assay, thresholds, idx))
    parts.append(_html_footer(pv))

    return ''.join(parts)


# ── File writing ──────────────────────────────────────────────────────────────

def _write_file(path: str, content: str, mode: int) -> None:
    """
    Write content to path with explicit permissions set atomically at creation.
    Uses os.open() rather than open() + chmod() to avoid a race window.
    """
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC | os.O_NOFOLLOW, mode)
    try:
        with os.fdopen(fd, 'w', encoding='utf-8') as fh:
            fh.write(content)
    except Exception:
        try:
            os.unlink(path)
        except OSError:
            pass
        raise
    log.debug('Wrote %s (mode %s, %d chars)', path, oct(mode), len(content))


# ── Public API ────────────────────────────────────────────────────────────────

def write_html_report(
        name:               str,
        out_path:           str,
        genome_file:        Optional[str]   = None,
        assay_file:         Optional[str]   = None,
        min_tm:             Optional[float] = None,
        variant_threshold:  Optional[float] = None,
        tool_version:       Optional[str]   = None,
        pass_threshold:     float           = INCLUSIVITY_PASS,
        caution_threshold:  float           = INCLUSIVITY_CAUTION,
        file_mode:          int             = _OUTPUT_FILE_MODE,
) -> None:
    """
    Generate and write the HTML report.

    Call this AFTER PCR_strainer's write_assay_report(), write_variant_report(),
    write_missed_seqs_report(), and write_tntblast_results() have run — this
    function reads the TSV files those functions produced.

    Parameters
    ----------
    name              : output name prefix (same value passed to PCR_strainer -o)
    out_path          : output directory
    genome_file       : path to reference FASTA (basename shown in provenance)
    assay_file        : path to assay CSV (basename shown in provenance)
    min_tm            : minimum Tm threshold used (shown in provenance)
    variant_threshold : variant prevalence threshold used (shown in provenance)
    tool_version      : PCR_strainer version string
    pass_threshold    : inclusivity %% >= this -> PASS
    caution_threshold : inclusivity %% >= this -> CAUTION (else ACTION)
    file_mode         : octal permission mode for the output file (default 0o640)

    Integration example in PCR_strainer main():
        write_assay_report(name, out_path, ...)
        write_variant_report(name, out_path, ...)
        write_missed_seqs_report(name, out_path, ...)
        write_tntblast_results(name, out_path, ...)
        from pcr_strainer_report import write_html_report
        write_html_report(
            name, out_path,
            genome_file=args['-g'], assay_file=args['-a'],
            min_tm=args['-m'], variant_threshold=args['-t'],
            tool_version=version,
        )
    """
    prefix = os.path.join(out_path, name)
    dfs    = load_outputs(prefix)

    config = ReportConfig(
        genome_file=genome_file,
        assay_file=assay_file,
        min_tm=min_tm,
        variant_threshold=variant_threshold,
        tool_version=tool_version or 'PCR_strainer',
        pass_threshold=pass_threshold,
        caution_threshold=caution_threshold,
    )

    report_data = build_report_data(dfs, config)
    html        = render_html(report_data)
    out_file    = os.path.join(out_path, f'{name}_report.html')

    _write_file(out_file, html, file_mode)
    log.info('HTML report written to %s', out_file)


# ── Standalone entry point ────────────────────────────────────────────────────

def main(argv: Optional[List[str]] = None) -> None:
    args = parse_args(argv)
    _configure_logging(args.log_file)

    if not args.output_prefix:
        log.error('--output-prefix is required.')
        sys.exit(1)

    try:
        resolved = _validate_output_prefix(args.output_prefix)
    except ValueError as exc:
        log.error('Invalid output prefix: %s', exc)
        sys.exit(1)

    log.info('PCR_strainer Report Generator')
    log.info('Output will be written to: %s_report.html', resolved)

    dfs = load_outputs(args.output_prefix)
    if all(df.empty for df in dfs.values()):
        log.error('No PCR_strainer output files found for prefix: %s  '
                  'Check --output-prefix.', args.output_prefix)
        sys.exit(1)

    config = ReportConfig(
        genome_file=args.genome_file,
        assay_file=args.assay_file,
        min_tm=args.min_tm,
        variant_threshold=args.variant_threshold,
        tool_version=args.tool_version or 'PCR_strainer',
        pass_threshold=args.pass_threshold,
        caution_threshold=args.caution_threshold,
    )

    report_data = build_report_data(dfs, config)
    html        = render_html(report_data)
    out_file    = f'{args.output_prefix}_report.html'

    _write_file(out_file, html, args.file_mode)
    log.info('Report written to: %s', out_file)
    log.info('Open in any web browser to view.')


if __name__ == '__main__':
    main()
