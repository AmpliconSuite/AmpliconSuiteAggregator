"""
asa_aggregator.py
Core implementation of the AmpliconSuiteAggregator pipeline.

Pipeline stages (implemented in separate segments):
  Stage 1  — CLI + scaffolding          (this file: constants, helpers, Aggregator skeleton)
  Stage 2  — Extraction engine
  Stage 3  — Discovery engine
  Stage 4  — result_table parser + sample registry
  Stage 5  — Output tree builder
  Stage 6  — run.json constructor
"""

from __future__ import annotations  # enables built-in generic hints (list[], dict[], etc.) on Python 3.8/3.9

import json
import os
import re
import shutil
import sys
import tarfile
import zipfile
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

__version__ = "7.2.0"

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Files/extensions we never want to include in output archives
EXCLUSION_SUFFIXES: Tuple[str, ...] = (
    ".bam", ".fq", ".fastq", ".cram",
    ".fq.gz", ".fastq.gz",
    "aln_stage.stderr", ".fai",
    ".call.cns.gz", ".cnr.gz",
)

# Legacy CNV BED suffix (pre-dates _CNV_CALLS.bed naming).
# Covers both [sname]_uncorr_CN.bed and [sname].cs.rmdup_uncorr_CN.bed —
# the .cs.rmdup prefix is handled by the same dot-split logic as _CNV_CALLS.bed.
LEGACY_CNV_BED_SUFFIX: str = "_uncorr_CN.bed"

# Whitelist of file suffixes kept when copying a sample's AA results directory.
# Anything not matching one of these suffixes is dropped (hidden files are always dropped).
AA_DIR_INCLUDE_SUFFIXES: Tuple[str, ...] = (
    "_graph.txt", "_cycles.txt", ".png", ".pdf", "_summary.txt", ".log",
)

# Recognised archive extensions (in priority order for nested extraction)
ARCHIVE_EXTENSIONS: Tuple[str, ...] = (".tar.gz", ".tar", ".zip")

# Suffixes that identify an AmpliconSuite-related directory
AMPLICON_DIR_SUFFIXES: Tuple[str, ...] = (
    "_AA_results",
    "_cnvkit_output",
    "_cnvkit_outputs",
    "_classification",
)

# Files that identify a directory as containing AC classification output.
# A classification dir must have BOTH a profiles file AND a result_table file
# (warn if profiles found without result_table).
AC_PROFILES_SUFFIX = "_amplicon_classification_profiles.tsv"
AC_RESULT_TABLE_SUFFIX = "_result_table.tsv"

# Miscellaneous per-sample AA files (found as siblings to _AA_results/)
MISC_AA_SUFFIXES: Tuple[str, ...] = (
    "_AA_CNV_SEEDS.bed",
    "_CNV_SEEDS.bed",       # CoRAL's naming (no "AA_" infix)
    "_CNV_CALLS_unfiltered_gains.bed",
    "_finish_flag.txt",
    "_run_metadata.json",
    "_sample_metadata.json",
    "_timing_log.txt",
    # NOTE: top-level [sample_name].log is handled separately because the
    # suffix alone (".log") is too generic — we match by exact stem.
)

# Content-sniffing signals used to tell CoRAL output from AmpliconArchitect
# output, since AC emits no tool-provenance column of its own. See
# is_aa_summary_content()/is_coral_summary_content() below and
# Aggregator._detect_reconstruction_tool() in asa_stages.py.
CORAL_SUMMARY_RE = re.compile(r"\d+/\d+ amplicons solved")
CORAL_BANNER_PREFIX = "CoRAL"          # e.g. "CoRAL v2.2.0" summary banner line
CORAL_HEADER_PREFIX = "#CoRAL"         # future graph.txt/cycles.txt header line

# Tool version banner lines, scanned out of AA/AmpliconSuite-pipeline/AC logs.
# "Amp\w*Suite" tolerates AmpliconSuite-pipeline's own logging typo
# ("AmpiconSuite-pipeline version ...", missing the "l") alongside the
# correctly-spelled form, in case it's ever fixed upstream.
AA_VERSION_RE = re.compile(r"AmpliconArchitect version\s+(\S+)")
AS_PIPELINE_VERSION_RE = re.compile(r"Amp\w*Suite-pipeline version\s+(\S+)")
# AC's own banner has no "version" keyword, e.g. "AmpliconClassifier 2.0.0".
AC_VERSION_RE = re.compile(r"AmpliconClassifier\s+(\S+)")

# Per-amplicon file kind -> canonical filename suffix. This is the naming
# convention every per-amplicon file is normalized to on output (Stage 5
# copy, Stage 6 resolve) — one source of truth for both.
AMPLICON_FILE_EXT_MAP: Dict[str, str] = {
    "cycles_png": "_cycles.png",
    "cycles_pdf": "_cycles.pdf",
    "cycles":     "_cycles.txt",
    "graph":      "_graph.txt",
    "pdf":        ".pdf",
    "png":        ".png",
}

# (suffix_to_match, key) pairs used to *recognise* a per-amplicon source
# file at discovery time. Ordered longest/most-specific-suffix-first so
# collision-prone lookups (_cycles.png before the generic .png) resolve
# correctly by construction. Deliberately broader than AMPLICON_FILE_EXT_MAP:
# CoRAL's plotting code names the graph-plot image "{sname}_amplicon{N}_graph.png"
# / "_graph.pdf" (infixed) rather than AA's bare "{sname}_amplicon{N}.png" /
# ".pdf" — both spellings map to the same "png"/"pdf" key here, and get
# normalized to the bare canonical form by Stage 5's copy.
AMPLICON_FILE_DISCOVERY_SUFFIXES: List[Tuple[str, str]] = [
    ("_cycles.png", "cycles_png"),
    ("_cycles.pdf", "cycles_pdf"),
    ("_cycles.txt", "cycles"),
    ("_graph.txt",  "graph"),
    ("_graph.png",  "png"),
    ("_graph.pdf",  "pdf"),
    (".pdf",        "pdf"),
    (".png",        "png"),
]

# AC output files/dirs to merge into consolidated_classification/
# Each entry: suffix, has_header (True/False/None=optional file), is_dir,
# and optionally prefix_output (default True) — set False for AC outputs
# whose on-disk name is NOT prefixed by the AC output name (e.g.
# bfbarchitect_outputs/, which AC always names literally regardless of -o) —
# and compress (default False) — set True to emit the merged directory as a
# single .tar.gz rather than a loose subdirectory (bfbarchitect_outputs can
# be large; compressing it keeps the aggregated tree tidy). A compressed
# target round-trips on reaggregation because Stage 2's nested-archive pass
# expands the emitted .tar.gz back into a plain dir before discovery.
AC_MERGE_TARGETS: List[Dict] = [
    {"suffix": "_amplicon_classification_profiles.tsv", "has_header": True,  "is_dir": False},
    {"suffix": "_annotated_cycles_files",               "has_header": None,  "is_dir": True,  "rescue_suffix": "_annotated_cycles.txt"},
    {"suffix": "bfbarchitect_outputs",                  "has_header": None,  "is_dir": True,  "prefix_output": False, "compress": True},
    {"suffix": "_classification_bed_files",             "has_header": None,  "is_dir": True,  "rescue_suffix": "_intervals.bed"},
    {"suffix": "_ecDNA_context_calls.tsv",              "has_header": False, "is_dir": False},  # no header, may not exist
    {"suffix": "_ecDNA_counts.tsv",                     "has_header": True,  "is_dir": False},
    {"suffix": "_fan_calls.tsv",                        "has_header": True,  "is_dir": False},  # AC2+, may not exist
    {"suffix": "_feature_basic_properties.tsv",         "has_header": True,  "is_dir": False},
    {"suffix": "_feature_complexity.tsv",               "has_header": True,  "is_dir": False},  # AC2+ rename of _feature_entropy.tsv
    {"suffix": "_feature_entropy.tsv",                  "has_header": True,  "is_dir": False},  # pre-AC2
    {"suffix": "_feature_similarity_scores.tsv",        "has_header": True,  "is_dir": False},  # AC2+, may not exist
    {"suffix": "_gene_list.tsv",                        "has_header": True,  "is_dir": False},
    {"suffix": "_lncRNA_list.tsv",                      "has_header": True,  "is_dir": False},  # AC2+, may not exist
    {"suffix": ".log",                                  "has_header": False, "is_dir": False},  # AC's own log (top-level only; scraped for AC version, see _scan_classification_log_version)
    {"suffix": "_result_table.tsv",                     "has_header": True,  "is_dir": False},
    {"suffix": "_SV_summaries",                         "has_header": None,  "is_dir": True,  "rescue_suffix": "_SV_summary.tsv"},
]

# Columns in result_table.tsv that hold file paths needing re-resolution
PATH_COLUMNS: Tuple[str, ...] = (
    "Feature BED file",
    "CNV BED file",
    "Graph PNG file",
    "Graph PDF file",
    "Cycles PNG file",
    "Cycles PDF file",
    "Cycles file",
    "Graph file",
    "Run metadata JSON",
    "Sample metadata JSON",
)

# List-valued columns in result_table / run.json
LIST_COLUMNS: Tuple[str, ...] = ("Location", "Oncogenes", "All genes")

# Columns for aggregated_results.tsv — specific order, excludes reconstruction/cnvkit directory.
# List fields are rendered as Python repr strings (e.g. ['EGFR', 'MYC']).
AGG_CSV_COLUMNS: Tuple[str, ...] = (
    "Sample name", "AA amplicon number", "Feature ID", "Classification",
    "Contains viral", "Reconstruction tool",
    "AmpliconArchitect version", "AmpliconSuite-pipeline version",
    "AmpliconClassifier version",
    "Location", "Oncogenes", "Complexity score", "Captured interval length",
    "Feature median copy number", "Feature maximum copy number", "Filter flag",
    "Reference version", "Tissue of origin", "Sample type",
    "Feature BED file", "CNV BED file",
    "Graph PNG file", "Graph PDF file", "Cycles PNG file", "Cycles PDF file",
    "Cycles file", "Graph file",
    "Run metadata JSON", "Sample metadata JSON", "All genes",
)

# Canonical column order for run.json feature dicts (mirrors result_table + extras)
RUN_JSON_COLUMNS: Tuple[str, ...] = (
    "Sample name",
    "AA amplicon number",
    "Feature ID",
    "Classification",
    "Contains viral",
    "Reconstruction tool",
    "AmpliconArchitect version",
    "AmpliconSuite-pipeline version",
    "AmpliconClassifier version",
    "Location",
    "Oncogenes",
    "All genes",
    "Complexity score",
    "Captured interval length",
    "Feature median copy number",
    "Feature maximum copy number",
    "Filter flag",
    "Reference version",
    "Tissue of origin",
    "Sample type",
    "Feature BED file",
    "CNV BED file",
    "Graph PNG file",
    "Graph PDF file",
    "Cycles PNG file",
    "Cycles PDF file",
    "Cycles file",
    "Graph file",
    "Run metadata JSON",
    "Sample metadata JSON",
    "Reconstruction directory",
    "cnvkit directory",
)

# Sentinel value used when a file cannot be located
NOT_PROVIDED = "Not Provided"

# Name of the working extraction directory
EXTRACTION_DIR = "extracted_from_zips"

# Name of the staging results directory
RESULTS_DIR = "results"

# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class SampleRecord:
    """
    All discovered filesystem resources associated with a single sample.
    Populated during Stage 3 (Discovery). Consumed in Stages 5 & 6.
    """
    name: str

    # Directories (absolute paths in the extraction tree)
    aa_results_dir: Optional[str] = None       # [name]_AA_results/  (or synthesised)
    cnvkit_dir: Optional[str] = None           # [name]_cnvkit_output/ or _outputs/

    # Individual files (absolute paths)
    cnv_calls_bed: Optional[str] = None        # [name]_CNV_CALLS.bed  (canonical copy)
    cnv_calls_unfiltered_gains_bed: Optional[str] = None  # [name]_CNV_CALLS_unfiltered_gains.bed
    run_metadata_json: Optional[str] = None
    sample_metadata_json: Optional[str] = None
    aa_cnv_seeds_bed: Optional[str] = None
    finish_flag: Optional[str] = None
    timing_log: Optional[str] = None
    pipeline_log: Optional[str] = None         # top-level [name].log

    # Tool versions. Primary source is the sample's own result_table.tsv row
    # ("AS-p version"/"AA version"/"AC version" columns, written by AC's
    # make_results_table.py from run_metadata.json) when populated; falls
    # back to log content scraping otherwise (see extract_tool_versions() and
    # Aggregator._extract_log_versions()/_scan_classification_log_version()).
    # The result_table columns are reliably populated only when AC was run
    # as part of a single combined AmpliconSuite-pipeline.py --run_AC
    # invocation — the common "reclassify an existing AA cohort" workflow
    # leaves them "NA", which is exactly what the log-scrape fallback covers.
    aa_version: Optional[str] = None
    amplicon_suite_pipeline_version: Optional[str] = None
    ac_version: Optional[str] = None

    # AA/CoRAL per-amplicon files discovered inside aa_results_dir
    # { amplicon_num (int) -> { 'pdf': path, 'png': path, 'cycles': path,
    #                           'graph': path, 'cycles_png': path, 'cycles_pdf': path } }
    amplicon_files: Dict = field(default_factory=dict)

    # Summary file inside aa_results_dir
    aa_summary_file: Optional[str] = None

    # Classification dir where this sample's AC output originated (set in
    # Stage 4 from the result_table.tsv's own location; also where
    # _scan_classification_log_version() looks for AC's own log)
    ac_source_dir: Optional[str] = None

    # Reconstruction tool that produced this sample's graph/cycles files.
    # Defaults to AmpliconArchitect; upgraded to "CoRAL" during discovery
    # when a content-based CoRAL signal is found (see
    # Aggregator._detect_reconstruction_tool in asa_stages.py).
    reconstruction_tool: str = "AmpliconArchitect"

    # Output-tree destination paths — set by Stage 5, consumed by Stage 6
    aa_dir_dest:    Optional[str] = None  # results/samples/[s]/[s]_reconstruction_results/ (uncompressed dir)
    cnvkit_tarball: Optional[str] = None  # results/samples/[s]/[s]_cnvkit_output.tar.gz
    cnv_bed_dest:   Optional[str] = None  # results/samples/[s]/[s]_CNV_CALLS.bed (uncompressed)


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def rchop(s: str, suffix: str) -> str:
    """Remove a known suffix from a string."""
    if suffix and s.endswith(suffix):
        return s[:-len(suffix)]
    return s


def not_provided(value: object) -> bool:
    """Return True if a value represents a missing/not-provided path."""
    if value is None:
        return True
    s = str(value).strip()
    return s.lower() in ("not provided", "not_provided", "na", "nan", "none", "")


def parse_list_field(value: object) -> List[str]:
    """
    Convert a stringified list field from result_table into an actual Python list.
    Handles both pipe-separated ('a|b|c') and bracket-formatted ("['a', 'b']") forms.
    Returns an empty list if the value is missing/not-provided.
    """
    if not_provided(value):
        return []
    s = str(value).strip()
    if "|" in s:
        return [x.strip() for x in s.split("|") if x.strip()]
    # Strip surrounding brackets and split on commas
    s = s.strip("[]")
    parts = [x.strip().strip("'\"") for x in s.split(",")]
    return [p for p in parts if p]


def read_name_map(name_map_file: Optional[str]) -> Dict[str, str]:
    """
    Read a two-column sample-rename file.
    Returns a dict mapping old_name -> new_name.
    Blank lines and lines starting with '#' are ignored.
    """
    remap: Dict[str, str] = {}
    if not name_map_file:
        return remap
    try:
        with open(name_map_file) as fh:
            for lineno, line in enumerate(fh, 1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 2:
                    print(f"Warning: name_map line {lineno} has fewer than 2 columns, skipping: {line!r}")
                    continue
                remap[parts[0]] = parts[1]
    except OSError as e:
        sys.stderr.write(f"Error reading name_map '{name_map_file}': {e}\n")
        sys.exit(1)
    print(f"Name map loaded: {len(remap)} entry/entries.")
    return remap


def extract_tool_versions(text: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """
    Scan AA/AmpliconSuite-pipeline/AmpliconClassifier log text for version
    banner lines, e.g. "AmpliconArchitect version 1.5.r1", "AmpiconSuite-
    pipeline version 1.3.6", and "AmpliconClassifier 2.0.0". Returns
    (aa_version, pipeline_version, ac_version); any may be None if its
    banner line isn't present in this particular log.
    """
    aa_m = AA_VERSION_RE.search(text)
    pipeline_m = AS_PIPELINE_VERSION_RE.search(text)
    ac_m = AC_VERSION_RE.search(text)
    return (aa_m.group(1) if aa_m else None,
            pipeline_m.group(1) if pipeline_m else None,
            ac_m.group(1) if ac_m else None)


def is_aa_summary_content(first_line: str) -> bool:
    """Check whether a summary file's first line matches AA's convention."""
    return first_line.startswith("#Amplicons")


def is_coral_summary_content(lines: List[str]) -> bool:
    """
    Check whether a summary file's leading lines match CoRAL's convention.
    CoRAL summaries may open with a "CoRAL v{version}" banner line before
    the "N/M amplicons solved." line, or omit the banner entirely depending
    on version — scan the first few lines rather than assuming line 1.
    """
    for line in lines[:5]:
        if line.startswith(CORAL_BANNER_PREFIX) or CORAL_SUMMARY_RE.search(line):
            return True
    return False


def is_valid_aa_results_dir(dirpath: str) -> bool:
    """
    Validate that a directory is a genuine AA or CoRAL results directory.

    Rules:
      - Must contain exactly one *_summary.txt file (this also matches
        CoRAL's older *_amplicon_summary.txt naming, which ends in the
        same suffix).
      - That file's leading content must match either AA's or CoRAL's
        summary convention (see is_aa_summary_content/is_coral_summary_content).
      - The directory must NOT be named 'files' (AC copy directory).
    """
    if os.path.basename(dirpath) == "files":
        return False
    try:
        entries = os.listdir(dirpath)
    except OSError:
        return False

    summaries = [f for f in entries if f.endswith("_summary.txt")]
    if len(summaries) != 1:
        if len(summaries) > 1:
            print(f"Warning: multiple _summary.txt files in {dirpath} — skipping as AA results dir.")
        return False

    summary_path = os.path.join(dirpath, summaries[0])
    try:
        with open(summary_path) as fh:
            lines = [fh.readline() for _ in range(5)]
        return is_aa_summary_content(lines[0]) or is_coral_summary_content(lines)
    except OSError as e:
        print(f"Warning: could not read {summary_path}: {e}")
        return False


def is_classification_dir(dirpath: str) -> Tuple[bool, Optional[str], Optional[str]]:
    """
    Check whether a directory contains AC classification output.

    Returns (is_classification_dir, profiles_file_path, result_table_path).
    Prints a warning if a profiles file is found without a matching result table.
    Profiles and result_table must share the same prefix.
    """
    try:
        entries = set(os.listdir(dirpath))
    except OSError:
        return False, None, None

    profiles_files = [f for f in entries if f.endswith(AC_PROFILES_SUFFIX) and not f.startswith(".")]
    if not profiles_files:
        return False, None, None

    for pf in profiles_files:
        prefix = rchop(pf, AC_PROFILES_SUFFIX)
        expected_rt = prefix + AC_RESULT_TABLE_SUFFIX
        profiles_path = os.path.join(dirpath, pf)
        if expected_rt in entries:
            result_table_path = os.path.join(dirpath, expected_rt)
            return True, profiles_path, result_table_path
        else:
            print(f"Warning: found '{pf}' in {dirpath} but no matching '{expected_rt}' — skipping this classification dir.")

    return False, None, None


def make_tarball(source_dir: str, dest_tar_path: str,
                 exclusions: Tuple[str, ...] = EXCLUSION_SUFFIXES) -> None:
    """
    Create a .tar.gz archive of source_dir at dest_tar_path.
    Files matching any suffix in exclusions are omitted.
    The archive is rooted at the basename of source_dir (i.e. extracting
    produces [basename]/ not a full absolute path tree).
    """
    parent = os.path.dirname(source_dir)
    with tarfile.open(dest_tar_path, "w:gz") as tar:
        for root, dirs, files in os.walk(source_dir):
            dirs[:] = [d for d in dirs if not d.startswith(".") and d != "__MACOSX"]
            for fname in files:
                if fname.startswith(".") or any(fname.endswith(excl) for excl in exclusions):
                    continue
                fpath = os.path.join(root, fname)
                arcname = os.path.relpath(fpath, parent)
                tar.add(fpath, arcname=arcname)


def convert_cnvkit_cns_to_bed(cns_path: str, dest_path: str, min_cn: float = 0.0) -> None:
    """
    Convert a plain CNVkit .cns segment file (not .call.cns/.bintest.cns)
    into a CNV_CALLS.bed, recovering absolute copy number from the log2
    ratio column. Mirrors PrepareAA/scripts/convert_cns_to_bed.py's
    conversion formula (2**(log2_ratio + 1)) — used as a fallback for
    tools like CoRAL that produce raw cnvkit segments but, unlike
    AmpliconSuite-pipeline, never run this conversion step themselves.
    """
    with open(cns_path) as infile, open(dest_path, "w") as outfile:
        next(infile)  # header
        for line in infile:
            fields = line.rstrip("\n").split("\t")
            s, e = int(fields[1]), int(fields[2])
            cn = 2 ** (float(fields[4]) + 1)
            if not cn < min_cn:
                outfile.write("\t".join(fields[0:3] + [f"size={e - s}", str(cn)]) + "\n")


def safe_copy_file(src: str, dest: str) -> bool:
    """
    Copy a file from src to dest, creating parent directories as needed.
    Returns True on success, False on failure.
    """
    try:
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        shutil.copy2(src, dest)
        return True
    except OSError as e:
        print(f"Warning: could not copy {src} -> {dest}: {e}")
        return False


def safe_copytree(src: str, dest: str,
                  exclusions: Tuple[str, ...] = (),
                  inclusions: Tuple[str, ...] = ()) -> bool:
    """
    Recursively copy a directory tree from src to dest.
    Hidden names (starting with '.') are always skipped.
    If inclusions is non-empty, only files whose names end with one of those
    suffixes are copied (directories are always descended into).
    Otherwise, files whose names end with any suffix in exclusions are skipped.
    Returns True on success, False on failure.
    """
    def _ignore(directory: str, contents: list) -> list:
        ignored = []
        for f in contents:
            if f.startswith("."):
                ignored.append(f)
            elif inclusions:
                full = os.path.join(directory, f)
                if os.path.isfile(full) and not any(f.endswith(s) for s in inclusions):
                    ignored.append(f)
            elif any(f.endswith(e) for e in exclusions):
                ignored.append(f)
        return ignored

    try:
        shutil.copytree(src, dest, dirs_exist_ok=True, ignore=_ignore)
        return True
    except Exception as e:
        print(f"Warning: could not copy tree {src} -> {dest}: {e}")
        return False


def relative_to_results(abs_path: str, results_dir: str) -> str:
    """
    Express abs_path as a path relative to results_dir.
    Returns NOT_PROVIDED if abs_path does not start with results_dir.
    """
    try:
        rel = os.path.relpath(abs_path, results_dir)
        if rel.startswith(".."):
            return NOT_PROVIDED
        return rel
    except ValueError:
        return NOT_PROVIDED