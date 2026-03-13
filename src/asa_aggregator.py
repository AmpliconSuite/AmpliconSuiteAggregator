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
import shutil
import sys
import tarfile
import zipfile
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

__version__ = "6.0.0"

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

# Additional exclusions applied only when copying a sample's AA results directory
AA_DIR_EXCLUSION_SUFFIXES: Tuple[str, ...] = EXCLUSION_SUFFIXES + (
    "_cnseg.txt", "_edges.txt", "_logs.txt", ".out",
)

# Additional exclusions applied only when copying AA results directories
AA_DIR_EXCLUSION_SUFFIXES: Tuple[str, ...] = EXCLUSION_SUFFIXES + (
    "_cnseg.txt", "_edges.txt", "_logs.txt", ".out",
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
    "_finish_flag.txt",
    "_run_metadata.json",
    "_sample_metadata.json",
    "_timing_log.txt",
    # NOTE: top-level [sample_name].log is handled separately because the
    # suffix alone (".log") is too generic — we match by exact stem.
)

# AC output files/dirs to merge into consolidated_classification/
# Each entry: suffix, has_header (True/False/None=optional file), is_dir
AC_MERGE_TARGETS: List[Dict] = [
    {"suffix": "_amplicon_classification_profiles.tsv", "has_header": True,  "is_dir": False},
    {"suffix": "_annotated_cycles_files",               "has_header": None,  "is_dir": True},
    {"suffix": "_classification_bed_files",             "has_header": None,  "is_dir": True},
    {"suffix": "_ecDNA_context_calls.tsv",              "has_header": False, "is_dir": False},  # no header, may not exist
    {"suffix": "_ecDNA_counts.tsv",                     "has_header": True,  "is_dir": False},
    {"suffix": "_feature_basic_properties.tsv",         "has_header": True,  "is_dir": False},
    {"suffix": "_feature_entropy.tsv",                  "has_header": True,  "is_dir": False},
    {"suffix": "_gene_list.tsv",                        "has_header": True,  "is_dir": False},
    {"suffix": "_result_table.tsv",                     "has_header": True,  "is_dir": False},
    {"suffix": "_SV_summaries",                         "has_header": None,  "is_dir": True},
]

# Columns in result_table.tsv that hold file paths needing re-resolution
PATH_COLUMNS: Tuple[str, ...] = (
    "Feature BED file",
    "CNV BED file",
    "AA PNG file",
    "AA PDF file",
    "AA cycles file",
    "AA graph file",
    "Run metadata JSON",
    "Sample metadata JSON",
)

# List-valued columns in result_table / run.json
LIST_COLUMNS: Tuple[str, ...] = ("Location", "Oncogenes", "All genes")

# Columns for aggregated_results.tsv — specific order, excludes AA/cnvkit directory.
# List fields are rendered as Python repr strings (e.g. ['EGFR', 'MYC']).
AGG_CSV_COLUMNS: Tuple[str, ...] = (
    "Sample name", "AA amplicon number", "Feature ID", "Classification",
    "Location", "Oncogenes", "Complexity score", "Captured interval length",
    "Feature median copy number", "Feature maximum copy number", "Filter flag",
    "Reference version", "Tissue of origin", "Sample type",
    "Feature BED file", "CNV BED file", "AA PNG file", "AA PDF file",
    "AA cycles file", "AA graph file",
    "Run metadata JSON", "Sample metadata JSON", "All genes",
)

# Canonical column order for run.json feature dicts (mirrors result_table + extras)
RUN_JSON_COLUMNS: Tuple[str, ...] = (
    "Sample name",
    "AA amplicon number",
    "Feature ID",
    "Classification",
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
    "AA PNG file",
    "AA PDF file",
    "AA cycles file",
    "AA graph file",
    "Run metadata JSON",
    "Sample metadata JSON",
    "AA directory",
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
    run_metadata_json: Optional[str] = None
    sample_metadata_json: Optional[str] = None
    aa_cnv_seeds_bed: Optional[str] = None
    finish_flag: Optional[str] = None
    timing_log: Optional[str] = None
    pipeline_log: Optional[str] = None         # top-level [name].log

    # AA per-amplicon files discovered inside aa_results_dir
    # { amplicon_num (int) -> { 'pdf': path, 'png': path, 'cycles': path, 'graph': path } }
    amplicon_files: Dict = field(default_factory=dict)

    # Summary file inside aa_results_dir
    aa_summary_file: Optional[str] = None

    # Classification dir where this sample's AC output originated
    ac_source_dir: Optional[str] = None

    # Output-tree destination paths — set by Stage 5, consumed by Stage 6
    aa_dir_dest:    Optional[str] = None  # results/samples/[s]/[s]_AA_results/ (uncompressed dir)
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


def is_valid_aa_results_dir(dirpath: str) -> bool:
    """
    Validate that a directory is a genuine AA results directory.

    Rules:
      - Must contain exactly one *_summary.txt file.
      - That file's first line must start with '#Amplicons'.
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
            first_line = fh.readline()
        return first_line.startswith("#Amplicons")
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

    profiles_files = [f for f in entries if f.endswith(AC_PROFILES_SUFFIX) and not f.startswith("._")]
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
            dirs[:] = [d for d in dirs if d not in ("__MACOSX", ".DS_Store")]
            for fname in files:
                if any(fname.endswith(excl) for excl in exclusions):
                    continue
                fpath = os.path.join(root, fname)
                arcname = os.path.relpath(fpath, parent)
                tar.add(fpath, arcname=arcname)


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
                  exclusions: Tuple[str, ...] = ()) -> bool:
    """
    Recursively copy a directory tree from src to dest,
    skipping files whose names end with any suffix in exclusions.
    Returns True on success, False on failure.
    """
    def _ignore(directory: str, contents: list) -> list:
        return [f for f in contents
                if any(f.endswith(e) for e in exclusions)]

    try:
        shutil.copytree(src, dest, dirs_exist_ok=True,
                        ignore=_ignore if exclusions else None)
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