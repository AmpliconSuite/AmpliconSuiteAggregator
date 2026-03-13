"""
asa_stages.py
The Aggregator class — pipeline orchestration and all stage implementations.

Imports shared constants, data structures, and utilities from asa_aggregator.py.

Stages:
  2  — Extraction engine
  3  — Discovery engine
  4  — result_table parser + sample registry reconciliation  (this segment)
  5  — Output tree builder
  6  — run.json constructor
  F  — Finalise (tarball creation)
"""

from __future__ import annotations

import json
import time
import os
import shutil
import sys
import tarfile
import zipfile
from typing import Dict, List, Optional, Tuple

import pandas as pd

from asa_aggregator import (
    # constants
    ARCHIVE_EXTENSIONS, EXCLUSION_SUFFIXES, AA_DIR_EXCLUSION_SUFFIXES, MISC_AA_SUFFIXES,
    AC_MERGE_TARGETS, RUN_JSON_COLUMNS, AGG_CSV_COLUMNS, LIST_COLUMNS,
    NOT_PROVIDED, EXTRACTION_DIR, RESULTS_DIR,
    AC_PROFILES_SUFFIX, AC_RESULT_TABLE_SUFFIX,
    # data structures
    SampleRecord,
    # utilities
    rchop, not_provided, parse_list_field, read_name_map,
    is_valid_aa_results_dir, is_classification_dir,
    make_tarball, safe_copy_file, safe_copytree, relative_to_results,
)

from asa_aggregator import __version__


# ---------------------------------------------------------------------------
# Aggregator class
# ---------------------------------------------------------------------------

class Aggregator:
    """
    Orchestrates the full aggregation pipeline.

    Instance attributes (set in __init__, used across stages):
      input_paths       — validated list of input archive/directory paths
      project_name      — prefix for consolidated classification files and output archive
      name_map          — { old_sname -> new_sname } rename dict
      no_cleanup        — when True, temp dirs are preserved after completion
      work_dir          — absolute cwd at construction time
      extract_dir       — <work_dir>/extracted_from_zips/
      results_dir       — <work_dir>/results/
      samples_dir       — <results_dir>/samples/
      classif_dir       — <results_dir>/consolidated_classification/
      other_dir         — <results_dir>/other_files/
      sample_registry   — { sname -> SampleRecord }  (populated by Stage 3)
      classification_dirs — [ dirpath ]               (populated by Stage 3)
      aux_dirs          — [ dirpath ]                 (populated by Stage 3)
      _files_dirs       — { cls_dir -> files_subdir } (populated by Stage 3)
      _floating_cnv_beds — { sname -> path }          (populated by Stage 3)
      completed         — True only after successful _finalise()
    """

    def __init__(
        self,
        input_paths: List[str],
        project_name: str,
        name_map_file: Optional[str] = None,
        no_cleanup: bool = False,
        work_dir: Optional[str] = None,
    ):
        self.input_paths = input_paths
        self.project_name = project_name
        self.name_map = read_name_map(name_map_file)
        self.no_cleanup = no_cleanup
        self.completed = False
        self._start_time: float = time.perf_counter()

        self.work_dir    = os.path.abspath(work_dir) if work_dir else os.path.abspath(".")
        self.extract_dir = os.path.join(self.work_dir, EXTRACTION_DIR)
        self.results_dir = os.path.join(self.work_dir, RESULTS_DIR)
        self.samples_dir = os.path.join(self.results_dir, "samples")
        self.classif_dir = os.path.join(self.results_dir, "consolidated_classification")
        self.other_dir   = os.path.join(self.results_dir, "other_files")

        self.sample_registry:    Dict[str, SampleRecord] = {}
        self.classification_dirs: List[str] = []
        self.aux_dirs:            List[str] = []
        self._files_dirs:         Dict[str, str] = {}
        self._floating_cnv_beds:  Dict[str, str] = {}

        self._run_pipeline()

    # ------------------------------------------------------------------
    # Pipeline orchestration
    # ------------------------------------------------------------------

    def _run_pipeline(self) -> None:
        try:
            self._setup_work_dirs()
            self._stage2_extract()
            self._stage3_discover()
            self._stage4_parse_result_tables()
            self._stage5_build_output_tree()
            self._stage6_build_run_json()
            self._finalise()
            self.completed = True
        except SystemExit:
            self._cleanup(failure=True)
            raise
        except Exception:
            import traceback
            sys.stderr.write("\nUnexpected error during aggregation:\n")
            traceback.print_exc()
            self._cleanup(failure=True)

    # ------------------------------------------------------------------
    # Cross-stage helpers
    # ------------------------------------------------------------------

    def _setup_work_dirs(self) -> None:
        for d in [self.extract_dir, self.results_dir]:
            if os.path.exists(d):
                print(f"Warning: '{d}' already exists — removing it now.")
                shutil.rmtree(d)
        for d in [self.extract_dir, self.results_dir,
                  self.samples_dir, self.classif_dir, self.other_dir]:
            os.makedirs(d, exist_ok=True)
        print("Working directories initialised.")
        print(f"  Extraction dir : {self.extract_dir}")
        print(f"  Results dir    : {self.results_dir}")

    def _cleanup(self, failure: bool = False) -> None:
        if self.no_cleanup:
            print("--no_cleanup set: leaving working directories in place.")
            return
        if os.path.exists(self.extract_dir):
            shutil.rmtree(self.extract_dir, ignore_errors=True)
        if os.path.exists(self.results_dir):
            shutil.rmtree(self.results_dir, ignore_errors=True)

    def _abort(self, message: str) -> None:
        sys.stderr.write(f"FATAL: {message}\n")
        self._cleanup(failure=True)
        sys.exit(1)

    # ==================================================================
    # Stage 2 — Extraction engine
    # ==================================================================

    def _stage2_extract(self) -> None:
        print("\n--- Stage 2: Extraction ---")
        self._input_size_bytes = 0
        for path in self.input_paths:
            abs_path_size = os.path.abspath(path)
            if os.path.isfile(abs_path_size):
                self._input_size_bytes += os.path.getsize(abs_path_size)
            elif os.path.isdir(abs_path_size):
                for dirpath, _, fnames in os.walk(abs_path_size):
                    for fn in fnames:
                        try:
                            self._input_size_bytes += os.path.getsize(os.path.join(dirpath, fn))
                        except OSError:
                            pass
        input_mb = self._input_size_bytes / (1024 * 1024)
        print(f"  Input size: {input_mb:.2f} MB ({len(self.input_paths)} source(s))")
        for path in self.input_paths:
            abs_path = os.path.abspath(path)
            if os.path.isdir(abs_path):
                self._copy_input_dir(abs_path)
            elif os.path.isfile(abs_path):
                self._extract_archive(abs_path, self.extract_dir)
            else:
                print(f"Warning: input path does not exist or is not accessible: {abs_path}")

        pass_num = 0
        while True:
            pass_num += 1
            nested = self._find_nested_archives()
            if not nested:
                break
            print(f"  Nested archive pass {pass_num}: found {len(nested)} archive(s) to expand.")
            for archive_path in nested:
                dest_dir = os.path.dirname(archive_path)
                self._extract_archive(archive_path, dest_dir)
                try:
                    os.remove(archive_path)
                except OSError as e:
                    print(f"  Warning: could not remove nested archive {archive_path}: {e}")

        print(f"  Extraction complete. Working tree: {self.extract_dir}")

    def _copy_input_dir(self, src_dir: str) -> None:
        basename = os.path.basename(src_dir.rstrip("/"))
        dest = self._unique_dest(self.extract_dir, basename)
        print(f"  Copying directory: {src_dir} -> {dest}")
        try:
            shutil.copytree(src_dir, dest,
                            ignore=shutil.ignore_patterns("__MACOSX", ".DS_Store"))
        except Exception as e:
            print(f"  Warning: error copying directory {src_dir}: {e}")

    def _extract_archive(self, archive_path: str, dest_parent: str) -> Optional[str]:
        fname = os.path.basename(archive_path)
        if fname.endswith(".tar.gz"):
            stem = fname[:-7]
        elif fname.endswith(".tar"):
            stem = fname[:-4]
        elif fname.endswith(".zip"):
            stem = fname[:-4]
        else:
            print(f"  Warning: unrecognised archive type, skipping: {archive_path}")
            return None

        dest = self._unique_dest(dest_parent, stem)
        os.makedirs(dest, exist_ok=True)
        try:
            if archive_path.endswith(".zip"):
                with zipfile.ZipFile(archive_path, "r") as zf:
                    members = [m for m in zf.namelist()
                               if not m.startswith("__MACOSX") and ".DS_Store" not in m]
                    zf.extractall(dest, members=members)
            else:
                with tarfile.open(archive_path, "r:*") as tf:
                    members = [m for m in tf.getmembers()
                               if not m.name.startswith("__MACOSX")
                               and ".DS_Store" not in m.name]
                    tf.extractall(dest, members=members)
            print(f"  Extracted: {archive_path} -> {dest}")
            return dest
        except Exception as e:
            print(f"  Warning: failed to extract {archive_path}: {e}")
            if os.path.exists(dest) and not os.listdir(dest):
                shutil.rmtree(dest, ignore_errors=True)
            return None

    def _find_nested_archives(self) -> List[str]:
        archives = []
        for root, dirs, files in os.walk(self.extract_dir):
            dirs[:] = [d for d in dirs
                       if d not in ("__MACOSX",) and ".DS_Store" not in d]
            for fname in files:
                if any(fname.endswith(ext) for ext in ARCHIVE_EXTENSIONS):
                    archives.append(os.path.join(root, fname))
        return archives

    @staticmethod
    def _unique_dest(parent: str, name: str) -> str:
        candidate = os.path.join(parent, name)
        if not os.path.exists(candidate):
            return candidate
        counter = 2
        while True:
            candidate = os.path.join(parent, f"{name}_{counter}")
            if not os.path.exists(candidate):
                return candidate
            counter += 1

    # ==================================================================
    # Stage 3 — Discovery engine
    # ==================================================================

    def _stage3_discover(self) -> None:
        print("\n--- Stage 3: Discovery ---")
        self._files_dirs = {}
        self._floating_cnv_beds = {}

        for root, dirs, files in os.walk(self.extract_dir, topdown=True):
            dirs[:] = [d for d in dirs
                       if d not in ("__MACOSX",) and ".DS_Store" not in d]

            # 3a. Classify sub-directories
            for dname in list(dirs):
                dpath = os.path.join(root, dname)
                try:
                    dcontents = set(os.listdir(dpath))
                except OSError:
                    dcontents = set()

                if "AUX_DIR" in dcontents:
                    print(f"  AUX dir found: {dpath}")
                    self.aux_dirs.append(dpath)
                    dirs.remove(dname)
                    continue

                if dname.endswith("_AA_results"):
                    if is_valid_aa_results_dir(dpath):
                        sname = rchop(dname, "_AA_results")
                        self._register_aa_results_dir(sname, dpath)
                        dirs.remove(dname)
                        continue
                    else:
                        print(f"  Warning: '{dpath}' looks like _AA_results "
                              f"but failed validation — will descend.")
                        continue

                if is_valid_aa_results_dir(dpath):
                    print(f"  Warning: '{dpath}' passes AA validation but "
                          f"is not named _AA_results — registering anyway.")
                    sname = self._sname_from_summary(dpath) or dname
                    self._register_aa_results_dir(sname, dpath)
                    dirs.remove(dname)
                    continue

                if dname.endswith("_cnvkit_output") or dname.endswith("_cnvkit_outputs"):
                    suffix = ("_cnvkit_outputs" if dname.endswith("_cnvkit_outputs")
                              else "_cnvkit_output")
                    sname = rchop(dname, suffix)
                    rec = self._get_or_create_record(sname)
                    if rec.cnvkit_dir:
                        print(f"  Warning: duplicate cnvkit dir for '{sname}', "
                              f"keeping first: {rec.cnvkit_dir}")
                    else:
                        rec.cnvkit_dir = dpath
                        print(f"  cnvkit dir   : {sname} -> {dpath}")
                        for fname in dcontents:
                            if fname.endswith("_CNV_CALLS.bed"):
                                rec.cnv_calls_bed = os.path.join(dpath, fname)
                                break
                    dirs.remove(dname)
                    continue

                is_cls, _profiles, _rt = is_classification_dir(dpath)
                if is_cls:
                    self.classification_dirs.append(dpath)
                    print(f"  Classif dir  : {dpath}")
                    files_sub = os.path.join(dpath, "files")
                    if os.path.isdir(files_sub):
                        self._files_dirs[dpath] = files_sub
                    dirs.remove(dname)
                    self._walk_classification_dir(dpath)
                    continue

            # 3b. Scan files at this level
            for fname in files:
                fpath = os.path.join(root, fname)
                if fname.startswith("._"):
                    continue

                matched = False
                for suffix in MISC_AA_SUFFIXES:
                    if fname.endswith(suffix):
                        sname = rchop(fname, suffix)
                        rec = self._get_or_create_record(sname)
                        self._register_misc_file(rec, suffix, fpath)
                        matched = True
                        break

                if not matched:
                    if fname.endswith(".log") and root != self.extract_dir:
                        sname = fname[:-4]
                        sibling_aa = os.path.join(root, f"{sname}_AA_results")
                        if os.path.isdir(sibling_aa):
                            rec = self._get_or_create_record(sname)
                            if not rec.pipeline_log:
                                rec.pipeline_log = fpath

                    if fname.endswith("_CNV_CALLS.bed"):
                        parent_name = os.path.basename(root)
                        if not (parent_name.endswith("_cnvkit_output")
                                or parent_name.endswith("_cnvkit_outputs")):
                            sname = self._cnv_bed_sample_name(
                                fname, parent_name,
                                known_snames=set(self.sample_registry.keys()))
                            if sname not in self._floating_cnv_beds:
                                self._floating_cnv_beds[sname] = fpath
                                print(f"  Floating CNV : {sname} -> {fpath}")

        # 3c. Attach floating CNV beds to records that lack one
        for sname, bed_path in self._floating_cnv_beds.items():
            rec = self._get_or_create_record(sname)
            if not rec.cnv_calls_bed:
                rec.cnv_calls_bed = bed_path

        print(f"  Discovery complete: {len(self.sample_registry)} sample(s) inferred, "
              f"{len(self.classification_dirs)} classification dir(s), "
              f"{len(self.aux_dirs)} AUX dir(s).")

    def _walk_classification_dir(self, cls_dir: str) -> None:
        for root, dirs, files in os.walk(cls_dir, topdown=True):
            dirs[:] = [d for d in dirs
                       if d not in ("__MACOSX",) and ".DS_Store" not in d]
            for dname in list(dirs):
                dpath = os.path.join(root, dname)
                if dname == "files":
                    dirs.remove(dname)
                    continue
                if dname.endswith("_AA_results") and is_valid_aa_results_dir(dpath):
                    sname = rchop(dname, "_AA_results")
                    self._register_aa_results_dir(sname, dpath)
                    dirs.remove(dname)
                    continue
                if (dname.endswith("_cnvkit_output")
                        or dname.endswith("_cnvkit_outputs")):
                    suffix = ("_cnvkit_outputs" if dname.endswith("_cnvkit_outputs")
                              else "_cnvkit_output")
                    sname = rchop(dname, suffix)
                    rec = self._get_or_create_record(sname)
                    if not rec.cnvkit_dir:
                        rec.cnvkit_dir = dpath
                        try:
                            for f in os.listdir(dpath):
                                if f.endswith("_CNV_CALLS.bed"):
                                    rec.cnv_calls_bed = os.path.join(dpath, f)
                                    break
                        except OSError:
                            pass
                    dirs.remove(dname)
                    continue
            for fname in files:
                if fname.startswith("._"):
                    continue
                fpath = os.path.join(root, fname)
                for suffix in MISC_AA_SUFFIXES:
                    if fname.endswith(suffix):
                        sname = rchop(fname, suffix)
                        rec = self._get_or_create_record(sname)
                        self._register_misc_file(rec, suffix, fpath)
                        break

    def _register_aa_results_dir(self, sname: str, dpath: str) -> None:
        rec = self._get_or_create_record(sname)
        if rec.aa_results_dir:
            print(f"  Warning: duplicate _AA_results dir for '{sname}', "
                  f"keeping first: {rec.aa_results_dir}")
            return
        rec.aa_results_dir = dpath
        print(f"  AA results   : {sname} -> {dpath}")
        try:
            for fname in os.listdir(dpath):
                fpath = os.path.join(dpath, fname)
                if fname.endswith("_summary.txt"):
                    rec.aa_summary_file = fpath
                    continue
                for ext, key in ((".pdf", "pdf"), (".png", "png"),
                                 ("_cycles.txt", "cycles"), ("_graph.txt", "graph")):
                    if fname.endswith(ext) and "_amplicon" in fname:
                        num = self._parse_amplicon_num(fname, sname)
                        if num is not None:
                            rec.amplicon_files.setdefault(num, {})[key] = fpath
                        break
        except OSError as e:
            print(f"  Warning: could not index files in {dpath}: {e}")

    def _register_misc_file(self, rec: SampleRecord, suffix: str, fpath: str) -> None:
        if suffix == "_AA_CNV_SEEDS.bed":
            if not rec.aa_cnv_seeds_bed:    rec.aa_cnv_seeds_bed = fpath
        elif suffix == "_finish_flag.txt":
            if not rec.finish_flag:         rec.finish_flag = fpath
        elif suffix == "_run_metadata.json":
            if not rec.run_metadata_json:   rec.run_metadata_json = fpath
        elif suffix == "_sample_metadata.json":
            if not rec.sample_metadata_json: rec.sample_metadata_json = fpath
        elif suffix == "_timing_log.txt":
            if not rec.timing_log:          rec.timing_log = fpath

    def _get_or_create_record(self, sname: str) -> SampleRecord:
        if sname not in self.sample_registry:
            self.sample_registry[sname] = SampleRecord(name=sname)
        return self.sample_registry[sname]

    @staticmethod
    def _sname_from_summary(dirpath: str) -> Optional[str]:
        try:
            for fname in os.listdir(dirpath):
                if fname.endswith("_summary.txt"):
                    return rchop(fname, "_summary.txt")
        except OSError:
            pass
        return None

    @staticmethod
    def _parse_amplicon_num(fname: str, sname: str) -> Optional[int]:
        stem = fname
        if stem.startswith(sname + "_amplicon"):
            stem = stem[len(sname) + len("_amplicon"):]
        elif "_amplicon" in stem:
            stem = stem.split("_amplicon", 1)[1]
        else:
            return None
        num_str = ""
        for ch in stem:
            if ch.isdigit():
                num_str += ch
            else:
                break
        try:
            return int(num_str)
        except ValueError:
            return None

    @staticmethod
    def _cnv_bed_sample_name(fname: str, parent_dir_name: str = "",
                              known_snames: Optional[set] = None) -> str:
        """
        Derive a sample name from a CNV BED filename.

        Handled patterns (in priority order):
          1. [sname]_CNV_CALLS.bed                    -> sname (no dots in prefix)
          2. [sname].cs.rmdup_CNV_CALLS.bed           -> sname if prefix matches a known sample
          3. [sname].[anything]_CNV_CALLS.bed         -> sname (truncate at first dot)
          4. [any_prefix].cs.rmdup_CNV_CALLS.bed
             where prefix is not a known sample name  -> parent_dir_name (directory fallback)
             e.g. /path/SAMPLE/ERR3345421.cs.rmdup_CNV_CALLS.bed -> "SAMPLE"

        Note: inside a [sname]_cnvkit_output/ dir, sname is already known from
        the directory name so _cnv_bed_sample_name is not used — any file ending
        in _CNV_CALLS.bed (including .cs.rmdup_CNV_CALLS.bed) is accepted directly.
        """
        prefix = rchop(fname, "_CNV_CALLS.bed")

        # Case 1: no dots — sample name is unambiguous
        if "." not in prefix:
            return prefix

        # Cases 2/3: prefix contains dots
        pre_dot = prefix.split(".")[0]

        # Case 2: check if the pre-dot part matches a known sample name
        if known_snames and pre_dot in known_snames:
            return pre_dot

        # Case 4: .cs.rmdup pattern with unrecognised prefix — use parent dir name
        if fname.endswith(".cs.rmdup_CNV_CALLS.bed") and parent_dir_name:
            return parent_dir_name

        # Case 3: fallback — truncate at first dot
        return pre_dot

    # ==================================================================
    # Stage 4 — result_table parser + sample registry reconciliation
    # ==================================================================

    def _stage4_parse_result_tables(self) -> None:
        """
        Find and parse all *_result_table.tsv files within the classification
        dirs discovered in Stage 3.  This is the authoritative source of:
          - the canonical sample name list
          - per-feature metadata (classification, location, oncogenes, etc.)

        After parsing, every result_table row is attached to the matching
        SampleRecord in self.sample_registry.  Records inferred by Stage 3
        that do NOT appear in any result_table trigger a warning.  Samples
        that appear in result_tables but have NO Stage 3 discovery data get a
        stub SampleRecord so later stages can still produce partial output.

        self.feature_rows is set to a flat list of dicts — one per feature row
        across all result tables — preserving the 'sample_N' grouping key used
        in run.json.
        """
        print("\n--- Stage 4: Result table parsing ---")

        # { 'sample_N' -> [ row_dict, ... ] }  preserves run.json index keys
        self.run_json_groups: Dict[str, List[dict]] = {}

        # Flat ordered list of (sample_key, sname, row_dict) for easy iteration
        self.all_feature_rows: List[Tuple[str, str, dict]] = []

        sample_counter = 1
        found_any = False

        # Search classification dirs first, then fall back to a broader walk
        # of the entire extraction tree (handles result tables that landed
        # outside a recognised classification dir).
        search_dirs = list(self.classification_dirs)
        if not search_dirs:
            search_dirs = [self.extract_dir]

        seen_result_tables: List[str] = []

        for search_root in search_dirs:
            for root, dirs, files in os.walk(search_root):
                dirs[:] = [d for d in dirs
                           if d not in ("__MACOSX",) and ".DS_Store" not in d]
                for fname in files:
                    if not fname.endswith("_result_table.tsv"):
                        continue
                    if fname.startswith("._"):
                        continue
                    rt_path = os.path.join(root, fname)
                    if rt_path in seen_result_tables:
                        continue
                    seen_result_tables.append(rt_path)
                    self._parse_single_result_table(rt_path, sample_counter)
                    # Advance counter by however many new sample keys were added
                    sample_counter = len(self.run_json_groups) + 1
                    found_any = True

        if not found_any:
            self._abort(
                "No *_result_table.tsv files found. "
                "Ensure AmpliconClassifier has been run before aggregating."
            )

        # Flatten into self.all_feature_rows for Stage 5 / 6 iteration
        self.all_feature_rows = []
        for skey, rows in self.run_json_groups.items():
            for row in rows:
                self.all_feature_rows.append((skey, row["Sample name"], row))

        # Reconcile: warn about Stage 3 discoveries not in any result table
        rt_snames = {sname for _, sname, _ in self.all_feature_rows}
        for inferred_sname in list(self.sample_registry.keys()):
            if inferred_sname not in rt_snames:
                # Only warn for records that have substantive data (AA or cnvkit)
                rec = self.sample_registry[inferred_sname]
                if rec.aa_results_dir or rec.cnvkit_dir:
                    print(f"  Warning: sample '{inferred_sname}' has AA/cnvkit data "
                          f"but does not appear in any result_table.tsv — "
                          f"it will not be included in the output.")

        # Reconcile: create stub records for result_table samples with no Stage 3 data
        for sname in rt_snames:
            if sname not in self.sample_registry:
                print(f"  Warning: sample '{sname}' found in result_table "
                      f"but no AA/cnvkit data was discovered for it.")
                self.sample_registry[sname] = SampleRecord(name=sname)

        total_features = len(self.all_feature_rows)
        total_samples  = len(rt_snames)
        print(f"  Parsed {total_features} feature row(s) across "
              f"{total_samples} sample(s) from "
              f"{len(seen_result_tables)} result_table file(s).")

    def _parse_single_result_table(self, rt_path: str, start_counter: int) -> None:
        """
        Read one *_result_table.tsv, group rows by 'Sample name', and merge
        them into self.run_json_groups.

        Duplicate sample names (across multiple result tables) are handled by
        keeping the last-seen set of rows and printing a light warning.
        """
        print(f"  Reading: {rt_path}")
        try:
            df = pd.read_csv(rt_path, sep="\t", dtype=str, keep_default_na=False)
        except Exception as e:
            print(f"  Warning: could not read {rt_path}: {e} — skipping.")
            return

        if "Sample name" not in df.columns:
            print(f"  Warning: 'Sample name' column missing in {rt_path} — skipping.")
            return

        # Normalise NaN-like values to our sentinel string
        df = df.where(df != "", other=NOT_PROVIDED)

        # Find the reverse mapping sname -> existing sample_key so we can
        # overwrite on duplicate rather than creating a second key.
        existing_key_for: Dict[str, str] = {}
        for skey, rows in self.run_json_groups.items():
            if rows:
                existing_key_for[rows[0]["Sample name"]] = skey

        counter = start_counter
        for sname, group in df.groupby("Sample name", sort=False):
            sname = str(sname)  # guard: ensure key is always a string (e.g. numeric sample names)
            rows = group.to_dict(orient="records")
            # Also normalise the Sample name field inside each row
            for row in rows:
                row["Sample name"] = str(row["Sample name"])
            if sname in existing_key_for:
                print(f"  Warning: duplicate sample name '{sname}' — "
                      f"overwriting earlier result_table rows.")
                self.run_json_groups[existing_key_for[sname]] = rows
            else:
                skey = f"sample_{counter}"
                self.run_json_groups[skey] = rows
                counter += 1
            print(f"    {sname}: {len(rows)} feature row(s)")

    # ==================================================================
    # Stage stubs — implemented in subsequent segments
    # ==================================================================


    # ==================================================================
    # Stage 5 — Output tree builder
    # ==================================================================

    def _stage5_build_output_tree(self) -> None:
        """
        Build the physical results/ tree:

          results/samples/[sname]/
              [sname]_AA_results/            (uncompressed dir)
              [sname]_cnvkit_output.tar.gz
              [sname]_CNV_CALLS.bed          (uncompressed copy)
              [sname]_run_metadata.json      (if found)
              [sname]_sample_metadata.json   (if found)
              [sname]_AA_CNV_SEEDS.bed       (if found)
              [sname]_finish_flag.txt        (if found)
              [sname]_timing_log.txt         (if found)
              [sname].log                    (if found)

          results/consolidated_classification/
              [project_name]_amplicon_classification_profiles.tsv
              [project_name]_annotated_cycles_files/
              [project_name]_classification_bed_files/
              [project_name]_ecDNA_context_calls.tsv
              [project_name]_ecDNA_counts.tsv
              [project_name]_feature_basic_properties.tsv
              [project_name]_feature_entropy.tsv
              [project_name]_gene_list.tsv
              [project_name]_result_table.tsv
              [project_name]_SV_summaries/

          results/other_files/
              [AUX dir contents copied here]
        """
        print("\n--- Stage 5: Building output tree ---")

        # Authoritative sample name list comes from result_table parsing
        rt_snames = {sname for _, sname, _ in self.all_feature_rows}

        for sname in sorted(rt_snames):
            self._build_sample_dir(sname)

        self._build_consolidated_classification()
        self._copy_aux_dirs()

        print("  Output tree construction complete.")

    # ------------------------------------------------------------------
    # Stage 5a — per-sample directory
    # ------------------------------------------------------------------

    def _build_sample_dir(self, sname: str) -> None:
        """
        Create results/samples/[sname]/ and populate it with tarballs,
        an uncompressed CNV BED copy, and all misc files.
        """
        rec = self.sample_registry.get(sname)
        if rec is None:
            # Should not happen — Stage 4 always creates a stub
            rec = SampleRecord(name=sname)
            self.sample_registry[sname] = rec

        sample_out = os.path.join(self.samples_dir, sname)
        os.makedirs(sample_out, exist_ok=True)

        # ---- AA results directory (uncompressed copy) ------------------
        rec.aa_dir_dest = self._copy_aa_results_dir(sname, rec, sample_out)

        # ---- cnvkit tarball ---------------------------------------------
        rec.cnvkit_tarball = self._build_cnvkit_tarball(sname, rec, sample_out)

        # ---- Uncompressed CNV BED copy ----------------------------------
        rec.cnv_bed_dest = self._copy_cnv_bed(sname, rec, sample_out)

        # ---- Miscellaneous files ----------------------------------------
        misc_map = [
            (rec.run_metadata_json,    f"{sname}_run_metadata.json"),
            (rec.sample_metadata_json, f"{sname}_sample_metadata.json"),
            (rec.aa_cnv_seeds_bed,     f"{sname}_AA_CNV_SEEDS.bed"),
            (rec.finish_flag,          f"{sname}_finish_flag.txt"),
            (rec.timing_log,           f"{sname}_timing_log.txt"),
            (rec.pipeline_log,         f"{sname}.log"),
        ]
        for src, dest_name in misc_map:
            if src and os.path.isfile(src):
                dest = os.path.join(sample_out, dest_name)
                safe_copy_file(src, dest)

        print(f"  Built sample dir: {sname}")

    def _copy_aa_results_dir(self, sname: str, rec: SampleRecord,
                              sample_out: str) -> Optional[str]:
        """
        Copy AA results into sample_out/[sname]_AA_results/ (uncompressed).

        Priority:
          1. Use rec.aa_results_dir (canonical dir) — copy tree directly.
          2. Synthesise a [sname]_AA_results/ dir from individual files:
             a. amplicon_files (pdf, png, cycles, graph) from rec
             b. aa_summary_file from rec
             c. Fall back to files/ subdir of any classification dir for this sample
          3. If nothing found, return None.
        """
        dest_dir = os.path.join(sample_out, f"{sname}_AA_results")

        if rec.aa_results_dir and os.path.isdir(rec.aa_results_dir):
            safe_copytree(rec.aa_results_dir, dest_dir,
                          exclusions=AA_DIR_EXCLUSION_SUFFIXES)
            return dest_dir

        # Synthesise from individual files
        os.makedirs(dest_dir, exist_ok=True)
        found_any = False

        for amp_num, file_dict in rec.amplicon_files.items():
            for key, src in file_dict.items():
                if src and os.path.isfile(src):
                    shutil.copy2(src, dest_dir)
                    found_any = True

        if rec.aa_summary_file and os.path.isfile(rec.aa_summary_file):
            shutil.copy2(rec.aa_summary_file, dest_dir)
            found_any = True

        if not found_any:
            found_any = self._pull_aa_files_from_files_dirs(sname, dest_dir)

        if found_any:
            return dest_dir

        shutil.rmtree(dest_dir, ignore_errors=True)
        print(f"  Warning: no AA results found for '{sname}' — AA directory: Not Provided")
        return None

    def _pull_aa_files_from_files_dirs(self, sname: str, staging: str) -> bool:
        """
        Search all known files/ subdirs inside classification dirs for AA files
        belonging to sname (matching [sname]_amplicon*, [sname]_summary.txt).
        Copy matches into staging/.  Returns True if anything was copied.
        """
        found = False
        for files_dir in self._files_dirs.values():
            if not os.path.isdir(files_dir):
                continue
            try:
                for fname in os.listdir(files_dir):
                    if not (fname.startswith(sname + "_amplicon")
                            or fname == f"{sname}_summary.txt"):
                        continue
                    src = os.path.join(files_dir, fname)
                    if os.path.isfile(src):
                        shutil.copy2(src, staging)
                        found = True
            except OSError as e:
                print(f"  Warning: could not read files dir {files_dir}: {e}")
        return found

    def _build_cnvkit_tarball(self, sname: str, rec: SampleRecord,
                               sample_out: str) -> Optional[str]:
        """
        Build [sname]_cnvkit_output.tar.gz in sample_out.

        If rec.cnvkit_dir exists, tar it (minus excluded suffixes).
        Otherwise, return None — cnvkit directory will be Not Provided.
        """
        tar_name = f"{sname}_cnvkit_output.tar.gz"
        tar_dest = os.path.join(sample_out, tar_name)

        if rec.cnvkit_dir and os.path.isdir(rec.cnvkit_dir):
            # Compress .cns files in place before tarring, consistent with
            # the old aggregator behaviour.
            self._gzip_files_in_dir(rec.cnvkit_dir, ".cns")
            make_tarball(rec.cnvkit_dir, tar_dest)
            return tar_dest

        print(f"  Warning: no cnvkit dir found for '{sname}' — cnvkit directory: Not Provided")
        return None

    def _copy_cnv_bed(self, sname: str, rec: SampleRecord,
                      sample_out: str) -> Optional[str]:
        """
        Copy the _CNV_CALLS.bed file uncompressed into sample_out.
        This uncompressed copy is what CNV BED file in run.json points to.
        Returns the destination path, or None if no BED was found.
        """
        if not rec.cnv_calls_bed or not os.path.isfile(rec.cnv_calls_bed):
            return None

        dest = os.path.join(sample_out, f"{sname}_CNV_CALLS.bed")
        safe_copy_file(rec.cnv_calls_bed, dest)
        return dest

    @staticmethod
    def _gzip_files_in_dir(dirpath: str, suffix: str) -> None:
        """gzip all files with the given suffix inside dirpath (in-place)."""
        import subprocess
        try:
            for fname in os.listdir(dirpath):
                if fname.endswith(suffix) and not fname.endswith(suffix + ".gz"):
                    fpath = os.path.join(dirpath, fname)
                    subprocess.call(["gzip", "-fq", fpath])
        except OSError:
            pass

    # ------------------------------------------------------------------
    # Stage 5b — consolidated classification
    # ------------------------------------------------------------------

    def _build_consolidated_classification(self) -> None:
        """
        Merge all AC output files from all classification dirs into
        results/consolidated_classification/ using [project_name] as prefix.
        """
        print(f"  Building consolidated classification for project '{self.project_name}'...")

        if not self.classification_dirs:
            print("  Warning: no classification dirs found — "
                  "consolidated_classification/ will be empty.")
            return

        for target in AC_MERGE_TARGETS:
            suffix   = target["suffix"]
            has_hdr  = target["has_header"]
            is_dir   = target["is_dir"]
            out_name = f"{self.project_name}{suffix}"

            if is_dir:
                out_path = os.path.join(self.classif_dir, out_name)
                os.makedirs(out_path, exist_ok=True)
                self._merge_ac_subdirs(suffix, out_path)
            else:
                out_path = os.path.join(self.classif_dir, out_name)
                self._merge_ac_tsvs(suffix, out_path, has_header=has_hdr)

    def _merge_ac_tsvs(self, suffix: str, out_path: str,
                        has_header: Optional[bool]) -> None:
        """
        Concatenate all TSV files with the given suffix found across all
        classification dirs into out_path.

        - has_header=True:  keep first file's header; drop subsequent headers;
                            warn if headers differ across files.
        - has_header=False: straight concatenation, no header logic.
        - has_header=None:  file is optional; skip silently if none found.
        """
        sources = self._find_ac_files(suffix)

        if not sources:
            if has_header is not None:
                # True = expected; warn only for expected files
                if has_header is True:
                    print(f"  Warning: no files found for '{suffix}' — "
                          f"output file will not be created.")
            return

        written_header: Optional[str] = None
        lines_written = 0

        with open(out_path, "w") as out_fh:
            for src_path in sources:
                try:
                    with open(src_path, "r") as in_fh:
                        raw_lines = in_fh.readlines()
                except OSError as e:
                    print(f"  Warning: could not read {src_path}: {e}")
                    continue

                if not raw_lines:
                    continue

                if has_header:
                    file_header = raw_lines[0]
                    data_lines  = raw_lines[1:]
                    if written_header is None:
                        written_header = file_header
                        out_fh.write(file_header)
                    elif file_header.strip() != written_header.strip():
                        print(f"  Warning: header mismatch in '{src_path}' "
                              f"for suffix '{suffix}' — using first header, "
                              f"skipping this file's header line.")
                    for line in data_lines:
                        out_fh.write(line)
                        lines_written += 1
                else:
                    # No header — straight concat
                    for line in raw_lines:
                        out_fh.write(line)
                        lines_written += 1

        print(f"  Merged {len(sources)} file(s) -> {os.path.basename(out_path)} "
              f"({lines_written} data line(s))")

    def _merge_ac_subdirs(self, suffix: str, out_dir: str) -> None:
        """
        Copy all files from every subdirectory matching suffix across all
        classification dirs into the flat output directory out_dir.
        Filename collisions are handled by appending a numeric suffix.
        """
        sources = self._find_ac_dirs(suffix)
        if not sources:
            return

        files_copied = 0
        for src_dir in sources:
            try:
                for fname in os.listdir(src_dir):
                    src_file = os.path.join(src_dir, fname)
                    if not os.path.isfile(src_file):
                        continue
                    dest_file = self._unique_dest(out_dir, fname)
                    shutil.copy2(src_file, dest_file)
                    files_copied += 1
            except OSError as e:
                print(f"  Warning: could not read subdir {src_dir}: {e}")

        print(f"  Merged {len(sources)} dir(s) -> {os.path.basename(out_dir)}/ "
              f"({files_copied} file(s))")

    def _find_ac_files(self, suffix: str) -> List[str]:
        """
        Return sorted list of file paths ending with suffix found directly
        inside any known classification dir (not recursing into subdirs).
        Also searches the full classif_dir tree to catch files at any depth.
        """
        found = []
        seen = set()
        for cls_dir in self.classification_dirs:
            try:
                for fname in os.listdir(cls_dir):
                    if fname.endswith(suffix) and not fname.startswith("._"):
                        fpath = os.path.join(cls_dir, fname)
                        if os.path.isfile(fpath) and fpath not in seen:
                            found.append(fpath)
                            seen.add(fpath)
            except OSError:
                continue
        return sorted(found)

    def _find_ac_dirs(self, suffix: str) -> List[str]:
        """
        Return sorted list of directory paths ending with suffix found directly
        inside any known classification dir.
        """
        found = []
        seen = set()
        for cls_dir in self.classification_dirs:
            try:
                for dname in os.listdir(cls_dir):
                    if dname.endswith(suffix) and not dname.startswith("._"):
                        dpath = os.path.join(cls_dir, dname)
                        if os.path.isdir(dpath) and dpath not in seen:
                            found.append(dpath)
                            seen.add(dpath)
            except OSError:
                continue
        return sorted(found)

    # ------------------------------------------------------------------
    # Stage 5c — AUX dirs
    # ------------------------------------------------------------------

    def _copy_aux_dirs(self) -> None:
        """Copy each AUX_DIR-marked directory wholesale into other_files/."""
        for aux_dir in self.aux_dirs:
            dname = os.path.basename(aux_dir)
            dest  = self._unique_dest(self.other_dir, dname)
            print(f"  Copying AUX dir: {aux_dir} -> {dest}")
            safe_copytree(aux_dir, dest)

    # ==================================================================
    # Stage stubs — implemented in subsequent segments
    # ==================================================================


    # ==================================================================
    # Stage 6 — run.json constructor
    # ==================================================================

    def _stage6_build_run_json(self) -> None:
        """
        Build results/run.json and results/aggregated_results.csv.

        For each feature row in self.run_json_groups:
          - Parse list-valued fields (Location, Oncogenes, All genes)
          - Re-resolve all stale file paths to new locations in results/
          - Apply sample name remapping if configured
          - Write the structured run.json
          - Write a flat aggregated_results.csv alongside
        """
        print("\n--- Stage 6: Building run.json ---")

        ref_genomes: set = set()
        processed = 0

        for skey, rows in self.run_json_groups.items():
            for row in rows:
                processed += 1
                sname = str(row["Sample name"])  # guard against numeric sample names
                row["Sample name"] = sname
                rec = self.sample_registry.get(sname)

                # ── Parse list-valued fields ──────────────────────────────
                for col in LIST_COLUMNS:
                    if col in row:
                        row[col] = parse_list_field(row[col])
                    else:
                        row[col] = []

                # ── Reference genome consistency check ───────────────────
                ref = row.get("Reference version", NOT_PROVIDED)
                if not not_provided(ref):
                    ref_genomes.add(ref)
                    if len(ref_genomes) > 1:
                        self._abort(
                            f"Multiple reference genomes detected: {ref_genomes}. "
                            "AmpliconRepository only supports single-reference projects."
                        )

                # ── Amplicon number (integer coercion) ────────────────────
                try:
                    row["AA amplicon number"] = int(row["AA amplicon number"])
                except (ValueError, TypeError, KeyError):
                    pass

                # ── Numeric field coercions ───────────────────────────────
                for num_col in ("Complexity score", "Captured interval length",
                                "Feature median copy number",
                                "Feature maximum copy number"):
                    try:
                        row[num_col] = float(row[num_col])
                    except (ValueError, TypeError, KeyError):
                        pass

                # ── Path re-resolution ────────────────────────────────────
                self._resolve_paths(row, rec)

                # ── Sample name remapping ─────────────────────────────────
                if sname in self.name_map:
                    row["Sample name"] = self.name_map[sname]
                elif self.name_map:
                    print(f"  Warning: sample '{sname}' not found in name_map.")

                # ── Ensure canonical column order ─────────────────────────
                ordered = {}
                for col in RUN_JSON_COLUMNS:
                    ordered[col] = row.get(col, NOT_PROVIDED)
                self.run_json_groups[skey][
                    self.run_json_groups[skey].index(row)] = ordered

                if processed % 100 == 0:
                    print(f"  Processed {processed} feature rows...")

        # ── Write run.json ────────────────────────────────────────────────
        run_json_path = os.path.join(self.results_dir, "run.json")
        with open(run_json_path, "w") as fh:
            json.dump({"runs": self.run_json_groups}, fh, indent=2, sort_keys=True)
        print(f"  Wrote run.json  ({processed} feature rows, "
              f"{len(self.run_json_groups)} sample key(s))")

        # ── Write aggregated_results.csv ──────────────────────────────────
        # List fields rendered as Python repr strings e.g. ['EGFR', 'MYC'].
        # Columns are AGG_CSV_COLUMNS subset in spec order (no AA/cnvkit dir).
        # Rows sorted by: Sample name, AA amplicon number, Feature ID.
        import csv as _csv
        all_flat_rows = [row for rows in self.run_json_groups.values() for row in rows]
        all_flat_rows.sort(key=lambda r: (
            str(r.get("Sample name", "")),
            int(r.get("AA amplicon number", 0)) if str(r.get("AA amplicon number", "")).isdigit() else 0,
            str(r.get("Feature ID", "")),
        ))
        csv_path = os.path.join(self.results_dir, "aggregated_results.csv")
        with open(csv_path, "w", newline="") as fh:
            writer = _csv.writer(fh)
            writer.writerow(AGG_CSV_COLUMNS)
            for row in all_flat_rows:
                cells = []
                for col in AGG_CSV_COLUMNS:
                    v = row.get(col, NOT_PROVIDED)
                    cells.append(str(v) if isinstance(v, list) else
                                 ("" if v is None else str(v)))
                writer.writerow(cells)
        print(f"  Wrote aggregated_results.csv")

    # ------------------------------------------------------------------
    # Stage 6 helpers
    # ------------------------------------------------------------------

    def _resolve_paths(self, row: dict, rec: Optional[SampleRecord]) -> None:
        """
        Re-resolve all stale file path fields in a feature row dict.
        Modifies row in-place.  All resolved paths are expressed relative
        to self.results_dir.
        """
        sname = row.get("Sample name", "")
        amp_num_raw = row.get("AA amplicon number")
        try:
            amp_num = int(amp_num_raw)
        except (ValueError, TypeError):
            amp_num = None

        # ── Feature BED file ─────────────────────────────────────────────
        # Resolve from consolidated_classification bed files dir.
        feat_bed_basename = os.path.basename(
            str(row.get("Feature BED file", "")))
        bed_files_dir = os.path.join(
            self.classif_dir,
            f"{self.project_name}_classification_bed_files")
        row["Feature BED file"] = self._resolve_from_dir(
            feat_bed_basename, bed_files_dir)

        # ── CNV BED file ─────────────────────────────────────────────────
        # Always use the uncompressed copy placed by Stage 5.
        if rec and rec.cnv_bed_dest and os.path.isfile(rec.cnv_bed_dest):
            row["CNV BED file"] = relative_to_results(
                rec.cnv_bed_dest, self.results_dir)
        else:
            row["CNV BED file"] = NOT_PROVIDED

        # ── AA PNG file ──────────────────────────────────────────────────
        row["AA PNG file"] = self._resolve_amplicon_file(
            rec, amp_num, "png", sname)

        # ── AA PDF file ──────────────────────────────────────────────────
        row["AA PDF file"] = self._resolve_amplicon_file(
            rec, amp_num, "pdf", sname)

        # ── AA cycles file ────────────────────────────────────────────────
        row["AA cycles file"] = self._resolve_amplicon_file(
            rec, amp_num, "cycles", sname)

        # ── AA graph file ─────────────────────────────────────────────────
        row["AA graph file"] = self._resolve_amplicon_file(
            rec, amp_num, "graph", sname)

        # ── Run metadata JSON ────────────────────────────────────────────
        row["Run metadata JSON"] = self._resolve_misc_path(
            rec, "run_metadata_json", sname, "_run_metadata.json")

        # ── Sample metadata JSON ─────────────────────────────────────────
        row["Sample metadata JSON"] = self._resolve_misc_path(
            rec, "sample_metadata_json", sname, "_sample_metadata.json")

        # ── AA directory ─────────────────────────────────────────────────
        if rec and rec.aa_dir_dest and os.path.isdir(rec.aa_dir_dest):
            row["AA directory"] = relative_to_results(
                rec.aa_dir_dest, self.results_dir)
        else:
            row["AA directory"] = NOT_PROVIDED

        # ── cnvkit directory ─────────────────────────────────────────────
        if rec and rec.cnvkit_tarball and os.path.isfile(rec.cnvkit_tarball):
            row["cnvkit directory"] = relative_to_results(
                rec.cnvkit_tarball, self.results_dir)
        else:
            row["cnvkit directory"] = NOT_PROVIDED

    def _resolve_from_dir(self, basename: str, search_dir: str) -> str:
        """
        Look for a file with the given basename inside search_dir.
        Returns a results-relative path, or NOT_PROVIDED.
        """
        if not basename or not_provided(basename):
            return NOT_PROVIDED
        candidate = os.path.join(search_dir, basename)
        if os.path.isfile(candidate):
            return relative_to_results(candidate, self.results_dir)
        return NOT_PROVIDED

    def _resolve_amplicon_file(self, rec: Optional[SampleRecord],
                                amp_num: Optional[int],
                                key: str, sname: str) -> str:
        """
        Resolve a per-amplicon file (pdf, png, cycles, graph) for a given amplicon number.
        Looks for the file inside rec.aa_dir_dest (the uncompressed AA results dir).
        Returns a results-relative path, or NOT_PROVIDED.
        """
        if rec is None or amp_num is None or not rec.aa_dir_dest:
            return NOT_PROVIDED
        ext_map = {"pdf": ".pdf", "png": ".png",
                   "cycles": "_cycles.txt", "graph": "_graph.txt"}
        ext = ext_map.get(key)
        if not ext:
            return NOT_PROVIDED
        path = os.path.join(rec.aa_dir_dest, f"{sname}_amplicon{amp_num}{ext}")
        if os.path.isfile(path):
            return relative_to_results(path, self.results_dir)
        return NOT_PROVIDED

    def _resolve_misc_path(self, rec: Optional[SampleRecord],
                            attr: str, sname: str, suffix: str) -> str:
        """
        Resolve a miscellaneous file (metadata JSON, etc.) that was copied
        placed into the sample dir by Stage 5.
        Falls back to checking the sample dir directly by expected filename.
        Returns a results-relative path, or NOT_PROVIDED.
        """
        if rec is None:
            return NOT_PROVIDED

        # Primary: use the path recorded on the SampleRecord
        src = getattr(rec, attr, None)
        if src and os.path.isfile(src):
            # File lives in the extraction tree; check if there's a copy in
            # the sample dir (placed by Stage 5 misc copy loop).
            sample_dir = os.path.join(self.samples_dir, sname)
            dest = os.path.join(sample_dir, f"{sname}{suffix}")
            if os.path.isfile(dest):
                return relative_to_results(dest, self.results_dir)
            # Fall back to the extraction tree path (unusual but safe)
            return relative_to_results(src, self.results_dir)

        return NOT_PROVIDED

    # ==================================================================
    # Finalise — create output archive and clean up
    # ==================================================================

    def _finalise(self) -> None:
        """
        Tar the entire results/ tree into [project_name].tar.gz and
        remove working directories (unless --no_cleanup).
        """
        print("\n--- Finalise: Creating output archive ---")
        output_archive = os.path.join(self.work_dir,
                                      f"{self.project_name}.tar.gz")
        print(f"  Writing: {output_archive}")
        make_tarball(self.results_dir, output_archive)
        output_bytes = os.path.getsize(output_archive)
        output_mb = output_bytes / (1024 * 1024)
        input_bytes = getattr(self, "_input_size_bytes", 0)
        print(f"  Output archive: {output_mb:.2f} MB")
        if input_bytes > 0:
            input_mb = input_bytes / (1024 * 1024)
            ratio = output_bytes / input_bytes
            print(f"  Size change: {input_mb:.2f} MB -> {output_mb:.2f} MB "
                  f"({ratio:.2f}x, {(1 - ratio) * 100:.1f}% reduction)")
        self._cleanup(failure=False)
        elapsed = time.perf_counter() - self._start_time
        minutes, seconds = divmod(elapsed, 60)
        if minutes >= 1:
            print(f"  Total time: {int(minutes)}m {seconds:.1f}s")
        else:
            print(f"  Total time: {elapsed:.1f}s")
        print(f"  Done.")