#!/usr/bin/env python3
"""
AmpliconSuiteAggregator.py
Entry point for the AmpliconSuite aggregation pipeline.

Takes one or more AmpliconSuite-pipeline output archives (or directories)
and produces a single structured archive ready for upload to AmpliconRepository.org.

Authors: Jens Luebeck et al.
"""

import argparse
import sys
import os

from asa_stages import Aggregator
from asa_aggregator import __version__


def get_paths_from_filelist(filelist_fp: str) -> list[str]:
    """Read input paths from a text file (one path per line, blank lines ignored)."""
    paths = []
    try:
        with open(filelist_fp) as fh:
            for line in fh:
                p = line.strip()
                if p:
                    paths.append(p)
    except OSError as e:
        sys.stderr.write(f"Error reading filelist '{filelist_fp}': {e}\n")
        sys.exit(1)

    if not paths:
        sys.stderr.write(f"Filelist '{filelist_fp}' is empty or contains no valid paths.\n")
        sys.exit(1)

    return paths


def validate_inputs(paths: list[str]) -> None:
    """Verify each input path exists and is a recognised type."""
    valid_extensions = (".tar.gz", ".tar", ".zip")
    errors = []
    for p in paths:
        if not os.path.exists(p):
            errors.append(f"  Input not found: {p}")
        elif os.path.isfile(p) and not any(p.endswith(ext) for ext in valid_extensions):
            errors.append(f"  Unrecognised file type (expected .tar.gz / .tar / .zip): {p}")

    if errors:
        sys.stderr.write("Input validation failed:\n" + "\n".join(errors) + "\n")
        sys.exit(1)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="AmpliconSuiteAggregator",
        description="Aggregate AmpliconSuite-pipeline outputs into a structured archive "
                    "for upload to AmpliconRepository.org.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # --- Input (mutually exclusive) ---
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "-flist", "--filelist",
        metavar="FILE",
        type=str,
        help="Text file listing input paths to aggregate (one per line). "
             "Each path may be a .tar.gz, .tar, .zip archive, or an uncompressed directory.",
    )
    input_group.add_argument(
        "--files",
        nargs="+",
        metavar="PATH",
        type=str,
        help="One or more input paths to aggregate directly "
             "(.tar.gz / .tar / .zip archives or uncompressed directories).",
    )

    # --- Output / project ---
    parser.add_argument(
        "-o", "--output_name",
        metavar="NAME",
        type=str,
        default="output",
        help="Output filename prefix and project name used for consolidated classification files. "
             "The final archive will be written to <NAME>.tar.gz. (default: 'output')",
    )

    # --- Optional features ---
    parser.add_argument(
        "--name_map",
        metavar="FILE",
        type=str,
        default=None,
        help="Two-column TSV/space-delimited file mapping current sample identifiers (col 1) "
             "to replacement names (col 2). Enables batch renaming of samples in run.json.",
    )
    parser.add_argument(
        "--no_cleanup",
        action="store_true",
        default=False,
        help="Skip removal of temporary working directories on completion. "
             "Useful for debugging.",
    )

    # --- Version ---
    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    # Resolve input paths
    if args.filelist:
        input_paths = get_paths_from_filelist(args.filelist)
    else:
        input_paths = args.files

    validate_inputs(input_paths)

    print(f"AmpliconSuiteAggregator {__version__}")
    print(f"Project name  : {args.output_name}")
    print(f"Input paths   : {len(input_paths)} item(s)")
    for p in input_paths:
        print(f"  {p}")
    print(f"Output archive: {args.output_name}.tar.gz")
    print(f"No-cleanup    : {args.no_cleanup}")
    if args.name_map:
        print(f"Name map      : {args.name_map}")
    print()

    aggregator = Aggregator(
        input_paths=input_paths,
        project_name=args.output_name,
        name_map_file=args.name_map,
        no_cleanup=args.no_cleanup,
    )

    if not aggregator.completed:
        sys.stderr.write(
            "Aggregation did not complete successfully. "
            "Check warnings above for details.\n"
        )
        sys.exit(1)

    print(f"\nDone. Output written to: {args.output_name}.tar.gz")


if __name__ == "__main__":
    main()