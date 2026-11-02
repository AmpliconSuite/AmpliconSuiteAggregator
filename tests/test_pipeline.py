"""Exercise the CLI with generated inputs, including complete archive round trips."""

import csv
import io
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import zipfile

import pytest

from fixtures import featureless, prefix_collision


CLI = Path(__file__).resolve().parents[1] / "src" / "AmpliconSuiteAggregator.py"


def run_cli(work, source, *options, success=True):
    work.mkdir()
    result = subprocess.run(
        [sys.executable, str(CLI), "--files", str(source), "-o", "project", *options],
        cwd=work, capture_output=True, text=True,
    )
    if not success:
        assert result.returncode != 0, result.stdout + result.stderr
        assert not (work / "project.tar.gz").exists()
        return result
    assert result.returncode == 0, result.stdout + result.stderr
    archive = work / "project.tar.gz"
    assert archive.is_file(), result.stdout + result.stderr
    return archive


def members(archive):
    with tarfile.open(archive) as tar:
        return {m.name: tar.extractfile(m).read() for m in tar if m.isfile()}


def runs(contents):
    return json.loads(contents["results/run.json"])["runs"]


def rows_by_sample(contents):
    return {rows[0]["Sample name"]: rows for rows in runs(contents).values()}


def assert_same_archive(first, second):
    # Compare contents first to identify which file caused a round-trip failure.
    assert members(first) == members(second)
    assert first.read_bytes() == second.read_bytes()


@pytest.mark.parametrize("rename", [False, True])
def test_round_trip_uses_final_sample_names(tmp_path, rename):
    source = tmp_path / "input"
    featureless.build_project(source)
    mapping = {
        "SYNTH_filtered": "Zulu",
        "SYNTH_focal": "Middle",
        "SYNTH_noamp": "Alpha",
    }
    options = []
    if rename:
        name_map = tmp_path / "names.tsv"
        name_map.write_text("".join(f"{old}\t{new}\n" for old, new in mapping.items()))
        options = ["--name_map", str(name_map)]
    first = run_cli(tmp_path / "first", source, *options)
    repeated = run_cli(tmp_path / "repeat", source, *options)
    # Renaming is a transformation of the original names; do not apply it twice.
    second = run_cli(tmp_path / "reaggregate", first)
    assert_same_archive(first, repeated)
    assert_same_archive(first, second)
    expected = sorted(mapping.values() if rename else mapping)
    assert {key: rows[0]["Sample name"] for key, rows in runs(members(first)).items()} == {
        f"sample_{i}": name for i, name in enumerate(expected, 1)
    }


@pytest.mark.parametrize("input_format", ["directory", "tar", "zip"])
def test_featureless_samples_retain_reconstruction_and_metadata(tmp_path, input_format):
    source = tmp_path / "input"
    featureless.build_project(source)
    if input_format == "tar":
        packed = tmp_path / "input.tar.gz"
        with tarfile.open(packed, "w:gz") as tar:
            tar.add(source, arcname="input")
        source = packed
    elif input_format == "zip":
        packed = tmp_path / "input.zip"
        with zipfile.ZipFile(packed, "w") as zf:
            for path in sorted(source.rglob("*")):
                if path.is_file():
                    zf.write(path, str(Path("input") / path.relative_to(source)))
        source = packed
    contents = members(run_cli(tmp_path / "output", source))
    rows = rows_by_sample(contents)
    assert set(rows) == {"SYNTH_focal", "SYNTH_noamp", "SYNTH_filtered"}
    assert rows["SYNTH_focal"][0]["Classification"] == "ecDNA"
    for name in ("SYNTH_noamp", "SYNTH_filtered"):
        assert len(rows[name]) == 1
        assert rows[name][0]["Classification"] == "NA"
        assert rows[name][0]["AmpliconArchitect version"] == "1.5.r5"
        assert rows[name][0]["AmpliconClassifier version"] == "1.5.1"
        for field in ("Run metadata JSON", "Sample metadata JSON"):
            assert "results/" + rows[name][0][field] in contents
    prefix = "results/samples/SYNTH_filtered/SYNTH_filtered_reconstruction_results/"
    for amp in (1, 2):
        for kind, expected in (("graph", featureless.graph_txt()),
                               ("cycles", featureless.cycles_txt())):
            assert contents[f"{prefix}SYNTH_filtered_amplicon{amp}_{kind}.txt"].decode() == expected
    csv_rows = list(csv.DictReader(io.StringIO(contents["results/aggregated_results.csv"].decode())))
    assert {row["Sample name"] for row in csv_rows} == set(rows)


@pytest.mark.parametrize("excluded,retained", [
    ("HCC827", "HCC827_res"), ("HCC827_res", "HCC827"),
])
def test_exclusion_preserves_sample_with_shared_prefix(tmp_path, excluded, retained):
    source = tmp_path / "input"
    prefix_collision.build_project(source)
    archive = run_cli(tmp_path / "output", source, "--exclude_samples", excluded)
    contents = members(archive)
    assert set(rows_by_sample(contents)) == {retained}
    assert rows_by_sample(contents)[retained][0]["Classification"] == "ecDNA"
    assert not any(name.startswith(f"results/samples/{excluded}/") for name in contents)
    for name, data in contents.items():
        assert f"{excluded}_amplicon" not in name
        if name.endswith((".tsv", ".csv", ".json")):
            assert f"{excluded}_amplicon".encode() not in data
    cls = "results/consolidated_classification/project"
    for suffix, filename in (
        ("classification_bed_files", f"{retained}_amplicon1_ecDNA_1_intervals.bed"),
        ("SV_summaries", f"{retained}_amplicon1_SV_summary.tsv"),
        ("annotated_cycles_files", f"{retained}_amplicon1_annotated_cycles.txt"),
    ):
        assert f"{cls}_{suffix}/{filename}" in contents
    assert contents[f"{cls}_ecDNA_context_calls.tsv"].decode() == (
        f"{retained}_amplicon1_ecDNA_1\tSimple circular\n"
    )
    for suffix, column in (
        ("result_table", "Sample name"),
        ("amplicon_classification_profiles", "sample_name"),
        ("gene_list", "sample_name"),
        ("ecDNA_counts", "#sample"),
    ):
        rows = csv.DictReader(io.StringIO(contents[f"{cls}_{suffix}.tsv"].decode()), delimiter="\t")
        assert {row[column] for row in rows} == {retained}
    assert_same_archive(archive, run_cli(tmp_path / "reaggregate", archive))


def test_excluding_all_samples_fails_without_output(tmp_path):
    source = tmp_path / "input"
    prefix_collision.build_project(source)
    result = run_cli(tmp_path / "output", source, "--exclude_samples",
                     "HCC827", "HCC827_res", success=False)
    assert "Every sample in the input was excluded" in result.stdout + result.stderr
