"""Check copy-number recovery from synthetic CNVkit segments."""

import pytest

from asa_aggregator import convert_cnvkit_cns_to_bed


@pytest.mark.parametrize("minimum,expected", [
    (0, "chr1\t0\t100\tsize=100\t1.0\n"
        "chr1\t100\t300\tsize=200\t2.0\n"
        "chr2\t20\t70\tsize=50\t8.0\n"),
    (2, "chr1\t100\t300\tsize=200\t2.0\n"
        "chr2\t20\t70\tsize=50\t8.0\n"),
    (9, ""),
])
def test_cns_conversion_preserves_intervals_and_filters_copy_number(tmp_path, minimum, expected):
    source = tmp_path / "sample.cns"
    source.write_text(
        "chromosome\tstart\tend\tgene\tlog2\tprobes\n"
        "chr1\t0\t100\t-\t-1\t10\n"
        "chr1\t100\t300\t-\t0\t20\n"
        "chr2\t20\t70\tEGFR\t2\t5\n"
    )
    destination = tmp_path / "calls.bed"
    convert_cnvkit_cns_to_bed(str(source), str(destination), min_cn=minimum)
    assert destination.read_text() == expected
