"""Synthetic featureless fixture, adapted from the local generator."""

import json
import os

AA_VERSION = "1.5.r5"
ASP_VERSION = "1.5.0"
AC_VERSION = "1.5.1"
REF = "GRCh38"

# Column order as AC's make_results_table.py emits it.
RT_COLUMNS = [
    "Sample name", "AA amplicon number", "Feature ID", "Classification",
    "Location", "Oncogenes", "All genes", "NCBI Gene IDs",
    "Complexity score", "ecDNA context", "Captured interval length",
    "Feature median copy number", "Feature maximum copy number",
    "Filter flag", "Reference version", "Tissue of origin", "Sample type",
    "Feature BED file", "CNV BED file",
    "AS-p version", "AA version", "AC version",
    "AA PNG file", "AA PDF file", "AA summary file",
    "Run metadata JSON", "Sample metadata JSON",
]

PROFILES_COLUMNS = [
    "sample_name", "amplicon_number", "amplicon_decomposition_class",
    "ecDNA+", "BFB+", "ecDNA_amplicons",
]


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def tsv(columns, rows):
    out = ["\t".join(columns)]
    for row in rows:
        out.append("\t".join(str(row.get(c, "NA")) for c in columns))
    return "\n".join(out) + "\n"


def aa_summary(sname, n_amplicons):
    """AA's _summary.txt. First line is what is_aa_summary_content() keys on."""
    lines = [f"#Amplicons = {n_amplicons}",
             "-" * 89]
    for i in range(1, n_amplicons + 1):
        lines += [
            f"[amplicon{i}] AmpliconID = {i}",
            f"[amplicon{i}] #Intervals = 2",
            f"[amplicon{i}] Intervals = chr7:54000000-55000000,chr7:55500000-56000000",
            f"[amplicon{i}] OncogenesAmplified = EGFR",
            "-" * 89,
        ]
    return "\n".join(lines) + "\n"


def graph_txt():
    return (
        "SequenceEdge: StartPosition, EndPosition, PredecessorPosition, "
        "SuccessorPosition, Coverage, Size, NumberReadsMapped\n"
        "sequence\tchr7:54000000-\tchr7:55000000+\t8.5\t1000000\t4200\n"
    )


def cycles_txt():
    return (
        "Interval\t1\tchr7\t54000000\t55000000\n"
        "Segment\t1\tchr7\t54000000\t55000000\n"
        "Cycle=1;Copy_count=8.5;Segments=1+,0+\n"
    )


def metadata_jsons(sample_dir, sname):
    write(os.path.join(sample_dir, f"{sname}_run_metadata.json"),
          json.dumps({
              "AA_version": AA_VERSION,
              "AS-p_version": ASP_VERSION,
              "AC_version": AC_VERSION,
              "ref_genome": REF,
          }, indent=2) + "\n")
    write(os.path.join(sample_dir, f"{sname}_sample_metadata.json"),
          json.dumps({"sample_type": "", "tissue_of_origin": ""}, indent=2) + "\n")


def base_row(sname):
    """Sample-level fields shared by every result_table row."""
    return {
        "Sample name": sname,
        "Reference version": REF,
        "AS-p version": ASP_VERSION,
        "AA version": AA_VERSION,
        "AC version": AC_VERSION,
        "Tissue of origin": "",
        "Sample type": "",
        "AA summary file": f"{sname}_AA_results/{sname}_summary.txt",
        "Run metadata JSON": f"{sname}_run_metadata.json",
        "Sample metadata JSON": f"{sname}_sample_metadata.json",
    }


def build_sample(root, sname, n_amplicons, rt_rows, profile_rows,
                 copy_summary_to_files):
    """
    Lay out one sample in AmpliconSuite-pipeline's on-disk shape.

    copy_summary_to_files mirrors AC's real behaviour: it populates
    _classification/files/ from the features it kept, so a sample with no
    features ends up with only the metadata JSONs there.
    """
    sample_dir = os.path.join(root, sname)
    aa_dir = os.path.join(sample_dir, f"{sname}_AA_results")
    cls_dir = os.path.join(sample_dir, f"{sname}_classification")
    files_dir = os.path.join(cls_dir, "files")

    # --- AA results ---
    write(os.path.join(aa_dir, f"{sname}_summary.txt"),
          aa_summary(sname, n_amplicons))
    write(os.path.join(aa_dir, f"{sname}.log"),
          f"AmpliconArchitect version {AA_VERSION}\nRunning AA on {sname}\n")
    for i in range(1, n_amplicons + 1):
        write(os.path.join(aa_dir, f"{sname}_amplicon{i}_graph.txt"), graph_txt())
        write(os.path.join(aa_dir, f"{sname}_amplicon{i}_cycles.txt"), cycles_txt())

    # --- pipeline-level files ---
    write(os.path.join(sample_dir, f"{sname}.log"),
          f"AmpliconSuite-pipeline version {ASP_VERSION}\n"
          f"AmpliconArchitect version {AA_VERSION}\n")
    write(os.path.join(sample_dir, f"{sname}_finish_flag.txt"), "All stages completed\n")
    metadata_jsons(sample_dir, sname)

    # --- AC classification ---
    write(os.path.join(cls_dir, f"{sname}_result_table.tsv"),
          tsv(RT_COLUMNS, rt_rows))
    write(os.path.join(cls_dir, f"{sname}_amplicon_classification_profiles.tsv"),
          tsv(PROFILES_COLUMNS, profile_rows))
    write(os.path.join(cls_dir, f"{sname}.log"),
          f"AmpliconClassifier {AC_VERSION}\n"
          f"Classification stage completed\n")

    # AC's copy dir — only ever holds what AC had features for.
    metadata_jsons(files_dir, sname)
    if copy_summary_to_files:
        write(os.path.join(files_dir, f"{sname}_summary.txt"),
              aa_summary(sname, n_amplicons))


def build_project(root):
    """Write synthetic inputs under the caller-provided temporary directory."""
    # 1. Ordinary sample: AA found an amplicon, AC kept it.
    focal = dict(base_row("SYNTH_focal"))
    focal.update({
        "AA amplicon number": "1",
        "Feature ID": "SYNTH_focal_amplicon1_ecDNA_1",
        "Classification": "ecDNA",
        "Location": "['chr7:54000000-55000000']",
        "Oncogenes": "['EGFR']",
        "All genes": "['EGFR']",
        "NCBI Gene IDs": "['1956']",
        "Complexity score": "0.61",
        "ecDNA context": "simple",
        "Captured interval length": "1000000",
        "Feature median copy number": "8.5",
        "Feature maximum copy number": "9.2",
        "Filter flag": "PASS",
        "Feature BED file": "SYNTH_focal_amplicon1_ecDNA_1_intervals.bed",
    })
    build_sample(
        root, "SYNTH_focal", n_amplicons=1,
        rt_rows=[focal],
        profile_rows=[{
            "sample_name": "SYNTH_focal", "amplicon_number": "amplicon1",
            "amplicon_decomposition_class": "Cyclic", "ecDNA+": "Positive",
            "BFB+": "None", "ecDNA_amplicons": "1",
        }],
        copy_summary_to_files=True,
    )

    # 2. AA found nothing at all — AC writes its own "NA" placeholder row.
    #    Worked before 8.2.0; present here to prove the fix does not disturb it.
    noamp = dict(base_row("SYNTH_noamp"))
    noamp.update({
        "AA amplicon number": "NA",
        "Feature ID": "SYNTH_noamp_NA",
        "Classification": "NA",
        "Location": "[]", "Oncogenes": "[]", "All genes": "[]",
        "NCBI Gene IDs": "[]",
    })
    build_sample(
        root, "SYNTH_noamp", n_amplicons=0,
        rt_rows=[noamp],
        profile_rows=[],
        copy_summary_to_files=True,
    )

    # 3. THE REGRESSION CASE: AA found 2 amplicons, AC classified none as
    #    focal and emitted a header-only result_table.tsv. No placeholder row,
    #    so before 8.2.0 the sample vanished from the project.
    build_sample(
        root, "SYNTH_filtered", n_amplicons=2,
        rt_rows=[],
        profile_rows=[],
        copy_summary_to_files=False,
    )

