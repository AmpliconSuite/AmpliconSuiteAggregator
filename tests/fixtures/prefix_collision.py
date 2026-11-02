"""Synthetic prefix collision fixture, adapted from the local generator."""

import json
import os

AA_VERSION = "1.5.r5"
ASP_VERSION = "1.5.0"
AC_VERSION = "1.5.1"
REF = "GRCh38"

COHORT = "cohort"

# The two names are the whole point: one is a strict prefix of the other.
SAMPLES = ["HCC827", "HCC827_res"]

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

GENE_LIST_COLUMNS = [
    "sample_name", "amplicon_number", "feature", "gene", "gene_cn",
    "truncated", "is_canonical_oncogene", "ncbi_id",
]

FEATURE_PROPS_COLUMNS = [
    "feature_ID", "captured_region_size_bp", "median_feature_CN",
    "max_feature_CN", "borderline_flag",
]

ECDNA_COUNTS_COLUMNS = ["#sample", "ecDNA_count"]

SV_SUMMARY_COLUMNS = [
    "chrom1", "pos1", "chrom2", "pos2", "sv_type", "read_support",
    "features", "orientation", "pos1_flanking_coordinate",
    "pos2_flanking_coordinate", "homology_length", "homology_sequence",
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


def aa_summary(sname):
    return "\n".join([
        "#Amplicons = 1",
        "-" * 89,
        "[amplicon1] AmpliconID = 1",
        "[amplicon1] #Intervals = 2",
        "[amplicon1] Intervals = chr7:54000000-55000000,chr7:55500000-56000000",
        "[amplicon1] OncogenesAmplified = EGFR",
        "-" * 89,
    ]) + "\n"


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


def build_sample_dir(root, sname):
    """One sample in AmpliconSuite-pipeline's on-disk shape (no AC dir)."""
    sample_dir = os.path.join(root, sname)
    aa_dir = os.path.join(sample_dir, f"{sname}_AA_results")

    write(os.path.join(aa_dir, f"{sname}_summary.txt"), aa_summary(sname))
    write(os.path.join(aa_dir, f"{sname}.log"),
          f"AmpliconArchitect version {AA_VERSION}\nRunning AA on {sname}\n")
    write(os.path.join(aa_dir, f"{sname}_amplicon1_graph.txt"), graph_txt())
    write(os.path.join(aa_dir, f"{sname}_amplicon1_cycles.txt"), cycles_txt())

    write(os.path.join(sample_dir, f"{sname}.log"),
          f"AmpliconSuite-pipeline version {ASP_VERSION}\n"
          f"AmpliconArchitect version {AA_VERSION}\n")
    write(os.path.join(sample_dir, f"{sname}_finish_flag.txt"),
          "All stages completed\n")
    metadata_jsons(sample_dir, sname)


def result_row(sname):
    return {
        "Sample name": sname,
        "AA amplicon number": "1",
        "Feature ID": f"{sname}_amplicon1_ecDNA_1",
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
        "Reference version": REF,
        "Tissue of origin": "",
        "Sample type": "",
        "Feature BED file": f"{sname}_amplicon1_ecDNA_1_intervals.bed",
        "AS-p version": ASP_VERSION,
        "AA version": AA_VERSION,
        "AC version": AC_VERSION,
        "AA summary file": f"{sname}_AA_results/{sname}_summary.txt",
        "Run metadata JSON": f"{sname}_run_metadata.json",
        "Sample metadata JSON": f"{sname}_sample_metadata.json",
    }


def build_cohort_classification(root):
    """
    One AC run covering both samples — tables hold rows for each, which is
    what forces exclusion to filter at merge time rather than simply skipping
    a directory.
    """
    cls_dir = os.path.join(root, f"{COHORT}_classification")

    write(os.path.join(cls_dir, f"{COHORT}_result_table.tsv"),
          tsv(RT_COLUMNS, [result_row(s) for s in SAMPLES]))

    write(os.path.join(cls_dir, f"{COHORT}_amplicon_classification_profiles.tsv"),
          tsv(PROFILES_COLUMNS, [{
              "sample_name": s, "amplicon_number": "amplicon1",
              "amplicon_decomposition_class": "Cyclic", "ecDNA+": "Positive",
              "BFB+": "None", "ecDNA_amplicons": "1",
          } for s in SAMPLES]))

    write(os.path.join(cls_dir, f"{COHORT}_gene_list.tsv"),
          tsv(GENE_LIST_COLUMNS, [{
              "sample_name": s, "amplicon_number": "amplicon1",
              "feature": "ecDNA_1", "gene": "EGFR", "gene_cn": "8.5",
              "truncated": "None", "is_canonical_oncogene": "TRUE",
              "ncbi_id": "1956",
          } for s in SAMPLES]))

    write(os.path.join(cls_dir, f"{COHORT}_feature_basic_properties.tsv"),
          tsv(FEATURE_PROPS_COLUMNS, [{
              "feature_ID": f"{s}_amplicon1_ecDNA_1",
              "captured_region_size_bp": "1000000",
              "median_feature_CN": "8.5", "max_feature_CN": "9.2",
              "borderline_flag": "None",
          } for s in SAMPLES]))

    write(os.path.join(cls_dir, f"{COHORT}_ecDNA_counts.tsv"),
          tsv(ECDNA_COUNTS_COLUMNS,
              [{"#sample": s, "ecDNA_count": "1"} for s in SAMPLES]))

    # Headerless: column 0 is the Feature ID. Nothing tells a reader which
    # column is which, so exclusion has to address it by index.
    write(os.path.join(cls_dir, f"{COHORT}_ecDNA_context_calls.tsv"),
          "".join(f"{s}_amplicon1_ecDNA_1\tSimple circular\n" for s in SAMPLES))

    write(os.path.join(cls_dir, f"{COHORT}.log"),
          f"AmpliconClassifier {AC_VERSION}\n"
          f"Classification stage completed for {len(SAMPLES)} samples\n")

    # The three flat per-amplicon subdirectories.
    for s in SAMPLES:
        write(os.path.join(cls_dir, f"{COHORT}_classification_bed_files",
                           f"{s}_amplicon1_ecDNA_1_intervals.bed"),
              "chr7\t54000000\t55000000\n")
        write(os.path.join(cls_dir, f"{COHORT}_SV_summaries",
                           f"{s}_amplicon1_SV_summary.tsv"),
              tsv(SV_SUMMARY_COLUMNS, [{
                  "chrom1": "chr7", "pos1": "54000000",
                  "chrom2": "chr7", "pos2": "55000000",
                  "sv_type": "interchromosomal", "read_support": "6",
                  "features": "ecDNA_1", "orientation": "++",
              }]))
        write(os.path.join(cls_dir, f"{COHORT}_annotated_cycles_files",
                           f"{s}_amplicon1_annotated_cycles.txt"),
              cycles_txt())

    # AC's copy of the per-sample inputs.
    files_dir = os.path.join(cls_dir, "files")
    for s in SAMPLES:
        metadata_jsons(files_dir, s)
        write(os.path.join(files_dir, f"{s}_summary.txt"), aa_summary(s))


def build_project(root):
    """Write synthetic inputs under the caller-provided temporary directory."""
    for sname in SAMPLES:
        build_sample_dir(root, sname)
    build_cohort_classification(root)

