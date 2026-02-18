###
# Authors: Edwin Huang, Jens Luebeck
# This script takes in completed AmpliconArchitect and
# AmpliconClassifier results and aggregates the results.
###
import sys
import os
import tarfile
import shutil
import json
import subprocess
import zipfile
from collections import defaultdict
import glob

import pandas as pd

__version__ = "5.2"

EXCLUSION_LIST = ['.bam', '.fq', '.fastq', '.cram', '.fq.gz', '.fastq.gz', 'aln_stage.stderr', '.fai']
# STRING BELOW CANNOT HAVE ANY SPACES IN IT
CLASSIFIER_FILES = "_classification,_amplicon_classification_profiles.tsv,_annotated_cycles_files,_classification_bed_files," \
                   "_ecDNA_counts.tsv,ecDNA_context_calls.tsv,_edge_classification_profiles.tsv,_feature_basic_properties.tsv,_feature_entropy.tsv," \
                   "_feature_similarity_scores.tsv,_features_to_graph.txt,_gene_list.tsv,.input,_SV_summaries,_result_data.json," \
                   "_result_table.tsv,files,index.html".split(',')


def rchop(s, suffix):
    if suffix and s.endswith(suffix):
        return s[:-len(suffix)]
    return s


def string_to_list(string):
    """
    Converts a string representation of a list to a list.
    """
    if "|" in string:
        return string.split("|")
    else:
        return string.strip('][').replace(" ", "").split(',')


def read_name_remap(name_remap_file):
    name_remap = {}
    if name_remap_file:
        with open(name_remap_file) as infile:
            for l in infile:
                fields = l.rstrip().rsplit()
                name_remap[fields[0]] = fields[1]

        print("Name remap has " + str(len(name_remap)) + " keys")

    else:
        print("No name remap file provided. Sample names will detected from filenames.")

    return name_remap


def unzip_file(fp, dest_root):
    """
    unzips file based on zip type
    """
    try:
        if fp.endswith(".tar.gz"):
            zip_name = os.path.basename(fp).replace(".tar.gz", "")
            destination = f'{dest_root}/{zip_name}'
            with tarfile.open(fp, 'r') as output_zip:
                output_zip.extractall(destination)
            output_zip.close()
        elif fp.endswith(".zip"):
            zip_name = os.path.basename(fp).replace(".zip", "")
            destination = f'{dest_root}/{zip_name}'
            with zipfile.ZipFile(fp, 'r') as zip_ref:
                zip_ref.extractall(destination)
            zip_ref.close()

        elif fp.endswith(".tar"):
            zip_name = os.path.basename(fp).replace(".tar", "")
            destination = f'{dest_root}/{zip_name}'
            with tarfile.open(fp, 'r') as output_zip:
                output_zip.extractall(destination)
            output_zip.close()

        else:
            print("File " + fp + " is not a zip or tar file. It may be ignored!")

    except Exception as e:
        print("ERROR WHILE EXTRACTING FILES!")
        print(e)
        if ":" not in str(e):  # not due to legacy ':' in AA files
            sys.exit(1)


def clean_dirs(dlist):
    for d in dlist:
        if d and not d == "/":
            try:
                shutil.rmtree(d)
            except Exception as e:
                print(e)


def check_run_json(dest_root, bname):
    runj_path = os.path.join(dest_root, bname, 'results', 'run.json')
    return os.path.exists(runj_path)


def is_valid_aa_dir(fp):
    """
    Checks if a directory is a valid AA results directory by verifying
    that a _summary.txt file exists and its first line begins with '#Amplicons'.
    """
    try:
        summary_files = [x for x in os.listdir(fp) if x.endswith("_summary.txt")]
        if not summary_files:
            return False
        if len(summary_files) > 1:
            print(f"Multiple summaries in {fp}, suggests mixed AA results")
            return False
        if fp.endswith("/files"):  # this is a classification localization directory
            return False
        if fp.endswith("_AA_results"):
            return True

        summary_path = os.path.join(fp, summary_files[0])
        with open(summary_path, 'r') as sf:
            first_line = sf.readline()
        return first_line.startswith("#Amplicons")
    except OSError as e:
        print(f"Warning: could not read summary in {fp}: {e}")
        return False


class Aggregator:
    def __init__(self, filelist, root, output_name, run_classifier, ref, py3_path, no_clean=False, name_remap_file=None, uuid=""):
        self.root = root
        if not os.path.exists(self.root):
            os.makedirs(self.root)

        self.completed = False
        self.DEST_ROOT = os.path.join(f"{self.root}/extracted_from_zips")
        self.OUTPUT_PATH = os.path.join(f"{self.root}/results/AA_outputs")
        self.OTHER_FILES = os.path.join(f"{self.root}/results/other_files")
        print(f"AmpliconSuiteAggregator {__version__}")
        print("DEST_ROOT: " + self.DEST_ROOT)

        self.zip_paths = filelist
        self.output_name = output_name
        if uuid:
            self.output_name += f"/{uuid}"

        self.aggregated_filename = f'{self.output_name}.tar.gz'

        self.run_classifier = run_classifier
        self.ref = ref
        self.py3_path = py3_path
        self.NO_CLEAN = no_clean
        print(f"No clean mode: {no_clean}")
        self.name_remap = read_name_remap(name_remap_file)

        self.unzip()
        if self.run_classifier == "Yes":
            self.run_amp_classifier()
        self.samp_AA_dct, self.samp_ckit_dct = defaultdict(str), defaultdict(str)
        self.samp_mdata_dct, self.run_mdata_dct = defaultdict(str), defaultdict(str)
        self.samp_cnv_calls_dct = defaultdict(str)
        self.locate_dirs_and_metadata_jsons()
        self.sample_to_ac_location_dct = self.aggregate_tables()
        self.json_modifications()

        # self.move_files()
        self.cleanup()

    def unzip(self):
        for temp_output_dir in [f'{self.root}/results', self.DEST_ROOT]:
            if os.path.exists(temp_output_dir):
                print("Warning: Directory " + temp_output_dir + " already exists! Cleaning now.")
                clean_dirs([temp_output_dir])

        for zip_fp in self.zip_paths:
            fp = os.path.join(self.root, zip_fp)
            try:
                unzip_file(fp, self.DEST_ROOT)
            except Exception as e:
                print("The zip: " + fp + " ran into an error when unzipping.")
                print(e)

        # Extract any nested zips — collect first, then extract (avoid walk/modify race)
        print("Finding compressed files...")
        nested_zips = []
        for root, dirs, files in os.walk(self.DEST_ROOT):
            dirs[:] = [d for d in dirs if "__MACOSX" not in d and "DS_Store" not in d]
            for file in files:
                if file.endswith(".tar.gz") or file.endswith(".zip") or file.endswith(".tar"):
                    nested_zips.append(os.path.join(root, file))
        for fp in nested_zips:
            unzip_file(fp, self.DEST_ROOT)

        if not os.path.exists(self.OUTPUT_PATH):
            os.makedirs(self.OUTPUT_PATH)
        if not os.path.exists(self.OTHER_FILES):
            os.makedirs(self.OTHER_FILES)

        # Pass 1: collect all valid AA results dirs and their parents.
        # Validity is determined by reading the first line of the _summary.txt file.
        # Ending in _AA_results strengthens confidence but is not required.
        print("Crawling files for AA, classification, and CN data...")
        aa_parents = {}      # parent_path -> set of child dir names found under it
        aa_results_paths = set()
        for root, dirs, files in os.walk(self.DEST_ROOT, topdown=True):
            dirs[:] = [d for d in dirs if "__MACOSX" not in d and "DS_Store" not in d]
            for d in dirs:
                fp = os.path.join(root, d)
                if is_valid_aa_dir(fp):
                    if not d.endswith("_AA_results"):
                        print(f"Warning: found valid AA results dir with unexpected name: {fp}")
                    parent = os.path.dirname(fp)
                    aa_parents.setdefault(parent, set()).add(d)
                    aa_results_paths.add(fp)
                elif d.endswith("_AA_results"):
                    print(f"Warning: {fp} looks like _AA_results but failed validation, skipping.")

        # Pass 2: copy each AA parent dir into OUTPUT_PATH.
        # Copying (not moving) keeps the source tree intact for Pass 3.
        copied_parents = set()
        for parent, children in aa_parents.items():
            dest = os.path.join(self.OUTPUT_PATH, os.path.basename(parent))
            try:
                shutil.copytree(parent, dest, dirs_exist_ok=True)
                copied_parents.add(parent)
                print(f"Copied AA parent {parent} -> {dest}")
            except Exception as e:
                print(f"Error copying {parent} to {dest}: {e}")

        # Pass 3: find dirs not under a copied AA parent that belong in OTHER_FILES.
        # Detection rules (direct children only, not recursive):
        #   Classification dirs: ends with _classification, OR directly contains
        #                        a _result_table.tsv, OR directly contains AUX_DIR
        #   CNV dirs: directly contains a _CNV_CALLS.bed file
        #             (ending in _cnvkit_output(s) strengthens confidence but is not required)
        for root, dirs, files in os.walk(self.DEST_ROOT, topdown=True):
            dirs[:] = [d for d in dirs if "__MACOSX" not in d and "DS_Store" not in d]
            for d in dirs:
                fp = os.path.join(root, d)
                if any(fp.startswith(p) for p in copied_parents):
                    continue

                try:
                    dir_contents = set(os.listdir(fp))
                except OSError as e:
                    print(f"Warning: could not list {fp}: {e}")
                    continue

                is_classification_dir = (
                    d.endswith("_classification") or
                    any(f.endswith("_result_table.tsv") for f in dir_contents) or
                    "AUX_DIR" in dir_contents
                )
                is_cnv_dir = any(f.endswith("_CNV_CALLS.bed") for f in dir_contents)

                if is_classification_dir:
                    dest = os.path.join(self.OTHER_FILES, d)
                    try:
                        shutil.copytree(fp, dest, dirs_exist_ok=True)
                        print(f"Copied classification dir {fp} -> {dest}")
                    except Exception as e:
                        print(f"Error copying {fp} to {dest}: {e}")

                elif is_cnv_dir:
                    dest = os.path.join(self.OTHER_FILES, d)
                    try:
                        shutil.copytree(fp, dest, dirs_exist_ok=True)
                        print(f"Copied CNV dir {fp} -> {dest}")
                    except Exception as e:
                        print(f"Error copying {fp} to {dest}: {e}")

    def locate_dirs_and_metadata_jsons(self):
        # Walk both OUTPUT_PATH and OTHER_FILES so that cnvkit dirs that landed
        # in other_files (e.g. per-collection runs) are still discovered and
        # registered in samp_ckit_dct / samp_cnv_calls_dct for the rescue logic
        # in json_modifications().
        print("Locating metadata and unlinked files...")
        for search_path in [self.OUTPUT_PATH, self.OTHER_FILES]:
            for root, dirs, files in os.walk(search_path, topdown=True):
                for dir in dirs:
                    fp = os.path.join(root, dir)
                    if not os.path.exists(fp):
                        continue

                    try:
                        dir_contents = os.listdir(fp)
                    except OSError as e:
                        print(f"Warning: could not list {fp}: {e}")
                        continue

                    if fp.endswith("_AA_results"):
                        implied_sname = rchop(fp, "_AA_results").rsplit("/")[-1]
                        self.samp_AA_dct[implied_sname] = fp

                    elif fp.endswith("_cnvkit_output") or fp.endswith("_cnvkit_outputs"):
                        # Use rchop with the longer suffix first so that the trailing
                        # 's' variant is handled correctly with a single call.
                        suffix = "_cnvkit_outputs" if fp.endswith("_cnvkit_outputs") else "_cnvkit_output"
                        implied_sname = rchop(os.path.basename(fp), suffix)
                        self.clean_by_suffix(".cnr.gz", fp)
                        self.clean_by_suffix(".cnr", fp)
                        self.clean_by_suffix("call.cns.gz", fp)
                        self.clean_by_suffix("call.cns", fp)
                        cmd = f'gzip -fq {fp}/*.cns 2> /dev/null'
                        subprocess.call(cmd, shell=True)
                        self.samp_ckit_dct[implied_sname] = fp

                        # Register CNV_CALLS.bed using the dir-derived sample name.
                        # Older AmpSuite-pipeline used the bam basename in the bed
                        # filename rather than the sample name, so we must not derive
                        # the sample name from the bed filename in this case.
                        for f in dir_contents:
                            if f.endswith("_CNV_CALLS.bed"):
                                self.samp_cnv_calls_dct[implied_sname] = fp + "/" + f

                    for f in dir_contents:
                        if f.endswith("_run_metadata.json"):
                            implied_sname = rchop(f, "_run_metadata.json")
                            self.run_mdata_dct[implied_sname] = fp + "/" + f
                        elif f.endswith("_sample_metadata.json"):
                            implied_sname = rchop(f, "_sample_metadata.json")
                            self.samp_mdata_dct[implied_sname] = fp + "/" + f
                        elif f.endswith("_CNV_CALLS.bed"):
                            # Only handle _CNV_CALLS.bed here for non-cnvkit dirs.
                            # cnvkit dirs handle their own _CNV_CALLS.bed above using
                            # the dir-derived sample name to avoid bam-basename confusion.
                            if not (fp.endswith("_cnvkit_output") or fp.endswith("_cnvkit_outputs")):
                                cnv_sname = rchop(f, "_CNV_CALLS.bed")
                                self.samp_cnv_calls_dct[cnv_sname] = fp + "/" + f

    def run_amp_classifier(self):
        """
        Goes into the OUTPUT_PATH to look for and delete amplicon classifier results.

        will run amplicon classifier on AA_results
        1. downloads genome based on
            - run metadata
            - log
            - default to GRCh38 or user input
        2. run make_inputs.sh on sample directories
        3. run Amplicon Classifier on each sample
        """

        print('************ STARTING AMPLICON CLASSIFIER **************')
        print(f'removing classifier files and directories ending with:\n {CLASSIFIER_FILES}')
        removed_files = 0
        for root, dirs, files in os.walk(self.OUTPUT_PATH):
            for name in files:
                if any([name.endswith(i) for i in CLASSIFIER_FILES]):
                    print(name, "to remove")
                    fp = os.path.join(root, name)
                    try:
                        os.remove(fp)
                        removed_files += 1
                    except Exception as e:
                        print(e)

            for name in dirs:
                if any([name.endswith(i) for i in CLASSIFIER_FILES]):
                    print(name, "to remove")
                    fp = os.path.join(root, name)
                    shutil.rmtree(fp, ignore_errors=True)
                    removed_files += 1

        print(f"Removed {removed_files} files and directories from previous classification results")

        ## run amplicon classifier on AA_outputs:

        ## 1. make inputs
        try:
            AC_SRC = os.environ['AC_SRC']
        except KeyError:
            sys.stderr.write("AC_SRC variable not found! AmpliconClassifier is not properly installed.\n")
            self.cleanup(failure=True)
            sys.exit(1)

        print(f"AC_SRC is set to {AC_SRC}")
        input_ec = os.system(f"{AC_SRC}/make_input.sh {self.OUTPUT_PATH} {self.OUTPUT_PATH}/{self.output_name}")
        if input_ec != 0:
            print("Failed to make input for AC!")
            self.cleanup(failure=True)
            sys.exit(1)

        ## if reference isn't downloaded already, then download appropriate reference genome
        try:
            local_data_repo = os.environ['AA_DATA_REPO']
        except KeyError:
            sys.stderr.write(
                "AA_DATA_REPO variable not found! The AA data repo directory is not properly configured.\n")
            self.cleanup(failure=True)
            sys.exit(1)

        if not os.path.exists(os.path.join(local_data_repo, self.ref)):
            os.system(
                f"wget -q -P {local_data_repo} https://datasets.genepattern.org/data/module_support_files/AmpliconArchitect/{self.ref}.tar.gz")
            os.system(
                f"wget -q -P {local_data_repo} https://datasets.genepattern.org/data/module_support_files/AmpliconArchitect/{self.ref}_md5sum.tar.gz")
            os.system(f"tar zxf {local_data_repo}/{self.ref}.tar.gz --directory {local_data_repo}")

        ## 2. run amplicon classifier.py
        os.system(
            f"{self.py3_path} {AC_SRC}/amplicon_classifier.py -i {self.OUTPUT_PATH}/{self.output_name}.input --ref {self.ref} -o {self.OUTPUT_PATH}/{self.output_name} ")

        ## 3. run make_results_table.py
        os.system(f"{self.py3_path} {AC_SRC}/make_results_table.py --input {self.OUTPUT_PATH}/{self.output_name}.input \
                  --classification_file {self.OUTPUT_PATH}/{self.output_name}_amplicon_classification_profiles.tsv \
                    --summary_map {self.OUTPUT_PATH}/{self.output_name}_summary_map.txt \
                        --ref {self.ref}")

        print('************ ENDING AMPLICON CLASSIFIER **************')

    def aggregate_tables(self):
        """
        From the unzipped paths, start aggregating results
        This function will put together the separate AA runs into one run.
        """

        output_head = ["Sample name", "AA amplicon number", "Feature ID", "Classification", "Location", "Oncogenes",
                       "Complexity score", "Captured interval length", "Feature median copy number",
                       "Feature maximum copy number", "Filter flag", "Reference version", "Tissue of origin",
                       "Sample type", "Feature BED file", "CNV BED file", "AA PNG file", "AA PDF file",
                       "AA summary file",
                       "Run metadata JSON", "Sample metadata JSON"]
        aggregate = pd.DataFrame(columns=output_head)
        runs = {}
        sample_to_ac_dir = {}
        sample_num = 1
        # aggregate results
        print("Aggregating results into tables")
        found_res_table = False
        for res_dir in [self.OUTPUT_PATH, self.OTHER_FILES]:
            for root, dirs, files in os.walk(res_dir, topdown=False):
                for name in files:
                    if name.endswith("_result_table.tsv") and not name.startswith("._"):
                        result_table_fp = os.path.join(root, name)
                        found_res_table = True
                        try:
                            df = pd.read_csv(result_table_fp, delimiter='\t')
                            aggregate = pd.concat([aggregate, df], ignore_index=True)
                        except Exception as e:
                            print(e)
                            continue

                        gdf = df.groupby('Sample name')
                        for sname, group in gdf:
                            print(f'{sname} has {str(len(group))} rows')
                            # sample_N key is retained as downstream code depends on it;
                            # the actual sample name is encoded inside each record dict.
                            runs[f'sample_{sample_num}'] = json.loads(group.to_json(orient='records'))
                            if sname in sample_to_ac_dir:
                                print(f"Warning: duplicate sample name {sname} found in multiple result tables.")
                            sample_to_ac_dir[f'sample_{sample_num}'] = os.path.dirname(result_table_fp)
                            sample_num += 1

        if not found_res_table:
            print(
                "Error: No results tables found! Aggregation will be empty or invalid. Please make sure to first run make_results_table.py from AmpliconClassifier.")
            self.cleanup(failure=True)
            self.completed = False
            return {}

        ## output the table
        with open(f'{self.root}/results/run.json', 'w') as run_file:
            json.dump({'runs': runs}, run_file)

        run_file.close()
        return sample_to_ac_dir

    def tardir(self, path, tar_name, keep_root=True):
        """
        compresses the path to a tar_name
        """
        with tarfile.open(tar_name, "w:gz") as tar_handle:
            base_path = os.path.dirname(path) if keep_root else path

            for root, dirs, files in os.walk(path):
                for file in files:
                    if not any([file.endswith(x) for x in EXCLUSION_LIST]):
                        file_path = os.path.join(root, file)

                        if keep_root:
                            # Create archive path relative to parent of target directory
                            arcname = os.path.relpath(file_path, base_path)
                        else:
                            # Just use the immediate parent directory
                            arcname = os.path.join(os.path.basename(root), file)

                        tar_handle.add(file_path, arcname=arcname)

    def cleanup(self, failure=False):
        """
        Zips the aggregate results, and deletes files for cleanup
        """
        if not self.NO_CLEAN:
            clean_dirs(self.samp_AA_dct.values())
        # self.clean_files(self.samp_ckit_dct.values())
        if not failure:
            print(f"Creating {self.aggregated_filename} and cleaning unpacked files...")
            self.tardir(f'{self.root}/results', self.aggregated_filename)

        if not self.NO_CLEAN:
            clean_dirs([f'{self.root}/results', self.DEST_ROOT])  # ./extracted_from_zips

    def json_modifications(self):
        """
        Modify outputs of the final run.json

        Modifications include:
        1. correcting string of list
        2. absolute pathing to relative pathing
        """
        with open(f'{self.root}/results/run.json') as json_file:
            dct = json.load(json_file)
        json_file.close()

        ref_genomes = set()
        processed_feat_count = 0
        print("Forming JSONs...")
        for sample in dct['runs'].keys():
            for sample_dct in dct['runs'][sample]:
                processed_feat_count += 1
                # updating string of lists to lists
                sample_name = sample_dct["Sample name"]
                ref_genomes.add(sample_dct["Reference version"])
                if sample_dct["Reference version"] is None:
                    sys.stderr.write("WARNING: " + sample_name + " had no reference genome build indicated!\n")
                if len(ref_genomes) > 1:
                    sys.stderr.write(str(ref_genomes) + "\n")
                    sys.stderr.write("ERROR! Multiple reference genomes detected in project.\n AmpliconRepository only "
                                     "supports single-reference projects currently. Exiting.\n")
                    self.cleanup(failure=True)
                    return

                potential_str_lsts = [
                    'Location',
                    'Oncogenes',
                    'All genes',
                ]
                for lname in potential_str_lsts:
                    try:
                        sample_dct[lname] = string_to_list(sample_dct[lname])
                    except:
                        print(f'Feature: {lname} doesnt exist for sample: {sample_name}')

                # update each path in run.json by finding them in outputs folder
                # separately for feature bed file because location is different
                feat_basename = os.path.basename(sample_dct['Feature BED file'])
                cfiles = os.listdir(self.sample_to_ac_location_dct[sample])
                cbf_hits = [x for x in cfiles if x.endswith("_classification_bed_files") and not x.startswith("._")]
                if cbf_hits:
                    cbf = cbf_hits[0]
                    feat_file = f'{self.sample_to_ac_location_dct[sample]}/{cbf}/{feat_basename}'
                    if os.path.exists(feat_file):
                        feat_file = feat_file.replace(f'{self.root}/results/', "")
                        sample_dct['Feature BED file'] = feat_file
                    else:
                        if not feat_file.endswith("/NA"):
                            print(
                                f'Feature: "Feature BED file" {feat_file} doesnt exist for sample {sample_dct["Sample name"]}')

                        sample_dct['Feature BED file'] = "Not Provided"

                else:
                    sample_dct['Feature BED file'] = "Not Provided"

                features_of_interest = [
                    'CNV BED file',
                    'AA PNG file',
                    'AA PDF file',
                    'AA summary file',
                    'Run metadata JSON',
                    'Sample metadata JSON',
                ]

                for feature in features_of_interest:
                    if feature in sample_dct and sample_dct[feature]:
                        feat_basename = os.path.basename(sample_dct[feature])
                        feat_file = f'{self.sample_to_ac_location_dct[sample]}/files/{feat_basename}'
                        if feature == "CNV BED file" and any([feat_file.endswith(x) for x in
                                                              ["AA_CNV_SEEDS.bed", "CNV_CALLS_pre_filtered.bed",
                                                               "Not provided", "Not Provided"]]):
                            if self.samp_cnv_calls_dct[sample_name]:
                                feat_file = self.samp_cnv_calls_dct[sample_name]

                            cnvkit_dir = self.samp_ckit_dct[sample_dct['Sample name']]
                            if cnvkit_dir:
                                for f in os.listdir(cnvkit_dir):
                                    if f.endswith("_CNV_CALLS.bed"):
                                        feat_file = cnvkit_dir + "/" + f

                        elif feature == "Run metadata JSON" and any(
                                [feat_file.endswith(x) for x in ["Not provided", "Not Provided"]]):
                            rmj = self.run_mdata_dct[sample_dct['Sample name']]
                            if rmj:
                                feat_file = rmj

                        elif feature == "Sample metadata JSON" and any(
                                [feat_file.endswith(x) for x in ["Not provided", "Not Provided"]]):
                            smj = self.samp_mdata_dct[sample_dct['Sample name']]
                            if smj:
                                feat_file = smj

                        if os.path.exists(feat_file):
                            feat_file = feat_file.replace(f'{self.root}/results/', "")
                            sample_dct[feature] = feat_file

                        else:
                            if feature not in ['Run metadata JSON', 'Sample metadata JSON']:
                                print(f'Feature: "{feature}" file doesn\'t exist for sample {sample_dct["Sample name"]}: {feat_file}')

                            sample_dct[feature] = "Not Provided"

                    else:
                        sample_dct[feature] = "Not Provided"

                dirs_of_interest = [
                    'AA directory',
                    'cnvkit directory'
                ]
                for dirfname, sdct in zip(dirs_of_interest, [self.samp_AA_dct, self.samp_ckit_dct]):
                    if sdct[sample_dct['Sample name']]:
                        orig_dir = sdct[sample_dct['Sample name']]
                        # clean .out and .cnr.gz files
                        self.clean_by_suffix("*.out", orig_dir)
                        self.clean_by_suffix("*.cnr.gz", orig_dir)
                        self.clean_by_suffix("*.cnr", orig_dir)

                        tarf = orig_dir + ".tar.gz"
                        self.tardir(orig_dir, tarf, keep_root=False)
                        sample_dct[dirfname] = tarf.replace(f'{self.root}/results/', "")

                    else:
                        sample_dct[dirfname] = "Not Provided"

                # rename the sample if user specifies
                if sample_name in self.name_remap:
                    sample_dct['Sample name'] = self.name_remap[sample_name]
                elif self.name_remap and sample_name not in self.name_remap:
                    sys.stderr.write("Could not rename sample: " + sample_name + "\n")

                if processed_feat_count % 100 == 0:
                    print("Processed " + str(processed_feat_count) + " features...")

        with open(f'{self.root}/results/run.json', 'w') as json_file:
            json.dump(dct, json_file, sort_keys=True, indent=2)
        json_file.close()

        flattened_samples = []
        for _, v in dct['runs'].items():
            flattened_samples.extend(v)

        aggregate = pd.DataFrame.from_records(flattened_samples)
        aggregate.to_csv(f'{self.root}/results/aggregated_results.csv')
        aggregate.to_html(f'{self.root}/results/aggregated_results.html')
        self.completed = True

    def clean_by_suffix(self, suffix, dir):
        # Validate inputs
        if suffix and dir and dir != "/" and suffix != "*":
            if not suffix.startswith("*"):
                suffix = "*" + suffix

            # Use glob to safely handle file listing
            file_pattern = os.path.join(dir, suffix)
            files_to_remove = glob.glob(file_pattern)
            for file_path in files_to_remove:
                try:
                    os.remove(file_path)
                except OSError as e:
                    print(f"Error removing {file_path}: {e}")


# TODO: VALIDATE IS NEVER USED!
def validate(root):
    """
    Validates that all the file locations exist.
    """
    ## check that everything is in the right place
    if os.path.exists(f'{root}/output'):
        print('output folder exists')
    else:
        raise Exception("OUTPUT FOLDER NOT CREATED...")

    with open(f'{root}/results/run.json') as json_file:
        dct = json.load(json_file)
    json_file.close()

    ## TODO: go through fields to see if each datatype is consistent
    ## and that paths exist