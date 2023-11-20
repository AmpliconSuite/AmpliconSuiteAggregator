###
# Author: Edwin Huang
# This script takes in completed AmpliconArchitect and
# AmpliconClassifier results and aggregates the results.
###
import sys
import tarfile
import os
import re
import pandas as pd
import tarfile
import shutil
import json
import subprocess
import ast
import zipfile
from collections import defaultdict
import requests


DEST_ROOT = os.path.join("./extracted_from_zips")
OUTPUT_PATH = os.path.join("./results/AA_outputs")
OTHER_FILES = os.path.join("./results/other_files")
EXCLUSION_LIST = ['.bam', '.fq', '.fastq', '.cram', '.fq.gz', '.fastq.gz']
# STRING BELOW CANNOT HAVE ANY SPACES IN IT
CLASSIFIER_FILES = "_classification,_amplicon_classification_profiles.tsv,_annotated_cycles_files,_classification_bed_files," \
                   "_ecDNA_counts.tsv,_edge_classification_profiles.tsv,_feature_basic_properties.tsv,_feature_entropy.tsv," \
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

    except Exception as e:
        print(e)


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


class Aggregator:
    def __init__(self, filelist, root, output_name, run_classifier, ref, py3_path, name_remap_file=None):
        self.zip_paths = filelist
        self.root = root
        self.output_name = output_name
        self.run_classifier = run_classifier
        self.ref = ref
        self.py3_path = py3_path
        self.name_remap = read_name_remap(name_remap_file)

        self.unzip()
        if self.run_classifier == "Yes":
            self.run_amp_classifier()
        self.samp_AA_dct, self.samp_ckit_dct = defaultdict(str), defaultdict(str)
        self.samp_mdata_dct, self.run_mdata_dct = defaultdict(str), defaultdict(str)
        self.locate_dirs_and_metadata_jsons()
        # print(self.samp_ckit_dct)
        # print(self.samp_AA_dct)
        # print(self.samp_mdata_dct)
        # print(self.run_mdata_dct)
        self.sample_to_ac_location_dct = self.aggregate_tables()
        self.json_modifications()
        # self.move_files()
        self.cleanup()

    def unzip(self):
        """
        Unzips the zip files, and get directories for files within

        """
        for zip_fp in self.zip_paths:
            fp = os.path.join(self.root, zip_fp)
            try:
                unzip_file(fp, DEST_ROOT)

            except Exception as e:
                print("The zip: " + fp + " ran into an error when unzipping.")
                print(e)

        ## find all nested zips and extract them.
        print("Finding compressed files...")
        for root, dirs, files in os.walk(DEST_ROOT):
            for file in files:
                if file.endswith(".tar.gz") or file.endswith(".zip"):
                    fp = os.path.join(root, file)
                    unzip_file(fp, DEST_ROOT)

        ## moving required AA files to correct destination
        if not os.path.exists(OUTPUT_PATH):
            os.makedirs(OUTPUT_PATH)

        if not os.path.exists(OTHER_FILES):
            os.makedirs(OTHER_FILES)

        ## find samples and move files
        # samples = []
        aa_samples_found = 0
        print("Crawling files for AA, classification, and CN data...")
        for root, dirs, files in os.walk(DEST_ROOT, topdown = True):
            for dir in dirs: 
                fp = os.path.join(root, dir)
                if "__MACOSX" in fp or "DS_Store" in fp or not os.path.exists(fp):
                    continue

                if fp.endswith("_AA_results"):
                    if any([x.endswith("_summary.txt") for x in os.listdir(fp)]):
                        print(f'Moving AA results to {OUTPUT_PATH}')
                        os.system(f'mv -vf {os.path.dirname(fp)} {OUTPUT_PATH}')
                        aa_samples_found+=1
                        if aa_samples_found % 100 == 0:
                            print("Crawled through " + str(aa_samples_found) + " AA results...")

                if fp.endswith("_cnvkit_output"):
                    print(f'Moving cnvkit_outputs to {OUTPUT_PATH}')
                    os.system(f'mv -vf {os.path.dirname(fp)} {OUTPUT_PATH}')

                if os.path.exists(fp) and any([x.endswith("_result_table.tsv") or x == "AUX_DIR" for x in os.listdir(fp)]):
                    pre = fp.rstrip("_classification")
                    # don't need to move things that will get brought already into AA_outputs
                    if not fp.endswith("_classification") and not os.path.exists(pre + "_AA_results"):
                        print(f'Moving file {fp} to {OTHER_FILES}')
                        os.system(f'mv -vf {fp} {OTHER_FILES}')

    def locate_dirs_and_metadata_jsons(self):
        # post-move identification of AA and CNVKit files:
        print("Locating metadata...")
        for root, dirs, files in os.walk(OUTPUT_PATH, topdown = True):
            for dir in dirs:
                fp = os.path.join(root, dir)
                if not os.path.exists(fp):
                    continue

                if fp.endswith("_AA_results"):
                    implied_sname = rchop(fp,"_AA_results").rsplit("/")[-1]
                    self.samp_AA_dct[implied_sname] = fp

                elif fp.endswith("_cnvkit_output"):
                    implied_sname = rchop(fp,"_cnvkit_output").rsplit("/")[-1]
                    self.clean_by_suffix(".cnr.gz", fp)
                    cmd = f'gzip -fq {fp}/*.cns 2> /dev/null'
                    subprocess.call(cmd, shell=True)
                    # print(fp.rstrip("_cnvkit_output"), implied_sname)
                    self.samp_ckit_dct[implied_sname] = fp

                for f in os.listdir(fp):
                    if f.endswith("_run_metadata.json"):
                        implied_sname = rchop(f, "_run_metadata.json")
                        self.run_mdata_dct[implied_sname] = fp + "/" + f

                    elif f.endswith("_sample_metadata.json"):
                        implied_sname = rchop(f, "_sample_metadata.json")
                        self.samp_mdata_dct[implied_sname] = fp + "/" + f

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

        # for each sample in AA_outputs

        # python3 /files/src/AmpliconSuiteAggregator.py -flist /files/gpunit/inputs/input_list.txt --run_classifier Yes -ref GRCh38

        print('************ STARTING AMPLICON CLASSIFIER **************')
        print(f'removing classifier files and directories ending with:\n {CLASSIFIER_FILES}')
        removed_files = 0
        for root, dirs, files in os.walk(OUTPUT_PATH):
            for name in files:
                if any([name.endswith(i) for i in CLASSIFIER_FILES]):
                    print(name, "to remove")
                    fp = os.path.join(root, name)
                    try:
                        os.remove(fp)
                        removed_files+=1
                    except Exception as e:
                        print(e)

            for name in dirs:
                if any([name.endswith(i) for i in CLASSIFIER_FILES]):
                    print(name, "to remove")
                    fp = os.path.join(root, name)
                    shutil.rmtree(fp, ignore_errors=True)
                    removed_files+=1

        print(f"Removed {removed_files} files and directories from previous classification results")

        ## run amplicon classifier on AA_outputs:

        ## 1. make inputs
        try:
            AC_SRC = os.environ['AC_SRC']
        except KeyError:
            sys.stderr.write("AC_SRC variable not found! AmpliconClassifier is not properly installed.\n")
            sys.exit(1)

        print(f"AC_SRC is set to {AC_SRC}")
        os.system(f"{AC_SRC}/make_input.sh {OUTPUT_PATH} {OUTPUT_PATH}/{self.output_name}" )

        ## if reference isn't downloaded already, then download appropriate reference genome
        try:
            local_data_repo = os.environ['AA_DATA_REPO']
        except KeyError:
            sys.stderr.write("AA_DATA_REPO variable not found! The AA data repo directory is not properly configured.\n")
            sys.exit(1)

        if not os.path.exists(os.path.join(local_data_repo, self.ref)):
            os.system(f"wget -q -P {local_data_repo} https://datasets.genepattern.org/data/module_support_files/AmpliconArchitect/{self.ref}.tar.gz")
            os.system(f"wget -q -P {local_data_repo} https://datasets.genepattern.org/data/module_support_files/AmpliconArchitect/{self.ref}_md5sum.tar.gz")
            os.system(f"tar zxf {local_data_repo}/{self.ref}.tar.gz --directory {local_data_repo}")

        ## 2. run amplicon classifier.py
        os.system(f"{self.py3_path} {AC_SRC}/amplicon_classifier.py -i {OUTPUT_PATH}/{self.output_name}.input --ref {self.ref} -o {OUTPUT_PATH}/{self.output_name} ")
        
        ## 3. run make_results_table.py
        os.system(f"{self.py3_path} {AC_SRC}/make_results_table.py --input {OUTPUT_PATH}/{self.output_name}.input \
                  --classification_file {OUTPUT_PATH}/{self.output_name}_amplicon_classification_profiles.tsv \
                    --summary_map {OUTPUT_PATH}/{self.output_name}_summary_map.txt \
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
                       "Sample type", "Feature BED file", "CNV BED file", "AA PNG file", "AA PDF file", "AA summary file",
                       "Run metadata JSON", "Sample metadata JSON"]
        aggregate = pd.DataFrame(columns = output_head)
        runs = {}
        sample_to_ac_dir = {}
        sample_num = 1
        # aggregate results
        print("Aggregating results into tables")
        for res_dir in [OUTPUT_PATH, OTHER_FILES]:
            for root, dirs, files in os.walk(res_dir, topdown = False):
                for name in files:
                    if name.endswith("_result_table.tsv") and not name.startswith("._"):
                        result_table_fp = os.path.join(root, name)
                        try:
                            df = pd.read_csv(result_table_fp, delimiter = '\t')
                            aggregate = pd.concat([aggregate, df], ignore_index = True)

                        except Exception as e:
                            print(e)
                            continue

                        # print(df)
                        gdf = df.groupby('Sample name')
                        for sname, group in gdf:
                            print(f'{sname} has {str(len(group))} rows')
                            runs[f'sample_{sample_num}'] = json.loads(group.to_json(orient = 'records'))
                            sample_to_ac_dir[f'sample_{sample_num}'] = os.path.dirname(result_table_fp)
                            sample_num += 1

        ## output the table
        with open('./results/run.json', 'w') as run_file:
            json.dump({'runs': runs}, run_file)

        run_file.close()
        return sample_to_ac_dir

    def tardir(self, path, tar_name, keep_root=True):
        """
        compresses the path to a tar_mname
        """
        with tarfile.open(tar_name, "w:gz") as tar_handle:
            for root, dirs, files in os.walk(path):
                for file in files:
                    if not any([file.endswith(x) for x in EXCLUSION_LIST]):
                        if not keep_root:
                            arcname = os.path.join(os.path.basename(root), file)

                        else:
                            arcname = None

                        tar_handle.add(os.path.join(root, file), arcname=arcname)

        tar_handle.close()

    def cleanup(self):
        """
        Zips the aggregate results, and deletes files for cleanup
        """
        clean_dirs(self.samp_AA_dct.values())
        # self.clean_files(self.samp_ckit_dct.values())
        print("Creating tar.gz...")
        self.tardir('./results', f'{self.output_name}.tar.gz')
        clean_dirs(['./results', DEST_ROOT]) # ./extracted_from_zips

    # def find_file(self, basename):
    #     """
    #     Finds a specific file given sample name and file suffix
    #     """
    #     for root, dirs, files in os.walk('./results/AA_outputs', topdown = True):
    #         for file in files:
    #             fp = os.path.join(root, file)
    #             if basename == file:
    #                 return fp.replace('./results/', "")
    #
    #     return 'Not Provided'

    def json_modifications(self):
        """
        Modify outputs of the final run.json

        Modifications include:
        1. correcting string of list
        2. absolute pathing to relative pathing
        """
        with open('./results/run.json') as json_file:
            dct = json.load(json_file)
        json_file.close()

        ref_genomes = set()
        processed_feat_count = 0
        print("Forming JSONs...")
        for sample in dct['runs'].keys():
            for sample_dct in dct['runs'][sample]:
                processed_feat_count+=1
                # updating string of lists to lists
                sample_name = sample_dct["Sample name"]
                ref_genomes.add(sample_dct["Reference version"])
                if len(ref_genomes) > 1:
                    sys.stderr.write(str(ref_genomes) + "\n")
                    sys.stderr.write("ERROR! Multiple reference genomes detected in project.\n AmpliconRepository only "
                                     "supports single-reference projects currently. Exiting.\n")
                    sys.exit(1)

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
                    # print("Feature BED file", feat_file)
                    if os.path.exists(feat_file):
                        feat_file = feat_file.replace('./results/', "")
                        sample_dct['Feature BED file'] = feat_file
                    else:
                        if not feat_file.endswith("/NA"):
                            print(f'Feature: "Feature BED file" {feat_file} doesnt exist for sample {sample_dct["Sample name"]}')

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
                        if feature == "CNV BED file" and any([feat_file.endswith(x) for x in ["AA_CNV_SEEDS.bed", "CNV_CALLS_pre_filtered.bed", "Not provided", "Not Provided"]]):
                            cnvkit_dir = self.samp_ckit_dct[sample_dct['Sample name']]
                            if cnvkit_dir:
                                for f in os.listdir(cnvkit_dir):
                                    if f.endswith("_CNV_CALLS.bed"):
                                        feat_file = cnvkit_dir + "/" + f

                        elif feature == "Run metadata JSON" and any([feat_file.endswith(x) for x in ["Not provided", "Not Provided"]]):
                            rmj = self.run_mdata_dct[sample_dct['Sample name']]
                            if rmj:
                                feat_file = rmj

                        elif feature == "Sample metadata JSON" and any([feat_file.endswith(x) for x in ["Not provided", "Not Provided"]]):
                            smj = self.samp_mdata_dct[sample_dct['Sample name']]
                            if smj:
                                feat_file = smj

                        if os.path.exists(feat_file):
                            feat_file = feat_file.replace('./results/', "")
                            sample_dct[feature] = feat_file

                        else:
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

                        tarf = orig_dir + ".tar.gz"
                        self.tardir(orig_dir, tarf, keep_root=False)
                        sample_dct[dirfname] = tarf.replace('./results/', "")

                    else:
                        sample_dct[dirfname] = "Not Provided"

                # rename the sample if user specifies
                if sample_name in self.name_remap:
                    # print("Renaming sample " + sample_name + " to " + self.name_remap[sample_name])
                    sample_dct['Sample name'] = self.name_remap[sample_name]
                elif self.name_remap and sample_name not in self.name_remap:
                    sys.stderr.write("Could not rename sample: " + sample_name + "\n")

                if processed_feat_count % 100 == 0:
                    print("Processed " + str(processed_feat_count) + " features...")

        with open('./results/run.json', 'w') as json_file:
            json.dump(dct, json_file, sort_keys=True, indent=2)
        json_file.close()

        flattened_samples = []
        for _,v in dct['runs'].items():
            flattened_samples.extend(v)

        aggregate = pd.DataFrame.from_records(flattened_samples)
        aggregate.to_csv('./results/aggregated_results.csv')
        aggregate.to_html('./results/aggregated_results.html')

    def clean_by_suffix(self, suffix, dir):
        if suffix and dir and not dir == "/" and not suffix == "*":
            if not suffix.startswith("*"):
                suffix = "*" + suffix

            cmd = f'rm -f {dir}/{suffix}'
            # try:
            subprocess.call(cmd, shell=True)
            # except FileNotFoundError:
            #     pass


# TODO: VALIDATE IS NEVER USED!
def validate():
    """
    Validates that all of the file locations exist.
    """
    ## check that everything is in the right place
    if os.path.exists('./output'):
        print('output folder exists')
    else:
        raise Exception("OUTPUT FOLDER NOT CREATED...")
    
    with open('./results/run.json') as json_file:
        dct = json.load(json_file)
    json_file.close()

    ## TODO: go through fields to see if each datatype is consistent
    ## and that paths exist
