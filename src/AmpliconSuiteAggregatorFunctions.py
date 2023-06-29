###
# Author: Edwin Huang
# This script takes in completed AmpliconArchitect and
# AmpliconClassifier results and aggregates the results.
###
import argparse
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
import intervaltree


DEST_ROOT = os.path.join("./extracted_from_zips")
OUTPUT_PATH = os.path.join("./results/AA_outputs")
OTHER_FILES = os.path.join("./results/other_files")
global EXCLUSION_LIST
EXCLUSION_LIST = ['.bam', '.gz', '.fq', '.fastq', '.cram']
CLASSIFIER_FILES = "_classification, _amplicon_classification_profiles.tsv,_annotated_cycles_files, _classification_bed_files, _ecDNA_counts.tsv, _edge_classification_profiles.tsv, _feature_basic_properties.tsv, _feature_entropy.tsv, _feature_similarity_scores.tsv, _features_to_graph.txt, _gene_list.tsv, .input, _SV_summaries,_result_data.json,_result_table.tsv,files,index.html".split(',')



class Aggregator():
    def __init__(self, filelist, root, output_name, run_classifier, ref):
        self.zip_paths = filelist
        self.root = root
        self.output_name = output_name
        self.run_classifier = run_classifier
        self.ref = ref
        self.unzip()
        if self.run_classifier == "Yes":
            self.run_amp_classifier()
        self.samp_AA_dct, self.samp_ckit_dct = defaultdict(str), defaultdict(str)
        self.samp_mdata_dct, self.run_mdata_dct = defaultdict(str), defaultdict(str)
        self.locate_dirs_and_metadata_jsons()
        print(self.samp_ckit_dct)
        print(self.samp_AA_dct)
        self.sample_to_ac_location_dct = self.aggregate_tables()
        self.json_modifications()
        self.move_files()
        self.cleanup()

    def unzip_file(self, fp, dest_root):
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

            # else:
            #     print(f'No need to extract: {fp}')

        except Exception as e:
            print(e)

    def unzip(self):
        """
        Unzips the zip files, and get directories for files within

        """

        for zip_fp in self.zip_paths:
            fp = os.path.join(self.root, zip_fp)
            try:
                self.unzip_file(fp, DEST_ROOT)

            except Exception as e:
                print("The zip: " + fp + " ran into an error when unzipping.")
                print(e)

        ## find all nested zips and extract them. 
        for root, dirs, files in os.walk(DEST_ROOT):
            for file in files:
                fp = os.path.join(root, file)
                self.unzip_file(fp, DEST_ROOT)

        ## moving required AA files to correct destination
        if not os.path.exists(OUTPUT_PATH):
            os.makedirs(OUTPUT_PATH)

        if not os.path.exists(OTHER_FILES):
            os.makedirs(OTHER_FILES)

        ## find samples and move files
        # samples = []
        for root, dirs, files in os.walk(DEST_ROOT, topdown = True):
            for dir in dirs: 
                fp = os.path.join(root, dir)
                if "__MACOSX" in fp or "DS_Store" in fp or not os.path.exists(fp):
                    continue

                if fp.endswith("_AA_results"):
                    if any([x.endswith("_summary.txt") for x in os.listdir(fp)]):
                        print(f'Moving AA results to {OUTPUT_PATH}')
                        os.system(f'mv -vf {os.path.dirname(fp)} {OUTPUT_PATH}')

                if fp.endswith("_cnvkit_output"):
                    # tarname = os.path.dirname(fp) + ".tar.gz"
                    # self.tardir(os.path.dirname(fp), tarname)
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
        for root, dirs, files in os.walk(OUTPUT_PATH, topdown = True):
            for dir in dirs:
                fp = os.path.join(root, dir)
                if not os.path.exists(fp):
                    continue

                if fp.endswith("_AA_results"):
                    implied_sname = fp.rstrip("_AA_results").rsplit("/")[-1]
                    self.samp_AA_dct[implied_sname] = fp

                elif fp.endswith("_cnvkit_output"):
                    implied_sname = fp.rstrip("_cnvkit_output").rsplit("/")[-1]
                    self.samp_ckit_dct[implied_sname] = fp

                for f in os.listdir(fp):
                    if f.endswith("_run_metadata.json"):
                        implied_sname = f.rstrip("_run_metadata.json")
                        self.run_mdata_dct[implied_sname] = fp + "/" + f

                    elif f.endswith("_sample_metadata.json"):
                        implied_sname = f.rstrip("_sample_metadata.json")
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


        print('************ STARTING AMPLICON CLASSIFIER NOW **************')
        print(f'removing classifier files and directories{CLASSIFIER_FILES}')
        for root, dirs, files in os.walk(OUTPUT_PATH):
            for name in files:
                ## paths to files
                fp = os.path.join(root, name)
                for i in CLASSIFIER_FILES:
                    if i in fp:
                        try:
                            os.remove(fp)
                        except Exception as e:
                            print(e)

            for name in dirs:
                fp = os.path.join(root, name)
                for i in CLASSIFIER_FILES:
                    if i in fp:
                        shutil.rmtree('fp', ignore_errors=True)

        ## run amplicon classifier on AA_outputs:

        ## 1. make inputs
        os.system(f"$AC_SRC/make_input.sh {OUTPUT_PATH} {OUTPUT_PATH}/{self.output_name}" )

        ## if reference isn't downloaded already, then download appropriate reference genome
        if not os.path.exists(os.path.join(os.environ['AA_DATA_REPO'], self.ref)):
            os.system(f"wget -q -P $AA_DATA_REPO https://datasets.genepattern.org/data/module_support_files/AmpliconArchitect/{self.ref}.tar.gz")
            os.system(f"wget -q -P $AA_DATA_REPO https://datasets.genepattern.org/data/module_support_files/AmpliconArchitect/{self.ref}_indexed_md5sum.tar.gz")
            os.system(f"tar zxf $AA_DATA_REPO/{self.ref}.tar.gz --directory $AA_DATA_REPO")

        ## 2. run amplicon classifier.py
        os.system(f"python3 $AC_SRC/amplicon_classifier.py -i {OUTPUT_PATH}/{self.output_name}.input --ref {self.ref} -o {OUTPUT_PATH}/{self.output_name} ")
        
        ## 3. run make_results_table.py
        os.system(f"python3 $AC_SRC/make_results_table.py --input {OUTPUT_PATH}/{self.output_name}.input \
                  --classification_file {OUTPUT_PATH}/{self.output_name}_amplicon_classification_profiles.tsv \
                    --summary_map {OUTPUT_PATH}/{self.output_name}_summary_map.txt \
                        --ref {self.ref}")


        print('************ ENDING AMPLICON CLASSIFIER NOW **************')



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
        aggregate.to_csv('./results/aggregated_results.csv')
        aggregate.to_html('./results/aggregated_results.html')
        with open('./results/run.json', 'w') as run_file:
            json.dump({'runs': runs}, run_file)

        run_file.close()
        return sample_to_ac_dir

    def move_files(self):
        """
        Move files to correct location for output
        """
        shutil.move(DEST_ROOT, OUTPUT_PATH)

    def tardir(self, path, tar_name):
        """
        compresses the path to a tar_mname
        """
        with tarfile.open(tar_name, "w:gz") as tar_handle:
            for root, dirs, files in os.walk(path):
                for file in files:
                    fp = os.path.join(root, file)
                    extension = os.path.splitext(fp)[-1]
                    if extension not in EXCLUSION_LIST:
                        tar_handle.add(fp)
        tar_handle.close()

    def cleanup(self):
        """
        Zips the aggregate results, and deletes files for cleanup

        """
        self.clean_cnr_gzs()
        self.clean_files(self.samp_AA_dct.values())
        self.clean_files(self.samp_ckit_dct.values())
        print("Creating tar.gz...")
        self.tardir('./results', f'{self.output_name}.tar.gz')
        shutil.rmtree('./results')

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
    
    def string_to_list(self, string):
        """
        Converts a string representation of a list to a list.
        """
        if "|" in string:
            return string.split("|")
        else:
            return string.strip('][').replace(" ", "").split(',')

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
        for sample in dct['runs'].keys():
            for sample_dct in dct['runs'][sample]:
                # updating string of lists to lists
                sample_name = sample_dct["Sample name"]
                ref_genomes.add(sample_dct["Reference version"])
                if len(ref_genomes) > 1:
                    sys.stderr.write(str(ref_genomes) + "\n")
                    sys.stderr.write("ERROR! Multiple reference genomes detected in project.\n AmpliconRepository only "
                                     "supports single-reference projects currently. Exiting.")
                    sys.exit(1)

                potential_str_lsts = [
                    'Location',
                    'Oncogenes',
                    'All genes',
                ]
                for lists in potential_str_lsts:
                    try:
                        sample_dct[lists] = self.string_to_list(sample_dct[lists])
                    except:
                        print(f'Feature: {lists} doesnt exist for sample: {sample_name}')

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
                        print(f'Feature: "Feature BED file" doesnt exist for sample: {sample_dct["Sample name"]}')
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
                            print(f'Feature: "{feature}" doesn\'t exist for sample: {sample_dct["Sample name"]}')
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
                        tarf = orig_dir + ".tar.gz"
                        self.tardir(orig_dir, tarf)
                        sample_dct[dirfname] = tarf.replace('./results/', "")


                    else:
                        sample_dct[dirfname] = "Not Provided"

        with open('./results/run.json', 'w') as json_file:
            json.dump(dct, json_file, sort_keys=True, indent=2)
        json_file.close()


    def clean_cnr_gzs(self):
        print("Removing raw copy number cnr.gz files to save space.")
        cmd = f'find {OUTPUT_PATH} -name "*.cnr.gz" -exec rm {{}} \;'
        print(cmd)
        subprocess.call(cmd, shell=True)


    def clean_files(self, flist):
        for f in flist:
            if not f == "/":
                cmd = f'rm -rf {f}'
                subprocess.call(cmd, shell=True)


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
