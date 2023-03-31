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

DEST_ROOT = os.path.join("./extracted_from_zips")
OUTPUT_PATH = os.path.join("./results/AA_outputs")
OTHER_FILES = os.path.join("./results/other_files")

class Aggregator():
    def __init__(self, filelist, root, output_name):
        self.filelist_fp = filelist
        self.root = root
        self.output_name = output_name
        self.get_zip_paths()
        self.unzip()
        self.sample_to_ac_location_dct = self.aggregate_tables()
        self.json_modifications()
        self.move_files()
        self.cleanup()



    def get_zip_paths(self):
        """
        Gets the individual file paths from the list of zip filepaths

        calls from self
        returns --> a list of filepaths to zip files
        """
        files = []
        with open(self.filelist_fp, 'r') as filelist:
            for line in filelist:
                parsed = line.strip()
                files.append(parsed)
        filelist.close()
        self.zip_paths = files
        return files

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
                        os.system(f'mv -vf {os.path.dirname(fp)} {OUTPUT_PATH}')

                elif any([x.endswith("_result_table.tsv") or x == "AUX_DIR" for x in os.listdir(fp)]):
                    pre = fp.rstrip("_classification")
                    if not fp.endswith("_classification") and not os.path.exists(pre + "_AA_results"):
                        print(f'Moving file {fp} to {OTHER_FILES}')
                        os.system(f'mv -vf {fp} {OTHER_FILES}')

                    # try: # check if this is an AA outputs folder
                    #     sample_name = re.findall(f'{DEST_ROOT}\/.*\/(.*)_AA_results$', fp)[0]
                    #     samples.append(sample_name)
                    #     print(os.path.dirname(fp))
                    #     print(os.path.exists(os.path.dirname(fp)))
                    #     os.system(f'mv -vf {os.path.dirname(fp)} {OUTPUT_PATH}')
                    #     # shutil.move(os.path.dirname(fp), output_dir)
                    # except:
                    #     continue

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
                        print(result_table_fp)
                        try:
                            df = pd.read_csv(result_table_fp, delimiter = '\t')
                            aggregate = pd.concat([aggregate, df], ignore_index = True)

                        except Exception as e:
                            print(e)
                            continue

                        # print(df)
                        gdf = df.groupby(['Sample name'])
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
                    tar_handle.add(os.path.join(root, file))
        tar_handle.close()

    def cleanup(self):
        """
        Zips the aggregate results, and deletes files for cleanup

        """

        self.tardir('./results', f'{self.output_name}.tar.gz')
        shutil.rmtree('./results')

    def find_file(self, basename):
        """
        Finds a specific file given sample name and file suffix
        """
        for root, dirs, files in os.walk('./results/AA_outputs', topdown = True):
            for file in files:
                fp = os.path.join(root, file)
                if basename == file:
                    return fp.replace('./results/', "")

        return 'Not Provided'
    
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

        for sample in dct['runs'].keys():
            for sample_dct in dct['runs'][sample]:
                # updating string of lists to lists
                sample_name = sample_dct["Sample name"]
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
                    'Sample metadata JSON'
                ]
                for feature in features_of_interest:
                    if feature in sample_dct and sample_dct[feature]:
                        feat_basename = os.path.basename(sample_dct[feature])
                        feat_file = f'{self.sample_to_ac_location_dct[sample]}/files/{feat_basename}'
                        if os.path.exists(feat_file):
                            feat_file = feat_file.replace('./results/', "")
                            sample_dct[feature] = feat_file

                        else:
                            print(f'Feature: "{feature}" doesnt exist for sample: {sample_dct["Sample name"]}')
                            sample_dct[feature] = "Not Provided"

                    else:
                        sample_dct[feature] = "Not Provided"

        with open('./results/run.json', 'w') as json_file:
            json.dump(dct, json_file, sort_keys=True, indent=2)
        json_file.close()

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

