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




class Aggregator():
    def __init__(self, filelist, root, output_name):
        self.filelist_fp = filelist
        self.root = root
        self.output_name = output_name
        print(self.get_zip_paths())
        self.unzip()
        self.aggregate_tables()
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


    def unzip(self):
        """
        Unzips the zip files, and get directories for files within

        returns -- None
        """

        # os.mkdir('all_AA_results')

        for zip_fp in self.zip_paths:
            fp = os.path.join(self.root, zip_fp)
            try:
                with tarfile.open(fp, 'r') as output_zip:
                    output_zip.extractall('./extracted')
                output_zip.close()
            except Exception as e:
                print("The zip: " + fp + " ran into an error when unzipping.")
                print(e)

        samples = []
        # parse contents of zip to get sample_name
        for root, dirs, files in os.walk("extracted", topdown=False):
            sample_name = ""
            for name in dirs:
                dir_name = os.path.join(root, name)
                if "_cnvkit_output" in dir_name:
                    dir_name = dir_name.split("/")[-1]
                    sample_name = dir_name.replace("_cnvkit_output", "")
                    sample_name = sample_name.replace('extracted/', "")
                    samples.append(sample_name)
        self.samples = samples
        print(self.samples)



    def aggregate_tables(self):
        """
        From the unzipped paths, start aggregating results

        """

        output_head = ["Sample name", "AA amplicon number", "Feature ID", "Classification", "Location", "Oncogenes",
                           "Complexity score",
                           "Captured interval length", "Feature median copy number", "Feature maximum copy number",
                           "Reference version", "Feature BED file", "CNV BED file", "AA PNG file", "AA PDF file",
                           "Run metadata JSON"]
        output_table_lines = [output_head, ]
        aggregate = pd.DataFrame(columns = output_head)
        cycle_files = []
        graph_files = []
        runs = {}
        sample_num = 1
        # aggregate results
        for root, dirs, files in os.walk("extracted", topdown = False):
            for name in files:
                if "_result_table.tsv" in name:
                    result_table_fp = os.path.join(root, name)
                    df = pd.read_csv(result_table_fp, delimiter = '\t')
                    aggregate = pd.concat([aggregate, df], ignore_index = True)

                    ## convert to json
                    runs[f'sample_{sample_num}'] = json.loads(df.to_json(orient = 'records'))
                    sample_num += 1

            for dir in dirs:
                if "files" == dir:
                    dir_fp = os.path.join(root, dir)
                    os.system(f'cp {dir_fp}/* /opt/genepatt/files')



        ## output the table
        aggregate.to_csv('aggregated_results.csv')
        aggregate.to_html('aggregated_results.html')
        with open('run.json', 'w') as run_file:
            json.dump({'runs': runs}, run_file)
        run_file.close()

    def aggregate_json(self, json_files):
        """
        Aggregates json files from separate samples

        """
        out_dict = {'run' : {}}
        sample_ct = 1
        for fp in json_files:
            sample_dict = json.load(open(fp))
            out_dict['run'][f'sample_{sample_ct}'] = sample_dict
            sample_ct += 1

        with open('run.json', 'w') as out_json:
            json.dump(out_dict, out_json)
        out_json.close()

    def move_files(self):
        """
        Move files to correct location for output
        """

        shutil.move('/opt/genepatt/files', '/opt/genepatt/extracted/files')
        shutil.copy('run.json', '/opt/genepatt/extracted/run.json')
        shutil.copy('aggregated_results.csv', '/opt/genepatt/extracted/aggregated_results.csv')
        shutil.copy('aggregated_results.html', '/opt/genepatt/extracted/aggregated_results.html')


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
        self.tardir('extracted', f'{self.output_name}.tar.gz')
        shutil.rmtree('extracted')



if __name__ == "__main__":
    # os.environ["AA_DATA_REPO"] = os.path.join('/opt/genepatt', '.data_repo')

    parser = argparse.ArgumentParser()
    parser.add_argument("-flist", "--filelist", type = str,
            help = "List of zip files to use")

    parser.add_argument("--output_name", type = str, help = "Output Prefix")

    args = parser.parse_args()


    # os.system(f'bash download_ref.sh {REFERENCE}')


    filelist = args.filelist
    root = '.'
    aggregate = Aggregator(filelist, root, args.output_name)
