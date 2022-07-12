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


class Aggregator():
    def __init__(self, filelist, root):
        self.filelist_fp = filelist
        self.root = root
        self.get_zip_paths()
        self.unzip()
        self.aggregate()
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



    def aggregate(self):
        """
        From the unzipped paths, start aggregating results

        """
        aggregate = pd.DataFrame(columns = ["Sample name", "AA amplicon number", "Feature ID", "Classification", "Location", "Oncogenes",
                   "Complexity score", "Feature BED file", "AA PNG file", "AA PDF file"])


        # aggregate results
        for root, dirs, files in os.walk("extracted", topdown = False):
            for name in files:
                if "_result_table.tsv" in name:
                    result_table_fp = os.path.join(root, name)
                    print(result_table_fp)
                    df = pd.read_csv(result_table_fp, delimiter = '\t')
                    aggregate = pd.concat([aggregate, df], ignore_index = True)
                    print(os.path.join(root, name))
        ## output the table
        aggregate.to_csv('aggregated_results.csv')
        aggregate.to_html('aggregated_results.html')


    def tardir(self, path, tar_name):
        """
        compresses the path to a tar_mname
        """
        with tarfile.open(tar_name, "w:gz") as tar_handle:
            for root, dirs, files in os.walk(path):
                for file in files:
                    tar_handle.add(os.path.join(root, file))

    def cleanup(self):
        """
        Zips the aggregate results, and deletes files for cleanup

        """

        self.tardir('./extracted', 'aggregated.tar.gz')

        shutil.rmtree('./extracted')



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-flist", "--filelist", type = str,
            help = "List of zip files to use")
    args = parser.parse_args()

    filelist = args.filelist
    root = '.'
    aggregate = Aggregator(filelist, root)
