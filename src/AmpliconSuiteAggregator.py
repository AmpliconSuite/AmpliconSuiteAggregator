from AmpliconSuiteAggregatorFunctions import *
from ASA_POST import * 
import numpy as np


def get_zip_paths(filelist_fp):
    """
    Gets the individual file paths from the list of zip filepaths

    calls from self
    returns --> a list of filepaths to zip files
    """
    files = []
    with open(filelist_fp) as filelist:
        for line in filelist:
            parsed = line.strip()
            if parsed:
                files.append(parsed)

    return files


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-flist", "--filelist", type=str, help="Text file with files to use (one per line)")
    group.add_argument("--files", nargs='+', type=str, help="List of files or directories to use")
    parser.add_argument("-o", "--output_name", type=str, help="Output Prefix", default="aggregated")
    parser.add_argument("-u", "--username", type = str, help = "Username for Amplicon Repository")
    parser.add_argument('-t', '--testing', action = 'store_true')



    args = parser.parse_args()

    if args.filelist:
        filelist = get_zip_paths(args.filelist)

    else:
        filelist = args.files

    root = '.'
    aggregate = Aggregator(filelist, root, args.output_name)
    output_fp = args.output_name + ".tar.gz"

    if args.username:
        current_path = os.getcwd()
        job_number = current_path.split('/')[-1]
        ## get the job number of the commandline, something like (from genepattern job number)
        ## get pwd and get the job number from there. 

        ## testrun: 
        # python3 /files/src/AmpliconSuiteAggregator.py -flist /files/gpunit/inputs/input_list.txt -u edh021@ucsd.edu

        user = args.username
        if args.testing:
            job_number = np.random.randint(0, 99999)
        data = {'project_name': f'Job_{job_number}',
                'description': f'Results transferred from GenePattern, job id: {job_number}',
                'private': True,
                'project_members': [args.username],}
        
        print(f'Creating project for user: {user}')
        post_package(output_fp, data, user)