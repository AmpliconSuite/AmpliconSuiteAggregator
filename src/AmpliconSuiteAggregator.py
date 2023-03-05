from AmpliconSuiteAggregatorFunctions import *



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-flist", "--filelist", type = str,
            help = "List of zip files to use")

    parser.add_argument("--output_name", type = str, help = "Output Prefix", default = "aggregated")

    args = parser.parse_args()


    filelist = args.filelist
    root = '.'
    aggregate = Aggregator(filelist, root, args.output_name)