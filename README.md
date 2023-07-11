# Amplicon Suite Aggregator

## Description
Aggregates the results from AmpliconSuite
  - Takes in zip files (completed results of individual or grouped Amplicon Suite runs)
  - Aggregates results
    - Packages results into a new file
    - Outputs an aggregated .html and .csv file of results. 
  - Result file (in .tar.gz) is the aggregated results of all individual AmpliconSutie runs. It can be directly loaded onto AmpliconRepository.
  - Can also take additional files along with the upload, provided the directory they are in contains a file named `AUX_DIR`.
    
## Parameters
  - *Amplicon Architect Results* (required)
    - Compressed (.tar.gz, .zip) files of the results from individual or grouped Amplicon Architect runs.
  - *project_name* (required)
    - Prefix for output .tar.gz. Result will be named: *output_prefix*.tar.gz
  - *Amplicon Repository Email*
    - If wanting to directly transfer results to AmpliconRepository.org, please enter your email.
  - *run amplicon classifier*
    - Option for users to re-run Amplicon Classifier.
  - *reference genome*
    - Reference genome used for Amplicon Architect results in the input. 

## Docker
  - Uses Docker image: genepattern/amplicon-suite-aggregator
  - To send results to Amplicon Repository directly from your local machine, follow these steps:
      1. Create an `input_list.txt` file in this format where each line is a filepath to compressed aa_results:
         <img width="269" alt="image" src="https://github.com/genepattern/AmpliconSuiteAggregator/assets/43209173/cca8c8af-a58f-4290-a8a8-4dbe82e5941e">


      2. `docker pull genepattern/amplicon-suite-aggregator`
      3. `docker run --rm -it -v PATH/TO/INPUTS/FOLDER:/inputs/ genepattern/amplicon-suite-aggregator:v2.3 python3 /opt/genepatt/AmpliconSuiteAggregator.py -flist /path/to/input_list.txt -u YOUR_AMPLICON_REPOSITORY_EMAIL -t`
      4. Log into Amplicon Repository, you should see a new project with the name "Job_*****". 

## Programming Language
  - Python
