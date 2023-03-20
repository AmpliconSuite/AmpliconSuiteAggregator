# Amplicon Suite Aggregator

## Description
Aggregates the results from AmpliconSuite
  - Takes in zip files (completed results of individual Amplicon Suite runs)
  - Aggregates results
    - Packages results into a new file
    - Outputs an aggregated .html and .csv file of results. 
  - Result file (in .tar.gz) is the aggregated results of all individual AmpliconSutie runs. It can be directly loaded onto AmpliconRepository.
    
## Parameters
  - *Amplicon Architect Results* (required)
    - Compressed (.tar.gz, .zip) files of the results from individual or grouped Amplicon Architect runs.
  - *output_preifx* (required)
    - Prefix for output .tar.gz. Result will be named: *output_prefix*.tar.gz

## Docker
  - Uses Docker image: genepattern/amplicon-suite-aggregator:v2.1

## Programming Language
  - Python
