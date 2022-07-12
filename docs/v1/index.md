# Amplicon Suite Aggregator
Author: Edwin Huang <br>
Email: edh021@cloud.ucsd.edu

### Summary: 
  - This module takes in compressed (tar.gz) files from the outputs of individual or grouped Amplicon Suite (Architect) runs, aggregates the outputs, and returns an aggregated table of results as well as a new compressed (tar.gz) file of all input files.
  - Results from: https://github.com/genepattern/AmpliconSuite

### Module Parameters:
  - Amplicon Architect Results: Compressed (tar.gz) files of the results from individual or grouped Amplicon Architect runs. Can take a list of zip files. 

### Outputs:
  - aggregated.tar.gz: A compressed (tar.gz) file of the aggregated inputs
  - aggregated_results.csv: An aggregated table of results
  - aggregated_results.html: Web format for aggregated_results.csv

Module was written in Python3.
