# Amplicon Suite Aggregator

## Description
Aggregates results from [AmpliconSuite-pipeline](https://github.com/AmpliconSuite/AmpliconSuite-pipeline) runs.
- Accepts compressed input files (`.tar.gz`, `.zip`) of individual or grouped AmpliconSuite results, or plain directories.
- Aggregates and packages results into a new `.tar.gz` file, along with an aggregated `.csv` and `.html` summary.
- The output `.tar.gz` can be directly uploaded to [AmpliconRepository.org](https://AmpliconRepository.org).
- Supports batch sample renaming via a name map file.
- Can include auxiliary files alongside the upload; place them in a directory containing a file named `AUX_DIR`.

## Installation

**Option 1 — pip**
```bash
pip install AmpliconSuiteAggregator
```

**Option 2 — Git clone**
```bash
git clone https://github.com/AmpliconSuite/AmpliconSuiteAggregator.git
cd AmpliconSuiteAggregator
pip install -r requirements.txt
```

## Dependencies
Python packages: `pandas`, `requests`

## Usage

```bash
python src/AmpliconSuiteAggregator.py -flist <input_list.txt> -o <project_name> [options]
```

`input_list.txt` is a plain text file with one input path per line (`.tar.gz`, `.zip`, or directory).

### Aggregation options

| Flag | Description |
|---|---|
| `-flist FILE` | Text file listing input paths, one per line |
| `--files PATH [PATH ...]` | Input files or directories directly on the command line |
| `-o NAME` | Output prefix / project name (required) |
| `--name_map FILE` | Two-column file: col 1 = current sample name, col 2 = replacement name. Applies a deep rename throughout all output files and tables. |
| `--exclude_samples NAME [NAME ...]` | Sample names to leave out of the output entirely. Matched against the names the samples arrive under, i.e. before `--name_map` is applied. |
| `--no_cleanup` | Keep the temporary working directories after finishing (for debugging) |
| `-v`, `--version` | Print the version and exit |

Aggregation no longer uploads anything, and no longer runs AmpliconClassifier
for you — run AC first, then aggregate its output. The `-u`/`-s`/`--upload_only`
upload flags and the `-c`/`--ref` classifier flags were removed in 8.0.0.

## Examples

**Aggregate a set of results:**
```bash
python src/AmpliconSuiteAggregator.py -flist input_list.txt -o MyProject
```

**Rename samples on the way through:**
```bash
python src/AmpliconSuiteAggregator.py -flist input_list.txt -o MyProject --name_map renames.tsv
```

**Drop samples from a project you already aggregated:**
```bash
python src/AmpliconSuiteAggregator.py --files MyProject.tar.gz -o MyProject_v2 \
    --exclude_samples SampleA SampleB
```
An aggregated archive is itself valid input, so this is how you remove a sample
from an existing project: everything belonging to `SampleA`/`SampleB` is left
out — their directories under `results/samples/`, their rows in `run.json`,
`aggregated_results.csv` and every consolidated classification table, and their
per-amplicon files. A name matching no sample is reported as a warning rather
than failing the run.

## Input files

For a breakdown of which per-sample files make up a classification project —
which are essential, which are strongly recommended, and the naming patterns
the Aggregator expects — see [INPUT_FILES.md](INPUT_FILES.md).

## Authors
- [Jens Luebeck](https://github.com/jluebeck) *(main contact)*
- Thorin Tabor
- Edwin Huang

## Issues
Bug reports and feature requests: [GitHub Issues](https://github.com/AmpliconSuite/AmpliconSuiteAggregator/issues)
