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
| `-c {Yes,No}` | Re-run Amplicon Classifier on inputs (`Yes`/`No`) |
| `--ref GENOME` | Reference genome: `hg19`, `GRCh37`, `GRCh38`, `GRCh38_viral`, or `mm10` |

### AmpliconRepository upload options

| Flag | Description |
|---|---|
| `-u EMAIL` | AmpliconRepository username (email). If provided, triggers upload after aggregation. |
| `--upload_only {Yes,No}` | Skip aggregation and upload an existing `.tar.gz` directly |
| `-s {prod,dev}` | Target server (`prod` for most users) |

## Examples

**Aggregate a set of results:**
```bash
python src/AmpliconSuiteAggregator.py -flist input_list.txt -o MyProject
```

**Aggregate and upload to AmpliconRepository:**
```bash
python src/AmpliconSuiteAggregator.py -flist input_list.txt -o MyProject -u you@email.com -s prod
```

**Upload an already-aggregated file without re-aggregating:**
```bash
python src/AmpliconSuiteAggregator.py --files MyProject.tar.gz -o MyProject -u you@email.com --upload_only Yes -s prod
```

## Authors
- [Jens Luebeck](https://github.com/jluebeck) *(main contact)*
- Thorin Tabor
- Edwin Huang

## Issues
Bug reports and feature requests: [GitHub Issues](https://github.com/AmpliconSuite/AmpliconSuiteAggregator/issues)
