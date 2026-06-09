
# CACalibrate

A calibration toolkit for Clover Array data.

## Overview

A set of R scripts for ingesting Cubix energy calibration outputs, aligning non-reference runs to a reference run, fitting a linear calibration, and correcting residual structure with a spline model.

## Dependencies

- R
- Cubix (the scripts expect Cubix calibration output files)
- R packages: ggplot2 and mgcv


## Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/tmkowalewski/CACalibrate.git
cd CACalibrate
```

## Usage
### 1) Confirm filename pattern parsing

The parser in energy_fit.r expects filenames shaped like:

RUN_TYPE_DETECTORCRYSTAL_calibration.ext

Example:

70Ge_11B-7.28_C1E1_calibration.dat

with:

- RUN_TYPE: 70Ge_11B-7.28
- DETECTOR: C1
- CRYSTAL: E1
- ext: dat, func, or res

If your files use a different naming pattern, update parse_filename_meta in energy_fit.r.

### 2) Configure align_and_merge.R

Set the Configuration block in align_and_merge.R:

- cubix_files: root directory containing your Cubix workspaces
- outdir: output directory for calibration files and plots
- reference_run_type: the parsed run_type used as the alignment reference (for example 70Ge_AllCal)
- match_run_type_pattern: regex selecting which non-reference run types to align
- gain_anchor_energies and gain_anchor_window: fixed gain-alignment anchors and matching window

Important:

- align_and_merge.R wipes the contents of outdir at startup before writing new results.

### 3) Run

From the repository directory:

```bash
Rscript align_and_merge.R
```

### 4) Outputs

For each detector/crystal pair, the script writes files in outdir including:

- DETECTORECRYSTAL.cal_params.txt
- DETECTORECRYSTAL.models.rds
- DETECTORECRYSTAL.lin_res.png
- DETECTORECRYSTAL.spline_res.png (when spline residual output is available)
- DETECTORECRYSTAL.alignment_diagnostics.tsv

The diagnostics table records, per run type, whether it was added, skipped, excluded by pattern, or missing the reference, along with match counts and reasons.

## Contributing

Contributions are welcome. Please submit pull requests with clear descriptions.

## License

MIT License
