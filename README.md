
# CACalibrate

A calibration toolkit for Clover Array data.

## Overview

A simple R script for taking in Cubix energy calibration data, performing a linear fit, and correcting for non-linearities with a penalized regression spline fit on the residuals. 

## Dependencies

- R (needed to run the script)
- Cubix (script is designed to run on outputs from Cubix's energy calibration tool)


## Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/tmkowalewski/CACalibrate.git
cd CACalibrate
```

## Usage
Modify the regex expression on line 9 of `energy_fit.r` to fit your naming scheme for your cubix calibration data files, as well as the Configuration section of `align_and_merge.R` to point to the right paths and use your preferred fitting options.

To run the script, use
```
source("/path/to/cloned/repository/align_and_merge.R", encoding = "UTF-8")
```
within an R interactive prompt,

## Contributing

Contributions are welcome. Please submit pull requests with clear descriptions.

## License

MIT License
