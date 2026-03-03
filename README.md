
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
git clone <repository-url>
cd CACalibrate
```

## Usage

IN an R interactive prompt, call the script with
```
source("/path/to/cloned/repository/align_and_merge.R", encoding = "UTF-8")
```

## Contributing

Contributions are welcome. Please submit pull requests with clear descriptions.

## License

MIT License
