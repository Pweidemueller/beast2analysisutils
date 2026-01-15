# Welcome to beast2analysisutils

A miscellaneous collection of analysis utilities that I frequently use in the BEAST2 ecosystem.

Broadly there are two types of analysis currently in this package:

1. Using ReMASTER output to construct BEAST2 XML files from templates by inserting simulated trees and/or sequences.
2. Calculating the effective sample sizes (ESS) for all parameters in a BEAST2 log file.

## Installation

If you already have a python environment set up, you can install this package via `pip` directly from GitHub.

**1. Conda Environment**
```bash
conda activate my_analysis_env
pip install git+https://github.com/Pweidemueller/beast2analysisutils.git
```

**2. Standard pip/venv**
```bash
source .venv/bin/activate
pip install git+https://github.com/Pweidemueller/beast2analysisutils.git
```

**3. uv Environment**
```bash
uv pip install git+https://github.com/Pweidemueller/beast2analysisutils.git
```

## Usage

See the [Reference](reference.md) page for API details.
