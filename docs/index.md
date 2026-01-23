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

## ReMASTER to BEAST2 XML Workflow

This package provides functionality to generate BEAST2 XML configuration files by combining a template XML with simulation data from ReMASTER (Nexus alignment and tree).

### XML Template Requirements
Your input template XML must be a valid BEAST2 XML file (except for the placeholders) and **must** contain the following three specific placeholders, which will be replaced by the extracted data:

1.  `INSERTSEQUENCES`: Replaced with the `<sequence>` blocks inside the alignment `<data>` element.
2.  `INSERTTRAITDATES`: Replaced with a comma-separated string of `taxon=date` pairs (e.g., `taxonA=2020/01/01,taxonB=2020/02/15`).
3.  `INSERTTRAITTYPES`: Replaced with a comma-separated string of `taxon=type` pairs (e.g., `taxonA=location1,taxonB=location2`).

**Example Template Structure:**
```xml
<data id="alignment" spec="Alignment">
    INSERTSEQUENCES
</data>

<trait id="dateTrait" spec="beast.base.evolution.tree.TraitSet" traitname="date" value="INSERTTRAITDATES">
    ...
</trait>

<trait id="typeTrait" spec="trait.TraitSet" traitname="type" value="INSERTTRAITTYPES">
    ...
</trait>
```

### Time and Date Interpretation
ReMASTER simulations typically output trees with a `time` attribute for each node. In this package, this `time` attribute is interpreted as the **time (in years) passed since the start of the simulation**.

To generate valid calendar dates for BEAST2:
1.  You provide an **artificial start date** (e.g., "2000/01/01").
2.  The package reads the `time` (years) from the ReMASTER tree for each leaf.
3.  It calculates the leaf date as: `StartDate + time (converted to days)`.
4.  These calculated dates are inserted into the `INSERTTRAITDATES` placeholder in the format `YYYY/MM/DD`.

*Note: The `INSERTTRAITDATES` replacement string assumes the template expects a date format of `YYYY/MM/DD`.*
