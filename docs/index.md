# Welcome to beast2analysisutils

A miscellaneous collection of analysis utilities that I frequently use in the BEAST2 ecosystem.

Broadly there are two types of analysis currently in this package:

1. **ReMASTER XML Generation**: Construct BEAST2 XML files from templates by inserting simulated trees and/or sequences from ReMASTER output.
2. **ESS Analysis**: Calculate the effective sample sizes (ESS) for all parameters in a BEAST2 log file.

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

## Quick Start

See the [Reference](reference.md) page for full API details.

### 1. BEAST2 XML Generation from ReMASTER Output

This module helps you turn ReMASTER simulation outputs (Nexus alignment + Tree) into runnable BEAST2 XMLs given a template XML.

#### Step 1: Construct Templates
You will need to construct your own BEAST2 templates based on the analysis that you would like to do. If you want to see two examples, the package provides two templates for fixed and estimated tree analyses. You can download them directly from GitHub:

- **Alignment based Template**: [Template.xml](https://raw.githubusercontent.com/Pweidemueller/beast2analysisutils/main/src/beast2analysisutils/data/templates/Template.xml)
- **Fixed Tree based Template**: [Template_fixedtree.xml](https://raw.githubusercontent.com/Pweidemueller/beast2analysisutils/main/src/beast2analysisutils/data/templates/Template_fixedtree.xml)

#### Step 2: Generate XML
Use the `generate_xml` wrapper function.

**Scenario A: Standard Analysis (Tree will be estimated)**
Requires both a Nexus alignment file and a Tree file from ReMASTER.

```python
from beast2analysisutils.remaster import generate_xml

generate_xml(
    template_path="./Template.xml",
    tree_path="simulation.trees",      # ReMASTER tree output
    output_path="output.xml",
    alignment_path="simulation.nexus", # ReMASTER alignment output
    fixed_tree=False                   # Default
)
```

**Scenario B: Fixed Tree Analysis**
Uses the simulated tree directly in the XML. Nexus file is optional (dummy sequences will be generated if missing).

```python
from beast2analysisutils.remaster import generate_xml

generate_xml(
    template_path="./Template_fixedtree.xml",
    tree_path="simulation.trees",
    output_path="output_fixed.xml",
    alignment_path=None,               # Optional for fixed tree
    fixed_tree=True
)
```

### 2. ESS Analysis

Calculate ESS for all parameters in a BEAST2 log file to check for convergence.

```python
from beast2analysisutils.ess import analyze_ess

# Calculate ESS for all numeric columns
ess_df = analyze_ess(
    log_source="beast_output.log",
    output_path="ess_results.csv",
    burnin=0.1,             # 10% burn-in
    check_threshold=True    # Check if ESS > 200 for key parameters
)

print(ess_df.head())
```

## Detailed Documentation

### ReMASTER to BEAST2 XML Workflow

This package provides functionality to generate BEAST2 XML configuration files by combining a template XML with simulation data.

#### Wrapper Function: `generate_xml`
The `generate_xml` function is the main entry point. It performs the following steps:
1.  **Validates Inputs**: Checks for required files based on `fixed_tree` mode.
2.  **Extracts Data**: Parses sequences, leaf dates, and types from input files.
3.  **Consistency Check**: Verifies that the number of sequences matches the number of tree leaves.
4.  **Prints Stats**: Displays sequence counts, leaf counts, trait distribution, and date ranges.
5.  **Generates XML**: Fills the template with the processed data.

#### XML Template Requirements
Your input template XML must be a valid BEAST2 XML file containing specific placeholders:

1.  `INSERTSEQUENCES`: Replaced with `<sequence>` blocks.
2.  `INSERTTRAITDATES`: (not needed for fixed tree) Replaced with `taxon=date` pairs (YYYY/MM/DD).
3.  `INSERTTRAITTYPES`: Replaced with `taxon=type` pairs.
4.  `INSERTNEWICKTREE`: (Fixed tree only) Replaced with the Newick tree string.

**Example Template Structure:**
```xml
<data id="alignment" spec="Alignment">
    INSERTSEQUENCES
</data>

<trait id="dateTrait" spec="beast.base.evolution.tree.TraitSet" traitname="date" value="INSERTTRAITDATES">
    ...
</trait>

<trait id="typeTrait" spec="mascot.util.InitializedTraitSet" traitname="type" value="INSERTTRAITTYPES">
    ...
</trait>

<!-- For Fixed Tree templates -->
<init id="NewickTree" spec="beast.base.evolution.tree.TreeParser" newick="INSERTNEWICKTREE" .../>
```

### Time and Date Interpretation
ReMASTER simulations typically output trees with a `time` attribute (years since simulation start).
The package converts these to calendar dates using a reference `start_date` (default "2000/01/01").
leaf_date = `start_date` + `time` (converted to days).
