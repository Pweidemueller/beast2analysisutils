# beast2analysisutils

Analysis utilities for BEAST2.

## Installation

## Installation

### Install into an existing environment

If you already have a python environment set up, you can install this package via `pip`.

**1. Conda Environment**
```bash
conda activate my_analysis_env
# Install from local source (if you cloned the repo)
pip install .

# OR install directly from GitHub
pip install git+https://github.com/Pweidemueller/beast2analysisutils.git
```

**2. Standard pip/venv**
```bash
source .venv/bin/activate
# Install from local source
pip install .

# OR install directly from GitHub
pip install git+https://github.com/Pweidemueller/beast2analysisutils.git
```

**3. uv Environment**
```bash
# Install from local source
uv pip install .

# OR install directly from GitHub
uv pip install git+https://github.com/Pweidemueller/beast2analysisutils.git
```

### Development Setup (from scratch)

To set up a fresh development environment using `uv`:

```bash
# Install uv if you haven't
curl -LsSf https://astral.sh/uv/install.sh | sh

# Initialize environment and install in editable mode with dev dependencies
uv venv
uv pip install -e ".[dev]"
```

## Usage

You can use the utilities in your Python scripts:

```python
from beast2analysisutils import analyze_ess, read_log_file

# Analyze a log file and save results to CSV
df = analyze_ess("path/to/logfile.log", "output.csv")

# Or read explicitly
header, data = read_log_file("path/to/logfile.log")
```

Run tests with:
```bash
uv run pytest
```
