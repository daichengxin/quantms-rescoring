# quantms-rescoring
    
[![Python package](https://github.com/bigbio/quantms-rescoring/actions/workflows/python-package.yml/badge.svg)](https://github.com/bigbio/quantms-rescoring/actions/workflows/python-package.yml)
[![codecov](https://codecov.io/gh/bigbio/quantms-rescoring/branch/main/graph/badge.svg?token=3ZQZQ2ZQ2D)](https://codecov.io/gh/bigbio/quantms-rescoring)
[![PyPI version](https://badge.fury.io/py/quantms-rescoring.svg)](https://badge.fury.io/py/quantms-rescoring)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

quantms-rescoring is a Python tool for rescoring peptide-spectrum matches (PSMs) in idXML files. It is part of the quantms ecosystem package and leverages the MS²Rescore framework to improve identification confidence in proteomics data analysis.

## Features

- Enhanced Rescoring: Utilizes advanced rescoring engines like Percolator to refine PSM scores.
- Flexible Feature Generators: Supports feature extraction using tools like MS²PIP, DeepLC, and custom generators.
- Metadata Retention: Preserves essential metadata from the input idXML files.
- Error Handling: Skips invalid PSMs and logs issues for transparent processing.
- Seamless Integration: Built to integrate into proteomics workflows.

## Installation

To use quantms-rescoring, ensure the following dependencies are installed:

- Python 3.8+
- click
- pyopenms
- ms2rescore
- psm_utils

### Issues and Contributions

For any issues or contributions, please open an issue in the [GitHub repository](https://github.com/bigbio/quantms/issues) - we use the quantms repo to control all issues—or PR in the [GitHub repository](https://github.com/bigbio/quantms-rescoring/pulls). 

