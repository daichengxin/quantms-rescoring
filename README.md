# quantms-ms2rescore
quantms-ms2rescore is a Python tool for rescoring peptide-spectrum matches (PSMs) in idXML files. It is part of the quantms ecosystem package and leverages the MS²Rescore framework to improve identification confidence in proteomics data analysis.

## Features

- Enhanced Rescoring: Utilizes advanced rescoring engines like Percolator to refine PSM scores.
- Flexible Feature Generators: Supports feature extraction using tools like MS²PIP, DeepLC, and custom generators.
- Metadata Retention: Preserves essential metadata from the input idXML files.
- Error Handling: Skips invalid PSMs and logs issues for transparent processing.
- Seamless Integration: Built to integrate into proteomics workflows.

## Installation
To use quantms-ms2rescore, ensure the following dependencies are installed:

- Python 3.8+
- click
- pyopenms
- ms2rescore
- psm_utils

