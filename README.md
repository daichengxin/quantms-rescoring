# quantms-rescoring
    
[![Python package](https://github.com/bigbio/quantms-rescoring/actions/workflows/python-package.yml/badge.svg)](https://github.com/bigbio/quantms-rescoring/actions/workflows/python-package.yml)
[![codecov](https://codecov.io/gh/bigbio/quantms-rescoring/branch/main/graph/badge.svg?token=3ZQZQ2ZQ2D)](https://codecov.io/gh/bigbio/quantms-rescoring)
[![PyPI version](https://badge.fury.io/py/quantms-rescoring.svg)](https://badge.fury.io/py/quantms-rescoring)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

quantms-rescoring is a Python tool that aims to add features to peptide-spectrum matches (PSMs) in idXML files using multiple tools including SAGE features, quantms spectrum features, MS2PIP and DeepLC. It is part of the quantms ecosystem package and leverages the MS²Rescore framework to improve identification confidence in proteomics data analysis.

### Core Components

- **Annotator Engine**: Integrates [MS2PIP](https://github.com/compomics/ms2pip) and [DeepLC](https://github.com/compomics/DeepLC) models to improve peptide-spectrum match (PSM) confidence. 
- **Feature Generation**: Extracts signal-to-noise ratios, spectrum metrics, SAGE extra features and add them to each PSM for posterior downstream with Percolator.
- **OpenMS Integration**: Processes idXML and mzML files with custom validation methods.

### CLI Tools

```sh
 quantms-rescoring msrescore2feature --help
```
Annotates PSMs with prediction-based features from MS2PIP and DeepLC

```sh
 quantms-rescoring add_sage_feature --help
```
Incorporates additional features from SAGE into idXML files. 

```sh
 quantms-rescoring spectrum2feature --help
```
Add additional spectrum feature like signal-to-noise to each PSM in the idXML.

### Technical Implementation Details

#### Model Selection and Optimization

- **MS2PIP Model Selection**: 
  - Automatically evaluate the quality of the MS2PIP model selected by the user. If the correlation between predicted and experiemtanl spectra is lower than a given threshold, we will try to find the best model to use (`annotator.py`)
- **DeepLC Model Selection**: 
  - Automatically select the best DeepLC model for each run based on the retention time calibration and prediction accuracy. Different to ms2rescore, the tool will try to use the best model from MS2PIP and benchmark it with the same model by using transfer learning (`annotator.py`). The best model is selected to be used to predict the retention time of PSMs.

#### Feature Engineering Pipeline

- **Retention Time Analysis**:
  - Calibrates DeepLC models per run to account for chromatographic variations.
  - Calculates delta RT (predicted vs. observed) as a discriminative feature
  - Normalizes RT differences for cross-run comparability

- **Spectral Feature Extraction**:
  - Computes signal-to-noise ratio using maximum intensity relative to background noise
  - Calculates spectral entropy to quantify peak distribution uniformity
  - Analyzes TIC (Total Ion Current) distribution across peaks for quality assessment
  - Determines weighted standard deviation of m/z values for spectral complexity estimation
- **Feature Selection**: The parameters `only_features` allows to select the features to be added to the idXML file. For example: `--only_features "DeepLC:RtDiff,DeepLC:PredictedRetentionTimeBest,Ms2pip:DotProd"` . 

#### Data Processing of idXML Files

- **Parallel Processing**: Implements multiprocessing capabilities for handling large datasets efficiently
- **OpenMS Compatibility Layer**: Custom helper classes that gather statistics of number of PSMs by MS levels / dissociation methods, etc.
- **Feature Validation**: Convert all Features from MS2PIP, DeepLC, and quantms into OpenMS features with well-established names (`constants.py`)
- **PSM Filtering and Validation**: 
  - Filter PSMs with **missing spectra information** or **empty peaks**.
  - Breaks the analysis of the input file contains more than one MS level or dissociation method, **only support for MS2 level** spectra. 
- **Output / Input files**: 
  - Only works for OpenMS formats idXML, and mzML as input and export to idXML with the annotated features. 

### Installation

Install quantms-rescoring using one of the following methods:

**Using `pip`**

```sh
❯ pip install quantms-rescoring
```

**Using `conda`** 

```sh
❯ conda install -c bioconda quantms-rescoring
```

**Build from source:**

1. Clone the quantms-rescoring repository:

   ```sh
   ❯ git clone https://github.com/bigbio/quantms-rescoring
   ```

2. Navigate to the project directory:

   ```sh
   ❯ cd quantms-rescoring
   ```

3. Install the project dependencies:

   - Using `pip`:

     ```sh
     ❯ pip install -r requirements.txt
     ```

   - Using `conda`:

     ```sh
     ❯ conda env create -f environment.yml
     ```
  
4. Install the package using `poetry`:

   ```sh
   ❯ poetry install
   ```

### TODO

- [ ] Add support for multiple Files combined idXML and mzML

### Issues and Contributions

For any issues or contributions, please open an issue in the [GitHub repository](https://github.com/bigbio/quantms/issues) - we use the quantms repo to control all issues—or PR in the [GitHub repository](https://github.com/bigbio/quantms-rescoring/pulls). 

