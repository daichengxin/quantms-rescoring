# quantms-rescoring
    
[![Python package](https://github.com/bigbio/quantms-rescoring/actions/workflows/python-package.yml/badge.svg)](https://github.com/bigbio/quantms-rescoring/actions/workflows/python-package.yml)
[![codecov](https://codecov.io/gh/bigbio/quantms-rescoring/branch/main/graph/badge.svg?token=3ZQZQ2ZQ2D)](https://codecov.io/gh/bigbio/quantms-rescoring)
[![PyPI version](https://badge.fury.io/py/quantms-rescoring.svg)](https://badge.fury.io/py/quantms-rescoring)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

quantms-rescoring is a Python tool that aims to add features to peptide-spectrum matches (PSMs) in idXML files using multiple tools including SAGE features, quantms spectrum features, MS2PIP and DeepLC. It is part of the quantms ecosystem package and leverages the MS²Rescore framework to improve identification confidence in proteomics data analysis.

## Core Components

- **Annotator Engine**: Integrates [MS2PIP](https://github.com/compomics/ms2pip) and [DeepLC](https://github.com/compomics/DeepLC) models to improve peptide-spectrum match (PSM) confidence. 
- **Feature Generation**: Extracts signal-to-noise ratios, spectrum metrics, SAGE extra features and add them to each PSM for posterior downstream with Percolator.
- **OpenMS Integration**: Processes idXML and mzML files with custom validation methods.

## CLI Tools

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

## Advanced Algorithms and Improvements

quantms-rescoring significantly enhances the capabilities of MS2PIP, DeepLC, and MS2Rescore through several innovative approaches:

### MS2PIP Integration Enhancements

- **Intelligent Model Selection**: Automatically evaluates and selects the optimal MS2PIP model for each dataset based on fragmentation type and correlation quality. If the user-selected model performs poorly, the system will intelligently search for a better alternative.
- **Adaptive MS2 Tolerance**: Dynamically adjusts MS2 tolerance based on the dataset characteristics, analyzing both reported and predicted tolerances to find the optimal setting.
- **Correlation Validation**: Implements a robust validation system that ensures the selected model achieves sufficient correlation with experimental spectra, preventing the use of inappropriate models.
- **Enhanced Spectrum Processing**: Uses OpenMS for spectrum file reading instead of ms2rescore_rs, providing better compatibility with a wider range of mzML files and formats.

### DeepLC Innovations

- **Model Optimization**: Automatically benchmarks pretrained vs. retrained DeepLC models for each dataset, selecting the one with the lowest Mean Absolute Error (MAE) for retention time prediction.
- **Per-Run Calibration**: Calibrates DeepLC models for each run to account for chromatographic variations between experiments, improving prediction accuracy.
- **Best Peptide Retention Time**: Tracks the best retention time prediction for each peptide across multiple PSMs, providing more reliable retention time features.
- **Transfer Learning**: Leverages transfer learning to adapt models to specific experimental conditions, improving prediction accuracy for challenging datasets.

### Spectrum Feature Analysis

Unlike traditional rescoring approaches, quantms-rescoring incorporates advanced spectrum quality metrics:

- **Signal-to-Noise Ratio (SNR)**: Calculates the ratio of maximum intensity to background noise, providing a robust measure of spectrum quality.
- **Spectral Entropy**: Quantifies the uniformity of peak distribution, helping to distinguish between high and low-quality spectra.
- **TIC Distribution Analysis**: Analyzes the distribution of Total Ion Current across peaks, identifying spectra with concentrated signal in top peaks.
- **Weighted m/z Standard Deviation**: Estimates spectral complexity by calculating the intensity-weighted standard deviation of m/z values.

### SAGE Feature Integration

- **Seamless Integration**: Incorporates additional features from SAGE (Spectrum Agnostic Generation of Embeddings) into the rescoring pipeline.
- **Feature Validation**: Ensures all features are properly validated and formatted for compatibility with OpenMS and downstream tools.

### Advantages Over Existing Tools

- **Compared to MS2PIP**: Adds automatic model selection, validation, and tolerance optimization, eliminating the need for manual parameter tuning.
- **Compared to DeepLC**: Provides automatic model selection between pretrained and retrained models, with per-run calibration for improved accuracy.
- **Compared to MS2Rescore**: Offers a more comprehensive feature set including spectrum quality metrics, better integration with OpenMS, and improved handling of different fragmentation methods and MS levels.

## Technical Implementation Details

#### Model Selection and Optimization

- **MS2PIP Model Selection**: 
  - Automatically evaluate the quality of the MS2PIP model selected by the user. If the correlation between predicted and experimental spectra is lower than a given threshold, we will try to find the best model to use (`annotator.py`). For example, if the user provides as model parameter HCD for a CI experiment, the tool will try to find the best model for this experiment within the CID models available. 
  - If the `ms_tolerance` is to restrictive for the data (e.g. 0.05 Da for a 0.5 Da dataset), the tool will try to find the annotated tolerances in the idXML file and use the best model for this tolerance.
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
- **Feature Selection**: The parameters `only_features` allows to select the features to be added to the idXML file. For example: `--only_features "DeepLC:RtDiff,DeepLC:PredictedRetentionTimeBest,Ms2pip:DotProd"`. 

##### Features

<details>
<summary>MS2PIP Feature Mapping Table</summary>

| MMS2Rescore MS2PIP Feature     | quantms-rescoring Name            |
|--------------------------------|-----------------------------------|
| spec_pearson                   | MS2PIP:SpecPearson                |
| cos_norm                       | MS2PIP:SpecCosineNorm             |
| spec_pearson_norm              | MS2PIP:SpecPearsonNorm            |
| dotprod                        | MS2PIP:DotProd                    |
| ionb_pearson_norm              | MS2PIP:IonBPearsonNorm            |
| iony_pearson_norm              | MS2PIP:IonYPearsonNorm            |
| spec_mse_norm                  | MS2PIP:SpecMseNorm                |
| ionb_mse_norm                  | MS2PIP:IonBMseNorm                |
| iony_mse_norm                  | MS2PIP:IonYMseNorm                |
| min_abs_diff_norm              | MS2PIP:MinAbsDiffNorm             |
| max_abs_diff_norm              | MS2PIP:MaxAbsDiffNorm             |
| abs_diff_Q1_norm               | MS2PIP:AbsDiffQ1Norm              |
| abs_diff_Q2_norm               | MS2PIP:AbsDiffQ2Norm              |
| abs_diff_Q3_norm               | MS2PIP:AbsDiffQ3Norm              |
| mean_abs_diff_norm             | MS2PIP:MeanAbsDiffNorm            |
| std_abs_diff_norm              | MS2PIP:StdAbsDiffNorm             |
| ionb_min_abs_diff_norm         | MS2PIP:IonBMinAbsDiffNorm         |
| ionb_max_abs_diff_norm         | MS2PIP:IonBMaxAbsDiffNorm         |
| ionb_abs_diff_Q1_norm          | MS2PIP:IonBAbsDiffQ1Norm          |
| ionb_abs_diff_Q2_norm          | MS2PIP:IonBAbsDiffQ2Norm          |
| ionb_abs_diff_Q3_norm          | MS2PIP:IonBAbsDiffQ3Norm          |
| ionb_mean_abs_diff_norm        | MS2PIP:IonBMeanAbsDiffNorm        |
| ionb_std_abs_diff_norm         | MS2PIP:IonBStdAbsDiffNorm         |
| iony_min_abs_diff_norm         | MS2PIP:IonYMinAbsDiffNorm         |
| iony_max_abs_diff_norm         | MS2PIP:IonYMaxAbsDiffNorm         |
| iony_abs_diff_Q1_norm          | MS2PIP:IonYAbsDiffQ1Norm          |
| iony_abs_diff_Q2_norm          | MS2PIP:IonYAbsDiffQ2Norm          |
| iony_abs_diff_Q3_norm          | MS2PIP:IonYAbsDiffQ3Norm          |
| iony_mean_abs_diff_norm        | MS2PIP:IonYMeanAbsDiffNorm        |
| iony_std_abs_diff_norm         | MS2PIP:IonYStdAbsDiffNorm         |
| dotprod_norm                   | MS2PIP:DotProdNorm                |
| dotprod_ionb_norm              | MS2PIP:DotProdIonBNorm            |
| dotprod_iony_norm              | MS2PIP:DotProdIonYNorm            |
| cos_ionb_norm                  | MS2PIP:CosIonBNorm                |
| cos_iony_norm                  | MS2PIP:CosIonYNorm                |
| ionb_pearson                   | MS2PIP:IonBPearson                |
| iony_pearson                   | MS2PIP:IonYPearson                |
| spec_spearman                  | MS2PIP:SpecSpearman               |
| ionb_spearman                  | MS2PIP:IonBSpearman               |
| iony_spearman                  | MS2PIP:IonYSpearman               |
| spec_mse                       | MS2PIP:SpecMse                    |
| ionb_mse                       | MS2PIP:IonBMse                    |
| iony_mse                       | MS2PIP:IonYMse                    |
| min_abs_diff_iontype           | MS2PIP:MinAbsDiffIonType          |
| max_abs_diff_iontype           | MS2PIP:MaxAbsDiffIonType          |
| min_abs_diff                   | MS2PIP:MinAbsDiff                 |
| max_abs_diff                   | MS2PIP:MaxAbsDiff                 |
| abs_diff_Q1                    | MS2PIP:AbsDiffQ1                  |
| abs_diff_Q2                    | MS2PIP:AbsDiffQ2                  |
| abs_diff_Q3                    | MS2PIP:AbsDiffQ3                  |
| mean_abs_diff                  | MS2PIP:MeanAbsDiff                |
| std_abs_diff                   | MS2PIP:StdAbsDiff                 |
| ionb_min_abs_diff              | MS2PIP:IonBMinAbsDiff             |
| ionb_max_abs_diff              | MS2PIP:IonBMaxAbsDiff             |
| ionb_abs_diff_Q1               | MS2PIP:IonBAbsDiffQ1              |
| ionb_abs_diff_Q2               | MS2PIP:IonBAbsDiffQ2              |
| ionb_abs_diff_Q3               | MS2PIP:IonBAbsDiffQ3              |
| ionb_mean_abs_diff             | MS2PIP:IonBMeanAbsDiff            |
| ionb_std_abs_diff              | MS2PIP:IonBStdAbsDiff             |
| iony_min_abs_diff              | MS2PIP:IonYMinAbsDiff             |
| iony_max_abs_diff              | MS2PIP:IonYMaxAbsDiff             |
| iony_abs_diff_Q1               | MS2PIP:IonYAbsDiffQ1              |
| iony_abs_diff_Q2               | MS2PIP:IonYAbsDiffQ2              |
| iony_abs_diff_Q3               | MS2PIP:IonYAbsDiffQ3              |
| iony_mean_abs_diff             | MS2PIP:IonYMeanAbsDiff            |
| iony_std_abs_diff              | MS2PIP:IonYStdAbsDiff             |
| dotprod_ionb                   | MS2PIP:DotProdIonB                |
| dotprod_iony                   | MS2PIP:DotProdIonY                |
| cos_ionb                       | MS2PIP:CosIonB                    |
| cos_iony                       | MS2PIP:CosIonY                    |

</details>

<details>
<summary>DeepLC Feature Mapping Table</summary>

| MMS2Rescore DeepLC Feature    | quantms-rescoring Name            |
|-------------------------------|-----------------------------------|
| observed_retention_time       | DeepLC:ObservedRetentionTime      |
| predicted_retention_time      | DeepLC:PredictedRetentionTime     |
| rt_diff                       | DeepLC:RtDiff                     |
| observed_retention_time_best  | DeepLC:ObservedRetentionTimeBest  |
| predicted_retention_time_best | DeepLC:PredictedRetentionTimeBest |
| rt_diff_best                  | DeepLC:RtDiffBest                 |

</details>

<details>
<summary>Spectrum Feature Mapping Table</summary>

| Spectrum Feature    | quantms-rescoring Name            |
|---------------------|-----------------------------------|
| snr                 | Quantms:Snr                       |
| spectral_entropy    | Quantms:SpectralEntropy           |
| fraction_tic_top_10 | Quantms:FracTICinTop10Peaks       |
| weighted_std_mz     | Quantms:WeightedStdMz             |

</details>

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
