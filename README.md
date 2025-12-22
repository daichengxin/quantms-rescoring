# quantms-rescoring
    
[![Python package](https://github.com/bigbio/quantms-rescoring/actions/workflows/python-package.yml/badge.svg)](https://github.com/bigbio/quantms-rescoring/actions/workflows/python-package.yml)
[![codecov](https://codecov.io/gh/bigbio/quantms-rescoring/branch/main/graph/badge.svg?token=3ZQZQ2ZQ2D)](https://codecov.io/gh/bigbio/quantms-rescoring)
[![PyPI version](https://badge.fury.io/py/quantms-rescoring.svg)](https://badge.fury.io/py/quantms-rescoring)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)

quantms-rescoring is a Python tool that aims to add features to peptide-spectrum matches (PSMs) in idXML files using multiple tools including SAGE features, quantms spectrum features, MS2PIP, AlphaPeptDeep and DeepLC. It is part of the quantms ecosystem package and leverages the MS²Rescore framework to improve identification confidence in proteomics data analysis.

## Core Components

- **Annotator Engine**: Integrates [MS2PIP](https://github.com/compomics/ms2pip), [AlphaPeptDeep](https://github.com/MannLabs/alphapeptdeep) and [DeepLC](https://github.com/compomics/DeepLC) models to improve peptide-spectrum match (PSM) confidence. 
- **Feature Generation**: Extracts signal-to-noise ratios, spectrum metrics, SAGE extra features and add them to each PSM for posterior downstream with Percolator.
- **OpenMS Integration**: Processes idXML and mzML files with custom validation methods.
- **Transfer learning**: Implemented train and fine-tune the ms2 model to generate project-specific model and pass the new model to quantms workflow for rescoring.

## CLI Tools

```sh
 rescoring msrescore2feature --help
```
Annotates PSMs with prediction-based features from MS2PIP and DeepLC

```sh
 rescoring add_sage_feature --help
```
Incorporates additional features from SAGE into idXML files.

```sh
 rescoring spectrum2feature --help
```
Add additional spectrum feature like signal-to-noise to each PSM in the idXML.

```sh
 rescoring psm_feature_clean --help
```
Check and clean invalid PSM with invalid features in the idXML.

```sh
 rescoring download_models --help
```
Download all required models (MS2PIP, AlphaPeptDeep) for offline use. This is useful for running quantms-rescoring in environments without internet access, such as HPC clusters.

```sh
 rescoring transfer_learning --help
```
Perform train and fine-tuning model for AlphaPeptDeep to generate project-specific model. This is useful for running quantms-rescoring in PTM datasets, which are unseen for AlphaPeptDeep pretrained model.

## Advanced Algorithms and Improvements

quantms-rescoring significantly enhances the capabilities of MS2PIP, AlphaPeptDeep, DeepLC, and MS2Rescore through several innovative approaches:

### MS2PIP Integration Enhancements

- **Intelligent Model Selection**: Automatically evaluates and selects the optimal MS2PIP model for each dataset based on fragmentation type and correlation quality. If the user-selected model performs poorly, the system will intelligently search for a better alternative.
- **Adaptive MS2 Tolerance**: Dynamically adjusts MS2 tolerance based on the dataset characteristics, analyzing both reported and predicted tolerances to find the optimal setting.
- **Correlation Validation**: Implements a robust validation system that ensures the selected model achieves sufficient correlation with experimental spectra, preventing the use of inappropriate models.
- **Enhanced Spectrum Processing**: Uses OpenMS for spectrum file reading instead of ms2rescore_rs, providing better compatibility with a wider range of mzML files and formats.

### AlphaPeptDeep Innovations

- **Fine-tuning**: Leverages fine-tuning to adapt models to specific experimental/project conditions based on identifications file (idXML) from quantms, improving prediction accuracy for challenging datasets, such as PTM.
- **Model Optimization**: Automatically benchmarks pretrained vs. retrained AlphaPeptDeep models for each dataset, selecting the one with the better Median PCC for MS2 intensity prediction.
- **Enhanced Spectrum Processing**: AlphaPeptDeep does not support idXML input, so we use OpenMS for spectrum file reading and pass it to AlphaPeptDeep for prediction and fine-tuning.
- **Correlation Validation**: Implements a robust validation system that ensures the pretrained and retrained models achieve sufficient correlation with experimental spectra, preventing the use of inappropriate models.

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

- **Compared to MS2PIP**: Adds automatic model selection, validation, features calculations and tolerance optimization, eliminating the need for manual parameter tuning.
- **Compared to DeepLC**: Provides automatic model selection between pretrained and retrained models, with per-run calibration for improved accuracy.
- **Compared to MS2Rescore**: Integrates a broader range of MS2 prediction models, including AlphaPeptDeep, and supports fine-tuning to generate project-specific models. It provides a richer feature set encompassing spectrum quality metrics, tighter integration with OpenMS, and more robust support for diverse fragmentation methods and MS levels. 
- **Compared to AlphaPeptDeep**: Seamlessly integrates into the quantms workflow and natively supports the quantms identification results format. It introduces automatic model selection and validation, delivers an expanded feature set, and offers improved handling of different fragmentation methods and MS levels.


## Technical Implementation Details

#### Model Selection and Optimization

- **MS2PIP Model Selection**: 
  - Automatically evaluate the quality of the MS2PIP model selected by the user. If the correlation between predicted and experimental spectra is lower than a given threshold, we will try to find the best model to use (`annotator.py`). For example, if the user provides as model parameter HCD for a CI experiment, the tool will try to find the best model for this experiment within the CID models available. 
  - If the `ms_tolerance` is to restrictive for the data (e.g. 0.05 Da for a 0.5 Da dataset), the tool will try to find the annotated tolerances in the idXML file and use the best model for this tolerance.
- **AlphaPeptDeep Model**:
  - Automatically evaluate the quality of the AlphaPeptDeep model weight passed by the user. If the correlation between predicted and experimental spectra is lower than a given threshold, we will skip MS2 features generation to avoid potential erroneous results.
  - When enabling `transfer_learning`, the tool will try to fine-tune the AlphaPeptDeep model on the given idXML and mzML files and compare it with the pretrained model, finally using the best model based on similarity metrics.
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

#### Features

<details>
<summary>MS2PIP and AlphaPeptDeep Feature Mapping Table</summary>

| MS2PIP and AlphaPeptDeep Feature | quantms-rescoring Name                   |
|----------------------------------|------------------------------------------|
| spec_pearson                     | MS2PIP/AlphaPeptDeep:SpecPearson         |
| cos_norm                         | MS2PIP/AlphaPeptDeep:SpecCosineNorm      |
| spec_pearson_norm                | MS2PIP/AlphaPeptDeep:SpecPearsonNorm     |
| dotprod                          | MS2PIP/AlphaPeptDeep:DotProd             |
| ionb_pearson_norm                | MS2PIP/AlphaPeptDeep:IonBPearsonNorm     |
| iony_pearson_norm                | MS2PIP/AlphaPeptDeep:IonYPearsonNorm     |
| spec_mse_norm                    | MS2PIP/AlphaPeptDeep:SpecMseNorm         |
| ionb_mse_norm                    | MS2PIP/AlphaPeptDeep:IonBMseNorm         |
| iony_mse_norm                    | MS2PIP/AlphaPeptDeep:IonYMseNorm         |
| min_abs_diff_norm                | MS2PIP/AlphaPeptDeep:MinAbsDiffNorm      |
| max_abs_diff_norm                | MS2PIP/AlphaPeptDeep:MaxAbsDiffNorm      |
| abs_diff_Q1_norm                 | MS2PIP/AlphaPeptDeep:AbsDiffQ1Norm       |
| abs_diff_Q2_norm                 | MS2PIP/AlphaPeptDeep:AbsDiffQ2Norm       |
| abs_diff_Q3_norm                 | MS2PIP/AlphaPeptDeep:AbsDiffQ3Norm       |
| mean_abs_diff_norm               | MS2PIP/AlphaPeptDeep:MeanAbsDiffNorm     |
| std_abs_diff_norm                | MS2PIP/AlphaPeptDeep:StdAbsDiffNorm      |
| ionb_min_abs_diff_norm           | MS2PIP/AlphaPeptDeep:IonBMinAbsDiffNorm  |
| ionb_max_abs_diff_norm           | MS2PIP/AlphaPeptDeep:IonBMaxAbsDiffNorm  |
| ionb_abs_diff_Q1_norm            | MS2PIP/AlphaPeptDeep:IonBAbsDiffQ1Norm   |
| ionb_abs_diff_Q2_norm            | MS2PIP/AlphaPeptDeep:IonBAbsDiffQ2Norm   |
| ionb_abs_diff_Q3_norm            | MS2PIP/AlphaPeptDeep:IonBAbsDiffQ3Norm   |
| ionb_mean_abs_diff_norm          | MS2PIP/AlphaPeptDeep:IonBMeanAbsDiffNorm |
| ionb_std_abs_diff_norm           | MS2PIP/AlphaPeptDeep:IonBStdAbsDiffNorm  |
| iony_min_abs_diff_norm           | MS2PIP/AlphaPeptDeep:IonYMinAbsDiffNorm  |
| iony_max_abs_diff_norm           | MS2PIP/AlphaPeptDeep:IonYMaxAbsDiffNorm  |
| iony_abs_diff_Q1_norm            | MS2PIP/AlphaPeptDeep:IonYAbsDiffQ1Norm   |
| iony_abs_diff_Q2_norm            | MS2PIP/AlphaPeptDeep:IonYAbsDiffQ2Norm   |
| iony_abs_diff_Q3_norm            | MS2PIP/AlphaPeptDeep:IonYAbsDiffQ3Norm   |
| iony_mean_abs_diff_norm          | MS2PIP/AlphaPeptDeep:IonYMeanAbsDiffNorm |
| iony_std_abs_diff_norm           | MS2PIP/AlphaPeptDeep:IonYStdAbsDiffNorm  |
| dotprod_norm                     | MS2PIP/AlphaPeptDeep:DotProdNorm         |
| dotprod_ionb_norm                | MS2PIP/AlphaPeptDeep:DotProdIonBNorm     |
| dotprod_iony_norm                | MS2PIP/AlphaPeptDeep:DotProdIonYNorm     |
| cos_ionb_norm                    | MS2PIP/AlphaPeptDeep:CosIonBNorm         |
| cos_iony_norm                    | MS2PIP/AlphaPeptDeep:CosIonYNorm         |
| ionb_pearson                     | MS2PIP/AlphaPeptDeep:IonBPearson         |
| iony_pearson                     | MS2PIP/AlphaPeptDeep:IonYPearson         |
| spec_spearman                    | MS2PIP/AlphaPeptDeep:SpecSpearman        |
| ionb_spearman                    | MS2PIP/AlphaPeptDeep:IonBSpearman        |
| iony_spearman                    | MS2PIP/AlphaPeptDeep:IonYSpearman        |
| spec_mse                         | MS2PIP/AlphaPeptDeep:SpecMse             |
| ionb_mse                         | MS2PIP/AlphaPeptDeep:IonBMse             |
| iony_mse                         | MS2PIP/AlphaPeptDeep:IonYMse             |
| min_abs_diff_iontype             | MS2PIP/AlphaPeptDeep:MinAbsDiffIonType   |
| max_abs_diff_iontype             | MS2PIP/AlphaPeptDeep:MaxAbsDiffIonType   |
| min_abs_diff                     | MS2PIP/AlphaPeptDeep:MinAbsDiff          |
| max_abs_diff                     | MS2PIP/AlphaPeptDeep:MaxAbsDiff          |
| abs_diff_Q1                      | MS2PIP/AlphaPeptDeep:AbsDiffQ1           |
| abs_diff_Q2                      | MS2PIP/AlphaPeptDeep:AbsDiffQ2           |
| abs_diff_Q3                      | MS2PIP/AlphaPeptDeep:AbsDiffQ3           |
| mean_abs_diff                    | MS2PIP/AlphaPeptDeep:MeanAbsDiff         |
| std_abs_diff                     | MS2PIP/AlphaPeptDeep:StdAbsDiff          |
| ionb_min_abs_diff                | MS2PIP/AlphaPeptDeep:IonBMinAbsDiff      |
| ionb_max_abs_diff                | MS2PIP/AlphaPeptDeep:IonBMaxAbsDiff      |
| ionb_abs_diff_Q1                 | MS2PIP/AlphaPeptDeep:IonBAbsDiffQ1       |
| ionb_abs_diff_Q2                 | MS2PIP/AlphaPeptDeep:IonBAbsDiffQ2       |
| ionb_abs_diff_Q3                 | MS2PIP/AlphaPeptDeep:IonBAbsDiffQ3       |
| ionb_mean_abs_diff               | MS2PIP/AlphaPeptDeep:IonBMeanAbsDiff     |
| ionb_std_abs_diff                | MS2PIP/AlphaPeptDeep:IonBStdAbsDiff      |
| iony_min_abs_diff                | MS2PIP/AlphaPeptDeep:IonYMinAbsDiff      |
| iony_max_abs_diff                | MS2PIP/AlphaPeptDeep:IonYMaxAbsDiff      |
| iony_abs_diff_Q1                 | MS2PIP/AlphaPeptDeep:IonYAbsDiffQ1       |
| iony_abs_diff_Q2                 | MS2PIP/AlphaPeptDeep:IonYAbsDiffQ2       |
| iony_abs_diff_Q3                 | MS2PIP/AlphaPeptDeep:IonYAbsDiffQ3       |
| iony_mean_abs_diff               | MS2PIP/AlphaPeptDeep:IonYMeanAbsDiff     |
| iony_std_abs_diff                | MS2PIP/AlphaPeptDeep:IonYStdAbsDiff      |
| dotprod_ionb                     | MS2PIP/AlphaPeptDeep:DotProdIonB         |
| dotprod_iony                     | MS2PIP/AlphaPeptDeep:DotProdIonY         |
| cos_ionb                         | MS2PIP/AlphaPeptDeep:CosIonB             |
| cos_iony                         | MS2PIP/AlphaPeptDeep:CosIonY             |

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

**Using Docker**

```sh
# Pull the latest image from GitHub Container Registry
❯ docker pull ghcr.io/bigbio/quantms-rescoring:latest

# Run the container
❯ docker run --rm ghcr.io/bigbio/quantms-rescoring:latest --help

# Run with data mounted
❯ docker run --rm -v /path/to/data:/data ghcr.io/bigbio/quantms-rescoring:latest rescoring msrescore2feature --help
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

**Build Docker image locally:**

```sh
# Build the Docker image
❯ docker build -t quantms-rescoring:latest .

# Test the image
❯ ./scripts/test-docker.sh latest
```

### Offline Model Download

For environments without internet access (e.g., HPC clusters), you can download all required models ahead of time using the `download_models` command:

```sh
# Download all models to default cache locations
❯ rescoring download_models

# Download models to a specific directory
❯ rescoring download_models --model_dir /path/to/models

# Download only specific models
❯ rescoring download_models --models deeplc,alphapeptdeep

# Get help
❯ rescoring download_models --help
```

This command downloads models for:

- **MS2PIP**: Fragment ion intensity prediction models (bundled with ms2pip package)
- **AlphaPeptDeep**: MS2 spectrum, retention time, and CCS prediction models

Once downloaded, you can transfer the models to your offline environment and use them with the processing commands. For AlphaPeptDeep models, use the `--ms2_model_dir` option when running `msrescore2feature`.

### TODO

- [ ] Add support for multiple Files combined idXML and mzML

### Issues and Contributions

For any issues or contributions, please open an issue in the [GitHub repository](https://github.com/bigbio/quantms/issues) - we use the quantms repo to control all issues—or PR in the [GitHub repository](https://github.com/bigbio/quantms-rescoring/pulls).
