MS2PIP_FEATURES = [
    {"MS2PIP:SpecPearsonNorm": "spec_pearson"},
    {"MS2PIP:SpecCosineNorm": "cos_norm"},
    {"MS2PIP:IonbPersonNorm", "ionb_pearson_norm"},
    {"MS2PIP:IonyPersonNorm", "iony_pearson_norm"},
]

DEEPLC_FEATURES = [
    {"DeepLC:RtPredicted": "predicted_retention_time"},
    {"DeepLC:RtPredictedBest": "predicted_retention_time_best"},
    {"DeepLC:RtDiff": "rt_diff"},  # Different retention time between observed and predicted
    {"DeepLC:RtDiffBest": "rt_diff_best"},  # Best retention time difference
    {"DeepLC:RtPredictedBest": "predicted_retention_time_best"},
]
