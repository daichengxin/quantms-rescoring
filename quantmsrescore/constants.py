MS2PIP_FEATURES = {
    "Ms2pip:SpecPearson": "spec_pearson",
    "Ms2pip:SpecCosineNorm": "cos_norm",
    "Ms2pip:SpecPearsonNorm": "spec_pearson_norm",
    "Ms2pip:DotProd": "dotprod",
    "Ms2pip:IonBPearsonNorm": "ionb_pearson_norm",
    "Ms2pip:IonYPearsonNorm": "iony_pearson_norm",
    "Ms2pip:SpecMseNorm": "spec_mse_norm",
    "Ms2pip:IonBMseNorm": "ionb_mse_norm",
    "Ms2pip:IonYMseNorm": "iony_mse_norm",
    "Ms2pip:MinAbsDiffNorm": "min_abs_diff_norm",
    "Ms2pip:MaxAbsDiffNorm": "max_abs_diff_norm",
    "Ms2pip:AbsDiffQ1Norm": "abs_diff_Q1_norm",
    "Ms2pip:AbsDiffQ2Norm": "abs_diff_Q2_norm",
    "Ms2pip:AbsDiffQ3Norm": "abs_diff_Q3_norm",
    "Ms2pip:MeanAbsDiffNorm": "mean_abs_diff_norm",
    "Ms2pip:StdAbsDiffNorm": "std_abs_diff_norm",
    "Ms2pip:IonBMinAbsDiffNorm": "ionb_min_abs_diff_norm",
    "Ms2pip:IonBMaxAbsDiffNorm": "ionb_max_abs_diff_norm",
    "Ms2pip:IonBAbsDiffQ1Norm": "ionb_abs_diff_Q1_norm",
    "Ms2pip:IonBAbsDiffQ2Norm": "ionb_abs_diff_Q2_norm",
    "Ms2pip:IonBAbsDiffQ3Norm": "ionb_abs_diff_Q3_norm",
    "Ms2pip:IonBMeanAbsDiffNorm": "ionb_mean_abs_diff_norm",
    "Ms2pip:IonBStdAbsDiffNorm": "ionb_std_abs_diff_norm",
    "Ms2pip:IonYMinAbsDiffNorm": "iony_min_abs_diff_norm",
    "Ms2pip:IonYMaxAbsDiffNorm": "iony_max_abs_diff_norm",
    "Ms2pip:IonYAbsDiffQ1Norm": "iony_abs_diff_Q1_norm",
    "Ms2pip:IonYAbsDiffQ2Norm": "iony_abs_diff_Q2_norm",
    "Ms2pip:IonYAbsDiffQ3Norm": "iony_abs_diff_Q3_norm",
    "Ms2pip:IonYMeanAbsDiffNorm": "iony_mean_abs_diff_norm",
    "Ms2pip:IonYStdAbsDiffNorm": "iony_std_abs_diff_norm",
    "Ms2pip:DotProdNorm": "dotprod_norm",
    "Ms2pip:DotProdIonBNorm": "dotprod_ionb_norm",
    "Ms2pip:DotProdIonYNorm": "dotprod_iony_norm",
    "Ms2pip:CosIonBNorm": "cos_ionb_norm",
    "Ms2pip:CosIonYNorm": "cos_iony_norm",
    "Ms2pip:IonBPearson": "ionb_pearson",
    "Ms2pip:IonYPearson": "iony_pearson",
    "Ms2pip:SpecSpearman": "spec_spearman",
    "Ms2pip:IonBSpearman": "ionb_spearman",
    "Ms2pip:IonYSpearman": "iony_spearman",
    "Ms2pip:SpecMse": "spec_mse",
    "Ms2pip:IonBMse": "ionb_mse",
    "Ms2pip:IonYMse": "iony_mse",
    "Ms2pip:MinAbsDiffIonType": "min_abs_diff_iontype",
    "Ms2pip:MaxAbsDiffIonType": "max_abs_diff_iontype",
    "Ms2pip:MinAbsDiff": "min_abs_diff",
    "Ms2pip:MaxAbsDiff": "max_abs_diff",
    "Ms2pip:AbsDiffQ1": "abs_diff_Q1",
    "Ms2pip:AbsDiffQ2": "abs_diff_Q2",
    "Ms2pip:AbsDiffQ3": "abs_diff_Q3",
    "Ms2pip:MeanAbsDiff": "mean_abs_diff",
    "Ms2pip:StdAbsDiff": "std_abs_diff",
    "Ms2pip:IonBMinAbsDiff": "ionb_min_abs_diff",
    "Ms2pip:IonBMaxAbsDiff": "ionb_max_abs_diff",
    "Ms2pip:IonBAbsDiffQ1": "ionb_abs_diff_Q1",
    "Ms2pip:IonBAbsDiffQ2": "ionb_abs_diff_Q2",
    "Ms2pip:IonBAbsDiffQ3": "ionb_abs_diff_Q3",
    "Ms2pip:IonBMeanAbsDiff": "ionb_mean_abs_diff",
    "Ms2pip:IonBStdAbsDiff": "ionb_std_abs_diff",
    "Ms2pip:IonYMinAbsDiff": "iony_min_abs_diff",
    "Ms2pip:IonYMaxAbsDiff": "iony_max_abs_diff",
    "Ms2pip:IonYAbsDiffQ1": "iony_abs_diff_Q1",
    "Ms2pip:IonYAbsDiffQ2": "iony_abs_diff_Q2",
    "Ms2pip:IonYAbsDiffQ3": "iony_abs_diff_Q3",
    "Ms2pip:IonYMeanAbsDiff": "iony_mean_abs_diff",
    "Ms2pip:IonYStdAbsDiff": "iony_std_abs_diff",
    "Ms2pip:DotProdIonB": "dotprod_ionb",
    "Ms2pip:DotProdIonY": "dotprod_iony",
    "Ms2pip:CosIonB": "cos_ionb",
    "Ms2pip:CosIonY": "cos_iony",
}

DEEPLC_FEATURES = {
    "DeepLC:ObservedRetentionTime": "observed_retention_time",
    "DeepLC:PredictedRetentionTime": "predicted_retention_time",
    "DeepLC:RtDiff": "rt_diff",
    "DeepLC:ObservedRetentionTimeBest": "observed_retention_time_best",
    "DeepLC:PredictedRetentionTimeBest": "predicted_retention_time_best",
    "DeepLC:RtDiffBest": "rt_diff_best",
}

QUANTMS_FEATURES = {
    "Quantms:Snr": "snr",
    "Quantms:SpectralEntropy": "spectral_entropy",
    "Quantms:FracTICinTop10Peaks": "fraction_tic_top_10",
    "Quantms:WeightedStdMz": "weighted_std_mz",
}

SUPPORTED_MODELS_MS2PIP = {
    "HCD": [
        "HCD2019",  # HCD from 2019
        "HCD2021",  # Default model
        "Immuno-HCD",  # Immuno-HCD
        "HCDch2",  # HCD with charge 2
        "TMT",  # TMT 10-plex
        "iTRAQ",  # iTRAQ
        "iTRAQphospho",  # iTRAQ phospho
    ],
    "CID": [
        "CID",  # Collision-induced dissociation
        "CIDch2",  # CID with charge 2
        "CID-TMT",  # CID-TMT
    ],
}

# This is the list of disassociation methods that are supported by OPENMS.
# This list is a path for release 3.3.0 of OpenMS.
OPENMS_DISSOCIATION_METHODS_PATCH = [
    {
        "CID": "Collision-induced dissociation (MS:1000133) (also CAD; parent term, but unless otherwise stated often used as synonym for trap-type CID)"
    },
    {"PSD": "Post-source decay."},
    {"PD": "Plasma desorption."},
    {"SID": "Surface-induced dissociation."},
    {"BIRD": "Blackbody infrared radiative dissociation."},
    {"ECD": "Electron capture dissociation (MS:1000250)"},
    {"IMD": "Infrared multiphoton dissociation."},
    {"SORI": "Sustained off-resonance irradiation."},
    {"HCID": "High-energy collision-induced dissociation."},
    {"LCID": "Low-energy collision-induced dissociation."},
    {"PHD": "Photodissociation."},
    {"ETD": "Electron transfer dissociation."},
    {"ETciD": "Electron transfer and collision-induced dissociation (MS:1003182)"},
    {"EThcD": "Electron transfer and higher-energy collision dissociation (MS:1002631)"},
    {"PQD": "Pulsed q dissociation (MS:1000599)"},
    {"TRAP": "trap-type collision-induced dissociation (MS:1002472)"},
    {"HCD": "beam-type collision-induced dissociation (MS:1000422)"},
    {"INSOURCE": "in-source collision-induced dissociation (MS:1001880)"},
    {"LIFT": "Bruker proprietary method (MS:1002000)"},
]
