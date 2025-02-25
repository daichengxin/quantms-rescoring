MS2PIP_FEATURES = [
    {"Ms2pip:SpecPearsonNorm": "spec_pearson"},
    {"Ms2pip:SpecCosineNorm": "cos_norm"},
]

DEEPLC_FEATURES = [
    {"DeepLC:RtPredicted": "predicted_retention_time"},
    {"DeepLC:RtPredictedBest": "predicted_retention_time_best"},
    {"DeepLC:RtDiff": "rt_diff"},  # Different retention time between observed and predicted
    {"DeepLC:RtDiffBest": "rt_diff_best"},  # Best retention time difference
]

QUANTMS_FEATURES = [
    {"Quantms:Snr": "snr"},
    {"Quantms:SpectralEntropy": "spectral_entropy"},
    {"Quantms:FracTICinTop10Peaks": "fraction_tic_top_10"},
    {"Quantms:WeightedStdMz": "weighted_std_mz"},
]

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
