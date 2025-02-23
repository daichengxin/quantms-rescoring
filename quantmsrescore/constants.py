MS2PIP_FEATURES = [
    {"MS2PIP:SpecPearsonNorm": "spec_pearson"},
    {"MS2PIP:SpecCosineNorm": "cos_norm"},
    {"MS2PIP:IonbPersonNorm": "ionb_pearson_norm"},
    {"MS2PIP:IonyPersonNorm": "iony_pearson_norm"},
]

DEEPLC_FEATURES = [
    {"DeepLC:RtPredicted": "predicted_retention_time"},
    {"DeepLC:RtPredictedBest": "predicted_retention_time_best"},
    {"DeepLC:RtDiff": "rt_diff"},  # Different retention time between observed and predicted
    {"DeepLC:RtDiffBest": "rt_diff_best"},  # Best retention time difference
    {"DeepLC:RtPredictedBest": "predicted_retention_time_best"},
]

SUPPORTED_MODELS_MS2PIP = [
    "HCD2019",  # HCD from 2019
    "HCD2021",  # Default model
    "HCD",  # Default model
    "CID",  # Collision-induced dissociation
    "iTRAQ",  # iTRAQ
    "iTRAQphospho",  # iTRAQ phospho
    "TMT",  # TMT 10-plex
    "TTOF5600",  # TTOF 5600
    "HCDch2",  # HCD with charge 2
    "CIDch2",  # CID with charge 2
    "Immuno-HCD",  # Immuno-HCD
    "CID-TMT",  # CID-TMT
    "timsTOF2023",  # TIMS-TOF 2023
    "timsTOF2024",  # TIMS-TOF 2024
]

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
