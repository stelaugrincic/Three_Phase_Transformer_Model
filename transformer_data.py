#---------------------------------------------------------------------------------------------------------------
# Author: Stela Ugrincic
# Affiliation: University of Maribor, Faculty of Electrical Engineering and Computer Science, Maribor, Slovenia
# Email: stela.ugrincic1@um.si
# Last modified: 16.03.2026
#---------------------------------------------------------------------------------------------------------------

# TRANSFORMER CONFIGURATION DATA
# Description: This file contains rated data, construction parameters, 
# and magnetic characteristic of a three-phase three-limb transformer

# Three-phase three-limb transformer:
# 750 VA, Upn = 400 V, Usn = 300 V, Ipn = 1,2 A, Isn = 1,49 A
# Np = 464, Ns = 355

import numpy as np

podatki_transformatorja = {
    # Rated data
    "Sn": 750,        #(VA) - Rated power
    "Upn_rms": 400,   #(V) - Primary rated voltage RMS
    "Usn_rms": 300,   #(V) - Secondary rated voltage RMS
    "Ipn": 1.2,       #(A) - Primary rated current
    "Isn": 1.49,      #(A) - Secondary rated current
    "f": 50.0,        #(Hz) - Rated frequency

    "PRIMAR": 'Y',    # Primary connection (Y)
    "SEKUNDAR": 'y',  # Secondary connection (y)
    "VEZALNA_SKUPINA": 'Yy0', # Vector group

    # Construction data
    "Np": 464,  # Number of primary turns
    "Ns": 355,  # Number of secondary turns
    "Am": np.array([0.001498, 0.001498, 0.001498]), # Core cross-section area (m^2)
    "lm": np.array([0.21, 0.120, 0.21]), # Magnetic path length (m)

    # Electrical parameters
    # Ohmic resistance
    "Rp": np.diag([4.852, 4.852, 4.852]),
    "Rs": np.diag([3.387, 3.387, 3.387]),
    # Leakage inductances
    "Lsp": np.diag([0.001127, 0.001127, 0.001127]),
    "Lss": np.diag([0.000742, 0.000742, 0.000742]),

    "Rg": np.diag([0, 0, 0]),      # Ω   
    "Lg": np.diag([3e-3, 3e-3, 3e-3]),   # H 

    # Magnetic Characteristic (B-H curve data)
    "B_orig_pr": np.array([0, 0.975, 1.000, 1.025, 1.050, 1.075, 1.100, 1.125, 1.150, 1.175, 1.201, 1.226, 1.251,
                           1.276, 1.300, 1.326, 1.351, 1.376, 1.401, 1.426, 1.451, 1.476, 1.502, 1.527, 1.552,
                           1.578, 1.604, 1.630, 1.655, 1.682, 1.709, 1.735, 1.762, 1.788, 1.816, 1.843]),  
    "H_orig_pr": np.array([0, 156.5, 162.8, 168.7, 175.4, 182.4, 190.0, 198.1, 206.9, 216.5, 226.8, 238.2, 250.8,
                           264.8, 280.6, 298.7, 319.6, 344.3, 373.9, 410.3, 456.4, 514.1, 589.6, 687.5, 814.5,
                           977.8, 1181, 1432, 1730, 2087, 2494, 2954, 3472, 4013, 4676, 5368]),                                            

    
    # Initial conditions - remanent flux density

    "B0": np.array([0*1.5, 0*1.5, 0*1.5]), 

    # Load data:
    # "Rb": np.diag([116.0, 116.0, 116.0]),
    # "Lb": np.diag([0, 0, 0]),
    # No load:
    "Rb": np.diag([0, 0, 0]),
    "Lb": np.diag([0, 0, 0])
}
