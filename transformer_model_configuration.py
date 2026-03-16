#---------------------------------------------------------------------------------------------------------------
# Author: Stela Ugrincic
# Affiliation: University of Maribor, Faculty of Electrical Engineering and Computer Science, Maribor, Slovenia
# Email: stela.ugrincic1@um.si
# Last modified: 16.03.2026
#---------------------------------------------------------------------------------------------------------------

# TRANSFORMER MODEL CONFIGURATION
# Description: Setup of simulation scenarios, transformation matrices, 
# and B-H curve modeling using pchip interpolation
import numpy as np
from scipy.interpolate import PchipInterpolator
from transformer_data import podatki_transformatorja as trafo_data

# SIMULATION SETTINGS
# Simulation scenarios: 'vklopni_pojav' (Inrush), 'stacionarno_stanje' (Steady-state)
                        #'kratek_stik' (Short-circuit), 'RL' (RL load), 'prosti_tek' (No load)

# Slovenian: Definiranje scenarija: 'prosti_tek', 'stacionarno_stanje', 'kratek_stik', 'RL'
NACIN_SIMULACIJE = 'vklopni_pojav'

if NACIN_SIMULACIJE == 'vklopni_pojav':
    tip_bremena = 'prosti_tek'
elif NACIN_SIMULACIJE == 'stacionarno_stanje':
    tip_bremena = 'RL'
elif NACIN_SIMULACIJE == 'kratek_stik':
    tip_bremena = 'kratek_stik'    

PRIMARNA_VEZAVA = trafo_data['PRIMAR']
SEKUNDARNA_VEZAVA = trafo_data['SEKUNDAR']

# PRIMARY SIDE
# SLOVENIAN: NASTAVITVE ZA PRIMARNO STRAN
if PRIMARNA_VEZAVA == 'D':
    m = 2 # Number of primary state variables
    Tp = (1/3) * np.array([[1, -1], [1, 2], [-2, -1]])
    T_omrezja = np.array([[1, 0], [0, 1], [-1, -1]])
elif PRIMARNA_VEZAVA == 'Y':
    m = 3
    Tp = np.identity(3)
    T_omrezja = np.identity(3)
elif PRIMARNA_VEZAVA == 'YN':
    m = 3
    Tp = np.identity(3)
    T_omrezja = np.identity(3)

# SECONDARY SIDE
# SLOVENIAN: NASTAVITVE ZA SEKUNDARNO STRAN 
if SEKUNDARNA_VEZAVA == 'd':
    n = 2 # Number of secondary state variables
    Ts = np.array([[1, 1], [1, 2], [0, 1]])
    Tb = np.array([[1, 0], [0, 1], [-1, -1]])
elif SEKUNDARNA_VEZAVA == 'y':
    n = 3
    Ts = np.identity(3)
    Tb = np.identity(3) 
elif SEKUNDARNA_VEZAVA == 'yn':
    n = 3
    Ts = np.identity(3)
    Tb = np.identity(3)        

# Incidence matrix Cm for three-limb transformer core
Cm = np.array([[1, -1, 0],
               [0, -1, 1]])

Np = trafo_data['Np']
Ns = trafo_data['Ns']

# Magnetic flux and Magnetomotive force matrices
Npsi = np.array([
    [Np, 0, 0],
    [0, Np, 0],
    [0, 0, Np],
    [Ns, 0, 0],
    [0, Ns, 0],
    [0, 0, Ns]
]) 

Ntps = np.array([[Np, 0, 0, -Ns, 0, 0],
                 [0, Np, 0, 0, -Ns, 0],
                 [0, 0, Np, 0, 0, -Ns]])

# Voltage calculations
omega = 2 * np.pi * trafo_data['f']
# Izracun U_peak glede na simulirano vezavo
if PRIMARNA_VEZAVA == 'D':
    U_peak = trafo_data['Upn_rms'] * np.sqrt(2)
elif PRIMARNA_VEZAVA in ['Y', 'YN']:
    U_peak = (trafo_data['Upn_rms'] / np.sqrt(3)) * np.sqrt(2)     

# Iskanje nicle ali temenske vrednosti
phi = np.array([0.0, -2*np.pi/3, +2*np.pi/3])

def iskanje_tv(t_ref, phi_k, tocka="peak", smer="down"):
    """
    ENGLISH
    Calculates the time of the next voltage event (zero crossing or peak)
    after the reference time t_ref.
    """
    """
    SLOVENIAN
    Iskanje casa naslednjega dogodka za fazo s faznim kotom phi_k po t_ref
    tocka -> "peak"
    tocka -> "zero"; nicelni prehod:
        smer -> "up" ... narascajoci (f'>0)
             -> "down" ... padajoci (f'<0)
             -> "any" ... katerikoli
    """
    w = omega
    if tocka == "peak":
        # ωt + φ = π/2 + 2πk
        k = np.ceil((w*t_ref + phi_k - np.pi/2) / (2*np.pi))
        return (np.pi/2 - phi_k + 2*np.pi*k) / w

    elif tocka == "zero":
        if smer == "up":
            # ωt + φ = 2πm  (naraščajoči prehod)
            m = np.ceil((w*t_ref + phi_k) / (2*np.pi))
            return (2*np.pi*m - phi_k) / w
        elif smer == "down":
            # ωt + φ = π + 2πm  (padajoči prehod)
            m = np.ceil((w*t_ref + phi_k - np.pi) / (2*np.pi))
            return (np.pi + 2*np.pi*m - phi_k) / w
        elif smer == "any":
            # ωt + φ = πk  (katerikoli ničelni prehod)
            k = np.ceil((w*t_ref + phi_k) / np.pi)
            return (np.pi*k - phi_k) / w
        else:
            raise ValueError("zero must be 'up', 'down', or 'any'")
    else:
        raise ValueError("kind must be 'peak' or 'zero'")



# Switching conditions | Controlled switching
# Slovenian: Definiranje trenutkov vklopa
t0 = 0.095 # base switching time

t_vklop_a = float(iskanje_tv(t0, phi[0], tocka="zero", smer="down")) 
t_vklop_b = t_vklop_a
t_vklop_c = t_vklop_a


def Uomrezja_func(t):
    return np.array([
        U_peak * np.sin(omega * t),
        U_peak * np.sin(omega * t - 2 * np.pi / 3),
        U_peak * np.sin(omega * t + 2 * np.pi / 3)
    ])

# Defining the grid voltage with switching logic
def Uomrezja_vklop(t):
    # Kontrolirani vklop
    Ua = U_peak * np.sin(omega * t) if t >= t_vklop_a else 0
    Ub = U_peak * np.sin(omega * t - 2 * np.pi / 3) if t >= t_vklop_b else 0
    Uc = U_peak * np.sin(omega * t + 2 * np.pi / 3) if t >= t_vklop_c else 0
    return np.array([Ua, Ub, Uc])


B0 = trafo_data['B0']


mu0 = 4 * np.pi * 1e-7

B_orig = trafo_data['B_orig_pr'].copy()
H_orig = trafo_data['H_orig_pr'].copy()
B_zadnji = B_orig[-1]
H_zadnji = H_orig[-1]

k_sat = 1.0 
B_dodatne = np.array([2.00, 2.10])
H_dodatne = H_zadnji + (B_dodatne - B_zadnji) / (k_sat*mu0)

B_razsirjen = np.concatenate([B_orig, B_dodatne])
H_razsirjen = np.concatenate([H_orig, H_dodatne])

B_pozitivne = B_razsirjen[1:]
H_pozitivne = H_razsirjen[1:]
# B_pozitivne = B_orig[1:]
# H_pozitivne = H_orig[1:]


B_negativne = -np.flip(B_pozitivne)
H_negativne = -np.flip(H_pozitivne)

B_simetricen = np.concatenate([B_negativne, [0], B_pozitivne])
H_simetricen = np.concatenate([H_negativne, [0], H_pozitivne])

idxH = np.argsort(H_simetricen)
H_sim_sortH = H_simetricen[idxH]
B_sim_sortH = B_simetricen[idxH]

H_unique, B_unique_indices = np.unique(H_sim_sortH, return_index=True)
B_unique = B_sim_sortH[B_unique_indices]

idxB = np.argsort(B_simetricen)
B_sim_sortB = B_simetricen[idxB]
H_sim_sortB = H_simetricen[idxB]

B_unique_inv, H_unique_indices_inv = np.unique(B_sim_sortB, return_index=True)
H_unique_inv = H_sim_sortB[H_unique_indices_inv]

B_H_pchip_spline = PchipInterpolator(H_unique, B_unique, extrapolate = False)
dB_dH_pchip_spline = B_H_pchip_spline.derivative(nu=1)

H_B_pchip_spline = PchipInterpolator(B_unique_inv, H_unique_inv, extrapolate = False)
dH_dB_pchip_spline = H_B_pchip_spline.derivative(nu=1)

H_min, H_max = H_unique_inv[0], H_unique_inv[-1]
B_min, B_max = B_unique_inv[0], B_unique_inv[-1]

H0 = H_B_pchip_spline(B0)

# Initial conditions
fi0 = B0*trafo_data['Am']
fiz0_1 = fi0[0]
fiz0_2 = fi0[2]
fiz0 = np.array([fiz0_1, fiz0_2])

parametri = {
    # Dictionary unpacking
    **trafo_data,
    'm': m,
    'n': n,
    'Tp': Tp,
    'Ts': Ts,
    'T_omrezja': T_omrezja,
    'Tb': Tb,
    'Cm': Cm,
    'Npsi': Npsi,
    'Ntps': Ntps,
    'tip_bremena': tip_bremena,
    'Uomrezja_func': Uomrezja_func,
    'Uomrezja_vklop': Uomrezja_vklop,
    'omega': omega,
    't0': t0,
    't_vklop_a': t_vklop_a,
    't_vklop_b': t_vklop_b,
    't_vklop_c': t_vklop_c,
    'B_H_pchip_spline': B_H_pchip_spline,
    'dB_dH_pchip_spline': dB_dH_pchip_spline,
    'H_B_pchip_spline': H_B_pchip_spline,
    'dH_dB_pchip_spline': dH_dB_pchip_spline,
    'H_min': H_min,
    'H_max': H_max,
    'B_min': B_min,
    'B_max': B_max,
    'B0': B0,
    'H0': H0,
    'fiz0': fiz0
}