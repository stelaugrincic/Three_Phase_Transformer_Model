#---------------------------------------------------------------------------------------------------------------
# Author: Stela Ugrincic
# Affiliation: University of Maribor, Faculty of Electrical Engineering and Computer Science, Maribor, Slovenia
# Email: stela.ugrincic1@um.si
# Last modified: 16.03.2026
#---------------------------------------------------------------------------------------------------------------

# THREE-PHASE DYNAMIC TRANSFORMER MODEL
# Description: This script defines the nonlinear system of differential equations 
# for a three-phase three-limb transformer

import numpy as np

def blok_matrika(Tp, Ts):
    """
    ENGLISH:
    Assembles a block diagonal matrix [Tp, 0; 0, Ts] for primary and secondary
    transformation matrices.

    SLOVENIAN:
    Sestavi matriko oblike [Tp, 0; 0, Ts] s pravilnimi dimenzijami,
    ne glede na število zančnih tokov (m in n).

    Args:
        Tp (np.ndarray): Primarna transformacijska matrika oblike (3, m).
        Ts (np.ndarray): Sekundarna transformacijska matrika oblike (3, n).

    Returns:
        np.ndarray: Sestavljena matrika oblike (6, m + n).
    """

    st_faz_p, m = Tp.shape
    st_faz_s, n = Ts.shape
    
    if st_faz_p != st_faz_s:
        raise ValueError(f"Število vrstic v Tp ({st_faz_p}) se ne ujema s številom vrstic v Ts ({st_faz_s}).")

    nic_zgoraj_desno = np.zeros((st_faz_p, n))
    nic_spodaj_levo = np.zeros((st_faz_s, m))
    
    T_blok = np.block([
        [Tp, nic_zgoraj_desno],
        [nic_spodaj_levo, Ts]
    ])

    return T_blok

def transformer_ode_extended(t, y, p):
    """
    ENGLISH: Extended dynamic model of a three-phase transformer including
    secondary currents.
    SLOVENIAN: Razširjeni model z integracijo magnetnih pretokov fiz
    """
    # Parameter extraction
    m, n = p['m'], p['n']
    Tp, Ts, Tb = p['Tp'], p['Ts'], p['Tb']
    lm, Am = p['lm'], p['Am']
    Cm = p['Cm']
    Ntps, Npsi = p['Ntps'], p['Npsi']
    Rp, Rs = p['Rp'], p['Rs']
    Lsp, Lss = p['Lsp'], p['Lss']
    tip_bremena, Rb, Lb = p['tip_bremena'], p['Rb'], p['Lb']
    Uomrezja_vklop = p['Uomrezja_vklop'](t)
    Tomrezja = p['T_omrezja']
    dH_dB_spline = p['dH_dB_pchip_spline']
    mu0 = 4 * np.pi * 1e-7

    # English: State vector unpacking
    # Slovenian: Razpakiranje vektorja stanja y v koraku integracije k:
    Ipz = y[0:m]
    Isz = y[m:m+n]
    Ips = np.concatenate([Tp @ Ipz, Ts @ Isz])
    Teta_ps = Ntps @ Ips

    fiz = y[m+n: m+n+2]

    fi = Cm.T @ fiz
    B = fi / Am

    dH_dB = dH_dB_spline(B)
        
    mask = (B < p['B_min']) | (B > p['B_max']) | np.isnan(dH_dB)
    if np.any(mask):
        dH_dB = np.where(mask, 1.0/mu0, dH_dB)

    mi_d = 1/dH_dB
    Rmd = np.diag(lm / (mi_d*Am))
    Rmdz = Cm @ Rmd @ Cm.T

    # English: Incremental Inductance Matrix
    # Slovenian: Izracun matrike inkrementalnih induktivnosti
    Ld = Npsi @ Cm.T @ np.linalg.inv(Rmdz) @ Cm @ Ntps
    Ldpp = Ld[0:3, 0:3]
    Ldps = Ld[0:3, 3:6]
    Ldsp = Ld[3:6, 0:3]
    Ldss = Ld[3:6, 3:6]

    # English: Load logic
    # Slovenian: Upostevanje tipa bremena
    if tip_bremena == 'prosti_tek':
        Rb = np.diag([1e9, 1e9, 1e9])
        Lb = np.zeros((3,3))
    elif tip_bremena == 'kratek_stik':
        Rb = np.zeros((3,3))
        Lb = np.zeros((3,3))
    elif tip_bremena == 'RL':
        Rb = p['Rb']
        Lb = p['Lb']
    else:
        Rb = np.diag([1e9, 1e9, 1e9])
        Lb = np.zeros((3,3))           
    
    # English: Calculation of the system matrix
    # Slovenian: Sestavljanje (izracun) matrik sistema
    Rpz = Tp.T @ (Rp + p['Rg']) @ Tp
    Rsz = Ts.T @ Rs @ Ts + Tb.T @ Rb @ Tb
    R_sistema = np.block([[Rpz, np.zeros((m, n))], [np.zeros((n, m)), Rsz]])

    Lpz = Tp.T @ (Lsp + Ldpp + p['Lg']) @ Tp
    Lpsz = Tp.T @ Ldps @ Ts
    Lspz = Ts.T @ Ldsp @ Tp
    Lsz = Ts.T @ (Lss + Ldss) @ Ts + Tb.T @ Lb @ Tb
    L_sistema = np.block([[Lpz, Lpsz], [Lspz, Lsz]])

    Upz = Tomrezja.T @ Uomrezja_vklop
    U_sistema = np.concatenate([Upz, np.zeros(n)])
    Iz = np.concatenate([Ipz, Isz])

    # English: Derivatives calculation
    # Slovenian: Izracun odvodov zancnih tokov
    # Stabilnejsi nacin racunanja kot iskanje inverzne matrike
    dIz_dt = np.linalg.solve(L_sistema, U_sistema - R_sistema @ Iz)
    T_ps = blok_matrika(Tp, Ts)
    dfiz_dt = np.linalg.inv(Rmdz) @ Cm @ Ntps @ T_ps @ dIz_dt
    
    return np.concatenate([dIz_dt, dfiz_dt])

def transformer_ode_simplified(t, y, p):
    """
    ENGLISH: Simplified model (Primary currents only), can be used for no-load energization
    SLOVENIAN: Poenostavljen model za vklopni pojav (Isz = 0)
    """

    # Parameter extraction
    # Pridobitev parametrov iz slovarja 'p'
    m = p['m']
    Tp, Cm = p['Tp'], p['Cm']
    lm, Am = p['lm'], p['Am']
    Rp, Lsp, Rg, Lg = p['Rp'], p['Lsp'], p['Rg'], p['Lg']
    Ntps, Npsi = p['Ntps'], p['Npsi']
    Uomrezja_vklop = p['Uomrezja_vklop'](t)
    Tomrezja = p['T_omrezja']
    H_B_spline = p['H_B_pchip_spline']
    dH_dB_spline = p['dH_dB_pchip_spline']
    mu0 = 4 * np.pi * 1e-7

    Ipz = y[0: m]
    Ip = Tp @ Ipz
    Teta_ps = Ntps[:, 0:3] @ Ip

    fiz = y[m: m+2]

    fi = Cm.T @ fiz
    B = fi / Am

    dH_dB = dH_dB_spline(B)

    mask = (B < p['B_min']) | (B > p['B_max']) | np.isnan(dH_dB)
    if np.any(mask):
        dH_dB = np.where(mask, 1.0/mu0, dH_dB)

    mi_d = 1/dH_dB

    Rmd = np.diag(lm / (mi_d*Am))
    Rmdz = Cm @ Rmd @ Cm.T

    Ld = Npsi @ Cm.T @ np.linalg.inv(Rmdz) @ Cm @ Ntps
    Ldpp = Ld[0:3, 0:3]

    Rpz = Tp.T @ (Rp + Rg) @ Tp
    Lpz = Tp.T @ (Lsp + Ldpp + Lg) @ Tp

    # English: Derivatives calculation
    # Izracun casovnega odvoda primarnih tokov zank
    dIpz_dt = np.linalg.solve(Lpz, Tomrezja.T @ Uomrezja_vklop - Rpz @ Ipz)
    dfiz_dt = np.linalg.inv(Rmdz) @ Cm @ Ntps[:, 0:3] @ Tp @ dIpz_dt


    # Diagnostic block (optional)
    # Ta blok periodično izpisuje ključne vrednosti za lažje iskanje napak.
    
    if not hasattr(transformer_ode_simplified, 'counter'):
        transformer_ode_simplified.counter = 0
    if transformer_ode_simplified.counter % 1000 == 0:
        L_magnetenja = Ldpp[0,0]
        Z_magnetenja = p['omega'] * L_magnetenja

        H_skupno = p['H_B_pchip_spline'](B)
        mi_d = p['lm'][0] / (Rmd[0,0] * p['Am'][0])
        
        print(f"\n--- DIAGNOSTIKA @ t = {t:.6f} (s) ---")

        with np.printoptions(precision=4, suppress=True):
            print(f"  Vhodni tok (stanje) Ipz (A): {y}")
            print(f"  Skupna gostota magnetnega pretoka B (T):         \n{B}")
            print(f"  Skupna jakost mag. polja (A/m):           \n{H_skupno}")
            print(f"  Inkrementalna mag. permeabilnost (Vs/Am): \n{mi_d:.4e}")
            print(f"  Izracunana L_mag (H):             {L_magnetenja:.4f}")
            print(f"  Impedanca wL_mag (Ohm):           {Z_magnetenja:.2f}")
            print(f"  Matrika Ldpp (magnetilna) (H):    \n{Ldpp}")
            print(f"  Matrika Lpz (zančna) (H):         \n{Lpz}")
            print(f"  Koncni odvod dIpz/dt:             {dIpz_dt}")
    
    transformer_ode_simplified.counter += 1

    return np.concatenate([dIpz_dt, dfiz_dt])