#---------------------------------------------------------------------------------------------------------------
# Author: Stela Ugrincic
# Affiliation: University of Maribor, Faculty of Electrical Engineering and Computer Science, Maribor, Slovenia
# Email: stela.ugrincic1@um.si
# Last modified: 16.03.2026
#---------------------------------------------------------------------------------------------------------------

# MAIN TRANSFORMER MODEL SIMULATION SCRIPT
# Description: Solves the nonlinear system of differenial equations 
# for the dynamic model of a three-phase three-limb transformer; performes sampling; generates plots; exports results

import matplotlib.pyplot as plt
import argparse
from pathlib import Path
import numpy as np
import os
from scipy.integrate import solve_ivp
from transformer_model import transformer_ode_extended, transformer_ode_simplified
from transformer_model_configuration import parametri

# Command line Arguments
parser = argparse.ArgumentParser(description="Simulacija + izvoz + izris (ločen od main)")
parser.add_argument("--output-dir", default="./results", help="Mapa za izvoz slik/datotek")
parser.add_argument("--basename", default="run", help="Osnovno ime slik/datotek")
parser.add_argument("--save-plots", action="store_true", help="Če podano: slike shrani")
parser.add_argument("--no-show", dest="show_plots", action="store_false", help="Ne prikazuj grafov na zaslon")
parser.add_argument("--export-excel", action="store_true", help="Izvozi rezultate v Excel (.xlsx)")
parser.add_argument("--export-csv", action="store_true", help="Izvozi rezultate v CSV datoteke")
parser.add_argument("--export-mat", action="store_true", help="Izvozi rezultate v MATLAB .mat datoteko")
parser.add_argument("--simple", action="store_true", help="Uporabi poenostavljen model (brez Isz)")
args = parser.parse_args()

OUTPUT_DIR = Path(args.output_dir).resolve()
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
BASENAME   = args.basename
SAVE_PLOTS = args.save_plots
SHOW_PLOTS = args.show_plots if hasattr(args, 'show_plots') else True
USE_EXTENDED = not args.simple


# Help Functions
def savefig_maybe(fig, filename):
    if SAVE_PLOTS:
        fig.savefig(os.path.join(OUTPUT_DIR, filename), dpi=300)
def show_or_close(fig):
    if SHOW_PLOTS:
        plt.show()
    else:
        plt.close(fig)

# English: Simulation setup
# Slovenian: Nastavitve simulacije
t_start = 0.0
t_end = 2.0

m = parametri['m']
n = parametri.get('n', 0)

# English: Initialization and model setup
# Slovenian: Inicializacija začetnih pogojev in izbira modela
if USE_EXTENDED:
    y0 = np.concatenate([np.zeros(m + n), parametri['fiz0']])
    ode_function = lambda t, y: transformer_ode_extended(t, y, parametri)
else:
    transformer_ode_simplified.fiz_prev = parametri['fiz0'].copy()    
    y0 = np.concatenate([np.zeros(m), parametri['fiz0']])
    ode_function = lambda t, y: transformer_ode_simplified(t, y, parametri)

# English: Events (Switching)
# Slovenian: Definiranje funkcije dogodka, ki zazna trenutek vklopa t0
def vklop(t, y):
    return t - parametri['t0']
vklop.terminal = False
vklop.direction = 1

def vklop_faze_a(t, y):
    return t - parametri['t_vklop_a']
vklop_faze_a.terminal = False
vklop_faze_a.direction = 1

def vklop_faze_b(t, y):
    return t - parametri['t_vklop_b']
vklop_faze_b.terminal = False
vklop_faze_b.direction = 1

def vklop_faze_c(t, y):
    return t - parametri['t_vklop_c']
vklop_faze_c.terminal = False
vklop_faze_c.direction = 1

def dogodek_vklop_bremena(t, y):
    return t - parametri['t_vklop_a'] 
dogodek_vklop_bremena.terminal = False

dogodki = [vklop_faze_a, vklop_faze_b, vklop_faze_c, dogodek_vklop_bremena]

# English: Integration 
# Slovenian: Integracija
print("Starting simulation...")
sol = solve_ivp(ode_function, t_span = (t_start, t_end), y0 = y0,
                method='Radau', events = dogodki,
                rtol=1e-5, atol=1e-7, max_step=1e-4, dense_output = True)
print("Simulation finished.")

# Results unpacking

Cm = parametri['Cm']
Tp = parametri['Tp']
T_omrezja = parametri['T_omrezja']
Ts = parametri['Ts']
Tb = parametri['Tb']

dt = 1e-4 
t_min = sol.t[0]
t_max = sol.t[-1]
time = np.arange(t_min, t_max, dt)  
y_uniform = sol.sol(time)   

if USE_EXTENDED:
    Ipz = y_uniform[0:m, :].T              
    Isz = y_uniform[m:m+n, :].T            
    fiz = y_uniform[m+n:(m+n+2), :].T 

    Is = (Ts @ Isz.T).T                   
    Ib = (Tb @ Isz.T).T
    Us = (parametri['Rb'] @ Is.T).T
else:
    Ipz = y_uniform[0:m, :].T            
    fiz = y_uniform[m:m+2, :].T            

fi = (Cm.T @ fiz.T).T    
Ip = (Tp @ Ipz.T).T
Iom = (T_omrezja @ Ipz.T).T

time_plot = time - time[0]
Up_vals = np.array([parametri['Uomrezja_vklop'](ti) for ti in time])

# PLOTTING

plt.rcParams['font.family'] = 'Calibri'
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 0.5
color_l1 = '#0072BD'  
color_l2 = '#D95319'  
color_l3 = '#EDB120' 

print("\nGenerating plots..")

# Figure 1: Primary Voltages and Currents (Full view)
fig1, (ax1_1, ax1_2) = plt.subplots(2, 1, figsize=(10, 8))
fig1.suptitle('') 
# Primary phase voltages
ax1_1.plot(time_plot, Up_vals[:, 0], linewidth=1.25, label=r'$\mathit{u}_{\mathrm{p1}}$', color=color_l1)
ax1_1.plot(time_plot, Up_vals[:, 1], linewidth=1.25, label=r'$\mathit{u}_{\mathrm{p2}}$', color=color_l2)
ax1_1.plot(time_plot, Up_vals[:, 2], linewidth=1.25, label=r'$\mathit{u}_{\mathrm{p3}}$', color=color_l3)
ax1_1.set_ylabel(r'$\mathit{u}_{\mathrm{p}}$ (V)')
ax1_1.set_xlim(0, 2)
ax1_1.grid(True)
lgd1_1 = ax1_1.legend(loc='upper right')
lgd1_1.get_frame().set_linewidth(0.5)

# Primary phase currents
ax1_2.plot(time_plot, Ip[:, 0], linewidth=1.25, label=r'$\mathit{i}_{\mathrm{p1}}$', color=color_l1)
ax1_2.plot(time_plot, Ip[:, 1], linewidth=1.25, label=r'$\mathit{i}_{\mathrm{p2}}$', color=color_l2)
ax1_2.plot(time_plot, Ip[:, 2], linewidth=1.25, label=r'$\mathit{i}_{\mathrm{p3}}$', color=color_l3)
ax1_2.set_xlabel(r'$\mathit{t}$ (s)')
ax1_2.set_ylabel(r'$\mathit{i}_{\mathrm{p}}$ (A)')
ax1_2.set_xlim(0, 2)
ax1_2.grid(True)
lgd1_2 = ax1_2.legend(loc='upper right')
lgd1_2.get_frame().set_linewidth(0.5)

fig1.tight_layout(); savefig_maybe(fig1, f"{BASENAME}_p_2s.png"); show_or_close(fig1)

# Figure 2: Primary Voltages and Currents (0-500ms)
fig3, (ax3_1, ax3_2) = plt.subplots(2, 1, figsize=(10, 8))
fig3.suptitle('')

# Primary phase voltages
# Zgornji graf: Primarne napetosti
ax3_1.plot(time_plot, Up_vals[:, 0], linewidth=1.25, label=r'$\mathit{u}_{\mathrm{p1}}$', color=color_l1)
ax3_1.plot(time_plot, Up_vals[:, 1], linewidth=1.25, label=r'$\mathit{u}_{\mathrm{p2}}$', color=color_l2)
ax3_1.plot(time_plot, Up_vals[:, 2], linewidth=1.25, label=r'$\mathit{u}_{\mathrm{p3}}$', color=color_l3)
ax3_1.set_ylabel(r'$\mathit{u}_{\mathrm{p}}$ (V)')
ax3_1.set_xlim(0.05, 0.5)
ax3_1.grid(True)
lgd3_1 = ax3_1.legend(loc='upper right')
lgd3_1.get_frame().set_linewidth(0.5)

# Primary phase currents
# Spodnji graf: Primarni tokovi
ax3_2.plot(time_plot, Ip[:, 0], linewidth=1.25, label=r'$\mathit{i}_{\mathrm{p1}}$', color=color_l1)
ax3_2.plot(time_plot, Ip[:, 1], linewidth=1.25, label=r'$\mathit{i}_{\mathrm{p2}}$', color=color_l2)
ax3_2.plot(time_plot, Ip[:, 2], linewidth=1.25, label=r'$\mathit{i}_{\mathrm{p3}}$', color=color_l3)
ax3_2.set_xlabel(r'$\mathit{t}$ (s)')
ax3_2.set_ylabel(r'$\mathit{i}_{\mathrm{p}}$ (A)')
ax3_2.set_xlim(0.05, 0.5)
ax3_2.grid(True)
lgd3_2 = ax3_2.legend(loc='upper right')
lgd3_2.get_frame().set_linewidth(0.5)

fig3.tight_layout(); savefig_maybe(fig3, f"{BASENAME}_p_500ms.png"); show_or_close(fig3)

# Figure 3: # Primary phase currents (0-200ms)
fig7, ax7 = plt.subplots(figsize=(10, 4))
fig7.suptitle('')
ax7.plot(time_plot, Ip[:, 0], linewidth=1.25, label=r'$\mathit{i}_{\mathrm{p1}}$', color=color_l1)
ax7.plot(time_plot, Ip[:, 1], linewidth=1.25, label=r'$\mathit{i}_{\mathrm{p2}}$', color=color_l2)
ax7.plot(time_plot, Ip[:, 2], linewidth=1.25, label=r'$\mathit{i}_{\mathrm{p3}}$', color=color_l3)
ax7.set_xlabel(r'$\mathit{t}$ (s)')
ax7.set_ylabel(r'$\mathit{i}_{\mathrm{p}}$ (A)')
ax7.set_xlim(0.105, 0.25)
ax7.grid(True)
lgd7 = ax7.legend(loc='upper right')
lgd7.get_frame().set_linewidth(0.5)
fig7.tight_layout(); savefig_maybe(fig7, f"{BASENAME}_p_tok_200ms.png"); show_or_close(fig7)

# Figure 4: # Primary phase currents (Steadly state 1.8-2s)
fig8, ax8 = plt.subplots(figsize=(10, 4))
fig8.suptitle('')
ax8.plot(time_plot, Ip[:, 0], linewidth=1.25, label=r'$\mathit{i}_{\mathrm{p1}}$', color=color_l1)
ax8.plot(time_plot, Ip[:, 1], linewidth=1.25, label=r'$\mathit{i}_{\mathrm{p2}}$', color=color_l2)
ax8.plot(time_plot, Ip[:, 2], linewidth=1.25, label=r'$\mathit{i}_{\mathrm{p3}}$', color=color_l3)
ax8.set_xlabel(r'$\mathit{t}$ (s)')
ax8.set_ylabel(r'$\mathit{i}_{\mathrm{p}}$ (A)')
ax8.set_xlim(1.8, 2.0)
ax8.set_ylim(-0.5, 0.5)
ax8.grid(True)
lgd8 = ax8.legend(loc='upper right')
lgd8.get_frame().set_linewidth(0.5)
fig8.tight_layout(); savefig_maybe(fig8, f"{BASENAME}_p_tok_z200ms.png"); show_or_close(fig8)

print("\nPlots generated.")

# DATA EXPORT
if args.export_excel or args.export_csv or args.export_mat:
    try:
        import pandas as pd
    except ImportError:
        raise SystemExit("Pandas ni nameščen. Namesti s: pip install pandas openpyxl")

    df_primary = pd.DataFrame(
        np.column_stack([time, Up_vals, Ip, Iom]),
        columns=['t_s','Up1_V','Up2_V','Up3_V','Ip1_A','Ip2_A','Ip3_A','Iom1_A','Iom2_A','Iom3_A']
    )

    df_flux = pd.DataFrame(
        np.column_stack([time, fi, fi*1e3]),
        columns=['t_s','fi1_Wb','fi2_Wb','fi3_Wb','fi1_mVs','fi2_mVs','fi3_mVs']
    )

    df_secondary = None
    if USE_EXTENDED:
        df_secondary = pd.DataFrame(
            np.column_stack([time, Us, Is]),
            columns=['t_s','Us1_V','Us2_V','Us3_V','Is1_A','Is2_A','Is3_A']
        )

    # Excel
    if args.export_excel:
        xlsx_path = os.path.join(OUTPUT_DIR, f"{BASENAME}_rezultati.xlsx")
        with pd.ExcelWriter(xlsx_path) as writer:
            df_primary.to_excel(writer, sheet_name='primary', index=False)
            df_flux.to_excel(writer, sheet_name='flux', index=False)
            if df_secondary is not None:
                df_secondary.to_excel(writer, sheet_name='secondary', index=False)
        print(f"Excel file exported: {xlsx_path}")

    # CSV
    if args.export_csv:
        df_primary.to_csv(os.path.join(OUTPUT_DIR, f"{BASENAME}_primary.csv"), index=False)
        df_flux.to_csv(os.path.join(OUTPUT_DIR, f"{BASENAME}_flux.csv"), index=False)
        if df_secondary is not None:
            df_secondary.to_csv(os.path.join(OUTPUT_DIR, f"{BASENAME}_secondary.csv"), index=False)
        print("CSV files exported.")

    # MAT
    if args.export_mat:
        try:
            from scipy.io import savemat
        except ImportError:
            raise SystemExit("scipy ni nameščen. Namesti z: pip install scipy")
        mdict = {'time': time, 'Up': Up_vals, 'Ip': Ip, 'Iom': Iom, 't_fi': time, 'fi': fi}
        if USE_EXTENDED:
            mdict.update({'Us': Us, 'Is': Is})
        mat_path = os.path.join(OUTPUT_DIR, f"{BASENAME}_rezultati.mat")
        savemat(mat_path, mdict)
        print(f"MAT file exported: {mat_path}")