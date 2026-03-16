#---------------------------------------------------------------------------------------------------------------
# Author: Stela Ugrincic
# Affiliation: University of Maribor, Faculty of Electrical Engineering and Computer Science, Maribor, Slovenia
# Email: stela.ugrincic1@um.si
# Last modified: 16.03.2026
#---------------------------------------------------------------------------------------------------------------

# Script for running the simulation and defining the export settings
import sys, subprocess
from pathlib import Path

# Get the directory where this script is located
ROOT = Path(__file__).resolve().parent
PY   = sys.executable  # Uses the current Python environment

# Define the results directory relative to the project root
OUTPUT_DIR = ROOT/"results"
OUTPUT_DIR.mkdir(exist_ok=True) # Create folder if it doesn't exist
MAIN_SCRIPT = "main_transformer_model.py"
BASE_NAME = "MOD_Transformer_Inrush"

# Execution
cmd = [
    PY, str(ROOT / MAIN_SCRIPT),
    "--output-dir", str(OUTPUT_DIR),
    "--basename", BASE_NAME,
    # "--simple"        # Uncomment to use simplified model, otherwise the
                        # extended (full) model is being used by default
    # "--export-mat",   # Uncomment to export data to mat file
    # "--save-plots",   # Uncomment to save figures
    # "--no-show",      # Uncomment to skip displaying plots on the screen
]
subprocess.run(cmd, check=True)