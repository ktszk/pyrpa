from pathlib import Path
import numpy as np

# Resolve lib path relative to this file so import works from any current working directory.
_LIB_DIR = Path(__file__).resolve().parents[1]
_lib = np.ctypeslib.load_library("libfmod.so", str(_LIB_DIR))
