import numpy as np
import os.path
import merl_nrec
import titopo_load_full

from utility_functions import *

brdf1 = titopo_load_full.LoadTitopoFromMERL("../Data/pvc.binary")
brdf2 = titopo_load_full.LoadTitopoFromMERL("../Data/pvc.D.binary")

# Computing distance between pvc and SGD approximated PVC, no cubic root applied, no angle cut-off
ComputeDistances(brdf1, brdf2, False, 90)

# Computing distance between pvc and SGD approximated PVC, cubic root applied, 80 degrees angle cut-off
ComputeDistances(brdf1, brdf2, True, 80)