# Python code for metrics in the paper *Perceptual Quality of BRDF approximations; dataset and metrics*

This python code compute all nine metrics used in the paper *Perceptual Quality of BRDF approximations: dataset and metrics*, for a single BRDF or for the whole dataset. The input BRDF must be in a proprietary format, Titopo, which lineal samples the BRDF along the Theta_In, Theta_out and Phi_out angles, with the same number of samples than the binary (MERL) format. Binary files may be converted using the provided BrdfConvert code, written in C++. Please refer to the corresponding README for compilation and usage instructions. We also provide a python script to do this, but the conversion is very slow in this case.

## Usage :
### Computing all 9 metrics for two single BRDFs : 
In a python script of your choice, on in a notebook : 

First, import Numpy and the BRDFs loading/utility functions.
```python
import numpy as np
import os.path
import merl_nrec
import titopo_load_full
```

Second, load the BRDFs data:
If BRDFs are in binary format (beware : BRDF resampling is really slow in python !):
```python 
brdf1 = titopo_load_full.LoadTitopoFromMERL("../Data/pvc.binary")
brdf2 = titopo_load_full.LoadTitopoFromMERL("../Data/pvc.D.binary")
```

 If BRDFs are in titopo format (much faster):
```python
brdf1 = titopo_load_full.LoadTitopo("../Data/pvc.titopo »)
brdf2 = titopo_load_full.LoadTitopo("../Data/pvc.D.titopo »)
```

Finally, call the conversion function:
```python
from utility_functions import *
ComputeDistances(brdf1, brdf2, False, 80)
```
The computeDistance fonction outputs an array with 10 values : 9 are the metrics from the paper, in the same order, and the tenth is the Mean Opinion Score value from the crowdsourcing experiment. The third and fourth parameter of the ComputeDistance function controls, respectively, the initial cubic root applied to the data, and the grazing angle cutoff.

### Computing metrics for the whole dataset 
To compute metrics values for the whole dataset, you will need the file with the results of the crowdsourcing experiment (*MOS_Grace.csv*). You may also use your own if you comply to the same format.

The code is in the *mos_dist_compare.py* script, which exists in two versions depending if you have a titopo or a binary dataset. Modify the paths at the beginning of the script according to your own environment, then run 
```shell
python mos_dist_compare_titopo.py
```
 or
```shell
python mos_dist_compare_merl.py
```
 Please note that the latter has a very long computing time…
