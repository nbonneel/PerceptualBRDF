# Python code for Optimal transport and energy metrics in the paper *Perceptual Quality of BRDF approximations; dataset and metrics*

The OT.py Python script computes all BRDF distances using Optimal Transport and Energy distances, for a single set of parameters.
This relies on the GeomLoss library. I could not make GeomLoss run under Windows despite significant effort (as of 2020), but this is easy to install on Linux.
Please edit BvqmFiles with paths to your BRDF directories.

Note that due to randomness in GeomLoss, results may slightly differ. Depending on random seed, the GeomLoss computation may occasionnally crash (e.g., due to an empty cluster in the K-means step). If that happens, please restart from the failed BRDF (uncomment line 161).
Displays correlation results in the console, and appends distances to a file ./Mosvalidation.csv