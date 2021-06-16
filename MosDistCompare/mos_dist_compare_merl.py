import numpy as np
import os.path
import titopo_load_full
from utility_functions import ComputeDistances

# Parameters : replace with your own
limitAngle = 80
distanceFile = "./AllResultsForLearningWithoutGPLVM.csv"
resultFile = "./MosvalidationClamp80Cbrt.csv"
MerlRefDirectory = "/Volumes/Nouveau nom/bvqm/merl/"
MerlApproxDirectory = "/Volumes/Nouveau nom/bvqm/merl_approx/"

# Main code  : parsing data file and computint all distances for all brdf couples  

liste_types = {'"DUPUY"':'T', '"GGX"':'X', '"BECKMAN"':'E', '"LAFORTUNE"':'L', '"WARD"':'Z', '"SGD"':'D', '"CHEBYCHEV"':'C', '"BAGHER"':'B',
'"LEGENDRE"':'G', '"ABC-OURS"':'F', '"LAFORTUNE-OURS"':'N', '"WARD-OURS"':'W', '"ABC"':'A', '"BLINN-PHONG"':'P', '"BLINN-PHONG-OURS"':'O',
'"GPLVM_APPROX.1"':'1', '"GPLVM_APPROX.2"':'2', '"GPLVM_APPROX.3"':'3', '"GPLVM_APPROX.4"':'4', '"GPLVM_APPROX.5"':'5',
'"GPLVM_APPROX.6"':'6'
}

print('using Distance File  : ' + distanceFile)
print('using Merl Directory : ' + MerlRefDirectory)
print('using Merl Approx    : ' + MerlApproxDirectory)

# Reading distance file and computing various L2 distances between BRDFs pairs
fd2 = open(distanceFile,'r')
fd3 = open(resultFile,'w')

for i,ligne in enumerate(fd2):
    approx, nom, distance = ligne.split(";")
    if approx != '"METHOD"': # Skipping the first line
        brdfName = nom[1:len(nom)-1] # Removing the quotes
        brdfRefPath = MerlRefDirectory + brdfName + '.binary'        
        brdfApproxPath = MerlApproxDirectory + brdfName + '.' + liste_types[approx] + '.binary'      
        brdfApproxName = brdfName + '.' + approx[1:len(approx)-1]

        # Loading files
        brdf1 = titopo_load_full.LoadTitopoFromMERL(brdfRefPath)
        brdf2 = titopo_load_full.LoadTitopoFromMERL(brdfApproxPath)
        
        # Computing distances
        dist = ComputeDistances(brdf1, brdf2, True, limitAngle)

        print("In progress : ",i, end='\r')
        line = approx[1:len(approx)-1] + ";" + brdfName + ";" + str(dist[0]) + ";" + str(dist[1]) + ";" + str(dist[2]) + ";" + str(dist[3]) + ";" + str(dist[4]) + ";"  + str(dist[5]) + ";"  + str(dist[6]) + ";" + str(dist[7]) + ";" + str(dist[8]) + ";" + distance
        fd3.write(line)
print("Test finished")
fd2.close()
fd3.close()