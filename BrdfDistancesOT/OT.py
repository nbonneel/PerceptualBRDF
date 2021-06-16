import torch
import time

from geomloss import SamplesLoss  # See also ImagesLoss, VolumesLoss
import scipy.ndimage
import scipy.ndimage as ndimage
import skimage

import numpy as np
import os.path
import titopoh_loadRaw
from scipy.stats import pearsonr

##############################################
#### PARAMETERS (see paper for notations) ####
##############################################

alpha = 0.1     # scaling factor for the domain vs. the BRDF values
gamma1 = 0.25   # we take the final (red+green+blue)/3 value at the power gamma
gamma2 = 0.35   # we can test several gammas without recomputing the optimal transport ; I thus test 4 values at the same time.
gamma3 = 0.5
gamma4 = 1
loss = "sinkhorn"  # could be "sinkhorn" (optimal transport) or "energy" (Maximum Mean Discrepancies)
p = 2           # p = 2 means taking the squared Euclidean distance as ground metric, p = 1 means taking the Euclidean distance.
epsilon = 0.05  # entropic regularization
s = 0.8         # scaling factor in the multi-scale approach
downratio = 2   # downscaling BRDFs by a factor of 2 for performance reasons

eta = 0.01      # see function phi below
powert = 1./4.  # see function phi below


def phi(x):     # function to transform BRDF values
    y = np.power(np.maximum(x,0.), powert)
    #y = np.log(eta+np.maximum(x,0.))
    return y

##############################################

angles = np.arange(0,90)
cosCoefs = np.cos(angles*np.pi/180)
plr = np.arange(0,int(90/downratio))
cosCoefslr = np.cos(plr*np.pi/180*downratio)

cosines = np.tile(cosCoefs.reshape([90, 1, 1, 1]), [1, 90, 360, 3])
cosines1 = np.tile(cosCoefs.reshape([90, 1, 1]), [1, 90, 360])
cosineslr = np.tile(cosCoefslr.reshape([int(90/downratio), 1, 1]), [1, int(90/downratio), int(360/downratio)])
coords = np.zeros((int((90*90*360)/(downratio*downratio*downratio)),6))
for i in range(0,int(90/downratio)):
  for j in range(0,int(90/downratio)):
    for k in range(0,int(360/downratio)):
      id = i*int((90*360)/(downratio*downratio))+j*int(360/downratio)+k
      coords[id, 0] = alpha*(i+0.5)/90.*downratio
      coords[id, 1] = alpha*(j+0.5)/90.*downratio
      coords[id, 2] = alpha*(k+0.5)/360.*downratio
      coords[id, 3] = 0   # will store BRDF values later



use_cuda = torch.cuda.is_available()
print(use_cuda)
dtype    = torch.cuda.FloatTensor if use_cuda else torch.FloatTensor
N = int(90*90*360/(downratio*downratio*downratio))
M = N

##############################

import BvqmFiles as bf
liste_types = {'"DUPUY"':'T', '"GGX"':'X', '"BECKMAN"':'E', '"LAFORTUNE"':'L', '"WARD"':'Z', '"SGD"':'D', '"CHEBYCHEV"':'C', '"BAGHER"':'B',
'"LEGENDRE"':'G', '"ABC-OURS"':'F', '"LAFORTUNE-OURS"':'N', '"WARD-OURS"':'W', '"ABC"':'A', '"BLINN-PHONG"':'P', '"BLINN-PHONG-OURS"':'O',
'"GPLVM_APPROX.1"':'1', '"GPLVM_APPROX.2"':'2', '"GPLVM_APPROX.3"':'3', '"GPLVM_APPROX.4"':'4', '"GPLVM_APPROX.5"':'5',
'"GPLVM_APPROX.6"':'6'
}

distanceFile = bf.distanceFile()
MerlRefDirectory = bf.merlDirectory()
MerlApproxDirectory = bf.merlApproxDirectory()

print('using Distance File  : ' + distanceFile)
print('using Merl Directory : ' + MerlRefDirectory)
print('using Merl Approx    : ' + MerlApproxDirectory)


######## Clustering -- see GeomLoss examples ############

from pykeops.torch import generic_argmin

def KMeans(x, K=10, Niter=30, verbose=True):
    N, D = x.shape  # Number of samples, dimension of the ambient space

    # Define our KeOps CUDA kernel:
    nn_search = generic_argmin(  # Argmin reduction for generic formulas:
        'SqDist(x,y)',           # A simple squared L2 distance
        'ind = Vi(1)',           # Output one index per "line" (reduction over "j")
        'x = Vi({})'.format(D),  # 1st arg: one point per "line"
        'y = Vj({})'.format(D))  # 2nd arg: one point per "column"

    # K-means loop:
    # - x  is the point cloud,
    # - cl is the vector of class labels
    # - c  is the cloud of cluster centroids
    start = time.time()

    # Simplistic random initialization for the cluster centroids:
    perm = torch.randperm(N)
    idx = perm[:K]
    c = x[idx, :].clone()

    for i in range(Niter):
        cl  = nn_search(x,c).view(-1)  # Points -> Nearest cluster
        Ncl = torch.bincount(cl).type(dtype)  # Class weights
        for d in range(D):  # Compute the cluster centroids with torch.bincount:
            c[:, d] = torch.bincount(cl, weights=x[:, d]) / Ncl
    if use_cuda: torch.cuda.synchronize()
    end = time.time()
    #if verbose: print("KMeans performed in {:.3f}s.".format(end-start))

    return cl, c


######### Correct downsampling of BRDF values ###########
#### See https://github.com/scipy/scipy/issues/7324 #####

def resize_indices(input_shape: tuple, output_shape: tuple) -> np.ndarray:
  """Produces indices for resampling input shape into output shape."""
  if len(input_shape) != len(output_shape):
    raise ValueError('incompatible shape dimensions')

  indices = np.indices(output_shape, dtype=np.float32)
  for d in range(len(output_shape)):
    s = input_shape[d] / output_shape[d]
    indices[d] *= np.float32(s)
    indices[d] += np.float32(s / 2 - 0.5)  # input vs output pixel center displacement

  return indices


def resize(im: np.ndarray, shape: tuple, optimize_map: bool=True, **kwargs) -> np.ndarray:
  """Somehow "proper" implementation, instead of buggy scipy.ndimage.zoom"""
  if len(im.shape) != len(shape):
    raise ValueError('incompatible shape dimensions')

  if 'output' in kwargs:
    output = kwargs['output']
    del kwargs['output']
  else:
    output = np.empty(shape, dtype=im.dtype)

  if im.shape == shape:
    np.copyto(output, im)
    return output

  # optimize doing the interpolation for N-1 dimensions
  if optimize_map and im.shape[0] == shape[0]:
    indices = resize_indices(im.shape[1:], shape[1:])
    for i in range(shape[0]):
      map_coordinates(im[i], indices, output=output[i], **kwargs)
    return output

  return ndimage.map_coordinates(im, resize_indices(im.shape, shape), output=output, **kwargs)


################ Looping through all BRDFs to compute distances ############################

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
fd2 = open(distanceFile,'r')
index = 0
allvals = []
allrefs = []
for i,ligne in enumerate(fd2):
    approx, nom, distance = ligne.split(";")
    if approx != '"METHOD"': # Skipping the first line
        index = index+1
        #if index<=XXX:  #if the code crashes, skips the first XXX correctly computed values
        #  continue

        nom2 = nom[1:len(nom)-1] # Removing the quotes
        brdfRefPath = MerlRefDirectory + nom2 + '.titopo'        
        brdfApproxPath = MerlApproxDirectory + nom2 + '.' + liste_types[approx] + '.titopo'
        brdfName = nom2 + '.titopo'        
        brdfApproxName = nom2 + '.' + liste_types[approx] + '.titopo'

        # Loading files
        brdf1 = titopoh_loadRaw.LoadTitopo(brdfRefPath).reshape(90,90, 360,3)
        brdf2 = titopoh_loadRaw.LoadTitopo(brdfApproxPath).reshape(90,90, 360,3)

        wdist = [np.nan, np.nan, np.nan]
        for col in range(0, 3):	
          while np.isnan(wdist[col]):
            if downratio==1:
              brdf1lr = brdf1[:,:,:,col]*cosines1
              brdf2lr = brdf2[:,:,:,col]*cosines1
            else:
              brdf1lr = resize(brdf1[:,:,:,col]*cosines1, (int(90/downratio),int(90/downratio),int(360/downratio)), order=3)
              brdf2lr = resize(brdf2[:,:,:,col]*cosines1, (int(90/downratio),int(90/downratio),int(360/downratio)), order=3)
            values1 = brdf1lr.reshape(M)
            values2 = brdf2lr.reshape(N)

            coords1 = coords.copy()
            coords2 = coords.copy()
            coords1[:,3] = phi(values1)
            coords2[:,3] = phi(values2)		

            x = torch.from_numpy(coords1).float().to(device)
            y = torch.from_numpy(coords2).float().to(device)


            lab_i, c_i = KMeans(x, K= int(np.sqrt(M)) if use_cuda else 10)
            lab_j, c_j = KMeans(y, K= int(np.sqrt(N)) if use_cuda else 10)
            std_i = (( x - c_i[lab_i, :] )**2).sum(1).mean().sqrt()
            std_j = (( y - c_j[lab_j, :] )**2).sum(1).mean().sqrt()

            # MAIN COMPUTATION HERE
            loss = SamplesLoss(loss=loss, p=p, blur=epsilon, scaling = s, cluster_scale = max(std_i, std_j))
            a_i = torch.ones(N).type(dtype) / N
            b_j = torch.ones(M).type(dtype) / M

            #L = loss(x, y)  # By default, use constant weights = 1/number of samples
            L = loss( lab_i, a_i, x, lab_j, b_j, y )
            torch.cuda.synchronize()
            wdist[col] = L.cpu().numpy()

        

        line = brdfName + ";" + brdfApproxName + ";" + str(dist1) + ";" + str(dist2) + ";" + str(dist3)+ ";" + str(wdist[0])+ ";" + str(wdist[1])+ ";" + str(wdist[2]) + ";" + distance    
        allvals.append((np.maximum(0,wdist[0])+np.maximum(0,wdist[1])+np.maximum(0,wdist[2]))/3);
        allrefs.append(distance)
        corr1, _ = pearsonr(np.power(allvals, gamma1).astype(np.float), np.array(allrefs).astype(np.float)) 
        corr2, _ = pearsonr(np.power(allvals, gamma2).astype(np.float), np.array(allrefs).astype(np.float)) 
        corr3, _ = pearsonr(np.power(allvals, gamma3).astype(np.float), np.array(allrefs).astype(np.float) )
        corr4, _ = pearsonr(np.power(allvals, gamma4).astype(np.float), np.array(allrefs).astype(np.float) )
        print("Current Pearsons at line ",i," : ",  corr1, " ", corr2, " ", corr3, " ", corr4) # one for each value of gamma
        fd3 = open("./Mosvalidation.csv",'a+')        
        fd3.write(line)
        fd3.close()
        print(line)
print("Test finished")
fd2.close()
fd3.close()




