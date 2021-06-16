import numpy as np

def fLAB(t):
    tempVal1 = np.power((24/116),3)
    tempVal2 = 84./108.
    tempVal3 = 16/116

    return np.less(t,tempVal1)*np.cbrt(t) + np.multiply(np.greater_equal(t,tempVal1),np.multiply(tempVal2,t)+tempVal3)

def RGB2LAB(brdf):
    XYZMatrix = np.array([[0.412453,0.357580,0.180423],[0.212671,0.715160,0.072169],[0.019334,0.119193,0.950227]])
    D65Point = np.array([95.047,100.0,108.883])
    
    brdfXYZ = np.matmul(brdf,XYZMatrix)
    
    labBRDF = np.empty((90,90,360,3))
    labBRDF[:,:,:,0] = 116 * (brdfXYZ[:,:,:,1]/D65Point[1])-16
    labBRDF[:,:,:,1] = 500 * (fLAB(brdfXYZ[:,:,:,0]/D65Point[0]) - fLAB(brdfXYZ[:,:,:,1]/D65Point[1]))
    labBRDF[:,:,:,2] = 200 * (fLAB(brdfXYZ[:,:,:,1]/D65Point[1]) - fLAB(brdfXYZ[:,:,:,2]/D65Point[2]))
    
    return(labBRDF)
    
# Distance functions
def ComputeDistances(brdf1, brdf2, applyCbrt, plimitAngle):
    dist = np.empty(9)

    p = np.arange(0,plimitAngle)
    cosCoefs = np.cos(p*np.pi/180)
    sinCoefs = np.sin(p*np.pi/180)
    cosThI = np.tile(cosCoefs.reshape([plimitAngle, 1, 1, 1]), [1, plimitAngle, 360, 3])
    cosThO = np.tile(cosCoefs.reshape([1, plimitAngle, 1, 1]), [plimitAngle, 1, 360, 3])
    sinThO = np.tile(sinCoefs.reshape([1, plimitAngle, 1, 1]), [plimitAngle, 1, 360, 3])

    brdf1Lab = RGB2LAB(brdf1)
    brdf2Lab = RGB2LAB(brdf2)
    
    # restricting angles if needed
    if(plimitAngle < 90):
        brdf1 = brdf1[0:plimitAngle,0:plimitAngle,:,:]
        brdf2 = brdf2[0:plimitAngle,0:plimitAngle,:,:]
        brdf1Lab = brdf1Lab[0:plimitAngle,0:plimitAngle,:,:]
        brdf2Lab = brdf2Lab[0:plimitAngle,0:plimitAngle,:,:]
    
    if(applyCbrt == True):
        brdf1 = np.cbrt(brdf1)
        brdf2 = np.cbrt(brdf2)
        brdf1Lab = np.cbrt(brdf1Lab)
        brdf2Lab = np.cbrt(brdf2Lab)

    myEpsilon = 10**(-3)
    
    cosMap = np.maximum(cosThI*cosThO, myEpsilon)
    brdf1Safe = np.maximum(brdf1*cosMap + myEpsilon, myEpsilon)
    brdf2Safe = np.maximum(brdf2*cosMap + myEpsilon, myEpsilon)
                    
    # Computing L2 Distances
    distA = np.sum(np.square(brdf2-brdf1)) 
    distB = np.sum(np.square((brdf2-brdf1)*cosThI)) # Same as NGan
    distC = np.sum(np.cbrt(np.square((brdf2-brdf1)*cosThI))) # Fores
    distD = np.sum(np.square(brdf2-brdf1)*cosThO*sinThO) # Our L2 distance
    # Additionnal distances
    distE = np.sum(np.abs(brdf2-brdf1)) # L1 Distance
    distF = np.sum(np.square(np.abs(brdf2Lab-brdf1Lab))*cosThO*cosThI) # From Ryman 2018
    # From Tianchen et al. 2018
    distG = np.sum(np.cbrt(np.square((brdf2-brdf1)*cosMap))) 
    distH = np.sum(np.abs(np.log(brdf1Safe/brdf2Safe)))
    distI = np.sum(np.square(np.log(brdf1Safe/brdf2Safe)))
    
    nbSamples = plimitAngle*plimitAngle*360
    
    dist[0] = np.sqrt(distA/(nbSamples))
    dist[1] = np.sqrt(distB/(nbSamples))
    dist[2] = np.sqrt(distC/(nbSamples))
    dist[3] = np.sqrt(distD/(nbSamples))
    dist[4] = distE/(nbSamples)
    dist[5] = np.sqrt(distF/(nbSamples))
    dist[6] = distG/(nbSamples)
    dist[7] = distH/(nbSamples)
    dist[8] = np.sqrt(distI/(nbSamples))

    return dist
