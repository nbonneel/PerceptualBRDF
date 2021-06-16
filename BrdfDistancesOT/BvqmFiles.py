import string
import os.path

# ---- Directories and data location : put your own here and add them to Paths array

DistanceFiles = {
        '../../Data/ExperimentResults/RealExperiment/AllResultsForLearningWithoutGPLVM.csv'

}
MerlDirectories = {
	'/media/nbonneel/Elements1/FINAL_BVQM/merl/',
	'D:/nbonneel/CNRS/csol/Data/titopoh/merl/',
        '../../Data/merl/',
        '/home/jpfarrug/Documents/Dev/bvqm/Data/merl/',
        '/Users/meros_k/Documents/Dev/MachineLearning/TensorFlowP3/Data/merl/',
        '/media/disc/csoler/MERL/',
}

MerlApproxDirectories = {
        '/media/nbonneel/Elements1/FINAL_BVQM/merl_approx/',
        'D:/nbonneel/CNRS/csol/Data/titopoh/merl_approx/',
        '../../Data/merl_approx/',
        '/home/jpfarrug/Documents/Dev/bvqm/Data/merl_approx/',
        '/Users/meros_k/Documents/Dev/MachineLearning/TensorFlowP3/Data/merl_approx/',
        '/media/disc2/csoler/Bvqm_Data/',    
}

BrdfListFiles = {
	'../../Data/MaterialType.csv'
}

def distanceFile():
	for file in DistanceFiles:
		if os.path.isfile(file) == True:
			return file
	
	return ""


def merlDirectory():
	for file in MerlDirectories:
		if os.path.isdir(file):
			return file
	
	return ""


def merlApproxDirectory():
	for file in MerlApproxDirectories:
		if os.path.isdir(file):
			return file
	
	return ""

def brdfListFile():
	for file in BrdfListFiles:
		if os.path.isfile(file) == True:
			return file
	
	return ""
