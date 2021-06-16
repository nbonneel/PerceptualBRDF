import string
import os.path

# ---- Directories and data location : put your own here and add them to Paths array

DistanceFiles = {
        '../Data/PerceptualDistancesWithoutGPLVM.csv'

}
MerlDirectories = {
	'../Data/merl/',
}

MerlApproxDirectories = {
        '../Data/merl_approx/',
}

BrdfListFiles = {
	'../Data/MaterialType.csv'
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
