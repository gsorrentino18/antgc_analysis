#!/usr/bin/env python3


# configure options
from argparse import ArgumentParser
parser.add_argument('--modelFilePath', type=str, help='BDT model file',
					default='/home/rusack/wadud/phoBDT/optimizedV0/aNTGC_photon_BDT.pkl', action='store')
parser.add_argument('--saveDir', type=str, default='/home/rusack/wadud/phoBDT/optimizedV0/',
					help='Save directory', action='store')

args = parser.parse_args()

from os import system
try:
	system('mkdir -p %s' % args.saveDir)
except OSError:
	print("\nCreation of results directory %s failed" % args.saveDir)
else:
	print("\nSuccessfully created directory %s " % args.saveDir)

import datetime
now = datetime.datetime.now()
logFileName = args.saveDir + '/predictBDT_' + now.strftime("%Y_%m_%d_%H_%M_%S") + ".log"

import sys



print("************************************************************************************************************************************************\n")
now = datetime.datetime.now()
print(now.strftime("%Y-%m-%d %H:%M:%S"))
print("Starting photon ID training BDT\n")
print("inFilePath\t=\t" + args.inFilePath)
print("modelFilePath\t=\t" + args.modelFilePath)
print("saveDir\t\t=\t" + args.saveDir)
print("inTreeName\t=\t" + args.inTreeName)
print("logFile\t\t=\t" + logFileName) 

BDTfeats = ["phoR9Full5x5", "phoS4Full5x5", "phoEmaxOESCrFull5x5", "phoE2ndOESCrFull5x5", "phoE2ndOEmaxFull5x5", "pho2x2OE3x3Full5x5", "phoE1x3OESCrFull5x5", "phoE2x5OESCrFull5x5", "phoE5x5OESCrFull5x5",
			"phoEmaxOE3x3Full5x5", "phoE2ndOE3x3Full5x5", "phoSigmaIEtaIEta", "phoSigmaIEtaIPhi", "phoSigmaIPhiIPhi", "phoSieieOSipipFull5x5", "phoEtaWidth", "phoPhiWidth", "phoEtaWOPhiWFull5x5"]
BDTfeats.sort()
print("\nBDT input features (" + str(len(BDTfeats)) + ") :")
print(BDTfeats)
print("\n")



phoMode = pickle.load(open(args.modelFilePath, "rb"))


‘weight’: the number of times a feature is used to split the data across all trees.

‘gain’: the average gain across all splits the feature is used in.

‘cover’: the average coverage across all splits the feature is used in.

‘total_gain’: the total gain across all splits the feature is used in.

‘total_cover’

#weights
print("\nFeature importances (" + str(len(phoModel.get_fscore())) + " features) with option weight:")
featImpDumpFilePath=args.saveDir + 'featImpsWeights.txt'
featImpDumpFile = open("%s" % featImpDumpFilePath, "w")
fScoreList = phoModel.get_score( importance_type='weight')
fScoreList = sorted(fScoreList.items(), key=lambda item: item[1], reverse=True)
usedFeatures = list()
for iFeat in fScoreList:
	print("\t\t" + iFeat[0] + "\t->\t" + str(iFeat[1]))
	featImpDumpFile.write(str(iFeat[0]) + ' , ' + str(iFeat[1]) + '\n')
	usedFeatures = usedFeatures + [str(iFeat[0])]
zeroImpFeats = list([str(x) for x in BDTfeats if x not in usedFeatures])
zeroImpFeats.sort()
for iFeat in zeroImpFeats:
	if iFeat is None:
		continue
	print("\t\t" + iFeat + "\t->\t0")
	featImpDumpFile.write(iFeat + ' , 0' + '\n')
featImpDumpFile.close()

#weights
print("\nFeature importances (" + str(len(phoModel.get_fscore())) + " features) with option ‘gain’:")
featImpDumpFilePath=args.saveDir + 'featImpsGain.txt'
featImpDumpFile = open("%s" % featImpDumpFilePath, "w")
fScoreList = phoModel.get_score( importance_type='‘gain’')
fScoreList = sorted(fScoreList.items(), key=lambda item: item[1], reverse=True)
usedFeatures = list()
for iFeat in fScoreList:
	print("\t\t" + iFeat[0] + "\t->\t" + str(iFeat[1]))
	featImpDumpFile.write(str(iFeat[0]) + ' , ' + str(iFeat[1]) + '\n')
	usedFeatures = usedFeatures + [str(iFeat[0])]
zeroImpFeats = list([str(x) for x in BDTfeats if x not in usedFeatures])
zeroImpFeats.sort()
for iFeat in zeroImpFeats:
	if iFeat is None:
		continue
	print("\t\t" + iFeat + "\t->\t0")
	featImpDumpFile.write(iFeat + ' , 0' + '\n')
featImpDumpFile.close()