#!/usr/bin/env python3

####################################################################
# Photon ID BDT for aNTGC search in Z(->nu nu) + Gamma channel     #
# Mohammad Abrar Wadud (BD), 04/02/2020                            #
####################################################################

# configure options
from argparse import ArgumentParser
parser = ArgumentParser(
	description='Photon ID BDT Training for CMS aNTGC search in Z(->nu nu) + gamma channel')
parser.add_argument('--inFilePath', type=str, help='Input root file',
					default='/hdfs/cms/user/wadud/anTGC/BDTdata/mergedSamplesShuffled.root', action='store')
parser.add_argument('--modelFilePath', type=str, help='BDT model file',
					default='/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/BDT/data/trainingEta0p02//aNTGC_photon_BDT.pkl', action='store')
parser.add_argument('--saveDir', type=str, default='/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/BDT/data/trainingEta0p02/',
					help='Save directory', action='store')
parser.add_argument('--inTreeName', type=str, default='fullEB_Tree',
					help='Input tree name', action='store')
parser.add_argument('--chunksize', default=2000000,
					type=int, help='Chunk size', action='store')
parser.add_argument('--outTreeName', default='fullEB_BDT_Tree',
					type=str, help='Out tree name', action='store')

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


class Logger(object):

	def __init__(self):
		self.terminal = sys.stdout
		self.log = open("%s" % logFileName, "aw")

	def write(self, message):
		self.terminal.write(message)
		self.log.write(message)

sys.stdout = Logger()

print("************************************************************************************************************************************************\n")
now = datetime.datetime.now()
chunksize = args.chunksize
print(now.strftime("%Y-%m-%d %H:%M:%S"))
print("Starting photon ID training BDT\n")
print("inFilePath\t=\t" + args.inFilePath)
print("modelFilePath\t=\t" + args.modelFilePath)
print("saveDir\t\t=\t" + args.saveDir)
print("inTreeName\t=\t" + args.inTreeName)
print("chunksize\t=\t" + str(args.chunksize))
print("logFile\t\t=\t" + logFileName)

# load model and apply to test data
import pickle
import xgboost as xg
from progressbar import ProgressBar
from root_pandas.readwrite import read_root, to_root
from custom_auc import aucW, sampleStats

BDTfeats = ["phoR9Full5x5", "phoS4Full5x5", "phoEmaxOESCrFull5x5", "phoE2ndOESCrFull5x5", "phoE2ndOEmaxFull5x5", "phoE1x3OESCrFull5x5", "phoE2x5OESCrFull5x5", "phoE5x5OESCrFull5x5",
			"phoEmaxOE3x3Full5x5", "phoE2ndOE3x3Full5x5", "phoSigmaIEtaIEta", "phoSigmaIEtaIPhi", "phoSigmaIPhiIPhi", "phoSieieOSipipFull5x5", "phoEtaWidth", "phoPhiWidth", "phoEtaWOPhiWFull5x5"]
BDTfeats.sort()
print("\nBDT input features (" + str(len(BDTfeats)) + ") :")
print(BDTfeats)
print("\n")

print "\nRunning BDT predictor... " + now.strftime("%Y-%m-%d %H:%M:%S")

phoModelIterated = pickle.load(open(args.modelFilePath, "rb"))

now = datetime.datetime.now()
print "\nFeature importances (" + str(len(phoModelIterated.get_fscore())) + " features) :"
featImpDumpFile = open("%s/featImps.txt" % args.saveDir, "w")
fScoreList = phoModelIterated.get_fscore()
fScoreList = sorted(fScoreList.items(), key=lambda item: item[1], reverse=True)
usedFeatures = list()

for iFeat in fScoreList:
	print "\t\t" + iFeat[0] + "\t->\t" + str(iFeat[1])
	featImpDumpFile.write(str(iFeat[0]) + ' , ' + str(iFeat[1]) + '\n')
	usedFeatures = usedFeatures + [str(iFeat[0])]


zeroImpFeats = list([str(x) for x in BDTfeats if x not in usedFeatures])
zeroImpFeats.sort()


for iFeat in zeroImpFeats:
	if iFeat is None:
		continue
	print "\t\t" + iFeat + "\t->\t0"
	featImpDumpFile.write(iFeat + ' , 0' + '\n')

featImpDumpFile.close()

print('\nGetting predictions @ ' + now.strftime("%Y-%m-%d %H:%M:%S") + "\n")

iChunk = int(0)

for totalData in  read_root(paths=args.inFilePath, key=args.inTreeName, columns=(BDTfeats + ['PtEtaRwBG', 'xSecW', 'isSignal', 'isTrain', 'isValidation', 'splitRand']), chunksize=chunksize):

	print("\nTotal data:")
	totalData['bdtWeight'] = totalData['xSecW'] * totalData['PtEtaRwBG']
	sampleStats(totalData)

	allPho = xg.DMatrix(totalData[BDTfeats].values, label=totalData['isSignal'].values, weight=totalData['bdtWeight'].values, feature_names=BDTfeats, nthread=-1)

	saveDF = totalData[['isSignal', 'isTrain', 'isValidation', 'splitRand']].copy()
	saveDF.columns = ['isSignalF', 'isTrainF', 'isValidationF', 'splitRandF']
	saveDF['bdtScore'] = phoModelIterated.predict(allPho)


	outFileName = args.saveDir + '/' + 'BDTresults_' + str(iChunk) + '.root'
	saveDF.to_root('%s' % outFileName, key=args.outTreeName)
	print ('Data frame saved in %s)' % outFileName)

	print('Chunk ' + str(iChunk) + ' processed!\n')

	iChunk+=1


print "\n************************************************************************************************************************************************"
now = datetime.datetime.now()
print "Done @ " + now.strftime("%Y-%m-%d %H:%M:%S")
print
"************************************************************************************************************************************************"