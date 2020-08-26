#!/usr/bin/env python3

####################################################################
# Photon ID BDT for aNTGC search in Z(->nu nu) + Gamma channel     #
# Mohammad Abrar Wadud (BD), 04/02/2020                            #
####################################################################

import numpy as np

# configure options
from argparse import ArgumentParser
parser = ArgumentParser(
	description='Photon ID BDT Training for CMS aNTGC search in Z(->nu nu) + gamma channel')
parser.add_argument('--inFilePath', type=str, help='Input root file',
					default='/hdfs/cms/user/wadud/anTGC/BDTdata/mergedSamplesShuffled.root', action='store')
parser.add_argument('--saveDir', type=str, default='/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/BDT/data/test/',
					help='Save directory', action='store')
parser.add_argument('--inTreeName', type=str, default='fullEB_Tree',
					help='Input tree name', action='store')
parser.add_argument('--chunksize', default=10000,
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
logFileName = args.saveDir + '/trainWithTune_' + now.strftime("%Y_%m_%d_%H_%M_%S") + ".log"

validationProgressFilename = args.saveDir + '/validationProgress.log'
validationProgressFile = open("%s" % validationProgressFilename,"w")

testProgressFilename = args.saveDir + '/testProgress.log'
testProgressFile = open("%s" % testProgressFilename,"w")


import sys


class Logger(object):

	def __init__(self):
		self.terminal = sys.stdout
		self.log = open("%s" % logFileName, "a")

	def write(self, message):
		self.terminal.write(message)
		self.log.write(message)

	def flush(self):
		pass

sys.stdout = Logger()


print("************************************************************************************************************************************************\n")
now = datetime.datetime.now()
chunksize = args.chunksize
print(now.strftime("%Y-%m-%d %H:%M:%S"))
print("Starting photon ID training BDT\n")
print("inFilePath\t=\t" + args.inFilePath)
print("saveDir\t\t=\t" + args.saveDir)
print("inTreeName\t=\t" + args.inTreeName)
print("chunksize\t=\t" + str(args.chunksize))
print("logFile\t\t=\t" + logFileName)


# using root_pandas - thanks to https://github.com/scikit-hep/root_pandas
import pandas as pd
import xgboost as xg
import pickle
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib as mpl
from matplotlib import rcParams
from matplotlib import rc
import matplotlib.pyplot as plt
from custom_auc import aucW, sampleStats
import json
import ROOT
import math


BDTfeats = ["phoR9Full5x5", "phoS4Full5x5", "phoEmaxOESCrFull5x5", "phoE2ndOESCrFull5x5", "phoE2ndOEmaxFull5x5", "phoE1x3OESCrFull5x5", "phoE2x5OESCrFull5x5", "phoE5x5OESCrFull5x5",
			"phoEmaxOE3x3Full5x5", "phoE2ndOE3x3Full5x5", "phoSigmaIEtaIEta", "phoSigmaIEtaIPhi", "phoSigmaIPhiIPhi", "phoSieieOSipipFull5x5", "phoEtaWidth", "phoPhiWidth", "phoEtaWOPhiWFull5x5"]
BDTfeats.sort()
print("\nBDT input features (" + str(len(BDTfeats)) + ") :")
print(BDTfeats)

params = {
			'objective': 'binary:logistic',
			'eta': 0.02,
			'eval_metric': 'auc',
			'booster': 'gbtree',
			'verbosity': 1,
			'nthread': -1,
			'gamma': 0.05,
			'max_depth': 22,
			'min_child_weight': 1.8,
			'subsample': 0.5,
			'colsample_bytree': 0.5,
			'alpha' : 0.01,
			'sampling_method': 'uniform',
			'tree_method': 'exact',
			'predictor': 'cpu_predictor'
		  }

print("Training with params:\n" + str(params))
print('\n')

phoModelIterated = None

totalData = ROOT.RDataFrame(args.inTreeName, args.inFilePath)
nEntries =  int(totalData.Count().GetValue())
nChunks = int(math.ceil(float(nEntries)/float(args.chunksize)))

print('\nChunks = %d/%d = %d\n' % (nEntries, chunksize, nChunks))

for iChunk in range(nChunks):

	iStart = iChunk * chunksize
	iEnd = (iChunk+1) * chunksize	
	iChunkData = totalData.Range(iStart, iEnd)
	print('Processing chunk '+str(iChunk)+': range '+str(iStart)+' to ' + str(iStart+iChunkData.Count().GetValue()))

	iChunkData = iChunkData.Define("isSignalB", "int(isSignal)")
	iChunkData = iChunkData.Define("isTrainB", "int(isTrain)")
	iChunkData = iChunkData.Define("isValidationB", "int(isValidation)")

	iChunkData = pd.DataFrame(iChunkData.AsNumpy(columns=BDTfeats + ['PtEtaRwBG', 'xSecW', 'isSignalB', 'isTrainB', 'isValidationB']))

	iChunkData.rename(columns={'isSignalB': 'isSignal', 'isTrainB': 'isTrain', 'isValidationB': 'isValidation'}, inplace=True)

	print("\nTotal data:")
	iChunkData['bdtWeight'] = iChunkData['xSecW'] * iChunkData['PtEtaRwBG']
	sampleStats(iChunkData)

	print("\nTraining data:")
	trainDF = iChunkData[(iChunkData.isTrain == 1) & (iChunkData.isValidation == 0)]
	sampleStats(trainDF)

	print("\nValidation data:")
	validationDF = iChunkData[(iChunkData.isTrain == 1) & (iChunkData.isValidation == 1)]
	sampleStats(validationDF)

	print("\nTest data:")
	testDF = iChunkData[(iChunkData.isTrain == 0)]
	sampleStats(testDF)

	print('\n')

	trainPho = xg.DMatrix(trainDF[BDTfeats].values, label=trainDF['isSignal'].values, weight=trainDF['bdtWeight'].values, feature_names=BDTfeats, nthread=-1)
	validationPho = xg.DMatrix(validationDF[BDTfeats].values, label=validationDF['isSignal'].values, weight=validationDF['bdtWeight'].values, feature_names=BDTfeats, nthread=-1)
	testPho = xg.DMatrix(testDF[BDTfeats].values, label=testDF['isSignal'].values, weight=testDF['bdtWeight'].values, feature_names=BDTfeats, nthread=-1)
	watchlist = [(trainPho, 'train'), (testPho, 'test'), (validationPho, "validation")]
	accuracyProgress = dict()
	phoModel = xg.train(params, dtrain=trainPho, evals=watchlist, evals_result=accuracyProgress, early_stopping_rounds=10, verbose_eval=True, num_boost_round=100000, xgb_model=phoModelIterated)

	del trainDF
	del validationDF
	del trainPho
	del validationPho

	print("\nValidation progress:\n")
	print(json.dumps(accuracyProgress, indent=4, sort_keys=True))
	validationProgressFile.write(json.dumps(accuracyProgress, indent=4, sort_keys=True))
	print('\n')

	trainDF = iChunkData[(iChunkData.isTrain == 1)]
	trainPho = xg.DMatrix(trainDF[BDTfeats].values, label=trainDF['isSignal'].values, weight=trainDF['bdtWeight'].values, feature_names=BDTfeats, nthread=-1)
	accuracyProgress = dict()
	watchlist = [(trainPho, 'train'), (testPho, 'test')]
	phoModelIterated = xg.train(params, dtrain=trainPho, evals=watchlist, evals_result=accuracyProgress, num_boost_round=phoModel.best_iteration + 1, xgb_model=phoModelIterated)

	print("\nPost-validation progress:\n")
	print(json.dumps(accuracyProgress, indent=4, sort_keys=True))
	testProgressFile.write(json.dumps(accuracyProgress, indent=4, sort_keys=True))
	print('\n')

	print('Chunk ' + str(iChunk) + ' processed!\n')

	if(iChunk == 2): break



now = datetime.datetime.now()
print('\nTraining complete @ ' + now.strftime("%Y-%m-%d %H:%M:%S") + "\n")

pickle.dump(phoModelIterated, open('%s/aNTGC_photon_BDT.pkl' % args.saveDir, "wb"))
phoModelIterated.save_model('%s/aNTGC_photon_BDT.model' % args.saveDir)
phoModelIterated.dump_model(fout='%s/modelDump.txt' % args.saveDir)

validationProgressFile.close()
testProgressFile.close()

print now.strftime("%Y-%m-%d %H:%M:%S") + '\nModel saved as %s/aNTGC_photon_BDT.(pkl/model/txt)' % args.saveDir
print('\nValidation and test progress saved in %s/validatiobProgress.log and testProgress.log/' % args.saveDir)

print("\nFeature importances (" + str(len(phoModelIterated.get_fscore())) + " features) :")
featImpDumpFile = open("%s/featImps.txt" % args.saveDir, "w")
fScoreList = phoModelIterated.get_fscore()
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

print('\nGetting predictions @ ' + now.strftime("%Y-%m-%d %H:%M:%S") + "\n")

for iChunk in range(nChunks):

	iStart = iChunk * chunksize
	iEnd = (iChunk+1) * chunksize	
	iChunkData = totalData.Range(iStart, iEnd)
	print('Processing chunk '+str(iChunk)+': range '+str(iStart)+' to ' + str(iChunkData.Count().GetValue()))

	iChunkData = iChunkData.Define("isSignalB", "int(isSignal)")
	iChunkData = iChunkData.Define("isTrainB", "int(isTrain)")
	iChunkData = iChunkData.Define("isValidationB", "int(isValidation)")

	iChunkData = pd.DataFrame(iChunkData.AsNumpy(columns=BDTfeats + ['PtEtaRwBG', 'xSecW', 'isSignalB', 'isTrainB', 'isValidationB', 'splitRand']))

	iChunkData['isSignalB'] = iChunkData['isSignalB'].astype(bool)
	iChunkData['isTrainB'] = iChunkData['isTrainB'].astype(bool)
	iChunkData['isValidationB'] = iChunkData['isValidationB'].astype(bool)
	iChunkData.rename(columns={'isSignalB': 'isSignal', 'isTrainB': 'isTrain', 'isValidationB': 'isValidation'}, inplace=True)

	iChunkData['bdtWeight'] = iChunkData['xSecW'] * iChunkData['PtEtaRwBG']
	sampleStats(iChunkData)

	iChunkDM = xg.DMatrix(iChunkData[BDTfeats].values, label=iChunkData['isSignal'].values, weight=iChunkData['bdtWeight'].values, feature_names=BDTfeats, nthread=-1)

	saveDF = iChunkData[['isSignal', 'isTrain', 'isValidation', 'bdtWeight', 'splitRand']].copy()
	saveDF.columns = ['isSignalF', 'isTrainF', 'isValidationF', 'bdtWeightF', 'splitRandF']
	saveDF['bdtScore'] = phoModelIterated.predict(iChunkDM)

	del iChunkData

	outFileName = args.saveDir + '/' + 'BDTresults_' + str(iChunk) + '.csv'
	saveDF.to_csv(outFileName, encoding='utf-8', index=False)
	print ('Data frame saved in %s' % outFileName)

	del saveDF

	saveDF = ROOT.RDF.MakeCsvDataFrame(outFileName);
	outFileName = args.saveDir + '/' + 'BDTresults_' + str(iChunk) + '.root'
	saveDF.Snapshot(args.outTreeName, outFileName)
	print ('Data frame saved in %s' % outFileName)

	print('Chunk ' + str(iChunk) + ' processed!\n')

	if(iChunk == 2): break


print("\n************************************************************************************************************************************************")
now = datetime.datetime.now()
print("Done @ " + now.strftime("%Y-%m-%d %H:%M:%S"))
print("************************************************************************************************************************************************")
