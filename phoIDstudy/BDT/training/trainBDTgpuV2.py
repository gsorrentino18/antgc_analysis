#!/usr/bin/env python3

####################################################################
# Photon ID BDT for aNTGC search in Z(->nu nu) + Gamma channel     #
# Mohammad Abrar Wadud (BD), 04/02/2020                            #
####################################################################

# configure options
from argparse import ArgumentParser
parser = ArgumentParser(
	description='Photon ID BDT parameter tuning for CMS aNTGC search in Z(->nu nu) + gamma channel')
parser.add_argument('--featFilePath', type=str, help='Features root file',
					default='dataV2/mergedSamples.root', action='store')
parser.add_argument('--featTreeName', type=str, default='fullEBTree',
					help='Input tree name', action='store')
parser.add_argument('--indexFilePath', type=str, help='Index root file',
					default='dataV2/mergedSamplesShuffledIndices.root', action='store')
parser.add_argument('--indexTreeName', type=str, default='randomized_split',
					help='Input tree name', action='store')
parser.add_argument('--saveDir', type=str, default='trainingV2/',
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
timeStamp=now.strftime("%Y_%m_%d_%H_%M_%S")
logFileName = args.saveDir + '/trainingLog_' + timeStamp + ".log"

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
print(now.strftime("%Y-%m-%d %H:%M:%S") +
	  "\nStarting photon ID training BDT\n")

scoreLogName = args.saveDir + '/trainingScoreLog_' +timeStamp + ".log"
scoreLogFile = open(scoreLogName, "w")

print("featFilePath\t=\t" + args.featFilePath)
print("featTreeName\t=\t" + args.featTreeName)
print("indexFilePath\t=\t" + args.indexFilePath)
print("indexTreeName\t=\t" + args.indexTreeName)
print("saveDir\t\t=\t" + args.saveDir)
print("logFile\t\t=\t" + logFileName)


# using root_pandas - thanks to https://github.com/scikit-hep/root_pandas
import numpy as np
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
from random import randint

BDTfeats = ["phoR9Full5x5", "phoS4Full5x5", "phoEmaxOESCrFull5x5", "phoE2ndOESCrFull5x5", "pho2x2OE3x3Full5x5", "phoE1x3OESCrFull5x5", "phoE2x5OESCrFull5x5", "phoE5x5OESCrFull5x5",
			"phoSigmaIEtaIEta", "phoSigmaIEtaIPhi", "phoSigmaIPhiIPhi", "phoSieieOSipipFull5x5", "phoEtaWidth", "phoPhiWidth", "phoEtaWOPhiWFull5x5"]
BDTfeats.sort()
print("\nBDT input features (" + str(len(BDTfeats)) + ") :")
print(BDTfeats)

params = {
			'objective': 'binary:logistic',
			'eval_metric': 'auc',
			'booster': 'gbtree',
			'sampling_method': 'gradient_based',
			'tree_method': 'gpu_hist',
			'predictor': 'gpu_predictor',
			'verbosity': 1,
			'nthread': -1, 
			'seed': randint(10000,10000000),

			'scale_pos_weight': 1.,

			'max_depth': 4, 
			'min_child_weight': 5140, 

			'alpha': 572.9, 
			'lambda': 18.0,
			'gamma': 87.45,

			'subsample': 0.28, 
			'colsample_bytree': 0.67, 

			'eta': 0.0395
}


print("\nTraining with params:\n" + str(params))
print('\n')

# Ntune = 1000000
totalData = ROOT.RDataFrame(args.featTreeName, args.featFilePath)
indexData = ROOT.RDataFrame(args.indexTreeName, args.indexFilePath)
# totalData = totalData.Range(0, Ntune)
# indexData = indexData.Range(0, Ntune)
totalData = pd.DataFrame(totalData.AsNumpy(columns=BDTfeats + ['flatPtEtaRwNoXsec', 'isSignal']))
totalData['isTrain'] = pd.DataFrame(indexData.AsNumpy(columns=['isTrain'])).isTrain.values
totalData['isValidation'] = pd.DataFrame(indexData.AsNumpy(columns=['isValidation'])).isValidation.values
totalData['isTest'] = pd.DataFrame(indexData.AsNumpy(columns=['isTest'])).isTest.values
del indexData

print("\nTotal data:")
sampleStats(totalData)

print("\nTraining data:")
trainDF = totalData[(totalData.isTrain == 1)]
trainDF = trainDF.sample(frac = 1)
sampleStats(trainDF)

print("\nValidation data:")
validationDF = totalData[(totalData.isValidation == 1)]
validationDF = validationDF.sample(frac = 1)
sampleStats(validationDF)

print("\nTest data:")
testDF = totalData[(totalData.isTest == 1)]
testDF = testDF.sample(frac = 1)
sampleStats(testDF)


trainPho = xg.DMatrix(trainDF[BDTfeats].values, label=trainDF['isSignal'].values, weight=trainDF['flatPtEtaRwNoXsec'].values, feature_names=BDTfeats, nthread=-1)
del trainDF
validationPho = xg.DMatrix(validationDF[BDTfeats].values, label=validationDF['isSignal'].values, weight=validationDF['flatPtEtaRwNoXsec'].values, feature_names=BDTfeats, nthread=-1)
del validationDF
testPho = xg.DMatrix(testDF[BDTfeats].values, label=testDF['isSignal'].values, weight=testDF['flatPtEtaRwNoXsec'].values, feature_names=BDTfeats, nthread=-1)
del testDF

watchlist = [(trainPho, 'train'), (testPho, 'test'), (validationPho, "validation")]
accuracyProgress = dict()

phoModel = xg.train(params, dtrain=trainPho, evals=watchlist, evals_result=accuracyProgress, early_stopping_rounds=20, verbose_eval=True, num_boost_round=100000)
print("\nLearning during validation:\n")
print(json.dumps(accuracyProgress, indent=4, sort_keys=True))
scoreLogFile.write(json.dumps(accuracyProgress, indent=4, sort_keys=True))

accuracyProgress = dict()
phoModel = xg.train(params, dtrain=trainPho, evals=watchlist, evals_result=accuracyProgress, verbose_eval=True, num_boost_round=phoModel.best_ntree_limit)
print("\nLearning during final training:\n")
print(json.dumps(accuracyProgress, indent=4, sort_keys=True))
scoreLogFile.write(json.dumps(accuracyProgress, indent=4, sort_keys=True))

del trainPho
del validationPho
del testPho

now = datetime.datetime.now()
print('\nTraining complete @ ' + now.strftime("%Y-%m-%d %H:%M:%S") + "\n")

modelFileName=args.saveDir+'/aNTGC_photon_BDT_'+timeStamp
pickle.dump(phoModel, open('%s.pkl' % modelFileName, "wb"))
phoModel.save_model('%s.model' % modelFileName)
phoModel.dump_model(fout='%s.txt' % modelFileName, with_stats=True)
scoreLogFile.close()

print(now.strftime("%Y-%m-%d %H:%M:%S") + '\nModel saved as %s.(pkl/model/txt)' % modelFileName)
print('\nLearning progress saved in %s' % scoreLogName)



fScores = pd.DataFrame.from_dict(phoModel.get_score(importance_type='weight'), orient='index', columns=['weight'])
fScores.index.name = 'Feature'
fScores['gain'] = fScores.index.map(phoModel.get_score(importance_type='gain'))
fScores['cover'] = fScores.index.map(phoModel.get_score(importance_type='cover'))
fScores['total_gain'] = fScores.index.map(phoModel.get_score(importance_type='total_gain'))
fScores['total_cover'] = fScores.index.map(phoModel.get_score(importance_type='total_cover'))
fScores=fScores.sort_values(by=['gain'], ascending=False)
featImpDumpFilePath=args.saveDir + 'featImps.txt'
fScores.to_csv(featImpDumpFilePath, index=True)

print("\nFeature importances (" + str(len(fScores.index)) + " features) saved to : %s " % featImpDumpFilePath)
print(fScores)


print('\nGetting predictions @ ' + now.strftime("%Y-%m-%d %H:%M:%S") + "\n")

totalDM = xg.DMatrix(totalData[BDTfeats].values, label=totalData['isSignal'].values, weight=totalData['flatPtEtaRwNoXsec'].values, feature_names=BDTfeats, nthread=-1)
saveDF = totalData[['isSignal', 'isTrain', 'isValidation', 'isTest', 'flatPtEtaRwNoXsec']].copy()
del totalData
saveDF.columns = ['isSignalF', 'isTrainF', 'isValidationF', 'isTestF', 'flatPtEtaRwNoXsecF']
saveDF['isSignalF'] = saveDF['isSignalF'].astype('int')
saveDF['isTrainF'] = saveDF['isTrainF'].astype('int')
saveDF['isValidationF'] = saveDF['isValidationF'].astype('int')
saveDF['isTestF'] = saveDF['isTestF'].astype('int')

saveDF['bdtV2Score'] = phoModel.predict(totalDM)
del totalDM

outFileName = args.saveDir + '/BDTresults_' + timeStamp+ '.csv'
saveDF.to_csv(outFileName, encoding='utf-8', index=False)
print ('Data frame saved in %s' % outFileName)

saveDF = ROOT.RDF.MakeCsvDataFrame(outFileName);
outFileName = args.saveDir + '/' + 'BDTresults_' + timeStamp+ '.root'
saveDF.Snapshot('BDTresults', outFileName)
print ('Data frame saved in %s' % outFileName)

print("\n************************************************************************************************************************************************")
now = datetime.datetime.now()
print("Done @ " + now.strftime("%Y-%m-%d %H:%M:%S"))
print("************************************************************************************************************************************************")