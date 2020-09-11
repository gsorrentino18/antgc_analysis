#!/usr/bin/env python3

####################################################################
# Photon ID BDT for aNTGC search in Z(->nu nu) + Gamma channel     #
# Mohammad Abrar Wadud (BD), 08/15/2020                            #
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
parser.add_argument('--saveDir', type=str, default='tuningCV/',
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
logFileName = args.saveDir + '/tuneCVLog_' + \
	now.strftime("%Y_%m_%d_%H_%M_%S") + ".log"

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

scoreLogName = args.saveDir + '/tuneCVScoreLog_' +now.strftime("%Y_%m_%d_%H_%M_%S") + ".log"
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
from custom_auc import aucW, sampleStats
import numpy as np
import xgboost as xg
from hyperopt import fmin, tpe, hp, STATUS_OK, STATUS_FAIL, Trials, space_eval
from hyperopt.pyll_utils import hp_quniform as quniform
from hyperopt.pyll_utils import hp_uniformint as uniformint
import ROOT
import json
from random import randint

pd.set_option("display.max_rows", None, "display.max_columns", None)

BDTfeats = ["phoR9Full5x5", "phoS4Full5x5", "phoEmaxOESCrFull5x5", "phoE2ndOESCrFull5x5", "phoE2ndOEmaxFull5x5", "pho2x2OE3x3Full5x5", "phoE1x3OESCrFull5x5", "phoE2x5OESCrFull5x5", "phoE5x5OESCrFull5x5",
			"phoEmaxOE3x3Full5x5", "phoE2ndOE3x3Full5x5", "phoSigmaIEtaIEta", "phoSigmaIEtaIPhi", "phoSigmaIPhiIPhi", "phoSieieOSipipFull5x5", "phoEtaWidth", "phoPhiWidth", "phoEtaWOPhiWFull5x5"]
BDTfeats.sort()
print("\nBDT input features (" + str(len(BDTfeats)) + ") :")
print(BDTfeats)
print("\n\n")

space = {
			'objective': 'binary:logistic',
			'eval_metric': 'auc',
			'booster': 'gbtree',
			'sampling_method': 'gradient_based',
			'tree_method': 'gpu_hist',
			'predictor': 'gpu_predictor',
			'verbosity': 1,
			'nthread': -1,
			'seed' : randint(100, 100000),

			'max_depth': 10,
			'min_child_weight': hp.choice('min_child_weight', [30,35,40,60,80,100,120,140,160,180,200,300,400,500]),

			'alpha' : 0.,
			'lambda' : 0.,
			'gamma': 0.,
						
			'subsample': 0.6,
			'colsample_bytree':	0.8,

			'eta': 0.1,

			'scale_pos_weight': 1.
}


# Ntune = 1000
totalData = ROOT.RDataFrame(args.featTreeName, args.featFilePath)
indexData = ROOT.RDataFrame(args.indexTreeName, args.indexFilePath)
# totalData = totalData.Range(0, Ntune)
# indexData = indexData.Range(0, Ntune)
totalData = pd.DataFrame(totalData.AsNumpy(columns=BDTfeats + ['flatPtEtaRwNoXsec', 'isSignal']))
totalData['isTrain'] = pd.DataFrame(indexData.AsNumpy(columns=['isTrain'])).isTrain.values
totalData['isValidation'] = pd.DataFrame(indexData.AsNumpy(columns=['isValidation'])).isValidation.values

del indexData

print("\nTotal data:")
sampleStats(totalData)

print("\nTraining data:")
trainDF = totalData[(totalData.isTrain == 1)]
trainDF = trainDF.sample(frac = 1)
sampleStats(trainDF)

trainPho = xg.DMatrix(trainDF[BDTfeats].values, label=trainDF['isSignal'].values, weight=trainDF['flatPtEtaRwNoXsec'].values, feature_names=BDTfeats, nthread=-1)

del totalData
del trainDF

trials = Trials()

print("Search space:")
print(str(space))

def score(params):
	global trials
	if len(trials.trials)>1:
		global space
		for x in trials.trials[:-1]:
			space_point_index = dict([(key,value[0]) for key,value in x['misc']['vals'].items() if len(value)>0])
			if params == space_eval(space,space_point_index):
				return {'loss': 100000000., 'status': STATUS_FAIL}

	print("\n\nTraining with params : " + str(params))
	
	cv_results = xg.cv(params, dtrain=trainPho, num_boost_round=10000, nfold=5, metrics=['auc'], early_stopping_rounds=20, verbose_eval=True, show_stdv=True, seed=randint(100, 100000))

	mean_auc = cv_results['test-auc-mean'].max()
	boost_rounds = cv_results['test-auc-mean'].argmax()

	print("\tAUC {} for {} rounds".format(mean_auc, boost_rounds))
	print(str(cv_results.loc[[boost_rounds]]))

	theScore = (1. - mean_auc)

	scoreLogFile.write("\n\n")
	scoreLogFile.write(str(params) + "\n")
	scoreLogFile.write(str(cv_results) + "\n\n")

	return {'loss': theScore, 'status': STATUS_OK}

def optimize():
	global trials
	trials = Trials()
	now = datetime.datetime.now()
	print('\n\n' +now.strftime("%Y-%m-%d %H:%M:%S") + '\nLaunching tuning module... \n\n')
	best = fmin(score, space, algo=tpe.suggest, trials=trials, max_evals=1000)
	print("Optimum max_depth and min_child_weight:")
	print(best)

	scoreLogFile.write("Optimum:\n")
	scoreLogFile.write(str(best))

optimize()

# space = {
# 			'objective': 'binary:logistic',
# 			'eval_metric': 'auc',
# 			'booster': 'gbtree',
# 			'sampling_method': 'gradient_based',
# 			'tree_method': 'gpu_hist',
# 			'predictor': 'gpu_predictor',
# 			'verbosity': 1,
# 			'nthread': -1,
# 			'seed' : randint(100, 100000),

# 			'max_depth': 10,
# 			'min_child_weight': hp.choice('min_child_weight', [30,35,40,60,80,100,120,140,160,180,200,300,400,500]),

# 			'alpha' : 0.,
# 			'lambda' : 0.,
# 			'gamma': 0.,
						
# 			'subsample': 0.6,
# 			'colsample_bytree':	0.8,

# 			'eta': 0.1,

# 			'scale_pos_weight': 1.
# }
