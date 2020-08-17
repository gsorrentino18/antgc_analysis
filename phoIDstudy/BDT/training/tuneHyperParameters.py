#!/usr/bin/env python3

####################################################################
# Photon ID BDT for aNTGC search in Z(->nu nu) + Gamma channel     #
# Mohammad Abrar Wadud (BD), 08/15/2020                            #
####################################################################


# configure options
from argparse import ArgumentParser
parser = ArgumentParser(
	description='Photon ID BDT parameter tuning for CMS aNTGC search in Z(->nu nu) + gamma channel')
parser.add_argument('--inFilePath', type=str, help='Input root file',
					default='/hdfs/cms/user/wadud/anTGC/BDTdata/mergedSamplesShuffled.root', action='store')
parser.add_argument('--saveDir', type=str, default='/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/BDT/data/training/',
					help='Save directory', action='store')
parser.add_argument('--testSize', default=0.2, type=float,
					help='Fraction of events for testing', action='store')
parser.add_argument('--validSize', default=0.25, type=float,
					help='Fraction of events for overtraing reduction', action='store')
parser.add_argument('--inTreeName', type=str, default='fullEB_Tree',
					help='Input tree name', action='store')
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
logFileName = args.saveDir + '/trainLog_' + \
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

scoreLogName = args.saveDir + '/tuneScoreLog_' +now.strftime("%Y_%m_%d_%H_%M_%S") + ".log"
scoreLogFile = open(scoreLogName, "wb")

testSize = args.testSize
validSize = args.validSize
print("inFilePath\t=\t" + args.inFilePath)
print("saveDir\t\t=\t" + args.saveDir)
print("inTreeName\t=\t" + args.inTreeName)
print("testSize\t=\t" + str(testSize))
print("validSize\t=\t" + str(validSize))
print("logFile\t\t=\t" + logFileName)


# using root_pandas - thanks to https://github.com/scikit-hep/root_pandas
import numpy as np
import pandas as pd
from custom_auc import aucW, sampleStats
from root_pandas import read_root, to_root
from progressbar import ProgressBar
import numpy as np
import xgboost as xg
from hyperopt import fmin, tpe, hp, STATUS_OK, Trials
from hyperopt.pyll_utils import hp_quniform as quniform
from hyperopt.pyll_utils import hp_uniformint as uniformint

BDTfeats = ["phoR9Full5x5", "phoS4Full5x5", "phoEmaxOESCrFull5x5", "phoE2ndOESCrFull5x5", "phoE2ndOEmaxFull5x5", "phoE1x3OESCrFull5x5", "phoE2x5OESCrFull5x5", "phoE5x5OESCrFull5x5",
			"phoEmaxOE3x3Full5x5", "phoE2ndOE3x3Full5x5", "phoSigmaIEtaIEta", "phoSigmaIEtaIPhi", "phoSigmaIPhiIPhi", "phoSieieOSipipFull5x5", "phoEtaWidth", "phoPhiWidth", "phoEtaWOPhiWFull5x5"]
BDTfeats.sort()
print("\nBDT input features (" + str(len(BDTfeats)) + ") :")
print(BDTfeats)
print("\n\n")

pbar = ProgressBar()
for totalData in pbar(read_root(paths=args.inFilePath, key=args.inTreeName, columns=(BDTfeats + ['PtEtaRwBG', 'xSecW', 'isSignal', 'isTrain', 'isValidation']), chunksize=500000)):

	# set weights : reweight background pT-eta to signal; use normal xsec
	# weights for signal
	print("\nTotal data:")
	totalData['bdtWeight'] = totalData['xSecW'] * totalData['PtEtaRwBG']
	sampleStats(totalData)

	print("\nTraining data:")
	trainDF = totalData[(totalData.isTrain == 1) & (totalData.isValidation == 0)]
	sampleStats(trainDF)

	print("\nValidation data:")
	validationDF = totalData[(totalData.isTrain == 1) & (totalData.isValidation == 1)]
	sampleStats(validationDF)

	trainPho = xg.DMatrix(trainDF[BDTfeats].values, label=trainDF['isSignal'].values, weight=trainDF['bdtWeight'].values, feature_names=BDTfeats, nthread=-1)
	validationPho = xg.DMatrix(validationDF[BDTfeats].values, label=validationDF['isSignal'].values, weight=validationDF['bdtWeight'].values, feature_names=BDTfeats, nthread=-1)

	def score(params):
		print("\n\nTraining with params : " + str(params))

		phoModel = xg.train(params, dtrain=trainPho, evals=[(validationPho, "val")], early_stopping_rounds=10, verbose_eval=True, num_boost_round=10000)

		trainPredictions = phoModel.predict(trainPho, ntree_limit=phoModel.best_iteration + 1)
		trainScore = aucW(trainDF['isSignal'].values, trainPredictions, trainDF['bdtWeight'].values)

		theScore = (1. - phoModel.best_score) * abs(trainScore - phoModel.best_score)

		params['bestValScore'] = phoModel.best_score
		params['bestIteration'] = phoModel.best_iteration
		params['bestNtreeLimit'] = phoModel.best_ntree_limit
		params['trainScore'] = trainScore
		params['theScore'] = theScore

		print("\tbestValScore: " + str(params['bestValScore']))
		print("\tbestNtreeLimit: " + str(params['bestNtreeLimit']))
		print("\ttrainScore : " + str(params['trainScore']))
		print("\ttheScore : " + str(params['theScore']))

		scoreLogFile.write(str(params) + "\n")

		return {'loss': theScore, 'status': STATUS_OK}

	trials = Trials()

	### step 1: min_child_weight and max_depth
	# def optimizeMD_MCW(trials):
	# 	now = datetime.datetime.now()
	# 	print('\n\n' +now.strftime("%Y-%m-%d %H:%M:%S") + '\nLaunching tuning module... \n\n')

	# 	space = {

	# 		'objective': 'binary:logistic',
	# 		'eta': 0.1,
	# 		'eval_metric': 'auc',
	# 		'booster': 'gbtree',
	# 		'verbosity': 1,
	# 		'nthread': -1,
	# 		'gamma': 0,
	# 		'max_depth': uniformint('max_depth', 6, 25),
	# 		'min_child_weight': quniform('min_child_weight', 0., 10., 0.1),
	# 		'subsample': 0.8,
	# 		'colsample_bytree': 0.8,

	# 		'sampling_method': 'uniform',
	# 		'tree_method': 'exact',

	# 		'predictor': 'cpu_predictor'
	# 	}
	# 	best = fmin(score, space, algo=tpe.suggest, trials=trials, max_evals=100)
	# 	print("Optimum max_depth and min_child_weight:")
	# 	print(best)

	# optimizeMD_MCW(trials)

	# ### step 2: gamma
	# def optimizeGamma(trials):
	# 	now = datetime.datetime.now()
	# 	print('\n\n' +now.strftime("%Y-%m-%d %H:%M:%S") + '\nLaunching tuning module... \n\n')

	# 	space = {

	# 		'objective': 'binary:logistic',
	# 		'eta': 0.1,
	# 		'eval_metric': 'auc',
	# 		'booster': 'gbtree',
	# 		'verbosity': 1,
	# 		'nthread': -1,
	# 		'gamma': quniform('gamma', 0., 5., 0.05),
	# 		'max_depth': 22,
	# 		'min_child_weight': 1.8,
	# 		'subsample': 0.8,
	# 		'colsample_bytree': 0.8,

	# 		'sampling_method': 'uniform',
	# 		'tree_method': 'exact',

	# 		'predictor': 'cpu_predictor'
	# 	}
	# 	best = fmin(score, space, algo=tpe.suggest, trials=trials, max_evals=100)
	# 	print("Optimum max_depth and min_child_weight:")
	# 	print(best)

	# optimizeGamma(trials)

	# ### step 3: subsample and colsample_bytree
	# def optimizeSS_CST(trials):
	# 	now = datetime.datetime.now()
	# 	print('\n\n' +now.strftime("%Y-%m-%d %H:%M:%S") + '\nLaunching tuning module... \n\n')

	# 	space = {

	# 		'objective': 'binary:logistic',
	# 		'eta': 0.1,
	# 		'eval_metric': 'auc',
	# 		'booster': 'gbtree',
	# 		'verbosity': 1,
	# 		'nthread': -1,
	# 		'gamma': 0.05,
	# 		'max_depth': 22,
	# 		'min_child_weight': 1.8,
	# 		'subsample': quniform('subsample', 0.5, 1., 0.1),
	# 		'colsample_bytree': quniform('colsample_bytree', 0.5, 1., 0.1),

	# 		'sampling_method': 'uniform',
	# 		'tree_method': 'exact',

	# 		'predictor': 'cpu_predictor'
	# 	}
	# 	best = fmin(score, space, algo=tpe.suggest, trials=trials, max_evals=100)
	# 	print("Optimum max_depth and min_child_weight:")
	# 	print(best)

	# optimizeSS_CST(trials)


	### step 4: alpha
	# def optimizeAlpha(trials):
	# 	now = datetime.datetime.now()
	# 	print('\n\n' +now.strftime("%Y-%m-%d %H:%M:%S") + '\nLaunching tuning module... \n\n')

	# 	space = {

	# 		'objective': 'binary:logistic',
	# 		'eta': 0.1,
	# 		'eval_metric': 'auc',
	# 		'booster': 'gbtree',
	# 		'verbosity': 1,
	# 		'nthread': -1,
	# 		'gamma': 0.05,
	# 		'max_depth': 22,
	# 		'min_child_weight': 1.8,
	# 		'subsample': 1.0,
	# 		'colsample_bytree': 0.8,
	# 		'sampling_method': 'uniform',
	# 		'tree_method': 'exact',
	# 		# 'alpha': hp.choice('alpha', [0, 0.001, 0.005, 0.01, 0.05]),
	# 		'alpha': hp.choice('alpha', [0]),
	# 		'predictor': 'cpu_predictor'
	# 	}
	# 	best = fmin(score, space, algo=tpe.suggest, trials=trials, max_evals=1)
	# 	print("Optimum max_depth and min_child_weight:")
	# 	print(best)

	# optimizeAlpha(trials)

	# ### step 5: lambda
	# def optimizeAlpha(trials):
	# 	now = datetime.datetime.now()
	# 	print('\n\n' +now.strftime("%Y-%m-%d %H:%M:%S") + '\nLaunching tuning module... \n\n')

	# 	space = {

	# 		'objective': 'binary:logistic',
	# 		'eta': 0.1,
	# 		'eval_metric': 'auc',
	# 		'booster': 'gbtree',
	# 		'verbosity': 1,
	# 		'nthread': -1,
	# 		'gamma': 0.05,
	# 		'max_depth': 22,
	# 		'min_child_weight': 1.8,
	# 		'subsample': 1.0,
	# 		'colsample_bytree': 0.8,
	# 		'sampling_method': 'uniform',
	# 		'tree_method': 'exact',
	# 		'alpha': 0.01,
	# 		'lambda': hp.choice('lambda', [0, 0.001, 0.01, 0.1, 1., 5., 10.]),
	# 		'predictor': 'cpu_predictor'
	# 	}
	# 	best = fmin(score, space, algo=tpe.suggest, trials=trials, max_evals=100)
	# 	print("Optimum max_depth and min_child_weight:")
	# 	print(best)

	# optimizeAlpha(trials)


	### step 6: eta
	def optimizeEta(trials):
		now = datetime.datetime.now()
		print('\n\n' +now.strftime("%Y-%m-%d %H:%M:%S") + '\nLaunching tuning module... \n\n')

		space = {

			'objective': 'binary:logistic',
			'eta': hp.choice('eta', [0.01, 0.02, 0.05, 0.08, 0.1, 0.12]),
			'eval_metric': 'auc',
			'booster': 'gbtree',
			'verbosity': 1,
			'nthread': -1,
			'gamma': 0.05,
			'max_depth': 22,
			'min_child_weight': 1.8,
			'subsample': 1.0,
			'colsample_bytree': 0.8,
			'sampling_method': 'uniform',
			'tree_method': 'exact',
			# 'alpha': hp.choice('alpha', [0, 0.001, 0.005, 0.01, 0.05]),
			'alpha': 0.,
			'predictor': 'cpu_predictor'
		}

		best = fmin(score, space, algo=tpe.suggest, trials=trials, max_evals=10)
		print("Optimum parameters:")
		print(best)

	optimizeEta(trials)


	# def optimize(trials):
	# 	now = datetime.datetime.now()
	# 	print(now.strftime("%Y-%m-%d %H:%M:%S") + '\nLaunching tuning module...')

	# 	space = {

	# 		'objective': 'binary:logistic',
	# 		'eta': quniform('eta', 0.02, 0.1, 0.01),
	# 		'eval_metric': 'auc',

	# 		'booster': 'gbtree',
	# 		'verbosity': 1,
	# 		'nthread': -1,

	# 		'gamma': quniform('gamma', 0., 5., 0.5),
			
	# 		'max_depth': uniformint('max_depth', 6, 25),

	# 		'min_child_weight': quniform('min_child_weight', 0., 10., 0.1),

	# 		'max_delta_step': quniform('max_delta_step', 0, 5, 1),

	# 		'subsample': quniform('subsample', 0.1, 1., 0.1),
	# 		'colsample_bytree': quniform('colsample_bytree', 0.5, 1., 0.05),

	# 		'sampling_method': 'uniform',

	# 		# 'scale_pos_weight': BoS,

	# 		# 'lambda': quniform('lambda', 0., 1., 0.05),
	# 		# 'alpha': quniform('alpha', 0., 0.5, 0.05),
	# 		'tree_method': 'exact',

	# 		'predictor': 'cpu_predictor'
	# 	}

	break

scoreLogFile.close()


print("\n************************************************************************************************************************************************")
now = datetime.datetime.now()
print("Done @ " + now.strftime("%Y-%m-%d %H:%M:%S"))
print("************************************************************************************************************************************************")
