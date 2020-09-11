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
parser.add_argument('--saveDir', type=str, default='tuningmGPUbays/',
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
global logFileName
logFileName = args.saveDir + '/mGPUbaysTuneLog_' + \
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

scoreLogName = args.saveDir + '/mGPUbaysTuneScoreLog_' +now.strftime("%Y_%m_%d_%H_%M_%S") + ".log"
scoreLogFile = open(scoreLogName, "w")

print("featFilePath\t=\t" + args.featFilePath)
print("featTreeName\t=\t" + args.featTreeName)
print("indexFilePath\t=\t" + args.indexFilePath)
print("indexTreeName\t=\t" + args.indexTreeName)
print("saveDir\t\t=\t" + args.saveDir)
print("logFile\t\t=\t" + logFileName)

import pandas as pd
from custom_auc import aucW, sampleStats
import xgboost as xg
from hyperopt import fmin, tpe, hp, STATUS_OK, STATUS_FAIL, Trials, space_eval
from hyperopt.pyll_utils import hp_quniform as quniform
from hyperopt.pyll_utils import hp_uniformint as uniformint
import ROOT
import json
from dask_cuda import LocalCUDACluster
from dask.distributed import Client
import dask.dataframe as dd
from random import randint
import numpy as np

global BDTfeats
BDTfeats = ["phoR9Full5x5", "phoS4Full5x5", "phoEmaxOESCrFull5x5", "phoE2ndOESCrFull5x5", "phoE2ndOEmaxFull5x5", "pho2x2OE3x3Full5x5", "phoE1x3OESCrFull5x5", "phoE2x5OESCrFull5x5", "phoE5x5OESCrFull5x5",
			"phoEmaxOE3x3Full5x5", "phoE2ndOE3x3Full5x5", "phoSigmaIEtaIEta", "phoSigmaIEtaIPhi", "phoSigmaIPhiIPhi", "phoSieieOSipipFull5x5", "phoEtaWidth", "phoPhiWidth", "phoEtaWOPhiWFull5x5"]
BDTfeats.sort()
print("\nBDT input features (" + str(len(BDTfeats)) + ") :")
print(BDTfeats)
print("\n")

space = {
			'objective': 'binary:logistic',
			'eval_metric': 'auc',
			'booster': 'gbtree',
			'sampling_method': 'gradient_based',
			'tree_method': 'gpu_hist',
			'predictor': 'gpu_predictor',
			'verbosity': 1,
			'nthread': -1, 
			'seed': randint(10000,10000000),

			'scale_pos_weight': 0.96,

			'max_depth': 6,
			'min_child_weight': 6521,

			'lambda': 41.6,
			'gamma': 28.0,
			'alpha': 546,

			'subsample': 0.29, 
			'colsample_bytree': 0.66, 


			'eta': hp.quniform('learning_rate', 0.01, 0.08, 0.005)
			# 'eta': 0.3
}

print("Search space:")
print(str(space))

# Ntune = 1000000
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

print("\nValidation data:")
validationDF = totalData[(totalData.isValidation == 1)]
validationDF = validationDF.sample(frac = 1)
sampleStats(validationDF)

del totalData

def main(client):
	global trainDF
	global validationDF
	global ddf_trainx
	global ddf_trainy
	global ddf_trainw
	ddf_trainx = dd.from_pandas(trainDF[BDTfeats], npartitions=8)
	ddf_trainy = dd.from_pandas(trainDF['isSignal'], npartitions=8)
	ddf_trainw = dd.from_pandas(trainDF['flatPtEtaRwNoXsec'], npartitions=8)
	print(ddf_trainx)
	del trainDF
	global trainPho
	trainPho = xg.dask.DaskDMatrix(client, ddf_trainx , label=ddf_trainy, weight=ddf_trainw, feature_names=BDTfeats)

	# ddf_trainx = dd.from_array(trainDF[BDTfeats].values)
	# ddf_trainy = dd.from_array(trainDF['isSignal'].values)
	# ddf_trainw = dd.from_array(trainDF['flatPtEtaRwNoXsec'].values)
	# trainPho = xg.dask.DaskDeviceQuantileDMatrix(client, ddf_trainx, label=ddf_trainy, weight=ddf_trainw, feature_names=BDTfeats, max_bin=512)

	global ddf_validationx
	global ddf_validationy
	global ddf_validationw
	ddf_validationx = dd.from_pandas(validationDF[BDTfeats], npartitions=8)
	ddf_validationy = dd.from_pandas(validationDF['isSignal'], npartitions=8)
	ddf_validationw = dd.from_pandas(validationDF['flatPtEtaRwNoXsec'], npartitions=8)
	print(ddf_validationx)
	del validationDF
	global validationPho
	validationPho = xg.dask.DaskDMatrix(client, ddf_validationx , label=ddf_validationy, weight=ddf_validationw, feature_names=BDTfeats)


def score(params):
	if len(trials.trials)>1:
		for x in trials.trials[:-1]:
			space_point_index = dict([(key,value[0]) for key,value in x['misc']['vals'].items() if len(value)>0])
			global space
			if params == space_eval(space,space_point_index):
				return {'loss': 100000000., 'status': STATUS_FAIL}

	scoreLogFile.write("\n\n\n\n")
	print("\n"+now.strftime("%Y_%m_%d_%H_%M_%S"))
	print("\nTraining with params : " + str(params))
	
	output = xg.dask.train(client, params, dtrain=trainPho, evals=[(trainPho, 'train'), (validationPho, "validation")], early_stopping_rounds=20, verbose_eval=True, num_boost_round=100000)
	accuracyProgress = output['history']
	phoModel = output['booster']

	theScore = (1. - phoModel.best_score)
	
	print("\tbestValScore: " + str(phoModel.best_score))
	print("\tbestNtreeLimit: " + str(phoModel.best_ntree_limit))
	print("\ttheScore : " + str(theScore))

	print(json.dumps(accuracyProgress, indent=4, sort_keys=True))
	scoreLogFile.write(str(params) + "\n")
	scoreLogFile.write(json.dumps(accuracyProgress, indent=4, sort_keys=True))

	return {'loss': theScore, 'status': STATUS_OK}

def optimize(nTrials=1000):
	global trials
	trials = Trials()
	now = datetime.datetime.now()
	print('\n\n' +now.strftime("%Y-%m-%d %H:%M:%S") + '\nLaunching tuning module... \n\n')
	best = fmin(score, space, algo=tpe.suggest, trials=trials, max_evals=nTrials)
	print("Optimum: "+str(best))
	scoreLogFile.write("Optimum:\n")
	scoreLogFile.write(str(best))

global cluster
global client
cluster=LocalCUDACluster()
client=Client(cluster)
main(client)

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
# 			'seed': randint(10000,10000000),

# 			'scale_pos_weight': 0.96,

# 			'max_depth': 6,
# 			'min_child_weight': hp.quniform('min_child_weight', 6100, 7000, 100),

# 			# 'alpha': hp.quniform('alpha', 300, 1000,1),
# 			# 'lambda': hp.quniform('lambda', 0.1, 50., 0.01),
# 			# 'gamma': hp.quniform('gamma', 5, 50, 0.01),

# 			# 'subsample': hp.quniform('subsample', 0.2, 0.5, 0.01),
# 			# 'colsample_bytree': hp.choice('colsample_bytree',[0.55, 0.56, 0.61, 0.66, 0.67, 0.72]),

# 			'alpha': 588,
# 			'lambda': 30.13,
# 			'gamma': 29.31,

# 			'subsample': 0.29,
# 			'colsample_bytree': 0.67,


# 			# 'eta': hp.qloguniform('learning_rate', np.log(0.01), np.log(0.08), 0.001)
# 			'eta': 0.3
# }

# Top 11
# 0.003126999999999991, 415
# {'alpha': 588.0, 'booster': 'gbtree', 'colsample_bytree': 0.67, 'eta': 0.3, 'eval_metric': 'auc', 'gamma': 29.310000000000002, 'lambda': 30.13, 'max_depth': 6, 'min_child_weight': 6998.0, 'nthread': -1, 'objective': 'binary:logistic', 'predictor': 'gpu_predictor', 'sampling_method': 'gradient_based', 'scale_pos_weight': 0.96, 'seed': 2299504, 'subsample': 0.29, 'tree_method': 'gpu_hist', 'verbosity': 1}

# 0.003129999999999966, 235
# {'alpha': 470.0, 'booster': 'gbtree', 'colsample_bytree': 0.56, 'eta': 0.3, 'eval_metric': 'auc', 'gamma': 39.08, 'lambda': 26.38, 'max_depth': 6, 'min_child_weight': 6193.0, 'nthread': -1, 'objective': 'binary:logistic', 'predictor': 'gpu_predictor', 'sampling_method': 'gradient_based', 'scale_pos_weight': 0.96, 'seed': 2299504, 'subsample': 0.31, 'tree_method': 'gpu_hist', 'verbosity': 1}

# 0.0031360000000000277, 399
# {'alpha': 801.0, 'booster': 'gbtree', 'colsample_bytree': 0.67, 'eta': 0.3, 'eval_metric': 'auc', 'gamma': 36.050000000000004, 'lambda': 32.57, 'max_depth': 6, 'min_child_weight': 7655.0, 'nthread': -1, 'objective': 'binary:logistic', 'predictor': 'gpu_predictor', 'sampling_method': 'gradient_based', 'scale_pos_weight': 0.96, 'seed': 2299504, 'subsample': 0.38, 'tree_method': 'gpu_hist', 'verbosity': 1}

# 0.0031480000000000397, 448
# {'alpha': 688.0, 'booster': 'gbtree', 'colsample_bytree': 0.67, 'eta': 0.3, 'eval_metric': 'auc', 'gamma': 31.34, 'lambda': 38.99, 'max_depth': 6, 'min_child_weight': 6070.0, 'nthread': -1, 'objective': 'binary:logistic', 'predictor': 'gpu_predictor', 'sampling_method': 'gradient_based', 'scale_pos_weight': 0.96, 'seed': 2299504, 'subsample': 0.3, 'tree_method': 'gpu_hist', 'verbosity': 1}

# 0.003149999999999986, 328
# {'alpha': 758.0, 'booster': 'gbtree', 'colsample_bytree': 0.67, 'eta': 0.3, 'eval_metric': 'auc', 'gamma': 34.02, 'lambda': 32.34, 'max_depth': 6, 'min_child_weight': 7946.0, 'nthread': -1, 'objective': 'binary:logistic', 'predictor': 'gpu_predictor', 'sampling_method': 'gradient_based', 'scale_pos_weight': 0.96, 'seed': 2299504, 'subsample': 0.35000000000000003, 'tree_method': 'gpu_hist', 'verbosity': 1}

# 0.003149999999999986, 354
# {'alpha': 506.0, 'booster': 'gbtree', 'colsample_bytree': 0.56, 'eta': 0.3, 'eval_metric': 'auc', 'gamma': 26.38, 'lambda': 45.61, 'max_depth': 6, 'min_child_weight': 7991.0, 'nthread': -1, 'objective': 'binary:logistic', 'predictor': 'gpu_predictor', 'sampling_method': 'gradient_based', 'scale_pos_weight': 0.96, 'seed': 2299504, 'subsample': 0.31, 'tree_method': 'gpu_hist', 'verbosity': 1}

# 0.0031400000000000317, 297
# {'alpha': 636.0, 'booster': 'gbtree', 'colsample_bytree': 0.56, 'eta': 0.3, 'eval_metric': 'auc', 'gamma': 12.92, 'lambda': 36.97, 'max_depth': 6, 'min_child_weight': 7668.0, 'nthread': -1, 'objective': 'binary:logistic', 'predictor': 'gpu_predictor', 'sampling_method': 'gradient_based', 'scale_pos_weight': 0.96, 'seed': 2299504, 'subsample': 0.25, 'tree_method': 'gpu_hist', 'verbosity': 1}

# 0.003147000000000011, 400
# {'alpha': 713.0, 'booster': 'gbtree', 'colsample_bytree': 0.67, 'eta': 0.3, 'eval_metric': 'auc', 'gamma': 40.88, 'lambda': 28.85, 'max_depth': 6, 'min_child_weight': 7291.0, 'nthread': -1, 'objective': 'binary:logistic', 'predictor': 'gpu_predictor', 'sampling_method': 'gradient_based', 'scale_pos_weight': 0.96, 'seed': 2299504, 'subsample': 0.37, 'tree_method': 'gpu_hist', 'verbosity': 1}

# 0.003141999999999978, 444
# {'alpha': 750.0, 'booster': 'gbtree', 'colsample_bytree': 0.67, 'eta': 0.3, 'eval_metric': 'auc', 'gamma': 42.38, 'lambda': 44.12, 'max_depth': 6, 'min_child_weight': 7428.0, 'nthread': -1, 'objective': 'binary:logistic', 'predictor': 'gpu_predictor', 'sampling_method': 'gradient_based', 'scale_pos_weight': 0.96, 'seed': 2299504, 'subsample': 0.29, 'tree_method': 'gpu_hist', 'verbosity': 1}

# 0.003147000000000011, 387
# {'alpha': 411.0, 'booster': 'gbtree', 'colsample_bytree': 0.56, 'eta': 0.3, 'eval_metric': 'auc', 'gamma': 26.810000000000002, 'lambda': 45.0, 'max_depth': 6, 'min_child_weight': 6966.0, 'nthread': -1, 'objective': 'binary:logistic', 'predictor': 'gpu_predictor', 'sampling_method': 'gradient_based', 'scale_pos_weight': 0.96, 'seed': 2299504, 'subsample': 0.32, 'tree_method': 'gpu_hist', 'verbosity': 1}
 
#  0.0031409999999999494, 368
#  {'alpha': 676.0, 'booster': 'gbtree', 'colsample_bytree': 0.67, 'eta': 0.3, 'eval_metric': 'auc', 'gamma': 29.66, 'lambda': 33.24, 'max_depth': 6, 'min_child_weight': 7202.0, 'nthread': -1, 'objective': 'binary:logistic', 'predictor': 'gpu_predictor', 'sampling_method': 'gradient_based', 'scale_pos_weight': 0.96, 'seed': 2299504, 'subsample': 0.32, 'tree_method': 'gpu_hist', 'verbosity': 1}
 
# Optimum: {'min_child_weight': 6500.0}
# >>> space['min_child_weight']=quniform('min_child_weight', 6460, 6540, 10)
# Optimum: {'min_child_weight': 6520.0}
# >>> space['min_child_weight']=quniform('min_child_weight', 6516, 6524, 1)
# Optimum: {'min_child_weight': 6521.0}

# >>> space['alpha'] = hp.quniform('alpha', 500, 800,1)
# >>> space['lambda'] = hp.quniform('lambda', 30., 46., 0.1)
# >>> space['gamma'] = hp.quniform('gamma', 26, 40, 0.1)
# >>> space['subsample'] = hp.quniform('subsample', 0.25, 0.35, 0.01)
# >>> space['colsample_bytree'] = hp.choice('colsample_bytree',[0.55, 0.56, 0.61, 0.66, 0.67, 0.72])
# >>> optimize(500)
# Optimum: {'alpha': 547.0, 'colsample_bytree': 66, 'gamma': 27.900000000000002, 'lambda': 41.6, 'subsample': 0.29}
# >>> space['alpha'] = hp.quniform('alpha', 545, 550,1)
# >>> space['colsample_bytree'] = hp.choice('colsample_bytree',[0.55, 0.61, 0.66, 0.72])
# >>> space['gamma'] = hp.quniform('gamma', 25, 32, 0.02)
# >>> space['lambda'] = hp.quniform('lambda', 38., 44., 0.05)
# >>> space['subsample'] = hp.quniform('subsample', 0.28, 0.30, 0.01)
# optimize(200)
# Optimum: {'alpha': 548.0, 'colsample_bytree': 0.55, 'gamma': 26.42, 'lambda': 42.300000000000004, 'subsample': 0.29}
# >>> space['alpha'] = hp.quniform('alpha', 545, 550,1)
# >>> space['colsample_bytree'] = 0.66
# >>> space['gamma'] = hp.quniform('gamma', 27.2, 28.4, 0.2)
# >>> space['lambda'] = hp.quniform('lambda', 41., 43., 0.2)
# >>> space['subsample'] = 0.29
# optimize(1000)
# Optimum: {'alpha': 547.0, 'gamma': 27.8, 'lambda': 41.400000000000006}
# >>> space['lambda'] = hp.quniform('lambda', 41.2, 41.6., 0.05)
# >>> space['gamma'] = hp.quniform('gamma', 27.6, 27.8, 0.05)
# >>> space['alpha'] = hp.quniform('alpha', 546, 548, 0.5)
# Optimum: {'alpha': 548.0, 'gamma': 27.8, 'lambda': 41.25}

# >>> space['lambda'] = hp.quniform('lambda', 41.1, 41.4, 0.1)
# >>> space['gamma'] = hp.quniform('gamma', 27.6, 27.9, 0.1)
# >>> space['alpha'] = hp.quniform('alpha', 547, 549, 0.5)
# Optimum: {'alpha': 547.0, 'gamma': 27.900000000000002, 'lambda': 41.400000000000006}
# >>> space['alpha'] = hp.quniform('alpha', 545, 548, 1)
# >>> space['gamma'] = hp.quniform('gamma', 27.8, 28.2, 0.2)
# >>> space['lambda'] = hp.quniform('lambda', 41.2, 41.8, 0.2)
# Optimum: {'alpha': 546.0, 'gamma': 28.0, 'lambda': 41.6}
# >>> space['lambda'] = 41.6
# >>> space['gamma'] = 28.0
# >>> space['alpha'] = 546
# >>> space['eta']= hp.quniform('learning_rate', 0.01, 0.08, 0.005)
# Optimum: {'learning_rate': 0.045}
# >>> space['eta']= hp.quniform('learning_rate', 0.042, 0.048, 0.001)
# Optimum: {'learning_rate': 0.048}
# >>> space['eta']= hp.quniform('learning_rate', 0.046, 0.049, 0.0002)
# Optimum: {'learning_rate': 0.046400000000000004}
# >>> space['eta']= hp.quniform('learning_rate', 0.0462, 0.0466, 0.0001)
# Optimum: {'learning_rate': 0.0466}
# >>> space['eta']= hp.quniform('learning_rate', 0.0465, 0.0468, 0.0001)
# Optimum: {'learning_rate': 0.0466}
# >>> space['eta']= hp.quniform('learning_rate', 0.0465, 0.0466, 0.0001)
# Optimum: {'learning_rate': 0.0466}