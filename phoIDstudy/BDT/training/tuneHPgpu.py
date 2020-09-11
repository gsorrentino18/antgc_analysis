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
parser.add_argument('--saveDir', type=str, default='tuning/',
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
logFileName = args.saveDir + '/tuneLog_' + \
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

			'max_depth': hp.choice('max_depth', [6,8,10]),
			'min_child_weight': quniform('min_child_weight', 0,40,10),

			'alpha' : 0.,
			'lambda' : 0.,
			'gamma': 0.,
						
			'subsample': 0.6,
			'colsample_bytree':	0.8,

			'eta': 0.05,

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

print("\nValidation data:")
validationDF = totalData[(totalData.isValidation == 1)]
validationDF = validationDF.sample(frac = 1)
sampleStats(validationDF)

trainPho = xg.DMatrix(trainDF[BDTfeats].values, label=trainDF['isSignal'].values, weight=trainDF['flatPtEtaRwNoXsec'].values, feature_names=BDTfeats, nthread=-1)
validationPho = xg.DMatrix(validationDF[BDTfeats].values, label=validationDF['isSignal'].values, weight=validationDF['flatPtEtaRwNoXsec'].values, feature_names=BDTfeats, nthread=-1)

del totalData
del trainDF
del validationDF

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
	accuracyProgress = dict()
	phoModel = xg.train(params, dtrain=trainPho, evals=[(trainPho, "train"), (validationPho, "val")], evals_result=accuracyProgress, early_stopping_rounds=20, verbose_eval=True, num_boost_round=10000)
	theScore = (1. - phoModel.best_score)

	scoreLogFile.write("\n\n\n\n")

	print("\tbestValScore: " + str(phoModel.best_score))
	print("\tbestNtreeLimit: " + str(phoModel.best_ntree_limit))
	print("\ttheScore : " + str(theScore))

	scoreLogFile.write(str(params) + "\n")
	scoreLogFile.write(json.dumps(accuracyProgress, indent=4, sort_keys=True))

	return {'loss': theScore, 'status': STATUS_OK}

def optimize():
	global trials
	trials = Trials()
	now = datetime.datetime.now()
	print('\n\n' +now.strftime("%Y-%m-%d %H:%M:%S") + '\nLaunching tuning module... \n\n')
	best = fmin(score, space, algo=tpe.suggest, trials=trials, max_evals=5000)
	print("Optimum max_depth and min_child_weight:")
	print(best)

	scoreLogFile.write("Optimum max_depth and min_child_weight:\n")
	scoreLogFile.write(str(best))

optimize()

























































# space = {
# 			'objective': 'binary:logistic',
# 			'eta': 0.3,
# 			'eval_metric': 'auc',
# 			'booster': 'gbtree',
# 			'verbosity': 1,
# 			'nthread': -1,
# 			'gamma': 0.,
# 			'max_depth': 20,
# 			'min_child_weight': 0.,
# 			'subsample': 1.,
# 			'colsample_bytree': 1.,

# 			# 'alpha' : 1.1,
# 			# 'lambda' : 0.5,
# 			# 'max_bin': 512,
# 			'sampling_method': 'gradient_based',
# 			'tree_method': 'gpu_hist',
# 			'predictor': 'gpu_predictor'
# }        
# 			bestValScore: 0.99801                           
# 			bestNtreeLimit: 63                              
# 			theScore : 0.0019900000000000473


# 1a 	|	'max_depth': uniformint('max_depth', 18,25)
# 		|	'min_child_weight': quniform('min_child_weight', 16., 30., 2.)                                               
# 		best loss: 0.0019019999999999593] 
# 		{'max_depth': 23.0, 'min_child_weight': 20.0}
#
#
# 1b 	|	'max_depth': uniformint('max_depth', 21,25),
# 		|	'min_child_weight': quniform('min_child_weight', 18., 22., 0.5),                                                
# 		best loss: 0.0019010000000000415] 
# 		{'max_depth': 22.0, 'min_child_weight': 21.0}
#
#
# 1c	|	space['min_child_weight'] = quniform('min_child_weight', 20., 22., 0.2),
# 		|	space['max_depth'] = uniformint('max_depth', 21,24)                                                 
# 		best loss: 0.0019010000000000415] 
# 		{'max_depth': 22.0, 'min_child_weight': 21.0}
#
#
# 2a 	|	'gamma': quniform('gamma', 0.1., 1., 0.2)
#		best loss: 0.0019169999999999465] 
# 		{'gamma': 0.8}
#
# 2b	|	'gamma': quniform('gamma', 0.7, 0.9, 0.04)
#		best loss: 0.0019160000000000288] 
# 			{'gamma': 0.88}
#
# 2c	|	'gamma': quniform('gamma', 0.82, 0.96, 0.01)
#		best loss: 0.0019120000000000248] 
#			{'gamma': 0.9400000000000001}
#
# 2d	|	'gamma': quniform('gamma', 0., 1.5, 0.1)
#		best loss: 0.0019010000000000415] 
#			{'gamma': 0.0}
#
# 2e	|	'gamma'=quniform('gamma', 0., 0.2, 0.02)
#		best loss: 0.0019010000000000415] 
#			{'gamma': 0.0}
#
# 3a	|	'subsample'=quniform('subsample', 0., 1., 0.2)
#		|	'colsample_bytree'=quniform('colsample_bytree', 0., 1., 0.2)
#		0.0018920000000000048]
# 			{'colsample_bytree': 0.8, 'subsample': 0.6000000000000001}
#
# 3b	|	'colsample_bytree'=quniform('colsample_bytree', 0.7, 0.9, 0.05)			-> different terminal
#		|	'subsample']=quniform('subsample', 0.5, 0.7, 0.05)
#			best loss: 0.0018930000000000335] 
#			{'colsample_bytree': 0.75, 'subsample': 0.5}
#
# 4a	|	'alpha'=hp.choice('alpha', [0,0.001,0.1,1,10])
#			best loss: 0.0018519999999999648] 
#			{'alpha': 10}
# 
# 4b	|	'alpha'=hp.choice('alpha', [10,20,50,100])
# 			best loss: 0.0018200000000000438] 
#			{'alpha': 100}

# 4c	|	'alpha'=hp.choice('alpha', [100,500,1000])		
# 			best loss: 0.0018089999999999495
#			{'alpha': 1000}
#
# 4d	|	'alpha'=hp.choice('alpha', [2000,10000,100000])
#			best loss: 0.0018089999999999495] 
#			{'alpha': 2000}
#
# 4e	|	'alpha'=hp.choice('alpha', [2000,5000,8000])
#			best loss:	0.0018179999999999863
#			{'alpha': 2000}
#
# 4d	|	'alpha'=hp.choice('alpha', [1200,1500,1800])
#			best loss:	0.0018000000000000238
#			{'alpha': 1200}
#
# 4d	|	'alpha'=hp.choice('alpha', [800,900,1000,1100,1200])
#			best loss:	0.0017970000000000486	
#			{'alpha': 800}******
#
# 4e	|	'alpha'=hp.choice('alpha', [700,750,800,850,900])
#			best loss:	0.0018000000000000238
#			{'alpha': 700}
#
# 5a	|	'eta'=quniform('eta',0.05,0.15,0.05)
#			best loss:	0.001735000000000042
#			{'eta': 0.05}
#
# 5b	|	'eta'=quniform('eta',0.02,0.04,0.02)
#			best loss:	0.0017329999999999846
#			{'eta': 0.04}
#
# 5c	|	'eta'=quniform('eta',0.03,0.04,0.005)
#			best loss:	0.0017329999999999846
#			{'eta': 0.04}
#
# 5d	|	'eta'=quniform('eta',0.03,0.04,0.005)
#			best loss:	0.0017369999999999886
#			{'eta': 0.04}
#
# 5e	|	'eta'=quniform('eta',0.039,0.043,0.003)
#			best loss:	0.0017329999999999846
#			{'eta': 0.042}
#
# 6a	|	'lambda'=hp.choice('lambda',[1,10,100,1000])
#		best loss: 0.0017319999999999558
#			{'lambda': 1}

# 6b	|	'lambda']=hp.choice('lambda',[0.1,0.01,0.001,0.0001,0.00001]
# 		best loss: 0.0017319999999999558
# 			{'lambda': 0}




# print("\n************************************************************************************************************************************************")
# now = datetime.datetime.now()
# print("Done @ " + now.strftime("%Y-%m-%d %H:%M:%S"))
# print("************************************************************************************************************************************************")



# 
# space = {
# 			'objective': 'binary:logistic',
# 			'eta': 0.08,
# 			'eval_metric': 'auc',
# 			'booster': 'gbtree',
# 			'verbosity': 1,
# 			'nthread': -1,
# 			'gamma': 0.88,
# 			'max_depth': uniformint('max_depth', 18,25),
# 			'min_child_weight': quniform('min_child_weight', 16., 18., 0.5),
# 			'subsample': 0.84,
# 			'colsample_bytree': 0.8,
# 			'alpha' : 0.6,
# 			'lambda' : 1.08,
# 			'max_bin': 512,
# 			'sampling_method': 'gradient_based',
# 			'grow_policy': 'lossguide',
# 			'tree_method': 'gpu_hist',
# 			'predictor': 'gpu_predictor'
# }
# 1a 			'max_depth': uniformint('max_depth', 18,25),
# 				'min_child_weight': quniform('min_child_weight', 16., 18., 0.5), 											-> 	{'max_depth': 22.0, 'min_child_weight': 17.5} 	0.0020459999999999923

# 1b			'max_depth': uniformint('max_depth', 21,23),
# 				'min_child_weight': quniform('min_child_weight', 17.3, 17.7, 0.1),											->	{'max_depth': 21.0, 'min_child_weight': 17.400000000000002} 0.002051000000000025

# 1c	z3		'max_depth': uniformint('max_depth', 20,23),
# 				'min_child_weight': 17.5,																					->	{'max_depth': 20.0}	0.0020930000000000115

# 1d	z3		'max_depth': uniformint('max_depth', 18,22),
# 				'min_child_weight': 17.5,																					->	{'max_depth': 21.0}	0.0020809999999999995

# 1d	z2		'max_depth': 21.0,
# 				'min_child_weight': quniform('min_child_weight', 17.3, 17.6, 0.05),											->	{'min_child_weight': 17.3} 0.0020599999999999508

#		z2		'max_depth': 21,
#				'min_child_weight': quniform('min_child_weight', 17., 17.6, 0.1),											->	{'min_child_weight': 17.5} 0.0020369999999999555

#		z2		'gamma': quniform('gamma', 0.7, 0.9,0.02),
#				'max_depth': 21,
#				'min_child_weight': 17.3,																					->	{'gamma': 0.74}	0.0020480000000000498

#		z2		'gamma': quniform('gamma', 0.72, 0.76,0.01),																->	{'gamma': 0.76} 0.0020499999999999963

#		z2		'gamma': quniform('gamma', 0.74, 0.80,0.01),																->	{'gamma': 0.78}	0.0020850000000000035

#		z2		'subsample': quniform('subsample', 0.1, 0.9,0.2),
#				'colsample_bytree': quniform('colsample_bytree', 0.1, 0.9,0.2),												->	{'colsample_bytree': 0.8, 'subsample': 0.4}	0.0020809999999999995

#		z2		'subsample': quniform('subsample', 0.5, 0.9,0.1),
#				'colsample_bytree': quniform('colsample_bytree', 0.3, 0.5,0.1),												->	{'colsample_bytree': 0.5, 'subsample': 0.5}	0.0020379999999999843

#		z2		'alpha' : quniform('alpha', 0., 1.,0.5),
#				'lambda' : quniform('lambda', 0., 1.5,0.5),																	->	{'alpha': 1.0, 'lambda': 0.5} 0.0020559999999999468

#		z2		'alpha' : quniform('alpha', 0.8, 1.2,0.1),
#				'lambda' : 0.5																								->	{'alpha': 1.1}	0.002047000000000021

#		z2		'eta': quniform('eta', 0.06, 0.10,0.01),																	->	{'eta': 0.09} 	0.0020540000000000003



# 2M events
# theScore = (1. - phoModel.best_score) * abs(trainScore - phoModel.best_score)
# 1a
# space = {
# 			'objective': 'binary:logistic',
# 			'eta': 0.1,
# 			'eval_metric': 'auc',
# 			'booster': 'gbtree',
# 			'verbosity': 1,
# 			'nthread': -1,
# 			'gamma': 0.0,
# 			'max_depth': 22,
# 			'min_child_weight': quniform('min_child_weight', 3., 10., 0.5),
# 			'subsample': 0.8,
# 			'colsample_bytree': 0.8,
# 			'alpha' : 0.0,
# 			'max_bin': 512,
# 			# 'num_parallel_tree': 3,
# 			'sampling_method': 'gradient_based',
# 			'grow_policy': 'lossguide',
# 			'tree_method': 'gpu_hist',
# 			'updater': 'grow_gpu_hist',
# 			'predictor': 'gpu_predictor',
# 			'n_gpus':1
# }																			-> 	9.5			3.950519024722641e-06

# 1b 	'min_child_weight': quniform('min_child_weight', 8.5, 10.5, 0.1) 	->	10.5		3.7946155815837854e-06
# 1c 	'min_child_weight': quniform('min_child_weight', 8., 20., 1.)		-> 	15.0		3.7602084610908994e-06
# 1d 	'min_child_weight': quniform('min_child_weight', 14., 16., 0.1)		->	15.7		3.7427294516689265e-06

# 2 	uniformint('max_depth', 18, 25) 									->	20			3.830058820493023e-06 -> Forgot to update max depth to 20. Used 22.

# 3a	'gamma': quniform('gamma', 0., 5, 0.5)								->	2.			3.7299658338633463e-06
# 3b	'gamma': quniform('gamma', 0., 3., 0.2)								->	0.8			3.7731156337502075e-06
# 3c	'gamma': quniform('gamma', 0.6, 1.2, 0.05)							->	1.0			3.963717688264565e-06
# 3d	'gamma': quniform('gamma', 0.76, 1.1, 0.02)							->	0.88		3.8073339538602274e-06

# 4a	'subsample': quniform('subsample', 0.1, 1., 0.1)					->	0.9			3.7468373262150536e-06
# 4b	'subsample': quniform('subsample', 0.7, 1., 0.02)					->	0.84		3.7922586640565815e-06

# 5a	'colsample_bytree': quniform('colsample_bytree', 0.1, 1., 0.1)		->	0.7			3.911737141024762e-06
# 5b	'colsample_bytree': quniform('colsample_bytree', 0.6, 0.8, 0.02)	->	0.8			3.985835999887746e-06

# 6a	'alpha' : quniform('alpha', 0., 1., 0.05)							->	0.6			3.7121822608154757e-06
# 6b	'alpha' : quniform('alpha', 0.4, 0.8, 0.02)							->	0.6			3.630041046576212e-06

# 7a	'lambda' : quniform('lambda', 0., 2., 0.1)							->	1.0			3.910859470950211e-06
# 7b	'lambda' : quniform('lambda', 0.8, 1.2, 0.02)						->	1.08		5.997737214074227e-06

# 8a	'eta'	: quniform('eta', 0.0, 0.3, 0.1)							-> 	0.1			3.94461577508743e-06
# 8b	'eta'	: quniform('eta', 0.0, 0.2, 0.02)							->	0.08		3.811139209552489e-06












# 2M events
# 	theScore = (1. - phoModel.best_score)
# 1a
# space = {
# 			'objective': 'binary:logistic',
# 			'eta': 0.08,
# 			'eval_metric': 'auc',
# 			'booster': 'gbtree',
# 			'verbosity': 1,
# 			'nthread': -1,
# 			'gamma': 0.88,
# 			'max_depth': 22,
# 			'min_child_weight': quniform('min_child_weight', 14., 17., 0.1),
# 			'subsample': 0.84,
# 			'colsample_bytree': 0.8,
# 			'alpha' : 0.6,
# 			'lambda' : 1.08,
# 			'max_bin': 512,
# 			'sampling_method': 'gradient_based',
# 			'grow_policy': 'lossguide',
# 			'tree_method': 'gpu_hist',
# 			'predictor': 'gpu_predictor'
# }																			->	16.0		0.001931000000000016


# 3M events
# 	theScore = (1. - phoModel.best_score)
# 1a
# space = {
# 			'objective': 'binary:logistic',
# 			'eta': 0.08,
# 			'eval_metric': 'auc',
# 			'booster': 'gbtree',
# 			'verbosity': 1,
# 			'nthread': -1,
# 			'gamma': 0.88,
# 			'max_depth': 20,
# 			'min_child_weight': quniform('min_child_weight', 15., 17., 0.05),
# 			'subsample': 0.84,
# 			'colsample_bytree': 0.8,
# 			'alpha' : 0.6,
# 			'lambda' : 1.08,
# 			'max_bin': 512,
# 			'sampling_method': 'gradient_based',
# 			'grow_policy': 'lossguide',
# 			'tree_method': 'gpu_hist',
# 			'predictor': 'gpu_predictor'
# }
# 1a		quniform('min_child_weight', 15., 17., 0.05)					->	17.0		0.002082999999999946
# 1b		'min_child_weight': quniform('min_child_weight', 16., 18., 0.5) ->	


