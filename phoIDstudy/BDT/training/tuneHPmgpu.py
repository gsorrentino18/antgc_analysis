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
parser.add_argument('--saveDir', type=str, default='tuningMGPU/',
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
logFileName = args.saveDir + '/mGPUTuneLog_' + \
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

scoreLogName = args.saveDir + '/mGPUTuneScoreLog_' +now.strftime("%Y_%m_%d_%H_%M_%S") + ".log"
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

			'max_depth': 8, 
			'min_child_weight': 9177, 

			'alpha': quniform('alpha',370,385,5),
			'lambda': hp.choice('lambda',[0.615,0.62,0.625]),
			'gamma': hp.choice('gamma',[6.4,6.5,6.6]),

			'subsample': 0.35, 
			'colsample_bytree': 0.61, 

			'eta': 0.2			
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

def optimize():
	global trials
	trials = Trials()
	now = datetime.datetime.now()
	print('\n\n' +now.strftime("%Y-%m-%d %H:%M:%S") + '\nLaunching tuning module... \n\n')
	best = fmin(score, space, algo=tpe.suggest, trials=trials, max_evals=1000)
	print("Optimum: "+str(best))
	scoreLogFile.write("Optimum:\n")
	scoreLogFile.write(str(best))

global cluster
global client
cluster=LocalCUDACluster()
client=Client(cluster)
main(client)

# optimize()

# if __name__ == '__main__':
# 	global cluster
# 	global client
# 	cluster=LocalCUDACluster()
# 	client=Client(cluster)
# 	# with LocalCUDACluster() as cluster:
# 	# 	with Client(cluster) as client:
# 	main(client)


# Total data:
# Total Events    =       61970811 (B/S= 51775921 / 10194890)
# BackgroundSumW/SignalSumW       =       269189226389.518951 / 53004456496.799980        =       5.078615

# Training data:
# Total Events    =       24782828 (B/S= 20704871 / 4077957)
# BackgroundSumW/SignalSumW       =       107648148835.787064 / 21429305093.572651        =       5.023408

# Validation data:
# Total Events    =       18597782 (B/S= 15539211 / 3058571)
# BackgroundSumW/SignalSumW       =       80802652594.335175 / 15699627885.225983 =       5.146788

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

# 			'scale_pos_weight': quniform('scale_pos_weight', 0.6,1.,0.05),

# 			'max_depth': 10, 
# 			'min_child_weight': 9303, 

# 			'alpha': 406, 
# 			'lambda': 0.66, 
# 			'gamma': 0.28, 

# 			'subsample': 0.35, 
# 			'colsample_bytree': 0.78, 

# 			'eta': 0.2			
# }
# Optimum: {'scale_pos_weight': 0.9500000000000001}
# space['scale_pos_weight']=quniform('scale_pos_weight', 0.95,1.0,0.01)
# Optimum: {'scale_pos_weight': 0.96}
# >>> space['scale_pos_weight'] = 0.96
# >>> space['max_depth'] = uniformint('max_depth',6,8)
# >>> space['min_child_weight']= quniform('min_child_weight', 7000,10000,1000)
# Optimum: {'max_depth': 8.0, 'min_child_weight': 10000.0}
# >>> space['max_depth'] = 8
# >>> space['min_child_weight']= quniform('min_child_weight', 10000,12000,1000)
# Optimum: {'min_child_weight': 10000.0}
# >>> space['min_child_weight']= quniform('min_child_weight', 9000,9500,100)
# Optimum: {'min_child_weight': 9100.0}
# >>> space['min_child_weight']= quniform('min_child_weight', 9050,9150,10)
# Optimum: {'min_child_weight': 9150.0}
# >>> space['min_child_weight']= quniform('min_child_weight', 9140,9200,20)
# Optimum: {'min_child_weight': 9180.0}
# >>> space['min_child_weight']= quniform('min_child_weight', 9170,9190,5)
# Optimum: {'min_child_weight': 9175.0}
# >>> space['min_child_weight']= quniform('min_child_weight', 9172,9178,1)
# Optimum: {'min_child_weight': 9177.0}
# >>> space['min_child_weight']= quniform('min_child_weight', 9177,9179,1)
# Optimum: {'min_child_weight': 9177.0}
# >>> space['min_child_weight']=9177
# >>> space['alpha']= quniform('alpha',400,420,10)
# >>> space['lambda']=quniform('lambda',0.62,0.66,0.01)
# >>> space['gamma']=quniform('gamma',0.3,0.36,0.02)
# Optimum: {'alpha': 410.0, 'gamma': 0.3, 'lambda': 0.64}
# >>> space['alpha']= quniform('alpha',406,414,2)
# >>> space['lambda']=0.64
# >>> space['gamma']=quniform('gamma',0.28,0.31,0.01)
# Optimum: {'alpha': 408.0, 'gamma': 0.3}
# >>> space['gamma']=quniform('gamma',0.28,0.31,0.01)
# >>> space['alpha']= quniform('alpha',407,409,1)
# >>> space['lambda']=quniform('lambda',0.63,0.65,0.01)
# Optimum: {'alpha': 409.0, 'gamma': 0.28, 'lambda': 0.65}
# >>> space['alpha']= quniform('alpha',408,410,1)
# >>> space['lambda']=quniform('lambda',0.64,0.66,0.01)
# >>> space['gamma']= quniform('gamma',0.27,0.29,0.02)
# Optimum: {'alpha': 409.0, 'gamma': 0.28, 'lambda': 0.65}
# >>> space['alpha']=409
# >>> space['gamma']=0.28
# >>> space['lambda']=0.65
# >>> space['subsample']=quniform('subsample',0.3,0.4,0.01)
# >>> space['colsample_bytree']=hp.choice('colsample_bytree', [0.56,0.62,0.67,0.73,0.78,0.84,0.89])
# Optimum: {'colsample_bytree': 0.62, 'subsample': 0.35000000000000003}
# >>> space['subsample']=0.35
# >>> space['colsample_bytree']=quniform('colsample_bytree',0.6,0.7,0.01)
# Optimum: {'colsample_bytree': 0.61}
# space['eta']=quniform('eta',0.05,0.1,0.01)
# Optimum: {'eta': 0.05}
# >>> space['eta']=quniform('eta',0.01,0.05,0.01)
# Optimum: {'eta': 0.03}
# space['eta']=quniform('eta',0.025,0.035,0.001)
# Optimum: {'eta': 0.031}
#=====================================================================================================
# space['eta']=0.2
# space['lambda']= hp.choice('lambda',[0.7,1.,2.,5.,10.])
# space['gamma']= hp.choice('gamma',[0.3,1.,2.,5.,10.])
# Optimum: {'gamma': 2., 'lambda': 1.}
# >>> space['gamma']=hp.choice('gamma',[0.28,0.5,1.5,2.,2.5])
# >>> space['lambda']=hp.choice('lambda',[0.65,0.8,0.9,1.,1.5])
# Optimum: {'gamma': 0.5, 'lambda': 0.9}
# >>> space['alpha']= quniform('alpha',400,600,100)
# >>> space['lambda']=hp.choice('lambda',[0.65,0.8,0.9,1.,1.5])
# >>> space['gamma']=hp.choice('gamma',[0.28,0.5,1.5,2.,2.5])
# Optimum: {'alpha': 400.0, 'gamma': 2.5, 'lambda': 0.65}
# >>> space['alpha']= quniform('alpha',390,410,10)
# >>> space['gamma']=hp.choice('gamma',[2.5,5.,7.5,10.])
# >>> space['lambda']=hp.choice('lambda',[0.5,0.6,0.65,0.8])
# Optimum: {'alpha': 390.0, 'gamma': 7.5, 'lambda': 0.65}
# >>> space['gamma']=hp.choice('gamma',[6.5,7.,7.5,8.])
# >>> space['alpha']= quniform('alpha',382,390,4)
# >>> space['lambda']=hp.choice('lambda',[0.63,0.65,0.67])
# Optimum: {'alpha': 388.0, 'gamma': 6.5, 'lambda': 0.63}
# >>> space['alpha']= quniform('alpha',386,390,2)
# >>> space['gamma']=hp.choice('gamma',[6.,6.5,7.])
# >>> space['lambda']=hp.choice('lambda',[0.61,0.62,0.63])
# Optimum: {'alpha': 386.0, 'gamma': 6.5, 'lambda': 0.62}
# >>> space['alpha']= quniform('alpha',370,385,5)
# >>> space['gamma']=hp.choice('gamma',[6.4,6.5,6.6])
# >>> space['lambda']=hp.choice('lambda',[0.615,0.62,0.625])
# Optimum: {'alpha': 380.0, 'gamma': 6.4, 'lambda': 0.625}
# >>> space['alpha']= quniform('alpha',376,384,2)
# >>> space['gamma']=quniform('gamma',6.2,6.5,0.1)
# >>> space['lambda']=quniform('lambda',0.62,0.63,0.02)
# Optimum: {'alpha': 384.0, 'gamma': 6.4, 'lambda': 0.62}
# >>> space['gamma']=6.4
# >>> space['alpha']= quniform('alpha',382,388,2)
# >>> space['lambda']=quniform('lambda',0.6,0.64,0.01)
# Optimum: {'alpha': 388.0, 'lambda': 0.63}
# >>> space['lambda']=quniform('lambda',0.62,0.63,0.005)
# >>> space['alpha']= quniform('alpha',387,390,1)
# Optimum: {'alpha': 387.0, 'lambda': 0.63}
# >>> space['lambda']=quniform('lambda',0.4,0.6,0.1)
# >>> space['alpha']= quniform('alpha',350,390,10)
# Optimum: {'alpha': 380.0, 'lambda': 0.5}
# >>> space['lambda']=quniform('lambda',0.46,0.54,0.02)
# >>> space['alpha']= quniform('alpha',376,384,2)
# Optimum: {'alpha': 382.0, 'lambda': 0.52}
# >>> space['lambda']=quniform('lambda',0.51,0.53,0.01)
# >>> space['alpha']= quniform('alpha',381,383,1)
# Optimum: {'alpha': 381.0, 'lambda': 0.51}
# >>> space['alpha']= quniform('alpha',380,384,1)
# >>> space['lambda']=quniform('lambda',0.50,0.54,0.01)
# Optimum: {'alpha': 384.0, 'lambda': 0.5}
# >>> space['lambda']=quniform('lambda',0.45,0.52,0.01)
# >>> space['alpha']= quniform('alpha',376,384,1)
# Optimum: {'alpha': 378.0, 'lambda': 0.52}
# >>> space['alpha']= 378
# >>> space['lambda']=quniform('lambda',0.51,0.58,0.01)
# Optimum: {'lambda': 0.58}
# >>> space['lambda']=quniform('lambda',0.6,0.8,0.02)
# Optimum: {'lambda': 0.8}
# >>> space['lambda']=quniform('lambda',0.6,1.5,0.1)
# Optimum: {'lambda': 1.0}
# >>> space['lambda']=quniform('lambda',0.9,1.1,0.02)
# >>> space['alpha']= quniform('alpha',377,379,1)
# Optimum: {'alpha': 379.0, 'lambda': 0.96}
# >>> space['alpha']= quniform('alpha',380,400,4)
# >>> space['lambda']=quniform('lambda',0.94,0.98,0.02)
# Optimum: {'alpha': 396.0, 'lambda': 0.96}
# >>> space['lambda']=quniform('lambda',0.95,0.97,0.01)
# >>> space['alpha']= quniform('alpha',394,398,1)
# Optimum: {'alpha': 394.0, 'lambda': 0.96}
# >>> space['lambda']=quniform('lambda',0.956,0.964,0.004)
# >>> space['alpha']= quniform('alpha',392,396,1)
# Optimum: {'alpha': 393.0, 'lambda': 0.964}
# >>> space['alpha']= quniform('alpha',392,394,1)
# >>> space['lambda']=quniform('lambda',0.956,0.964,0.002)
# >>> space['gamma']=quniform('gamma',6.3,6.5,0.1)
# Optimum: {'alpha': 393.0, 'gamma': 6.4, 'lambda': 0.96}
# >>> space['gamma']=6.4
# >>> space['lambda']=0.96
# >>> space['alpha']=393
# >>> space['subsample']=quniform('subsample',0.3,0.4,0.01)
# >>> space['colsample_bytree']=hp.choice('colsample_bytree',[0.56,0.61,0.67,0.72])
# Optimum: {'colsample_bytree': 0.61, 'subsample': 0.35000000000000003}
# >>> space['colsample_bytree']=0.61
# >>> space['subsample']=0.35
# >>> space['eta']=quniform('eta',0.02,0.04,0.01)
# Optimum: {'eta': 0.04}
# >>> space['eta']=quniform('eta',0.036,0.044,0.001)
# Optimum: {'eta': 0.043000000000000003}
















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

# 			'max_depth': hp.choice('max_depth', [6,8,10]),
# 			'min_child_weight': quniform('min_child_weight', 0,40,10),

# 			'alpha' : 0.,
# 			'lambda' : 0.,
# 			'gamma': 0.,
						
# 			'subsample': 0.6,
# 			'colsample_bytree':	0.8,

# 			'eta': 0.05,

# 			'scale_pos_weight': 1.
# }
# {'max_depth': 2 -> 10, 'min_child_weight': 40.0}

# >>> space['min_child_weight']=quniform('min_child_weight', 40,80,20)
# >>> space['max_depth']=uniformint('max_depth',6,8)
# {'max_depth': 8.0, 'min_child_weight': 80.0}

# 'max_depth': uniformint('max_depth',6,10),
# 'min_child_weight': quniform('min_child_weight', 40,160,20),
# -> {'max_depth': 10.0, 'min_child_weight': 140.0}


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

# 			'eta': 0.2,

# }
# 											{'min_child_weight': 13} -> 500
#
# 'min_child_weight': hp.choice('min_child_weight', [500,800,1000,1500])
# Optimum: {'min_child_weight': 1500}
# space['min_child_weight']= quniform('min_child_weight', 1500,10000,500)
# Optimum: {'min_child_weight': 9500.0}
# space['min_child_weight']= quniform('min_child_weight', 9000,10200,400)
# Optimum: {'min_child_weight': 9200.0}
# space['min_child_weight']= quniform('min_child_weight', 9120,9300,30)
# Optimum: {'min_child_weight': 9300.0}
# space['min_child_weight']= quniform('min_child_weight', 9250,9400,50)
# Optimum: {'min_child_weight': 9300.0}
# space['min_child_weight']= quniform('min_child_weight', 9280,9320,10)
# Optimum: {'min_child_weight': 9310.0}
# space['min_child_weight']= quniform('min_child_weight', 9300,9320,4)
# Optimum: {'min_child_weight': 9304.0}
# space['min_child_weight']= quniform('min_child_weight', 9301,9307,1)
# Optimum: {'min_child_weight': 9303.0}
# space['min_child_weight']=9303
# space['gamma']= hp.choice('gamma',[0.,0.2,0.5])
# space['lambda']= hp.choice('lambda',[0.,0.5,0.1])
# space['alpha']= hp.choice('alpha',[0.,0.1,1,100,1000])
# {'alpha': 100, 'gamma': 0.2, 'lambda': 0.5}
# space['gamma']=quniform('gamma',0.1,0.3,0.1)
# space['alpha']= quniform('alpha',80,240,80)
# space['lambda']=quniform('lambda',0.3,0.7,0.2)
# Optimum: {'alpha': 240.0, 'gamma': 0.30000000000000004, 'lambda': 0.6000000000000001}
# space['alpha']= quniform('alpha',200,300,50)
# space['gamma']=quniform('gamma',0.25,0.45,0.1)
# space['lambda']=quniform('lambda',0.4,0.9,0.1)
# Optimum: {'alpha': 300.0, 'gamma': 0.4, 'lambda': 0.7000000000000001}
# space['alpha']= quniform('alpha',300,500,100)
# space['gamma']=quniform('gamma',0.3,0.6,0.1)
# space['lambda']=quniform('lambda',0.6,0.8,0.1)
# Optimum: {'alpha': 400.0, 'gamma': 0.30000000000000004, 'lambda': 0.8}
# space['alpha']= quniform('alpha',340,460,30)
# space['gamma']=quniform('gamma',0.25,0.35,0.05)
# space['lambda']=quniform('lambda',0.6,0.9,0.1)
# Optimum: {'alpha': 390.0, 'gamma': 0.30000000000000004, 'lambda': 0.7000000000000001}
# space['gamma']=quniform('gamma',0.28,0.32,0.02)
# space['lambda']=quniform('lambda',0.6,0.8,0.04)
# space['alpha']= quniform('alpha',380,400,5)
# {'alpha': 400.0, 'gamma': 0.28, 'lambda': 0.64}
#  space['gamma']=quniform('gamma',0.25,0.29,0.01)
#  space['lambda']=quniform('lambda',0.62,0.66,0.01)
#  space['alpha']= quniform('alpha',396,404,2)
# Optimum: {'alpha': 398.0, 'gamma': 0.28, 'lambda': 0.64}
# space['gamma']=0.28
# space['lambda']=0.64
# space['alpha']= quniform('alpha',397,399,1)
# Optimum: {'alpha': 399.0}
# space['alpha']= quniform('alpha',399,401,1)
# Optimum: {'alpha': 399.0}
# space['alpha']=399
# space['colsample_bytree']=quniform('colsample_bytree',0.111,0.999,0.1111)
# space['subsample']=quniform('subsample',0.1,0.9,0.1)
# Optimum: {'colsample_bytree': 0.8888, 'subsample': 0.4}
# space['subsample']=quniform('subsample',0.36,0.44,0.02)
# space['colsample_bytree']=hp.choice('colsample_bytree',[0.84,0.89,0.94])
# Optimum: {'colsample_bytree': 0.84, 'subsample': 0.36}
# space['colsample_bytree']=hp.choice('colsample_bytree',[0.72,0.78,0.84,0.89])
# space['subsample']=quniform('subsample',0.34,0.38,0.01)
# Optimum: {'colsample_bytree': 0.72, 'subsample': 0.35000000000000003}
# space['subsample']=0.35
# space['colsample_bytree']=hp.choice('colsample_bytree',[0.67,0.73,0.78])
# Optimum: {'colsample_bytree': 78}
# space['colsample_bytree']=hp.choice('colsample_bytree',[0.73,0.78,0.84])
# Optimum: {'colsample_bytree': 0.78}
# space['colsample_bytree']=0.78
# space['alpha']= quniform('alpha',398,401,1)
# space['lambda']=quniform('lambda',0.62,0.66,0.01)
# Optimum: {'alpha': 400.0, 'lambda': 0.66}
# space['lambda']=quniform('lambda',0.65,0.67,0.01)
# space['alpha']= quniform('alpha',400,402,1)
# Optimum: {'alpha': 402.0, 'lambda': 0.66}
# space['lambda']=0.66
# space['alpha']= quniform('alpha',402,410,2)
# Optimum: {'alpha': 410.0}
# space['alpha']= hp.choice('alpha',[400,500,800,1000])
# Optimum: {'alpha': 400}
# space['alpha']= quniform('alpha',400,410,1)
# Optimum: {'alpha': 406.0}
# space['eta']=hp.choice('eta', [0.01,0.02,0.04,0.06])
# Optimum: {'eta': 0.06}
# space['scale_pos_weight']=hp.choice('scale_pos_weight', [1,5,10,100])
# Optimum: {'scale_pos_weight': 1}
# space['scale_pos_weight']=hp.choice('scale_pos_weight', [0.1,0.5,1])
# Optimum: {'scale_pos_weight': 0.5}
# space['scale_pos_weight']=hp.choice('scale_pos_weight', [0.5,0.6,0.7,0.8,0.9,1])
# Optimum: {'scale_pos_weight': 0.6}
# space['scale_pos_weight']=quniform('scale_pos_weight', 0.5,0.7,0.02)
# Optimum: {'scale_pos_weight': 0.7000000000000001}
