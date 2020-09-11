#!/usr/bin/env python3

####################################################################
# Photon ID BDT for aNTGC search in Z(->nu nu) + Gamma channel     #
# Mohammad Abrar Wadud (BD), 04/02/2020                            #
###################################################################

# configure options
from argparse import ArgumentParser
parser = ArgumentParse

	description='Photon ID BDT Training for CMS aNTGC search in Z(->nu nu) + gamma channel')
parser.add_argument('--inFilePath', type=str, help='Input root file',
					default='/home/rusack/wadud/phoBDT/mergedSamplesShuffled.root', action='store')
parser.add_argument('--modelFilePath', type=str, help='BDT model file',
					default='/home/rusack/wadud/phoBDT/optimizedV0/aNTGC_photon_BDT.pkl', action='store')
parser.add_argument('--saveDir', type=str, default='/home/rusack/wadud/phoBDT/optimizedV0/',
					help='Save directory', action='store')
parser.add_argument('--inTreeName', type=str, default='fullEB_Tree',
					help='Input tree name', action='store')
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
		self.log = open("%s" % logFileName, "a")

	def write(self, message):
		self.terminal.write(message)
		self.log.write(message)

	def flush(self):
		pass

sys.stdout = Logger()

print("************************************************************************************************************************************************\n")
now = datetime.datetime.now()
print(now.strftime("%Y-%m-%d %H:%M:%S"))
print("Starting photon ID training BDT\n")
print("inFilePath\t=\t" + args.inFilePath)
print("modelFilePath\t=\t" + args.modelFilePath)
print("saveDir\t\t=\t" + args.saveDir)
print("inTreeName\t=\t" + args.inTreeName)
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


BDTfeats = ["phoR9Full5x5", "phoS4Full5x5", "phoEmaxOESCrFull5x5", "phoE2ndOESCrFull5x5", "phoE2ndOEmaxFull5x5", "phoE1x3OESCrFull5x5", "phoE2x5OESCrFull5x5", "phoE5x5OESCrFull5x5",
			"phoEmaxOE3x3Full5x5", "phoE2ndOE3x3Full5x5", "phoSigmaIEtaIEta", "phoSigmaIEtaIPhi", "phoSigmaIPhiIPhi", "phoSieieOSipipFull5x5", "phoEtaWidth", "phoPhiWidth", "phoEtaWOPhiWFull5x5"]
BDTfeats.sort()
print("\nBDT input features (" + str(len(BDTfeats)) + ") :")
print(BDTfeats)

phoModelIterated = pickle.load(open(args.modelFilePath, "rb"))

totalData = ROOT.RDataFrame(args.inTreeName, args.inFilePath)
print('\nGetting predictions @ ' + now.strftime("%Y-%m-%d %H:%M:%S") + "\n")


# totalData =	totalData.Range(0, 1000)
totalData = pd.DataFrame(totalData.AsNumpy(columns=BDTfeats + ['PtEtaRwBG', 'xSecW', 'isSignal', 'isTrain', 'isValidation', 'splitRand']))
sampleStats(totalData)

totalData['bdtWeight'] = totalData['xSecW'] * totalData['PtEtaRwBG']

totalDataDM = xg.DMatrix(totalData[BDTfeats].values, label=totalData['isSignal'].values, feature_names=BDTfeats, nthread=-1)

saveDF = totalData[['isSignal', 'isTrain', 'isValidation', 'bdtWeight', 'splitRand']].copy()

booleanDictionary = {True: 1, False: 0}
totalData = totalData.replace(booleanDictionary)

saveDF['isTrain'] = saveDF['isTrain'].astype('int')
saveDF['isSignal'] = saveDF['isSignal'].astype('int')
saveDF['isValidation'] = saveDF['isValidation'].astype('int')
saveDF['bdtWeight'] = saveDF['bdtWeight'].astype('int')

saveDF['bdtScore'] = phoModelIterated.predict(totalDataDM)

del totalData

outFileName = args.saveDir + '/' + 'BDTresults' + '.csv'
saveDF.to_csv(outFileName, encoding='utf-8', index=False)
print ('Data frame saved in %s' % outFileName)

del saveDF

saveDF = ROOT.RDF.MakeCsvDataFrame(outFileName);
outFileName = args.saveDir + '/' + 'BDTresults' + '.root'
saveDF.Snapshot(args.outTreeName, outFileName)
print ('Data frame saved in %s' % outFileName)


print("\n************************************************************************************************************************************************")
now = datetime.datetime.now()
print("Done @ " + now.strftime("%Y-%m-%d %H:%M:%S"))
print("************************************************************************************************************************************************")
