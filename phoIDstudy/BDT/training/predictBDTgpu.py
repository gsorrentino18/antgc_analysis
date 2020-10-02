#!/usr/bin/env python3

####################################################################
# Photon ID BDT for aNTGC search in Z(->nu nu) + Gamma channel     #
# Mohammad Abrar Wadud (BD), 04/02/2020                            #
###################################################################

# configure options
from argparse import ArgumentParser
parser = ArgumentParser(description='Photon ID BDT Training for CMS aNTGC search in Z(->nu nu) + gamma channel')
parser.add_argument('--inFilePath', type=str, help='Input root file',
					default='data/2020_09_15_Prompt_AlternateSamples/ZGTo2NuGPtG130TuneCP513TeVamcatnloFXFXpythia8.root', action='store')
parser.add_argument('--modelFilePath', type=str, help='BDT model file',
					default='/home/rusack/wadud/phoBDT/training/trainingV2/aNTGC_photon_BDT_2020_09_12_19_08_19.pkl', action='store')
parser.add_argument('--saveDir', type=str, default='data/2020_09_15_Prompt_AlternateSamples/predictions/',
					help='Save directory', action='store')
parser.add_argument('--inTreeName', type=str, default='fullEB/fullEBTree',
					help='Input tree name', action='store')
parser.add_argument('--outTreeName', default='BDTresults',
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
import json
import ROOT
import math
import ntpath

BDTfeats = ["phoR9Full5x5", "phoS4Full5x5", "phoEmaxOESCrFull5x5", "phoE2ndOESCrFull5x5", "pho2x2OE3x3Full5x5", "phoE1x3OESCrFull5x5", "phoE2x5OESCrFull5x5", "phoE5x5OESCrFull5x5",
			"phoSigmaIEtaIEta", "phoSigmaIEtaIPhi", "phoSigmaIPhiIPhi", "phoSieieOSipipFull5x5", "phoEtaWidth", "phoPhiWidth", "phoEtaWOPhiWFull5x5"]
BDTfeats.sort()
print("\nBDT input features (" + str(len(BDTfeats)) + ") :")
print(BDTfeats)

phoModel = pickle.load(open(args.modelFilePath, "rb"))

totalData = ROOT.RDataFrame(args.inTreeName, args.inFilePath)
print('\nGetting predictions @ ' + now.strftime("%Y-%m-%d %H:%M:%S") + "\n")


# totalData =	totalData.Range(0, 1000)
totalData = pd.DataFrame(totalData.AsNumpy(columns=BDTfeats))
# sampleStats(totalData)

totalDataDM = xg.DMatrix(totalData[BDTfeats].values, feature_names=BDTfeats, nthread=-1)

saveDF = pd.DataFrame()

saveDF['bdtScore'] = phoModel.predict(totalDataDM)

del totalData

outFileName = args.saveDir + ntpath.basename(args.inFilePath)
outFileName.replace('.root', '.csv')
saveDF.to_csv(outFileName, encoding='utf-8', index=False)
print ('Data frame saved in %s' % outFileName)

del saveDF

saveDF = ROOT.RDF.MakeCsvDataFrame(outFileName);


outFileName = args.saveDir + ntpath.basename(args.inFilePath)
saveDF.Snapshot(args.outTreeName, outFileName)
print ('Data frame saved in %s' % outFileName)


print("\n************************************************************************************************************************************************")
now = datetime.datetime.now()
print("Done @ " + now.strftime("%Y-%m-%d %H:%M:%S"))
print("************************************************************************************************************************************************")
