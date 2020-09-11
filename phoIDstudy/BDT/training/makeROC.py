#!/usr/bin/env python3

####################################################################
# Photon ID BDT for aNTGC search in Z(->nu nu) + Gamma channel     #
# Mohammad Abrar Wadud (BD), 05/06/2020                            #
####################################################################

import sys
from sklearn.metrics import roc_curve
import numpy as np
from custom_auc import aucW
import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd
# from root_pandas.readwrite import read_root
import ROOT

saveDir='dataV2/training/'
bdtFile='dataV2/training/BDTresults_2020_09_10_23_04_13.root'

columns=['isSignalF', 'isTrainF', 'isValidationF', 'isTestF', 'flatPtEtaRwNoXsecF', 'bdtV2Score']
data = ROOT.RDataFrame('BDTresults', bdtFile)
# data = read_root(paths=bdtFile, key='BDTresults', columns=columns)
data = pd.DataFrame(data.AsNumpy(columns=columns))
data = data.sample(frac = 1)

print("\nTotal data:")
sigCount = len(data[data['isSignalF'] == 1].index)
bkgCount = len(data[data['isSignalF'] == 0].index)
print('Total Events\t=\t%d (B/S= %d / %d)' %  ((sigCount + bkgCount), bkgCount, sigCount))

print("Data loaded!")

data['flatPtEtaRwNoXsecF'] = data['flatPtEtaRwNoXsecF'].astype('float')
data['bdtV2Score'] = data['bdtV2Score'].astype('float')
data['isTrainF'] = data['isTrainF'].astype('bool')
data['isSignalF'] = data['isSignalF'].astype('bool')
data['isValidationF'] = data['isValidationF'].astype('bool')


trainY = data[(data.isTrainF == 1)]['isSignalF']
trainPredY = data[(data.isTrainF == 1)]['bdtV2Score']
trainW = data[(data.isTrainF == 1)]['flatPtEtaRwNoXsecF']
trainAUC = aucW(trainY, trainPredY, trainW)
print("Train AUC computed!")

validationY = data[(data.isValidationF==1)]['isSignalF']
validationPredY = data[(data.isValidationF==1)]['bdtV2Score']
validationW = data[(data.isValidationF==1)]['flatPtEtaRwNoXsecF']
validationAUC = aucW(validationY, validationPredY, validationW)
print("Validation AUC computed!")

testY = data[data.isTestF == 1]['isSignalF']
testPredY = data[data.isTestF == 1]['bdtV2Score']
testW = data[data.isTestF == 1]['flatPtEtaRwNoXsecF']
testAUC = aucW(testY, testPredY, testW)
print("Test AUC computed!")


print("Making plot!")

mpl.use('Agg')
rcParams.update({'figure.autolayout': True})
plt.rc('legend', fontsize='medium')

plt.figure(1)

bkgEff, sigEff, thres = roc_curve(trainY, trainPredY, sample_weight=trainW)
plt.plot(bkgEff, sigEff, label='Train: AUC=%1.5f, S=%.0f(%.0f%%), B=%.0f(%.0f%%)' %
		 (trainAUC, np.count_nonzero(trainY == 1), 100. * np.count_nonzero(trainY == 1) / len(trainY), np.count_nonzero(trainY == 0), 100. * np.count_nonzero(trainY == 0) / len(trainY)), color='#034e7b')


bkgEff, sigEff, thres = roc_curve(validationY, validationPredY, sample_weight=validationW)
plt.plot(bkgEff, sigEff, label='Validation: AUC=%1.5f, S=%.0f(%.0f%%), B=%.0f(%.0f%%)' %
		 (validationAUC, np.count_nonzero(validationY == 1), 100. * np.count_nonzero(validationY == 1) / len(validationY), np.count_nonzero(validationY == 0), 100. * np.count_nonzero(validationY == 0) / len(validationY)), color='#1b7837')


bkgEff, sigEff, thres = roc_curve(testY, testPredY, sample_weight=testW)
plt.plot(bkgEff, sigEff, label='Test: AUC=%1.5f, S=%.0f(%.0f%%), B=%.0f(%.0f%%)' %
		 (testAUC, np.count_nonzero(testY == 1), 100. * np.count_nonzero(testY == 1) / len(testY), np.count_nonzero(testY == 0), 100. * np.count_nonzero(testY == 0) / len(testY)), color='#d73027')

plt.legend(loc="lower right")
plt.grid()
plt.xlabel('Background Efficiency')
plt.ylabel('Signal Efficiency')

plt.ylim(0.65, 1.01)
plt.xlim(0., 0.05)
plt.savefig('%s/ROC1.pdf' % saveDir, dpi=600)
plt.savefig('%s/ROC1.png' % saveDir, dpi=600)


plt.ylim(0.9, 1.01)
plt.xlim(0., 0.05)
plt.savefig('%s/ROC2.pdf' % saveDir, dpi=600)
plt.savefig('%s/ROC2.png' % saveDir, dpi=600)

plt.ylim(0., 1.05)
plt.xlim(0., 1.05)
plt.savefig('%s/ROCfull1.pdf' % saveDir, dpi=600)
plt.savefig('%s/ROCfull1.png' % saveDir, dpi=600)


print("Done!" )