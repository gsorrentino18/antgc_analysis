#!/usr/bin/env python3

####################################################################
# Photon ID BDT for aNTGC search in Z(->nu nu) + Gamma channel     #
# Mohammad Abrar Wadud (BD), 05/06/2020                            #
####################################################################

import sys
sys.path.insert(1, 'training/')

from sklearn.metrics import roc_auc_score, roc_curve
import numpy as np
from custom_auc import aucW
import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib import rc
from root_pandas.readwrite import read_root, to_root
import pandas as pd

saveDir='/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/BDT/plots/'
featFile='/hdfs/cms/user/wadud/anTGC/BDTdata/mergedSamplesShuffled.root'
bdtFile='/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/BDT/data/optimizedV1/BDTresults_0.root'

data = read_root(paths=featFile, key='fullEB_Tree', columns=['isSignal', 'isTrain', 'isValidation', 'phoBDTpredictionHoE', 'phoIDbit', 'splitRand'])
bdtData = read_root(paths=bdtFile, key='fullEB_BDT_Tree', columns=['bdtWeightF', 'bdtScore'])
data['bdtWeightF'] = bdtData['bdtWeightF'].to_numpy()
data['bdtScore'] = bdtData['bdtScore'].to_numpy()
del bdtData

print("Data loaded!")

data['bdtWeightF'] = data['bdtWeightF'].astype('float')
data['bdtScore'] = data['bdtScore'].astype('float')
data['isTrain'] = data['isTrain'].astype('bool')
data['isSignal'] = data['isSignal'].astype('bool')
data['isValidation'] = data['isValidation'].astype('bool')

sumSig = data[(data.isSignal==1) & (data.isTrain == 0) & (data.isValidation==0)]['bdtWeightF'].sum()
sumBg = data[(data.isSignal==0) & (data.isTrain == 0)& (data.isValidation==0)]['bdtWeightF'].sum()

trainY = data[(data.isTrain == 1) & (data.isValidation==0)]['isSignal']
trainPredY = data[(data.isTrain == 1) & (data.isValidation==0)]['bdtScore']
trainW = data[(data.isTrain == 1) & (data.isValidation==0)]['bdtWeightF']
trainAUC = aucW(trainY, trainPredY, trainW)
print("Train AUC computed!")

validationY = data[(data.isValidation==1)]['isSignal']
validationPredY = data[(data.isValidation==1)]['bdtScore']
validationW = data[(data.isValidation==1)]['bdtWeightF']
validationAUC = aucW(validationY, validationPredY, validationW)
print("Validation AUC computed!")

testY = data[data.isTrain == 0]['isSignal']
testPredY = data[data.isTrain == 0]['bdtScore']
testW = data[data.isTrain == 0]['bdtWeightF']
testAUC = aucW(testY, testPredY, testW)
print("Test AUC computed!")

# oldBDTPredY = data[data.isTrain == 0]['phoBDTpredictionHoE']
# oldBDTAUC = aucW(testY, oldBDTPredY, testW)
# print("Old AUC computed!")


# data['looseID'] = data['phoIDbit'] & (1 << 0)
# data['looseID'] = data['looseID'].astype('bool')
# lSigEff	= data[(data.isSignal==1) & (data.looseID==1) & (data.isTrain == 0)]['bdtWeightF'].sum()/sumSig
# lBgEff	= data[(data.isSignal==0) & (data.looseID==1) & (data.isTrain == 0)]['bdtWeightF'].sum()/sumBg
# print("Loose efficiencies calculated!")

# data['mediumID'] = data['phoIDbit'] & (1 << 1)
# data['mediumID'] = data['mediumID'].astype('bool')
# mSigEff	= data[(data.isSignal==1) & (data.mediumID==1) & (data.isTrain == 0)]['bdtWeightF'].sum()/sumSig
# mBgEff	= data[(data.isSignal==0) & (data.mediumID==1) & (data.isTrain == 0)]['bdtWeightF'].sum()/sumBg
# print("Medium efficiencies calculated!")

# data['tightID'] = data['phoIDbit'] & (1 << 2)
# data['tightID'] = data['tightID'].astype('bool')
# tSigEff	= data[(data.isSignal==1) & (data.tightID==1) & (data.isTrain == 0)]['bdtWeightF'].sum()/sumSig
# tBgEff	= data[(data.isSignal==0) & (data.tightID==1) & (data.isTrain == 0)]['bdtWeightF'].sum()/sumBg
# print("Tight efficiencies calculated!")

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

# bkgEff, sigEff, thres = roc_curve(testY, oldBDTPredY, sample_weight=testW)
# plt.plot(bkgEff, sigEff, label='H/E BDT: AUC=%1.5f, S=%.0f(%.0f%%), B=%.0f(%.0f%%)' %
# 		 (oldBDTAUC, np.count_nonzero(testY == 1), 100. * np.count_nonzero(testY == 1) / len(testY), np.count_nonzero(testY == 0), 100. * np.count_nonzero(testY == 0) / len(testY)), color='#238443')

# plt.plot([lBgEff], [lSigEff], marker='o', markersize=3, color="black")
# plt.text(lBgEff*0.98,lSigEff,'Loose(%1.3f,%1.3f)' %(lBgEff, lSigEff),horizontalalignment='left', va='bottom', fontsize=9)

# plt.plot([mBgEff], [mSigEff], marker='o', markersize=3, color="black")
# plt.text(mBgEff*0.98,mSigEff,'Medium(%1.3f,%1.3f)' %(mBgEff, mSigEff),horizontalalignment='left', va='bottom', fontsize=9)

# plt.plot([tBgEff], [tSigEff], marker='o', markersize=3, color="black")
# plt.text(tBgEff*0.98,tSigEff,'Tight(%1.3f,%1.3f)' %(tBgEff, tSigEff),horizontalalignment='left', va='bottom', fontsize=9)

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