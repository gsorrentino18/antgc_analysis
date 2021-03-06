#!/usr/bin/env python3

####################################################################
# Photon ID BDT for aNTGC search in Z(->nu nu) + Gamma channel     #
# Mohammad Abrar Wadud (BD), 05/06/2020                            #
####################################################################

import sys
sys.path.insert(1, 'training/')

from sklearn.metrics import roc_auc_score, roc_curve
import numpy as np
import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib import rc
from root_pandas.readwrite import read_root, to_root
import pandas as pd

saveDir='/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/efficiency/plots/'
featFile='/hdfs/cms/user/wadud/anTGC/BDTdata/mergedSamplesShuffled.root'
bdtFile='/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/BDT/data/optimizedV1/BDTresults_0.root'
idFile='/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/efficiency/data/idShuffledTree.root'

data = read_root(paths=featFile, key='fullEB_Tree', columns=['isSignal', 'isTrain', 'isValidation', 'phoIDbit', 'phoSCeta'])
bdtData = read_root(paths=bdtFile, key='fullEB_BDT_Tree', columns=['bdtWeightF', 'bdtScore'])
idData = read_root(paths=idFile, key='fullEB_IdTree', columns=['pass95', 'pass90', 'pass80', 'pass70', 'cutbased95', 'cutbased90', 'cutbased80', 'cutbased70'])
data['bdtWeightF'] = bdtData['bdtWeightF'].to_numpy()
data['bdtScore'] = bdtData['bdtScore'].to_numpy()
data['pass95'] = idData['pass95'].to_numpy()
data['pass90'] = idData['pass90'].to_numpy()
data['pass80'] = idData['pass80'].to_numpy()
data['pass70'] = idData['pass70'].to_numpy()

data['cutbased95'] = idData['cutbased95'].to_numpy()
data['cutbased90'] = idData['cutbased90'].to_numpy()
data['cutbased80'] = idData['cutbased80'].to_numpy()
data['cutbased70'] = idData['cutbased70'].to_numpy()

del bdtData
del idData

print('Data loaded!')

data['bdtWeightF'] = data['bdtWeightF'].astype('float')
data['bdtScore'] = data['bdtScore'].astype('float')
data['isTrain'] = data['isTrain'].astype('bool')
data['isSignal'] = data['isSignal'].astype('bool')
data['isValidation'] = data['isValidation'].astype('bool')
data['pass95'] = data['pass95'].astype('bool')
data['pass90'] = data['pass90'].astype('bool')
data['pass80'] = data['pass80'].astype('bool')
data['pass70'] = data['pass70'].astype('bool')
data['cutbased95'] = data['cutbased95'].astype('bool')
data['cutbased90'] = data['cutbased90'].astype('bool')
data['cutbased80'] = data['cutbased80'].astype('bool')
data['cutbased70'] = data['cutbased70'].astype('bool')

sumSig = data[(data.isSignal==1) & (data.isTrain == 0) & (abs(data.phoSCeta) <= 1.4442)]['bdtWeightF'].sum()
sumBg = data[(data.isSignal==0) & (data.isTrain == 0) & (abs(data.phoSCeta) <= 1.4442)]['bdtWeightF'].sum()

llSigEff = 100.*data[(data.isSignal==1) & (data.pass95==1) & (data.isTrain == 0) & (abs(data.phoSCeta) <= 1.4442)]['bdtWeightF'].sum()/sumSig
llBgEff	= 100.*data[(data.isSignal==0) & (data.pass95==1) & (data.isTrain == 0) & (abs(data.phoSCeta) <= 1.4442)]['bdtWeightF'].sum()/sumBg
print('95 %% efficiencies calculated!')

lSigEff	= 100.*data[(data.isSignal==1) & (data.pass90==1) & (data.isTrain == 0) & (abs(data.phoSCeta) <= 1.4442)]['bdtWeightF'].sum()/sumSig
lBgEff	= 100.*data[(data.isSignal==0) & (data.pass90==1) & (data.isTrain == 0) & (abs(data.phoSCeta) <= 1.4442)]['bdtWeightF'].sum()/sumBg
print('90 %% efficiencies calculated!')


mSigEff	= 100.*data[(data.isSignal==1) & (data.pass80==1) & (data.isTrain == 0) & (abs(data.phoSCeta) <= 1.4442) ]['bdtWeightF'].sum()/sumSig
mBgEff	= 100.*data[(data.isSignal==0) & (data.pass80==1) & (data.isTrain == 0) & (abs(data.phoSCeta) <= 1.4442) ]['bdtWeightF'].sum()/sumBg
print('80 %% efficiencies calculated!')

tSigEff	= 100.*data[(data.isSignal==1) & (data.pass70==1) & (data.isTrain == 0) & (abs(data.phoSCeta) <= 1.4442) ]['bdtWeightF'].sum()/sumSig
tBgEff	= 100.*data[(data.isSignal==0) & (data.pass70==1) & (data.isTrain == 0) & (abs(data.phoSCeta) <= 1.4442)]['bdtWeightF'].sum()/sumBg
print('70 %% efficiencies calculated!')



llSigCBEff	= 100.*data[(data.isSignal==1) & (data.cutbased95==1) & (data.isTrain == 0) & (abs(data.phoSCeta) <= 1.4442)]['bdtWeightF'].sum()/sumSig
llBgCBEff	= 100.*data[(data.isSignal==0) & (data.cutbased95==1) & (data.isTrain == 0) & (abs(data.phoSCeta) <= 1.4442)]['bdtWeightF'].sum()/sumBg
print('Cut based 95 %% Efficiencies calculated!')

lSigCBEff	= 100.*data[(data.isSignal==1) & (data.cutbased90==1) & (data.isTrain == 0) & (abs(data.phoSCeta) <= 1.4442)]['bdtWeightF'].sum()/sumSig
lBgCBEff	= 100.*data[(data.isSignal==0) & (data.cutbased90==1) & (data.isTrain == 0) & (abs(data.phoSCeta) <= 1.4442)]['bdtWeightF'].sum()/sumBg
print('Cut based 90 %% Efficiencies calculated!')


mSigCBEff	= 100.*data[(data.isSignal==1) & (data.cutbased80==1) & (data.isTrain == 0) & (abs(data.phoSCeta) <= 1.4442) ]['bdtWeightF'].sum()/sumSig
mBgCBEff	= 100.*data[(data.isSignal==0) & (data.cutbased80==1) & (data.isTrain == 0) & (abs(data.phoSCeta) <= 1.4442) ]['bdtWeightF'].sum()/sumBg
print('Cut based 80 %% Efficiencies calculated!')

tSigCBEff	= 100.*data[(data.isSignal==1) & (data.cutbased70==1) & (data.isTrain == 0) & (abs(data.phoSCeta) <= 1.4442) ]['bdtWeightF'].sum()/sumSig
tBgCBEff	= 100.*data[(data.isSignal==0) & (data.cutbased70==1) & (data.isTrain == 0) & (abs(data.phoSCeta) <= 1.4442)]['bdtWeightF'].sum()/sumBg
print('Cut based 70 %% Efficiencies calculated!')


data['looseID'] = data['phoIDbit'] & (1 << 0)
data['looseID'] = data['looseID'].astype('bool')
lEGMSigEff	= 100.*data[(data.isSignal==1) & (data.looseID==1) & (data.isTrain == 0) & (abs(data.phoSCeta) <= 1.4442)]['bdtWeightF'].sum()/sumSig
lEGMBgEff	= 100.*data[(data.isSignal==0) & (data.looseID==1) & (data.isTrain == 0) & (abs(data.phoSCeta) <= 1.4442)]['bdtWeightF'].sum()/sumBg
print('EGM:L efficiencies calculated!')

data['mediumID'] = data['phoIDbit'] & (1 << 1)
data['mediumID'] = data['mediumID'].astype('bool')
mEGMSigEff	= 100.*data[(data.isSignal==1) & (data.mediumID==1) & (data.isTrain == 0) & (abs(data.phoSCeta) <= 1.4442)]['bdtWeightF'].sum()/sumSig
mEGMBgEff	= 100.*data[(data.isSignal==0) & (data.mediumID==1) & (data.isTrain == 0) & (abs(data.phoSCeta) <= 1.4442)]['bdtWeightF'].sum()/sumBg
print('EGM:M efficiencies calculated!')

data['tightID'] = data['phoIDbit'] & (1 << 2)
data['tightID'] = data['tightID'].astype('bool')
tEGMSigEff	= 100.*data[(data.isSignal==1) & (data.tightID==1) & (data.isTrain == 0) & (abs(data.phoSCeta) <= 1.4442)]['bdtWeightF'].sum()/sumSig
tEGMBgEff	= 100.*data[(data.isSignal==0) & (data.tightID==1) & (data.isTrain == 0) & (abs(data.phoSCeta) <= 1.4442)]['bdtWeightF'].sum()/sumBg
print('EGM:T efficiencies calculated!')

print('Making plot!')

mpl.use('Agg')
rcParams.update({'figure.autolayout': True})
plt.rc('legend', fontsize='medium')

plt.figure(1)

plt.plot([llBgEff], [llSigEff], marker='H', markersize=4, color='#41ab5d')
plt.text(llBgEff*0.98,llSigEff*1.002,'VL(%1.3f,%1.1f)' %(llBgEff, llSigEff),horizontalalignment='left', va='bottom', fontsize=9)

plt.plot([lBgEff], [lSigEff], marker='H', markersize=4, color='#41ab5d')
plt.text(lBgEff*0.98,lSigEff*1.002,'L(%1.3f,%1.1f)' %(lBgEff, lSigEff),horizontalalignment='left', va='bottom', fontsize=9)

plt.plot([mBgEff], [mSigEff], marker='H', markersize=4, color='#41ab5d')
plt.text(mBgEff*0.98,mSigEff*1.002,'M(%1.3f,%1.1f)' %(mBgEff, mSigEff),horizontalalignment='left', va='bottom', fontsize=9)

plt.plot([tBgEff], [tSigEff], marker='H', markersize=4, color='#41ab5d')
plt.text(tBgEff*0.98,tSigEff*1.002,'T(%1.3f,%1.1f)' %(tBgEff, tSigEff),horizontalalignment='left', va='bottom', fontsize=9)



plt.plot([llBgCBEff], [llSigCBEff], marker="s", markersize=4, color='#2171b5')
plt.text(llBgCBEff*0.96,llSigCBEff*0.985,'CB:VL(%1.3f,%1.1f)' %(llBgCBEff, llSigCBEff),horizontalalignment='left', va='bottom', fontsize=9)

plt.plot([lBgCBEff], [lSigCBEff], marker='s', markersize=4, color='#2171b5')
plt.text(lBgCBEff*0.96,lSigCBEff*0.985,'CB:L(%1.3f,%1.1f)' %(lBgCBEff, lSigCBEff),horizontalalignment='left', va='bottom', fontsize=9)

plt.plot([mBgCBEff], [mSigCBEff], marker='s', markersize=4, color='#2171b5')
plt.text(mBgCBEff*0.96,mSigCBEff*0.985,'CB:M(%1.3f,%1.1f)' %(mBgCBEff, mSigCBEff),horizontalalignment='left', va='bottom', fontsize=9)

plt.plot([tBgCBEff], [tSigCBEff], marker='s', markersize=4, color='#2171b5')
plt.text(tBgCBEff*0.96,tSigCBEff*0.983,'CB:T(%1.3f,%1.1f)' %(tBgCBEff, tSigCBEff),horizontalalignment='left', va='bottom', fontsize=9)



plt.plot([lEGMBgEff], [lEGMSigEff], marker='o', markersize=4, color='#cc4c02')
plt.text(lEGMBgEff*0.982,lEGMSigEff*1.002,'EGM:L(%1.3f,%1.1f)' %(lEGMBgEff, lEGMSigEff),horizontalalignment='right', va='bottom', fontsize=9)

plt.plot([mEGMBgEff], [mEGMSigEff], marker='o', markersize=4, color='#cc4c02')
plt.text(mEGMBgEff*0.98,mEGMSigEff*1.002,'EGM:M(%1.3f,%1.1f)' %(mEGMBgEff, mEGMSigEff),horizontalalignment='left', va='bottom', fontsize=9)

plt.plot([tEGMBgEff], [tEGMSigEff], marker='o', markersize=4, color='#cc4c02')
plt.text(tEGMBgEff*0.98,tEGMSigEff*1.002,'EGM:T(%1.3f,%1.1f)' %(tEGMBgEff, tEGMSigEff),horizontalalignment='left', va='bottom', fontsize=9)

plt.grid(alpha=0.5)
plt.xlabel('Background Efficiency (%)')
plt.ylabel('Signal Efficiency (%)')

plt.savefig('%s/ROC.pdf' % saveDir, dpi=600)
plt.savefig('%s/ROC.png' % saveDir, dpi=600)

print('Done!' )
