#!/usr/bin/env python3

####################################################################
# Photon ID BDT for aNTGC search in Z(->nu nu) + Gamma channel     #
# Mohammad Abrar Wadud (BD), 05/06/2020                            #
####################################################################

import json
import matplotlib.pyplot as plt

saveDir='dataV2/trainingV2/'
f =open('dataV2/trainingV2/trainingScoreLog_2020_09_12_19_08_19_f.log')
data = json.load(f) 

# data = [json.loads(line) for line in f]

trainAUC=data['train']['auc']
testAUC=data['test']['auc']
validationAUC=data['validation']['auc']


plt.plot(trainAUC, label='Train', color='#034e7b')
plt.plot(testAUC, label='Test', color='#d73027')
plt.plot(validationAUC, label='Validation', color='#238443')

plt.plot(0,trainAUC[0],'ro', color='#034e7b', markersize=2)
plt.plot(0,testAUC[0],'ro', color='#d73027', markersize=2)
plt.plot(0,validationAUC[0],'ro', color='#238443', markersize=2)

plt.legend(loc="lower right")
plt.grid()
plt.xlabel('Number of Boosting Rounds')
plt.ylabel('AUC')

plt.xlim(left=-50)

plt.savefig('%s AUCvsTrees.png' % saveDir, dpi=600)
plt.savefig('%s AUCvsTrees.pdf' % saveDir, dpi=600)

yMin=min([trainAUC[0], testAUC[0], validationAUC[0]])
plt.ylim(bottom=0.999*yMin)
plt.yscale(value='log')



plt.savefig('%s AUCvsTrees_logY.png' % saveDir, dpi=600)
plt.savefig('%s AUCvsTrees_logY.pdf' % saveDir, dpi=600)
