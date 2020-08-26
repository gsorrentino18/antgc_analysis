#!/usr/bin/env python3

####################################################################
# Photon ID BDT for aNTGC search in Z(->nu nu) + Gamma channel     #
# Mohammad Abrar Wadud (BD), 05/06/2020                            #
####################################################################

import json
import matplotlib.pyplot as plt

f =open('/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/BDT/data/optimizedV1/validationProgress.log')
data = json.load(f) 

trainAUC=data['train']['auc']
testAUC=data['test']['auc']
validationAUC=data['validation']['auc']

plt.plot(trainAUC, label='Train', color='#034e7b')
plt.plot(testAUC, label='Test', color='#d73027')
plt.plot(validationAUC, label='Validation', color='#238443')

plt.legend(loc="lower right")
plt.grid()
plt.xlabel('Number of Boosting Rounds')
plt.ylabel('AUC')

plt.savefig('AUCvsTrees.png', dpi=600)
