#!/usr/bin/env python3

import matplotlib as mpl
from matplotlib import rcParams
from matplotlib import rc
from xml.dom import minidom    
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np

rcParams.update({'figure.autolayout': True, 'lines.markersize': 2, 'axes.axisbelow': True})
plt.rc('legend', fontsize='medium')

def getROCpts(xmlFile):
	xmldoc = minidom.parse(xmlFile)
	binlist = xmldoc.getElementsByTagName('Bin')

	sEff=[]
	bEff=[]

	for iBin in binlist:

		sEff.append(100.*float(iBin.attributes['effS'].value))
		bEff.append(100.*float(iBin.attributes['effB'].value))

	return bEff, sEff

saveDir='/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/cutOptimize/trainingV2/plots/'

plt.figure(1)

plt.xticks(np.arange(-10, 120, step=5))
plt.yticks(np.arange(-10, 120, step=5))

plt.grid(which='major', linestyle=':', linewidth='0.4', color='#4d4d4d')

plt.xlabel('Background Efficiency (%)')
plt.ylabel('Signal Efficiency (%)')

plt.ylim(-1., 102)
plt.xlim(-1, 50)

bdtCutsFile='/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/cutOptimize/trainingV2/dataset/weights/optimizeCuts_CutsGA.weights.xml' 
bkgEff, sigEff = getROCpts(bdtCutsFile)
plt.scatter(bkgEff, sigEff, label='BDT + H/E + Iso ID', color='#034e7b')

sietaCutsFile='/local/cms/user/wadud/aNTGCmet/aNTGC_analysis/phoIDstudy/cutOptimize/trainingV2/dataset_sietaieta/weights/optimizeCuts_CutsGA.weights.xml'
bkgEff, sigEff = getROCpts(sietaCutsFile)
plt.scatter(bkgEff, sigEff, label='Cut-based ID', color='#d73027')

plt.legend(loc="lower right")

plt.savefig('%s/ROC.pdf' % saveDir, dpi=600)
plt.savefig('%s/ROC.png' % saveDir, dpi=600)
