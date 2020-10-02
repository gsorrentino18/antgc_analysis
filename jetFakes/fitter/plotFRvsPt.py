#!/usr/local/bin/python3

import numpy as np
import pandas as pd
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.use('Agg')
rcParams.update({'figure.autolayout': True})
plt.rc('legend', fontsize='medium')

dat=pd.read_csv('fits/vsPt/fr.txt') 

dat['pTc'] = (dat['pTmax'] + dat['pTmin'])/2.
dat['pTerr'] = (dat['pTmax'] - dat['pTmin'])/2.


plt.figure(1)

plt.errorbar(x=dat['pTc'].values, y=dat['Fr'].values, yerr=dat['FrErr'].values, xerr=dat['pTerr'].values, color='#d73027',marker="H", markersize=4, markeredgewidth=0)

plt.grid()
plt.xlabel(r"$p_{T}$ (GeV)", fontsize=24)
plt.ylabel('Fake Ratio', fontsize=24)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.savefig('fits/fr_vs_pt.pdf', dpi=600)
plt.savefig('fits/fr_vs_pt.png', dpi=600)