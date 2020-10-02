#!/usr/local/bin/python3

import numpy as np
import pandas as pd
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.use('Agg')
rcParams.update({'figure.autolayout': True})
plt.rc('legend', fontsize='medium')

dat=pd.read_csv('fits/vsEta/fr.txt') 

dat['etaC'] = (dat['etaMin'] + dat['etaMax'])/2.
dat['etaErr'] = (dat['etaMin'] - dat['etaMax'])/2.


plt.figure(1)

plt.errorbar(x=dat['etaC'].values, y=dat['Fr'].values, yerr=dat['FrErr'].values, xerr=dat['etaErr'].values, color='#d73027',marker="H", markersize=4, markeredgewidth=0)

plt.grid()
plt.xlabel(r"$|\eta|$", fontsize=24)
plt.ylabel('Fake Ratio', fontsize=24)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.savefig('fits/fr_vs_eta.pdf', dpi=600)
plt.savefig('fits/fr_vs_eta.png', dpi=600)