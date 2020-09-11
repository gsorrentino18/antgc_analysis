#!/usr/bin/env python3


# configure options
from argparse import ArgumentParser
parser = ArgumentParser(description='Photon ID BDT parameter tuning for CMS aNTGC search in Z(->nu nu) + gamma channel')
parser.add_argument('--modelFilePath', type=str, help='BDT model file',
					default='training/aNTGC_photon_BDT_2020_09_10_23_04_13.pkl', action='store')
parser.add_argument('--saveDir', type=str, default='training/',
					help='Save directory', action='store')

args = parser.parse_args()

from os import system
system('mkdir -p %s' % args.saveDir)
import pandas as pd
import xgboost as xg
import pickle

pd.set_option("display.max_rows", None, "display.max_columns", None)

print("modelFilePath\t=\t" + args.modelFilePath)
print("saveDir\t\t=\t" + args.saveDir)

phoModel = pickle.load(open(args.modelFilePath, "rb"))

fScores = pd.DataFrame.from_dict(phoModel.get_score(importance_type='weight'), orient='index', columns=['weight'])
fScores.index.name = 'Feature'
fScores['gain'] = fScores.index.map(phoModel.get_score(importance_type='gain'))
fScores['cover'] = fScores.index.map(phoModel.get_score(importance_type='cover'))
fScores['total_gain'] = fScores.index.map(phoModel.get_score(importance_type='total_gain'))
fScores['total_cover'] = fScores.index.map(phoModel.get_score(importance_type='total_cover'))

fScores=fScores.sort_values(by=['gain'], ascending=False)

featImpDumpFilePath=args.saveDir + 'featImps.txt'
fScores.to_csv(featImpDumpFilePath, index=True)

print("\nFeature importances (" + str(len(fScores.index)) + " features) saved to : %s " % featImpDumpFilePath)
print(fScores)