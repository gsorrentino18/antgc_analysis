def aucW(classes, predictions, weights):
	"""
	Calculating ROC AUC score as the probability of correct ordering
	"""

	if weights is None:
		weights = np.ones_like(predictions)

	assert len(classes) == len(predictions) == len(weights)
	assert classes.ndim == predictions.ndim == weights.ndim == 1
	class0, class1 = sorted(np.unique(classes))

	data = np.empty(
		shape=len(classes),
		dtype=[('c', classes.dtype),
			   ('p', predictions.dtype),
			   ('w', weights.dtype)]
	)
	data['c'], data['p'], data['w'] = classes, predictions, weights

	data = data[np.argsort(data['c'])]
	# here we're relying on stability as we need class orders preserved
	data = data[np.argsort(data['p'], kind='mergesort')]

	correction = 0.
	# mask1 - bool mask to highlight collision areas
	# mask2 - bool mask with collision areas' start points
	mask1 = np.empty(len(data), dtype=bool)
	mask2 = np.empty(len(data), dtype=bool)
	mask1[0] = mask2[-1] = False
	mask1[1:] = data['p'][1:] == data['p'][:-1]
	if mask1.any():
		mask2[:-1] = ~mask1[:-1] & mask1[1:]
		mask1[:-1] |= mask1[1:]
		ids, = mask2.nonzero()
		correction = sum([((dsplit['c'] == class0) * dsplit['w'] * msplit).sum() *
						  ((dsplit['c'] == class1) *
						   dsplit['w'] * msplit).sum()
						  for dsplit, msplit in zip(np.split(data, ids), np.split(mask1, ids))]) * 0.5

	weights_0 = data['w'] * (data['c'] == class0)
	weights_1 = data['w'] * (data['c'] == class1)
	cumsum_0 = weights_0.cumsum()

	return ((cumsum_0 * weights_1).sum() - correction) / (weights_1.sum() * cumsum_0[-1])

import numpy as np
import pandas as pd

def sampleStats(df):

	sigCount = len(df[df['isSignal'] == 1].index)
	bkgCount = len(df[df['isSignal'] == 0].index)
	print('Total Events\t=\t%d (B/S= %d / %d)' %
		  ((sigCount + bkgCount), bkgCount, sigCount))


	# print('Total Events\t=\t%d (B/S= %d / %d = %1.2f)' %
	# 	  ((sigCount + bkgCount), bkgCount, sigCount, float(bkgCount) / float(sigCount)))

	# sigSumW = np.sum(df[df.isSignal == 1]['xSecW'].values)
	# bkgSumW = np.sum(df[df.isSignal == 0]['xSecW'].values)
	# BoS = bkgSumW / sigSumW
	# print('BackgroundSumW/SignalSumW \t=\t%.6f / %.6f\t=\t%.6f' %
	# 	  (bkgSumW, sigSumW, BoS))

	# sigSumPtEtaW = np.sum(df[df.isSignal == 1]['bdtWeight'].values)
	# bkgSumPtEtaW = np.sum(df[df.isSignal == 0]['bdtWeight'].values)
	# BoSPtEtaW = bkgSumPtEtaW / sigSumPtEtaW
	# print('BackgroundSumPtEtaW/SignalSumPtEtaW \t=\t%.6f / %.6f\t=\t%.6f' %
	# 	  (bkgSumPtEtaW, sigSumPtEtaW, BoSPtEtaW))
