import numpy as np
import scipy.cluster.hierarchy as sch
import pandas as pd 


df = pd.read_csv("plots/1_featCorrelations.txt") 
df=df.drop(['Unnamed: 0'], axis=1)

cluster_th = 3

X = df.values
d = sch.distance.pdist(X)
L = sch.linkage(d, method='complete')
ind = sch.fcluster(L, 0.5*d.max(), 'distance')

columns = [df.columns.tolist()[i] for i in list(np.argsort(ind))]
df = df.reindex_axis(columns, axis=1)

unique, counts = np.unique(ind, return_counts=True)
counts = dict(zip(unique, counts))

i = 0
j = 0
columns = []
for cluster_l1 in set(sorted(ind)):
	j += counts[cluster_l1]
	sub = df[df.columns.values[i:j]]
	if counts[cluster_l1]>cluster_th:        
		X = sub.corr().values
		d = sch.distance.pdist(X)
		L = sch.linkage(d, method='complete')
		ind = sch.fcluster(L, 0.5*d.max(), 'distance')
		col = [sub.columns.tolist()[i] for i in list((np.argsort(ind)))]
		sub = sub.reindex_axis(col, axis=1)
	cols = sub.columns.tolist()
	columns.extend(cols)
	i = j
df = df.reindex_axis(columns, axis=1) 

print(df.columns)
