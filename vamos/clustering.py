    import logging

    logger = logging.getLogger(__name__)

    def run_clustering(**kwargs):
        """
        Auto-wrapped pipeline step. The original notebook code has been
        converted to a callable function. You can pass configuration via kwargs.
        """
        #!/usr/bin/env python
# coding: utf-8

# In[ ]:


#import all necessary packages
import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from collections import Counter
import seaborn as sns
import statistics as st

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import ListedColormap
from sklearn.metrics.pairwise import pairwise_distances
import sklearn

import itertools
import os

from kneed import KneeLocator, DataGenerator as dg

import numpy


# In[ ]:


pd.set_option('display.max_columns', None)


# In[ ]:


#read in data
df1 = pd.read_csv('alphafold_mapped_confidence_above70.csv')


# In[ ]:


#first do one gene, here VPS13D
df2 = df1[df1['gene_name']=='VPS13D']


# In[ ]:


df = pd.DataFrame().assign(x =df2['x_coord'], y =df2['y_coord'], z =df2['z_coord'])

nbrs = NearestNeighbors(n_neighbors=5).fit(df)
neigh_dist, neigh_ind = nbrs.kneighbors(df)
sort_neigh_dist = np.sort(neigh_dist, axis=0)
k_dist = sort_neigh_dist[:, 4]
k_dist_list = k_dist.tolist()
length = list(range(len(k_dist_list)))


x, y = length, k_dist_list
kl = KneeLocator(x, y, curve="convex")
epsilon = int(k_dist_list[kl.knee])/2

clusters = DBSCAN(eps= epsilon, min_samples=6).fit(df)
cluster_list = clusters.labels_.tolist()
df.insert(0, 'cluster', cluster_list)
df['gene'] = 'VPS13d'
df = df.reset_index()


# In[ ]:


cluster_list2 = [*set(cluster_list)]
cluster_list2


# In[ ]:


# do one cluster

df_cluster_loop = df[df['cluster'] == 0]
dist_matrix = pd.DataFrame().assign(x=df_cluster_loop['x'], y=df_cluster_loop['y'], z=df_cluster_loop['z'])
pw_dist = pairwise_distances(dist_matrix)
avg_dist = numpy.average(pw_dist)
df_cluster_loop["dist"] = avg_dist
df_cluster_loop['num_members'] = len(df_cluster_loop['gene'])
df_cluster = df_cluster_loop.copy()


# In[ ]:


#do the rest
cluster_list2.remove(0)

for i in cluster_list2:
    df_cluster_loop = df[df['cluster'] == i]
    dist_matrix = pd.DataFrame().assign(x=df_cluster_loop['x'], y=df_cluster_loop['y'], z=df_cluster_loop['z'])
    pw_dist = pairwise_distances(dist_matrix)
    avg_dist = numpy.average(pw_dist)
    df_cluster_loop["dist"] = avg_dist
    df_cluster_loop['num_members'] = len(df_cluster_loop['gene'])

    dfs = [df_cluster, df_cluster_loop]
    df_cluster = pd.concat(dfs)


# In[ ]:


#get all genes that have 5 or more mutations
genes = df1['gene_name'].to_list()
genes.remove('VPS13D')

res  = []
for x in set(genes):
    if genes.count(x) >= 5:
        res.append(x)


# In[ ]:


#do all genes 
#possible error: small number of genes may not have a clear "knee", epsilon for those has to be picked manually
for i in res:
    subset = df1[df1["gene_name"] == i]
    xyz = pd.DataFrame().assign(x =subset['x_coord'], y =subset['y_coord'], z =subset['z_coord'])
    nbrs = NearestNeighbors(n_neighbors=5).fit(xyz)
    neigh_dist, neigh_ind = nbrs.kneighbors(xyz)
    sort_neigh_dist = np.sort(neigh_dist, axis=0)
    k_dist = sort_neigh_dist[:, 4]
    k_dist_list = k_dist.tolist()
    length = list(range(len(k_dist_list)))

    x, y = length, k_dist_list
    kl = KneeLocator(x, y, curve="convex")
    epsilon = int(k_dist_list[kl.knee])/2

    clusters = DBSCAN(eps= epsilon, min_samples=6).fit(xyz)
    cluster_list = clusters.labels_.tolist()
    xyz.insert(0, 'cluster', cluster_list)
    xyz['gene'] = i
    xyz = xyz.reset_index()

    cluster_list2 = [*set(cluster_list)]

    for x in cluster_list2:
        df_cluster_loop = xyz[xyz['cluster'] == x]
        dist_matrix = pd.DataFrame().assign(x=df_cluster_loop['x'], y=df_cluster_loop['y'], z=df_cluster_loop['z'])
        pw_dist = pairwise_distances(dist_matrix)
        avg_dist = numpy.average(pw_dist)
        df_cluster_loop["dist"] = avg_dist
        df_cluster_loop['num_members'] = len(df_cluster_loop['gene'])

        dfs = [df_cluster, df_cluster_loop]
        df_cluster = pd.concat(dfs)


# In[ ]:


#prep og df to merge with cluster info
df1 = df1.rename(columns={"x_coord": "x", "y_coord": "y","z_coord":'z'})
df1['index1'] = df1.index
df1 = df1.rename(columns={"gene_name": "gene"})


# In[ ]:


#prep cluster info df to merge with og df
df_cluster = df_cluster.rename(columns={"index": "index1"})


# In[ ]:


#merge both dfs to get all info
df_final = pd.merge(df1,df_cluster,on=['x','y','z','gene','index1'])


# In[ ]:


#save
df_final.to_csv('ready_for_ML_processed.csv')

        logger.info("Completed run_clustering.")
        return True
