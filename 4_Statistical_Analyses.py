#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#import necessary packages
import pandas as pd 
import numpy as np


# In[ ]:


#load in data
df = pd.read_csv('NRF2_CCLE_data_afterML.csv')


# In[ ]:


genes = [*set(df['gene'].to_list())]


# In[ ]:


#calculate log_odds
gene_list = []
cluster_list = []
log_odds_list = []

for x in genes:
    sub = df[df['gene']==x]
    N = len(sub)

    array_clust = sub['cluster'].to_numpy()
    array_class1 = sub['class'].to_numpy()

    joint_counts = pd.crosstab(array_clust,array_class1, rownames=['cluster'], colnames=['class'])
    joint_prob = joint_counts/N

    sub_cluster_list = [*set(sub['cluster'].to_list())]

    for i in sub_cluster_list:
        sub_cluster_df = sub[sub['cluster']==i]

        if len([*set(sub_cluster_df['label'].to_list())]) >= 3:

            if 0 in set(sub['class']) and 1 in set(sub_cluster_df['class']) and 0 in set(sub_cluster_df['class']):
                d = sum(joint_counts[0])-joint_counts[0][i]
                b = joint_counts[0][i]
                c = sum(joint_counts[1])-joint_counts[1][i]
                a = joint_counts[1][i]

                odds = np.true_divide(a,c+0.00000001)
                OR = np.true_divide(np.multiply(a,d),np.multiply(b,c))
                LOR = np.log(OR)

                gene_list.append(x)
                cluster_list.append(i)
                log_odds_list.append(LOR)
            elif 0 in set(sub['class']) and 1 in set(sub_cluster_df['class']):
                d = sum(joint_counts[0])
                b = 0
                c = sum(joint_counts[1])-joint_counts[1][i]
                a = joint_counts[1][i]

                odds = np.true_divide(a,c)
                OR = np.true_divide(np.multiply(a,d),(np.multiply(b,c)+0.001))
                LOR = np.log(OR)

                gene_list.append(x)
                cluster_list.append(i)
                log_odds_list.append(LOR)

            elif 1 in set(sub_cluster_df['class']):
                d = 0
                b = 0
                c = sum(joint_counts[1])-joint_counts[1][i]
                a = joint_counts[1][i]

                odds = np.true_divide(a,c)
                OR = np.true_divide(np.multiply(a,d),(np.multiply(b,c)+0.001))
                LOR = np.log(OR)

                gene_list.append(x)
                cluster_list.append(i)
                log_odds_list.append(LOR)

            else:
                continue

        else:
            continue


# In[ ]:


df2= pd.DataFrame(list(zip(gene_list, cluster_list,log_odds_list)),
                  columns =['gene', 'cluster','log_odds'])


# In[ ]:


#get rid of outliers
df3 = df2[df2['log_odds']<100]
df3 = df3[df3['log_odds']>-100]
df3 = df3.dropna()
df3


# In[ ]:


#work on getting size of each cluster next
#assign cluster id to each df
df3['id'] = df3['gene'] + '_' + df3['cluster'].astype(str)
df['id'] = df['gene'] + '_' + df['cluster'].astype(str)


# In[ ]:


cluster_id = df3['id'].to_list()


# In[ ]:


#get length of each cluster
id_list = []
num_list = []

for i in cluster_id:
    sub = df[df['id'] == i]
    id_list.append(i)
    num_list.append(len(sub))


# In[ ]:


df_num= pd.DataFrame(list(zip(id_list, num_list)),
                  columns =['id', 'num'])
df_num = df_num.set_index('id')


# In[ ]:


df3['num'] = 'none'
df3 = df3.reset_index()


# In[ ]:


#add lengths to df
for i in range(len(df3["id"])):
    if df3["id"][i] in df_num.index:
        df3['num'][i] = df_num.loc[df3["id"][i]][0]


# In[ ]:


#get norm diff of class1 to class0 for each cluster and gene
id_list2 = []
diff_list = []

for i in cluster_id:
    sub = df[df['id'] == i]
    sub1 = sub[sub['class']==1]
    sub0 = sub[sub['class']==0]
    id_list2.append(i)
    diff_list.append(sub1['label'].mean()-sub0['label'].mean())


# In[ ]:


df_diff= pd.DataFrame(list(zip(id_list2, diff_list)),
                  columns =['id', 'diff'])


# In[ ]:


df_diff = df_diff.set_index('id')


# In[ ]:


df3['diff'] = 'none'


# In[ ]:


for i in range(len(df3["id"])):
    if df3["id"][i] in df_diff.index:
        df3['diff'][i] = df_diff.loc[df3["id"][i]][0]


# In[ ]:


print(df['label'].mean(), df['label'].std())


# In[ ]:


df3['norm_diff'] = 'none'


# In[ ]:


for i in range(len(df3)):
    df3['norm_diff'][i] = (df3['diff'][i] - 2279.635620930903)/1907.116650152018


# In[ ]:


df3


# In[ ]:


df3.to_csv('filtered_sorted_NRF2_CCLE_RNAseq.csv')
#can now filter and sort this dataframe to find genes with high density, log_odds, norm_diff


# In[ ]:


df['filter'] = 'no'


# In[ ]:


clusters_list = df4['id'].to_list()


# In[ ]:


df4 = df3[df3['num'] > 12]
df4 = df4[df4['norm_diff'] > 1]


# In[ ]:


df4


# In[ ]:


for i in range(len(df['id'])):
    if df["id"][i] in clusters_list:
        df['filter'][i] = 'yes'


# In[ ]:


df_filtered = df[df['filter']=='yes']


# In[ ]:


df_filtered.to_csv('NRF2_CCLE_RNAseq_postfilter.csv')


# In[ ]:




