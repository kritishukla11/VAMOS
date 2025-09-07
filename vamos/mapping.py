    import logging

    logger = logging.getLogger(__name__)

    def run_mapping(**kwargs):
        """
        Auto-wrapped pipeline step. The original notebook code has been
        converted to a callable function. You can pass configuration via kwargs.
        """
        #!/usr/bin/env python
# coding: utf-8

# In[ ]:


#importing all necessary packages

#this is for mapping
from unipressed import IdMappingClient
import time

#pandas
import pandas as pd

#cif reading
import sys
from gemmi import cif

#plotting
import matplotlib.pyplot as plt

#joining lists, opening files
import itertools
import os


# In[ ]:


pd.set_option('display.max_columns', None)


# In[ ]:


#this is the data from CCLE
depmap = pd.read_csv("CCLE_mutations.csv")
depmap = depmap.rename(columns={"Hugo_Symbol": "gene_name"})


# In[ ]:


#this is a file with the gene name to uniprot mapping info
info = pd.read_csv('mut_confidence_mapped.csv')


# In[ ]:


#merging the two files and processing the columns
metadata = pd.merge(depmap, info, on='gene_name')
metadata = metadata.loc[metadata['Variant_Type'] == 'SNP']
metadata = metadata.reset_index()
metadata[['c.','info']] = metadata['cDNA_Change'].str.split('.',expand=True)
metadata.drop('c.', inplace=True, axis=1)
metadata[['mut','change']] = metadata['info'].str.split('>',expand=True)
metadata.drop('change', inplace=True, axis=1)
metadata.drop('info', inplace=True, axis=1)
metadata[['atom_mut', 'mut_id']] = metadata['mut'].str.extract('(\d+\.?\d*)([A-Za-z]*)', expand = True)
metadata.drop('mut', inplace=True, axis=1)
metadata.drop('mut_id', inplace=True, axis=1)
metadata['atom_mut'] = pd.to_numeric(metadata['atom_mut'])


# In[ ]:


#do one one gene first, here VPS13D
VPS = metadata[metadata['gene_name'] == 'VPS13D']
VPS


# In[ ]:


#map to Alphafold 
doc = cif.read('AF-Q5THJ4-F1-model_v3.cif')
block = doc.sole_block()
atom_id = block.find_loop('_atom_site.id')
x_coord = block.find_loop('_atom_site.Cartn_x')
y_coord = block.find_loop('_atom_site.Cartn_y')
z_coord = block.find_loop('_atom_site.Cartn_z')

atom_id = [eval(i) for i in atom_id]
atom_id = pd.to_numeric(atom_id)
x_coord = [eval(i) for i in x_coord] 
y_coord = [eval(i) for i in y_coord] 
z_coord = [eval(i) for i in z_coord] 

xyz_coord_dic = {'atom_mut':atom_id,'x_coord':x_coord, 'y_coord':y_coord, 'z_coord':z_coord}
xyz_coord_df = pd.DataFrame(xyz_coord_dic)


# In[ ]:


#merge xyz with gene info
tot = pd.merge(xyz_coord_df, VPS, on='atom_mut')


# In[ ]:


#get all genes
gene_list = metadata['gene_name'].to_list()
gene_list = [*set(gene_list)]


# In[ ]:


#map all genes
def bar(s):
    return 'AF-' + s + '-F1-model_v3.cif'

problem_genes = []

for x in gene_list:
    sub = metadata[metadata['gene_name'] == x]
    sub = sub.reset_index()
    y = sub['uniprot_id'][0]

    if os.path.isfile(bar(y)):
        with open(bar(y), 'r') as f1:

            doc = cif.read('AF-Q5THJ4-F1-model_v3.cif')
            block = doc.sole_block()
            atom_id = block.find_loop('_atom_site.id')
            x_coord = block.find_loop('_atom_site.Cartn_x')
            y_coord = block.find_loop('_atom_site.Cartn_y')
            z_coord = block.find_loop('_atom_site.Cartn_z')

            atom_id = [eval(i) for i in atom_id]
            atom_id = pd.to_numeric(atom_id)
            x_coord = [eval(i) for i in x_coord] 
            y_coord = [eval(i) for i in y_coord] 
            z_coord = [eval(i) for i in z_coord] 

            xyz_coord_dic = {'atom_mut':atom_id,'x_coord':x_coord, 'y_coord':y_coord, 'z_coord':z_coord}
            xyz_coord_df = pd.DataFrame(xyz_coord_dic)

            df_new = pd.merge(xyz_coord_df, sub, on='atom_mut')

            dfs = [tot, df_new]

            tot = pd.concat(dfs)

    else:
        problem_genes.append(x)
        continue


# In[ ]:


tot = tot.reset_index(drop=True)


# In[ ]:


#save this df

tot.to_csv('alphafold_mapped.csv')


# In[ ]:


#add structural info, first for VPS13D

doc = cif.read('AF-Q5THJ4-F1-model_v3.cif')
block = doc.sole_block()
pos1 = block.find_loop('_struct_conf.beg_auth_seq_id')
struc1 = block.find_loop('_struct_conf.conf_type_id')
pos = [eval(i) for i in pos1]
pos = pd.to_numeric(pos)
struc = [i for i in struc1]
struc_dic = {'position': pos,'structure':struc}
struc_df = pd.DataFrame(struc_dic)
struc_df


# In[ ]:


pos2 = block.find_loop('_struct_conf.end_auth_seq_id')
struc2 = block.find_loop('_struct_conf.id')
pos_ = [eval(i) for i in pos2]
pos_ = pd.to_numeric(pos_)
struc_ = [i for i in struc2]
struc2_dic = {'position': pos_,'structure':struc_}
struc2_df = pd.DataFrame(struc2_dic)
struc2_df


# In[ ]:


dfs = [struc_df,struc2_df]
df = pd.concat(dfs)
df = df.reset_index(drop = True)
df['structure_pos'] = df['structure'].str.replace('([A-Z]+)', '')
df['structure_type'] = df['structure'].str.extract('([A-Z]+)')
df = df.drop(columns=['structure_pos'])
df = df.drop(columns=['structure'])
df['disorder'] = df['structure_type'].map(lambda x: x == "STRN")
df['position'] = df['position'].astype(int)
df


# In[ ]:


df['position'] = df['position'].astype(int)


# In[ ]:


tot['Protein_Change'] = tot['Protein_Change'].str.replace('p.', '')
tot[['original', 'position', 'changed_mut']] = tot['Protein_Change'].str.extract('([A-Za-z]+)(\d+\.?\d*)([A-Za-z]*)', expand = True)


# In[ ]:


tot['position'] = pd.to_numeric(tot['position'])
tot = tot.dropna(subset=['position'])
tot['position'] = tot['position'].astype(int)


# In[ ]:


sub = tot[tot["gene_name"] == 'VPS13D']


# In[ ]:


merged = pd.merge(df, sub, on=['position'])


# In[ ]:


#do it for all genes now
genes = tot['gene_name'].tolist()
genes = [*set(genes)]


# In[ ]:


genes.remove('VPS13D')


# In[ ]:


def bar(s):
    return 'AF-' + s + '-F1-model_v3.cif'

for i in genes:
    sub = tot[tot['gene_name'] == i]
    sub = sub.reset_index(drop=True)

    y = sub['uniprot_id'][0]

    if os.path.isfile(bar(y)):
        with open(bar(y), 'r') as f1:

            doc = cif.read(bar(y))
            block = doc.sole_block()
            pos1 = block.find_loop('_struct_conf.beg_auth_seq_id')
            struc1 = block.find_loop('_struct_conf.conf_type_id')
            pos = [eval(i) for i in pos1]
            pos = pd.to_numeric(pos)
            struc = [i for i in struc1]
            struc_dic = {'position': pos,'structure':struc}
            struc_df = pd.DataFrame(struc_dic)

            pos2 = block.find_loop('_struct_conf.end_auth_seq_id')
            struc2 = block.find_loop('_struct_conf.id')
            pos_ = [eval(i) for i in pos2]
            pos_ = pd.to_numeric(pos_)
            struc_ = [i for i in struc2]
            struc2_dic = {'position': pos_,'structure':struc_}
            struc2_df = pd.DataFrame(struc2_dic)

            dfs = [struc_df,struc2_df]
            structure_df = pd.concat(dfs)

            structure_df['structure_pos'] = structure_df['structure'].str.replace('([A-Z]+)', '')
            structure_df['structure_type'] = structure_df['structure'].str.extract('([A-Z]+)')
            structure_df = structure_df.drop(columns=['structure_pos'])
            structure_df = structure_df.drop(columns=['structure'])
            structure_df['disorder'] = structure_df['structure_type'].map(lambda x: x == "STRN")
            structure_df['position'] = structure_df['position'].astype(int)

            merged2 = pd.merge(structure_df, sub, on=['position'])

            dfs = [merged2,merged]
            merged = pd.concat(dfs)


# In[ ]:


#save
merged.to_csv('alphafold_mapped_struc.csv')


# In[ ]:


#get rid of disordered regions
xyz = merged[merged['disorder'] == False]


# In[ ]:


#formatting dataframe to add confidence scores by isolating protein residue
xyz[['firstcol', 'mut_residue', 'thirdcol']] = xyz['Protein_Change'].str.extract('([A-Za-z]+)(\d+\.?\d*)([A-Za-z]*)', expand = True)
xyz = xyz.drop(columns=['firstcol','thirdcol'])
xyz = xyz.dropna(subset=['mut_residue'])
xyz = xyz.reset_index(drop=True)


# In[ ]:


#add scores
def bar(s):
    return 'AF-' + s + '-F1-model_v3.cif'

scores_list = []

try:
    for i in range(len(xyz)):
        x = merged['uniprot_id'][i] 
        if os.path.isfile(bar(x)):
            with open(bar(x), 'r') as f1:
                doc = cif.read(bar(x))
                block = doc.sole_block()
                residue = block.find_loop('_ma_qa_metric_local.label_seq_id')
                score = block.find_loop('_ma_qa_metric_local.metric_value')
                scorelist = list(score)
                residuelist = list(residue)
                residuelist = [eval(i) for i in residuelist]
                scorelist = list(score)
                y = int(xyz['mut_residue'][i])
                if y in residuelist:
                    scores_list.append(scorelist[y-1])
                else:
                    scores_list.append('no') 

except:
    print(i)


# In[ ]:


xyz['score'] = scores_list


# In[ ]:


xyz.to_csv('alphafold_nodisorder_confidence_mapped.csv')


# In[ ]:


#can do the same for the whole mapped dataset without structural info too

def bar(s):
    return 'AF-' + s + '-F1-model_v3.cif'

scores_list = []

try:
    for i in range(len(tot)):
        x = tot['uniprot_id'][i] 
        if os.path.isfile(bar(x)):
            with open(bar(x), 'r') as f1:
                doc = cif.read(bar(x))
                block = doc.sole_block()
                residue = block.find_loop('_ma_qa_metric_local.label_seq_id')
                score = block.find_loop('_ma_qa_metric_local.metric_value')
                scorelist = list(score)
                residuelist = list(residue)
                residuelist = [eval(i) for i in residuelist]
                scorelist = list(score)
                y = int(tot['position'][i])
                if y in residuelist:
                    scores_list.append(scorelist[y-1])
                else:
                    scores_list.append('no') 

except:
    print(i)


# In[ ]:


#add to df
tot['score'] = scores_list


# In[ ]:


#get rid of anything with no info
tot = tot[tot['score']!= 'no']
tot = tot.reset_index(drop=True)


# In[ ]:


#save
tot.to_csv('alphafold_mapped_confidence.csv')


# In[ ]:


#can filter for confidence level, say above 70, if needed

xyz2 = xyz[xyz['score'].astype(float)>70]


# In[ ]:


xyz2.to_csv('alphafold_mapped_confidence_above70.csv')

        logger.info("Completed run_mapping.")
        return True
