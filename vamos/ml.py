    import logging

    logger = logging.getLogger(__name__)

    def run_ml(**kwargs):
        """
        Auto-wrapped pipeline step. The original notebook code has been
        converted to a callable function. You can pass configuration via kwargs.
        """
        #!/usr/bin/env python
# coding: utf-8

# In[ ]:


#import all packages
import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE

from sklearn.neighbors import NearestNeighbors
from matplotlib import pyplot as plt
import xgboost as xg
import seaborn as sns

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import VotingClassifier, RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
from sklearn.utils.class_weight import compute_class_weight, compute_sample_weight

from collections import defaultdict

from feature_engine.wrappers import SklearnTransformerWrapper

import warnings
import pandas as pd
from sklearn.cluster import DBSCAN
import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE

from sklearn.neighbors import NearestNeighbors
from matplotlib import pyplot as plt
import xgboost as xg
import seaborn as sns

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import VotingClassifier, RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
from sklearn.utils.class_weight import compute_class_weight, compute_sample_weight
from sklearn.model_selection import cross_validate
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import make_scorer, accuracy_score, precision_score, recall_score, f1_score
from sklearn.metrics import plot_roc_curve

from collections import defaultdict


# In[ ]:


pd.set_option('display.max_columns', None)


# In[ ]:


#read in data
new_data = pd.read_csv('ready_for_ML_processed.csv')


# In[ ]:


#need to add tumor sample barcodes from depmap
#read in file with tumor sample barcode info
df1 = pd.read_csv('sample_info.csv')
df1 = df1.set_index(['DepMap_ID'])


# In[ ]:


#add blank column to df
new_data['Tumor_Sample_Barcode'] = 'none'


# In[ ]:


#add barcode info to df
for i in range(len(new_data["DepMap_ID"])):
    if new_data["DepMap_ID"][i] in df1.index:
        new_data['Tumor_Sample_Barcode'][i] = df1.loc[new_data["DepMap_ID"][i]][0]


# In[ ]:


#load gsea expression or crispr sensitivity data  
data_label = pd.read_csv('DF_NRF2v2_GSEA.csv')
data_label = data_label.set_index('index')


# In[ ]:


#create empty column in df
new_data['label'] = 0


# In[ ]:


#add expression/sensitivity data to df
for i in range(len(new_data["Tumor_Sample_Barcode"])):
    if new_data["Tumor_Sample_Barcode"][i] in data_label.index:
        new_data['label'][i] = data_label.loc[new_data["Tumor_Sample_Barcode"][i]][0]


# In[ ]:


new_data['label'].describe()


# In[ ]:


#change SD and mean values based on earlier output
SD =1907.116650
Mean =2279.635621
up_thresh = Mean + SD
down_thresh = Mean - SD


# In[ ]:


#make function to sep class1 from class 0
def categorize(df, threshold):
    row = df['label']
    df['class'] = np.where( row > threshold, 1, 0)


# In[ ]:


#assign preliminary class1 or class0
categorize(new_data, up_thresh)


# In[ ]:


#balance data as needed by oversampling
class1 = new_data[new_data['class']==1]
class0 = new_data[new_data['class']==0]

print(len(class1),len(class0))


# In[ ]:


class1_balanced = pd.concat([class1,class1]) #repeat this line as many times as needed to have equal numbers
new_data_b = pd.concat([class1_balanced,class0])


# In[ ]:


#perform machine learning
X = new_data_b[['x', 'y', 'z','cluster','dist','num_members']]
y = new_data_b['class']

rf_r = RandomForestClassifier(n_estimators=300, random_state=42, class_weight= 'balanced')

scoring_choice = {'accuracy' : make_scorer(accuracy_score), 
                          'precision' : make_scorer(precision_score),
                          'recall' : make_scorer(recall_score), 
                          'f1_score' : make_scorer(f1_score),
                          'mcc:' : make_scorer(matthews_corrcoef)}

rf_r_scores = cross_validate(estimator = rf_r, X=X, y=y, cv=5,
                                      scoring = scoring_choice,
                                      return_train_score=True)
rf_r_results = cross_val_predict(rf_r, X, y, cv=5)


# In[ ]:


#get scores
rf_r_scores


# In[ ]:


#merge machine learning output with previous df
df_results = X.copy()
df_results['class'] = pd.Series(rf_r_results)
new_data_2 = pd.merge(df_results,new_data_b,on=['x','y','z'])


# In[ ]:


new_data_2.to_csv('NRF2_CCLE_data_afterML.csv')

        logger.info("Completed run_ml.")
        return True
