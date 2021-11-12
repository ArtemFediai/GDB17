#!/usr/bin/env python
# coding: utf-8

# # Example of Kmeans for SAMPLE_SIZE molecules and molecules with size MOL_SIZE_OF_INTEREST
# use positional arguments
# In[47]:

from os import uname
import os.path
from rdkit import Chem
from rdkit import RDLogger
import numpy as np
import sys
import matplotlib.pyplot as plt
import seaborn
import pickle
import os
import sklearn
RDLogger.DisableLog('rdApp.*')
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin

# ## PARAMETERS
MOL_SIZE_OF_INTEREST = 15  # number of atoms
SAMPLE_SIZE = 10000  #<-- data set size

# ## POSITIONAL ARGS
if len(sys.argv) > 1:  # 0-th argv is a file name
    MOL_SIZE_OF_INTEREST = sys.argv[1]  # 1-st arg
    print(f'According to positional argument, MOL_SIZE_OF_INTEREST is {MOL_SIZE_OF_INTEREST}')
if len(sys.argv) > 2:
    SAMPLE_SIZE = sys.argv[2]
    print(f'According to positional argument, SAMPLE SIZE is {SAMPLE_SIZE}')

# ## SET CONSTANTS
# In[48]:

# add the folder where you store the data set -->
host_name = uname().nodename
if host_name == 'artem-pc':
    DATA_FOLDER = '/home/artem/dataset_gdb17'
elif host_name == 'int-nano':
    DATA_FOLDER = '/home/ws/bh5670/dataset_gdb17'
else:
    print(f'I  do not know, where DATA_FOLDER is at {host_name}? exiting...')
    exit()
# < -- add the folder where you store the sata set

#DATA_FILE_NAME = 'GDB17.50000000.smi'
DATA_FILE_NAME = str(MOL_SIZE_OF_INTEREST) + '.smi'
DATA_PATH =  DATA_FOLDER + '/' + str(SAMPLE_SIZE) + '/' + DATA_FILE_NAME

_N_CLUSTERS = 100  # we want 100 mols. Not more, not fewer --> 'protected'
METRICS = 'euclidean'

# In[49]:


#estimates exe time
def timeit(fn):
    """
    measures execution time. used as decorator
    @rtype: object
    """
    from time import perf_counter
    from functools import wraps

    @wraps(fn)
    def inner(*args, **kwargs):
        start = perf_counter()
        result = fn(*args, **kwargs)
        end = perf_counter()
        elapsed = end - start

        args_ = [str(a) for a in args]
        kwargs_ = ['{0}={1}'.format(k, v) for (k, v) in kwargs.items()]
        all_args = args_ + kwargs_
        args_str = ','.join(all_args)
        # print('{0}({1}) took {2:.6f}s to run.'.format(fn.__name__,
        #                                               args_str,
        #                                               elapsed))
        print('{0} took {1:.6f}s to run.'.format(fn.__name__, elapsed))
        return result

    return inner


# ## DATA

# In[50]:


DATA_PATH


# ## LOAD DATA
# 

# In[51]:


suppl = Chem.SmilesMolSupplier(DATA_PATH)
smiles_list = [Chem.MolToSmiles(m) for m in suppl]
smiles_list[0:3]
# smiles_list is our raw data. now we need to convert them into some representation


# ## GET FEATURES FOR EVERY MOLECULE
# 

# In[52]:


features = [Chem.RDKFingerprint(m) for m in suppl]


# In[53]:


fp_array = np.array([[int(i) for i in Chem.RDKFingerprint(m).ToBitString()] for m in suppl], dtype=int)


# In[54]:


fp_array


# ## KMEANS

# In[55]:


@timeit
def returns_kmeans(n_clusters, fp_array):
    return KMeans(n_clusters=_N_CLUSTERS).fit(X=fp_array)

kmeans = returns_kmeans(n_clusters=_N_CLUSTERS, fp_array=fp_array)


# In[56]:


kmeans.labels_

kmeans.cluster_centers_


# In[57]:


central_mol_ids = pairwise_distances_argmin(X=kmeans.cluster_centers_, Y=fp_array, metric=METRICS)


# ## SAVE

# In[58]:


# - as a list of mol idx
fname_wo_extension =f'{DATA_FOLDER}/{str(SAMPLE_SIZE)}/{MOL_SIZE_OF_INTEREST}_{METRICS}'
np.savetxt(f'{fname_wo_extension}.txt', central_mol_ids, fmt='%i')
# - as smiles
list_of_central_smiles = [Chem.MolToSmiles(suppl[int(mol_id)]) for mol_id in central_mol_ids]
print(list_of_central_smiles[0:2])
list_of_central_smiles

current_path = f'{fname_wo_extension}.smi'
with Chem.SmilesWriter(current_path) as w:
    for mol_id in central_mol_ids:
        w.write(suppl[int(mol_id)])
print('I am done')


# In[58]:




