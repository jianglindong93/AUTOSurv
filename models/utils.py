#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import math
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
import lifelines
from lifelines.utils import concordance_index

dtype = torch.FloatTensor


# In[ ]:


"""
Dataloaders

Note:
While applying the following data loaders, Please pay attention to the difference in covariate names between the two cancer types.
"""

def sort_data(path):
    data = pd.read_csv(path)
    data.sort_values("OS.time", ascending = False, inplace = True)
    patient_id = data.loc[:, ["patient_id"]]
    patient_id.index = range(0, patient_id.shape[0], 1)
    x = data.drop(["patient_id", "OS", "OS.time", "age", "race_white", "stage_i", "stage_ii"], axis = 1).values
    ytime = data.loc[:, ["OS.time"]].values
    yevent = data.loc[:, ["OS"]].values
    age = data.loc[:, ["age"]].values
    stage_i = data.loc[:, ["stage_i"]].values
    stage_ii = data.loc[:, ["stage_ii"]].values
    race_white = data.loc[:, ["race_white"]].values
    return(patient_id, x, ytime, yevent, age, stage_i, stage_ii, race_white)

def load_data(path, dtype):
    patient_id, x, ytime, yevent, age, stage_i, stage_ii, race_white = sort_data(path)
    X = torch.from_numpy(x).type(dtype)
    YTIME = torch.from_numpy(ytime).type(dtype)
    YEVENT = torch.from_numpy(yevent).type(dtype)
    AGE = torch.from_numpy(age).type(dtype)
    STAGE_I = torch.from_numpy(stage_i).type(dtype)
    STAGE_II = torch.from_numpy(stage_ii).type(dtype)
    RACE_WHITE = torch.from_numpy(race_white).type(dtype)
    if torch.cuda.is_available():
        X = X.cuda()
        YTIME = YTIME.cuda()
        YEVENT = YEVENT.cuda()
        AGE = AGE.cuda()
        STAGE_I = STAGE_I.cuda()
        STAGE_II = STAGE_II.cuda()
        RACE_WHITE = RACE_WHITE.cuda()
    return(patient_id, X, YTIME, YEVENT, AGE, STAGE_I, STAGE_II, RACE_WHITE)

def load_pathway(path, dtype):
    '''Load a bi-adjacency matrix of pathways and genes, and then covert it to a Pytorch tensor.
    Input:
        path: path to input dataset (which is expected to be a csv file).
        dtype: define the data type of tensor (i.e. dtype=torch.FloatTensor)
    Output:
        PATHWAY_MASK: a Pytorch tensor of the bi-adjacency matrix of pathways & genes.
    '''
    pathway_mask = pd.read_csv(path, index_col = 0).to_numpy()

    PATHWAY_MASK = torch.from_numpy(pathway_mask).type(dtype)
    PATHWAY_MASK = torch.transpose(PATHWAY_MASK, 0, 1)
    ###if gpu is being used
    if torch.cuda.is_available():
        PATHWAY_MASK = PATHWAY_MASK.cuda()
    ###
    return(PATHWAY_MASK)


# In[ ]:


"""
Loss function for KL-PMVAE
"""

def bce_recon_loss(recon_x, x):
    batch_size = x.size(0)
    assert batch_size != 0
    bce_loss = F.binary_cross_entropy(recon_x, x, reduction='sum').div(batch_size)
    return bce_loss

def kl_divergence(mu, logvar):
    batch_size = mu.size(0)
    assert batch_size != 0
    
    klds = -0.5*(1 + logvar - mu.pow(2) - logvar.exp())
    total_kld = klds.sum(1).mean(0, True)
    dimension_wise_kld = klds.mean(0)
    mean_kld = klds.mean(1).mean(0, True)

    return total_kld, dimension_wise_kld, mean_kld


# In[ ]:


"""
Loss function for LFSurv & C-index calculation
"""

def R_set(x):
    n_sample = x.size(0)
    matrix_ones = torch.ones(n_sample, n_sample)
    indicator_matrix = torch.tril(matrix_ones)
    return(indicator_matrix)

def neg_par_log_likelihood(pred, ytime, yevent):
    n_observed = yevent.sum(0)
    ytime_indicator = R_set(ytime)
    if torch.cuda.is_available():
        ytime_indicator = ytime_indicator.cuda()
    risk_set_sum = ytime_indicator.mm(torch.exp(pred)) 
    diff = pred - torch.log(risk_set_sum)
    sum_diff_in_observed = torch.transpose(diff, 0, 1).mm(yevent)
    cost = (- (sum_diff_in_observed / n_observed)).reshape((-1,))
    return(cost)

def c_index(true_T, true_E, pred_risk, include_ties=True):
    """
    Calculate c-index for survival prediction downstream task
    """
    # Ordering true_T, true_E and pred_score in descending order according to true_T
    order = np.argsort(-true_T[:,0].detach().cpu().numpy())

    true_T = true_T.detach().cpu().numpy()[order]
    true_E = true_E.detach().cpu().numpy()[order]
    pred_risk = pred_risk.detach().cpu().numpy()[order]

    # Calculating the c-index
    result = concordance_index(true_T, -pred_risk, true_E)

    return result


# In[ ]:


"""
Early stopping scheme when training LFSurv
"""

class EarlyStopping:
    def __init__(self, patience, verbose=False, delta=0, path="saved_model/sup_checkpoint.pt"):
        
        self.patience = patience
        self.verbose = verbose
        self.counter = 0
        self.epoch_count = 0
        self.best_epoch_num = 1
        self.early_stop = False
        self.max_acc = None
        self.delta = delta
        self.path = path

    def __call__(self, acc, model):
        if self.max_acc is None:
            self.epoch_count += 1
            self.best_epoch_num = self.epoch_count
            self.max_acc = acc
            self.save_checkpoint(model)
        elif acc < self.max_acc + self.delta:
            self.epoch_count += 1
            self.counter += 1
            if self.counter % 20 == 0:
                print(f'EarlyStopping counter: {self.counter} out of {self.patience}')
            if self.counter >= self.patience:
                self.early_stop = True
        else:
            self.epoch_count += 1
            self.best_epoch_num = self.epoch_count
            self.max_acc = acc
            if self.verbose:
                print(f'Validation accuracy increased ({self.max_acc:.6f} --> {acc:.6f}).  Saving model ...')
            self.save_checkpoint(model)
            self.counter = 0

    def save_checkpoint(self, model):
        torch.save(model.state_dict(), self.path)


# In[ ]:


"""
Match patient IDs
"""

def get_match_id(id_1, id_2):
    match_id_cache = []
    for i in id_2["patient_id"]:
        match_id = id_1.index[id_1["patient_id"]==i].tolist()
        match_id_cache += match_id
    return(match_id_cache)


# In[ ]:


"""
Sample from the groups
"""

def splitExprandSample(condition, sampleSize, expr):
    split_expr = expr[condition].T
    split_expr = split_expr.sample(n=sampleSize, axis=1).T
    return split_expr


# In[ ]:




