#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import math
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.nn.init as init
import torch.optim as optim
from torch.autograd import Variable
import time
from tqdm import tqdm
import lifelines
from lifelines.utils import concordance_index
from Surv_OmiVAE import Surv_OmiVAE_2omics, Surv_OmiVAE_genes, Surv_OmiVAE_mirnas
from utils import sort_data, load_data, load_pathway, bce_recon_loss, kl_divergence, R_set, neg_par_log_likelihood, c_index, EarlyStopping
from train_Surv_OmiVAE import trainSurv_OmiVAE

dtype = torch.FloatTensor


# In[ ]:


"""
This is an example using the breast cancer dataset

Please note that z_dim, Dropout_Rate_1, Dropout_Rate_2, Cutting_Ratio, p1_epoch_num, num_cycles are also hyperparameters 
that can be tuned, in this example we only show the tuning process for some of the hyperparameters
i.e., cox_level2_dim, L2_Lambda, Initial_Learning_Rate
"""

input_n1 = 5147
input_n2 = 865
z_dim = 16
cox_level2_dim = [8, 16, 32]
Dropout_Rate_1 = 0.3 
Dropout_Rate_2 = 0.1
Cutting_Ratio = 0.9
p1_epoch_num = 2400
p2_epoch_num = 500
patience2 = 200
num_cycles = 5
L2_Lambda = [0.1, 0.05, 0.005, 0.0005]
Initial_Learning_Rate = [0.05, 0.01, 0.0075, 0.005, 0.0025]

best_epoch_num = 0

patient_id_train, x_train_gene, ytime_train, yevent_train, age_train, stage_i_train, stage_ii_train, race_white_train = load_data("tune/minmax_normalized/data_train_gene_minmax_tune.csv", dtype)
patient_id_valid, x_valid_gene, ytime_valid, yevent_valid, age_valid, stage_i_valid, stage_ii_valid, race_white_valid = load_data("tune/minmax_normalized/data_valid_gene_minmax_tune.csv", dtype)
pathway_mask_tune = load_pathway("tune/minmax_normalized/pathway_mask.csv", dtype)

_, x_train_mirna, _, _, _, _, _, _ = load_data("tune/minmax_normalized/data_train_mirna_minmax_tune.csv", dtype)
_, x_valid_mirna, _, _, _, _, _, _ = load_data("tune/minmax_normalized/data_valid_mirna_minmax_tune.csv", dtype)

patient_id_train_overall, x_train_gene_overall, ytime_train_overall, yevent_train_overall, age_train_overall, stage_i_train_overall, stage_ii_train_overall, race_white_train_overall = load_data("train_test_split/minmax_normalized/data_train_gene_minmax_overall.csv", dtype)
patient_id_test_overall, x_test_gene_overall, ytime_test_overall, yevent_test_overall, age_test_overall, stage_i_test_overall, stage_ii_test_overall, race_white_test_overall = load_data("train_test_split/minmax_normalized/data_test_gene_minmax_overall.csv", dtype)
pathway_mask_test = load_pathway("train_test_split/minmax_normalized/pathway_mask.csv", dtype)

_, x_train_mirna_overall, _, _, _, _, _, _ = load_data("train_test_split/minmax_normalized/data_train_mirna_minmax_overall.csv", dtype)
_, x_test_mirna_overall, _, _, _, _, _, _ = load_data("train_test_split/minmax_normalized/data_test_mirna_minmax_overall.csv", dtype)

opt_l2 = 0.
opt_lr = 0.
opt_cox_level2_dim = 0
opt_cindex_va = np.float(0)
opt_cindex_tr = np.float(0)
for l2 in L2_Lambda:
    for lr in Initial_Learning_Rate:
        for Coxlevel2_dim in cox_level2_dim:
            _, _, cindex_train, cindex_valid, best_epoch_num_tune = trainSurv_OmiVAE(x_train_gene, x_train_mirna, age_train, stage_i_train, stage_ii_train, race_white_train, ytime_train, yevent_train,
                                                                                     x_valid_gene, x_valid_mirna, age_valid, stage_i_valid, stage_ii_valid, race_white_valid, ytime_valid, yevent_valid,
                                                                                     z_dim, input_n1, input_n2, pathway_mask_tune, Coxlevel2_dim, Dropout_Rate_1, Dropout_Rate_2,
                                                                                     lr, l2, Cutting_Ratio, p1_epoch_num, p2_epoch_num, patience2, num_cycles, dtype,
                                                                                     path = "saved_models/omivae_sup_checkpoint_tune.pt")
            
            if cindex_valid > opt_cindex_va:
                opt_l2 = l2
                opt_lr = lr
                opt_cox_level2_dim = Coxlevel2_dim
                opt_cindex_tr = cindex_train
                opt_cindex_va = cindex_valid
                best_epoch_num = best_epoch_num_tune
            print("L2: %s," %l2, "LR: %s." %lr, "cox_level2_dim: %s." %Coxlevel2_dim)
            print("Training C-index: %s," %cindex_train.round(4), "validation C-index: %s." %cindex_valid.round(4))

train_y_pred, test_y_pred, cindex_train, cindex_test, best_epoch_num_overall = trainSurv_OmiVAE(x_train_gene_overall, x_train_mirna_overall, age_train_overall, stage_i_train_overall, stage_ii_train_overall, race_white_train_overall, ytime_train_overall, yevent_train_overall,
                                                                                                x_test_gene_overall, x_test_mirna_overall, age_test_overall, stage_i_test_overall, stage_ii_test_overall, race_white_test_overall, ytime_test_overall, yevent_test_overall,
                                                                                                z_dim, input_n1, input_n2, pathway_mask_test, opt_cox_level2_dim, Dropout_Rate_1, Dropout_Rate_2,
                                                                                                opt_lr, opt_l2, Cutting_Ratio, p1_epoch_num, p2_epoch_num, patience2, num_cycles, dtype,
                                                                                                path = "saved_models/omivae_sup_checkpoint_overall.pt")
print("Optimal L2: %s," %opt_l2, "optimal LR: %s," %opt_lr, "optimal cox_level2_dim: %s," %opt_cox_level2_dim, "best epoch number in tuning: %s." %best_epoch_num)
print("Optimal training C-index: %s," %opt_cindex_tr.round(4), "optimal validation C-index: %s." %opt_cindex_va.round(4))
print("Training C-index: %s," %cindex_train.round(4), "testing C-index: %s." %cindex_test.round(4))
print("Best epoch num testing: %s." %best_epoch_num_overall)


# In[ ]:




