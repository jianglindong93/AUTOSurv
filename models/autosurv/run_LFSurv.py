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
import time
import lifelines
from lifelines.utils import concordance_index
from LFSurv import LFSurv
from utils import sort_data, load_data, R_set, neg_par_log_likelihood, c_index, EarlyStopping
from train_LFSurv import train_LFSurv

dtype = torch.FloatTensor


# In[ ]:


input_n = 16
level_2_dim = [8, 16, 32]
epoch_num = 500
patience = 200
Initial_Learning_Rate = [0.05, 0.01, 0.0075, 0.005, 0.0025]
L2_Lambda = [0.001, 0.00075, 0.0005, 0.00025, 0.0001]
Dropout_rate_1 = [0.1, 0.3, 0.5]
Dropout_rate_2 = [0.1, 0.3, 0.5]

best_epoch_num = 0

patient_id_train, x_train, ytime_train, yevent_train, age_train, stage_i_train, stage_ii_train, race_white_train = load_data("tune_tr_z_2omics.csv", dtype)
patient_id_valid, x_valid, ytime_valid, yevent_valid, age_valid, stage_i_valid, stage_ii_valid, race_white_valid = load_data("tune_val_z_2omics.csv", dtype)

patient_id_train_overall, x_train_overall, ytime_train_overall, yevent_train_overall, age_train_overall, stage_i_train_overall, stage_ii_train_overall, race_white_train_overall = load_data("tr_z_2omics.csv", dtype)
patient_id_test_overall, x_test_overall, ytime_test_overall, yevent_test_overall, age_test_overall, stage_i_test_overall, stage_ii_test_overall, race_white_test_overall = load_data("tes_z_2omics.csv", dtype)
opt_l2 = 0
opt_lr = 0
opt_dim = 0
opt_dr1 = 0
opt_dr2 = 0

opt_cindex_va = np.float(0)
opt_cindex_tr = np.float(0)
for l2 in L2_Lambda:
    for lr in Initial_Learning_Rate:
        for dim in level_2_dim:
            for dr1 in Dropout_rate_1:
                for dr2 in Dropout_rate_2:
                    _, _, cindex_train, cindex_valid, best_epoch_num_tune = train_LFSurv(x_train, age_train, stage_i_train, stage_ii_train, race_white_train, ytime_train, yevent_train,
                                                                                         x_valid, age_valid, stage_i_valid, stage_ii_valid, race_white_valid, ytime_valid, yevent_valid,
                                                                                         input_n, dim, dr1, dr2, lr, l2, epoch_num, patience, dtype,
                                                                                         path = "saved_models/sup_checkpoint_tune.pt")
                    
                    if cindex_valid > opt_cindex_va:
                        opt_l2 = l2
                        opt_lr = lr
                        opt_dim = dim
                        opt_dr1 = dr1
                        opt_dr2 = dr2
                        opt_cindex_tr = cindex_train
                        opt_cindex_va = cindex_valid
                        best_epoch_num = best_epoch_num_tune
                    print("L2: %s," %l2, "LR: %s." %lr, "dim: %s," %dim, "dr1: %s," %dr1, "dr2: %s." %dr2)
                    print("Training C-index: %s," %cindex_train.round(4), "validation C-index: %s." %cindex_valid.round(4))
train_y_pred, test_y_pred, cindex_train, cindex_test, best_epoch_num_overall = train_LFSurv(x_train_overall, age_train_overall, stage_i_train_overall, stage_ii_train_overall, race_white_train_overall, ytime_train_overall, yevent_train_overall,
                                                                                            x_test_overall, age_test_overall, stage_i_test_overall, stage_ii_test_overall, race_white_test_overall, ytime_test_overall, yevent_test_overall,
                                                                                            input_n, opt_dim, opt_dr1, opt_dr2, opt_lr, opt_l2, epoch_num, patience, dtype,
                                                                                            path = "saved_models/sup_checkpoint_overall.pt")
print("Optimal L2: %s," %opt_l2, "optimal LR: %s," %opt_lr, "optimal dim: %s," %opt_dim, "optimal dr1: %s," %opt_dr1, "optimal dr2: %s," %opt_dr2, "best epoch number in tuning: %s." %best_epoch_num)
print("Optimal training C-index: %s," %opt_cindex_tr.round(4), "optimal validation C-index: %s." %opt_cindex_va.round(4))
print("Testing phase: training C-index: %s," %cindex_train.round(4), "testing C-index: %s." %cindex_test.round(4))

