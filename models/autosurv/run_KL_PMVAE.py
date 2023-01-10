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
from KL_PMVAE import KL_PMVAE_2omics, KL_PMVAE_genes, KL_PMVAE_mirnas
from utils import sort_data, load_data, load_pathway, bce_recon_loss, kl_divergence, get_match_id
from train_KL_PMVAE import train_KL_PMVAE

dtype = torch.FloatTensor


# In[ ]:


"""
Example listed here is for TCGA-BRCA data. If using TCGA-OV data, please adjust the parameters/covariates accordingly
"""

input_n1 = 2699
input_n2 = 516
z_dim = [8, 16, 32, 64]
EPOCH_NUM = [800, 1200, 1600, 2000, 2400]
NUM_CYCLES = [2, 4, 5, 10]
CUTTING_RATIO = [0.3, 0.5, 0.7, 0.9]
Initial_Learning_Rate = [0.1, 0.05, 0.005, 0.001]
L2_Lambda = [0.1, 0.05, 0.005, 0.0005]

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
opt_dim = 0
opt_epoch_num = 0.
opt_num_cycle = 0.
opt_cr = 0.
opt_loss = torch.Tensor([float("Inf")])
if torch.cuda.is_available():
    opt_loss = opt_loss.cuda()
for l2 in L2_Lambda:
    for lr in Initial_Learning_Rate:
        for Z in z_dim:
            for Epoch_num in EPOCH_NUM:
                for Num_cycles in NUM_CYCLES:
                    for cutting_ratio in CUTTING_RATIO:
                        _, _, _, _, train_loss_unsup, eval_loss_unsup = train_KL_PMVAE(x_train_gene, x_train_mirna, x_valid_gene, x_valid_mirna,
                                                                                       Z, input_n1, input_n2, pathway_mask_tune,
                                                                                       lr, l2, cutting_ratio, Epoch_num, Num_cycles, dtype,
                                                                                       path = "saved_models/unsup_checkpoint_tune.pt")
            
                        if eval_loss_unsup < opt_loss:
                            opt_l2 = l2
                            opt_lr = lr
                            opt_dim = Z
                            opt_epoch_num = Epoch_num
                            opt_num_cycle = Num_cycles
                            opt_cr = cutting_ratio
                            opt_loss = eval_loss_unsup
                        print("num_epoch: %s," %Epoch_num, "num_cycles: %s," %Num_cycles, "cutting_ratio: %s." %cutting_ratio)
                        print("L2: %s," %l2, "LR: %s," %lr, "z_dim: %s," %Z, "loss in validation: %s," %np.array(eval_loss_unsup.detach().cpu().numpy()).round(4), "loss in training: %s." %np.array(train_loss_unsup.detach().cpu().numpy()).round(4))

print("Optimal num epoch: %s," %opt_epoch_num, "optimal num cycles: %s," %opt_num_cycle, "optimal cutting ratio: %s." %opt_cr)
print("Optimal L2: %s," %opt_l2, "optimal LR: %s," %opt_lr, "optimal z_dim: %s." %opt_dim)

train_mean, train_logvar, test_mean, test_logvar, train_loss_unsup, test_loss_unsup = train_KL_PMVAE(x_train_gene_overall, x_train_mirna_overall, x_test_gene_overall, x_test_mirna_overall,
                                                                                                     opt_dim, input_n1, input_n2, pathway_mask_test,
                                                                                                     opt_lr, opt_l2, opt_cr, opt_epoch_num, opt_num_cycle, dtype, save_model = False,
                                                                                                     path = "saved_models/unsup_checkpoint_overall.pt")
print("Loss in testing: %s," %np.array(test_loss_unsup.detach().cpu().numpy()).round(4), "loss in training: %s." %np.array(train_loss_unsup.detach().cpu().numpy()).round(4))


# In[ ]:


tr_z = train_mean
tes_z = test_mean

print("Training sample size: %s," %tr_z.size()[0], "testing sample size: %s." %tes_z.size()[0])

processed_tr_pre = torch.cat((ytime_train_overall, yevent_train_overall, age_train_overall, stage_i_train_overall, stage_ii_train_overall, race_white_train_overall, tr_z), 1)
processed_tes_pre = torch.cat((ytime_test_overall, yevent_test_overall, age_test_overall, stage_i_test_overall, stage_ii_test_overall, race_white_test_overall, tes_z), 1)

z_count = np.array(list(range(1, tr_z.size()[1]+1, 1))).astype('str')
z_names = np.char.add('Z_', z_count).tolist()

processed_tr = pd.DataFrame(processed_tr_pre, columns = ['OS.time', 'OS', 'age', 'stage_i', 'stage_ii', 'race_white'] + z_names)
processed_tr = processed_tr.astype(float)
processed_tr = pd.concat([patient_id_train_overall, processed_tr], axis=1)

processed_tes = pd.DataFrame(processed_tes_pre, columns = ['OS.time', 'OS', 'age', 'stage_i', 'stage_ii', 'race_white'] + z_names)
processed_tes = processed_tes.astype(float)
processed_tes = pd.concat([patient_id_test_overall, processed_tes], axis=1)

tune_id_train = get_match_id(patient_id_train_overall, patient_id_train)
tune_id_valid = get_match_id(patient_id_train_overall, patient_id_valid)
processed_tr_tune = processed_tr.iloc[tune_id_train]
processed_val_tune = processed_tr.iloc[tune_id_valid]

processed_tr.to_csv("tr_z_2omics.csv", index=False)
processed_tr_tune.to_csv("tune_tr_z_2omics.csv", index=False)
processed_val_tune.to_csv("tune_val_z_2omics.csv", index=False)
processed_tes.to_csv("tes_z_2omics.csv", index=False)


# In[ ]:




