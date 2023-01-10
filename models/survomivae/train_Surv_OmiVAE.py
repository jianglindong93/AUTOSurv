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
from utils import bce_recon_loss, kl_divergence, R_set, neg_par_log_likelihood, c_index, EarlyStopping

dtype = torch.FloatTensor


# In[ ]:


def trainSurv_OmiVAE(train_x1, train_x2, train_age, train_stage_i, train_stage_ii, train_race_white, train_ytime, train_yevent,
                     eval_x1, eval_x2, eval_age, eval_stage_i, eval_stage_ii, eval_race_white, eval_ytime, eval_yevent,
                     z_dim, input_n1, input_n2, Pathway_Mask, cox_level2_dim, Dropout_Rate_1, Dropout_Rate_2,
                     Learning_Rate, L2, Cutting_Ratio, p1_epoch_num, p2_epoch_num, patience2, num_cycles, dtype,
                     path = "saved_model/surv_omivae_checkpoint.pt"):
    net = Surv_OmiVAE_2omics(z_dim, input_n1, input_n2, Pathway_Mask, cox_level2_dim, Dropout_Rate_1, Dropout_Rate_2)
    train_real = torch.cat((train_x1, train_x2), 1)
    eval_real = torch.cat((eval_x1, eval_x2), 1)
    
    early_stopping_sup = EarlyStopping(patience = patience2, verbose = False, path = path)
    
    if torch.cuda.is_available():
        net.cuda()
    opt = optim.Adam(net.parameters(), lr = Learning_Rate, weight_decay = L2)
    cycle_iter = p1_epoch_num // num_cycles
    start_time = time.time()
    print("First stage training in progress...")
    for epoch in tqdm(range(p1_epoch_num)):
        tmp = float(epoch%cycle_iter)/cycle_iter
        
        if tmp == 0:
            beta = 0.1
        elif tmp <= Cutting_Ratio:
            beta = tmp/Cutting_Ratio
        else:
            beta = 1
        
        net.train()
        opt.zero_grad()
        
        _, recon_x1, recon_x2, mean, logvar, _ = net(train_x1, train_x2, train_age, train_stage_i, train_stage_ii, train_race_white, s_dropout = False, phase = "unsup")
        recon_x = torch.cat((recon_x1, recon_x2), 1)
        recon_loss = bce_recon_loss(recon_x, train_real)
        total_kld, _, _ = kl_divergence(mean, logvar)
        loss_unsup = recon_loss + beta*total_kld
        
        loss_unsup.backward()
        opt.step()
        
        """if (epoch+1) % 100 == 0:
            net.eval()
            _, train_recon1, train_recon2, train_mean, train_logvar, _ = net(train_x1, train_x2, train_age, train_stage_i, train_stage_ii, train_race_white, s_dropout = False, phase = "unsup")
            train_recon = torch.cat((train_recon1, train_recon2), 1)
            train_recon_loss = bce_recon_loss(train_recon, train_real)
            train_total_kld, _, _ = kl_divergence(train_mean, train_logvar)
            train_loss_unsup = train_recon_loss + beta*train_total_kld
            
            net.eval()
            _, eval_recon1, eval_recon2, eval_mean, eval_logvar, _ = net(eval_x1, eval_x2, eval_age, eval_stage_i, eval_stage_ii, eval_race_white, s_dropout = False, phase = "unsup")
            eval_recon = torch.cat((eval_recon1, eval_recon2), 1)
            eval_recon_loss = bce_recon_loss(eval_recon, eval_real)
            eval_total_kld, _, _ = kl_divergence(eval_mean, eval_logvar)
            eval_loss_unsup = eval_recon_loss + beta*eval_total_kld
            
            temp_epoch = epoch +1
            print("Stage I Epoch: %s," %temp_epoch, "loss in training: %s," %np.array(train_loss_unsup.detach().cpu().numpy()).round(4), "loss in validation: %s." %np.array(eval_loss_unsup.detach().cpu().numpy()).round(4))"""
    print("First stage training complete.")
    
    for epoch in range(p1_epoch_num, p1_epoch_num+p2_epoch_num+1):
        net.train()
        opt.zero_grad()
        
        _, _, _, _, _, y_pred = net(train_x1, train_x2, train_age, train_stage_i, train_stage_ii, train_race_white, s_dropout = True, phase = 'sup')
        loss_sup = neg_par_log_likelihood(y_pred, train_ytime, train_yevent)
        
        loss_sup.backward()
        opt.step()
        
        net.eval()
        _, _, _, _, _, eval_y_pred = net(eval_x1, eval_x2, eval_age, eval_stage_i, eval_stage_ii, eval_race_white, s_dropout = False, phase = 'sup')
        #eval_loss_sup = neg_par_log_likelihood(eval_y_pred, eval_ytime, eval_yevent)
        eval_cindex = c_index(eval_ytime, eval_yevent, eval_y_pred)
        
        early_stopping_sup(eval_cindex, net)
        if early_stopping_sup.early_stop:
            temp_epoch = epoch - p1_epoch_num + 1
            print("Early stopping, number of epochs: %s." %temp_epoch)
            print('Save model of Epoch {:d}'.format(early_stopping_sup.best_epoch_num))
            break
        if (epoch+1) % 100 == 0:
            net.eval()
            _, _, _, _, _, train_y_pred = net(train_x1, train_x2, train_age, train_stage_i, train_stage_ii, train_race_white, s_dropout = False, phase = 'sup')
            #train_loss_sup = neg_par_log_likelihood(train_y_pred, train_ytime, train_yevent)
            train_cindex = c_index(train_ytime, train_yevent, train_y_pred)
            print ("Training C-index: %s," %train_cindex.round(4), "validation C-index: %s." %eval_cindex.round(4))
    
    print("Loading model, best epoch: %s." %early_stopping_sup.best_epoch_num)
    net.load_state_dict(torch.load(path, map_location=torch.device('cpu')))
    
    net.eval()
    _, _, _, _, _, train_y_pred = net(train_x1, train_x2, train_age, train_stage_i, train_stage_ii, train_race_white, s_dropout = False, phase = "sup")
    train_cindex = c_index(train_ytime, train_yevent, train_y_pred)
    
    net.eval()
    _, _, _, _, _, eval_y_pred = net(eval_x1, eval_x2, eval_age, eval_stage_i, eval_stage_ii, eval_race_white, s_dropout = False, phase = "sup")
    eval_cindex = c_index(eval_ytime, eval_yevent, eval_y_pred)
    
    print ("Final training C-index: %s," %train_cindex.round(4), "final validation C-index: %s." %eval_cindex.round(4))
    time_elapse = np.array(time.time() - start_time).round(2)
    print("Total time elapse: %s." %time_elapse)
    return (train_y_pred, eval_y_pred, train_cindex, eval_cindex, early_stopping_sup.best_epoch_num)


# In[ ]:




