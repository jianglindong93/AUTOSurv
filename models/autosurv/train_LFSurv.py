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
from utils import R_set, neg_par_log_likelihood, c_index, EarlyStopping

dtype = torch.FloatTensor


# In[ ]:


"""
Please note the difference in covariate names between the two cancer types when inputing the variables
"""

def train_LFSurv(train_x1, train_age, train_stage_i, train_stage_ii, train_race_white, train_ytime, train_yevent,
                 eval_x1, eval_age, eval_stage_i, eval_stage_ii, eval_race_white, eval_ytime, eval_yevent,
                 input_n, level_2_dim, Dropout_Rate_1, Dropout_Rate_2, Learning_Rate, L2, epoch_num, patience, dtype,
                 path = "saved_model/sup_checkpoint.pt"):
    net = LFSurv(input_n, level_2_dim, Dropout_Rate_1, Dropout_Rate_2)
    
    early_stopping_sup = EarlyStopping(patience = patience, verbose = False, path = path)
    
    if torch.cuda.is_available():
        net.cuda()
    opt = optim.Adam(net.parameters(), lr = Learning_Rate, weight_decay = L2)
    
    start_time = time.time()
    for epoch in range(epoch_num):
        net.train()
        opt.zero_grad()
        
        y_pred = net(train_x1, train_age, train_stage_i, train_stage_ii, train_race_white, s_dropout = True)
        loss_sup = neg_par_log_likelihood(y_pred, train_ytime, train_yevent)
        
        loss_sup.backward()
        opt.step()
        
        net.eval()
        eval_y_pred = net(eval_x1, eval_age, eval_stage_i, eval_stage_ii, eval_race_white, s_dropout = False)
        eval_cindex = c_index(eval_ytime, eval_yevent, eval_y_pred)
        
        early_stopping_sup(eval_cindex, net)
        if early_stopping_sup.early_stop:
            print("Early stopping, number of epochs: ", epoch)
            print('Save model of Epoch {:d}'.format(early_stopping_sup.best_epoch_num))
            break
        if (epoch+1) % 100 == 0:
            net.eval()
            train_y_pred = net(train_x1, train_age, train_stage_i, train_stage_ii, train_race_white, s_dropout = False)
            train_cindex = c_index(train_ytime, train_yevent, train_y_pred)
            print ("Training C-index: %s," % train_cindex.round(4), "validation C-index: %s." % eval_cindex.round(4))
    
    print("Loading model, best epoch: %s." %early_stopping_sup.best_epoch_num)
    net.load_state_dict(torch.load(path, map_location=torch.device('cpu')))
    
    net.eval()
    train_y_pred = net(train_x1, train_age, train_stage_i, train_stage_ii, train_race_white, s_dropout = False)
    train_cindex = c_index(train_ytime, train_yevent, train_y_pred)
    
    net.eval()
    eval_y_pred = net(eval_x1, eval_age, eval_stage_i, eval_stage_ii, eval_race_white, s_dropout = False)
    eval_cindex = c_index(eval_ytime, eval_yevent, eval_y_pred)
    
    print ("Final training C-index: %s," % train_cindex.round(4), "final validation C-index: %s." % eval_cindex.round(4))
    time_elapse = np.array(time.time() - start_time).round(2)
    print("Total time elapse: %s." %time_elapse)
    
    return (train_y_pred, eval_y_pred, train_cindex, eval_cindex, early_stopping_sup.best_epoch_num)

