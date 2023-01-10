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
from utils import bce_recon_loss, kl_divergence

dtype = torch.FloatTensor


# In[ ]:


def train_KL_PMVAE(train_x1, train_x2, eval_x1, eval_x2,
                   z_dim, input_n1, input_n2, Pathway_Mask,
                   Learning_Rate, L2, Cutting_Ratio, p1_epoch_num, num_cycles, dtype, save_model = False,
                   path = "saved_model/unsup_checkpoint.pt"):
    net = KL_PMVAE_2omics(z_dim, input_n1, input_n2, Pathway_Mask)
    train_real = torch.cat((train_x1, train_x2), 1)
    eval_real = torch.cat((eval_x1, eval_x2), 1)
    
    if torch.cuda.is_available():
        net.cuda()
    opt = optim.Adam(net.parameters(), lr = Learning_Rate, weight_decay = L2)
    cycle_iter = p1_epoch_num // num_cycles
    start_time = time.time()
    for epoch in range(p1_epoch_num):
        tmp = float(epoch%cycle_iter)/cycle_iter
        
        if tmp == 0:
            beta = 0.1
        elif tmp <= Cutting_Ratio:
            beta = tmp/Cutting_Ratio
        else:
            beta = 1
        
        net.train()
        opt.zero_grad()
        
        mean, logvar,  _, recon_x1, recon_x2 = net(train_x1, train_x2, s_dropout = True)
        recon_x = torch.cat((recon_x1, recon_x2), 1)
        recon_loss = bce_recon_loss(recon_x, train_real)
        total_kld, _, _ = kl_divergence(mean, logvar)
        loss_unsup = recon_loss + beta*total_kld
        
        loss_unsup.backward()
        opt.step()
        
        if (epoch+1) % 100 == 0:
            net.eval()
            train_mean, train_logvar, _, train_recon1, train_recon2 = net(train_x1, train_x2, s_dropout = False)
            train_recon = torch.cat((train_recon1, train_recon2), 1)
            train_recon_loss = bce_recon_loss(train_recon, train_real)
            train_total_kld, _, _ = kl_divergence(train_mean, train_logvar)
            train_loss_unsup = train_recon_loss + beta*train_total_kld
            
            net.eval()
            eval_mean, eval_logvar, _, eval_recon1, eval_recon2 = net(eval_x1, eval_x2, s_dropout = False)
            eval_recon = torch.cat((eval_recon1, eval_recon2), 1)
            eval_recon_loss = bce_recon_loss(eval_recon, eval_real)
            eval_total_kld, _, _ = kl_divergence(eval_mean, eval_logvar)
            eval_loss_unsup = eval_recon_loss + beta*eval_total_kld
            
            temp_epoch = epoch +1
            print("Epoch: %s," %temp_epoch, "Loss in training: %s," %np.array(train_loss_unsup.detach().cpu().numpy()).round(4), "loss in validation: %s." %np.array(eval_loss_unsup.detach().cpu().numpy()).round(4))
    
    if save_model:
        print("Saving model...")
        torch.save(net.state_dict(), path)
        print("Model saved.")
    
    print(np.array(time.time() - start_time).round(2))
    return (train_mean, train_logvar, eval_mean, eval_logvar, train_loss_unsup, eval_loss_unsup)


# In[ ]:




