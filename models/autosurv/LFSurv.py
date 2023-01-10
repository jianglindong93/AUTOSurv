#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import torch
import torch.nn as nn

dtype = torch.FloatTensor


# In[ ]:


class LFSurv(nn.Module):
    """
    Feed-forward Cox proportional hazard regression network
    """

    def __init__(self, input_n, level_2_dim, Dropout_Rate_1, Dropout_Rate_2):
        super(LFSurv, self).__init__()
        self.input_n = input_n
        self.level_2_dim = level_2_dim
        self.tanh = nn.Tanh()
        
        # Cox fc layers
        self.c_bn_input = nn.BatchNorm1d(self.input_n)
        self.c_fc1 = nn.Linear(self.input_n+4, self.level_2_dim)
        self.c_bn2 = nn.BatchNorm1d(self.level_2_dim)
        self.c_fc2 = nn.Linear(self.level_2_dim, 1, bias = False)
        self.c_fc2.weight.data.uniform_(-0.001, 0.001)
        
        # Dropout
        self.dropout_1 = nn.Dropout(Dropout_Rate_1)
        self.dropout_2 = nn.Dropout(Dropout_Rate_2)
    
    def coxlayer(self, latent_features, c1, c2, c3, c4, s_dropout):
        if s_dropout:
            latent_features = self.dropout_1(latent_features)
        latent_features = self.c_bn_input(latent_features)
        clinical_layer = torch.cat((latent_features, c1, c2, c3, c4), 1)
        hidden_layer = self.tanh(self.c_fc1(clinical_layer))
        if s_dropout:
            hidden_layer = self.dropout_2(hidden_layer)
        hidden_layer = self.c_bn2(hidden_layer)
        y_pred = self.c_fc2(hidden_layer)
        
        return y_pred
    
    def forward(self, latent_features, c1, c2, c3, c4, s_dropout = False):
        y_pred = self.coxlayer(latent_features, c1, c2, c3, c4, s_dropout)
        return y_pred

