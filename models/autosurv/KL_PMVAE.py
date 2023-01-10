#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import math
import torch
import torch.nn as nn
from torch.autograd import Variable

dtype = torch.FloatTensor


# In[ ]:


class KL_PMVAE_2omics(nn.Module):

    def __init__(self, z_dim, input_n1, input_n2, Pathway_Mask):
        super(KL_PMVAE_2omics, self).__init__()
        self.pathway_mask = Pathway_Mask
        level_2_dim = Pathway_Mask.size(0)
        self.level_2_dim = level_2_dim
        self.input_n2 = input_n2
        self.sigm = nn.Sigmoid()
        self.relu = nn.ReLU()
        
        self.z_dim = z_dim
        # ENCODER fc layers
        # level 1
        self.e_fc1 = nn.Linear(input_n1, level_2_dim)
        self.e_bn1 = nn.BatchNorm1d(level_2_dim)
        
        # level 2
        self.e_fc2_mean = nn.Linear(level_2_dim + input_n2, z_dim)
        self.e_fc2_logvar = nn.Linear(level_2_dim + input_n2, z_dim)
        self.e_bn2_mean = nn.BatchNorm1d(z_dim)
        self.e_bn2_logvar = nn.BatchNorm1d(z_dim)
        
        # DECODER fc layers
        # level 2
        self.d_fc2 = nn.Linear(z_dim, level_2_dim + input_n2)
        self.d_bn2 = nn.BatchNorm1d(level_2_dim + input_n2)
        
        # level 1
        self.d_fc1 = nn.Linear(level_2_dim, input_n1)
        self.d_bn1 = nn.BatchNorm1d(input_n1)
        
        # Dropout
        self.dropout_1 = nn.Dropout(0.)
        self.dropout_2 = nn.Dropout(0.)
    
    def _reparameterize(self, mean, logvar, z = None):
        std = logvar.mul(0.5).exp()    
        if z is None:
            if torch.cuda.is_available():
                z = Variable(torch.cuda.FloatTensor(std.size()).normal_(0, 1))
            else:
                z = Variable(torch.FloatTensor(std.size()).normal_(0, 1))
        return z.mul(std) + mean
    
    def encode(self, x1, x2, s_dropout):
        if s_dropout:
            x1 = self.dropout_1(x1)
        self.e_fc1.weight.data = self.e_fc1.weight.data.mul(self.pathway_mask)
        level2_layer = self.relu(self.e_bn1(self.e_fc1(x1)))
        
        conc_layer = torch.cat((level2_layer, x2), 1)
        if s_dropout:
            conc_layer = self.dropout_2(conc_layer)
        latent_mean = self.e_bn2_mean(self.e_fc2_mean(conc_layer))
        latent_logvar = self.e_bn2_logvar(self.e_fc2_logvar(conc_layer))
        
        return latent_mean, latent_logvar
    
    def decode(self, z):
        level2_layer = self.d_bn2(self.d_fc2(z))
        level2_x1 = level2_layer.narrow(1, 0, self.level_2_dim)
        recon_x2 = self.sigm(level2_layer.narrow(1, self.level_2_dim, self.input_n2))
        
        self.d_fc1.weight.data = self.d_fc1.weight.data.mul(torch.transpose(self.pathway_mask, 0, 1))
        recon_x1 = self.sigm(self.d_bn1(self.d_fc1(level2_x1)))
        
        return recon_x1, recon_x2
        
    def forward(self, x1, x2, s_dropout = False):
        mean, logvar = self.encode(x1, x2, s_dropout)
        z = self._reparameterize(mean, logvar)
        recon_x1, recon_x2 = self.decode(z)
                    
        return mean, logvar, z, recon_x1, recon_x2


# In[ ]:


class KL_PMVAE_genes(nn.Module):

    def __init__(self, z_dim, input_n, Pathway_Mask):
        super(KL_PMVAE_genes, self).__init__()
        self.pathway_mask = Pathway_Mask
        level_2_dim = Pathway_Mask.size(0)
        self.sigm = nn.Sigmoid()
        self.relu = nn.ReLU()
        
        self.z_dim = z_dim
        # ENCODER fc layers
        # level 1
        self.e_fc1 = nn.Linear(input_n, level_2_dim)
        self.e_bn1 = nn.BatchNorm1d(level_2_dim)
        
        # level 2
        self.e_fc2_mean = nn.Linear(level_2_dim, z_dim)
        self.e_fc2_logvar = nn.Linear(level_2_dim, z_dim)
        self.e_bn2_mean = nn.BatchNorm1d(z_dim)
        self.e_bn2_logvar = nn.BatchNorm1d(z_dim)
        
        # DECODER fc layers
        # level 2
        self.d_fc2 = nn.Linear(z_dim, level_2_dim)
        self.d_bn2 = nn.BatchNorm1d(level_2_dim)
        
        # level 1
        self.d_fc1 = nn.Linear(level_2_dim, input_n)
        self.d_bn1 = nn.BatchNorm1d(input_n)
        
        # Dropout
        self.dropout_1 = nn.Dropout(0.)
        self.dropout_2 = nn.Dropout(0.)
    
    def _reparameterize(self, mean, logvar, z = None):
        std = logvar.mul(0.5).exp()    
        if z is None:
            if torch.cuda.is_available():
                z = Variable(torch.cuda.FloatTensor(std.size()).normal_(0, 1))
            else:
                z = Variable(torch.FloatTensor(std.size()).normal_(0, 1))
        return z.mul(std) + mean
    
    def encode(self, x1, s_dropout):
        if s_dropout:
            x1 = self.dropout_1(x1)
        self.e_fc1.weight.data = self.e_fc1.weight.data.mul(self.pathway_mask)
        level2_layer = self.relu(self.e_bn1(self.e_fc1(x1)))
        
        if s_dropout:
            level2_layer = self.dropout_2(level2_layer)
        latent_mean = self.e_bn2_mean(self.e_fc2_mean(level2_layer))
        latent_logvar = self.e_bn2_logvar(self.e_fc2_logvar(level2_layer))
        
        return latent_mean, latent_logvar
    
    def decode(self, z):
        level2_layer = self.relu(self.d_bn2(self.d_fc2(z)))
        
        self.d_fc1.weight.data = self.d_fc1.weight.data.mul(torch.transpose(self.pathway_mask, 0, 1))
        
        recon_x = self.sigm(self.d_bn1(self.d_fc1(level2_layer)))
        
        return recon_x
        
    def forward(self, x1, s_dropout = False):
        mean, logvar = self.encode(x1, s_dropout)
        z = self._reparameterize(mean, logvar)
        recon_x = self.decode(z)
                    
        return mean, logvar, z, recon_x


# In[ ]:


class KL_PMVAE_mirnas(nn.Module):

    def __init__(self, z_dim, input_n):
        super(KL_PMVAE_mirnas, self).__init__()
        self.sigm = nn.Sigmoid()
        self.relu = nn.ReLU()
        self.z_dim = z_dim
        
        # ENCODER fc layers
        # level 1
        self.e_fc1_mean = nn.Linear(input_n, z_dim)
        self.e_fc1_logvar = nn.Linear(input_n, z_dim)
        self.e_bn1_mean = nn.BatchNorm1d(z_dim)
        self.e_bn1_logvar = nn.BatchNorm1d(z_dim)
        
        # DECODER fc layers
        # level 1
        self.d_fc1 = nn.Linear(z_dim, input_n)
        self.d_bn1 = nn.BatchNorm1d(input_n)
        
        # Dropout
        self.dropout_1 = nn.Dropout(0.)
    
    def _reparameterize(self, mean, logvar, z = None):
        std = logvar.mul(0.5).exp()    
        if z is None:
            if torch.cuda.is_available():
                z = Variable(torch.cuda.FloatTensor(std.size()).normal_(0, 1))
            else:
                z = Variable(torch.FloatTensor(std.size()).normal_(0, 1))
        return z.mul(std) + mean
    
    def encode(self, x1, s_dropout):
        if s_dropout:
            x1 = self.dropout_1(x1)
        latent_mean = self.e_bn1_mean(self.e_fc1_mean(x1))
        latent_logvar = self.e_bn1_logvar(self.e_fc1_logvar(x1))
        
        return latent_mean, latent_logvar
    
    def decode(self, z):
        recon_x = self.sigm(self.d_bn1(self.d_fc1(z)))
        
        return recon_x
        
    def forward(self, x1, s_dropout = False):
        mean, logvar = self.encode(x1, s_dropout)
        z = self._reparameterize(mean, logvar)
        recon_x = self.decode(z)
                    
        return mean, logvar, z, recon_x

