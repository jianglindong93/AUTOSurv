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
from torch.autograd import Variable
from KL_PMVAE import KL_PMVAE_2omics
from utils import get_match_id

dtype = torch.FloatTensor


# In[ ]:


def load_data_deepshap(data1, data2, dtype):
    EXPR1 = data1.drop(["patient_id", "OS", "OS.time", "age", "race_white", "stage_i", "stage_ii"], axis = 1).values.astype(np.float)
    EXPR2 = data2.drop(["patient_id", "OS", "OS.time", "age", "race_white", "stage_i", "stage_ii"], axis = 1).values.astype(np.float)
    EXPR1 = torch.from_numpy(EXPR1).type(dtype)
    EXPR2 = torch.from_numpy(EXPR2).type(dtype)
        
    return (EXPR1, EXPR2)


# In[ ]:


def load_pathway(path, dtype):
    '''Load a bi-adjacency matrix of pathways, and then covert it to a Pytorch tensor.
    Input:
        path: path to input dataset (which is expected to be a csv file).
        dtype: define the data type of tensor (i.e. dtype=torch.FloatTensor)
    Output:
        PATHWAY_MASK: a Pytorch tensor of the bi-adjacency matrix of pathways & genes.
    '''
    pathway_mask = pd.read_csv(path, index_col = 0).to_numpy()

    PATHWAY_MASK = torch.from_numpy(pathway_mask).type(dtype)
    PATHWAY_MASK = torch.transpose(PATHWAY_MASK, 0, 1)
    
    return(PATHWAY_MASK)


# In[ ]:


train_data_gene = pd.read_csv("train_test_split/minmax_normalized/data_train_gene_minmax_overall.csv")
train_data_mirna = pd.read_csv("train_test_split/minmax_normalized/data_train_mirna_minmax_overall.csv")

gene_names = list(train_data_gene.columns)[7:]
mirna_names = list(train_data_mirna.columns)[7:]
feature_names = gene_names + mirna_names


# In[ ]:


upper_data = pd.read_csv("saved_models/higher_PI_train.csv")
lower_data = pd.read_csv("saved_models/lower_PI_train.csv")


# In[ ]:


pathway_mask = pd.read_csv("train_test_split/minmax_normalized/pathway_mask.csv", index_col = 0)
pathway_names = list(pathway_mask.columns)
mask_names = pathway_names + mirna_names
pathway_mask_int = load_pathway("train_test_split/minmax_normalized/pathway_mask.csv", dtype)


# In[ ]:


"""
Breast cancer dataset used as example here
"""
# number of genes
input_n1 = 2699
# number of miRNAs
input_n2 = 516
# number of latent features: one of the hyperparameters, determined via grid search during training/validation process
z_dim = 16


# In[ ]:


"""
DeepSHAP

These codes were adapted from the work of Withnell et al. https://academic.oup.com/bib/article/22/6/bbab315/6353242
Please check their original implementation of DeepSHAP for more details at https://github.com/zhangxiaoyu11/XOmiVAE
"""

class Explainer(object):
    """ This is the superclass of all explainers.
    """

    def shap_values(self, X):
        raise Exception("SHAP values not implemented for this explainer!")

    def attributions(self, X):
        return self.shap_values(X)

class PyTorchDeepExplainer(Explainer):
    """
    This class has been adapted to explain OmiVAE. It is important that the correct output of the model
    (0=z dimension, 2=mean) as well as the dimension within this latent space. We allow a dimension to be chosen to explain from, and adjustable outputs to be explained. 
    Lundberg et al., 2017: http://papers.nips.cc/paper/7062-a-unified-approach-to-interpreting-model-predictions.pdf
    """

    def __init__(self, model, data1, data2, outputNumber, dim, explainLatentSpace):
        
        #data = list(load_data_deepshap(data, dtype))
        data = list(load_data_deepshap(data1, data2, dtype))
        
        # check if we have multiple inputs
        self.multi_input = False
        if type(data) == list:
            self.multi_input = True
        else:
            data = [data]
        self.data = data
        self.layer = None
        self.input_handle = None
        self.interim = False
        self.interim_inputs_shape = None
        self.expected_value = None  # to keep the DeepExplainer base happy
        if type(model) == tuple:

            self.interim = True
            model, layer = model
            model = model.eval()
            self.layer = layer
            self.add_target_handle(self.layer)

            # if we are taking an interim layer, the 'data' is going to be the input
            # of the interim layer; we will capture this using a forward hook
            with torch.no_grad():
                _ = model(*data)
                #model.e_fc1.weight.data = model.e_fc1.weight.data.mul(model.pathway_mask)
                #_ = model.relu(model.e_bn1(model.e_fc1(data[0])))
                interim_inputs = self.layer.target_input
                #print("Interim input type: %s." %type(interim_inputs))
                if type(interim_inputs) is tuple:
                    # this should always be true, but just to be safe
                    self.interim_inputs_shape = [i.shape for i in interim_inputs]
                else:
                    self.interim_inputs_shape = [interim_inputs.shape]
            self.target_handle.remove()
            del self.layer.target_input
        self.model = model.eval()
        self.multi_output = False
        self.num_outputs = 1
        with torch.no_grad():
            outputs = model(*data)
            #model.e_fc1.weight.data = model.e_fc1.weight.data.mul(model.pathway_mask)
            #outputs = model.relu(model.e_bn1(model.e_fc1(data[0])))
            #This is where specifies whether we want to explain the mean or z output
            if type(outputs) != tuple:
                output = outputs
            else:
                output = outputs[outputNumber]
            self.outputNum=outputNumber
            # Chosen dimension
            self.dim=None
            self.explainLatent = False
            if explainLatentSpace:
                self.explainLatent=True
                self.dimension=dim
                output = output[:, dim]
                output = output.reshape(output.shape[0], 1)
            # also get the device everything is running on
            self.device = output.device
            if output.shape[1] > 1:
                self.multi_output = True
                self.num_outputs = output.shape[1]
            self.expected_value = output.mean(0).cpu().numpy()

    def add_target_handle(self, layer):

        input_handle = layer.register_forward_hook(get_target_input)
        self.target_handle = input_handle

    def add_handles(self, model, forward_handle, backward_handle):
        """
        Add handles to all non-container layers in the model.
        Recursively for non-container layers
        """
        handles_list = []
        model_children = list(model.children())
        if model_children:
            for child in model_children:
                handles_list.extend(self.add_handles(child, forward_handle, backward_handle))
        else:  # leaves
            handles_list.append(model.register_forward_hook(forward_handle))
            handles_list.append(model.register_backward_hook(backward_handle))

        return handles_list

    def remove_attributes(self, model):
        """
        Removes the x and y attributes which were added by the forward handles
        Recursively searches for non-container layers
        """
        for child in model.children():
            if 'nn.modules.container' in str(type(child)):
                self.remove_attributes(child)
            else:
                try:
                    del child.x
                except AttributeError:
                    pass
                try:
                    del child.y
                except AttributeError:
                    pass

    def gradient(self, idx, inputs):

        self.model.zero_grad()
        X = [x.requires_grad_() for x in inputs]
        
        output = self.model(*X)
        #self.model.e_fc1.weight.data = self.model.e_fc1.weight.data.mul(self.model.pathway_mask)
        #output = self.model.relu(self.model.e_bn1(self.model.e_fc1(X[0])))

        #Specify the output to change
        if type(output) != tuple:
            outputs = output
        else:
            outputs = output[self.outputNum]

        #Specify the dimension to explain
        if self.explainLatent==True:

            outputs = outputs[:, self.dimension]
            outputs = outputs.reshape(outputs.shape[0], 1)


        selected = [val for val in outputs[:, idx]]

        grads = []
        if self.interim:
            interim_inputs = self.layer.target_input
            for idx, input in enumerate(interim_inputs):
                grad = torch.autograd.grad(selected, input,
                                           retain_graph=True if idx + 1 < len(interim_inputs) else None,
                                           allow_unused=True)[0]
                if grad is not None:
                    grad = grad.cpu().numpy()
                else:
                    grad = torch.zeros_like(interim_inputs[idx]).cpu().numpy()
                grads.append(grad)
            del self.layer.target_input
            return grads, [i.detach().cpu().numpy() for i in interim_inputs]
        else:
            for idx, x in enumerate(X):
                grad = torch.autograd.grad(selected, x,
                                           retain_graph=True if idx + 1 < len(X) else None,
                                           allow_unused=True)[0]
                if grad is not None:
                    grad = grad.cpu().numpy()
                else:
                    grad = torch.zeros_like(X[idx]).cpu().numpy()
                grads.append(grad)
            return grads

    def shap_values(self, X1, X2, ranked_outputs=None, output_rank_order="max", check_additivity=False):

        # X ~ self.model_input
        # X_data ~ self.data
        
        #X = list(load_data_deepshap(X, dtype))
        X = list(load_data_deepshap(X1, X2, dtype))

        # check if we have multiple inputs
        if not self.multi_input:
            assert type(X) != list, "Expected a single tensor model input!"
            X = [X]
        else:
            assert type(X) == list, "Expected a list of model inputs!"


        X = [x.detach().to(self.device) for x in X]

        # if ranked output is given then this code is run and only the 'max' value given is explained
        if ranked_outputs is not None and self.multi_output:
            with torch.no_grad():
                model_output_values = self.model(*X)
                #self.model.e_fc1.weight.data = self.model.e_fc1.weight.data.mul(self.model.pathway_mask)
                #model_output_values = self.model.relu(self.model.e_bn1(self.model.e_fc1(X[0])))
                
                #Whithnell's change to adjust for the additional outputs in VAE model
                model_output_values = model_output_values[self.outputNum]

            # rank and determine the model outputs that we will explain

            if output_rank_order == "max":
                _, model_output_ranks = torch.sort(model_output_values, descending=True)
            elif output_rank_order == "min":
                _, model_output_ranks = torch.sort(model_output_values, descending=False)
            elif output_rank_order == "max_abs":
                _, model_output_ranks = torch.sort(torch.abs(model_output_values), descending=True)
            else:
                assert False, "output_rank_order must be max, min, or max_abs!"

        else:
            # outputs an array of 0s so we know we are explaining the first value
            model_output_ranks = (torch.ones((X[0].shape[0], self.num_outputs)).int() *
                                  torch.arange(0, self.num_outputs).int())

        # add the gradient handles

        handles = self.add_handles(self.model, add_interim_values, deeplift_grad)
        if self.interim:
            self.add_target_handle(self.layer)

        # compute the attributions
        output_phis = []

        for i in range(model_output_ranks.shape[1]):

            phis = []
            #phis are shapLundberg values

            if self.interim:
                for k in range(len(self.interim_inputs_shape)):
                    phis.append(np.zeros((X[0].shape[0], ) + self.interim_inputs_shape[k][1: ]))
            else:
                for k in range(len(X)):
                    phis.append(np.zeros(X[k].shape))
            #shape is 5 as testing 5 samples
            for j in range(X[0].shape[0]):

                # tile the inputs to line up with the background data samples
                tiled_X = [X[l][j:j + 1].repeat(
                                   (self.data[l].shape[0],) + tuple([1 for k in range(len(X[l].shape) - 1)])) for l
                           in range(len(X))]
                joint_x = [torch.cat((tiled_X[l], self.data[l]), dim=0) for l in range(len(X))]
                # run attribution computation graph
                feature_ind = model_output_ranks[j, i]
                sample_phis = self.gradient(feature_ind, joint_x)
                # assign the attributions to the right part of the output arrays
                if self.interim:
                    sample_phis, output = sample_phis
                    x, data = [], []
                    for i in range(len(output)):
                        x_temp, data_temp = np.split(output[i], 2)
                        x.append(x_temp)
                        data.append(data_temp)
                    for l in range(len(self.interim_inputs_shape)):
                        phis[l][j] = (sample_phis[l][self.data[l].shape[0]:] * (x[l] - data[l])).mean(0)
                else:
                    for l in range(len(X)):
                        phis[l][j] = (torch.from_numpy(sample_phis[l][self.data[l].shape[0]:]).to(self.device) * (X[l][j: j + 1] - self.data[l])).cpu().numpy().mean(0)
            output_phis.append(phis[0] if not self.multi_input else phis)


        # cleanup; remove all gradient handles
        for handle in handles:
            handle.remove()
        self.remove_attributes(self.model)
        if self.interim:
            self.target_handle.remove()

        if not self.multi_output:
            return output_phis[0]
        elif ranked_outputs is not None:
            # EW: returns a list... only want first value
            return output_phis, model_output_ranks
        else:
            return output_phis

# Module hooks


def deeplift_grad(module, grad_input, grad_output):
    """The backward hook which computes the deeplift
    gradient for an nn.Module
    """
    # first, get the module type
    module_type = module.__class__.__name__

    # first, check the module is supported
    if module_type in op_handler:

        if op_handler[module_type].__name__ not in ['passthrough', 'linear_1d']:
            return op_handler[module_type](module, grad_input, grad_output)
    else:
        print('Warning: unrecognized nn.Module: {}'.format(module_type))
        return grad_input


def add_interim_values(module, input, output):
    """The forward hook used to save interim tensors, detached
    from the graph. Used to calculate the multipliers
    """
    try:
        del module.x
    except AttributeError:
        pass
    try:
        del module.y
    except AttributeError:
        pass
    module_type = module.__class__.__name__

    if module_type in op_handler:
        func_name = op_handler[module_type].__name__

        # First, check for cases where we don't need to save the x and y tensors
        if func_name == 'passthrough':
            pass
        else:
            # check only the 0th input varies
            for i in range(len(input)):
                if i != 0 and type(output) is tuple:
                    assert input[i] == output[i], "Only the 0th input may vary!"
            # if a new method is added, it must be added here too. This ensures tensors
            # are only saved if necessary
            if func_name in ['maxpool', 'nonlinear_1d']:
                # only save tensors if necessary
                if type(input) is tuple:
                    setattr(module, 'x', torch.nn.Parameter(input[0].detach()))
                else:
                    setattr(module, 'x', torch.nn.Parameter(input.detach()))
                if type(output) is tuple:
                    setattr(module, 'y', torch.nn.Parameter(output[0].detach()))
                else:
                    setattr(module, 'y', torch.nn.Parameter(output.detach()))
            if module_type in failure_case_modules:
                input[0].register_hook(deeplift_tensor_grad)


def get_target_input(module, input, output):
    """A forward hook which saves the tensor - attached to its graph.
    Used if we want to explain the interim outputs of a model
    """
    try:
        del module.target_input
    except AttributeError:
        pass
    setattr(module, 'target_input', input)

# Whithnell:
# From the documentation: "The current implementation will not have the presented behavior for
# complex Module that perform many operations. In some failure cases, grad_input and grad_output
# will only contain the gradients for a subset of the inputs and outputs.
# The tensor hook below handles such failure cases (currently, MaxPool1d). In such cases, the deeplift
# grad should still be computed, and then appended to the complex_model_gradients list. The tensor hook
# will then retrieve the proper gradient from this list.


failure_case_modules = ['MaxPool1d']


def deeplift_tensor_grad(grad):
    return_grad = complex_module_gradients[-1]
    del complex_module_gradients[-1]
    return return_grad


complex_module_gradients = []


def passthrough(module, grad_input, grad_output):
    """No change made to gradients"""
    return None


def maxpool(module, grad_input, grad_output):
    pool_to_unpool = {
        'MaxPool1d': torch.nn.functional.max_unpool1d,
        'MaxPool2d': torch.nn.functional.max_unpool2d,
        'MaxPool3d': torch.nn.functional.max_unpool3d
    }
    pool_to_function = {
        'MaxPool1d': torch.nn.functional.max_pool1d,
        'MaxPool2d': torch.nn.functional.max_pool2d,
        'MaxPool3d': torch.nn.functional.max_pool3d
    }
    delta_in = module.x[: int(module.x.shape[0] / 2)] - module.x[int(module.x.shape[0] / 2):]
    dup0 = [2] + [1 for i in delta_in.shape[1:]]
    # we also need to check if the output is a tuple
    y, ref_output = torch.chunk(module.y, 2)
    cross_max = torch.max(y, ref_output)
    diffs = torch.cat([cross_max - ref_output, y - cross_max], 0)

    # all of this just to unpool the outputs
    with torch.no_grad():
        _, indices = pool_to_function[module.__class__.__name__](
            module.x, module.kernel_size, module.stride, module.padding,
            module.dilation, module.ceil_mode, True)
        xmax_pos, rmax_pos = torch.chunk(pool_to_unpool[module.__class__.__name__](
            grad_output[0] * diffs, indices, module.kernel_size, module.stride,
            module.padding, list(module.x.shape)), 2)
    org_input_shape = grad_input[0].shape  # for the maxpool 1d
    grad_input = [None for _ in grad_input]
    grad_input[0] = torch.where(torch.abs(delta_in) < 1e-7, torch.zeros_like(delta_in),
                           (xmax_pos + rmax_pos) / delta_in).repeat(dup0)
    if module.__class__.__name__ == 'MaxPool1d':
        complex_module_gradients.append(grad_input[0])
        # the grad input that is returned doesn't matter, since it will immediately be
        # be overridden by the grad in the complex_module_gradient
        grad_input[0] = torch.ones(org_input_shape)
    return tuple(grad_input)


def linear_1d(module, grad_input, grad_output):
    """No change made to gradients."""
    return None


def nonlinear_1d(module, grad_input, grad_output):
    delta_out = module.y[: int(module.y.shape[0] / 2)] - module.y[int(module.y.shape[0] / 2):]

    delta_in = module.x[: int(module.x.shape[0] / 2)] - module.x[int(module.x.shape[0] / 2):]
    dup0 = [2] + [1 for i in delta_in.shape[1:]]
    # handles numerical instabilities where delta_in is very small by
    # just taking the gradient in those cases
    grads = [None for _ in grad_input]
    grads[0] = torch.where(torch.abs(delta_in.repeat(dup0)) < 1e-6, grad_input[0],
                           grad_output[0] * (delta_out / delta_in).repeat(dup0))
    return tuple(grads)


op_handler = {}

# passthrough ops, where we make no change to the gradient
op_handler['Dropout3d'] = passthrough
op_handler['Dropout2d'] = passthrough
op_handler['Dropout'] = passthrough
op_handler['AlphaDropout'] = passthrough

op_handler['Conv1d'] = linear_1d
op_handler['Conv2d'] = linear_1d
op_handler['Conv3d'] = linear_1d
op_handler['ConvTranspose1d'] = linear_1d
op_handler['ConvTranspose2d'] = linear_1d
op_handler['ConvTranspose3d'] = linear_1d
op_handler['Linear'] = linear_1d
op_handler['AvgPool1d'] = linear_1d
op_handler['AvgPool2d'] = linear_1d
op_handler['AvgPool3d'] = linear_1d
op_handler['AdaptiveAvgPool1d'] = linear_1d
op_handler['AdaptiveAvgPool2d'] = linear_1d
op_handler['AdaptiveAvgPool3d'] = linear_1d
op_handler['BatchNorm1d'] = linear_1d
op_handler['BatchNorm2d'] = linear_1d
op_handler['BatchNorm3d'] = linear_1d

op_handler['LeakyReLU'] = nonlinear_1d
op_handler['ReLU'] = nonlinear_1d
op_handler['ELU'] = nonlinear_1d
op_handler['Sigmoid'] = nonlinear_1d
op_handler["Tanh"] = nonlinear_1d
op_handler["Softplus"] = nonlinear_1d
op_handler['Softmax'] = nonlinear_1d

op_handler['MaxPool1d'] = maxpool
op_handler['MaxPool2d'] = maxpool
op_handler['MaxPool3d'] = maxpool


# In[ ]:


"""
The sampling procedure is modified for interpreting KL_PMVAE
"""

def splitExprandSample(condition_data, sampleSize, expr):
    match_id = get_match_id(expr, condition_data)
    split_expr = expr.iloc[match_id]
    split_expr.index = range(0, split_expr.shape[0], 1)
    split_expr = split_expr.T
    split_expr = split_expr.sample(n=sampleSize, axis=1).T
    return split_expr


# In[ ]:


"""
Specify if you want to obtain importance scores for genes/miRNA (i.e., input factors) (imp_feature = True, imp_mask = False) 
or pathways (imp_feature = False, imp_mask = True)
"""

def getTopShapValues(shap_vals, numberOfTopFeatures, path, file_count, featureNames, maskNames = None, imp_feature = True, imp_mask = False, absolute=True):
    multiple_input = False
    if type(shap_vals) == list:
        multiple_input = True
        shap_values = None
        for l in range(len(shap_vals)):
            if shap_values is not None:
                shap_values = np.concatenate((shap_values, shap_vals[l]), axis=1)
            else:
                shap_values = shap_vals[l]
        shap_vals = shap_values
    
    if absolute:
        vals = np.abs(shap_vals).mean(0)
    else:
        vals = shap_vals.mean(0)
    
    if imp_feature:
        feature_names = featureNames
    if imp_mask:
        feature_names = maskNames

    feature_importance = pd.DataFrame(list(zip(feature_names, vals)),
                                      columns=['features', 'importance_vals'])
    feature_importance.sort_values(by=['importance_vals'], ascending=False, inplace=True)

    mostImp_shap_values = feature_importance.head(numberOfTopFeatures)
    print(mostImp_shap_values)

    feature_importance.to_csv(path + "/unsup_feature_imp_" + file_count + "_input_all_top6.csv")
    #feature_importance.to_csv(path + "/unsup_feature_imp_" + file_count + "_1st_hidden_all_top6.csv")
    #feature_importance.to_csv(path + "/unsup_feature_imp_specific_pathway_" + file_count + "_top6.csv")
    """
    print(mostImp_shap_values)
    print("least importance absolute values")
    feature_importance.sort_values(by=['feature_importance_vals'], ascending=True, inplace=True)
    leastImp_shap_values = feature_importance.head(numberOfTopFeatures)
    print(leastImp_shap_values)
    """
    return mostImp_shap_values


# In[ ]:


def UnSupShapExplainer(train_overall_df_gene, train_overall_df_mirna, condition_data1, condition_data2, path, file_count, dimension, featureNames, maskNames, imp_feature, imp_mask):
    KL_PMVAE_model = KL_PMVAE_2omics(z_dim, input_n1, input_n2, pathway_mask_int)
    
    #load trained model
    KL_PMVAE_model.load_state_dict(torch.load('saved_models/unsup_checkpoint_overall.pt', map_location=torch.device('cpu')))
    
    upper_group_sample_gene = splitExprandSample(condition_data=condition_data1, sampleSize=100, expr=train_overall_df_gene)
    upper_mirna_match_id = get_match_id(train_overall_df_mirna, upper_group_sample_gene)
    upper_group_sample_mirna = train_overall_df_mirna.iloc[upper_mirna_match_id]
    
    lower_group_sample_gene = splitExprandSample(condition_data=condition_data2, sampleSize=100, expr=train_overall_df_gene)
    lower_mirna_match_id = get_match_id(train_overall_df_mirna, lower_group_sample_gene)
    lower_group_sample_mirna = train_overall_df_mirna.iloc[lower_mirna_match_id]
    
    if imp_feature:
        e = PyTorchDeepExplainer(KL_PMVAE_model, lower_group_sample_gene, lower_group_sample_mirna, outputNumber=0, dim=dimension, explainLatentSpace=True)
    if imp_mask:
        e = PyTorchDeepExplainer((KL_PMVAE_model, KL_PMVAE_model.e_fc2_mean), lower_group_sample_gene, lower_group_sample_mirna, outputNumber=0, dim=dimension, explainLatentSpace=True)
    print("calculating shap values")
    shap_values_obtained = e.shap_values(upper_group_sample_gene, upper_group_sample_mirna)
    
    print("calculated shap values")
    most_imp  = getTopShapValues(shap_vals=shap_values_obtained, numberOfTopFeatures=10, path=path, file_count=file_count, featureNames=featureNames, maskNames=maskNames, imp_feature=imp_feature, imp_mask=imp_mask, absolute=True)
    
    return most_imp


# In[ ]:


"""
identify key input factors (KIFs) (i.e., genes/miRNAs) for the important latent features;
important latent features were found by applying DeepSHAP to trained LFSurv
"""


important_latent = ["Z_5", "Z_2", "Z_10", "Z_8", "Z_1", "Z_9"]
#z_count = np.array(list(range(1, 17, 1))).astype('str')
#important_latent = np.char.add('Z_', z_count).tolist()
imp_features_cache = []
for i in important_latent:
    print("Latent dimension: " + i + ".")
    temp_dim = int(i.split("_")[1])-1
    temp_most_imp_fea = UnSupShapExplainer(train_overall_df_gene=train_data_gene, train_overall_df_mirna=train_data_mirna, condition_data1=upper_data, condition_data2=lower_data, 
                                           path='saved_models', file_count = i,
                                           dimension=temp_dim, featureNames=feature_names, maskNames=mask_names, imp_feature=True, imp_mask=False)
    imp_features_cache = imp_features_cache + list(temp_most_imp_fea.features)


# In[ ]:


freq = {} 
for item in imp_features_cache: 
    if (item in freq): 
        freq[item] += 1
    else: 
        freq[item] = 1


# In[ ]:


freq_summary = pd.DataFrame.from_dict(freq, orient ='index')
freq_summary.columns = ['Frequency']
freq_summary.sort_values(by=['Frequency'], ascending=False, inplace=True)
big_freq_summary = freq_summary[freq_summary['Frequency'] > 1]

big_freq_summary.to_csv("saved_models/top6Z_over1_impinput.csv")
#big_freq_summary.to_csv("saved_models/top6Z_over1_imp1sthidden.csv")


# In[ ]:


"""
Identify important genes for important pathways
"""

most_important_features = UnSupShapExplainer(train_overall_df_gene=train_data_gene, train_overall_df_mirna=train_data_mirna, condition_data1=upper_data, condition_data2=lower_data, 
                                             path='saved_models', file_count = "R.HSA.163560",
                                             dimension=276, featureNames=feature_names, maskNames=mask_names, imp_feature=True, imp_mask=False)


# In[ ]:




