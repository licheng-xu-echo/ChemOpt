# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 16:41:17 2021

@author: Li-Cheng Xu
"""
import pandas as pd
import numpy as np
from itertools import product

def get_domain(defined_chemical_space):
    '''
    The defined parameters are arranged and combined to generate a chemical space
    '''
    domain_list = [tmp_combine for tmp_combine in product(*[defined_chemical_space[tmp_key] for tmp_key in defined_chemical_space])]
    domain = pd.DataFrame.from_dict({tmp_category:[domain_list[i][idx] for i in range(len(domain_list))] \
                                     for idx,tmp_category in enumerate(defined_chemical_space)})
    return domain
def process_desc(array):
    '''
    Process descriptors, removing descriptors that contain incorrect values (such as NaN), 
    and removing descriptors that have the same value on all compounds
    '''
    array = np.array(array,dtype=np.float32)
    desc_len = array.shape[1]
    rig_idx = []
    for i in range(desc_len):
        try:
            desc_range = array[:,i].max() - array[:,i].min()
            if desc_range != 0 and not np.isnan(desc_range):
                rig_idx.append(i)
        except:
            continue
    array = array[:,rig_idx]
    return array
def maxminscale(array):
    return (array - array.min(axis=0))/(array.max(axis=0)-array.min(axis=0))
def zscorescale(array):
    return (array - array.mean(axis=0))/array.std(axis=0)