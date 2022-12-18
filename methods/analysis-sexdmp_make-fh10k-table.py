#!/usr/bin/env python3

"""
Authors: Sean Maden, Abhinav Nellore

Get hashed features from some tabular data.

"""

import mmh3
import numpy as np
import pandas as pd
import hnswlib, sys, os, re, pickle, random
from time import time
import faulthandler
faulthandler.enable()

def feature_hash(arr, target_dim=10000):
    """ Perform feature hashing on an array of data
    
    Perform feature hashing on the data in arr, into a vector of target_dim 
    total hashed features.

    Arguments:
        * arr: An array of values to be hashed.
        * target_dim: The target number of hashed values.
    Returns:
        * low_d_rep, or an array of hashed values of len == target_dim

    """ 
    low_d_rep = [0.0 for _ in range(target_dim)]
    for i, el in enumerate(arr):
        hashed = mmh3.hash(str(i))
        if hashed > 0.0:
            low_d_rep[hashed % target_dim] += arr[i]
        else:
            low_d_rep[hashed % target_dim] -= arr[i]
    return low_d_rep

def make_fhmatrix_autolabel(wf_name, of_name, lnotfloat = ['','NA','NaN'], 
    ndim = 10000, lstart = 0):
    """
    
    Get the hashed features table from an input data table. This function 
    automatically sets row labels as the first column of data in the input data
    table at `of_name.`

    Arguments:
        * wf_name: Name/path of hashed features table to write (required, 
            string, rows = samples, cols = hashed features).
        * of_name: Name/path of table to hash (required, tring, rows =  samples, 
            cols = probes).
        * lnotfloat: List of expected missing value symbols (required, ['','NA',
            'NaN']). These are replaced by the row-wise meidans of non-missing 
            values.
        * ndim: Number of hashed features (integer, 1000).
        * lstart: Line to start reading (required, int, 0).
    Returns:
        * None, saves new hashed features table to `wf_name`.

    """
    with open(of_name, "r") as fr:
        with open(wf_name, "w") as fw:
            for li, line in enumerate(fr):
                line_format = line.replace('\n', '').split(',')
                newrow = line_format[0]; lli = line_format[1::]
                if li >= lstart:
                    # replace NAs with median values
                    lli_median = np.median([float(ii) for ii in lli 
                        if not ii in lnotfloat])
                    lli_format = [float(ii) if not ii in lnotfloat
                                    else lli_median for ii in lli]
                    lli_fh = feature_hash(lli_format, target_dim = ndim)
                    newrow = newrow + ',' + ','.join([str(ii) for ii in lli_fh])
                    newrow = newrow + '\n'
                    print('Found new sample: '+newrow[0:100])
                    fw.write(newrow)    
                print("Finished with line number "+str(li))
    return None

import mmh3
import numpy as np

ndim = 10000
lstart = 1
lnotfloat = ['','NA','NaN']
of_name = "mval-sexdmp_pbmc_2platforms.csv"
wf_name = "mval-fh10k_pbmc.csv"
with open(of_name, "r") as fr:
    with open(wf_name, "w") as fw:
        for li, line in enumerate(fr):
            line_format = line.replace('\n', '').split(',')
            lli = line_format[0::]
            if li >= lstart:
                # replace NAs with median values
                lli_median = np.median([float(ii) for ii in lli 
                    if not ii in lnotfloat])
                lli_format = [float(ii) if not ii in lnotfloat
                                else lli_median for ii in lli]
                lli_fh = feature_hash(lli_format, target_dim = ndim)
                newrow = ','.join([str(ii) for ii in lli_fh])
                newrow = newrow + '\n'
                print('Found new sample: '+newrow[0:100])
                fw.write(newrow)    
            print("Finished with line number "+str(li))

with open(of_name, "r") as fr:
    for li, line in enumerate(fr):
        line_format = line.replace('\n', '').split(',')
        lli = line_format[0::]
        line_format2 = line.replace('\n', '').split(',')
        lli2 = line_format2[0::]