#!/usr/bin/env python3

"""

Authors: Sean Maden, Abhinav Nellore

Produce tables of hashed features from DNAm beta-value matrices, for samples 
and CpG probes.

Output tables include feature labels (first column) and hashed data (second to 
last column). 


"""

import mmh3; import numpy as np; import pandas as pd

def feature_hash(arr, target_dim=10000):
    low_d_rep = [0 for _ in range(target_dim)]
    for i, el in enumerate(arr):
        hashed = mmh3.hash(str(i))
        if hashed > 0:
            low_d_rep[hashed % target_dim] += arr[i]
        else:
            low_d_rep[hashed % target_dim] -= arr[i]
    return low_d_rep

def make_fhmatrix_autolabel(wf_name, of_name, 
    lnotfloat = ['','NA','NaN'], ndim = 10000):
    """
    Automatically assigns labels as first column of data in wf_name
    """
    with open(of_name, "r") as fr:
        with open(wf_name, "w") as fw:
            for li, line in enumerate(fr):
                line_format = line.replace('\n', '').split(',')
                newrow = line_format[0]; lli = line_format[1::]
                if li > 0:
                    # replace NAs with median values
                    lli_median = np.median([float(ii) for ii in lli 
                        if not ii in lnotfloat])
                    lli_format = [float(ii) if not ii in lnotfloat
                                    else lli_median for ii in lli]
                    lli_fh = feature_hash(lli_format, target_dim = ndim)
                    newrow = newrow + ',' + ','.join([str(ii) for ii in lli_fh])
                    newrow = newrow + '\n'
                    print('new data sample: '+newrow[0:100])
                    fw.write(newrow)    
                print("Finished with line number "+str(li))
    return None

def make_fhmatrix_specifylabels(labels_list, wf_name, of_name, 
    lnotfloat = ['','NA','NaN'], ndim = 10000):
    """ Generated matrix of hashed features from a table

    Generates a matrix of hashed features, of dim ndim

    Arguments
    * labels_list : Column names of the output table
    * wf_name : Name of the output table
    * of_name : Name of table to be hashed. Rows will be hashed. First line 
                is colnames and thus skipped.
    Output
    * None, produces a new hashed features table of dim nrow_of_name by 
    ndim+1 (columns), where first column has feature labels


    """
    with open(of_name, "r") as fr:
        with open(wf_name, "w") as fw:
            for li, line in enumerate(fr):
                lli = line.replace('\n', '').split(',')[1::]
                newrow = labels_list[li] # append label to new row
                # replace NAs with median values
                lli_median = np.median([float(ii) for ii in lli 
                    if not ii in lnotfloat])
                lli_format = [float(ii) if not ii in lnotfloat
                                else lli_median for ii in lli]
                lli_fh = feature_hash(lli_format, target_dim = ndim)
                newrow = newrow + ',' + ','.join([str(ii) for ii in lli_fh])
                newrow = newrow + '\n'
                print('Writing new row: '+newrow[0:100])
                fw.write(newrow)
                print("Finished with line number "+str(li))
    return None










def __main__:
    """
    """

# HM450k
# load the labels
samples_labels_fname = "samples_hm450k.csv"
probes_labels_fname = "probes_hm450k.csv"

samples_labels = []
with open(samples_labels_fname, "r") as of:
    for li, line in enumerate(of):
        if li > 0:
            line_format = line.split(',')[1]
            line_format = line_format.replace('\n', '').replace('"', '')
            samples_labels.append(line_format)

probes_labels = []
with open(probes_labels_fname, "r") as of:
    for li, line in enumerate(of):
        if li > 0:
            line_format = line.split(',')[1]
            line_format = line_format.replace('\n', '').replace('"', '')
            probes_labels.append(line_format)

of_name = "bval_hm450k_cols-cgids_rows-samples.csv"
wf_name = "bval-fh10k-samples_hm450k.csv"

make_fhmatrix_specifylabels(labels_list = samples_labels, 
    wf_name = wf_name, of_name = of_name, ndim = 10000)

of_name = "bval_hm450k_cols-samples_rows-cgids.csv"
wf_name = "bval-fh1k-cgprobes_hm450k.csv"

make_fhmatrix_specifylabels(labels_list = probes_labels, 
    wf_name = wf_name, of_name = of_name, ndim = 1000)

# EPIC
samples_labels_fname = "samples_epic.csv"
probes_labels_fname = "probes_epic.csv"

samples_labels = []
with open(samples_labels_fname, "r") as of:
    for li, line in enumerate(of):
        if li > 0:
            line_format = line.split(',')[1]
            line_format = line_format.replace('\n', '').replace('"', '')
            samples_labels.append(line_format)

probes_labels = []
with open(probes_labels_fname, "r") as of:
    for li, line in enumerate(of):
        if li > 0:
            line_format = line.split(',')[1]
            line_format = line_format.replace('\n', '').replace('"', '')
            probes_labels.append(line_format)

of_name = "bval_epic_cols-cgids_rows-samples.csv"
wf_name = "bval-fh10k-samples_epic.csv"

make_fhmatrix_specifylabels(labels_list = samples_labels, 
    wf_name = wf_name, of_name = of_name, ndim = 10000)

of_name = "bval_epic_cols-samples_rows-cgids.csv"
wf_name = "bval-fh1k-cgprobes_epic.csv"

make_fhmatrix_specifylabels(labels_list = probes_labels, 
    wf_name = wf_name, of_name = of_name, ndim = 1000)
