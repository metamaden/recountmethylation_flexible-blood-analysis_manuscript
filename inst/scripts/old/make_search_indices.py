#!/usr/bin/env python3

"""
    Authors: Sean Maden & Abhinav Nellore
    
    Make the search indices for samples and probes from h5se objects. This 
    script shows how to read in data from the Beta-value .csv tables produced
    as describedin the `makesi_bvalcsv.R` script. 

    The Beta-value tables are feature hashed into either 10,000 features 
    (sample tables) or 1,000 features (probe tables) before being indexed using 
    the Hierarchical Navigable Small Worlds method implemented in the `hnswlib` 
    Python library.

"""

import mmh3; import numpy as np; import pandas as pd
import numpy as np; import pandas as pd
import hnswlib, sys, os, h5py, re, time, hnswlib, pickle, time, random
import faulthandler; faulthandler.enable()
random.seed(0)

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
    low_d_rep = [0 for _ in range(target_dim)]
    for i, el in enumerate(arr):
        hashed = mmh3.hash(str(i))
        if hashed > 0:
            low_d_rep[hashed % target_dim] += arr[i]
        else:
            low_d_rep[hashed % target_dim] -= arr[i]
    return low_d_rep

def make_fhmatrix_specifylabels(labels_list, wf_name, of_name, 
    lnotfloat = ['','NA','NaN'], ndim = 10000):
    """ Make the hashed feaures matrix a Beta-values table

    Saves a .csv table of ndim total hashed features (columns) by elements, 
    where the number of elements is equal to the row count in the of_name table.
    Missing/NA values are replaced by the row median. New rows are processively
    written to the file wf_name.

    Arguments:
        * labels_list : Column names of the output table
        * wf_name : Name of the output table
        * of_name : Name of table to be hashed. Rows will be hashed. First line 
                    is colnames and thus skipped.
    Returns:
        * None, produces a new hashed features .csv table of dim nrow_of_name by 
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

def make_si(fname, index_name, dindex_name, space_val = 'l2', threads = 10, 
    efc_val = 2000, m_val = 1000, ef_val = 2000):
    """ Makes a new search index from a csv table

    Arguments:
        * fname: Name of CSV file containing index data.
        * index_name: Name of new search index file to save.
        * dindex_name: Name of new index dict, with str labels, to save.
        * space_val: Space value for the index.
        * threads: Number of threads for processing the index.
        * efc_val: EFC value for the index.
        * m_val: M value for the index.
        * ef_val: EF value for the index.

    Returns
        * None, saves the new index and index dict
    """
    print("Loading the dataset...")
    df=pd.read_csv(fname, sep=',',header=None)
    df_str_labels = df.iloc[0::,0]
    df = df.iloc[0::,1::] # subset to exclude string labels
    num_elements = df.shape[0]
    dim = df.shape[1]
    print("Making the new index...")
    p_iter = hnswlib.Index(space = space_val, dim = dim)
    p_iter.set_num_threads(threads)
    df_labels = np.arange(num_elements)
    p_iter.init_index(max_elements = num_elements, 
        ef_construction = efc_val, M = m_val)
    p_iter.set_ef(ef_val)
    print("Adding new data to the search index...")
    p_iter.add_items(df, df_labels)
    print("Saving the search index...")
    p_iter.save_index(index_name)
    print("Saving the index dict with str labels, as pickle object...")
    dindex = {'index' : p_iter, 'strlabels' : df_str_labels}
    file_open = open(dindex_name, "wb")
    pickle.dump(obj = dindex, file = file_open)
    return None

def siquery(data, siobject, si_labels = False, kval = 100):
    """ Perform a lookup on an hnswlib-saved search index

    Arguments:
        * data : Vector of feature hashed data, of dims R x C, where C is the 
                same as in siobject/search index object.
        * siobject : Search index object previously saved using hnswlib.
        * si_labels : Element labels corresponding to the index positions in 
                siobject.
        * kval : K number of nearest neighbors to return from siobject lookup.

    Returns:
        * Vector of elements, where elements are indices if si_labels = False, 
                or else element labels if si_labels is provided.

    """
    # rowdat = pd.DataFrame([random.uniform(-10, 10) for ii in range(10000)]).T
    # rowdat.shape
    kval_labels, kval_distances = siobject.knn_query(data, k = kval)
    element_list  = kval_labels
    return element_list

def main():
    """ Make the search indices from hash tables of Beta-values

    Makes 2 search indices per platform, one each for probes and samples. Takes
    as input the Beta-values CSV tables, makes the hash tables of 1k columns for 
    probes and 10k columns for samples, and finally makes the index objects 
    using the Hierarchical Navigable Small Worlds algorithm implemented in the
    `hnswlib` package.

    """
    # Make fh tables from beta-values
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
    # Write new fh tables
    of_name = "bval_hm450k_cols-cgids_rows-samples.csv"
    wf_name = "bval-fh10k-samples_hm450k.csv"
    make_fhmatrix_specifylabels(labels_list = samples_labels, 
        wf_name = wf_name, of_name = of_name, ndim = 10000)
    of_name = "bval_hm450k_cols-samples_rows-cgids.csv"
    wf_name = "bval-fh1k-cgprobes_hm450k.csv"
    make_fhmatrix_specifylabels(labels_list = probes_labels, 
        wf_name = wf_name, of_name = of_name, ndim = 1000)
    # EPIC
    # load labels
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
    # write new fh tables
    of_name = "bval_epic_cols-cgids_rows-samples.csv"
    wf_name = "bval-fh10k-samples_epic.csv"
    make_fhmatrix_specifylabels(labels_list = samples_labels, 
        wf_name = wf_name, of_name = of_name, ndim = 10000)
    of_name = "bval_epic_cols-samples_rows-cgids.csv"
    wf_name = "bval-fh1k-cgprobes_epic.csv"
    make_fhmatrix_specifylabels(labels_list = probes_labels, 
        wf_name = wf_name, of_name = of_name, ndim = 1000)
    
    # make the search index objects
    # HM450K samples si
    make_si(fname = "bval-fh10k-samples_hm450k.csv", 
        index_name = "hnsw-index_samples_hm450k.pickle", 
        dindex_name = "dindex_samples_hm450k.pickle")
    # HM450K probes si
    make_si(fname = "bval-fh1k-cgprobes_hm450k.csv", 
        index_name = "hnsw-index_probes_hm450k.pickle", 
        dindex_name = "dindex_probes_hm450k.pickle")
    # EPIC samples si
    make_si(fname = "bval-fh10k-samples_epic.csv", 
        index_name = "hnsw-index_samples_epic.pickle", 
        dindex_name = "dindex_samples_epic.pickle")
        # EPIC probes si
    make_si(fname = "bval-fh1k-cgprobes_epic.csv", 
        index_name = "hnsw-index_probes_epic.pickle", 
        dindex_name = "dindex_probes_epic.pickle")
        return None

if __name__ == "__main__":
    """ Make the search index objects
    
    Make search index objects from h5se database files.

    This script makes the search index objects from beta-value tables previously
    saved from h5se database files. Refer to the script "make_bval_h5.R" for
    details.

    """
    main()