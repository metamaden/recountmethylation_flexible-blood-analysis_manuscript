#!/usr/bin/env python3

"""

Author: Sean Maden

Test the probe and sample search indices for each platform using 
pre-identified blood subgroups.

Start by making the group labels from existing .csv files

QUERY SPEEDS:

    * hm450k, all samples, 8123 elements:
        - k = 10, time = 49.7sec user
        - k = 100, time = 50.4sec user
        - k = 500, time = 50sec user
        - k = 1000, time = 50sec user

"""

import pandas as pd; import numpy as np; from time import time
import os, pickle

def siquery(data, siobject, si_labels = [], kval = 100):
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
    print("Querying "+str(data.shape[0])+" elements in data with k = "+
        str(kval)+" nearest neighbors...")
    time_start = time()
    kval_labels, kval_distances = siobject.knn_query(data, k = kval)
    time_stop = time()
    print("Query completed in, time elapsed: "+str(time_stop-time_start))
    if len(si_labels) > 0:
        print("Applying labels to query results..."); kval_str_labels = []
        for ll in kval_labels:
            kval_str_labels.append([si_labels[ii] for ii in ll])
        return kval_str_labels, kval_labels, kval_distances
    return kval_labels, kval_distances

def make_dfk(groupcsv_fpath = "si2_blood-md-2platforms.csv", lgroup = ["all", 
    "cord_blood", "whole_blood", "peripheral_blood_mononuclear_cells"],
    lplatform = ["hm450k", "epic"], lk = [10, 100, 500, 1000],
    ldsi_fpath = ["dindex_samples_hm450k.pickle", "dindex_samples_epic.pickle"],
    lfhtable_fpath = ["bval-fh10k-samples_hm450k.csv", 
    "bval-fh10k-samples_epic.csv"]):
    """ Make a data frame from queries of ki members in lk

    Get the sample labels for the k nearest neighbors on a series of queries.

    Arguments:
        * groupcsv_fpath: Path to the CSV table defining the sample group
            members.
        * lgroup: The group terms to consider. If "all", select all non-NA 
            samples.
        * lplatform: The platform string names. Members correspond to memebers 
            provided in ldsi_fpath and lfhtable_fpath arguments.
        * lk: List of the k nearest neighbors to query.
        * ldsi_fpath: List of the pickled dictionaries containing the search 
            indices and string labels of index elements. Their ordering 
            matches members provided in lplatform and lfhtable_fpath args.
        * lfhtable_fpath: List of paths to hashed feature tables. These members
            correspond to members of lplatform and ldsi_fpath args.
        
    Returns:
        * dfk_final, a DataFrame object concatenating the query results across 
            groups.

    """
    print("Reading in the metadata from file at path '"+str(bloodsubfname)+"'...")
    dfall = pd.read_csv(bloodsubfname, sep=','); ldfk = []
    for group in lgroup:
        for pi, platform in enumerate(lplatform):
            dsi_fname = ldsi_fname[pi]; fht_fname = lfhtable_fpath[pi]
            print("Beginning group '"+group+"'' on platform '"+platform+"'...")
            print("Subsetting the sample metadata...")
            if group == "all":
                dfsub = dfall[dfall["blood_subgroup"].isnull() == False]
            else:
                dfsub = dfall[dfall["blood_subgroup"] == group]
                dfsub = dfsub[dfsub["platform"] == platform]
            print("After subsetting, retained "+str(dfsub.shape[0])+" samples.")
            lgsm = [gsm for gsm in dfsub["gsm"]]
            print("Getting fh data for samples in group..."); fhdict = {}
            with open(fht_fname, "r") as of:
                for line in of:
                    lline = line.split(","); gsm = lline[0].split(".")[0]
                    if gsm in lgsm:
                        fhdict[gsm] = [float(fhi.replace("\n", "")) for fhi in lline[1::]]
            dfi = pd.DataFrame.from_dict(fhdict).T
            of_sidict = open(dsi_fname, "rb"); sidict = pickle.load(of_sidict)
            dfk = pd.DataFrame(fhdict.keys())
            print("Beginning queries on k values in lk...")
            for ii, ki in enumerate(lk):
                kval_str_labels, kval_labels, kval_distances = siquery(dfi, 
                    siobject = sidict["index"], si_labels = sidict["strlabels"], 
                    kval = ki)
                kli_format = [";".join(ii) for ii in kval_str_labels]
                dfki = pd.DataFrame(kli_format)
                dfk[str(ki)] = [ii for ii in dfki[0]]
            dfk["subgroup"] = group; dfk["platform"] = platform; ldfk.append(dfk)
            print("Finished platform '"+platform+"' for group '"+group+"'.")
        print("Finished all platforms for group '"+group+"'.")
    dfk_final = pd.concat(ldfk)
    print("Finished all queries. Returning results DataFrame of "+
        str(dfk_final.shape[0]) + " rows and " + 
        str(dfk_final.shape[1]) + " columns.")
    return dfk_final

def parse_dfk():
    """
    """
    return rdfk

def main():
    """
    """
    dfk_final = make_dfk()
new_fname = "sitest_results_blood-groups-2platforms"
new_fname = new_fname + ".csv"
dfk_final.to_csv(new_fname)

    return dfk_final

if __name__ == "__main__":
    """
    """
    main()


