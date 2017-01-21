import sys
import os
import csv
import sqlite3

import numpy as np
import scipy as sp
import scipy.cluster

import multiprocessing as mp
import itertools

import Levenshtein
import collections

#generate distance matrix
def pdist(X):
    m = len(X)
    dm = np.zeros((m * (m - 1) / 2,), dtype=np.double)
    k = 0
    for i in xrange(0, m - 1):
        for j in xrange(i+1, m):
            dm[k] = Levenshtein.distance(X[i], X[j])
            k += 1
    return dm

#cluster seqs
def cluster_seqs(seqs,cutoff,linkage='single'):
    if len(seqs) == 0:
        return (np.array([]),{})

    #checks if there is only 1 unique seq
    if len(seqs) == 1:
        T = np.array([1]*len(seqs))
        return T

    #compute distance matrix
    Y = pdist(seqs)
    
    #compute linkage
    Z = sp.cluster.hierarchy.linkage(Y,method=linkage)

    # determine the clusters at level cutoff
    #returns array where array[i] is the cluster number that seqs[i] belongs to
    T = sp.cluster.hierarchy.fcluster(Z,cutoff,criterion='distance')

    return T
    
def main(db, subject, outfile):

    connection = sqlite3.connect(db)
    c = connection.cursor()
    connection.text_factory = str

    #get data from database
    query = "SELECT * FROM " + subject + ";"
    results = c.execute(query).fetchall()

    connection.close()

    #iterate through each subgroup in the list
    out = []
    subgroup_list = list(set([x[-1] for x in results]))
    subgroup_list = ['IGHV3-43*01_IGHJ1*01_16']
    
    
    for subgroup in subgroup_list:
        subgroup_data = []
        
        #get cutoff distance for clustering
        max_edit =  int(subgroup.split('_')[-1])*.15
        
        #get all the sequences in the subgroup
        for row in results:
            if row[-1] == subgroup:
                subgroup_data.append(list(row))

        #cluster sequences
        seqs = [x[7] for x in subgroup_data]        
        clusters = cluster_seqs(seqs, int(max_edit))
        clusters = [x for x in clusters]

        #format data for output
        for i in range(len(subgroup_data)):
            j = list(subgroup_data[i])
            k = clusters[i]
            j.append(k)

            out.append(j)
    
    #write to csv
    of = open(outfile, 'wb')
    csv_out = csv.writer(of)
    csv_out.writerows(out)
    
    
if __name__ == "__main__":
    #db = '/home/djlau/RepSeq3/IMGT_parsed.sqlite'
    #subject = 'IMGT_012'
    #outfile = '/home/djlau/RepSeq3/testshit.csv'

    db = sys.argv[1]
    subject = sys.argv[2]
    outfile = sys.argv[3]
    
    main(db, subject, outfile)

