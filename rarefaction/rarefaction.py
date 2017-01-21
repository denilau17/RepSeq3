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
    T = sp.cluster.hierarchy.fcluster(Z,cutoff,criterion='distance')

    return T
    
def main(db, subject, outfile, iters):

    connection = sqlite3.connect(db)
    c = connection.cursor()
    connection.text_factory = str

    query = "SELECT * FROM " + subject + ";"
    results = c.execute(query).fetchall()

    connection.close()
    
    np.random.seed(10)

    clone_counts = []
    for i in range(iters):
        random = np.random.choice(len(results), 54320)
    
        subset = [results[j] for j in random]
        subgroup_list = set([x[-1] for x in subset])

        clones = []
        for subgroup in subgroup_list:
            seqs = []
            max_edit =  int(subgroup.split('_')[-1])*.15
            print max_edit
            for row in subset:
                if row[-1] == subgroup:
                    seqs.append(row[7])
            clusters = cluster_seqs(seqs, int(max_edit))
            for cluster in clusters:
                name = subgroup + str(cluster)
                clones.append(name)
        counter=collections.Counter(clones)
        clone_counts.append(counter.values())
    print clone_counts

    out = open(outfile, 'wb')
    csv_out = csv.writer(out)
    csv_out.writerows(clone_counts)

if __name__ == "__main__":
    #db = '/Users/denise/Documents/RepSeq2/IMGT_parsed3.sqlite'
    #subject = 'IMGT_011'
    #outfile = '/Users/denise/Documents/RepSeq2/RepSeq/vdj/testshit.csv'

    db = sys.argv[1]
    subject = sys.argv[2]
    outfile = sys.argv[3]
    
    main(db, subject, outfile, 100)

