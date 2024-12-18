""" Initial clustering for repeater detection
"""
import os, shutil, glob
import numpy as np
from obspy import UTCDateTime
from dataset_cc import read_fpha_dict
import config

# i/o paths
cfg = config.Config()
fpha = cfg.fpha
event_dict = read_fpha_dict(fpha)
flink = 'output/cc_dd_link_rep.npy'
fout = open('output/%s_rep.clust'%cfg.ctlg_code,'w')
# linkage criteria
num_nbr_thres = cfg.num_nbr_thres[0]

# 1. get evid list
link_list = np.load(flink)
evid_list = np.unique(link_list)
evid_dict = {}
for ii,evid in enumerate(evid_list): evid_dict[evid] = ii

# 2. build link mat
num_events = len(evid_list)
link_mat = np.zeros([num_events, num_events])
for link in link_list:
    ii, jj = np.sort([evid_dict[evid] for evid in link])
    link_mat[ii,jj] = 1

def n_link_clustering(link_mat, num_link_thres):
    # select by num_nbr
    for ii in range(len(link_mat)):
        num_nbr = sum(link_mat[ii] + link_mat[:,ii])
        if num_nbr<num_link_thres:
            link_mat[ii, link_mat[ii]==1] = 0
            link_mat[link_mat[:,ii]==1, ii] = 0
    # single-link clustering
    clusters = []
    for i in range(len(link_mat)-1):
        nbrs  = list(np.where(link_mat[i]==1)[0])
        nbrs += list(np.where(link_mat[:,i]==1)[0])
        link_mat[i, link_mat[i]==1] = 0
        link_mat[link_mat[:,i]==1, i] = 0
        if len(nbrs)>0: clusters.append([i]+nbrs) # save the evid
        while len(nbrs)>0:
            new_nbrs = []
            for nbr in nbrs:
                new_nbrs += list(np.where(link_mat[nbr]==1)[0])
                new_nbrs += list(np.where(link_mat[:,nbr]==1)[0])
                link_mat[nbr, link_mat[nbr]==1] = 0
                link_mat[link_mat[:,nbr]==1, nbr] = 0
            clusters[-1] += new_nbrs
            nbrs = new_nbrs
    return [np.unique(cluster) for cluster in clusters]

# 3. n-link clustering
clusters = n_link_clustering(link_mat, num_nbr_thres)
print('%s clusters found'%len(clusters))

for i,cluster in enumerate(clusters):
    print('write %sth cluster'%i)
    fout.write('# cluster %s \n'%i)
    for idx in cluster: 
        fout.write(event_dict[evid_list[idx]][1])
fout.close()
