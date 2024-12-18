""" Relink with lower CC threshold for quasi-repeater detection
"""
import numpy as np
from dataset_cc import read_fsta, read_fpha_dict, calc_dist_km
import config

# i/o paths
cfg = config.Config()
ctlg_code = cfg.ctlg_code
fpha = cfg.fpha
event_dict = read_fpha_dict(fpha)
fsta = cfg.fsta
sta_dict = read_fsta(fsta)
flink_qrep = 'output/cc_dd_link_qrep.npy'
link_list = np.load(flink_qrep)
fclust_rep = 'output/%s_rep.clust'%cfg.ctlg_code
fout = open('output/%s_rep-qrep.clust'%cfg.ctlg_code,'w')
# linkage criteria
num_nbr_thres = cfg.num_nbr_thres[0]

# read fclust repeater
def get_recall(ref_list, links):
    det_dict = {}
    for link in links:
        if len(np.intersect1d(link, ref_list))>0:
            if link[0] not in det_dict: det_dict[link[0]] = [link[1]]
            else: det_dict[link[0]].append(link[1])
            if link[1] not in det_dict: det_dict[link[1]] = [link[0]]
            else: det_dict[link[1]].append(link[0])
    # apply num_nbr criteria
    recall_list, new_list = [],[]
    for evid, linked_ids in det_dict.items():
        if len(np.unique(linked_ids))<num_nbr_thres: continue
        if evid in ref_list: recall_list.append(evid)
        else: new_list.append(evid)
    #miss_list = np.setdiff1d(ref_list, recall_list)
    return recall_list, new_list

print('reading fclust_rep')
clusts = []
f=open(fclust_rep); lines=f.readlines(); f.close()
for line in lines:
    if line[0]=='#': clusts.append([]); continue
    evid = line.split(',')[-1][:-1]
    clusts[-1].append(evid)

print('detecting quasi-repeaters within repeater sequences')
for ii,clust in enumerate(clusts):
    recall_list, new_list = get_recall(clust, link_list)
    fout.write('# cluster %s\n'%ii)
    for evid in recall_list: fout.write(event_dict[evid][1])
    fout.write('## quasi-repeaters\n')
    for evid in new_list: fout.write(event_dict[evid][1])
fout.close()
