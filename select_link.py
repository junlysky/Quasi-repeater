""" Select cc_dd links with different CC threshold
"""
import numpy as np
from dataset_cc import read_fsta, read_fpha_dict, calc_dist_km
import config

cfg = config.Config()
# i/o paths
flink = 'output/cc_dd_link.csv' 
fout_rep_link = 'output/cc_dd_link_rep.npy'
fout_qrep_link = 'output/cc_dd_link_qrep.npy'
fpha = cfg.fpha
fsta = cfg.fsta
event_dict = read_fpha_dict(fpha)
sta_dict = read_fsta(fsta)
# thres for linking event pairs
dt_thres = cfg.dt_thres[1][2] # max dt
num_sta_thres = cfg.num_sta_thres[1] # min sta
loc_dev_thres = cfg.loc_dev_thres[1] # max dev loc
dep_dev_thres = cfg.dep_dev_thres[1] # max dev dep
dist_thres = cfg.dist_thres[1] # max epi-dist

def select_link(flink, cc_thres):
    link_list = []
    f=open(flink); lines=f.readlines(); f.close()
    for i,line in enumerate(lines):
      codes = line.split(',')
      if len(codes)==2:
        to_add = True
        evid1, evid2 = [str(int(code)) for code in codes]
        if evid1 not in event_dict or evid2 not in event_dict:
            to_add = False; continue
        lat1, lon1, dep1 = event_dict[evid1][0][0:3]
        lat2, lon2, dep2 = event_dict[evid2][0][0:3]
        # 1. select loc dev
        loc_dev = calc_dist_km([lat1,lat2], [lon1,lon2])
        dep_dev = abs(dep1 - dep2)
        if not (loc_dev<loc_dev_thres and dep_dev<dep_dev_thres):
            to_add = False; continue
        link_list.append([[evid1, evid2],[]])
      else:
        if not to_add: continue
        # 2. select by epicentral distance
        sta = codes[0]
        sta_lat, sta_lon = sta_dict[sta]
        dist1 = calc_dist_km([sta_lat,lat1], [sta_lon,lon1])
        dist2 = calc_dist_km([sta_lat,lat2], [sta_lon,lon2])
        if min(dist1, dist2)>dist_thres: continue
        # select by CC
        dt_sp, _,_, cc_det, cc_p, cc_s = [abs(float(code)) for code in codes[1:7]]
        is_self = 0
        if cc_det==1 and cc_p==1 and cc_s==1 and dt_sp==0: is_self=1
        if dt_sp>dt_thres: continue
        if min(cc_p, cc_s, cc_det)<cc_thres: continue
        link_list[-1][-1].append([line,is_self])
    # select links
    link_list_sel = []
    for [[evid1,evid2],cd_list] in link_list:
        if len(cd_list)<num_sta_thres: continue
        if sum([is_self for [_,is_self] in cd_list])>=len(cd_list)-1: continue
        link_list_sel.append([evid1,evid2])
    return link_list_sel

# save refined flink
print('selecting quasi-repeater links')
qrep_link_list = select_link(flink, cfg.cc_thres[0])
print('selecting repeater links')
rep_link_list = select_link(flink, cfg.cc_thres[1])
np.save(fout_qrep_link, qrep_link_list)
np.save(fout_rep_link, rep_link_list)
