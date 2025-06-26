'''
1. Add arrival times of stations for events in the cluster
2. Only include common stations shared by all events within the same cluster
'''
import os,glob,time
from obspy.core import read, UTCDateTime
import numpy as np

def dtime2str(dtime):
    date = ''.join(str(dtime).split('T')[0].split('-'))
    time = ''.join(str(dtime).split('T')[1].split(':'))[0:9]
    return date + time

# read phase file (full format)
def read_fpha(fpha):
    f=open(fpha); lines=f.readlines(); f.close()
    event_list = []
    dtype = [('evid','int'),('evt_name','O'),('loc','O'),('picks','O')]
    for line in lines:
        codes = line.split(',')
        if len(codes[0])>=14:
            ot = UTCDateTime(codes[0])
            event_name = dtime2str(UTCDateTime(ot))
            lat, lon, dep, mag = [float(code) for code in codes[1:5]]
            event_loc = [ot, lat, lon, dep, mag]
            evid = int(codes[-1][:-1])
            event_list.append((evid, event_name, event_loc, {}))
        else:
            net_sta = codes[0]
            tp = UTCDateTime(codes[1]) if codes[1]!='-1' else -1
            ts = UTCDateTime(codes[2]) if codes[2][:-1]!='-1' else -1
            event_list[-1][-1][net_sta] = [tp, ts]
    return np.array(event_list, dtype=dtype)

def read_cluster_file(cluster_file):
    cluster_dict = {}
    with open(cluster_file) as f:
        for line in f:
            codes = line.strip().split()
            if not codes:
                continue
            if codes[0] == '#':
                cluster_id = codes[-1]
                if cluster_id not in cluster_dict:
                    cluster_dict[cluster_id] = []
            else:
                parts = line.strip().split(',')
                evid = int(parts[-1])
                cluster_dict[cluster_id].append(evid)
    return cluster_dict

def build_evid_index(pha_list):
    evid_dict = {}
    for ev in pha_list:
        evid_dict[ev['evid']] = ev
    return evid_dict

def write_cluster_phase_file(cluster_dict, evid_dict, output_file):
    with open(output_file, 'w') as f:
        for clust_id, evids in cluster_dict.items():
            f.write(f"# cluster {clust_id}\n")

            all_station_sets = []
            for evid in evids:
                if evid not in evid_dict:
                    continue
                picks = evid_dict[evid]['picks']
                station_list = list(picks.keys())
                all_station_sets.append(set(station_list))

            #search common_stations
            if not all_station_sets:
                continue
            common_stations = set.intersection(*all_station_sets)

            for evid in evids:
                if evid not in evid_dict:
                    continue
                ot, lat, lon, dep, mag = evid_dict[evid]['loc']
                f.write(f"{ot.isoformat()},{lat},{lon},{dep},{mag},{evid}\n")
                # sort 
                sort_list = []
                for net_sta in common_stations:
                    tp, ts = evid_dict[evid]['picks'][net_sta]
                    if tp == -1 or ts == -1:
                        continue
                    delta = (ts - tp)
                    sort_list.append((delta, net_sta, tp, ts))
                # sort by (ts - tp) 
                sort_list.sort()
                for delta, net_sta, tp, ts in sort_list:
                    tp_str = tp.isoformat() if tp != -1 else '-1'
                    ts_str = ts.isoformat() if ts != -1 else '-1'
                    f.write(f"{net_sta},{tp_str},{ts_str}\n")


fpha = 'input/eg_cc_full.pha'
cluster_file = 'output/eg_rep.clust'
output_file = 'output/eg_clustered.pha'

pha_list = read_fpha(fpha)
cluster_dict = read_cluster_file(cluster_file)
evid_dict = build_evid_index(pha_list)
write_cluster_phase_file(cluster_dict, evid_dict, output_file)
print(f"output：{output_file}")
