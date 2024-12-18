import os, glob
import numpy as np
from obspy import read, UTCDateTime
import warnings
warnings.filterwarnings("ignore")

# i/o paths
fclust = '../output/eg_rep-qrep.clust'
fout = open('../output/eg_rep-qrep_rel-mag.clust','w')
event_root = '/data/bigdata/eg_events_data'
fpha = '../output/eg_rep-qrep_refine-mag_full.pha'
fsta = '../input/station_eg.csv'
# selection
win_amp = [1,10]
dist_max = 100
win_sta, win_lta = 0.8, 2.5
snr_thres = 8

def dtime2str(dtime):
    date = ''.join(str(dtime).split('T')[0].split('-'))
    time = ''.join(str(dtime).split('T')[1].split(':'))[0:9]
    return date + time

def read_fpha(fpha):
    event_dict = {}
    f=open(fpha); lines=f.readlines(); f.close()
    for line in lines:
        codes = line.split(',')
        if len(codes[0])>10:
            ot = UTCDateTime(codes[0])
            lat, lon, dep, mag = [float(code) for code in codes[1:5]]
            evid = codes[-1][:-1]
            event_loc = [ot, lat, lon, dep, mag]
            event_dict[evid] = [event_loc, {}]
        else:
            net_sta = codes[0]
            tp, ts = [UTCDateTime(code) for code in codes[1:3]]
            event_dict[evid][-1][net_sta] = [tp, ts, line]
    return event_dict

def read_fsta(sta_file):
    sta_dict = {}
    f=open(sta_file); lines=f.readlines(); f.close()
    for line in lines:
        codes = line.split(',')
        net_sta = codes[0]
        lat, lon, ele = [float(code) for code in codes[1:4]]
        # format 1: same gain for 3-chn & time invariant
        if len(codes[4:])==1: gain = float(codes[4])
        # format 2: different gain for 3-chn & time invariant
        elif len(codes[4:])==3: gain = [float(code) for code in codes[4:]]
        if net_sta not in sta_dict:
            sta_dict[net_sta] = [lat,lon,ele,gain]
        else:
            sta_dict[net_sta][-1].append(gain[0]) # if format 3
    return sta_dict

# read fclust
clusts = []
f=open(fclust); lines=f.readlines(); f.close()
for line in lines:
    if line[0:2]=='# ': 
        clusts.append([[line,[]],[],[]]) 
        idx=0; continue
    if line[0:4]=='## C': 
        clusts[-1][1] += [line,[]]
        idx=1; continue
    if line[0:4]=='## a': 
        clusts[-1][2] += [line,[]]
        idx=2; continue
    evid = line.split(',')[-1][:-1]
    clusts[-1][idx][1].append(evid)

event_dict = read_fpha(fpha)
sta_dict = read_fsta(fsta)

def calc_sta_lta(data, samp_rate):
    npts = len(data)
    win_lta_npts, win_sta_npts = int(win_lta*samp_rate), int(win_sta*samp_rate)
    if npts < win_lta_npts + win_sta_npts:
        print('input data too short!')
        return np.zeros(1)
    sta = np.zeros(npts)
    lta = np.ones(npts)
    data_cum = np.cumsum(data)
    sta[:-win_sta_npts] = data_cum[win_sta_npts:] - data_cum[:-win_sta_npts]
    sta /= win_sta_npts
    lta[win_lta_npts:]  = data_cum[win_lta_npts:] - data_cum[:-win_lta_npts]
    lta /= win_lta_npts
    sta_lta = sta/lta
    sta_lta[0:win_lta_npts] = 0.
    sta_lta[np.isinf(sta_lta)] = 0.
    sta_lta[np.isnan(sta_lta)] = 0.
    return sta_lta

def get_s_amp(velo, samp_rate):
    # remove mean
    velo -= np.reshape(np.mean(velo, axis=1), [velo.shape[0],1])
    # velocity to displacement
    disp = np.cumsum(velo, axis=1)
    disp /= samp_rate
    return np.amax(abs(np.sum(disp**2, axis=0)))**0.5

def calc_dist_km(lat, lon, dep):
    cos_lat = np.cos(np.mean(lat) * np.pi / 180)
    dx = cos_lat * (lon[1] - lon[0]) * 111
    dy = (lat[1] - lat[0]) * 111
    dz = dep[1] - dep[0]
    return (dx**2 + dy**2)**0.5, (dx**2 + dy**2 + dz**2)**0.5

def get_amp_dict(evid):
    amp_dict = {}
    event_loc, pha_dict = event_dict[evid]
    ot, lat, lon, dep, mag = event_loc
    event_name = dtime2str(ot)
    event_dir = os.path.join(event_root, event_name)
    for net_sta, [tp, ts, pha_line] in pha_dict.items():
        sta_lat, sta_lon, sta_ele, gain = sta_dict[net_sta]
        dist_epi, dist_hyp = calc_dist_km([lat,sta_lat],[lon,sta_lon],[dep,-sta_ele/1e3])
        if dist_epi>dist_max: continue
        st_paths = sorted(glob.glob(os.path.join(event_dir, "%s.*"%net_sta)))
        if len(st_paths)<3: continue
        st  = read(st_paths[0])
        st += read(st_paths[1])
        st += read(st_paths[2])
        samp_rate = st[0].stats.sampling_rate
        # check SNR
        sta_lta = calc_sta_lta(st[2].data**2, samp_rate)
        if np.amax(sta_lta)<snr_thres: continue
        # calc amp
        st_amp = st.slice(tp-win_amp[0], ts+win_amp[1])
        if len(st_amp)<3: continue
        data_amp = np.array([tr.data / gain[ii] * 1e6 for ii,tr in enumerate(st_amp)])
        s_amp = get_s_amp(data_amp, samp_rate)
        amp_dict[net_sta] = s_amp
    return amp_dict

def get_amp_ratio(amp_dict1, amp_dict2):
    amp_ratio = []
    sta_list = [sta for sta in amp_dict1.keys() if sta in list(amp_dict2.keys())]
    for sta in sta_list:
        amp_ratio.append(amp_dict1[sta] / amp_dict2[sta])
    return amp_ratio

for [evid_list1, evid_list2, evid_list3] in clusts:
    print('process %s'%evid_list1[0][:-1])
    print(evid_list1, evid_list2, evid_list3)
    # get repeater mag
    head_line1, rep_evid_list = evid_list1
    mag_list = np.array([event_dict[evid][0][-1] for evid in rep_evid_list])
    mag_med = np.sort(mag_list)[int(len(mag_list)/2)]
    med_evid = rep_evid_list[np.where(mag_list==mag_med)[0][0]]
    # get ref amp
    ref_amp_dict = get_amp_dict(med_evid)
    # calc rel mag
    fout.write(head_line1)
    for evid in rep_evid_list:
        ot, lat, lon, dep, mag = event_dict[evid][0]
        if evid==med_evid: 
            fout.write('%s,%s,%s,%s,%s,%s\n'%(ot,lat,lon,dep,mag,evid)); continue
        amp_dict_i = get_amp_dict(evid)
        amp_ratio = get_amp_ratio(amp_dict_i, ref_amp_dict)
        mag = round(np.median([np.log10(ai) + mag_med for ai in amp_ratio]),2)
        fout.write('%s,%s,%s,%s,%s,%s\n'%(ot,lat,lon,dep,mag,evid))
    # for additional events
    if len(evid_list2)>0:
        head_line2, add_evid_list2 = evid_list2
        fout.write(head_line2)
        for evid in add_evid_list2:
            ot, lat, lon, dep, mag = event_dict[evid][0]
            amp_dict_i = get_amp_dict(evid)
            amp_ratio = get_amp_ratio(amp_dict_i, ref_amp_dict)
            mag = round(np.median([np.log10(ai) + mag_med for ai in amp_ratio]),2)
            fout.write('%s,%s,%s,%s,%s,%s\n'%(ot,lat,lon,dep,mag,evid))
    if len(evid_list3)>0:
        head_line3, add_evid_list3 = evid_list3
        fout.write(head_line3)
        for evid in add_evid_list3:
            ot, lat, lon, dep, mag = event_dict[evid][0]
            amp_dict_i = get_amp_dict(evid)
            amp_ratio = get_amp_ratio(amp_dict_i, ref_amp_dict)
            mag = round(np.median([np.log10(ai) + mag_med for ai in amp_ratio]),2)
            fout.write('%s,%s,%s,%s,%s,%s\n'%(ot,lat,lon,dep,mag,evid))
fout.close()
