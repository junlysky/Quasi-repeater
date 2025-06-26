"""
plot waves for common stations in all events from clusters
"""
import os,sys,glob
import numpy as np
from obspy.core import read, UTCDateTime
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

#may be not to use,or to verfity cc 
def calc_cc(data, temp):
    ntemp, ndata = len(temp), len(data)
    if ntemp>ndata: return [0]
    nndata = len(data)
    norm_temp = np.sqrt(np.sum(temp**2))
    norm_data = np.sqrt(np.sum(data**2))
    cc= np.correlate(data, temp,mode='same')
    if not ((np.sqrt(np.sum(data**2)))*(np.sqrt(np.sum(temp**2))) == 0):
       cc = cc / ((np.sqrt(np.sum(data**2)))*(np.sqrt(np.sum(temp**2))))
       time_shift = np.argmax(cc) - int((len(temp) - 1) / 2) # mode = same
    cc[np.isinf(cc)] = 0.
    cc[np.isnan(cc)] = 0.
    cc_max = cc.max()
    return round(cc_max, 2), time_shift 

def preprocess_stream(stream, freqmin=1, freqmax=16):
    stream = stream.detrend('demean').detrend('linear') #.taper(max_percentage=0.05, max_length=5.)
    #stream = stream.filter('highpass', freq=1.)
    stream = stream.filter('bandpass', freqmin=freqmin,freqmax=freqmax)
    return stream.normalize()

#keep event name keep same
def dtime2str(dtime):
    date = ''.join(str(dtime).split('T')[0].split('-'))
    time = ''.join(str(dtime).split('T')[1].split(':'))[0:9]
    return date + time


def read_cluster_phase_file(pha_file):
    clusters = {}
    current_cluster = None
    current_event = None
    with open(pha_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("# cluster"):
                current_cluster = line.split()[-1]
                clusters[current_cluster] = []
            elif line.startswith("20"):
                parts = line.split(',')
                ot,lat,lon,dep,mag,evid = parts
                ot = UTCDateTime(ot)
                event_name = dtime2str(ot)
                ev_info = {
                    'evid': evid,
                    'event_name': event_name,
                    'origin_time': ot,
                    'lat':lat,
                    'lon':lon,
                    'dep':dep,
                    'mag':mag,
                    'picks': {}
                }
                current_event = ev_info
                clusters[current_cluster].append(current_event)
            else:
                parts = line.split(',')
                net_sta = parts[0]
                tp = UTCDateTime(parts[1]) if parts[1] != '-1' else -1
                ts = UTCDateTime(parts[2]) if parts[2] != '-1' else -1
                current_event['picks'][net_sta] = (tp, ts)
    return clusters

def load_and_preprocess_waveform(event, stations, data_path):
    out = {}
    event_path = os.path.join(data_path, event["event_name"])
    sac_files = glob.glob(os.path.join(event_path, "*.*"))
    for sta in stations:
        sta_files = [f for f in sac_files if sta in f]
        st = read(sta_files[0])
        st+= read(sta_files[1])
        st+= read(sta_files[2])
        out[sta] = st
    return out


def plot_cluster_waveforms(cluster_id, cluster_events, data_path):
    all_stations = [set(ev['picks'].keys()) for ev in cluster_events]
    common_stations = sorted(set.intersection(*all_stations))
    #验证波形路径中确实存在这些台站
    valid_common_stations = []
    for sta in common_stations:
        missing = False
        for ev in cluster_events:
            event_path = os.path.join(data_path, ev["event_name"])
            sac_files = glob.glob(os.path.join(event_path, "*.*"))
            # 检查该台站是否存在于当前事件的 SAC 文件中
            match = any(sta in os.path.basename(f) for f in sac_files)
            if not match:
               missing = True
               break
        if not missing:
           valid_common_stations.append(sta)
    common_stations = sorted(valid_common_stations)

    lw, hp = 1, 16
    template_event = cluster_events[0]
    st_template_all = load_and_preprocess_waveform(template_event, common_stations, data_path)
    new_samp_rate = 500
 
    fig, axes = plt.subplots(len(cluster_events) * len(common_stations),1,figsize=(12,10),sharex=True)
    for i_sta, sta in enumerate(common_stations):
        for i_ev, ev in enumerate(cluster_events):
            st_event_all = load_and_preprocess_waveform(ev, common_stations, data_path)
            ax_index = i_sta * len(cluster_events) + i_ev  #  索引顺序
            ax = axes[ax_index]

            st_event = st_event_all.get(sta)
            st_event = preprocess_stream(st_event,lw, hp)
            st_temp = st_template_all.get(sta)
            st_temp  = preprocess_stream(st_temp,lw, hp)
            samp_rate = st_event[0].stats.sampling_rate

            data = st_event[0].data[900:3000]
            temp = st_temp[0].data[900:3000]
            time = np.arange(0, len(data) / samp_rate, 1 / samp_rate)

            #interp_data = interp1d(time, data, kind="linear", fill_value="extrapolate")(new_time)
            #interp_temp = interp1d(time, temp, kind="linear", fill_value="extrapolate")(new_time)
            #new_time = np.arange(0, len(data) / samp_rate, 1 / new_samp_rate)

            cc, time_shift = calc_cc(data, temp)
            #cc, time_shift = calc_cc(interp_data, interp_temp)
            ax.plot(time, data,color='red',alpha=1.0, linewidth=0.5)
            #ax.plot(new_time, interp_data, color="black")
            ax.text(0.75, 0.85, f"{ev['event_name']} {sta}", transform=ax.transAxes, fontsize=8)
            #ax.text(0.8, 0.9, f"CC={cc:.2f}", transform=ax.transAxes, fontsize=10)
            ax.set_ylim([-1.3, 1.3])
            ax.set_yticks([])
            ax.set_ylabel('')         # 去除纵轴标签

    for j in range(len(axes)):
      for spine in axes[j].spines:
        if spine in ['left', 'right']:
            axes[j].spines[spine].set_visible(True)
        elif spine == 'top':
            axes[j].spines[spine].set_visible(j == 0)  # 仅首子图保留上边框
        elif spine == 'bottom':
            axes[j].spines[spine].set_visible(j == len(axes) - 1)  # 仅末子图保留下边框


    axes[-1].set_xlabel("Time (s)",fontsize=14)
    #plt.suptitle("CC>0.95, $\Delta$(S-P)<0.01 s", fontsize=14)
    axes[0].set_title("2024 Noto Repeaters, CC>0.95, $\Delta$(S-P)<0.01 s", fontsize=14)
    plt.subplots_adjust(hspace=0)  # 设置子图间距为0       
    plt.savefig(f"Cluster_{cluster_id}_waveforms.png")
    plt.show()
    plt.close()

data_path = '/data_path/eg_events_repeat'
cluster_dict = read_cluster_phase_file('output/eg_clustered.pha') 
for cluster_id in cluster_dict:
    #cluster_id = str(0)
    cluster_events = cluster_dict[cluster_id]




    plot_cluster_waveforms(cluster_id, cluster_events, data_path)    
