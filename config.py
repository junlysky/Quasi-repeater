""" Configure file for repeater and quasi-repeater detection
"""

class Config(object):
  def __init__(self):

    # i/o paths
    self.hypo_root = '/home/zhouyj/bin'
    self.ctlg_code = 'eg'
    self.fsta = 'input/station_eg.csv'
    self.fpha = 'input/phase_eg.csv'
    self.event_root = '/data/bigdata/eg_events_data'
    self.data_dir = '/data/Continuous_data_eg'
    # thresholds for event pair link
    self.cc_thres = [0.5,0.9] # CC thres for event pair
    self.dt_thres = [[0.5,0.8,0.04],[0.5,0.8,0.01]] # dt_p, dt_s, d(S-P)
    self.loc_dev_thres = [5,5] # km, maximum epicenter separation
    self.dep_dev_thres = [10,10] # km, maximum dep separation
    self.dist_thres = [120,100] # km, max epicentral dist
    self.num_sta_thres = [3,4] # min sta_num for one event pair
    self.num_nbr_thres = [2,200]
    self.temp_mag = 1.    # min mag for templates
    self.temp_sta = 4    # min sta_num for templates
    # data prep
    self.num_workers = 10
    self.freq_band = [1,16]
    self.samp_rate = 100
    self.chn_p = [[2],[0,1,2]][1] # chn for P picking
    self.chn_s = [[0,1],[0,1,2]][1] # chn for S picking
    self.win_event = [10, 30]    # event data cutting, just long enough
    self.win_temp_det = [1.,11.]
    self.win_temp_p = [0.5,2.]
    self.win_temp_s = [0.2,3.8]
    self.ot_range = '20200101-20230501'
    self.lat_range = [35.4,40]
    self.lon_range = [34.5,42]
