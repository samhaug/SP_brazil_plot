#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : sanity_check.py
Purpose : make sure amplitudes match after rotating
Creation Date : 05-04-2017
Last Modified : Thu 06 Apr 2017 11:06:10 AM EDT
Created By : Samuel M. Haugland

==============================================================================
'''

import numpy as np
from matplotlib import pyplot as plt
from subprocess import call
from os import listdir
import h5py
import obspy
from obspy.taup import TauPyModel
import seispy

def main():
    svr = svaxi_read('/home/samhaug/work1/SP_brazil_sims/SVaxi/0327_shareseismos/shareseismos/smslab_a-10_h10_dVs10/')
    svr = seispy.filter.range_filter(svr,(56,72))
    stz,str = data_read()
    svr,svz,str,stz = normalize_all(svr,str,stz)
    #seispy.plot.simple_section(svr)
    #seispy.plot.simple_section(str)
    #seispy.plot.simple_section(stz)
    print svr[6].stats.sac['gcarc']
    print str[3].stats.sac['gcarc']
    #plt.plot(svr[6].data)
    plt.plot(svz[6].data)
    #plt.plot(str[1].data)
    plt.plot(stz[3].data)
    #seispy.plot.simple_section(stz)
    plt.show()

def svaxi_read(dir):
    str = obspy.read(dir+'*0.R')
    str.interpolate(50)
    for idx,tr in enumerate(str):
        tr.stats.sac['o'] += -202
        str[idx] = seispy.data.phase_window(str[idx],phase=['P'],window=(-10,600))
    return str

def data_read():
    stz = obspy.read('/home/samhaug/work1/SP_brazil_data/sparse_Z.pk')
    str = obspy.read('/home/samhaug/work1/SP_brazil_data/sparse_R.pk')
    #stz.filter('highpass',freq=1/15.)
    #str.filter('highpass',freq=1/15.)
    stz.filter('lowpass',freq=0.5)
    str.filter('lowpass',freq=0.5)
    str.interpolate(50)
    stz.interpolate(50)
    for idx,tr in enumerate(str):
        str[idx] = seispy.data.phase_window(str[idx],phase=['P'],window=(-10,600))
        stz[idx] = seispy.data.phase_window(stz[idx],phase=['P'],window=(-10,600))
    return stz,str

def normalize_all(svr,str,stz):
    model = TauPyModel(model='prem50')
    svr.normalize()
    for idx,tr in enumerate(str):
        #a = seispy.data.phase_window(tr,phase=['S'],window=(-30,30))
        #plt.plot(a.data)
        #plt.show()
        norm = np.abs(str[idx].data).max()
        str[idx].data *= 1./norm
        stz[idx].data *= 1./norm
    svz = svr.copy()
    for idx,tr in enumerate(svz):
        arrivals = model.get_travel_times(source_depth_in_km=svz[idx].stats.sac['evdp'],
                                          distance_in_degree=svz[idx].stats.sac['gcarc'],
                                          phase_list=['S'])
        s_ang = arrivals[0].incident_angle
        svz[idx].data *= 1./np.tan(np.radians(s_ang))

    return svr,svz,str,stz

def radial_rotate(stn,ste):
    model = TauPyModel(model='prem50')
    r,t = seispy.rotate.rotate_ne_rt(stn,ste)
    seispy.plot.plot(r[5],phase_list=['P','S'])
    seispy.plot.plot(t[5],phase_list=['P','S'])
    plt.plot(r[5].data,alpha=0.5)
    plt.plot(t[5].data,alpha=0.5)
    plt.show()
    s = r.copy()
    p = r.copy()
    for idx, tr in enumerate(r):
        arrivals = model.get_travel_times(source_depth_in_km=r[idx].stats.sac['evdp'],
                                          distance_in_degree=r[idx].stats.sac['gcarc'],
                                          phase_list=['S','S1800P'])
        s_ang = arrivals[0].incident_angle
        conv_ang = arrivals[1].incident_angle

        s[idx].data *= 1./np.cos(np.radians(s_ang))
        a = seispy.data.phase_window(s[idx],phase=['S'],window=(-30,30))
        p[idx].data *= 1./np.max(np.abs(a.data))
        p[idx].data *= 1./np.sin(np.radians(conv_ang))
    return p

def rotate_conversion(stz,ste,stn):
    l,qf,tf = seispy.data.rotate_phase(stz,stn,ste,['S1800P'])
    lf,q,t = seispy.data.rotate_phase(stz,stn,ste,['S'])
    for idx,tr in enumerate(l):
        q[idx].data = obspy.signal.filter.envelope(q[idx].data)
        s_tr = seispy.data.phase_window(q[idx],['S'],window=(-50,50))
        l[idx].data *= 1./np.max(np.abs(s_tr.data))
    return l

main()





