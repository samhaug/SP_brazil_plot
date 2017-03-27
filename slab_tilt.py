#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : slab_tilt.py
Purpose : Compare SxP amplitudes from different tilting slabs.
Creation Date : 27-03-2017
Last Modified : Mon 27 Mar 2017 06:39:09 PM EDT
Created By : Samuel M. Haugland

==============================================================================
'''

import numpy as np
from matplotlib import pyplot as plt
from subprocess import call
from os import listdir
import h5py
import obspy
import seispy

def main():
    svdir = '/home/samhaug/work1/SP_brazil_sims/SVaxi/0327_shareseismos/shareseismos/'
    full = listdir(svdir)
    part = [i for i in full if 'h10_dVs10' in i]

    for dir in part[::5]:
        print dir
        st = setup_stream(svdir+dir)

def setup_stream(dir):
    st = obspy.read(dir+'/*0.Z')
    for tr in st:
        tr.stats.sac['o'] += -202.
    st = seispy.data.normalize_on_phase(st,phase=['S'])
    return st

def strip_conversion(st):
    for idx,tr in enumerate(st):
        st[idx] = seispy.data.phase_window(st[idx],phase=['S'],window=(-20,20))
    return st_strip

def stats_conversion(st_strip):
    conv = []
    for tr in st_strip:
        conv.append(tr.data)
    std = np.std(conv,axis=0)
    mean = np.mean(conv,axis=0)
    return mean,std

def plot_conversion(angle,mean,std,ax):
    ax.errorbar(angle,mean,yerr=std,color='k')

def plot_amplitude(amp,ax):
    print 'fuck'

def setup_figure():
    fig,ax_list = plt.subplots(2,1,figsize=(3,4.5))
    return fig,ax_list


main()
