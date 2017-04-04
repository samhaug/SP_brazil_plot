#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : slab_tilt.py
Purpose : Compare SxP amplitudes from different tilting slabs.
Creation Date : 27-03-2017
Last Modified : Tue 28 Mar 2017 02:06:48 PM EDT
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
from obspy.taup import TauPyModel
model = TauPyModel(model='prem')

def main():
    svdir = '/home/samhaug/work1/SP_brazil_sims/SVaxi/0327_shareseismos/shareseismos/'
    full = listdir(svdir)
    part = [i for i in full if 'h10_dVs10' in i]

    for dir in part[::5]:
        print dir
        stp = obspy.read(svdir+dir+'/stp.pk')
        sts = obspy.read(svdir+dir+'/sts.pk')
        for idx,tr in enumerate(sts):
            a = seispy.data.phase_window(tr,['S'],window=(-30,30))
            stp[idx].data *= 1./np.abs(a).max()

        #stp = strip_conversion(stp)
        for tr in stp:
            plt.plot(tr.data,color='k',alpha=0.5)
        plt.show()

def strip_conversion(st):
    for idx,tr in enumerate(st):
        st[idx] = seispy.data.phase_window(st[idx],phase=['S1800P'],window=(-30,10))
    return st

def plot_wave_conversion(st_strip,ax):
    for tr in st_strip:
        ax.plot(tr.data,color='k',alpha=0.5)
    plt.show()

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
