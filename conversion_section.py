#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : converion_section.py
Purpose : show record section of the manually aligned S1750P
Creation Date : 20-06-2017
Last Modified : Thu 17 Aug 2017 10:11:05 AM EDT
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
import matplotlib.gridspec as gridspec

def main():
    st = obspy.read('/home/samhaug/work1/SP_brazil_setup/Z_align.pk')
    st.normalize()
    fig,ax_sec,ax_over,ax_stack = setup_figure()
    plot_section(st,ax_sec)
    plot_overlap(st,ax_over)
    plot_stack(st,ax_stack)
    #plot_filt(st,ax_filt)
    #plt.tight_layout()
    plt.savefig('conversion_section.pdf')
    call('evince conversion_section.pdf',shell=True)

def plot_section(st,ax):
    len_time = st[0].stats.npts/st[0].stats.sampling_rate
    t = np.linspace(-len_time/2.,len_time/2.,num=st[0].stats.npts)
    for tr in st:
        ax.plot(t,tr.data+tr.stats.sac['gcarc'],color='k',lw=0.5)
    #ax.set_xlim(int(t[0]),int(t[-1]))
    ax.set_xlim(-25,25)
    ax.set_ylim(55,71)

def plot_overlap(st,ax):
    len_time = st[0].stats.npts/st[0].stats.sampling_rate
    t = np.linspace(-len_time/2.,len_time/2.,num=st[0].stats.npts)
    for tr in st:
        ax.plot(t,tr.data,color='k',lw=0.5,alpha=0.3)
    ax.set_xlim(-25,25)

def plot_stack(st,ax):
    len_time = st[0].stats.npts/st[0].stats.sampling_rate
    t = np.linspace(-len_time/2.,len_time/2.,num=st[0].stats.npts)
    a = []
    for tr in st:
        a.append(tr.data)
    std = np.std(a,axis=0)
    mean = np.mean(a,axis=0)
    ax.plot(t,mean,color='k',lw=0.5)
    #ax.text(-10,1.5*mean.max(),'BB',size=6)
    ax.fill_between(t,mean+std,mean-std,color='k',lw=0,alpha=0.3)
    ax.set_xlim(-25,25)

def plot_filt(st,ax):
    len_time = st[0].stats.npts/st[0].stats.sampling_rate
    sr = st[0].stats.sampling_rate
    st.filter('lowpass',freq=1/4.)
    t = np.linspace(-len_time/2.,len_time/2.,num=st[0].stats.npts)
    a = []
    for tr in st:
        a.append(tr.data)
    std = np.roll(np.std(a,axis=0),int(-2*sr))
    mean = np.roll(np.mean(a,axis=0),int(-2*sr))
    ax.plot(t,mean,color='k',lw=0.5)
    ax.text(-10,1.5*mean.max(),'f < 0.25 Hz',size=6)
    ax.fill_between(t,mean+std,mean-std,color='k',lw=0,alpha=0.3)
    ax.set_xlim(-25,25)

def setup_figure():
    fig = plt.figure(figsize=(5.5,8))
    gs = gridspec.GridSpec(70,100)
    ax_sec = plt.subplot(gs[0:48,0:])
    ax_over = plt.subplot(gs[50:59,0:])
    ax_stack = plt.subplot(gs[60:69,0:])
    #ax_filt = plt.subplot(gs[70:80,0:])
    for b in [ax_over,ax_stack]:
        b.spines['top'].set_visible(False)
        b.spines['right'].set_visible(False)
        b.spines['left'].set_visible(False)
        b.yaxis.set_ticklabels([])
        b.yaxis.set_ticks([])
        b.tick_params(axis='both', which='major', labelsize=6)
        b.xaxis.set_ticks_position('bottom')
    for b in [ax_sec,ax_over]:
        b.spines['bottom'].set_visible(False)
        b.xaxis.set_ticklabels([])
        b.xaxis.set_ticks([])
    ax_sec.spines['top'].set_visible(False)
    ax_sec.spines['right'].set_visible(False)
    ax_sec.yaxis.set_ticks_position('left')
    ax_sec.tick_params(axis='both', which='major', labelsize=6)
    ax_sec.set_ylabel(r'Epicentral distance (degrees)',size=8)
    ax_stack.set_xlabel('Time (s)',size=8)
    #ax_stack.xaxis.set_ticks(np.arange(-20,20,2))
    major_ticks = np.arange(-25,30,5)
    minor_ticks = np.arange(-25,26,1)
    ax_stack.set_xticks(major_ticks)
    ax_stack.set_xticks(minor_ticks,minor=True)
    ax_stack.spines['bottom'].set_bounds(-25,25)
    plt.figtext(0.07,0.93,'(a)',size=7)
    plt.figtext(0.07,0.30,'(b)',size=7)
    plt.figtext(0.07,0.20,'(c)',size=7)
    #plt.figtext(0.07,0.18,'(d)',size=8)
    return fig,ax_sec,ax_over,ax_stack

main()
