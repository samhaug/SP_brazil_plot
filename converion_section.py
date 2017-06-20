#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : converion_section.py
Purpose : show record section of the manually aligned S1750P
Creation Date : 20-06-2017
Last Modified : Tue 20 Jun 2017 12:24:20 PM EDT
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
    plt.tight_layout()
    plt.savefig('conversion_section.pdf')
    call('evince conversion_section.pdf',shell=True)

def plot_section(st,ax):
    len_time = st[0].stats.npts/st[0].stats.sampling_rate
    t = np.linspace(-len_time/2.,len_time/2.,num=st[0].stats.npts)
    for tr in st:
        ax.plot(t,tr.data+tr.stats.sac['gcarc'],color='k',lw=0.5)
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
    ax.fill_between(t,mean+std,mean-std,color='k',lw=0,alpha=0.3)
    ax.set_xlim(-25,25)

def setup_figure():
    fig = plt.figure(figsize=(4.5,8))
    gs = gridspec.GridSpec(100,100)
    ax_sec = plt.subplot(gs[0:60,0:])
    ax_over = plt.subplot(gs[61:80,0:])
    ax_stack = plt.subplot(gs[81:100,0:])
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
    ax_sec.set_ylabel(r'Epicentral distance ($\circ$)',size=8)
    ax_stack.set_xlabel('Time (s)',size=8)
    plt.figtext(0.05,0.98,'(a)',size=8)
    plt.figtext(0.05,0.38,'(b)',size=8)
    plt.figtext(0.05,0.18,'(c)',size=8)
    return fig,ax_sec,ax_over,ax_stack

main()
