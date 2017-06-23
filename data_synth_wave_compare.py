#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : data_synth_wave_compare.py
Purpose : figure for directly comparing synth/data waveforms.
Creation Date : 22-06-2017
Last Modified : Fri 23 Jun 2017 01:20:49 PM EDT
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
from matplotlib import rcParams

rcParams['axes.color_cycle'] = ['#1f77b4','#ff7f0e','#2ca02c',
                                '#d62728','#9467bd','#8c564b',
                                '#e377c2','#7f7f7f','#bcbd22','#17becf']

def main():
    dirname = '/home/samhaug/work1/SP_brazil_sims/SVaxi/full_shareseismos/shareseismos/'

    t_dat,mean_dat,std_dat = setup_data()
    mean_dat *= 2.0/mean_dat.max()
    mean_dat = np.roll(mean_dat,-100)
    t_dat,mean_dat = remove_trend(t_dat,mean_dat)
    fig,ax = setup_figure()
    plt.tight_layout()
    ax[0][0].plot(t_dat,mean_dat,color='k',alpha=0.5,lw=1.5)
    ax[1][0].plot(t_dat,mean_dat,color='k',alpha=0.5,lw=1.5)
    ax[1][1].plot(t_dat,mean_dat,color='k',alpha=0.5,lw=1.5)

    t,mean = setup_sv(dirname+'smslab_a0_h10_dVs10/stv_strip.pk')
    ax[0][0].plot(t,-1.*mean/np.abs(mean).max(),lw=0.5,label='dVs=-10%')
    t,mean = setup_sv(dirname+'smslab_a0_h10_dVs5/stv_strip.pk')
    wave = mean/np.abs(mean).max()/2.
    wave = wave-wave.mean()
    ax[0][0].plot(t,wave*-1,lw=0.5,label='dVs=-5%')
    ax[0][0].plot(t,0.1+wave,lw=0.5,label='dVs=5%')
    ax[0][0].legend(loc='upper left',prop={'size':6},frameon=False)
    ax[0][0].set_xlabel('Time (s)',size=6)

    t,mean = setup_sv(dirname+'smslab_a0_h2_dVs5/stv_strip.pk')
    ax[1][0].plot(t,-1.*mean/mean.max(),lw=0.5,label='2km')
    t,mean = setup_sv(dirname+'smslab_a0_h5_dVs5/stv_strip.pk')
    ax[1][0].plot(t,-1.*mean/mean.max(),lw=0.5,label='5km')
    t,mean = setup_sv(dirname+'smslab_a0_h10_dVs5/stv_strip.pk')
    ax[1][0].plot(t,-0.05+-1*mean/mean.max(),lw=0.5,label='10km')
    #t,mean = setup_sv(dirname+'smslab_a0_h20_dVs5/stv_strip.pk')
    #ax[1][0].plot(t,-1*mean/mean.max(),lw=0.5,label='20km')
    ax[1][0].legend(loc='upper left',prop={'size':6},frameon=False)
    ax[1][0].set_xlabel('Time (s)',size=6)

    t,mean = setup_sv(dirname+'smslab_a0_h10_dVs5/stv_strip.pk')
    ax[1][1].plot(t,-1.*mean/mean.max(),lw=0.5,label='0')
    t,mean = setup_sv(dirname+'smslab_a-10_h10_dVs5/stv_strip.pk')
    mean = np.roll(mean,-70)
    ax[1][1].plot(t,-1.*mean/mean.max(),lw=0.5,label='-10')
    t,mean = setup_sv(dirname+'smslab_a10_h10_dVs5/stv_strip.pk')
    mean = np.roll(mean,90)
    ax[1][1].plot(t,-1*mean/mean.max(),lw=0.5,label='10')
    #t,mean = setup_sv(dirname+'smslab_a-20_h10_dVs5/stv_strip.pk')
    #ax[1][1].plot(t,-1*mean/mean.max(),lw=0.5,label='-20')
    #t,mean = setup_sv(dirname+'smslab_a20_h10_dVs5/stv_strip.pk')
    #ax[1][1].plot(t,-1*mean/mean.max(),lw=0.5,label='20')
    ax[1][1].legend(loc='upper left',prop={'size':6},frameon=False)
    ax[1][1].set_xlabel('Time (s)',size=6)

    plt.savefig('data_synth_wave_compare.pdf')
    call('evince data_synth_wave_compare.pdf',shell=True)

def setup_data():
    std = obspy.read('/home/samhaug/work1/SP_brazil_setup/Z_align.pk')
    std.filter('lowpass',freq=1./4)
    a = []
    full_time = std[0].stats.npts/std[0].stats.sampling_rate
    t = np.linspace(-full_time/2.,full_time/2.,num=std[0].stats.npts)
    for tr in std:
        a.append(tr.data)
    return t,np.mean(a,axis=0),np.std(a,axis=0)

def remove_trend(t,mean):
    neg = np.argmin(np.abs(t+10))
    pos = np.argmin(np.abs(t-10))
    t = t[neg:pos]
    mean = mean[neg:pos]
    p = np.poly1d(np.polyfit(t,mean,3))
    #plt.plot(t,mean)
    #plt.plot(t,p(t))
    #plt.plot(t,mean-p(t))
    #plt.show()
    return t,mean-p(t)

def setup_sv(stdir):
    st = obspy.read(stdir)
    len_time = st[0].stats.npts/st[0].stats.sampling_rate
    sr = st[0].stats.sampling_rate
    t = np.linspace(-len_time/2.,len_time/2.,num=st[0].stats.npts)
    a = []
    for tr in st:
        a.append(tr.data[0:1796])
    std = np.roll(np.std(a,axis=0),int(-2*sr))
    mean = np.roll(np.mean(a,axis=0),int(-2*sr))
    return t[0:1796],np.roll(mean,-40)

def setup_figure():
    fig,ax = plt.subplots(2,2,figsize=(5.5,5.5))

    for b in ax.reshape(ax.size):
        b.spines['top'].set_visible(False)
        b.spines['right'].set_visible(False)
        b.spines['left'].set_visible(False)
        b.yaxis.set_ticklabels([])
        b.yaxis.set_ticks([])
        b.tick_params(axis='both', which='major', labelsize=6)
        b.xaxis.set_ticks_position('bottom')
        b.set_xlim(-10,10)
        b.set_ylim(-1.4,1.4)
        b.spines['bottom'].set_bounds(-5,5)
        b.get_xticklabels()[0].set_visible(False)
        b.get_xticklabels()[-1].set_visible(False)
        b.get_xticklines()[0].set_visible(False)
        b.get_xticklines()[-2].set_visible(False)

    ax[0][1].spines['bottom'].set_visible(False)
    ax[0][1].xaxis.set_ticklabels([])
    ax[0][1].xaxis.set_ticks([])
    return fig,ax

main()



