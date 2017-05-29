#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : frequency_compare.py
Purpose : compare conversion amplitude as fcn of frequency
Creation Date : 29-05-2017
Last Modified : Mon 29 May 2017 12:16:32 PM EDT
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
from subprocess import call
rcParams['axes.color_cycle'] = ['k','#1f77b4','#ff7f0e','#2ca02c',
                                '#d62728','#9467bd','#8c564b',
                                '#e377c2','#7f7f7f','#bcbd22','#17becf']

def main():
    fig,ax = setup_figure()
    stz,strr = read_stream('/home/samhaug/work1/SP_brazil_setup/')
    freq_list = [100,1,1/2.,1/3.,1/4.]

    for freq in freq_list:
        stz_c,str_c = filter_stream(stz,strr,freq)
        mean,std = normalize_streams(stz_c,str_c)
        t = np.linspace(-30,30,num=len(mean))
        if freq == 100:
            ax[0].plot(t,mean,lw=1,label='BB')
            ax[1].scatter(0,mean[1200:1800].max()-mean[1200:1800].min(),
                          color='k',marker='D')
            i = np.argmax(mean[1200:1800])
            err = std[1200+i]
            amp = mean[1200:1800].max()-mean[1200:1800].min()
            ax[1].errorbar(0,amp,yerr=err,color='k')
        else:
            ax[0].plot(t,mean,lw=1,label=str(int(1./freq))+' s')
            ax[1].scatter(int(1./freq),
                         mean[1200:1800].max()-mean[1200:1800].min(),
                         color='k',marker='D')
            i = np.argmax(mean[1200:1800])
            err = std[1200+i]
            amp = mean[1200:1800].max()-mean[1200:1800].min()
            ax[1].errorbar(int(1./freq),amp,yerr=err,color='k')
            print freq,amp,err
    ax[0].legend(loc='upper right',prop={'size':6},frameon=False)
    plt.savefig('frequency_compare.pdf')
    call('evince frequency_compare.pdf',shell=True)

def normalize_streams(stz,str):
    a = []
    for idx,tr in enumerate(stz):
        stz[idx].data *= 1./str[idx].data.max()
        a.append(stz[idx].data)
    mean = np.mean(a,axis=0)
    std = np.std(a,axis=0)
    return mean,std

def filter_stream(stz,str,freq):
    stz_c = stz.copy()
    str_c = str.copy()
    stz_c.filter('lowpass',freq=freq)
    str_c.filter('lowpass',freq=freq)
    return stz_c,str_c

def setup_figure():
    fig,ax_list = plt.subplots(2,1,figsize=(4,8))
    fig.patch.set_visible(False)
    for b in ax_list:
        b.set_xlim(-10,10)
        b.spines['top'].set_visible(False)
        b.spines['right'].set_visible(False)
        b.yaxis.set_ticks_position('left')
        b.tick_params(axis='both', which='major', labelsize=6)
        b.xaxis.set_ticks_position('bottom')
    ax_list[0].set_ylim(-0.1,0.1)
    ax_list[1].set_ylim(0.0,0.20)
    ax_list[0].set_xlim(-10,10)
    ax_list[1].set_xlim(-1,5)
    ax_list[0].set_xlabel('Time (s)',size=8)
    ax_list[0].set_ylabel('Amp.',size=8)
    ax_list[1].set_xlabel('Corner Period (s)',size=8)
    ax_list[1].set_ylabel('P2P Amp.',size=8)
    #labels = [item.get_text() for item in ax_list[1].get_xticklabels()]
    #print labels
    #labels[1] = 'BB'
    labels = ['','BB','1','2','3','4','']
    ax_list[1].set_xticklabels(labels)
    return fig,ax_list

def read_stream(streamdir):
    stz = obspy.read(streamdir+'Z_align.pk')
    str = obspy.read(streamdir+'R.pk')
    stz.sort(['station'])
    str.sort(['station'])
    return stz,str

main()
