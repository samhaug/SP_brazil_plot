#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : svaxi_waveform_param.py
Purpose : plot change of waveform shape on SVaxi simulation.
Creation Date : 17-04-2017
Last Modified : Mon 17 Apr 2017 02:17:05 PM EDT
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
model = TauPyModel(model='prem50')

def main():
    homedir = '/home/samhaug/work1/SP_brazil_sims/SVaxi/full_shareseismos/shareseismos/'
    prem_st = read_prem(homedir)

    fig,ax = setup_figure()
    thickness_compare(homedir,prem_st,ax[0][0])
    dvs_compare(homedir,prem_st,ax[1][0])
    dvp_compare(homedir,prem_st,ax[1][1])
    drho_compare(homedir,prem_st,ax[0][1])
    angle_compare(homedir,prem_st,ax[0][2])
    dist_compare(homedir,prem_st,ax[1][2])

    ax[0][0].set_title('Thickness',size=7)
    ax[1][0].set_title(r'$\delta V_{S}$',size=7)
    ax[0][1].set_title(r'$\delta \rho$',size=7)
    ax[1][1].set_title(r'$\delta V_{P}$',size=7)
    ax[0][2].set_title(r'$\alpha$',size=7)
    ax[1][2].set_title(r'$Distance$',size=7)

    ax[1][0].set_xlabel('Time (s)',size=7)
    ax[1][1].set_xlabel('Time (s)',size=7)
    ax[1][2].set_xlabel('Time (s)',size=7)

    for ax in ax.reshape(ax.size):
        ax.legend(loc='upper left',prop={'size':6},frameon=False)

    plt.savefig('svaxi_waveform_param.pdf')
    call('evince svaxi_waveform_param.pdf',shell=True)

def thickness_compare(homedir,prem_st,ax):
    st5 = prepare_stream(homedir,'smslab_a0_h5_dVs5',prem_st)
    a = seispy.data.phase_window(st5[10],['S1800P'],window=(-10,10))
    t = np.linspace(-10,10,num=a.stats.npts)
    ax.plot(t,a.data+-1*a.data[0],label='5km',lw=0.5)

    st10 = prepare_stream(homedir,'smslab_a0_h10_dVs5',prem_st)
    a = seispy.data.phase_window(st10[10],['S1800P'],window=(-10,10))
    t = np.linspace(-10,10,num=a.stats.npts)
    ax.plot(t,a.data+-1*a.data[0],label='10km',lw=0.5)

    st20 = prepare_stream(homedir,'smslab_a0_h20_dVs5',prem_st)
    a = seispy.data.phase_window(st20[10],['S1800P'],window=(-10,10))
    t = np.linspace(-10,10,num=a.stats.npts)
    ax.plot(t,a.data+-1*a.data[0],label='20km',lw=0.5)

def dvs_compare(homedir,prem_st,ax):
    st5 = prepare_stream(homedir,'smslab_a0_h10_dVs5',prem_st)
    a = seispy.data.phase_window(st5[10],['S1800P'],window=(-10,10))
    t = np.linspace(-10,10,num=a.stats.npts)
    ax.plot(t,a.data+-1*a.data[0],label='5',lw=0.5)

    st10 = prepare_stream(homedir,'smslab_a0_h10_dVs10',prem_st)
    a = seispy.data.phase_window(st10[10],['S1800P'],window=(-10,10))
    t = np.linspace(-10,10,num=a.stats.npts)
    ax.plot(t,a.data+-1*a.data[0],label='10',lw=0.5)

def dvp_compare(homedir,prem_st,ax):
    st0 = prepare_stream(homedir,'smslab_a0_h10_dVs5',prem_st)
    a = seispy.data.phase_window(st0[10],['S1800P'],window=(-10,10))
    t = np.linspace(-10,10,num=a.stats.npts)
    ax.plot(t,a.data+-1*a.data[0],label='0',lw=0.5)

    st5 = prepare_stream(homedir,'smslab_a0_h10_dVs5_dVp5',prem_st)
    a = seispy.data.phase_window(st5[10],['S1800P'],window=(-10,10))
    t = np.linspace(-10,10,num=a.stats.npts)
    ax.plot(t,a.data+-1*a.data[0],label='5',lw=0.5)

    stn5 = prepare_stream(homedir,'smslab_a0_h10_dVs5_dVp-5',prem_st)
    a = seispy.data.phase_window(stn5[10],['S1800P'],window=(-10,10))
    t = np.linspace(-10,10,num=a.stats.npts)
    ax.plot(t,a.data+-1*a.data[0],label='-5',lw=0.5)

def drho_compare(homedir,prem_st,ax):
    st0 = prepare_stream(homedir,'smslab_a0_h10_dVs5',prem_st)
    a = seispy.data.phase_window(st0[10],['S1800P'],window=(-10,10))
    t = np.linspace(-10,10,num=a.stats.npts)
    ax.plot(t,a.data+-1*a.data[0],label='0',lw=0.5)

    st1 = prepare_stream(homedir,'smslab_a0_h10_dVs5_drho1',prem_st)
    a = seispy.data.phase_window(st1[10],['S1800P'],window=(-10,10))
    t = np.linspace(-10,10,num=a.stats.npts)
    ax.plot(t,a.data+-1*a.data[0],label='1',lw=0.5)

    st3 = prepare_stream(homedir,'smslab_a0_h10_dVs5_drho3',prem_st)
    a = seispy.data.phase_window(st3[10],['S1800P'],window=(-10,10))
    t = np.linspace(-10,10,num=a.stats.npts)
    ax.plot(t,a.data+-1*a.data[0],label='3',lw=0.5)

    st5 = prepare_stream(homedir,'smslab_a0_h10_dVs5_drho5',prem_st)
    a = seispy.data.phase_window(st5[10],['S1800P'],window=(-10,10))
    t = np.linspace(-10,10,num=a.stats.npts)
    ax.plot(t,a.data+-1*a.data[0],label='5',lw=0.5)

def angle_compare(homedir,prem_st,ax):
    st = prepare_stream(homedir,'smslab_a20_h10_dVs10',prem_st)
    a = seispy.data.phase_window(st[10],['S1800P'],window=(-20,20))
    t = np.linspace(-20,20,num=a.stats.npts)
    ax.plot(t,a.data+-1*a.data[0],label='20',lw=0.5)

    st = prepare_stream(homedir,'smslab_a15_h10_dVs10',prem_st)
    a = seispy.data.phase_window(st[10],['S1800P'],window=(-20,20))
    t = np.linspace(-20,20,num=a.stats.npts)
    ax.plot(t,a.data+-1*a.data[0],label='15',lw=0.5)

    st = prepare_stream(homedir,'smslab_a10_h10_dVs10',prem_st)
    a = seispy.data.phase_window(st[10],['S1800P'],window=(-20,20))
    t = np.linspace(-20,20,num=a.stats.npts)
    ax.plot(t,a.data+-1*a.data[0],label='10',lw=0.5)

    st = prepare_stream(homedir,'smslab_a5_h10_dVs10',prem_st)
    a = seispy.data.phase_window(st[10],['S1800P'],window=(-20,20))
    t = np.linspace(-20,20,num=a.stats.npts)
    ax.plot(t,a.data+-1*a.data[0],label='5',lw=0.5)

    st = prepare_stream(homedir,'smslab_a0_h10_dVs10',prem_st)
    a = seispy.data.phase_window(st[10],['S1800P'],window=(-20,20))
    t = np.linspace(-20,20,num=a.stats.npts)
    ax.plot(t,a.data+-1*a.data[0],label='0',lw=0.5)

    st = prepare_stream(homedir,'smslab_a-10_h10_dVs10',prem_st)
    a = seispy.data.phase_window(st[10],['S1800P'],window=(-20,20))
    t = np.linspace(-20,20,num=a.stats.npts)
    ax.plot(t,a.data+-1*a.data[0],label='-10',lw=0.5)

    ax.set_xlim((-20,10))

def dist_compare(homedir,prem_st,ax):
    st = prepare_stream(homedir,'smslab_a0_h10_dVs10',prem_st)
    for tr in st[::4]:
        a = seispy.data.phase_window(tr,['S1800P'],window=(-10,10))
        t = np.linspace(-10,10,num=a.stats.npts)
        label = int(tr.stats.sac['gcarc'])
        ax.plot(t,a.data+-1*a.data[0],label=label,lw=0.5)

def setup_figure():
    fig,ax = plt.subplots(2,3,figsize=(7,5))
    fig.patch.set_visible(False)

    for b in ax.reshape(ax.size):
        b.spines['top'].set_visible(False)
        b.spines['right'].set_visible(False)
        b.spines['left'].set_visible(False)
        b.yaxis.set_ticklabels([])
        b.yaxis.set_ticks([])
        b.tick_params(axis='both', which='major', labelsize=7)
        b.xaxis.set_ticks_position('bottom')
        b.set_ylim((-0.05,0.05))
    '''
    ax[0].spines['bottom'].set_visible(False)
    ax[0].xaxis.set_ticks([])
    ax[1].spines['bottom'].set_visible(False)
    ax[1].xaxis.set_ticks([])
    #ax[2].xaxis.set_ticks([])
    ax[2].set_xlabel('Time (s)',size=10)
    '''
    return fig,ax
    plt.tight_layout()
    plt.show()

def prepare_stream(homedir,dir,prem_st):
    st = obspy.read(homedir+dir+'/stv.pk')
    for idx, tr in enumerate(st):
        st[idx].data = prem_st[idx].data-st[idx].data
    return st

def read_prem(homedir):
    st = obspy.read(homedir+'smoothprem/stv.pk')
    return st


main()
