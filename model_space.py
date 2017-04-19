#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : model_space.py
Purpose : Plot model space search
Creation Date : 18-04-2017
Last Modified : Wed 19 Apr 2017 10:22:14 AM EDT
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
    #homedir = '/home/samhaug/work1/SP_brazil_sims/SVaxi/0327_shareseismos/shareseismos/'
    homedir = '/home/samhaug/work1/SP_brazil_sims/SVaxi/full_shareseismos/shareseismos/'
    d_rat,s_rat = read_ratio('/home/samhaug/work1/SP_brazil_setup/ratio.dat')
    read_bootstrap('../SP_brazil_setup',d_rat)
    fig,ax_list = setup_figure(d_rat)
    stprem = obspy.read(homedir+'smoothprem/stv.pk')
    for dir in listdir(homedir):
        if dir == 'smoothprem':
            continue
        ang = float(dir.split('_')[1][1::])
        if ang > 15 or ang < -15:
            continue
        print dir
        st = read_synth(homedir+dir)

        for idx,tr in enumerate(st):
            st[idx].data = st[idx].data-stprem[idx].data

        if dir.split('_')[2]+'_'+dir.split('_')[3] == 'h5_dVs5':
            ax_list[0][0].set_title('h5_dVs5',size=10)
            ax = plot_amp(st,ax_list[0][0],ang)
        if dir.split('_')[2]+'_'+dir.split('_')[3] == 'h10_dVs10':
            ax_list[1][1].set_title('h10_dVs10',size=10)
            ax = plot_amp(st,ax_list[1][1],ang)
        if dir.split('_')[2]+'_'+dir.split('_')[3] == 'h20_dVs5':
            ax_list[0][1].set_title('h5_dVs10',size=10)
            ax = plot_amp(st,ax_list[0][1],ang)
        if dir.split('_')[2]+'_'+dir.split('_')[3] == 'h10_dVs5':
            ax_list[1][0].set_title('h10_dVs5',size=10)
            ax = plot_amp(st,ax_list[1][0],ang)
    plt.savefig('normalize_svaxi.pdf')
    call('evince normalize_svaxi.pdf',shell=True)

def read_synth(dir):
    st = obspy.read(dir+'/stv.pk')
    #st = seispy.filter.range_fiter(st,(56,72))
    #for tr in st:
    #    tr.stats.sac['o'] += -202
    #st = seispy.data.normalize_on_phase(st,phase=['S'],min=True)
    #for tr in st:
    #    tr.data *= 1./np.tan(np.radians(17))
    return st

def read_bootstrap(dir,d_rat):
    f = h5py.File(dir+'/align_bootstrap_old.h5','r')
    align_array = f['boot'][...]
    plt.plot(np.mean(align_array,axis=0)*np.mean(d_rat))
    plt.show()

def setup_figure(d_rat):
    fig,ax_list = plt.subplots(2,2,figsize=(5,5))

    amp = 0.059*d_rat
    std = 0.01
    for ax in ax_list.reshape(ax_list.size):
        #ax.yaxis.set_ticklabels([])
        #ax.yaxis.set_ticks([])
        ax.tick_params(axis='both', which='major', labelsize=6)
        ax.xaxis.set_ticks_position('bottom')
        ax.set_xlim((-20,20))
        ax.set_ylim((0,0.1))
        ax.axhline(amp,color='k',lw=0.5)
        ax.fill_between(np.linspace(-30,30),amp-std,amp+std,color='b',alpha=0.3)
        #ax.set_ylim((-1,14))

    #ax_list[0][0].set_xticklabels([0])
    #ax_list[0][1].set_xticklabels([0])
    ax_list[1][0].set_xlabel('Slab angle',size=8)
    ax_list[1][1].set_xlabel('Slab angle',size=8)
    ax_list[0][0].set_ylabel('Amplitude',size=8)
    ax_list[1][0].set_ylabel('Amplitude',size=8)


    return fig,ax_list

def plot_amp(st,ax,ang):
    a = []
    for tr in st:
        d = seispy.data.phase_window(tr,['S1800P'],window=(-12,8)).data
        a.append(d.max()-d.min())
    #for ii in b:
    #    ax1.plot(ii,color='k',alpha=0.5)
    #plt.show()
    amax = np.max(a)
    amin = np.min(a)
    mid = (amax+amin)/2.
    ax.errorbar(ang, mid,yerr=amax-mid,color='k')

def read_ratio(ratio_file):
    a = np.genfromtxt(ratio_file)
    data_rat = np.mean(a[:,1])
    synth_rat = np.mean(a[:,2])
    return data_rat,synth_rat

main()


