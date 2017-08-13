#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : model_space.py
Purpose : Plot model space search
Creation Date : 18-04-2017
Last Modified : Fri 11 Aug 2017 02:59:53 PM EDT
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
    homedir = '/home/samhaug/work1/SP_brazil_sims/SVaxi/full_shareseismos/shareseismos/'
    prem_st = read_prem(homedir)
    d_rat,s_rat = read_ratio('/home/samhaug/work1/SP_brazil_setup/ratio.dat')
    #read_bootstrap('../SP_brazil_setup',d_rat)
    dirty_amp,dirty_std = read_dirty()
    fig,ax_list = setup_figure(d_rat,dirty_amp,dirty_std)

    vs_plot(ax_list[0][0],homedir,prem_st)
    vp_plot(ax_list[0][1],homedir,prem_st)
    rho_plot(ax_list[0][2],homedir,prem_st)
    thick_plot(ax_list[1][0],homedir,prem_st)
    angle_plot(ax_list[1][2],homedir,prem_st)
    distance_plot(ax_list[1][1],homedir,prem_st)

    plt.savefig('model_space.pdf')
    call('evince model_space.pdf',shell=True)

def vs_plot(ax,homedir,prem_st):
    print('Vs plot')
    st = prepare_stream(homedir,'smslab_a0_h10_dVs5',prem_st)
    plot_amp(st,ax,-5)
    st = prepare_stream(homedir,'smslab_a0_h10_dVs10',prem_st)
    plot_amp(st,ax,-10)
    props = dict(boxstyle='square',facecolor='white',alpha=1.0,lw=0.5)
    textstr=r'$h=10km$, $\alpha=0^{\circ}$, $\delta V_{P}=0\%$, $\delta \rho=0\%$'
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=5,
            verticalalignment='top', bbox=props)

def vp_plot(ax,homedir,prem_st):
    print('Vp plot')
    st = prepare_stream(homedir,'smslab_a0_h10_dVs5',prem_st)
    plot_amp(st,ax,0,multiply=2.)
    st = prepare_stream(homedir,'smslab_a0_h10_dVs5_dVp-5',prem_st)
    plot_amp(st,ax,-5,multiply=2.)
    st = prepare_stream(homedir,'smslab_a0_h10_dVs5_dVp5',prem_st)
    plot_amp(st,ax,5,multiply=2.)
    props = dict(boxstyle='square',facecolor='white',alpha=1.0,lw=0.5)
    textstr=r'$h=10km$, $\alpha=0^{\circ}$, $\delta V_{S}=-10\%$, $\delta \rho=0\%$'
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=5,
            verticalalignment='top', bbox=props)

def rho_plot(ax,homedir,prem_st):
    print('Rho plot')
    st = prepare_stream(homedir,'smslab_a0_h10_dVs5',prem_st)
    plot_amp(st,ax,0,multiply=2.)
    st = prepare_stream(homedir,'smslab_a0_h10_dVs5_drho1',prem_st)
    plot_amp(st,ax,1,multiply=2.)
    st = prepare_stream(homedir,'smslab_a0_h10_dVs5_drho3',prem_st)
    plot_amp(st,ax,3,multiply=2.)
    st = prepare_stream(homedir,'smslab_a0_h10_dVs5_drho5',prem_st)
    plot_amp(st,ax,5,multiply=2.)
    props = dict(boxstyle='square',facecolor='white',alpha=1.0,lw=0.5)
    textstr=r'$h=10km$, $\alpha=0^{\circ}$, $\delta V_{S}=--10\%$, $\delta V_{P}=0\%$'
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=5,
            verticalalignment='top', bbox=props)

def thick_plot(ax,homedir,prem_st):
    print('Thickness plot')
    st = prepare_stream(homedir,'smslab_a0_h2_dVs5',prem_st)
    plot_amp(st,ax,2,multiply=2.)
    st = prepare_stream(homedir,'smslab_a0_h5_dVs5',prem_st)
    plot_amp(st,ax,5,multiply=2.)
    st = prepare_stream(homedir,'smslab_a0_h10_dVs5',prem_st)
    plot_amp(st,ax,10,multiply=2.)
    st = prepare_stream(homedir,'smslab_a0_h20_dVs5',prem_st)
    plot_amp(st,ax,20,multiply=2.)
    props = dict(boxstyle='square',facecolor='white',alpha=1.0,lw=0.5)
    textstr=r'$\alpha=0^{\circ}$, $\delta V_{S}=-10\%$, $\delta V_{P}=0\%$, $\delta \rho = 0\%$'
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=5,
            verticalalignment='top', bbox=props)

def angle_plot(ax,homedir,prem_st):
    print('Angle plot')
    st = prepare_stream(homedir,'smslab_a-20_h10_dVs5',prem_st)
    plot_amp(st,ax,-20,multiply=2.)
    st = prepare_stream(homedir,'smslab_a-10_h10_dVs5',prem_st)
    plot_amp(st,ax,-10,multiply=2.)
    st = prepare_stream(homedir,'smslab_a0_h10_dVs5',prem_st)
    plot_amp(st,ax,0,multiply=2.)
    st = prepare_stream(homedir,'smslab_a10_h10_dVs5',prem_st)
    plot_amp(st,ax,10,multiply=2.)
    st = prepare_stream(homedir,'smslab_a15_h10_dVs5',prem_st)
    plot_amp(st,ax,15,multiply=2.)
    st = prepare_stream(homedir,'smslab_a20_h10_dVs5',prem_st)
    plot_amp(st,ax,20,multiply=2.)
    props = dict(boxstyle='square',facecolor='white',alpha=1.0,lw=0.5)
    textstr=r'$h=10km$, $\delta V_{S}=-10\%$, $\delta V_{P}=0\%$, $\delta \rho = 0\%$'
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=5,
            verticalalignment='top', bbox=props)

def distance_plot(ax,homedir,prem_st):
    print('Distance plot')
    st = prepare_stream(homedir,'smslab_a0_h10_dVs10',prem_st)
    amp = []
    dist = []
    for tr in st[::1]:
        d = seispy.data.phase_window(tr,['S1800P'],window=(-12,8)).data
        amp.append(d.max()-d.min())
        dist.append(tr.stats.sac['gcarc'])
    ax.scatter(dist,amp,marker='D',color='k',s=3)
    ax.set_xlim((np.min(dist)-5.,np.max(dist)+5))
    props = dict(boxstyle='square',facecolor='white',alpha=1.0,lw=0.5)
    textstr=r'$h=10km$, $\alpha=0^{\circ}$, $\delta V_{S}=-10\%$, $\delta V_{P}=0\%$, $\delta \rho = 0\%$'
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=5,
            verticalalignment='top', bbox=props)

def setup_figure(d_rat,dirty_amp,dirty_std):
    fig,ax_list = plt.subplots(2,3,figsize=(7,5))

    #amp = 0.059*d_rat
    #std = 0.01
    #amp = 0.10*d_rat
    #std = 0.0425
    amp = 0.0417370943722*d_rat
    std = 0.032767920367*d_rat
    for ax in ax_list.reshape(ax_list.size):
        ax.tick_params(axis='both', which='major', labelsize=6)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_ylim((0,0.1))
        ax.axhline(amp,color='k',lw=0.5)
        #ax.axhline(dirty_amp,color='r',lw=0.5)
        ax.fill_between(np.linspace(-30,80),amp-std,amp+std,color='gray',alpha=0.2,lw=0)
        #ax.fill_between(np.linspace(-300,300),dirty_amp-dirty_std,dirty_amp+dirty_std,color='red',alpha=0.4,lw=0)
        #ax.set_ylim((-1,14))

    ax_list[0][1].yaxis.set_ticklabels([])
    ax_list[0][2].yaxis.set_ticklabels([])
    ax_list[1][1].yaxis.set_ticklabels([])
    ax_list[1][2].yaxis.set_ticklabels([])
    ax_list[1][2].yaxis.set_ticks([])
    ax_list[0][2].yaxis.set_ticks([])
    ax_list[0][1].yaxis.set_ticks([])
    ax_list[1][1].yaxis.set_ticks([])

    ax_list[0][0].set_xlabel(r'$\delta V_{S} (\%)$',size=7)
    ax_list[0][0].set_xlim((0,-15))

    ax_list[0][1].set_xlabel(r'$\delta V_{P} (\%)$',size=7)
    ax_list[0][1].set_xlim((-10,10))

    ax_list[0][2].set_xlabel(r'$\delta \rho (\%)$',size=7)
    ax_list[0][2].set_xlim((-2,8))

    ax_list[1][0].set_xlabel('Thickness (km)',size=7)
    ax_list[1][0].set_xlim((0,25))

    ax_list[1][1].set_xlabel(r'Distance $(^{\circ})$',size=7)

    ax_list[1][2].set_xlabel(r'Angle $(^{\circ})$',size=7)
    ax_list[1][2].set_xlim((-25,25))


    plt.tight_layout()

    return fig,ax_list

def plot_amp(st,ax,ang,**kwargs):
    multiply = kwargs.get('multiply',1.)
    a = []
    for tr in st:
        d = seispy.data.phase_window(tr,['S1800P'],window=(-12,8)).data
        d *= multiply
        a.append(d.max()-d.min())
    #for ii in b:
    #    ax1.plot(ii,color='k',alpha=0.5)
    #plt.show()
    amax = np.max(a)
    amin = np.min(a)
    mid = ((amax+amin)/2.)
    ax.errorbar(ang, mid,yerr=(amax-mid),color='k')

def read_ratio(ratio_file):
    a = np.genfromtxt(ratio_file)
    data_rat = np.mean(a[:,1])
    synth_rat = np.mean(a[:,2])
    return data_rat,synth_rat

def prepare_stream(homedir,dir,prem_st):
    st = obspy.read(homedir+dir+'/stv.pk')
    for idx, tr in enumerate(st):
        st[idx].data = prem_st[idx].data-st[idx].data
    return st

def read_prem(homedir):
    st = obspy.read(homedir+'smoothprem/stv.pk')
    return st

def read_dirty():
     f = h5py.File('/home/samhaug/work1/SP_brazil_data/2007-07-21-mw60-western-brazil-4/dirty_stacks.h5','r')
     b = f['dirty'][...]
     mean,std = np.mean(b,axis=0),np.std(b,axis=0)
     amp = std[900:1100].max()-std[900:1100].min()
     amp_std = std[900:1100].max()-std[900:1100].min()
     return amp,amp_std

def read_bootstrap(dir,d_rat):
    f = h5py.File(dir+'/align_bootstrap_old.h5','r')
    align_array = f['boot'][...]
    '''
    plt.plot(np.mean(align_array,axis=0)*np.mean(d_rat))
    plt.show()
    '''


main()


