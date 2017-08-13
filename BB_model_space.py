#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : model_space.py
Purpose : Plot model space search with high frequency runs
Creation Date : 18-04-2017
Last Modified : Fri 11 Aug 2017 04:08:56 PM EDT
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
    homedir = '/home/samhaug/work1/SP_brazil_sims/SVaxi/high_freq/'
    #prem_st = read_prem(homedir)
    #d_rat,s_rat = read_ratio('/home/samhaug/work1/SP_brazil_setup/ratio.dat')
    #read_bootstrap('../SP_brazil_setup',d_rat)
    #dirty_amp,dirty_std = read_dirty()
    st_data = obspy.read('/home/samhaug/work1/SP_brazil_setup/stz_amplitude_processed.pk')
    fig,ax_list = setup_figure(st_data)

    vs_plot(ax_list[0][0],homedir)
    #vp_plot(ax_list[0][1],homedir,prem_st)
    #rho_plot(ax_list[0][2],homedir,prem_st)
    thick_plot(ax_list[1][0],homedir)
    angle_plot(ax_list[0][1],homedir)
    #distance_plot(ax_list[1][1],homedir)
    lowf_vp_vs_plot(ax_list[1][1],'/home/samhaug/work1/SP_brazil_sims/SVaxi/full_shareseismos/shareseismos/')

    plt.savefig('BB_model_space.pdf')
    call('evince BB_model_space.pdf',shell=True)

def vs_plot(ax,homedir):
    print('Vs plot')
    y_5 = [0]
    x_5 = [0]
    y_10 = [0]
    x_10 = [0]
    amps = prepare_stream(homedir,'smslab_a0_h5_dVs5')
    plot_amp(amps,ax,-5,multiply=1.)
    x_5.append(-5)
    y_5.append(np.mean(amps))
    plot_amp(amps,ax,-5,multiply=1.2916,color='green')
    x_10.append(-5)
    y_10.append(np.mean(amps)*1.2916)

    amps = prepare_stream(homedir,'smslab_a0_h5_dVs10')
    plot_amp(amps,ax,-10,multiply=1.)
    x_5.append(-10)
    y_5.append(np.mean(amps))
    plot_amp(amps,ax,-10,multiply=1.2916,color='green')
    x_10.append(-10)
    y_10.append(np.mean(amps)*1.2916)

    amps = prepare_stream(homedir,'smslab_a0_h5_dVs15')
    plot_amp(amps,ax,-15,multiply=1.)
    x_5.append(-15)
    y_5.append(np.mean(amps))
    plot_amp(amps,ax,-15,multiply=1.2916,color='green')
    x_10.append(-15)
    y_10.append(np.mean(amps)*1.2916)

    props = dict(boxstyle='square',facecolor='white',alpha=1.0,lw=0.5)
    textstr=r'$\alpha=0^{\circ}$, $\delta V_{P}=0\%$, $\delta \rho=0\%$'
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=5,
            verticalalignment='top', bbox=props)
    ax.plot(np.unique(x_5), np.poly1d(np.polyfit(x_5, y_5, 1))(np.unique(x_5)),
            color='k',alpha=0.5,ls='--')
    ax.plot(np.unique(x_10), np.poly1d(np.polyfit(x_10, y_10, 1))(np.unique(x_10)),color='g',alpha=0.5,ls='--')
    ax.text(-12,0.064,'h=5km',color='k',rotation=45,size=6)
    ax.text(-12,0.084,'h=10km',color='g',rotation=50,size=6)
    #ax.axvspan(0.25,-2,alpha=0.2,color='green',lw=0)

def vp_plot(ax,homedir,prem_st):
    print('Vp plot')
    st = prepare_stream(homedir,'smslab_a0_h10_dVs5',prem_st)
    plot_amp(st,ax,0,multiply=1.)
    st = prepare_stream(homedir,'smslab_a0_h10_dVs5_dVp-5',prem_st)
    plot_amp(st,ax,-5,multiply=1.)
    st = prepare_stream(homedir,'smslab_a0_h10_dVs5_dVp5',prem_st)
    plot_amp(st,ax,5,multiply=1.)
    props = dict(boxstyle='square',facecolor='white',alpha=1.0,lw=0.5)
    textstr=r'$h=10km$, $\alpha=0^{\circ}$, $\delta V_{S}=-10\%$, $\delta \rho=0\%$'
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=5,
            verticalalignment='top', bbox=props)

def rho_plot(ax,homedir,prem_st):
    print('Rho plot')
    st = prepare_stream(homedir,'smslab_a0_h10_dVs5',prem_st)
    plot_amp(st,ax,0,multiply=1.)
    st = prepare_stream(homedir,'smslab_a0_h10_dVs5_drho1',prem_st)
    plot_amp(st,ax,1,multiply=1.)
    st = prepare_stream(homedir,'smslab_a0_h10_dVs5_drho3',prem_st)
    plot_amp(st,ax,3,multiply=1.)
    st = prepare_stream(homedir,'smslab_a0_h10_dVs5_drho5',prem_st)
    plot_amp(st,ax,5,multiply=1.)
    props = dict(boxstyle='square',facecolor='white',alpha=1.0,lw=0.5)
    textstr=r'$h=10km$, $\alpha=0^{\circ}$, $\delta V_{S}=--10\%$, $\delta V_{P}=0\%$'
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=5,
            verticalalignment='top', bbox=props)

def thick_plot(ax,homedir):
    print('Thickness plot')
    amps = prepare_stream(homedir,'smslab_a0_h2_dVs10')
    plot_amp(amps,ax,2)
    amps = prepare_stream(homedir,'smslab_a0_h5_dVs10')
    plot_amp(amps,ax,5)
    amps = prepare_stream(homedir,'smslab_a0_h10_dVs10')
    plot_amp(amps,ax,10)
    amps = prepare_stream(homedir,'smslab_a0_h20_dVs10')
    plot_amp(amps,ax,20)
    props = dict(boxstyle='square',facecolor='white',alpha=1.0,lw=0.5)
    textstr=r'$\alpha=0^{\circ}$, $\delta V_{S}=-10\%$, $\delta V_{P}=0\%$, $\delta \rho = 0\%$'
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=5,
            verticalalignment='top', bbox=props)

def angle_plot(ax,homedir):
    print('Angle plot')
    amps = prepare_stream(homedir,'smslab_a-10_h10_dVs10')
    plot_amp(amps,ax,-10,multiply=1)
    amps = prepare_stream(homedir,'smslab_a10_h10_dVs10')
    plot_amp(amps,ax,10,multiply=1)
    amps = prepare_stream(homedir,'smslab_a0_h10_dVs10')
    plot_amp(amps,ax,0,multiply=1)
    amps = prepare_stream(homedir,'smslab_a20_h10_dVs10')
    plot_amp(amps,ax,20,multiply=1)

    props = dict(boxstyle='square',facecolor='white',alpha=1.0,lw=0.5)
    textstr=r'$h=10km$, $\delta V_{S}=-10\%$, $\delta V_{P}=0\%$, $\delta \rho = 0\%$'
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=5,
            verticalalignment='top', bbox=props)

def distance_plot(ax,homedir):
    print('Distance plot')
    st = obspy.read(homedir+'smslab_a0_h10_dVs10/stv_strip.pk')
    amp = []
    dist = []
    for tr in st[7::]:
        d = tr.data[2000:4500]
        amp.append(d.max()-d.min())
        dist.append(tr.stats.sac['gcarc'])
    ax.scatter(dist,amp,marker='D',color='k',s=3)
    ax.set_xlim((np.min(dist)-5.,np.max(dist)+5))
    props = dict(boxstyle='square',facecolor='white',alpha=1.0,lw=0.5)
    textstr=r'$h=10km$, $\alpha=0^{\circ}$, $\delta V_{S}=-10\%$, $\delta V_{P}=0\%$, $\delta \rho = 0\%$'
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=5,
            verticalalignment='top', bbox=props)

def lowf_vp_vs_plot(ax,homedir):
    print('lowf_vp_vs')
    amps = lowf_prepare_stream(homedir,'smslab_a0_h10_dVs5_drho1')
    plot_amp(amps,ax,1,multiply=1)
    amps = lowf_prepare_stream(homedir,'smslab_a0_h10_dVs5_drho3')
    plot_amp(amps,ax,3,multiply=1)
    amps = lowf_prepare_stream(homedir,'smslab_a0_h10_dVs5_drho5')
    plot_amp(amps,ax,5,multiply=1)
    amps = lowf_prepare_stream(homedir,'smslab_a0_h10_dVs5_dVp5')
    plot_amp(amps,ax,5,multiply=1,color='g')
    amps = lowf_prepare_stream(homedir,'smslab_a0_h10_dVs5_dVp-5')
    plot_amp(amps,ax,-5,multiply=1,color='g')
    amps = lowf_prepare_stream(homedir,'smslab_a0_h10_dVs5')
    plot_amp(amps,ax,0,multiply=1,color='g')
    ax.set_xlim(-10,10)

    props = dict(boxstyle='square',facecolor='white',alpha=1.0,lw=0.5)
    textstr=r'$h=10km$, $\delta V_{S}=-5\%$, $\alpha=0^{\circ}$'
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=5,
            verticalalignment='top', bbox=props)

def setup_figure(stream):
    #fig,ax_list = plt.subplots(2,3,figsize=(7,5))
    fig,ax_list = plt.subplots(2,2,figsize=(5,5))
    a = []
    for tr in stream:
        a.append(tr.data)
    m = np.mean(a,axis=0)[1200:1800]
    s = np.std(a,axis=0)[1200:1800]
    amp = m.max()-m.min()
    stdmax = s[np.argmax(m)]
    stdmin = s[np.argmin(m)]

    for ax in ax_list.reshape(ax_list.size):
        ax.tick_params(axis='both', which='major', labelsize=6)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_ylim((0,0.1))
        ax.axhline(amp,color='k',lw=0.5)
        #ax.axhline(dirty_amp,color='r',lw=0.5)
        ax.fill_between(np.linspace(-40,80),amp-stdmin,amp+stdmax,color='gray',alpha=0.2,lw=0)
        #ax.fill_between(np.linspace(-300,300),dirty_amp-dirty_std,dirty_amp+dirty_std,color='red',alpha=0.4,lw=0)
        #ax.set_ylim((-1,14))

    ax_list[0][1].yaxis.set_ticklabels([])
    #ax_list[0][2].yaxis.set_ticklabels([])
    ax_list[1][1].yaxis.set_ticklabels([])
    #ax_list[1][2].yaxis.set_ticklabels([])
    #ax_list[1][2].yaxis.set_ticks([])
    #ax_list[0][2].yaxis.set_ticks([])
    ax_list[0][1].yaxis.set_ticks([])
    ax_list[1][1].yaxis.set_ticks([])
    ax_list[0][0].set_ylabel('Peak-to-peak Amplitude',size=7)
    ax_list[1][0].set_ylabel('Peak-to-peak Amplitude',size=7)

    ax_list[0][0].set_xlabel(r'$\delta V_{S} (\%)$',size=7)
    ax_list[0][0].set_xlim((0,-20))

    #ax_list[0][1].set_xlabel(r'$\delta V_{P} (\%)$',size=7)
    #ax_list[0][1].set_xlim((-10,10))

    #ax_list[0][2].set_xlabel(r'$\delta \rho (\%)$',size=7)
    #ax_list[0][2].set_xlim((-2,8))

    ax_list[1][0].set_xlabel('Thickness (km)',size=7)
    ax_list[1][0].set_xlim((0,25))

    ax_list[1][1].set_xlabel(r'dVp/drho $\%$',size=7)

    ax_list[0][1].set_xlabel(r'Angle $(^{\circ})$',size=7)
    ax_list[0][1].set_xlim((-15,25))
    plt.figtext(0.01,0.97,'(a)',size=10)
    plt.figtext(0.50,0.97,'(b)',size=10)
    plt.figtext(0.01,0.50,'(c)',size=10)
    plt.figtext(0.50,0.50,'(d)',size=10)

    plt.tight_layout()

    return fig,ax_list

def plot_amp(amps,ax,ang,**kwargs):
    m = kwargs.get('multiply',1.)
    color = kwargs.get('color','k')
    print ang,m*np.mean(amps)
    #ax.errorbar(ang,m*np.mean(amps),
    #            yerr=(m*np.max(amps)-m*np.mean(amps)),color=color)
    ax.errorbar(ang,m*np.mean(amps),
                yerr=(np.max(amps)-np.mean(amps)),color=color)

def read_ratio(ratio_file):
    a = np.genfromtxt(ratio_file)
    data_rat = np.mean(a[:,1])
    synth_rat = np.mean(a[:,2])
    return data_rat,synth_rat

def prepare_stream(homedir,dir):
    st = obspy.read(homedir+dir+'/stv_strip.pk')
    amps = []
    for tr in st[7::]:
        amps.append(tr.data[2000:4500].max()-tr.data[2000:4500].min())
    return amps

def lowf_prepare_stream(homedir,dir):
    st = obspy.read(homedir+dir+'/stv_strip.pk')
    amps = []
    for tr in st[7::]:
        amps.append(tr.data[850:1100].max()-tr.data[850:1100].min())
    return amps

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


