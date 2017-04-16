#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : waves_radiation.py
Purpose : Combination of S1750P_L and radiation_plot.py
Creation Date : 15-04-2017
Last Modified : Sun 16 Apr 2017 11:50:46 AM EDT
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
import matplotlib.patches as patches
import mplstereonet
model = obspy.taup.TauPyModel(model="prem50")
from obspy.taup import TauPyModel
import matplotlib.gridspec as gridspec
import itertools

def main():
    fig,axl,axsh,axsv,axlrad,axshrad,axsvrad = setup_figure()
    wave_plot(axl,axsh,axsv)
    beachball('/home/samhaug/work1/SP_brazil_data/2007-07-21-mw60-western-brazil-4/',
              axlrad,axshrad,axsvrad)
    plt.savefig('waves_radiation.pdf')
    call('evince waves_radiation.pdf',shell=True)

def setup_figure():
    gs = gridspec.GridSpec(72,115)
    fig = plt.figure(figsize=(6,5.85))
    fig.patch.set_visible(False)

    axl = plt.subplot(gs[0:45,0:35])
    axsh = plt.subplot(gs[0:45,40:75])
    axsv = plt.subplot(gs[0:45,80:115])

    axlrad = plt.subplot(gs[52:70,0:35],projection='stereonet')
    axshrad = plt.subplot(gs[52:70,40:75],projection='stereonet')
    axsvrad = plt.subplot(gs[52:70,80:115],projection='stereonet')

    axl.set_ylim(-2,14)
    axsh.set_ylim(-2,14)
    axsv.set_ylim(-2,14)
    axl.set_title(r'$S_{1750}P$'+' displacement',size=8)
    axsh.set_title(r'$S_H$'+' displacement',size=8)
    axsv.set_title(r'$S_V$'+' displacement',size=8)
    axl.set_xlabel('Time (s)',size=8)
    axsh.set_xlabel('Time after '+r'$S_H$ (s)',size=8)
    axsv.set_xlabel('Time after '+r'$S_V$ (s)',size=8)

    for b in [axl,axsh,axsv]:
        b.spines['top'].set_visible(False)
        b.spines['right'].set_visible(False)
        b.spines['left'].set_visible(False)
        b.yaxis.set_ticklabels([])
        b.yaxis.set_ticks([])
        b.tick_params(axis='both', which='major', labelsize=6)
        b.xaxis.set_ticks_position('bottom')

    plt.figtext(0.08,0.94,'(a)',size=7)
    plt.figtext(0.38,0.94,'(b)',size=7)
    plt.figtext(0.65,0.94,'(c)',size=7)

    plt.figtext(0.08,0.30,'(d)',size=7)
    plt.figtext(0.38,0.30,'(e)',size=7)
    plt.figtext(0.65,0.30,'(f)',size=7)

    plt.figtext(0.23,0.25,r'$S_{1750}P$',size=6,color='black')
    plt.figtext(0.46,0.26,'$S$',size=6,color='black')
    plt.figtext(0.45,0.23,'$S_{1750}P$',size=6,color='black')
    plt.figtext(0.73,0.26,'$S$',size=6,color='white')
    plt.figtext(0.72,0.23,'$S_{1750}P$',size=6,color='white')

    return fig,axl,axsh,axsv,axlrad,axshrad,axsvrad

def wave_plot(axl,axsh,axsv):
    '''
    Plot data of L component S1800P conversion alongside SV and SH
    data.
    '''
    stz = obspy.read('/home/samhaug/work1/SP_brazil_data/BHZ.pk')
    stn = obspy.read('/home/samhaug/work1/SP_brazil_data/BHN.pk')
    ste = obspy.read('/home/samhaug/work1/SP_brazil_data/BHE.pk')

    stz = seispy.filter.gimp_filter(stz)
    stz.filter('highpass',freq=1/50.)
    ste = seispy.filter.gimp_filter(ste)
    ste.filter('highpass',freq=1/50.)
    stn = seispy.filter.gimp_filter(stn)
    stn.filter('highpass',freq=1/50.)

    for tr in stz:
        tr.stats.az = tr.stats.sac['az']
    for tr in stn:
        tr.stats.az = tr.stats.sac['az']
    for tr in ste:
        tr.stats.az = tr.stats.sac['az']

    l,fq,ft = seispy.data.rotate_phase(stz,stn,ste,['S1850P'])
    fuck,q,t = seispy.data.rotate_phase(stz,stn,ste,['S'])


    l = seispy.data.align_on_phase(l,phase =['P'],window_tuple=(-20,20))
    q = seispy.data.align_on_phase(q,phase =['S'],window_tuple=(-20,20),min=False)
    t = seispy.data.align_on_phase(t,phase =['S'],window_tuple=(-20,20))

    for tr in l:
        tr.stats.az = tr.stats.sac['az']
    for tr in q:
        tr.stats.az = tr.stats.sac['az']
    for tr in t:
        tr.stats.az = tr.stats.sac['az']

    l.sort(['az'])
    q.sort(['az'])
    t.sort(['az'])

    skip_list = ['OR074','V19A']

    for idx, tr in enumerate(l):
        name = tr.stats.station
        if name in skip_list:
            l.remove(tr)
            q.remove(q[idx])
            t.remove(t[idx])

    l.sort(['az'])

    lmax = []
    for idx,tr in enumerate(l):
        name = tr.stats.station
        #print name
        #This is the same as windowing on S1775P
        a = seispy.data.phase_window(tr,['S1750P'],window=(-8,12))
        tl = np.linspace(-10,10,num=len(a.data))
        z = np.polyfit(tl,a.data,4)
        p = np.poly1d(z)
        a.data += -1*p(tl)
        #print np.mean(a.data)
        a.data += -1*np.mean(a.data)
        a_max = np.abs(a.data).max()
        a.normalize()
        a.data *= 0.7
        gcarc = tr.stats.sac['gcarc']
        az = tr.stats.sac['az']
        h = tr.stats.sac['evdp']

        lmax.append(a_max)
        axl.plot(tl,2.0*idx+a.data,color='k',lw=0.5)
        axl.text(-9,2.0*idx+0.48,r'$\Delta =$'+str(np.round(gcarc,1))+r'$^{\circ}$',size=6)
        axl.text(5,2.0*idx+0.48,r'$az =$'+str(np.round(az,2))+r'$^{\circ}$',size=6)
        axl.text(-14,2*idx+0.48,name,size=6)
        #axl.text(1.5,2*idx+0.44,str(int(a_max))+'$nm$',size=6)

    for idx,tr in enumerate(q):
        name = tr.stats.station
        a = seispy.data.phase_window(tr,['S'],window=(-30,30))
        env = obspy.signal.filter.envelope(a.data)
        a_max = env.max()
        #a_max = np.abs(a.data).max()
        a.normalize()
        a.data *= 0.9
        gcarc = tr.stats.sac['gcarc']
        az = tr.stats.sac['az']
        h = tr.stats.sac['evdp']
        tl = np.linspace(-30,30,num=len(a.data))
        axsv.plot(tl,2.0*idx+a.data,color='k',lw=0.5)
        #axsv.text(-9,2*idx+0.44,r'$\Delta =$'+str(np.round(gcarc,1))+r'$^{\circ}$',size=6)
        #axsv.text(-14,2*idx+0.44,name,size=6)
        axsv.set_xlim(-30,30)
        axsv.text(-20,2*idx+0.44,str(int(a_max))+'$nm$',size=6)
        axsv.text(20,2*idx+0.44,r'$\alpha=$'+str(np.round(int(a_max)/lmax[idx],1)),size=6)

    for idx,tr in enumerate(t):
        name = tr.stats.station
        a = seispy.data.phase_window(tr,['S'],window=(-30,30))
        env = obspy.signal.filter.envelope(a.data)
        a_max = env.max()
        #a_max = np.abs(a.data).max()
        a.normalize()
        a.data *= 0.9
        gcarc = tr.stats.sac['gcarc']
        az = tr.stats.sac['az']
        h = tr.stats.sac['evdp']
        tl = np.linspace(-30,30,num=len(a.data))
        axsh.plot(tl,2.0*idx+a.data,color='k',lw=0.5)
        #axsh.text(-9,2*idx+0.44,r'$\Delta =$'+str(np.round(gcarc,1))+r'$^{\circ}$',size=6)
        #axsh.text(7,2*idx+0.44,r'$az =$'+str(np.round(az,1))+r'$^{\circ}$',size=6)
        #axsh.text(-14,2*idx+0.44,name,size=6)
        axsh.set_xlim(-30,30)
        axsh.text(-20,2*idx+0.44,str(int(a_max))+'$nm$',size=6)
        axsh.text(20,2*idx+0.44,r'$\alpha=$'+str(np.round(int(a_max)/lmax[idx],1)),size=6)

    p = patches.Rectangle((-2.0, -1.0), 5.0, 14.0,alpha=0.1,color='k',linewidth=None,lw=0)
    axl.add_patch(p)

def beachball(homedir,ax1,ax2,ax3):
    def main(homedir):
        st = obspy.read(homedir+'sparse_R.pk')
        Mxyz = cmt2mxyz(homedir+'CMTSOLUTION')
        p_coords = stereonet_coords(Mxyz,find_p_amp)
        sv_coords = stereonet_coords(Mxyz,find_sv_amp)
        sh_coords = stereonet_coords(Mxyz,find_sh_amp)
        plot_coords(p_coords,ax1)
        plot_coords(sh_coords,ax2)
        plot_coords(sv_coords,ax3)
        ray_coords_S1750P = get_ray_coordinates(st,['S1750P'])
        ray_coords_S = get_ray_coordinates(st,['S'])
        for ii in ray_coords_S1750P:
            ax1.pole(90+ii[0],ii[1],markersize=2.0,color='black',mew=0.)
            ax2.pole(90+ii[0],ii[1],markersize=2.0,color='black',mew=0.)
            ax3.pole(90+ii[0],ii[1],markersize=2.0,color='white',mew=0.)
        for ii in ray_coords_S:
            ax2.pole(90+ii[0],ii[1],markersize=2.0,color='black',mew=0.)
            ax3.pole(90+ii[0],ii[1],markersize=2.0,color='white',mew=0.)


    def stereonet_coords(Mxyz,func):
        r = np.linspace(0,90,num=190)
        t = np.linspace(0,360,num=360)
        coords = list(itertools.product(r,t))
        mag = []
        for ii in range(0,len(coords)):
            mag.append(func(Mxyz,np.radians(coords[ii][1]),np.radians(coords[ii][0])))
        for ii in range(len(mag)):
            if mag[ii] <= 0:
                coords[ii] = 0
        coords = [i for i in coords if i != 0]
        coords = np.array(coords)
        return coords

    def plot_coords(coords,ax):
        ax.pole(90+coords[:,1],coords[:,0],color='k',marker='o',alpha=0.2,rasterized=True)
        ax.set_azimuth_ticklabels([])

    def get_ray_coordinates(st,phase):
        model = TauPyModel(model='prem50')
        ray_coord = []
        for tr in st:
            evdp = tr.stats.sac['evdp']
            gcarc = tr.stats.sac['gcarc']
            arrivals = model.get_travel_times(source_depth_in_km=evdp,
                                              distance_in_degree=gcarc,
                                              phase_list=phase)
            ang = arrivals[0].takeoff_angle
            ray_coord.append([tr.stats.sac['az'],ang])
        return ray_coord

    main(homedir)

def cmt2mxyz(file):
    '''modification on what jeroen did in orig_cmt2mxyz'''
    '''See box 4.4, aki richards page 113'''
    a = np.genfromtxt(file,skip_header=7)[:,1]
    Mrr = a[0]
    Mtt = a[1]
    Mpp = a[2]
    Mrt = a[3]
    Mrp = a[4]
    Mtp = a[5]
    Mxyz = np.array([[Mrr,Mrt,Mrp],[Mrt,Mtt,Mtp],[Mrp,Mtp,Mpp]])
    new_Mxyz = Mxyz.copy()
    new_Mxyz[0,0] = Mxyz[1,1]
    new_Mxyz[0,1] = -Mxyz[1,2]
    new_Mxyz[0,2] = Mxyz[0,1]

    new_Mxyz[1,0] = -Mxyz[1,2]
    new_Mxyz[1,1] = Mxyz[2,2]
    new_Mxyz[1,2] = -Mxyz[0,2]

    new_Mxyz[2,0] = Mxyz[0,1]
    new_Mxyz[2,1] = -Mxyz[0,2]
    new_Mxyz[2,2] = Mxyz[0,0]
    return new_Mxyz

def find_p_amp(Mxyz,az,toa):
    '''az:azimuthal direction, toa:takeoff angle'''
    Ap = 0
    si,ci = np.sin(toa),np.cos(toa)
    sf,cf = np.sin(az),np.cos(az)
    g = [si*cf, si*sf, ci]
    v = [ci*cf, ci*sf,-si]
    h = [-sf, cf, 0.]
    #for k1 in xrange(3):
    #    for k2 in xrange(3):
    #        Ap += Mxyz[k1][k2]*g[k1]*g[k2]
    Ap = np.dot(np.array(g),np.dot(np.array(Mxyz),np.array(g)))
    return Ap

def find_sv_amp(Mxyz,az,toa):
    Asv = 0
    si,ci = np.sin(toa),np.cos(toa)
    sf,cf = np.sin(az),np.cos(az)
    g = [si*cf, si*sf, ci]
    v = [ci*cf, ci*sf,-si]
    h = [-sf, cf, 0.]
    #for k1 in xrange(3):
    #    for k2 in xrange(3):
    #        Asv += Mxyz[k1][k2]*g[k1]*v[k2]
    Asv = np.dot(np.array(v),np.dot(np.array(Mxyz),np.array(g)))
    return Asv

def find_sh_amp(Mxyz,az,toa):
    Ash = 0
    si,ci = np.sin(toa),np.cos(toa)
    sf,cf = np.sin(az),np.cos(az)
    g = [si*cf, si*sf, ci]
    v = [ci*cf, ci*sf,-si]
    h = [-sf, cf, 0.]
    #for k1 in xrange(3):
    #    for k2 in xrange(3):
    #        Ash += Mxyz[k1][k2]*g[k1]*h[k2]
    Ash = np.dot(np.array(h),np.dot(np.array(Mxyz),np.array(g)))
    return Ash


#wave_plot()
#plt.savefig('S1800P_L.pdf')
#call('evince S1800P_L.pdf',shell=True)


main()
