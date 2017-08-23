#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : radiation_plot.py
Purpose : Find P,SV,SH amplitudes for different takeoff angles/azimuths.
          page 108-109. Aki/richards. equation (4.88). chapter : Elastic waves
          from a point dislocation source
Creation Date : 11-04-2017
Last Modified : Tue 11 Apr 2017 11:53:39 AM EDT
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
from sys import stdout
from obspy.taup import TauPyModel
import mplstereonet
import itertools
from deco import concurrent,synchronized

def main():
    beachball('/home/samhaug/work1/SP_brazil_data/2007-07-21-mw60-western-brazil-4/')

def beachball(homedir):
    def main(homedir):
        st = obspy.read(homedir+'sparse_R.pk')
        Mxyz = cmt2mxyz(homedir+'CMTSOLUTION')
        fig = plt.figure(figsize=(5,3))
        ax1 = fig.add_subplot(131, projection='stereonet')
        ax1.set_title('P',size=7)
        ax2 = fig.add_subplot(132, projection='stereonet')
        ax2.set_title('SV',size=7)
        ax3 = fig.add_subplot(133, projection='stereonet')
        ax3.set_title('SH',size=7)
        p_coords = stereonet_coords(Mxyz,find_p_amp)
        sv_coords = stereonet_coords(Mxyz,find_sv_amp)
        sh_coords = stereonet_coords(Mxyz,find_sh_amp)
        plot_coords(p_coords,ax1)
        plot_coords(sv_coords,ax2)
        plot_coords(sh_coords,ax3)

        ray_coords = get_ray_coordinates(st,['S1800P'])
        ax1.pole(90+ray_coords[1][0],ray_coords[1][1],markersize=3.0,
                 color='limegreen',mew=0.)
        ax2.pole(90+ray_coords[1][0],ray_coords[1][1],markersize=3.0,
                 color='limegreen',mew=0.)
        ax3.pole(90+ray_coords[1][0],ray_coords[1][1],markersize=3.0,
                 color='limegreen',mew=0.)

        ray_coords = get_ray_coordinates(st,['S'])
        ax1.pole(90+ray_coords[1][0],ray_coords[1][1],markersize=3.0,
                 color='red',mew=0.)
        ax2.pole(90+ray_coords[1][0],ray_coords[1][1],markersize=3.0,
                 color='red',mew=0.)
        ax3.pole(90+ray_coords[1][0],ray_coords[1][1],markersize=3.0,
                 color='red',mew=0.)

        ray_coords = get_ray_coordinates(st,['P'])
        ax1.pole(90+ray_coords[1][0],ray_coords[1][1],markersize=3.0,
                 color='yellow',mew=0.)
        ax2.pole(90+ray_coords[1][0],ray_coords[1][1],markersize=3.0,
                 color='yellow',mew=0.)
        ax3.pole(90+ray_coords[1][0],ray_coords[1][1],markersize=3.0,
                 color='yellow',mew=0.)

        plt.savefig('radiation_plot.pdf')
        call('evince radiation_plot.pdf',shell=True)

    def stereonet_coords(Mxyz,func):
        r = np.linspace(0,90,num=190)
        t = np.linspace(0,360,num=460)
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
        ax.pole(90+coords[:,1],coords[:,0],color='k',marker='o',rasterized=True)
        ax.set_azimuth_ticklabels([])

    def get_ray_coordinates(st,phase_list):
        model = TauPyModel(model='prem50')
        ray_coord = []
        for tr in st:
            evdp = tr.stats.sac['evdp']
            gcarc = tr.stats.sac['gcarc']
            arrivals = model.get_travel_times(source_depth_in_km=evdp,
                                              distance_in_degree=gcarc,
                                              phase_list=phase_list)
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

main()
