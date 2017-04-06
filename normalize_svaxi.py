#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : normalize_svaxi.py
Purpose : normalize SVaxi radial components to mimic Z component.
Creation Date : 06-04-2017
Last Modified : Thu 06 Apr 2017 11:37:59 AM EDT
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

def read_data():
    sts = obspy.read('/home/samhaug/work1/SP_brazil_data/q_S_sparse.pk')
    stp = obspy.read('/home/samhaug/work1/SP_brazil_data/l_S1800P_sparse.pk')
    for idx,tr in enumerate(sts):
        sts[idx].data = obspy.signal.filter.envelope(sts[idx].data)

    return sts

def read_synth(dir):
    st = obspy.read(dir+'/*0.R')
    st = seispy.filter.range_fiter(st,(56,72))
    for tr in st:
        tr.stats.sac['o'] += -202
    st = seispy.data.normalize_on_phase(st,phase=['S'],min=True)
    for tr in st:
        tr.data *= 1./np.tan(np.radians(17))

    return st


main()
