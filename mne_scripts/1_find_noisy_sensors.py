# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 12:44:37 2021

@author: Katharina

find noisy and flat sensors before applying maxfilter
"""

import os.path as op
import os
import sys
import numpy as np
import matplotlib

import mne
import matplotlib.pyplot as plt
import scipy
import glob

comp_chan = ['MEG2113']

# working directory
pth = r'Z:\Visual Search RFT'

os.chdir(pth+'/MNE scripts')

subj = os.listdir(pth+'/data') # subj id's

subfolders = [ f.path +'/meg/' for f in os.scandir(pth+'/data') if f.is_dir() ]

for i,n in enumerate(subfolders):
    #filenames = [f for f in os.listdir('n') if re.match(, f)]
    filenames = glob.glob(n + 'part*.fif')
    subjdir = os.path.join(pth,'results','meg','1 maxfilter',subj[i])
    
    try:
        os.makedirs(subjdir)
    except FileExistsError:
        print("subject directory exists") 
    #
    if os.path.isfile(os.path.join(subjdir,"artef_sens.txt")):
         continue
    else:
     raw = mne.io.read_raw_fif(filenames[1],allow_maxshield=True,preload=True,verbose=True)
     rawgs = raw.copy()
     raw_check = raw.copy()
     auto_noisy_chs, auto_flat_chs, auto_scores = mne.preprocessing.find_bad_channels_maxwell(raw_check, cross_talk='ct_sparse_SA.fif', calibration='sss_cal_SA.dat', return_scores=True, verbose=True)
     cols = ["red","orangered","darkred","coral","lightcoral","red","orangered","darkred","coral","lightcoral"]
     leg = []
     
     if auto_noisy_chs:
           for i,a in enumerate(auto_noisy_chs):
               rawbs = raw.copy()
               rawbs.pick_channels([a])
               # matplotlib.pyplot.psd(np.concatenate(rawbs._data),NFFT=2048,Fc=1,Fs=1000,color=cols[i],scale_by_freq=False)
               leg += [a]
              
     cols = ["blue","navy","purple","magenta","yellow"]
         
     if auto_flat_chs:
            for i,a in enumerate(auto_flat_chs):
                rawfs = raw.copy()
                rawfs.pick_channels([a])
                # matplotlib.pyplot.psd(np.concatenate(rawbs._data),NFFT=2048,Fc=1,Fs=1000,color=cols[i],scale_by_freq=False)
                leg += [a]
     
     rawgs.pick_channels(comp_chan)
     # matplotlib.pyplot.psd(np.concatenate(rawgs._data),NFFT=2048,Fc=1,Fs=1000,color="k",scale_by_freq=False)
     # matplotlib.pyplot.legend(leg+comp_chan)
     # matplotlib.pyplot.savefig(os.path.join(subjdir,'noisy_sens'))
     # matplotlib.pyplot.close()
     txtfl = open(os.path.join(subjdir,"artef_sens.txt"),"w")
     
     for sens in auto_noisy_chs:
         txtfl.write(sens + "\n")
         
     for sens in auto_flat_chs:
         txtfl.write(sens + "\n")
         txtfl.close()
         
     if auto_noisy_chs:
         del rawbs
         
     if auto_flat_chs:
         del rawfs
         
    del raw, rawgs, auto_flat_chs, auto_noisy_chs, auto_scores, subjdir, filenames
   
             

                    

                
        