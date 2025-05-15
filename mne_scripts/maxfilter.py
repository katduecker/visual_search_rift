# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 11:38:05 2021

@author: Katharina

Maxfilter on MEG data
"""

import os.path as op
import os
import sys
import numpy as np
import mne
import matplotlib.pyplot as plt
import scipy
import glob

# working directory
os.chdir('W:/Visual Search RFT/MNE scripts')

pth = 'W:/Visual Search RFT'
subj = os.listdir(pth+'/data') # subj id's

subfolders = [ f.path +'/meg/' for f in os.scandir(pth+'/data') if f.is_dir() ]

for i,n in enumerate(subfolders):
    #filenames = [f for f in os.listdir('n') if re.match(, f)]
    filenames = glob.glob(n + 'part*.fif')
    subjdir = os.path.join(pth,'results','meg','1 maxfilter',subj[i])
    
    # don't run again!
    if os.path.isfile(os.path.join(subjdir,'part'+str(1) + '_sss.fif')):
        continue
    else:
    
        # read in artefactual sensors
        sensfile = open(os.path.join(subjdir,'artef_sens.txt'), 'r')
        readsens = sensfile.readlines()
        
        # Strips the newline character
        count = 0
        artefsens = []
        for line in readsens:
            count += 1
            print(line.strip())
            artefsens += [line.strip()]
            
        del readsens
        
        sensfile.close()
    
        # maxfilter each data set
        for j,f in enumerate(filenames):
            # read in raw data
            if os.path.isfile(os.path.join(subjdir,'part'+str(j+1) + '_sss.fif')):
                continue
            else:
                raw = mne.io.read_raw_fif(f,allow_maxshield=True,preload=True,verbose=True)
        
                raw_sss = raw.copy()
                # mark bad channels
                raw_sss.info['bads'] = []
                if artefsens:
                    raw_sss.info['bads'].extend(artefsens)
                
                # fix magnetometer coil types
                raw_sss.fix_mag_coil_types()
                
                # maxfilter
                raw_sss = mne.preprocessing.maxwell_filter(raw_sss, cross_talk='ct_sparse_SA.fif', calibration='sss_cal_SA.dat')
                raw_sss.save(os.path.join(subjdir,'part'+str(j+1) + '_sss.fif'),overwrite=True)
                # plot to check
                raw.plot_psd(fmax = 100, n_fft = 1000)
                plt.savefig(os.path.join(subjdir,'part_'+str(j+1)+'_pre_max'))
                plt.close()
                raw_sss.plot_psd(fmax = 100, n_fft = 1000)
                plt.savefig(os.path.join(subjdir,'part_'+str(j+1)+'_post_max'))
                plt.close()
        
                del raw, raw_sss
            
        del artefsens, subjdir
