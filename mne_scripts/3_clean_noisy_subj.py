# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 09:04:22 2021

@author: Katharina
"""

import os.path as op
import os
import sys
import numpy as np
import matplotlib
import msvcrt as m
import mne
import matplotlib.pyplot as plt
import scipy
import glob

noisub = '20210604_b4bc'
pth = 'E:\\UoB\\Proj2_Visual Search\\data'

fold = os.path.join(pth,noisub,'meg\\')
filenames = glob.glob(fold + 'part*.fif')

f = 2
raw = mne.io.read_raw_fif(filenames[f])
fig = raw.plot()
fig.canvas.key_press_event('a')

raw.save(os.path.join(fold,'part'+str(f+1)+'_raw.fif'),overwrite=True)


