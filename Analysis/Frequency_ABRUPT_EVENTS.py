#!/usr/bin/python
# -*- coding: utf-8 -*-
from matplotlib import *
use('Agg')
from scipy import *
from numpy import *
from pylab import *
import matplotlib.pyplot as plt
import math
import random
from random import uniform
import os
import glob
from scipy import signal
from scipy.signal import argrelextrema
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from peak_detection import *
from scipy import stats
from collections import Counter

    


####################################
######### CALCULATIONS #############
####################################
Fluo_Dir= os.getenv('Fluo_Dir')

MUA = np.loadtxt('%s/Multi_Unit_Activity.dat' %Fluo_Dir,unpack=True)# 10 percent neurons pruning at a distance > average length


fout_single = open('%s/Frequency_abrupt_events.txt' %(Fluo_Dir),'w')


Avg_timing = []
Avg_maxima = []
maxima = []
minima = []
timing = []

maxtab,mintab = peakdet(MUA[1],.2) # maxtab,mintab = peakdet(MUA[1],.2)
for i,j in maxtab:
	timing.append(list(MUA[0]).index(i))
	maxima.append(j)
maxima = array(maxima)
timing = array(timing)
Ratio = float(len(timing[np.where(maxima>np.average(maxima)+1*np.std(maxima))] ) )/(float(MUA[0][-1])/(20*60))
			

Avg_timing.append(Ratio)
print >>fout_single, np.average(Avg_timing),np.std(Avg_timing)
