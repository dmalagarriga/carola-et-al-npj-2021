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


fout_single = open('%s/Ratio_abrupt_events.txt' %(Fluo_Dir),'w')



maxima = []
minima = []
timing = []
Avg_maxima = []
maxtab,mintab = peakdet(MUA[1],.2) #maxtab,mintab = peakdet(MUA[1],.2)
for i,j in maxtab:
	i = int(i)
	timing.append(MUA[0][i])
	maxima.append(j)
maxima = array(maxima)

Ratio = float(len(maxima[np.where(maxima>np.average(maxima)+1*np.std(maxima))] ) )/(float(len(maxima)))
			

Avg_maxima.append(Ratio)
print >>fout_single, np.average(Avg_maxima),np.std(Avg_maxima)
