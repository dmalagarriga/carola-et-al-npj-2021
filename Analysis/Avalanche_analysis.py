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

    
############ EXPLANATION #########
# Script to calculate peaks of the MULTI-UNIT ACTIVITY
##################################


Fluo_Dir= os.getenv('Fluo_Dir')


####################################
######### CALCULATIONS #############
####################################

MUA = np.loadtxt('%s/Multi_Unit_Activity.dat' %Fluo_Dir,unpack=True)# 10 percent neurons pruning at a distance > average length

fout = open('%s/Avalanches.txt' %Fluo_Dir,'w')

maxtab,mintab = peakdet(MUA[1],0.01)

for i,j in maxtab:
	i = int(i)
	print >> fout,MUA[0][i],j
