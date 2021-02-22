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
from Load_MUA_data import *
from collections import Counter

def namestr(obj, namespace):
	return [name for name in namespace if namespace[name] is obj]
    


First_set = [Exp_10_percent_Output_e0, Exp_10_percent_Output_e1, Exp_10_percent_Output_e2]
Second_set = [Exp_30_percent_Output_e0,Exp_30_percent_Output_e1,Exp_30_percent_Output_e2]
Third_set = [Exp_50_percent_Output_e0,Exp_50_percent_Output_e1,Exp_50_percent_Output_e2]

Exp_10_percent = [First_set]
Exp_30_percent = [Second_set]
Exp_50_percent = [Third_set]


####################################
######### CALCULATIONS #############
####################################

Thresh_prominence = 1.5

fout_all = open('Neurodegeneration_TH_neurite_retraction/Results_Analysis/Frequency_abrupt_events_all_SIMULATIONS.txt','w')
print >> fout_all,'CTR 1.59656677335 0.243234205985'

for MUA in Exp_10_percent:
	fout_single = open('Neurodegeneration_TH_neurite_retraction/Results_Analysis/Frequency_abrupt_events_%s.txt' %(namestr(Exp_10_percent, globals())[0]) ,'w')
	Avg_maxima = []
	for d in range(len(MUA)):
		maxima = []
		minima = []
		timing = []
		
		maxtab,mintab = peakdet(MUA[d][1],Thresh_prominence) #maxtab,mintab = peakdet(MUA[1],.2)
		for i,j in maxtab:
			i = int(i)
			timing.append(timing.append(list(MUA[d][0]).index(i)))
			maxima.append(j)
		maxima = array(maxima)
		timing = array(timing)
		#Ratio = float(len(maxima[np.where(maxima>np.average(MUA[d][1])+1*np.std(MUA[d][1]))] ) )/(float(len(maxima)))
		#Ratio = float(len(timing[np.where(maxima>np.average(MUA[d][1])+1*np.std(MUA[d][1]))] ) )/(float(MUA[d][0][-1])/(20*60))	
		
		Ratio = float(len(timing[np.where(maxima>np.average(maxima)+1*np.std(maxima))] ) )/(float(MUA[d][0][-1])/(20*60))
		Avg_maxima.append(Ratio)
		#print >>fout_single, np.average(Avg_maxima),np.std(Avg_maxima)
	
		print >>fout_single, '10%DOPA',Ratio
	print >>fout_all, '10%DOPA',np.average(Avg_maxima),np.std(Avg_maxima)
	
		
for MUA in Exp_30_percent:
	fout_single = open('Neurodegeneration_TH_neurite_retraction/Results_Analysis/Frequency_abrupt_events_%s.txt' %(namestr(Exp_30_percent, globals())[0]) ,'w')
	Avg_maxima = []
	for d in range(len(MUA)):
		maxima = []
		minima = []
		timing = []
		
		maxtab,mintab = peakdet(MUA[d][1],Thresh_prominence) #maxtab,mintab = peakdet(MUA[1],.2)
		for i,j in maxtab:
			i = int(i)
			timing.append(i)
			maxima.append(j)
		maxima = array(maxima)
		timing = array(timing)
		#Ratio = float(len(maxima[np.where(maxima>np.average(MUA[d][1])+1*np.std(MUA[d][1]))] ) )/(float(len(maxima)))
		#Ratio = float(len(timing[np.where(maxima>np.average(MUA[d][1])+1*np.std(MUA[d][1]))] ) )/(float(MUA[d][0][-1])/(20*60))		
		Ratio = float(len(timing[np.where(maxima>np.average(maxima)+1*np.std(maxima))] ) )/(float(MUA[d][0][-1])/(20*60))
		Avg_maxima.append(Ratio)
		#print >>fout_single, np.average(Avg_maxima),np.std(Avg_maxima)
		print >>fout_single,'30%DOPA',Ratio
	print >>fout_all, '30%DOPA',np.average(Avg_maxima),np.std(Avg_maxima)
	

for MUA in Exp_50_percent:
	fout_single = open('Neurodegeneration_TH_neurite_retraction/Results_Analysis/Frequency_abrupt_events_%s.txt' %(namestr(Exp_50_percent, globals())[0]) ,'w')
	Avg_maxima = []
	for d in range(len(MUA)):
		maxima = []
		minima = []
		timing = []
		
		maxtab,mintab = peakdet(MUA[d][1],Thresh_prominence) #maxtab,mintab = peakdet(MUA[1],.2)
		for i,j in maxtab:
			i = int(i)
			timing.append(i)
			maxima.append(j)
		maxima = array(maxima)
		timing = array(timing)
		#Ratio = float(len(maxima[np.where(maxima>np.average(MUA[d][1])+1*np.std(MUA[d][1]))] ) )/(float(len(maxima)))
		#Ratio = float(len(timing[np.where(maxima>np.average(MUA[d][1])+1*np.std(MUA[d][1]))] ) )/(float(MUA[d][0][-1])/(20*60))		
		Ratio = float(len(timing[np.where(maxima>np.average(maxima)+1*np.std(maxima))] ) )/(float(MUA[d][0][-1])/(20*60))
		Avg_maxima.append(Ratio)
		#print >>fout_single, np.average(Avg_maxima),np.std(Avg_maxima)
		print >>fout_single,'50%DOPA',Ratio
	print >>fout_all, '50%DOPA',np.average(Avg_maxima),np.std(Avg_maxima)