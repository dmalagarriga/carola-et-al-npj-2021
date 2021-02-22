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
from filtering import *
from peak_detection import *
from scipy import stats
from Load_all_data import *
from collections import Counter

def namestr(obj, namespace):
    return [name for name in namespace if namespace[name] is obj]
    
    
################# EXPLANATION ####################
# Script to compute the Kolmogorov-Smirnov test on the probability distribution functions 
# of Multi-Unit activity and Single-Unit activity of the experiments. This method is used
# to classify the experiments by similarity of their distribution functions. The pooling 
# of the experiments will be done between those who are more similar in MUA and SUA distributions.
# In order to show prominent differences between days, experiments showing larger KS difference
# will be taken as standard experiments.
#################################################

Pooled_1 = 0
Pooled = 0
Not_Pooled = 1

Day_35 = 0
Day_50 = 1
Day_80 = 0

if Not_Pooled == 1:

	# Ouput data Not Pooled 
	if Day_35 == 1:
		fout_D35_Not_Pooled = open('All_Experiments_Data/MUA_Not_Pooled/MUA_D35_Final_2.txt','w')
		#fout_D35_Not_Pooled = open('All_Experiments_Data/MUA_Not_Pooled/MUA_D35_Final_1.txt','w')
		#fout_D35_Not_Pooled = open('All_Experiments_Data/MUA_Not_Pooled/MUA_D35_MORE_EXP.txt','w')
		#fout_D35_Not_Pooled = open('All_Experiments_Data/MUA_Not_Pooled/MUA_D35_LESS_EXP.txt','w')
		#fout_D35_Not_Pooled = open('All_Experiments_Data/MUA_Not_Pooled/MUA_D35.txt','w')
		#fout_HIST_D35_Not_Pooled = open('All_Experiments_Data/MUA_Not_Pooled/HIST_MUA_D35.txt','w')
	elif Day_50 == 1:
		fout_D50_Not_Pooled = open('All_Experiments_Data/MUA_Not_Pooled/MUA_D50_Final_2.txt','w')
		#fout_D50_Not_Pooled = open('All_Experiments_Data/MUA_Not_Pooled/MUA_D50_Final_1.txt','w')
		#fout_D50_Not_Pooled = open('All_Experiments_Data/MUA_Not_Pooled/MUA_D50_MORE_EXP.txt','w')
		#fout_D50_Not_Pooled = open('All_Experiments_Data/MUA_Not_Pooled/MUA_D50_LESS_EXP.txt','w')
		#fout_D50_Not_Pooled = open('All_Experiments_Data/MUA_Not_Pooled/MUA_D50.txt','w')
		#fout_HIST_D50_Not_Pooled = open('All_Experiments_Data/MUA_Not_Pooled/HIST_MUA_D50.txt','w')
	elif Day_80 == 1:
		fout_D80_Not_Pooled = open('All_Experiments_Data/MUA_Not_Pooled/MUA_D80_Final_2.txt','w')
		#fout_D80_Not_Pooled = open('All_Experiments_Data/MUA_Not_Pooled/MUA_D80_Final_1.txt','w')
		#fout_D80_Not_Pooled = open('All_Experiments_Data/MUA_Not_Pooled/MUA_D80_MORE_EXP.txt','w')
		#fout_D80_Not_Pooled = open('All_Experiments_Data/MUA_Not_Pooled/MUA_D80_LESS_EXP.txt','w')
		#fout_D80_Not_Pooled = open('All_Experiments_Data/MUA_Not_Pooled/MUA_D80.txt','w')
		#fout_HIST_D80_Not_Pooled = open('All_Experiments_Data/MUA_Not_Pooled/HIST_MUA_D80.txt','w')
if Pooled_1 == 1:
	# Ouput data Pooled 
	if Day_35 == 1:
		fout_D35_Pooled = open('All_Experiments_Data/MUA_Pooled/MUA_D35_Final_2.txt','w')
		#fout_D35_Pooled = open('All_Experiments_Data/MUA_Pooled/MUA_D35_Final_1.txt','w')
		#fout_D35_Pooled = open('All_Experiments_Data/MUA_Pooled/MUA_D35_MORE_EXP.txt','w')
		#fout_D35_Pooled = open('All_Experiments_Data/MUA_Pooled/MUA_D35.txt','w')
		#fout_HIST_D35_Pooled = open('All_Experiments_Data/MUA_Pooled/HIST_MUA_D35.txt','w')
	elif Day_50 == 1:
		fout_D50_Pooled = open('All_Experiments_Data/MUA_Pooled/MUA_D50_Final_2.txt','w')
		#fout_D50_Pooled = open('All_Experiments_Data/MUA_Pooled/MUA_D50_Final_1.txt','w')
		#fout_D50_Pooled = open('All_Experiments_Data/MUA_Pooled/MUA_D50_MORE_EXP.txt','w')
		#fout_D50_Pooled = open('All_Experiments_Data/MUA_Pooled/MUA_D50.txt','w')
		#fout_HIST_D50_Pooled = open('All_Experiments_Data/MUA_Pooled/HIST_MUA_D50.txt','w')
	elif Day_80 == 1:
		fout_D80_Pooled = open('All_Experiments_Data/MUA_Pooled/MUA_D80_Final_2.txt','w')
		#fout_D80_Pooled = open('All_Experiments_Data/MUA_Pooled/MUA_D80_Final_1.txt','w')
		#fout_D80_Pooled = open('All_Experiments_Data/MUA_Pooled/MUA_D80_MORE_EXP.txt','w')
		#fout_D80_Pooled = open('All_Experiments_Data/MUA_Pooled/MUA_D80.txt','w')
		#fout_HIST_D80_Pooled = open('All_Experiments_Data/MUA_Pooled/HIST_MUA_D80.txt','w')

if Pooled  == 1:
	# Ouput data Pooled 
	if Day_35 == 1:
		fout_D35_Pooled = open('All_Experiments_Data/MUA_Pooled/MUA_D35_Final_2.txt','w')
		#fout_D35_Pooled = open('All_Experiments_Data/MUA_Pooled/MUA_D35_Final_1.txt','w')
		#fout_D35_Pooled = open('All_Experiments_Data/MUA_Pooled/MUA_D35_MORE_EXP_2.txt','w')
		#fout_D35_Pooled = open('All_Experiments_Data/MUA_Pooled/MUA_D35.txt','w')
		#fout_HIST_D35_Pooled = open('All_Experiments_Data/MUA_Pooled/HIST_MUA_D35.txt','w')
	elif Day_50 == 1:
		fout_D50_Pooled = open('All_Experiments_Data/MUA_Pooled/MUA_D50_Final_2.txt','w')
		#fout_D50_Pooled = open('All_Experiments_Data/MUA_Pooled/MUA_D50_Final_1.txt','w')
		#fout_D50_Pooled = open('All_Experiments_Data/MUA_Pooled/MUA_D50_MORE_EXP_2.txt','w')
		#fout_D50_Pooled = open('All_Experiments_Data/MUA_Pooled/MUA_D50.txt','w')
		#fout_HIST_D50_Pooled = open('All_Experiments_Data/MUA_Pooled/HIST_MUA_D50.txt','w')
	elif Day_80 == 1:
		fout_D80_Pooled = open('All_Experiments_Data/MUA_Pooled/MUA_D80_Final_2.txt','w')
		#fout_D80_Pooled = open('All_Experiments_Data/MUA_Pooled/MUA_D80_Final_1.txt','w')
		#fout_D80_Pooled = open('All_Experiments_Data/MUA_Pooled/MUA_D80_MORE_EXP_2.txt','w')
		#fout_D80_Pooled = open('All_Experiments_Data/MUA_Pooled/MUA_D80.txt','w')
		#fout_HIST_D80_Pooled = open('All_Experiments_Data/MUA_Pooled/HIST_MUA_D80.txt','w')

####################################
######### CALCULATIONS #############
####################################


# Selected data ('Cooked results') (Final_2)
SP11_D35 = [SP11_D35_MUA_1,SP11_D35_MUA_2]
SP12_D35 = [SP12_D35_MUA_2]
SP12_D35_ED = [SP12_ED_D35_MUA_1_n1,SP12_ED_D35_MUA_2_n1,SP12_ED_D35_MUA_3_n1,SP12_ED_D35_MUA_1_n2,SP12_ED_D35_MUA_2_n2,SP12_ED_D35_MUA_3_n2] 
SP13_D35 = [SP13_D35_MUA_1,SP13_D35_MUA_2]

SP11_D50 = [SP11_D50_MUA_1,SP11_D50_MUA_2]
SP12_D50 = [SP12_D50_MUA_1,SP12_D50_MUA_2]
SP12_D50_ED = [SP12_ED_D50_MUA_2_n1,SP12_ED_D50_MUA_1_n2,SP12_ED_D50_MUA_2_n2,SP12_ED_D50_MUA_3_n2]
SP13_D50 = [SP13_D50_MUA_1,SP13_D50_MUA_2]

SP11_D80 = [SP11_D80_MUA_2,SP11_D80_MUA_3]
SP12_D80 = [SP12_D80_MUA_1,SP12_D80_MUA_2,SP12_D80_MUA_3,SP12_D80_MUA_4,SP12_D80_MUA_E1,SP12_D80_MUA_E2]
SP12_D80_ED = [SP12_ED_D80_MUA_1_n1,SP12_ED_D80_MUA_2_n1]
SP13_D80 = [SP13_D80_MUA_1,SP13_D80_MUA_2]

'''
# Selected data ('Cooked results') (Final_1)
SP11_D35 = [SP11_D35_MUA_1,SP11_D35_MUA_2]
SP12_D35 = [SP12_D35_MUA_1,SP12_D35_MUA_2]
SP12_D35_ED = [SP12_ED_D35_MUA_1_n1,SP12_ED_D35_MUA_2_n1,SP12_ED_D35_MUA_3_n1,SP12_ED_D35_MUA_1_n2,SP12_ED_D35_MUA_2_n2,SP12_ED_D35_MUA_3_n2] 
SP13_D35 = [SP13_D35_MUA_1,SP13_D35_MUA_2]

SP11_D50 = [SP11_D50_MUA_1,SP11_D50_MUA_2]
SP12_D50 = [SP12_D50_MUA_1,SP12_D50_MUA_2,SP12_D50_MUA_E1,SP12_D50_MUA_E1_CNQX,SP12_D50_MUA_E1_l_dopa,SP12_D50_MUA_E2]
SP12_D50_ED = [SP12_ED_D50_MUA_2_n1,SP12_ED_D50_MUA_1_n2,SP12_ED_D50_MUA_2_n2,SP12_ED_D50_MUA_3_n2]
SP13_D50 = [SP13_D50_MUA_1,SP13_D50_MUA_2]

SP11_D80 = [SP11_D80_MUA_2,SP11_D80_MUA_3]
SP12_D80 = [SP12_D80_MUA_1,SP12_D80_MUA_2,SP12_D80_MUA_3,SP12_D80_MUA_4,SP12_D80_MUA_E1,SP12_D80_MUA_E2]
SP12_D80_ED = [SP12_ED_D80_MUA_1_n1,SP12_ED_D80_MUA_1_n1]
SP13_D80 = [SP13_D80_MUA_1,SP13_D80_MUA_2]


# Adding Drug_TestResults (_MORE_EXP)
SP11_D35 = [SP11_D35_MUA_1,SP11_D35_MUA_2]
SP12_D35 = [SP12_D35_MUA_1,SP12_D35_MUA_2]
SP12_D35_ED = [SP12_ED_D35_MUA_1_n1,SP12_ED_D35_MUA_2_n1,SP12_ED_D35_MUA_3_n1,SP12_ED_D35_MUA_1_n2,SP12_ED_D35_MUA_2_n2,SP12_ED_D35_MUA_3_n2] 
SP13_D35 = [SP13_D35_MUA_1,SP13_D35_MUA_2]

SP11_D50 = [SP11_D50_MUA_1,SP11_D50_MUA_2,SP11_D50_MUA_3,SP11_D50_MUA_E1_CNQX,SP11_D50_MUA_E1_l_dopa,SP11_D50_MUA_E2,SP11_D50_MUA_E3]
SP12_D50 = [SP12_D50_MUA_1,SP12_D50_MUA_2,SP12_D50_MUA_E1,SP12_D50_MUA_E1_CNQX,SP12_D50_MUA_E1_l_dopa,SP12_D50_MUA_E2]
SP12_D50_ED = [SP12_ED_D50_MUA_1_n1,SP12_ED_D50_MUA_2_n1,SP12_ED_D50_MUA_1_n2,SP12_ED_D50_MUA_2_n2,SP12_ED_D50_MUA_3_n2]
SP13_D50 = [SP13_D50_MUA_1,SP13_D50_MUA_2]

SP11_D80 = [SP11_D80_MUA_2,SP11_D80_MUA_3,SP11_D80_MUA_4,SP11_D80_MUA_E1,SP11_D80_MUA_E1_CNQX,SP11_D80_MUA_E1_l_dopa,SP11_D80_MUA_E2,SP11_D80_MUA_E2_CNQX,SP11_D80_MUA_E2_l_dopa]
SP12_D80 = [SP12_D80_MUA_1,SP12_D80_MUA_2,SP12_D80_MUA_3,SP12_D80_MUA_4,SP12_D80_MUA_E1,SP12_D80_MUA_E2]
SP12_D80_ED = [SP12_ED_D80_MUA_1_n1,SP12_ED_D80_MUA_1_n1]
SP13_D80 = [SP13_D80_MUA_1,SP13_D80_MUA_2]



# Removing 1st set of experiments from the set below (_LESS_EXP)
SP11_D35 = [SP11_D35_MUA_1]
SP12_D35 = [SP12_D35_MUA_2]
SP12_D35_ED = [SP12_ED_D35_MUA_1_n1,SP12_ED_D35_MUA_2_n1,SP12_ED_D35_MUA_3_n1,SP12_ED_D35_MUA_1_n2,SP12_ED_D35_MUA_3_n2] 
SP13_D35 = [SP13_D35_MUA_1,SP13_D35_MUA_2]

SP11_D50 = [SP11_D50_MUA_1,SP11_D50_MUA_2]
SP12_D50 = [SP12_D50_MUA_1,SP12_D50_MUA_2]
SP12_D50_ED = [SP12_ED_D50_MUA_2_n1,SP12_ED_D50_MUA_1_n2,SP12_ED_D50_MUA_2_n2,SP12_ED_D50_MUA_3_n2]
SP13_D50 = [SP13_D50_MUA_1,SP13_D50_MUA_2]

SP11_D80 = [SP11_D80_MUA_2,SP11_D80_MUA_3]
SP12_D80 = [SP12_D80_MUA_1,SP12_D80_MUA_2,SP12_D80_MUA_3,SP12_D80_MUA_4]
SP12_D80_ED = [SP12_ED_D80_MUA_1_n1,SP12_ED_D80_MUA_1_n1]
SP13_D80 = [SP13_D80_MUA_1,SP13_D80_MUA_2]


#  Prior to adding Drug_TestResults
SP11_D35 = [SP11_D35_MUA_1,SP11_D35_MUA_2]
SP12_D35 = [SP12_D35_MUA_1,SP12_D35_MUA_2]
SP12_D35_ED = [SP12_ED_D35_MUA_1_n1,SP12_ED_D35_MUA_2_n1,SP12_ED_D35_MUA_3_n1,SP12_ED_D35_MUA_1_n2,SP12_ED_D35_MUA_2_n2,SP12_ED_D35_MUA_3_n2] 
SP13_D35 = [SP13_D35_MUA_1,SP13_D35_MUA_2]

SP11_D50 = [SP11_D50_MUA_1,SP11_D50_MUA_2,SP11_D50_MUA_3]
SP12_D50 = [SP12_D50_MUA_1,SP12_D50_MUA_2]
SP12_D50_ED = [SP12_ED_D50_MUA_1_n1,SP12_ED_D50_MUA_2_n1,SP12_ED_D50_MUA_1_n2,SP12_ED_D50_MUA_2_n2,SP12_ED_D50_MUA_3_n2]
SP13_D50 = [SP13_D50_MUA_1,SP13_D50_MUA_2]

SP11_D80 = [SP11_D80_MUA_2,SP11_D80_MUA_3,SP11_D80_MUA_4]
SP12_D80 = [SP12_D80_MUA_1,SP12_D80_MUA_2,SP12_D80_MUA_3,SP12_D80_MUA_4]
SP12_D80_ED = [SP12_ED_D80_MUA_1_n1,SP12_ED_D80_MUA_1_n1]
SP13_D80 = [SP13_D80_MUA_1,SP13_D80_MUA_2]
'''

D35 = [SP11_D35,SP12_D35,SP12_D35_ED,SP13_D35]

D50 = [SP11_D50,SP12_D50,SP12_D50_ED,SP13_D50]

D80 = [SP11_D80,SP12_D80,SP12_D80_ED,SP13_D80]

#print len(SP11_D35_SUA_1[0]),len(SP11_D35_SUA_2[1])

if Not_Pooled == 1:
	################################################################################################
	######################## CALCULATION WITHOUT POOLING TOGETHER SIMILAR EXPERIMENTS ##############
	################################################################################################
	if Day_35 == 1:
		for day in D35:
			Pooled_Experiments = []
			if len(day)>1:
				for i in range(len(day)):
					for j in range(i):
				
						weights1 = np.ones_like(day[i][1])/float(len(day[i][1]))
						weights2 = np.ones_like(day[j][1])/float(len(day[j][1]))
				
						n1, bins1, patches1 = plt.hist(day[i][1], weights=weights1,bins=200)#,50, normed=1, facecolor='green', alpha=0.75)
						n2, bins2, patches2 = plt.hist(day[j][1], weights=weights2,bins=200)#,50, normed=1, facecolor='green', alpha=0.75)
				
						#print stats.ks_2samp(n1,n2)[1],namestr(day[i], globals()),namestr(day[j], globals()),namestr(day, globals())
						if stats.ks_2samp(n1,n2)[1]>= 0.0: # NO POOLING!!
							Pooled_Experiments.append(day[i][1].mean()/2.0)
							Pooled_Experiments.append(day[j][1].mean()/2.0)
						
						#	print i,j,day[i][1].mean(),day[j][1].mean(),namestr(day[i], globals()),namestr(day[j], globals())
		
				if Pooled_Experiments != []:
					if namestr(day, globals())[0] == 'day':
						print >> fout_D35_Not_Pooled,namestr(day, globals())[1],np.average(Pooled_Experiments), np.std(Pooled_Experiments)
					else:
						print >> fout_D35_Not_Pooled,namestr(day, globals())[0],np.average(Pooled_Experiments), np.std(Pooled_Experiments)
				
			elif len(day)<=1:
				for i in range(len(day)):
					Pooled_Experiments.append(day[i][1].mean())
				
				if Pooled_Experiments != []:
					if namestr(day, globals())[0] == 'day':
						print >> fout_D35_Not_Pooled,namestr(day, globals())[1],np.average(Pooled_Experiments), np.std(Pooled_Experiments)
					else:
						print >> fout_D35_Not_Pooled,namestr(day, globals())[0],np.average(Pooled_Experiments), np.std(Pooled_Experiments)
	if Day_50 == 1:
		for day in D50:
			Pooled_Experiments = []
			
			if len(day)>1:
				
				for i in range(len(day)):
					for j in range(i):
				
						weights1 = np.ones_like(day[i][1])/float(len(day[i][1]))
						weights2 = np.ones_like(day[j][1])/float(len(day[j][1]))
				
						n1, bins1, patches1 = plt.hist(day[i][1], weights=weights1,bins=200)#,50, normed=1, facecolor='green', alpha=0.75)
						n2, bins2, patches2 = plt.hist(day[j][1], weights=weights2,bins=200)#,50, normed=1, facecolor='green', alpha=0.75)
				
						#print stats.ks_2samp(n1,n2)[1],namestr(day[i], globals()),namestr(day[j], globals()),namestr(day, globals())
						if stats.ks_2samp(n1,n2)[1]>= 0.0: # NO POOLING!!
							Pooled_Experiments.append(day[i][1].mean()/2.0)
							Pooled_Experiments.append(day[j][1].mean()/2.0)
						
						
						#print i,j,day[i][1].mean(),day[j][1].mean(),namestr(day[i], globals()),namestr(day[j], globals())
		
				if Pooled_Experiments != []:
					if namestr(day, globals())[0] == 'day':
						print >> fout_D50_Not_Pooled,namestr(day, globals())[1],np.average(Pooled_Experiments), np.std(Pooled_Experiments)
					else:
						print >> fout_D50_Not_Pooled,namestr(day, globals())[0],np.average(Pooled_Experiments), np.std(Pooled_Experiments)
				
			elif len(day)<=1:
				for i in range(len(day)):
					Pooled_Experiments.append(day[i][1].mean())
				
				if Pooled_Experiments != []:
					if namestr(day, globals())[0] == 'day':
						print >> fout_D50_Not_Pooled,namestr(day, globals())[1],np.average(Pooled_Experiments), np.std(Pooled_Experiments)
					else:
						print >> fout_D50_Not_Pooled,namestr(day, globals())[0],np.average(Pooled_Experiments), np.std(Pooled_Experiments)
	
	if Day_80 == 1:
		for day in D80:
			Pooled_Experiments = []
			
			if len(day)>1:
				
				for i in range(len(day)):
					for j in range(i):
				
						weights1 = np.ones_like(day[i][1])/float(len(day[i][1]))
						weights2 = np.ones_like(day[j][1])/float(len(day[j][1]))
				
						n1, bins1, patches1 = plt.hist(day[i][1], weights=weights1,bins=200)#,50, normed=1, facecolor='green', alpha=0.75)
						n2, bins2, patches2 = plt.hist(day[j][1], weights=weights2,bins=200)#,50, normed=1, facecolor='green', alpha=0.75)
				
						#print stats.ks_2samp(n1,n2)[1],namestr(day[i], globals()),namestr(day[j], globals()),namestr(day, globals())
						if stats.ks_2samp(n1,n2)[1]>= 0.0: # NO POOLING!!
							Pooled_Experiments.append(day[i][1].mean()/2.0)
							Pooled_Experiments.append(day[j][1].mean()/2.0)
						
						#	print i,j,day[i][1].mean(),day[j][1].mean(),namestr(day[i], globals()),namestr(day[j], globals())
		
				if Pooled_Experiments != []:
					if namestr(day, globals())[0] == 'day':
						print >> fout_D80_Not_Pooled,namestr(day, globals())[1],np.average(Pooled_Experiments), np.std(Pooled_Experiments)
					else:
						print >> fout_D80_Not_Pooled,namestr(day, globals())[0],np.average(Pooled_Experiments), np.std(Pooled_Experiments)
				
			elif len(day)<=1:
				for i in range(len(day)):
					Pooled_Experiments.append(day[i][1].mean())
				
				if Pooled_Experiments != []:
					if namestr(day, globals())[0] == 'day':
						print >> fout_D80_Not_Pooled,namestr(day, globals())[1],np.average(Pooled_Experiments), np.std(Pooled_Experiments)
					else:
						print >> fout_D80_Not_Pooled,namestr(day, globals())[0],np.average(Pooled_Experiments), np.std(Pooled_Experiments)

if Pooled_1 == 1:
	################################################################################################
	######################## CALCULATION POOLING TOGETHER SIMILAR EXPERIMENTS ##############
	################################################################################################
	if Day_35 == 1:
		for day in D35:
			Pooled_Experiments = []
			Not_Pooled_Experiments = []
			if len(day)>1:
				
				for i in range(len(day)):
					for j in range(i):
				
						weights1 = np.ones_like(day[i][1])/float(len(day[i][1]))
						weights2 = np.ones_like(day[j][1])/float(len(day[j][1]))
				
						n1, bins1, patches1 = plt.hist(day[i][1], weights=weights1,bins=200)#,50, normed=1, facecolor='green', alpha=0.75)
						n2, bins2, patches2 = plt.hist(day[j][1], weights=weights2,bins=200)#,50, normed=1, facecolor='green', alpha=0.75)
				
						#print stats.ks_2samp(n1,n2)[1],namestr(day[i], globals()),namestr(day[j], globals()),namestr(day, globals())
						if stats.ks_2samp(n1,n2)[1]>= 0.01: # NO POOLING!!
							Pooled_Experiments.append(day[i][1].mean()/2.0)
							Pooled_Experiments.append(day[j][1].mean()/2.0)
						elif stats.ks_2samp(n1,n2)[1]< 0.01:
							Not_Pooled_Experiments.append(day[i][1].mean())
						#	print i,j,day[i][1].mean(),day[j][1].mean(),namestr(day[i], globals()),namestr(day[j], globals())
		
				if Pooled_Experiments != []:
					if namestr(day, globals())[0] == 'day':
						print >> fout_D35_Pooled,namestr(day, globals())[1],np.average(Pooled_Experiments), np.std(Pooled_Experiments)
					else:
						print >> fout_D35_Pooled,namestr(day, globals())[0],np.average(Pooled_Experiments), np.std(Pooled_Experiments)
				elif Not_Pooled_Experiments != []:
					if namestr(day, globals())[0] == 'day':
						print >> fout_D35_Pooled,namestr(day, globals())[1],np.average(Not_Pooled_Experiments), np.std(Not_Pooled_Experiments)
					else:
						print >> fout_D35_Pooled,namestr(day, globals())[0],np.average(Not_Pooled_Experiments), np.std(Not_Pooled_Experiments)
			elif len(day)<=1:
				for i in range(len(day)):
					Pooled_Experiments.append(day[i][1].mean())
				
				if Pooled_Experiments != []:
					if namestr(day, globals())[0] == 'day':
						print >> fout_D35_Pooled,namestr(day, globals())[1],np.average(Pooled_Experiments), np.std(Pooled_Experiments)
					else:
						print >> fout_D35_Pooled,namestr(day, globals())[0],np.average(Pooled_Experiments), np.std(Pooled_Experiments)
	if Day_50 == 1:
		for day in D50:
			Pooled_Experiments = []
			Not_Pooled_Experiments = []
			if len(day)>1:
				
				for i in range(len(day)):
					for j in range(i):
				
						weights1 = np.ones_like(day[i][1])/float(len(day[i][1]))
						weights2 = np.ones_like(day[j][1])/float(len(day[j][1]))
				
						n1, bins1, patches1 = plt.hist(day[i][1], weights=weights1,bins=200)#,50, normed=1, facecolor='green', alpha=0.75)
						n2, bins2, patches2 = plt.hist(day[j][1], weights=weights2,bins=200)#,50, normed=1, facecolor='green', alpha=0.75)
				
						#print stats.ks_2samp(n1,n2)[1],namestr(day[i], globals()),namestr(day[j], globals()),namestr(day, globals())
						if stats.ks_2samp(n1,n2)[1]>= 0.01: # NO POOLING!!
							Pooled_Experiments.append(day[i][1].mean()/2.0)
							Pooled_Experiments.append(day[j][1].mean()/2.0)
						elif stats.ks_2samp(n1,n2)[1]< 0.01:
							Not_Pooled_Experiments.append(day[i][1].mean())
						
						#print i,j,day[i][1].mean(),day[j][1].mean(),namestr(day[i], globals()),namestr(day[j], globals())
		
				if Pooled_Experiments != []:
					if namestr(day, globals())[0] == 'day':
						print >> fout_D50_Pooled,namestr(day, globals())[1],np.average(Pooled_Experiments), np.std(Pooled_Experiments)
					else:
						print >> fout_D50_Pooled,namestr(day, globals())[0],np.average(Pooled_Experiments), np.std(Pooled_Experiments)
				elif Not_Pooled_Experiments != []:
					if namestr(day, globals())[0] == 'day':
						print >> fout_D50_Pooled,namestr(day, globals())[1],np.average(Not_Pooled_Experiments), np.std(Not_Pooled_Experiments)
					else:
						print >> fout_D50_Pooled,namestr(day, globals())[0],np.average(Not_Pooled_Experiments), np.std(Not_Pooled_Experiments)
				
			elif len(day)<=1:
				for i in range(len(day)):
					Pooled_Experiments.append(day[i][1].mean())
				
				if Pooled_Experiments != []:
					if namestr(day, globals())[0] == 'day':
						print >> fout_D50_Pooled,namestr(day, globals())[1],np.average(Pooled_Experiments), np.std(Pooled_Experiments)
					else:
						print >> fout_D50_Pooled,namestr(day, globals())[0],np.average(Pooled_Experiments), np.std(Pooled_Experiments)
	
	if Day_80 == 1:
		for day in D80:
			Pooled_Experiments = []
			Not_Pooled_Experiments = []
			if len(day)>1:
				
				for i in range(len(day)):
					for j in range(i):
				
						weights1 = np.ones_like(day[i][1])/float(len(day[i][1]))
						weights2 = np.ones_like(day[j][1])/float(len(day[j][1]))
				
						n1, bins1, patches1 = plt.hist(day[i][1], weights=weights1,bins=200)#,50, normed=1, facecolor='green', alpha=0.75)
						n2, bins2, patches2 = plt.hist(day[j][1], weights=weights2,bins=200)#,50, normed=1, facecolor='green', alpha=0.75)
				
						#print stats.ks_2samp(n1,n2)[1],namestr(day[i], globals()),namestr(day[j], globals()),namestr(day, globals())
						if stats.ks_2samp(n1,n2)[1]>= 0.01: # NO POOLING!!
							Pooled_Experiments.append(day[i][1].mean()/2.0)
							Pooled_Experiments.append(day[j][1].mean()/2.0)
						elif stats.ks_2samp(n1,n2)[1]< 0.01:
							Not_Pooled_Experiments.append(day[i][1].mean())
						#	print i,j,day[i][1].mean(),day[j][1].mean(),namestr(day[i], globals()),namestr(day[j], globals())
		
				if Pooled_Experiments != []:
					if namestr(day, globals())[0] == 'day':
						print >> fout_D80_Pooled,namestr(day, globals())[1],np.average(Pooled_Experiments), np.std(Pooled_Experiments)
					else:
						print >> fout_D80_Pooled,namestr(day, globals())[0],np.average(Pooled_Experiments), np.std(Pooled_Experiments)
				elif Not_Pooled_Experiments != []:
					if namestr(day, globals())[0] == 'day':
						print >> fout_D80_Pooled,namestr(day, globals())[1],np.average(Not_Pooled_Experiments), np.std(Not_Pooled_Experiments)
					else:
						print >> fout_D80_Pooled,namestr(day, globals())[0],np.average(Not_Pooled_Experiments), np.std(Not_Pooled_Experiments)
			elif len(day)<=1:
				for i in range(len(day)):
					Pooled_Experiments.append(day[i][1].mean())
				
				if Pooled_Experiments != []:
					if namestr(day, globals())[0] == 'day':
						print >> fout_D80_Pooled,namestr(day, globals())[1],np.average(Pooled_Experiments), np.std(Pooled_Experiments)
					else:
						print >> fout_D80_Pooled,namestr(day, globals())[0],np.average(Pooled_Experiments), np.std(Pooled_Experiments)



if Pooled == 1:
	################################################################################################
	######################## CALCULATION POOLING TOGETHER SIMILAR EXPERIMENTS ##############
	################################################################################################
	# Distance measurement (in terms of KS test) and selection of the experiment with largest number of neurons
	
	
	
	if Day_35 == 1:
		All_Pooled =[]
		All_Non_Pooled = []
		for day in D35:
			Pooled_Experiments = []
			Non_Pooled_Experiments = []
		
			if len(day)>1:
				for i in range(len(day)):
					for j in range(i):
						
						#fig = plt.figure(1)         # create a figure instance
						#ax = fig.add_subplot(111)   # and axes
						weights1 = np.ones_like(day[i][1])/float(len(day[i][1]))
						weights2 = np.ones_like(day[j][1])/float(len(day[j][1]))
				
						n1, bins1, patches1 = plt.hist(day[i][1], weights=weights1,bins=200)#,50, normed=1, facecolor='green', alpha=0.75)
						n2, bins2, patches2 = plt.hist(day[j][1], weights=weights2,bins=200)#,50, normed=1, facecolor='green', alpha=0.75)
						#plt.xlim([0,.2])
						#plt.ylim([0,1.0])
						
						                
						#if stats.ks_2samp(day[i][1],day[j][1])[1]>= 0.01:
						if stats.ks_2samp(n1,n2)[1]>= 0.01: # POOLING SIMILAR EXPERIMENTS!!
							if namestr(day[i], globals()) and namestr(day[j], globals()) not in Pooled_Experiments:
								#fig.savefig('Histogram_%s_%s_P.png' %(namestr(day[i], globals()),namestr(day[j], globals()))) 
								#plt.close()
								Pooled_Experiments.append(namestr(day[i], globals()))
								Pooled_Experiments.append(namestr(day[j], globals()))
							
						#if stats.ks_2samp(day[i][1],day[j][1])[1]< 0.01:
						elif stats.ks_2samp(n1,n2)[1] < 0.01:
							if namestr(day[i], globals()) and namestr(day[j], globals()) not in Pooled_Experiments and namestr(day[i], globals()) and namestr(day[j], globals()) not in Non_Pooled_Experiments:
								#fig.savefig('Histogram_%s_%s_NP.png' %(namestr(day[i], globals()),namestr(day[j], globals()))) 
								#plt.close()
								Non_Pooled_Experiments.append(namestr(day[i], globals()))
								Non_Pooled_Experiments.append(namestr(day[j], globals()))
							
		
				for elem in Non_Pooled_Experiments:
					if elem in Pooled_Experiments:
						Non_Pooled_Experiments.remove(elem)
				#print Non_Pooled_Experiments
				#print Pooled_Experiments
			All_Pooled.append(Pooled_Experiments)
			All_Non_Pooled.append(Non_Pooled_Experiments)
	 
		#print All_Non_Pooled
		for i in range(len(All_Non_Pooled)):
			if All_Non_Pooled[i] == []:
				for j in range(len(All_Pooled[i])):
				
					All_Non_Pooled[i].append(All_Pooled[i][j])

		Selected_Experiments = []
		for i in range(len(All_Non_Pooled)):
			for j in range(len(All_Non_Pooled)):
				Index=[]
				for k in range(len(All_Non_Pooled[i])):
				
					for m in range(len(All_Non_Pooled[j])):
						if j!=i:
						
							#print (i,k),(j,m)	
							weights1 = np.ones_like(D35[i][k][1])/float(len(D35[i][k][1]))
							weights2 = np.ones_like(D35[j][m][1])/float(len(D35[j][m][1]))
						
							n1, bins1, patches1 = plt.hist(D35[i][k][1], weights=weights1,bins=200)#,50, normed=1, facecolor='green', alpha=0.75)
							n2, bins2, patches2 = plt.hist(D35[j][m][1], weights=weights2,bins=200)#,50, normed=1, facecolor='green', alpha=0.75)
						
							#print All_Non_Pooled[i][k],All_Non_Pooled[j][m],stats.ks_2samp(n1,n2)[1]
							Index.append((All_Non_Pooled[i][k],All_Non_Pooled[j][m],stats.ks_2samp(n1,n2)[1]))
							#print All_Non_Pooled[i][k],All_Non_Pooled[j][m],stats.ks_2samp(n1,n2)[1]
				if Index != []:
						Selected_Experiments.append(min(Index,key=lambda x: x[2]))
						#print min(Index,key=lambda x: x[2])	
	
		List_Experiments = []
		for i in range(len(Selected_Experiments)):
			List_Experiments.append(str(Selected_Experiments[i][0]))
		
		#print set(List_Experiments)
		x = Counter(List_Experiments)
		#print x.most_common()[0]
	
		


		for day in D35:
			for i in range(len(x.most_common()[:4])):
				Pooled_Experiments = []
				for exp in day:
					#print str(str(x.most_common()[i][0]).strip('[]')).strip("''")==str(namestr(exp, globals())[1])
					if str(str(x.most_common()[i][0]).strip('[]')).strip("''") == namestr(exp, globals())[0] or str(str(x.most_common()[i][0]).strip('[]')).strip("''") == namestr(exp, globals())[1]:
					
							
						
						 
						
						#print >> fout_HIST_D35_Pooled,bins1,n1
						
						if namestr(exp, globals())[0] == 'exp':
							print >> fout_D35_Pooled,namestr(exp, globals())[1],np.average(exp[1]), np.std(exp[1])
							
						else:
							print >> fout_D35_Pooled,namestr(exp, globals())[0],np.average(exp[1]), np.std(exp[1])
							
						plt.close()
	if Day_50 == 1:
		All_Pooled =[]
		All_Non_Pooled = []
		for day in D50:
			Pooled_Experiments = []
			Non_Pooled_Experiments = []
		
			if len(day)>1:
				for i in range(len(day)):
					for j in range(i):
				
						weights1 = np.ones_like(day[i][1])/float(len(day[i][1]))
						weights2 = np.ones_like(day[j][1])/float(len(day[j][1]))
				
						n1, bins1, patches1 = plt.hist(day[i][1], weights=weights1,bins=200)#,50, normed=1, facecolor='green', alpha=0.75)
						n2, bins2, patches2 = plt.hist(day[j][1], weights=weights2,bins=200)#,50, normed=1, facecolor='green', alpha=0.75)
				
						#if stats.ks_2samp(day[i][1],day[j][1])[1]>= 0.01:
						if stats.ks_2samp(n1,n2)[1]>= 0.01: # POOLING SIMILAR EXPERIMENTS!!
							if namestr(day[i], globals()) and namestr(day[j], globals()) not in Pooled_Experiments:
							
								Pooled_Experiments.append(namestr(day[i], globals()))
								Pooled_Experiments.append(namestr(day[j], globals()))
							
						#elif stats.ks_2samp(day[i][1],day[j][1])[1] < 0.01:
						elif stats.ks_2samp(n1,n2)[1] < 0.01:
							if namestr(day[i], globals()) and namestr(day[j], globals()) not in Pooled_Experiments and namestr(day[i], globals()) and namestr(day[j], globals()) not in Non_Pooled_Experiments:
						
								Non_Pooled_Experiments.append(namestr(day[i], globals()))
								Non_Pooled_Experiments.append(namestr(day[j], globals()))
							
		
				for elem in Non_Pooled_Experiments:
					if elem in Pooled_Experiments:
						Non_Pooled_Experiments.remove(elem)
				#print Non_Pooled_Experiments
				#print Pooled_Experiments
			All_Pooled.append(Pooled_Experiments)
			All_Non_Pooled.append(Non_Pooled_Experiments)
	 
		for i in range(len(All_Non_Pooled)):
			if All_Non_Pooled[i] == []:
				for j in range(len(All_Pooled[i])):
				
					All_Non_Pooled[i].append(All_Pooled[i][j])
		

		#print All_Non_Pooled[0][2]
				
		Selected_Experiments = []
		for i in range(len(All_Non_Pooled)):
			for j in range(len(All_Non_Pooled)):
				Index=[]
				for k in range(len(All_Non_Pooled[i])):
				
					for m in range(len(All_Non_Pooled[j])):
						if j!=i:
						
							#print (i,k),(j,m)	
							weights1 = np.ones_like(D50[i][k][1])/float(len(D50[i][k][1]))
							weights2 = np.ones_like(D50[j][m][1])/float(len(D50[j][m][1]))
						
							n1, bins1, patches1 = plt.hist(D50[i][k][1], weights=weights1,bins=200)#,50, normed=1, facecolor='green', alpha=0.75)
							n2, bins2, patches2 = plt.hist(D50[j][m][1], weights=weights2,bins=200)#,50, normed=1, facecolor='green', alpha=0.75)
						
						
							Index.append((All_Non_Pooled[i][k],All_Non_Pooled[j][m],stats.ks_2samp(n1,n2)[1]))
							
							#print All_Non_Pooled[i][k],All_Non_Pooled[j][m],stats.ks_2samp(n1,n2)[1]
				if Index != []:
						Selected_Experiments.append(min(Index,key=lambda x: x[2]))
						#print min(Index,key=lambda x: x[2])	
	
		List_Experiments = []
		for i in range(len(Selected_Experiments)):
			List_Experiments.append(str(Selected_Experiments[i][0]))
		
		#print set(List_Experiments)
		x = Counter(List_Experiments)
		#print x.most_common()[0]
	
		


		for day in D50:
			for i in range(len(x.most_common()[:4])):
				Pooled_Experiments = []
				for exp in day:
					#print str(str(x.most_common()[i][0]).strip('[]')).strip("''")==str(namestr(exp, globals())[1])
					if str(str(x.most_common()[i][0]).strip('[]')).strip("''") == namestr(exp, globals())[0] or str(str(x.most_common()[i][0]).strip('[]')).strip("''") == namestr(exp, globals())[1]:
					
							

						#weights1 = np.ones_like(exp[1])/float(len(exp[1]))
						
						#n1, bins1, patches1 = plt.hist(exp[1], weights=weights1,bins=200)#,50, normed=1, facecolor='green', alpha=0.75)
						
						#print >> fout_HIST_D50_Pooled,bins1,n1
	

						if namestr(exp, globals())[0] == 'exp':
							print >> fout_D50_Pooled,namestr(exp, globals())[1],np.average(exp[1]), np.std(exp[1])
						else:
							print >> fout_D50_Pooled,namestr(exp, globals())[0],np.average(exp[1]), np.std(exp[1])
	
	if Day_80 == 1:
		All_Pooled =[]
		All_Non_Pooled = []
		for day in D80:
			Pooled_Experiments = []
			Non_Pooled_Experiments = []
		
			if len(day)>1:
				for i in range(len(day)):
					for j in range(i):
				
						weights1 = np.ones_like(day[i][1])/float(len(day[i][1]))
						weights2 = np.ones_like(day[j][1])/float(len(day[j][1]))
				
						n1, bins1, patches1 = plt.hist(day[i][1], weights=weights1,bins=200)#,80, normed=1, facecolor='green', alpha=0.75)
						n2, bins2, patches2 = plt.hist(day[j][1], weights=weights2,bins=200)#,80, normed=1, facecolor='green', alpha=0.75)
				
						#if stats.ks_2samp(day[i][1],day[j][1])[1]>= 0.01:
						if stats.ks_2samp(n1,n2)[1]>= 0.01: # POOLING SIMILAR EXPERIMENTS!!
							if namestr(day[i], globals()) and namestr(day[j], globals()) not in Pooled_Experiments:
							
								Pooled_Experiments.append(namestr(day[i], globals()))
								Pooled_Experiments.append(namestr(day[j], globals()))
							
						#elif stats.ks_2samp(day[i][1],day[j][1])[1] < 0.01:
						elif stats.ks_2samp(n1,n2)[1] < 0.01:
							if namestr(day[i], globals()) and namestr(day[j], globals()) not in Pooled_Experiments and namestr(day[i], globals()) and namestr(day[j], globals()) not in Non_Pooled_Experiments:
						
								Non_Pooled_Experiments.append(namestr(day[i], globals()))
								Non_Pooled_Experiments.append(namestr(day[j], globals()))
							
		
				for elem in Non_Pooled_Experiments:
					if elem in Pooled_Experiments:
						Non_Pooled_Experiments.remove(elem)
				#print Non_Pooled_Experiments
				#print Pooled_Experiments
			All_Pooled.append(Pooled_Experiments)
			All_Non_Pooled.append(Non_Pooled_Experiments)
	 
		for i in range(len(All_Non_Pooled)):
			if All_Non_Pooled[i] == []:
				for j in range(len(All_Pooled[i])):
				
					All_Non_Pooled[i].append(All_Pooled[i][j])
		

		#print len(All_Non_Pooled[0])
				
		Selected_Experiments = []
		for i in range(len(All_Non_Pooled)):
			for j in range(len(All_Non_Pooled)):
				Index=[]
				for k in range(len(All_Non_Pooled[i])):
				
					for m in range(len(All_Non_Pooled[j])):
						if j!=i:
						
							#print (i,k),(j,m)	
							weights1 = np.ones_like(D80[i][k][1])/float(len(D80[i][k][1]))
							weights2 = np.ones_like(D80[j][m][1])/float(len(D80[j][m][1]))
						
							n1, bins1, patches1 = plt.hist(D80[i][k][1], weights=weights1,bins=200)#,80, normed=1, facecolor='green', alpha=0.75)
							n2, bins2, patches2 = plt.hist(D80[j][m][1], weights=weights2,bins=200)#,80, normed=1, facecolor='green', alpha=0.75)
						
						
							Index.append((All_Non_Pooled[i][k],All_Non_Pooled[j][m],stats.ks_2samp(n1,n2)[1]))
							
							#print All_Non_Pooled[i][k],All_Non_Pooled[j][m],stats.ks_2samp(n1,n2)[1]
				if Index != []:
						Selected_Experiments.append(min(Index,key=lambda x: x[2]))
						#print min(Index,key=lambda x: x[2])	
	
		List_Experiments = []
		for i in range(len(Selected_Experiments)):
			List_Experiments.append(str(Selected_Experiments[i][0]))
		
		#print set(List_Experiments)
		x = Counter(List_Experiments)
		#print x.most_common()[0]
	
		


		for day in D80:
			for i in range(len(x.most_common()[:4])):
				Pooled_Experiments = []
				for exp in day:
					#print str(str(x.most_common()[i][0]).strip('[]')).strip("''")==str(namestr(exp, globals())[1])
					if str(str(x.most_common()[i][0]).strip('[]')).strip("''") == namestr(exp, globals())[0] or str(str(x.most_common()[i][0]).strip('[]')).strip("''") == namestr(exp, globals())[1]:
					
							

						#weights1 = np.ones_like(exp[1])/float(len(exp[1]))
						
						#n1, bins1, patches1 = plt.hist(exp[1], weights=weights1,bins=200)#,50, normed=1, facecolor='green', alpha=0.75)
						
						#print >> fout_HIST_D80_Pooled,bins1,n1

	

						if namestr(exp, globals())[0] == 'exp':
							print >> fout_D80_Pooled,namestr(exp, globals())[1],np.average(exp[1]), np.std(exp[1])
						else:
							print >> fout_D80_Pooled,namestr(exp, globals())[0],np.average(exp[1]), np.std(exp[1])
