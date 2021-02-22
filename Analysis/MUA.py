#!/usr/bin/python
# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Agg')
import io
import base64
import numpy as np
import commands
from itertools import tee, izip
import matplotlib.pylab as plt
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerPatch
import os
from peak_detection import *

def window(iterable, size):
    iters = tee(iterable, size)
    for i in xrange(1, size):
        for each in iters[i:]:
            next(each, None)
    return izip(*iters)

Fluo_Dir= os.getenv('Fluo_Dir')
Dir_Input = os.getenv('Dir_Input')
Percentage_retracted =float(os.getenv('Percentage_retracted'))
Percentage_retraction=float(os.getenv('Percentage_retraction'))

Calculation = 0
Figures = 1
Draw_peaks = 0

# Calculation of the multiunit activity #
#Raster_calc = np.loadtxt('../EXPD50SP11_SP12/SP11#1/+fact/E2/New_Project/sp11_d50_e2/data/Raster_fluorescence.dat',delimiter=',',unpack=True)
Raster_calc = np.loadtxt('%s/Raster_fluorescence.dat' %Fluo_Dir,delimiter=',',unpack=True)
Dopaminergic_neurons=np.loadtxt('%s/Retracted_neurons.txt' %Dir_Input,unpack=True)
Exc_Inh_list = np.loadtxt('%s/Exc_Inh_List.dat' %Dir_Input,unpack=True)


if Calculation == 1:
	Exc_neurons = [i for i in Exc_Inh_list[1:int(0.8*len(Exc_Inh_list))] if i not in Dopaminergic_neurons]
	Inh_neurons = Exc_Inh_list[int(0.8*len(Exc_Inh_list))+1:]
	Dopa_neurons = [i for i in Dopaminergic_neurons]
	
	fout_MUA = open('%s/Multi_Unit_Activity.dat' %Fluo_Dir,'w')
	fout_SUA = open('%s/Single_Unit_Activity.dat' %Fluo_Dir,'w')

	fout_SUA_exc = open('%s/Single_Unit_Activity_exc.dat' %Fluo_Dir,'w')
	fout_MUA_exc = open('%s/Multi_Unit_Activity_exc.dat' %Fluo_Dir,'w')

	fout_SUA_inh = open('%s/Single_Unit_Activity_inh.dat' %Fluo_Dir,'w')
	fout_MUA_inh = open('%s/Multi_Unit_Activity_inh.dat' %Fluo_Dir,'w')


	fout_SUA_dopa = open('%s/Single_Unit_Activity_dopa.dat' %Fluo_Dir,'w')
	fout_MUA_dopa = open('%s/Multi_Unit_Activity_dopa.dat' %Fluo_Dir,'w')

	Bursts = []
	Bursts_exc = []
	Bursts_inh = []
	Bursts_dopa = []
	MUA = [] # multi-unit activity
	SUA = [] # single unit activity
	MUA_exc = [] # multi-unit activity
	SUA_exc = [] # single unit activity
	MUA_inh = [] # multi-unit activity
	SUA_inh = [] # single unit activity
	MUA_dopa = [] # multi-unit activity
	SUA_dopa = [] # single unit activity

	Raster_exc = [i for i in Raster_calc[1] if Raster_calc[0][list(Raster_calc[1]).index(i)] in Exc_neurons]
	Raster_inh = [i for i in Raster_calc[1] if Raster_calc[0][list(Raster_calc[1]).index(i)] in Inh_neurons]
	Raster_dopa = [i for i in Raster_calc[1] if Raster_calc[0][list(Raster_calc[1]).index(i)] in Dopa_neurons]
	for i in range(int(max(Raster_calc[1]))):
		Bursts_exc.append(float(list(Raster_exc).count(i))/float(int(max(Raster_calc[0]))))
		Bursts_inh.append(float(list(Raster_inh).count(i))/float(int(max(Raster_calc[0]))))  
		Bursts_dopa.append(float(list(Raster_dopa).count(i))/float(int(max(Raster_calc[0]))))
		Bursts.append(float(list(Raster_calc[1]).count(i))/float(int(max(Raster_calc[0]))))  
	    
	# 20.0 is the number of frames that form 1second recording             
	for i,each in enumerate(window(Bursts,20)):
		MUA.append((i,np.sum(each)*100))
		print >> fout_MUA,i,float(np.sum(each))*100
		
	for i,each in enumerate(window(Bursts_exc,20)):
		MUA_exc.append((i,np.sum(each)*100))
		print >> fout_MUA_exc,i,float(np.sum(each))*100

	for i,each in enumerate(window(Bursts_inh,20)):
		MUA_inh.append((i,np.sum(each)*100))
		print >> fout_MUA_inh,i,float(np.sum(each))*100

	for i,each in enumerate(window(Bursts_dopa,20)):
		MUA_dopa.append((i,np.sum(each)*100))
		print >> fout_MUA_dopa,i,float(np.sum(each))*100
    	
	'''
	for i in range(int(max(Raster_calc[0]))):

		SUA_exc.append((i,float(list(Raster_calc[0]).count(i))/(float(int(max(Raster_calc[1])))/20.0)))
		print >> fout_SUA_exc,i,float(list(Raster_calc[0]).count(i))/(float(int(max(Raster_calc[1])))/20.0)
   
		SUA_inh.append((i,float(list(Raster_calc[0]).count(i))/(float(int(max(Raster_calc[1])))/20.0)))
		print >> fout_SUA_inh,i,float(list(Raster_calc[0]).count(i))/(float(int(max(Raster_calc[1])))/20.0)
	
	for i in range(int(max(Raster_calc[0]))):
		if i in Dopa_neurons:
			SUA_dopa.append((i,float(list(Raster_calc[0]).count(i))/(float(int(max(Raster_calc[1])))/20.0)))
			print >> fout_SUA_dopa,i,float(list(Raster_calc[0]).count(i))/(float(int(max(Raster_calc[1])))/20.0)
	# 20.0 is the number of frames that form 1second recording             
	'''
	
if Figures == 1:

	Exc_neurons = [i for i in Exc_Inh_list[1:int(0.8*len(Exc_Inh_list))] if i not in Dopaminergic_neurons]
	Inh_neurons = Exc_Inh_list[int(0.8*len(Exc_Inh_list))+1:]
	Dopa_neurons = [i for i in Dopaminergic_neurons]
	Raster_exc = [i for i in Raster_calc[1] if Raster_calc[0][list(Raster_calc[1]).index(i)] in Exc_neurons]
	Raster_inh = [i for i in Raster_calc[1] if Raster_calc[0][list(Raster_calc[1]).index(i)] in Inh_neurons]
	Raster_dopa = [i for i in Raster_calc[1] if Raster_calc[0][list(Raster_calc[1]).index(i)] in Dopa_neurons]

	Exc_neurons = [i for i in Exc_Inh_list[1:int(0.8*len(Exc_Inh_list))] if i not in Dopaminergic_neurons]
	#Exc_Inh_list[1:int(0.8*len(Exc_Inh_list))]
	Inh_neurons = Exc_Inh_list[int(0.8*len(Exc_Inh_list))+1:]
	Dopa_neurons = [i for i in Dopaminergic_neurons]

	nNeurons = len(Exc_neurons)+len(Inh_neurons)+len(Dopa_neurons)
	#nDopaminergic = len(Dopaminergic_list[1:])
	times_exc = []
	times_inh = []
	times_dop = []
	spikes_exc = []
	spikes_inh = []
	spikes_dop = []


	for i in range(len(Raster_calc[1])):
		if int(Raster_calc[0][i]) in Exc_neurons:
			times_exc.append(Raster_calc[1][i])
			spikes_exc.append(int(Raster_calc[0][i]))
		elif int(Raster_calc[0][i]) in Inh_neurons:
			times_inh.append(Raster_calc[1][i])
			spikes_inh.append(int(Raster_calc[0][i]))
		elif int(Raster_calc[0][i]) in Dopa_neurons:
			times_dop.append(Raster_calc[1][i])
			spikes_dop.append(int(Raster_calc[0][i]))
	MUA_excitation = np.genfromtxt('%s/Multi_Unit_Activity_exc.dat' %Fluo_Dir,unpack=True,dtype=np.float64)
	MUA_inhibition = np.genfromtxt('%s/Multi_Unit_Activity_inh.dat' %Fluo_Dir,unpack=True,dtype=np.float64)
	MUA_dopa = np.genfromtxt('%s/Multi_Unit_Activity_dopa.dat' %Fluo_Dir,unpack=True,dtype=np.float64)
	MUA = np.genfromtxt('%s/Multi_Unit_Activity.dat' %Fluo_Dir,unpack=True,dtype=np.float64)
	
	###### CALCULATION OF PEAKS ####
	timing = []
	maxima = []
	maxtab,mintab = peakdet(MUA[1],1.) #maxtab,mintab = peakdet(MUA[1],.2)
	for i,j in maxtab:
		i = float(i)
		timing.append(i)
		maxima.append(j)

	
	
	
	# Plotting
	nullfmt = plt.NullFormatter()         # no labels
	label_size = 20
	plt.rcParams['xtick.labelsize'] = label_size 
	plt.rcParams['ytick.labelsize'] = label_size 

	# definitions for the axes
	left, width = 0.1, 0.65
	bottom, height = 0.1, 0.65
	bottom_h = left_h = left + width + 0.02

	rect_scatter = [left, bottom, width, height]
	rect_histx = [left, bottom_h, width, 0.2]
	rect_histy = [left_h, bottom, 0.2, height]

	# start with a rectangular Figure
	plt.figure(1, figsize=(12, 12))
	axScatter = plt.axes(rect_scatter)
	axHistx = plt.axes(rect_histx)
	#axHisty = plt.axes(rect_histy)

	# no labels
	axHistx.xaxis.set_major_formatter(nullfmt)
	#axHisty.yaxis.set_major_formatter(nullfmt)

	axScatter.plot(np.array(times_exc)/(20*60),spikes_exc,'|',markeredgewidth=3,ms=10,c='b') 
	axScatter.set_xlim([0.1,2.0])
	axScatter.set_ylim([.0,300])

	axScatter.yaxis.set_ticks(np.arange(axScatter.get_ylim()[0], axScatter.get_ylim()[1]+1, axScatter.get_ylim()[1]))
	axScatter.xaxis.set_ticks(np.arange(axScatter.get_xlim()[0], axScatter.get_xlim()[1]+1, axScatter.get_xlim()[1]-axScatter.get_xlim()[0]))
 
	axHistx.plot(np.array(MUA_excitation)[0]/(20*60),np.array(MUA_excitation)[1],'-',c='b',linewidth=2)
	axHistx.fill_between(np.array(MUA_excitation)[0]/(20*60), 0, np.array(MUA_excitation)[1],color='b')
	#axHisty.plot(np.array(SP11_D80_SUA)[1],np.array(SP11_D80_SUA)[0],'-',c='k')

	axHistx.yaxis.set_ticks(np.arange(axHistx.get_ylim()[0], axHistx.get_ylim()[1]+.1, axHistx.get_ylim()[1]))

	axHistx.set_xlim(axScatter.get_xlim())
	#axHisty.set_ylim(axScatter.get_ylim())
	#axHisty.set_ylim([0,1])

	axScatter.set_xlabel("Time (min)",fontsize=30,labelpad=-10)
	axScatter.set_ylabel("Neuron",fontsize=30,labelpad=-30)

	axHistx.set_ylabel("% neurons firing/s",fontsize=15,labelpad=-10)

	plt.savefig('%s/Raster_SP11_%02i_percent_ramdom_retraction_excitation.eps' %(Fluo_Dir,Percentage_retracted),transparent=True)
	plt.close()

	nullfmt = plt.NullFormatter()         # no labels
	label_size = 20
	plt.rcParams['xtick.labelsize'] = label_size 
	plt.rcParams['ytick.labelsize'] = label_size 

	# definitions for the axes
	left, width = 0.1, 0.65
	bottom, height = 0.1, 0.65
	bottom_h = left_h = left + width + 0.02

	rect_scatter = [left, bottom, width, height]
	rect_histx = [left, bottom_h, width, 0.2]
	rect_histy = [left_h, bottom, 0.2, height]

	# start with a rectangular Figure
	plt.figure(1, figsize=(12, 12))
	axScatter = plt.axes(rect_scatter)
	axHistx = plt.axes(rect_histx)
	#axHisty = plt.axes(rect_histy)

	# no labels
	axHistx.xaxis.set_major_formatter(nullfmt)
	#axHisty.yaxis.set_major_formatter(nullfmt)

	axScatter.plot(np.array(times_inh)/(20*60),spikes_inh,'|',markeredgewidth=3,ms=10,c='g') 
	axScatter.set_xlim([0.1,2.0])
	axScatter.set_ylim([.0,300])

	axScatter.yaxis.set_ticks(np.arange(axScatter.get_ylim()[0], axScatter.get_ylim()[1]+1, axScatter.get_ylim()[1]))
	axScatter.xaxis.set_ticks(np.arange(axScatter.get_xlim()[0], axScatter.get_xlim()[1]+1, axScatter.get_xlim()[1]-axScatter.get_xlim()[0]))
 
	axHistx.plot(np.array(MUA_inhibition)[0]/(20*60),np.array(MUA_inhibition)[1],'-',c='g',linewidth=2)
	axHistx.fill_between(np.array(MUA_inhibition)[0]/(20*60), 0, np.array(MUA_inhibition)[1],color='g')
	#axHisty.plot(np.array(SP11_D80_SUA)[1],np.array(SP11_D80_SUA)[0],'-',c='k')

	axHistx.yaxis.set_ticks(np.arange(axHistx.get_ylim()[0], axHistx.get_ylim()[1]+.1, axHistx.get_ylim()[1]))

	axHistx.set_xlim(axScatter.get_xlim())
	#axHisty.set_ylim(axScatter.get_ylim())
	#axHisty.set_ylim([0,1])

	axScatter.set_xlabel("Time (min)",fontsize=30,labelpad=-10)
	axScatter.set_ylabel("Neuron",fontsize=30,labelpad=-30)

	axHistx.set_ylabel("% neurons firing/s",fontsize=15,labelpad=-10)

	plt.savefig('%s/Raster_SP11_%02i_percent_ramdom_retraction_inhibition.eps' %(Fluo_Dir,Percentage_retracted),transparent=True)
	plt.close()

	# Plotting
	nullfmt = plt.NullFormatter()         # no labels
	label_size = 20
	plt.rcParams['xtick.labelsize'] = label_size 
	plt.rcParams['ytick.labelsize'] = label_size 

	# definitions for the axes
	left, width = 0.1, 0.65
	bottom, height = 0.1, 0.65
	bottom_h = left_h = left + width + 0.02

	rect_scatter = [left, bottom, width, height]
	rect_histx = [left, bottom_h, width, 0.2]
	rect_histy = [left_h, bottom, 0.2, height]

	# start with a rectangular Figure
	plt.figure(1, figsize=(12, 12))
	axScatter = plt.axes(rect_scatter)
	axHistx = plt.axes(rect_histx)
	#axHisty = plt.axes(rect_histy)

	# no labels
	axHistx.xaxis.set_major_formatter(nullfmt)
	#axHisty.yaxis.set_major_formatter(nullfmt)

	axScatter.plot(np.array(times_dop)/(20*60),spikes_dop,'|',markeredgewidth=3,ms=10,c='r') 
	axScatter.set_xlim([0.1,2.0])
	axScatter.set_ylim([.0,300])

	axScatter.yaxis.set_ticks(np.arange(axScatter.get_ylim()[0], axScatter.get_ylim()[1]+1, axScatter.get_ylim()[1]))
	axScatter.xaxis.set_ticks(np.arange(axScatter.get_xlim()[0], axScatter.get_xlim()[1]+1, axScatter.get_xlim()[1]-axScatter.get_xlim()[0]))
 
	axHistx.plot(np.array(MUA_dopa)[0]/(20*60),np.array(MUA_dopa)[1],'-',c='r',linewidth=2)
	axHistx.fill_between(np.array(MUA_dopa)[0]/(20*60), 0, np.array(MUA_dopa)[1],color='r')
	#axHisty.plot(np.array(SP11_D80_SUA)[1],np.array(SP11_D80_SUA)[0],'-',c='k')

	axHistx.yaxis.set_ticks(np.arange(axHistx.get_ylim()[0], axHistx.get_ylim()[1]+.1, axHistx.get_ylim()[1]))

	axHistx.set_xlim(axScatter.get_xlim())
	#axHisty.set_ylim(axScatter.get_ylim())
	#axHisty.set_ylim([0,1])

	axScatter.set_xlabel("Time (min)",fontsize=30,labelpad=-10)
	axScatter.set_ylabel("Neuron",fontsize=30,labelpad=-30)

	axHistx.set_ylabel("% neurons firing/s",fontsize=15,labelpad=-10)

	plt.savefig('%s/Raster_SP11_%02i_percent_ramdom_retraction_dopaminergic.eps' %(Fluo_Dir,Percentage_retracted),transparent=True)
	plt.close()
	
	# Plotting
	nullfmt = plt.NullFormatter()         # no labels
	label_size = 20
	plt.rcParams['xtick.labelsize'] = label_size 
	plt.rcParams['ytick.labelsize'] = label_size 

	# definitions for the axes
	left, width = 0.1, 0.65
	bottom, height = 0.1, 0.65
	bottom_h = left_h = left + width + 0.02

	rect_scatter = [left, bottom, width, height]
	rect_histx = [left, bottom_h, width, 0.2]
	rect_histy = [left_h, bottom, 0.2, height]

	# start with a rectangular Figure
	plt.figure(1, figsize=(12, 12))
	axScatter = plt.axes(rect_scatter)
	axHistx = plt.axes(rect_histx)
	#axHisty = plt.axes(rect_histy)

	# no labels
	axHistx.xaxis.set_major_formatter(nullfmt)
	#axHisty.yaxis.set_major_formatter(nullfmt)

	axScatter.plot(Raster_calc[1]/(20*60),Raster_calc[0],'|',markeredgewidth=3,ms=10,c='k') 
	axScatter.set_xlim([0.0,2.0])
	axScatter.set_ylim([.0,300])

	axScatter.yaxis.set_ticks(np.arange(axScatter.get_ylim()[0], axScatter.get_ylim()[1]+1, axScatter.get_ylim()[1]))
	axScatter.xaxis.set_ticks(np.arange(axScatter.get_xlim()[0], axScatter.get_xlim()[1]+1, axScatter.get_xlim()[1]-axScatter.get_xlim()[0]))
 
	axHistx.plot(np.array(MUA)[0]/(20*60),np.array(MUA)[1],'-',c='k',linewidth=2)
	axHistx.fill_between(np.array(MUA)[0]/(20*60), 0, np.array(MUA)[1],color='k')
	
	if Draw_peaks = 1:
		markers_on_1 = list([np.where(maxima>np.average(MUA[1])+1.0*np.std(MUA[1]))][0][0])
		markers_on_2 = list([np.where(maxima<np.average(MUA[1])+1.0*np.std(MUA[1]))][0][0])
		axHistx.plot(np.array(timing)/(20*60),np.array(maxima),'*',c='r',markersize=10,markevery=markers_on_1)
		axHistx.plot(np.array(timing)/(20*60),np.array(maxima),'o',c='g',markersize=7,markevery=markers_on_2)
	
	#print list([np.where(maxima>np.average(maxima)+1*np.std(maxima))][0][0])
	axHistx.yaxis.set_ticks(np.arange(axHistx.get_ylim()[0], axHistx.get_ylim()[1]+.1, axHistx.get_ylim()[1]))

	axHistx.set_xlim(axScatter.get_xlim())
	#axHisty.set_ylim(axScatter.get_ylim())

	axScatter.set_xlabel("Time (min)",fontsize=30,labelpad=-10)
	axScatter.set_ylabel("Neuron",fontsize=30,labelpad=-30)

	axHistx.set_ylabel("% neurons firing/s",fontsize=15,labelpad=-10)

	plt.savefig('%s/Raster_SP11_D80_Random_pruning_%02i_percent_DOPA.eps' %(Fluo_Dir,Percentage_retraction),transparent=True)
	#axHisty.set_xlabel("Freq. (Burst(s))/neuron")

	#fig,ax = plt.subplots(2,1,sharex=True)
	#gs = gridspec.GridSpec(2, 1, height_ratios=[.5, 3])

	#fig.set_figheight(10)
	#fig.set_figwidth(10)
	#ax[0] = plt.subplot(gs[0])
	#ax[0].plot(np.array(zip(*MUA))[0]/(20*60),np.array(zip(*MUA))[1],'-')
	#ax[1] = plt.subplot(gs[1])
	#ax[0].set_ylabel('Bursts/s')
	#ax[1].set_xlabel('Time(min)')
	#ax[1].set_ylabel('Neuron')
	#ax[1].plot(Raster_calc[1]/(20*60),Raster_calc[0],'|',ms=5)
	#ax[0].set_xlim([0.0,15.0])
	#ax[1].set_xlim([0.0,15.0])
	#ax[1].set_ylim([.0,819.0])
	#ax[0].set_ylim([.0,7.0])
	#plt.subplots_adjust(hspace=.0)

