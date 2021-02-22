#!/usr/bin/python
# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Agg')
from matplotlib.lines import Line2D   
import networkx as nx
import matplotlib.pyplot as plt
import pylab
import os
import numpy as np
import matplotlib.cm as cm
import matplotlib.animation as manimation
from matplotlib import collections  as mc
import random

SEED = int(os.getenv('SEED'))

random.seed(SEED)

Dir_Input = os.getenv('Dir_Input')

Write_Dopa_Neurons = 1
Prune_Input_Connections = 0
Prune_Output_Connections = 0
Excitatory = 0
Inhibitory = 0

Exc_Inh_List = np.loadtxt('%s/Exc_Inh_List.dat' %Dir_Input,skiprows=1,unpack=True)

Percentage_retracted =float(os.getenv('Percentage_retracted'))
Percentage_retraction=float(os.getenv('Percentage_retraction'))
Prop_exc_inh=float(os.getenv('Prop_exc_inh'))

Connections = np.loadtxt('%s/cons10_original.txt' %Dir_Input, skiprows=9, unpack=False)

#Retracted=random.sample(Exc_Inh_List[:int((1-Prop_exc_inh)*(len(Exc_Inh_List)))],int(Percentage_retracted*len(Exc_Inh_List)/100))
Retracted=random.sample(Exc_Inh_List[:int((Prop_exc_inh)*(len(Exc_Inh_List)))],int(Percentage_retracted*len(Exc_Inh_List[:int((Prop_exc_inh)*(len(Exc_Inh_List)))])/100))
if Write_Dopa_Neurons == 1:
	fout_retracted = open('%s/Retracted_neurons.txt' %Dir_Input,'w')
	for i in Retracted:
		print>>fout_retracted,int(i)
Inhibitory_neurons = Exc_Inh_List[int((Prop_exc_inh)*(len(Exc_Inh_List))):]
nNeurons = int((len(Exc_Inh_List)))
Kex1 = np.zeros([nNeurons,nNeurons])

if Excitatory == 1:
	fout = open('%s/cons10_new.txt' %Dir_Input,'w')
	if Prune_Input_Connections == 1:
		for j,i in Connections:
			i = int(i)
			j = int(j)
			#print i,j
			if i in Retracted:

				Kex1[i,j] = 1
	
			else:
				Kex1[i,j] = 2
	if Prune_Output_Connections == 1:
		for i,j in Connections:
			i = int(i)
			j = int(j)
			#print i,j
			if i in Retracted:

				Kex1[j,i] = 1
	
			else:
				Kex1[j,i] = 2
	print >> fout,'%-----------------------------------------------------------------'
	print >> fout,'% Neuron11 '
	print >> fout,'% Copyright (c) 2010 Javier G. Orlandi <javiergorlandi@gmail.com> '
	print >> fout,'%-----------------------------------------------------------------'
	print >> fout,'% Connection List'
	print >> fout,'% Format: Input | Output'
	print >> fout,'%-----------------------------------------------------------------'
	print >> fout,'% Seed: 0'
	print >> fout,'%-----------------------------------------------------------------'
	Total_initial_Length = len(zip(*np.where( Kex1==1)))

	if Percentage_retraction != 0:
		
		for i in range( int( Total_initial_Length)-int( Total_initial_Length-Total_initial_Length*Percentage_retraction/100 ) ):
			
			Kex1[random.choice(zip(*np.where( Kex1==1)))] = 0

	for j in range(len(Kex1)):
		for i in range(len(Kex1)):
			if Kex1[i,j] !=0:
				print >> fout,j,i

elif Inhibitory == 1:
	fout = open('%s/cons10_less_inh.txt' %Dir_Input,'w')
	for i,j in Connections:
		i = int(i)
		j = int(j)
		#print i,j
		if i in Inhibitory_neurons:

			Kex1[j,i] = 1
	
		else:
			Kex1[j,i] = 2

	print >> fout,'%-----------------------------------------------------------------'
	print >> fout,'% Neuron11 '
	print >> fout,'% Copyright (c) 2010 Javier G. Orlandi <javiergorlandi@gmail.com> '
	print >> fout,'%-----------------------------------------------------------------'
	print >> fout,'% Connection List'
	print >> fout,'% Format: Input | Output'
	print >> fout,'%-----------------------------------------------------------------'
	print >> fout,'% Seed: 0'
	print >> fout,'%-----------------------------------------------------------------'
	Total_initial_Length = len(zip(*np.where( Kex1==1)))

	if Percentage_retraction != 0:
		for i in range( int( Total_initial_Length)-int( Total_initial_Length-Total_initial_Length*Percentage_retraction/100 ) ):
		
			Kex1[random.choice(zip(*np.where( Kex1==1)))] = 0

	for j in range(len(Kex1)):
		for i in range(len(Kex1)):
			if Kex1[i,j] !=0:
				print >> fout,j,i