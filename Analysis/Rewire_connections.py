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

random.seed(0)

Dir_Input = os.getenv('Dir_Input')
Prop_exc_inh=float(os.getenv('Prop_exc_inh'))

Exc_Inh_List = np.loadtxt('%s/Exc_Inh_List.dat' %Dir_Input,skiprows=1,unpack=True)
fout = open('%s/cons10_new.txt' %Dir_Input,'w')
Percentage_retracted =float(os.getenv('Percentage_retracted'))
Number_std=float(os.getenv('Number_std')) # Number of std out of which rewire connections (0, 1 or 2)

Connections = np.loadtxt('%s/cons10_original.txt' %Dir_Input, skiprows=9, unpack=False)

Positions = np.loadtxt('%s/map10.txt' %Dir_Input, skiprows=7,unpack=False)

Retracted=random.sample(Exc_Inh_List[:int((Prop_exc_inh)*(len(Exc_Inh_List)))],int(Percentage_retracted*len(Exc_Inh_List)/100))

nNeurons = int((len(Exc_Inh_List)))
Kex1 = np.zeros([nNeurons,nNeurons])

ind_old = 0
length_each_neuron = []
New_Conns = []

for i,j in Connections:
	i = int(i)
	j = int(j)

	Dist = np.sqrt( (Positions[int(j)][1]-Positions[int(i)][1])**2 + ( Positions[int(j)][2] - Positions[int(i)][2] )**2 )
	if ind_old in Retracted:
		if i == ind_old:
			length_each_neuron.append((j,Dist))
	
		elif i!=ind_old:
			
			#if ind_old!=int(Connections[-1][0]):
			#if ind_old in Retracted:
			# Rewire connections that have a length above average+1std
			indices = [list(zip(*length_each_neuron)[1]).index(k) for k in zip(*length_each_neuron)[1] if k>np.average(zip(*length_each_neuron)[1])+Number_std*np.std(zip(*length_each_neuron)[1])]
			Retracted_Connections = [zip(*length_each_neuron)[0][l] for l in indices  ]

			D = [np.sqrt( (Positions[int(m)][1]-Positions[int(i)][1])**2 + ( Positions[int(m)][2] - Positions[int(i)][2] )**2 ) for m in range(nNeurons)]
			Less = np.where(np.array(D)<np.average(zip(*length_each_neuron)[1]))
			# Comprovations
			#print ind_old,Retracted_Connections
			#print ind_old,zip(*length_each_neuron)[0]
			#print ind_old, [o for o in zip(*length_each_neuron)[0] if o not in Retracted_Connections]
			#print ind_old,[random.choice(zip(*Less))[0] for n in range(len(Retracted_Connections))]
			New_wiring = [o for o in zip(*length_each_neuron)[0] if o not in Retracted_Connections]
			New_wiring.extend([random.choice(zip(*Less))[0] for n in range(len(Retracted_Connections))])
			#print ind_old,New_wiring
			for p in New_wiring:
				New_Conns.append((ind_old,p))
				Kex1[p,ind_old] = 1
		
			length_each_neuron = []
			length_each_neuron.append((j,Dist))
			
	else:
		if i == ind_old:
			length_each_neuron.append((j,Dist))
		elif i!=ind_old:
			
			
			for s in zip(*length_each_neuron)[0]:
				New_Conns.append((ind_old,s))
				Kex1[s,ind_old] = 1
	
			length_each_neuron = []
			length_each_neuron.append((j,Dist))
			
	ind_old = i
	if i==int(Connections[-1][0]):
		
		New_Conns.append((i,j))
		
print >> fout,'%-----------------------------------------------------------------'
print >> fout,'% Neuron11 '
print >> fout,'% Copyright (c) 2010 Javier G. Orlandi <javiergorlandi@gmail.com> '
print >> fout,'%-----------------------------------------------------------------'
print >> fout,'% Connection List'
print >> fout,'% Format: Input | Output'
print >> fout,'%-----------------------------------------------------------------'
print >> fout,'% Seed: 0'
print >> fout,'%-----------------------------------------------------------------'	


for q,r in sorted(New_Conns,key=lambda x:x[0]):
	print >>fout,q,r


