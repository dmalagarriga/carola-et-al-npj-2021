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
from scipy.stats import norm
from scipy.spatial import distance
import random

random.seed(0)

Create_Picture = 1

Dir_Input = os.getenv('P_Dir') 
Percentage_retracted =float(os.getenv('Percentage_retracted'))
Prop_exc_inh=float(os.getenv('Prop_exc_inh'))

Exc_Inh_List = np.loadtxt('%s/Exc_Inh_List.dat' %Dir_Input,skiprows=1,unpack=True)
Retracted=random.sample(Exc_Inh_List[:int((Prop_exc_inh)*(len(Exc_Inh_List)))],int(Percentage_retracted*len(Exc_Inh_List)/100))

Excitatory_list = [ i for i in Exc_Inh_List[:int((Prop_exc_inh)*(len(Exc_Inh_List)))] if i not in Retracted ] 
Inhibitory_list = Exc_Inh_List[int((Prop_exc_inh)*(len(Exc_Inh_List))):len(Exc_Inh_List)]


Real_pos_neurons = np.loadtxt('%s/Regions.dat' %Dir_Input, unpack=True, delimiter=',')
Real_pos_TH = []
#Real_pos_neurons = np.loadtxt('%s/Regions2nd.txt' %Data_Dir, unpack=True, delimiter=',')
#Real_pos_TH = np.loadtxt('%s/Regions2ndTH.txt' %Data_Dir, unpack=True, delimiter=',')
positions_neurons = []
positions_neurons_x = []
positions_neurons_y = []
if Real_pos_TH!=[]:
	positions_neurons_x_TH = []
	positions_neurons_y_TH = []
count=0
indices = []




Pos_N = zip(Real_pos_neurons[1],Real_pos_neurons[2])
if Real_pos_TH !=[]:
	Pos_TH = zip(Real_pos_TH[1],Real_pos_TH[2])

	d = distance.cdist(Pos_TH,Pos_N)
	#Len_Coin_Points.append(len(set(np.where(d<0.01)[1])))


	nNeuronsTH = len(Real_pos_TH[0])


nNeurons = len(Real_pos_neurons[0])

if Real_pos_TH !=[]:
	for i in list(set(np.where(d<13.9)[1])):
		positions_neurons_x_TH.append(Real_pos_neurons[1][i])
		positions_neurons_y_TH.append(Real_pos_neurons[2][i])



New_list= [int(x) for x in Real_pos_neurons[0]]




for i in range(len(Real_pos_neurons[0])):
	positions_neurons.append( ( Real_pos_neurons[1][i]-int( max(Real_pos_neurons[1])/2),Real_pos_neurons[2][i]-int( max(Real_pos_neurons[2])/2)) )
	#positions_neurons.append( ( Real_pos_neurons[1][i],Real_pos_neurons[2][i]) )
	positions_neurons_x.append(Real_pos_neurons[1][i]-int( max(Real_pos_neurons[1])/2))
	positions_neurons_y.append(-Real_pos_neurons[2][i]+int( max(Real_pos_neurons[2])/2))

if Create_Picture ==1:

	fig = plt.figure()
	ax = fig.add_subplot(111,axisbg='black')
	ax.set_aspect(1)
	plt.axis('off')	
	X=nx.Graph()
	X.add_nodes_from(range(nNeurons))
	
	pos_exc = {}
	pos_inh = {}
	pos_retracted = {}
	
	for i in Excitatory_list:
		i = int(i)
		pos_exc[i] = (positions_neurons_x[i],positions_neurons_y[i])
	for j in Inhibitory_list:
		j = int(j)
		pos_inh[j] = (positions_neurons_x[j],positions_neurons_y[j])
	for k in Retracted:
		k = int(k)
		pos_retracted[k] = (positions_neurons_x[k],positions_neurons_y[k])
	
	
	for i in Excitatory_list:
		a = nx.draw_networkx_nodes(X,pos_exc,node_size=100,nodelist=[i],alpha=1.0,node_color='blue')

	for j in Inhibitory_list:	
		b = nx.draw_networkx_nodes(X,pos_inh,node_size=100,nodelist=[j],alpha=1.0,node_color='green',node_shape='d')
	
	for k in Retracted:
		c = nx.draw_networkx_nodes(X,pos_retracted,node_size=150,nodelist=[k],alpha=1.0,node_color='red',node_shape='h')
	
	ax.set_xlim([min(positions_neurons_x)-15,max(positions_neurons_x)+15])
	ax.set_ylim([min(positions_neurons_y)-15,max(positions_neurons_y)+15])
	
	plt.savefig('%s/Network.png' %Dir_Input)
	plt.close()
	'''
	pos_neurons = {}
	for i in range(nNeurons):
		#pos_neurons[i] = (Real_pos_neurons[1][i],Real_pos_neurons[2][i])
		pos_neurons[i] = (positions_neurons_x[i],positions_neurons_y[i])
	for i in range(nNeurons):
		a = nx.draw_networkx_nodes(X,pos_neurons,node_size=50,nodelist=[i],node_color='darkgray',linewidths=0)
	if Real_pos_TH !=[]:
		pos_TH_1 = {}
		for i in range(nNeuronsTH):
			pos_TH_1[i] = (Real_pos_TH[1][i],Real_pos_TH[2][i])
		pos_TH_2 = {}
		for i in range(nNeuronsTH):
			pos_TH_2[i] = (positions_neurons_x_TH[i],positions_neurons_y_TH[i])
	
		

		for j in range(nNeuronsTH):
			b = nx.draw_networkx_nodes(X,pos_TH_1,node_size=30,nodelist=[j],node_color='red',linewidths=0)
	
		for k in range(nNeuronsTH):
			c = nx.draw_networkx_nodes(X,pos_TH_2,node_size=15,nodelist=[k],node_color='green',linewidths=0)
	#ax.set_xlim([min(Real_pos_neurons[1])-5,max(Real_pos_neurons[1])+5])
	#ax.set_ylim([min(Real_pos_neurons[2])-5,max(Real_pos_neurons[2])+5])
	ax.set_xlim([min(positions_neurons_x)-5,max(positions_neurons_x)+5])
	ax.set_ylim([min(positions_neurons_y)-5,max(positions_neurons_y)+5])
	plt.savefig('./Preprocessing/Network.png')
	'''
