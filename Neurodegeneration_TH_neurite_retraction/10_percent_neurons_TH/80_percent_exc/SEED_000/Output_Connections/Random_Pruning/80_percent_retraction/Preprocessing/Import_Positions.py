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


Create_Picture = 1

Data_Dir = os.getenv('Dir_Input') 

Output_Dir=os.getenv('Output_Map_Dir') 

fout = open('%s/map.txt' %Output_Dir, 'w')
fout2 = open('%s/Neurons_list.dat' %Data_Dir,'w')
fout3 = open('%s/Dopaminergic_List.dat' %Data_Dir,'w')
'''
print >> fout, '%-----------------------------------------------------------------'
print >> fout, '% Neuron11' 
print >> fout, '% Copyright (c) 2010 Javier G. Orlandi <javiergorlandi@gmail.com> '
print >> fout, '%-----------------------------------------------------------------'
print >> fout, '% Generated Positional Map'
print >> fout, '% Format: Neuron # | X | Y'
print >> fout, '%-----------------------------------------------------------------'
'''
Real_pos_neurons = np.loadtxt('%s/Regions.dat' %Data_Dir, unpack=True, delimiter=',')
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

	print len(set(np.where(d<13.9)[1]))
	nNeuronsTH = len(Real_pos_TH[0])


nNeurons = len(Real_pos_neurons[0])

if Real_pos_TH !=[]:
	for i in list(set(np.where(d<13.9)[1])):
		positions_neurons_x_TH.append(Real_pos_neurons[1][i])
		positions_neurons_y_TH.append(Real_pos_neurons[2][i])

		print >> fout3,int(Real_pos_neurons[0][i]-1)

#New_list= [int(x) for x in Real_pos_neurons[0] if x not in list(set(np.where(d<13.9)[1]))]
New_list= [int(x) for x in Real_pos_neurons[0]]
for k in New_list:
	print >> fout2, k


for i in range(len(Real_pos_neurons[0])):
	positions_neurons.append( ( Real_pos_neurons[1][i]-int( max(Real_pos_neurons[1])/2),Real_pos_neurons[2][i]-int( max(Real_pos_neurons[2])/2)) )
	#positions_neurons.append( ( Real_pos_neurons[1][i],Real_pos_neurons[2][i]) )
	positions_neurons_x.append(Real_pos_neurons[1][i]-int( max(Real_pos_neurons[1])/2))
	positions_neurons_y.append(-Real_pos_neurons[2][i]+int( max(Real_pos_neurons[2])/2))
for i,(j,k) in enumerate(zip(positions_neurons_x, positions_neurons_y)):
	print >> fout, i,j/100,k/100
if Create_Picture ==1:

	fig = plt.figure()
	ax = fig.add_subplot(111,axisbg='black')
	ax.set_aspect(1)
	plt.axis('off')	
	X=nx.Graph()
	X.add_nodes_from(range(nNeurons))
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
for i in range(len(Real_pos_TH[0])):
	positions_neurons.append( ( Real_pos_TH[1][i]-int( max(Real_pos_TH[1])/2),-Real_pos_TH[2][i]+int( max(Real_pos_TH[2])/2)) )
	#positions_neurons.append( ( Real_pos_TH[1][i],Real_pos_TH[2][i]) )
	positions_neurons_x.append(Real_pos_neurons[1][i]-int( max(Real_pos_neurons[1])/2))
	positions_neurons_y.append(-Real_pos_neurons[2][i]+int( max(Real_pos_neurons[2])/2))


for i,(j,k) in enumerate(zip(positions_neurons_x, positions_neurons_y)):
	print >> fout, i,j/100,k/100
'''

#for i, (j,k) in enumerate(positions_neurons):
#	print >>fout,i,j/100,k/100
#	print i,j/100,k/100
#plt.scatter(positions_neurons_x,positions_neurons_y)
#plt.scatter(Real_pos_neurons[1],-Real_pos_neurons[2])
#plt.show()