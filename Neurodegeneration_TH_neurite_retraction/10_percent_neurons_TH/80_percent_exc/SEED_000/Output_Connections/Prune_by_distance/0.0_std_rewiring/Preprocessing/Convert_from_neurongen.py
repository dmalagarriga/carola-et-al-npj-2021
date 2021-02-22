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
#FFMpegWriter = manimation.writers['ffmpeg']
#metadata = dict(title='Movie Test', artist='Matplotlib',
#                comment='Movie support!')
#writer = FFMpegWriter(fps=10, metadata=metadata)
#writer = FFMpegWriter(fps=100, codec = "libx264",bitrate=-1,extra_args=['-pix_fmt', 'yuv420p'],metadata=metadata)

from matplotlib import style
style.use('dark_background')

random.seed(0)


def chunker(seq, size):
	return (tuple(seq[pos:pos+size]) for pos in xrange(0, len(seq), size))

File_Dat1="./Input/Kex1.dat"
File_Dat2="./Input/Kex2.dat"
File_Dat3="./Input/Kin.dat"
File_Dat4="./Input/Exc_Inh_List.dat"
File_Dat5="./Input/Dopaminergic_List.dat"
File_Dat6="./Input/Delay_Matrix.dat"
File_Dat7="./Input/Connections_all.dat"
File_Dat8="./Input/Number_of_Connections_all.dat"

#nNeurons=int(os.getenv('Nodes'))

fout1=open(File_Dat1,'w')
fout2=open(File_Dat2,'w')
fout3=open(File_Dat3,'w')
fout4=open(File_Dat4,'w')
fout5=open(File_Dat5,'w')
fout6=open(File_Dat6,'w')
fout7=open(File_Dat7,'w')
fout8=open(File_Dat8,'w')

Positions = np.loadtxt('Input/map10.txt', skiprows=7,unpack=False)
Axons = open('Input/axons10.txt','r')
Connections = np.loadtxt('Input/cons10.txt', skiprows=9, unpack=False)

Prop_exc_inh=float(os.getenv('Prop_exc_inh'))

nDopaminergic = 0
nNeurons = int(len(Positions))
Kex1 = np.zeros([nNeurons,nNeurons])
Kex2 = np.zeros([nNeurons,nNeurons])
Kin = np.zeros([nNeurons,nNeurons])
Delay_Matrix = np.zeros([nNeurons,nNeurons]) 

Excitatory_list = random.sample(range(nNeurons-nDopaminergic),int((1-Prop_exc_inh)*(nNeurons-nDopaminergic)))
Inhibitory_list = [i for i in range(nNeurons-nDopaminergic) if i not in Excitatory_list]
Dopaminergic_list=list(np.arange(nNeurons-nDopaminergic,nNeurons))

dt = 0.05 # ms
velocity = 0.5 #mm/ms 
l_segment = 0.01 # mm 
'''
line_collection=[[] for _ in range(nNeurons)]
for i,line in enumerate(Axons.readlines()[7:]): # 7 is the first line with data
	
	x = range(3,len(line.split()))
	running = True

	while running:
		for idx,elem in enumerate( list(chunker(x,2)) ):
			thiselem = elem
			if idx == len(list(chunker(x,2)))-1:
				nextelem = thiselem
			else:
				nextelem = list(chunker(x,2))[(idx + 1) % len(list(chunker(x,2)))]
			
			
			line_collection[i].append([(float(line.split()[ thiselem[0] ]),float(line.split()[ thiselem[1] ])),(float(line.split()[ nextelem[0] ]),float(line.split()[ nextelem[1] ]))])
			
			
		running = False
		
	
	
Axons.close()

for i,j in Connections:
	i=int(i)
	j=int(j)
	for k in range(len(line_collection[int(i)])):
		if np.sqrt( (line_collection[int(i)][k][0][0]-Positions[int(j)][1])**2 + (line_collection[int(i)][k][0][1]-Positions[int(j)][2])**2 ) < 0.15:
			Delay_Matrix[j,i] = int(k*l_segment/(velocity*dt))
			#print i,j
			
'''
for i,j in Connections:
	i=int(i)
	j=int(j)
	Dist = np.sqrt( (Positions[int(j)][1]-Positions[int(i)][1])**2 + ( Positions[int(j)][2] - Positions[int(i)][2] )**2 )
	
	Delay_Matrix[j,i] = int(float(Dist)/velocity)






print >> fout4,nNeurons
for i in Excitatory_list:
	print >> fout4,i
for i in Inhibitory_list:
	print >> fout4,i
print >> fout5,nDopaminergic
for i in Dopaminergic_list:
	print >> fout5,i



counter=0
for i,j in Connections:
	i = int(i)
	j = int(j)
	
	if i in Excitatory_list[:int(len(Excitatory_list)/2)]:
		counter+=1
		
		Kex1[j,i] = 1
		
	elif i in Excitatory_list[int(len(Excitatory_list)/2):]:
		Kex1[j,i] = 2
		counter+=1
	elif int(i)> nNeurons-nDopaminergic:
		Kex2[j,i] = 3
		counter+=1
	elif i in Inhibitory_list:
		Kin[j,i] = -1
		counter+=1
print>>fout7, counter
#print np.count_nonzero(Kex1)+np.count_nonzero(Kin)

for i in Excitatory_list:
	counter=0
	for j in range(nNeurons): # Half of the excitatory are AMPA
		if Kex1[i,j]==1:
			
			
			counter+=1
			print>>fout7,j, "1" # AMPA 
			print>>fout6,int(Delay_Matrix[i,j])
	
		elif Kex1[i,j]==2:
			
			counter+=1
	
			print>>fout7, j, "2" # NMDA 
			print>>fout6,int(Delay_Matrix[i,j])
		elif Kin[i,j]==-1:
			counter+=1
			print>>fout7, j, "0" # GABA 
			print>>fout6,int(Delay_Matrix[i,j])
		elif Kex2[i,j]==3:
			counter+=1
			print>>fout7, j, "2" # As if DOPA was exc NMDA 
			print>>fout6,int(Delay_Matrix[i,j])
	print>>fout8,counter

for i in Inhibitory_list:
	counter=0
	for j in range(nNeurons):
		if Kex1[i,j]==1:
			
			
			counter+=1
			print>>fout7,j, "1" # AMPA 
			print>>fout6,int(Delay_Matrix[i,j])
	
		elif Kex1[i,j]==2:
			
			counter+=1
	
			print>>fout7, j, "2" # NMDA 
			print>>fout6,int(Delay_Matrix[i,j])
		elif Kin[i,j]==-1:
			counter+=1
			print>>fout7, j, "0" # GABA 
			print>>fout6,int(Delay_Matrix[i,j])
		elif Kex2[i,j]==3:
			counter+=1
			print>>fout7, j, "2" # As if DOPA was exc NMDA 
			print>>fout6,int(Delay_Matrix[i,j])
	
	print>>fout8,counter
for i in Dopaminergic_list:
	counter=0
	for j in range(nNeurons):
		if Kex1[i,j]==1:
			
			
			counter+=1
			print>>fout7,j, "1" # AMPA 
			print>>fout6,int(Delay_Matrix[i,j])
	
		elif Kex1[i,j]==2:
			
			counter+=1
	
			print>>fout7, j, "2" # NMDA 
			print>>fout6,int(Delay_Matrix[i,j])
		elif Kin[i,j]==-1:
			counter+=1
			print>>fout7, j, "0" # GABA 
			print>>fout6,int(Delay_Matrix[i,j])
		elif Kex2[i,j]==3:
			counter+=1
			print>>fout7, j, "2" # As if DOPA was exc NMDA 
			print>>fout6,int(Delay_Matrix[i,j])
	print>>fout8,counter

print >> fout1,nNeurons



np.savetxt(fout1,Kex1,fmt='%03.2f')
np.savetxt(fout2,Kex2,fmt='%03.2f')
np.savetxt(fout3,Kin,fmt='%03.2f')
#np.savetxt(fout6,Delay_Matrix,fmt='%03i')



