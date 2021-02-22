#!/usr/bin/python

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
FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='Movie Test', artist='Matplotlib',
                comment='Movie support!')
#writer = FFMpegWriter(fps=10, metadata=metadata)
writer = FFMpegWriter(fps=10, codec = "libx264",bitrate=-1,extra_args=['-pix_fmt', 'yuv420p'],metadata=metadata)

from matplotlib import style
style.use('dark_background')

def chunker(seq, size):
	return (tuple(seq[pos:pos+size]) for pos in xrange(0, len(seq), size))

#File_Dat1="./Input/Kex1.dat"
#File_Dat2="./Input/Kex2.dat"
#File_Dat3="./Input/Kin.dat"

#nNodes=int(os.getenv('Nodes'))

#fout1=open(File_Dat1,'w')
#fout2=open(File_Dat2,'w')
#fout3=open(File_Dat3,'w')

#print >>fout1,nNodes


Create_movie = 1

#Constructed_network = np.loadtxt('Input/cons10.txt', skiprows=9,unpack=True)
Positions = np.loadtxt('Input/map10.txt', skiprows=7,unpack=True)
Axons = open('Input/axons10.txt','r')

Neurons = np.loadtxt('./Input/Exc_Inh_List.dat',skiprows=1,unpack=True)
Dopa = np.loadtxt('./Input/Dopaminergic_list.dat',skiprows=1,unpack=True)
Exc = Neurons[:int(0.8*len(Neurons))]
Inh = Neurons[int(0.8*len(Neurons)):]


Calcium = np.loadtxt('./Analysis/Fluorescence/fluorescence_000_0m_25hz__all.dat',delimiter=' ',unpack=True)
#Calcium = np.loadtxt('Output/Kth0000/Data_Calcium_0000.dat',unpack=True)
Total_time = len(Calcium[0])
pos_exc = {}
pos_inh = {}
pos_dop = {}
for i in Positions[0]:
	if int(i) in Exc:
		pos_exc[int(i)] = (Positions[1][int(i)],Positions[2][int(i)])
	elif int(i) in Inh:
		pos_inh[int(i)] = (Positions[1][int(i)],Positions[2][int(i)])
	elif int(i) in Dopa:
		pos_dop[int(i)] = (Positions[1][int(i)],Positions[2][int(i)])

nNeurons = len(Positions[0])
X=nx.Graph()
X.add_nodes_from(range(nNeurons))




Center = [0,0]
Radius_assay = 5


if Create_movie==1:
	
	an = np.linspace(0, 2*np.pi, 100)
	Arena=[Center[0]+Radius_assay*np.cos(an),Center[1]+Radius_assay*np.sin(an)]
	
	fig = plt.figure()
	ax = fig.add_subplot(111,axisbg='black')
	fig.tight_layout()
	#nx.draw_networkx_nodes(X,pos,node_size=10,nodelist=range(0,nNeurons),node_color='gray')
	#ax.plot(Arena[0],Arena[1],lw=3,color='darkmagenta')
	#ax.set_aspect(1)
	#plt.axis('off')
	
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
	#Calcium.close()
	
	#nx.draw_networkx_nodes(X,pos,node_size=10,nodelist=range(0,nNeurons),node_color='k')
	
	
	norm = matplotlib.colors.Normalize(vmin=0.0, vmax=0.6)
	cmap1 = cm.hot #cm.Reds
	cmap2 = cm.hot #cm.Blues
	cmap3 = cm.hot #cm.Greens
	m1 = cm.ScalarMappable(norm=norm, cmap=cmap1)
	m2 = cm.ScalarMappable(norm=norm, cmap=cmap2)
	m3 = cm.ScalarMappable(norm=norm, cmap=cmap3)
	
	
	
	print 'Starting movie!'
	
	with writer.saving(fig, "./Movies/fluorescence_000_0m_25hz.mpg", 150):
	
		for i in range(Total_time):
			
			if i%1==0:
				
				ax.set_aspect(1)
				plt.axis('off')
				for j in Dopa:
					j=int(j)
					a = nx.draw_networkx_nodes(X,pos_dop,node_size=10,nodelist=[j],node_color=m1.to_rgba(Calcium[j][i]),linewidths=0)
					lc_1 = mc.LineCollection(line_collection[j], colors=m1.to_rgba(Calcium[j][i]), linewidths=0.1)
					ax.add_collection(lc_1)
					ax.collections[0].set_edgecolor("k") 
				
				for k in Inh:
					k = int(k)
					b = nx.draw_networkx_nodes(X,pos_inh,node_size=10,nodelist=[k],node_color=m2.to_rgba(Calcium[k][i]),linewidths=0)
					lc_2 = mc.LineCollection(line_collection[k], colors=m2.to_rgba(Calcium[k][i]), linewidths=0.1)
					ax.add_collection(lc_2)
					ax.collections[0].set_edgecolor("k") 
				
				for l in Exc:
					l = int(l)
					c = nx.draw_networkx_nodes(X,pos_exc,node_size=10,nodelist=[l],node_color=m3.to_rgba(Calcium[l][i]),linewidths=0)
					lc_3 = mc.LineCollection(line_collection[l], colors=m3.to_rgba(Calcium[l][i]), linewidths=0.1)
					ax.add_collection(lc_3)
					ax.collections[0].set_edgecolor("k") 
				writer.grab_frame()
				#if a!=[]:
				#	a.remove()
				b.remove()
				c.remove()
				ax.cla()
				
	print '...finished!'
	
	
