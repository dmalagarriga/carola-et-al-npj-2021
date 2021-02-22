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
writer = FFMpegWriter(fps=40, codec = "libx264",bitrate=-1,extra_args=['-pix_fmt', 'yuv420p'],metadata=metadata)

from matplotlib import style
style.use('dark_background')
#style.use('fivethirtyeight')

random.seed(0)


def chunker(seq, size):
	return (tuple(seq[pos:pos+size]) for pos in xrange(0, len(seq), size))


Create_Picture = 0
Dir_Input = os.getenv('Dir_Input')

Positions = np.loadtxt('%s/map10.txt' %Dir_Input, skiprows=7,unpack=True)
Axons = open('%s/axons10.txt' %Dir_Input,'r')
Exc_Inh_list = np.loadtxt('%s/Exc_Inh_List.dat' %Dir_Input,skip_rows=1,unpack=True)

Percentage_retracted = 10.0
Percentage_retraction = 25.0

Retracted=random.sample(Exc_inh_list[:int(0.8*(len(Exc_Inh_list)))],int(Percentage_retracted*len(Exc_Inh_list)/100))


nNeurons = len(Positions[0])

pos = {}
for i in range(nNeurons):
	pos[i] = (Positions[1][i],Positions[2][i])



X=nx.Graph()
X.add_nodes_from(range(nNeurons))





if Create_Picture ==1:
	fig = plt.figure()
	
	ax = fig.add_subplot(111,axisbg='black')
	ax.set_aspect(1)
	plt.axis('off')		
	fig.tight_layout()
	line_collection=[[] for _ in range(nNeurons)]
	for i,line in enumerate(Axons.readlines()[7:]): # 7 is the first line with data
		
		x = range(3,len(line.split()))
		running = True

		while running:
			for idx,elem in enumerate( list(chunker(x,2)) ):
				thiselem = elem
				if idx == (len(list(chunker(x,2)))-1)*Percentage_retraction/100#len(list(chunker(x,2)))-1:
					nextelem = thiselem
				else:
					nextelem = list(chunker(x,2))[(idx + 1) % len(list(chunker(x,2)))]
				
				
				line_collection[i].append([(float(line.split()[ thiselem[0] ]),float(line.split()[ thiselem[1] ])),(float(line.split()[ nextelem[0] ]),float(line.split()[ nextelem[1] ]))])
				
				
			running = False
	for j in range(nNeurons):
		a = nx.draw_networkx_nodes(X,pos,node_size=30,nodelist=[j],node_color='darkgray',linewidths=0)
		lc = mc.LineCollection(line_collection[j], colors='gray', linewidths=0.8)
		ax.add_collection(lc)
		ax.collections[0].set_edgecolor("k") 
	ax.set_xlim([min(Positions[1])-.1,max(Positions[1])+.1])
	ax.set_ylim([min(Positions[2])-.1,max(Positions[2])+.1])
	
	plt.savefig('Network.png',bbox_inches='tight')