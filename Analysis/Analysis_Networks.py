#!/usr/bin/python
# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use('Agg')
from matplotlib.lines import Line2D   
import networkx as nx
import matplotlib.pyplot as plt
import pylab
import os
import glob
import numpy as np
import matplotlib.cm as cm
import matplotlib.animation as manimation
from matplotlib import collections  as mc
from matplotlib.patches import FancyArrowPatch, Circle

from matplotlib import style
#style.use('dark_background')
style.use('fivethirtyeight')

def chunker(seq, size):
	return (tuple(seq[pos:pos+size]) for pos in xrange(0, len(seq), size))

def check_symmetric(a, tol=1e-8):
    return np.allclose(a, a.T, atol=tol)

def draw_network(G,pos,ax,sg=None):

    for n in G:
        c=Circle(pos[n],radius=10.2,alpha=0.5)
        #ax.add_patch(c)
        G.node[n]['patch']=c
        x,y=pos[n]
    seen={}
    for (u,v,d) in G.edges(data=True):
        n1=G.node[u]['patch']
        n2=G.node[v]['patch']
        rad=0.1
        if (u,v) in seen:
            rad=seen.get((u,v))
            rad=(rad+np.sign(rad)*0.1)*-1
        alpha=0.2
        color='k'

        e = FancyArrowPatch(n1.center,n2.center,patchA=n1,patchB=n2,
                            arrowstyle='-|>',
                            connectionstyle='arc3,rad=%s'%rad,
                            mutation_scale=10.0,
                            lw=0.7,
                            alpha=alpha,
                            color=color)
        seen[(u,v)]=rad
        ax.add_patch(e)
    for n in G:
        c=Circle(pos[n],radius=10.2,alpha=0.5,color='darkgreen')
        ax.add_patch(c)
    return e

Create_Picture = 0 # Not simultaneously with Analysis=1, please
Create_Picture_Overlay = 0
Analysis = 1
Exp_Dir = os.getenv('P_Dir') # Where regions.dat is
Data_Dir = os.getenv('Data_Dir') # Where adjacency matrix is
Percentage = 50



for file in glob.glob('%s/*_W_Matrix.dat' %Data_Dir):
	Adj_Matrix = np.loadtxt(file,unpack=True)


# Only taking into account connections stronger than
# a 80% of the maximum strength (Adj_Matrix.max()-Adj_Matrix.max()*20/100)

if Analysis == 0:
	Adj_Matrix[Adj_Matrix<=(Adj_Matrix.max()-Adj_Matrix.max()*Percentage/100)] = 0.0
	#Adj_Matrix[Adj_Matrix<3.0] = 0.0
elif Analysis == 1:
	Adj_Matrix[Adj_Matrix<2.0] = 0.0
	#Adj_Matrix[Adj_Matrix<=(Adj_Matrix.max()-Adj_Matrix.max()*50/100)] = 0.0
Conns = nx.from_numpy_matrix(Adj_Matrix)



G=nx.DiGraph(Adj_Matrix)
if Create_Picture_Overlay == 1:
	im = plt.imread('%s/bf1.BMP' %Exp_Dir)
	fig = plt.figure()
	#fig.set_size_inches(18.5, 10.5)
	
	ax = fig.add_subplot(111,axisbg='white')
	plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
	ax.set_aspect(1)
	implot = plt.imshow(im,zorder=0)
	Positions = np.loadtxt('%s/Regions.dat' %Exp_Dir,delimiter=',',unpack=True)
	nNeurons = len(Positions[0])
	pos = {}
	for i in range(nNeurons):
		pos[i] = (Positions[1][i],-Positions[2][i])
	
	draw_network(Conns, pos,ax)
	ax.set_xlim([min(Positions[1])-10,max(Positions[1])+10])
	ax.set_ylim([-min(Positions[2])-10,-max(Positions[2])+10])
	ax.grid(False)
	plt.axis('off') 
	#fig.tight_layout(pad=0)
	plt.gca().invert_yaxis()
	
	plt.savefig('%s/Functional_Network_Overlay.png' %Exp_Dir)
	plt.close()
	
if Create_Picture ==1:
	
	Positions = np.loadtxt('%s/Regions.dat' %Exp_Dir,delimiter=',',unpack=True)
	nNeurons = len(Positions[0])
	pos = {}
	for i in range(nNeurons):
		pos[i] = (Positions[1][i],-Positions[2][i])
		
	fig = plt.figure()
	ax = fig.add_subplot(111,axisbg='white')
	ax.set_aspect(1)
	plt.axis('off')		
	#a = nx.draw_networkx_nodes(X,pos,node_size=30,node_color='darkgray',linewidths=0)
	#a = nx.draw_networkx_nodes(Conns,pos,node_size=30,node_color='darkgray',linewidths=0)
	#b = nx.draw_networkx_edges(Conns,pos,width=.01)
	draw_network(Conns, pos,ax)
	ax.set_xlim([min(Positions[1])-10,max(Positions[1])+10])
	#ax.set_ylim([min(Positions[2])-10,max(Positions[2])+10])
	ax.set_ylim([min(-Positions[2])-10,max(-Positions[2])+10])
	
	plt.savefig('%s/Functional_Network_%s.png' %(Exp_Dir,Percentage))
	plt.close()
	
if Analysis == 1:
	###### IN AND OUT DEGREES ##### 
	fout_deg = open('%s/In_and_Out_Degrees.dat' %Exp_Dir,'w')
	In_degrees = dict(G.in_degree())
	In_values = sorted(set(In_degrees.values()))
	In_hist = [In_degrees.values().count(x) for x in In_values]
	
	Out_degrees = dict(G.out_degree())
	Out_values = sorted(set(Out_degrees.values()))
	Out_hist = [Out_degrees.values().count(x) for x in Out_values]
	
	for (i,j,k,l) in zip(*(In_values,In_hist,Out_values,Out_hist)):
		print >>fout_deg,i,j,k,l
	
	##### AVERAGE CLUSTERING COEFFICIENT ####
	fout_clust = open('%s/Avg_clustering.dat' %Exp_Dir,'w')
	
	G_ud = G.to_undirected()
	ccs = nx.clustering(G_ud)
	avg_clust = sum(ccs.values())/len(ccs)
	print >> fout_clust, avg_clust
	
	fout_all_data = open('%s/Nodes_characteristics.dat' %Exp_Dir,'a')
	print >> fout_all_data, 'Node','Clustering','Betweenness_centrality','Closeness_centrality','Eigenvector_centrality'
	
	#### Clustering ####
	clust_coefficients = nx.clustering(G_ud)
		
	#### Betweenness centrality ####
	G_components = list(nx.connected_component_subgraphs(G_ud))
	
	G_mc = G_components[0]
	
	
	bet_cen = nx.betweenness_centrality(G_mc)
	
	#### Closeness centrality ####
	clo_cen = nx.closeness_centrality(G_mc)
	
	
	#### Eigenvector centrality ####
	eig_cen = nx.eigenvector_centrality(G_mc)
	
	# PRINT TO FILE 
	results = [(k,clust_coefficients[k],bet_cen[k],clo_cen[k],eig_cen[k]) for k in G_mc]
	for item in results:
		fout_all_data.write(','.join(map(str,item)))
		fout_all_data.write('\n')
	fout_all_data.close()
	