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

#def draw_network(G,pos,ax,color_nodes,color_edges,sg=None):
def draw_network(G,pos,ax,color_nodes,sg=None):
    for n in G:
        c=Circle(pos[n],color=color_nodes,radius=10.2,alpha=1.0)
        
        G.node[n]['patch']=c
       	ax.add_patch(c)
       	#nx.draw_networkx_nodes(G,pos)
    	x,y=pos[n]
		
	
	#return c
def draw_edges(G1,G,pos,ax,color_edges,sg=None):
	for n in G:
		if n<max(list(G1.nodes())):
			c=Circle(pos[n],radius=10.2,alpha=1.0)
		#ax.add_patch(c)
		G.node[n]['patch']=c
		#x,y=pos[n]
	seen={}
	for (u,v,d) in G.edges(data=True):
		n1=G.node[u]['patch']
		n2=G.node[v]['patch']
		rad=0.1
		if (u,v) in seen:
			rad=seen.get((u,v))
			rad=(rad+np.sign(rad)*0.1)*-1
		alpha=0.5
		color=color_edges

		e = FancyArrowPatch(n1.center,n2.center,patchA=n1,patchB=n2,
							arrowstyle='-|>',
							connectionstyle='arc3,rad=%s'%rad,
							mutation_scale=10.0,
							lw=1.,
							alpha=alpha,
							color=color_edges)
		seen[(u,v)]=rad
		ax.add_patch(e)
	return e
	
Entropy = 0
Create_Picture = 1
Create_Picture_Overlay = 0
Analysis = 0
Exp_Dir = os.getenv('P_Dir')
Data_Dir = os.getenv('Data_Dir')
Percentage=40.0

if Entropy == 1:
	for file in glob.glob('%s/GTE.dat' %Data_Dir):
		Adj_Matrix = np.loadtxt(file,unpack=True)
else:
	for file in glob.glob('%s/*_W_Matrix.dat' %Data_Dir):
		Adj_Matrix = np.loadtxt(file,unpack=True)


# Only taking into account connections stronger than
# a 80% of the maximum strength (Adj_Matrix.max()-Adj_Matrix.max()*20/100)

Adj_Matrix[Adj_Matrix<=(Adj_Matrix.max()-Adj_Matrix.max()*Percentage/100)] = 0.0
Conns = nx.from_numpy_matrix(Adj_Matrix)


#print np.shape(Adj_Matrix)

G=nx.DiGraph(Adj_Matrix)
if Create_Picture_Overlay == 1:
	im = plt.imread('%s/1.BMP' %Exp_Dir)
	fig = plt.figure()
	#fig.set_size_inches(18.5, 10.5)
	
	ax = fig.add_subplot(111,axisbg='white')
	plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
	ax.set_aspect(1)
	implot = plt.imshow(im,zorder=0,cmap='gray')
	
	Positions_TH = np.loadtxt('%s/ROI_TH.txt' %Exp_Dir,unpack=True)
	Positions_nonTH = np.loadtxt('%s/ROI_nonTH.txt' %Exp_Dir,unpack=True)
	
	nNeurons = len(Positions_TH[0])+len(Positions_nonTH[0])
	Positions = np.hstack([Positions_nonTH,Positions_TH])
	
	Conns_nonTH = nx.DiGraph() 
	
	Conns_nonTH.add_edges_from(Conns.edges(range(len(Positions_nonTH[0]))))	
	
	
	Conns_TH = nx.DiGraph() 
	
	Conns_TH.add_edges_from(Conns.edges(range(len(Positions_nonTH[0]),len(Positions[0]))))
	
	
	pos = {}
	for i in range(nNeurons):
		pos[i] = (Positions[1][i],Positions[2][i])
		
	#fig = plt.figure()
	#ax = fig.add_subplot(111,axisbg='white')
	#ax.set_aspect(1)
	#plt.axis('off')		
	
	#draw_edges(Conns_nonTH,pos, ax,'b')
	
	TH = nx.DiGraph()
	TH.add_nodes_from(list(range(len(Positions_nonTH[0]),len(Positions[0]))))
	NonTH = nx.DiGraph()
	NonTH.add_nodes_from(list(range(len(Positions_nonTH[0]))))
	
	
	draw_edges(Conns_nonTH, pos,ax,'b')
	draw_edges(Conns_TH, pos,ax,'r')
	draw_network(NonTH, pos,ax,'b')
	draw_network(TH, pos,ax,'r')
	
	#ax.set_xlim([min(Positions[1])-10,max(Positions[1])+10])
	#ax.set_ylim([min(Positions[2])-10,max(Positions[2])+10])
	
	ax.grid(False)
	plt.axis('off') 
	#fig.tight_layout(pad=0)
	#plt.gca().invert_yaxis()
	
	plt.savefig('%s/Functional_Network_Overlay.png' %Exp_Dir)
	plt.close()
	
if Create_Picture ==1:
	
	Positions_TH = np.loadtxt('%s/ROI_TH.txt' %Exp_Dir,unpack=True)
	Positions_nonTH = np.loadtxt('%s/ROI_nonTH.txt' %Exp_Dir,unpack=True)
	
	nNeurons = len(Positions_TH[0])+len(Positions_nonTH[0])
	Positions = np.hstack([Positions_nonTH,Positions_TH])
	
	Conns_nonTH = nx.DiGraph() 
	
	Conns_nonTH.add_edges_from(Conns.edges(range(len(Positions_nonTH[0]))))	
	
	
	Conns_TH = nx.DiGraph() 
	
	Conns_TH.add_edges_from(Conns.edges(range(len(Positions_nonTH[0]),len(Positions[0]))))
	
	
	pos = {}
	
	for i in range(nNeurons):
		
		x1 = Positions[1][i] - max(Positions[1])/2
		y1 = Positions[2][i] - max(Positions[2])/2
		
		x2 = y1
		y2 = -x1
		
		x3 = x2 + max(Positions[1])/2
		y3 = y2 + max(Positions[2])/2
		#pos[i] = (Positions[1][i],Positions[2][i])
		pos[i] = (x3,y3)
	
	fig = plt.figure()
	ax = fig.add_subplot(111,axisbg='white')
	ax.set_aspect(1)
	plt.axis('off')		
	
	
	
	TH = nx.DiGraph()
	TH.add_nodes_from(list(range(len(Positions_nonTH[0]),len(Positions[0]))))
	NonTH = nx.DiGraph()
	NonTH.add_nodes_from(list(range(len(Positions_nonTH[0]))))
	
	draw_network(TH, pos,ax,'r')
	draw_network(NonTH, pos,ax,'b')
	draw_edges(NonTH,Conns_nonTH, pos,ax,'b')
	draw_edges(TH,Conns_TH, pos,ax,'r')
	
	
	ax.set_xlim([min(Positions[1])-150,max(Positions[1])+150])
	ax.set_ylim([min(Positions[2])-100,max(Positions[2])+100])
	plt.savefig('%s/Functional_Network_Thresh_%s.png' %(Exp_Dir,Percentage))
	plt.close()
	
	
if Analysis == 1:
	NonTH = 1
	TH = 0
	
	Positions_TH = np.loadtxt('%s/ROI_TH.txt' %Exp_Dir,unpack=True)
	Positions_nonTH = np.loadtxt('%s/ROI_nonTH.txt' %Exp_Dir,unpack=True)
	
	nNeurons = len(Positions_TH[0])+len(Positions_nonTH[0])
	Positions = np.hstack([Positions_nonTH,Positions_TH])
	
	Conns_nonTH = nx.DiGraph() 
	
	Conns_nonTH.add_edges_from(Conns.edges(range(len(Positions_nonTH[0]))))	
	
	
	Conns_TH = nx.DiGraph() 
	
	Conns_TH.add_edges_from(Conns.edges(range(len(Positions_nonTH[0]),len(Positions[0]))))
	if NonTH == 1:
		###############################################
		################## NON TH #####################
		###############################################
	
		###### IN AND OUT DEGREES NonTH ##### 
		fout_deg = open('%s/In_and_Out_Degrees_nonTH.dat' %Exp_Dir,'w')
		In_degrees = Conns_nonTH.in_degree()
		In_values = sorted(set(In_degrees.values()))
		In_hist = [In_degrees.values().count(x) for x in In_values]
	
		Out_degrees = Conns_nonTH.out_degree()
		Out_values = sorted(set(Out_degrees.values()))
		Out_hist = [Out_degrees.values().count(x) for x in Out_values]
	
		for (i,j,k,l) in zip(*(In_values,In_hist,Out_values,Out_hist)):
			print >>fout_deg,i,j,k,l
	
		##### AVERAGE CLUSTERING COEFFICIENT ####
		fout_clust = open('%s/Avg_clustering_nonTH.dat' %Exp_Dir,'w')
	
		Conns_nonTH_ud = Conns_nonTH.to_undirected()
		ccs = nx.clustering(Conns_nonTH_ud)
		avg_clust = sum(ccs.values())/len(ccs)
		print >> fout_clust, avg_clust
	
		fout_all_data = open('%s/Nodes_characteristics_nonTH.dat' %Exp_Dir,'a')
		print >> fout_all_data, 'Node','Clustering','Betweenness_centrality','Closeness_centrality','Eigenvector_centrality'
	
		#### Clustering ####
		clust_coefficients = nx.clustering(Conns_nonTH_ud)
		
		#### Betweenness centrality ####
		Conns_nonTH_components = list(nx.connected_component_subgraphs(Conns_nonTH_ud))
	
		Conns_nonTH_mc = Conns_nonTH_components[0]
	
	
		bet_cen = nx.betweenness_centrality(Conns_nonTH_mc)
	
		#### Closeness centrality ####
		clo_cen = nx.closeness_centrality(Conns_nonTH_mc)
	
	
		#### Eigenvector centrality ####
		eig_cen = nx.eigenvector_centrality(Conns_nonTH_mc)
	
		# PRINT TO FILE 
		results = [(k,clust_coefficients[k],bet_cen[k],clo_cen[k],eig_cen[k]) for k in Conns_nonTH_mc]
		for item in results:
			fout_all_data.write(','.join(map(str,item)))
			fout_all_data.write('\n')
		fout_all_data.close()
	elif TH == 1:
		###############################################
		################## TH #####################
		###############################################
	
		###### IN AND OUT DEGREES NonTH ##### 
		fout_deg = open('%s/In_and_Out_Degrees_TH.dat' %Exp_Dir,'w')
		In_degrees = Conns_TH.in_degree()
		In_values = sorted(set(In_degrees.values()))
		In_hist = [In_degrees.values().count(x) for x in In_values]
	
		Out_degrees = Conns_TH.out_degree()
		Out_values = sorted(set(Out_degrees.values()))
		Out_hist = [Out_degrees.values().count(x) for x in Out_values]
	
		for (i,j,k,l) in zip(*(In_values,In_hist,Out_values,Out_hist)):
			print >>fout_deg,i,j,k,l
	
		##### AVERAGE CLUSTERING COEFFICIENT ####
		fout_clust = open('%s/Avg_clustering_TH.dat' %Exp_Dir,'w')
	
		Conns_TH_ud = Conns_TH.to_undirected()
		ccs = nx.clustering(Conns_TH_ud)
		avg_clust = sum(ccs.values())/len(ccs)
		print >> fout_clust, avg_clust
	
		fout_all_data = open('%s/Nodes_characteristics_TH.dat' %Exp_Dir,'a')
		print >> fout_all_data, 'Node','Clustering','Betweenness_centrality','Closeness_centrality','Eigenvector_centrality'
	
		#### Clustering ####
		clust_coefficients = nx.clustering(Conns_TH_ud)
		
		#### Betweenness centrality ####
		Conns_TH_components = list(nx.connected_component_subgraphs(Conns_TH_ud))
	
		Conns_TH_mc = Conns_TH_components[0]
	
	
		bet_cen = nx.betweenness_centrality(Conns_TH_mc)
	
		#### Closeness centrality ####
		clo_cen = nx.closeness_centrality(Conns_TH_mc)
	
	
		#### Eigenvector centrality ####
		eig_cen = nx.eigenvector_centrality(Conns_TH_mc)
	
		# PRINT TO FILE 
		results = [(k,clust_coefficients[k],bet_cen[k],clo_cen[k],eig_cen[k]) for k in Conns_TH_mc]
		for item in results:
			fout_all_data.write(','.join(map(str,item)))
			fout_all_data.write('\n')
		fout_all_data.close()
