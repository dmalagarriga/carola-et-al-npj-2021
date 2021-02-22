#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
import networkx as nx
import matplotlib.pyplot as plt
import pylab
import os
import numpy as np


File_Dat1="./Input/Kex.dat"
File_Dat2="./Input/Kin.dat"

nNodes=int(os.getenv('Nodes'))
initial_nodes=int(os.getenv('Initial_Nodes'))
Seed=int(os.getenv('Seed'))
#d=int(os.getenv('Degree'))
#BARABASI=>Scale-free from Barabasi-Albert algorithm
#SCALE_FREE=>Scale_free directed graph
#SMALL_WORLD=>Watts Strogatz small world network
#SMALL_WORLD_CONNECTED=>small world network watts strogratz connected
#SMALL_WORLD_NEWMAN=>Return a Newman-Watts-Strogatz small world graph.

Which_network=str(os.getenv('Network'))
File_Dat3="./Input/Degree_distribution_%s_%s.dat" %(Which_network,initial_nodes)
File_Dat4="./Input/Edge_distribution_%s_%s.dat" %(Which_network,initial_nodes)
fout1=open(File_Dat1,'w')
fout2=open(File_Dat2,'w')
fout3=open(File_Dat3,'w')
fout4=open(File_Dat4,'w')
print >>fout1,nNodes
print >>fout2,nNodes 

if(Which_network=='BARABASI'):
	G=nx.barabasi_albert_graph(nNodes,initial_nodes,seed=Seed)
	H=G.to_directed()
	for u,v,d in sorted(H.edges(data=True)): #u,v -> pairs, d=weight
		d['weight']=1./H.in_degree(u)
		#print H.in_degree(u)
		#print >> fout3,(H.in_degree(v))
	np.savetxt(fout1,nx.adjacency_matrix(H),'%5.3f')
	np.savetxt(fout2,nx.adjacency_matrix(H),'%5.3f')
	nx.draw(H)
	pylab.savefig('./Input/Network_Barabasi_init_nodes_%s.ps' %initial_nodes)
	pylab.close()
	degree_sequence=sorted(nx.degree(H).values(),reverse=True)
	for i in degree_sequence:	
		print >> fout3,i
	for edges in G.edges():
		print >> fout4,edges
	plt.loglog(degree_sequence,'bo-')
	plt.xlabel('Node')
	plt.ylabel('Degree')
	pylab.savefig('./Input/Degree_Dist_Barabasi_%s.ps' %initial_nodes)
	pylab.close()
	fout1.close()
	fout2.close()

if(Which_network=='SCALE_FREE'):
	G=nx.scale_free_graph(nNodes,alpha=0.41, beta=0.54, gamma=0.05, delta_in=0.2, delta_out=0, create_using=None, seed=None)
	H=G.to_directed()
	for u,v,d in H.edges(data=True):
		d['weight']=1./H.in_degree(u)
	np.savetxt(fout1,nx.adjacency_matrix(H),'%5.3f')
	np.savetxt(fout2,nx.adjacency_matrix(H),'%5.3f')
	nx.draw(H,'b')
	pylab.savefig('./Input/Network_Scale_Free_Directed.ps')
	fout1.close()
	fout2.close()
	
if(Which_network=='SMALL_WORLD'):
	G=nx.watts_strogatz_graph(nNodes,2,p=0.0,seed=1)
	H=G.to_directed()
	for u,v,d in H.edges(data=True):
		d['weight']=1./H.in_degree(u)
	np.savetxt(fout1,nx.adjacency_matrix(H),'%3.1f')
	np.savetxt(fout2,nx.adjacency_matrix(H),'%3.1f')
	nx.draw_circular(G)
	#nx.draw(H)
	pylab.savefig('./Input/Network_Small_World_%s.ps' %initial_nodes)
	fout1.close()
	fout2.close()
	pylab.close()

if(Which_network=='SMALL_WORLD_CONNECTED'):
	G=nx.connected_watts_strogatz_graph(nNodes,2,p=0.2,seed=None)
	H=G.to_directed()
	for u,v,d in H.edges(data=True):
		d['weight']=1./H.in_degree(u)
	np.savetxt(fout1,nx.adjacency_matrix(H),'%5.3f')
	np.savetxt(fout2,nx.adjacency_matrix(H),'%5.3f')
	nx.draw(H)
	pylab.savefig('./Input/Network_Small_World_Connected.ps')
	fout1.close()
	fout2.close()

if(Which_network=='SMALL_WORLD_NEWMAN'):
	G=nx.newman_watts_strogatz_graph(nNodes,2,p=1,seed=1)
	H=G.to_directed()
	for u,v,d in H.edges(data=True):
		d['weight']=1./H.in_degree(u)
	np.savetxt(fout1,nx.adjacency_matrix(H),'%5.3f')
	np.savetxt(fout2,nx.adjacency_matrix(H),'%5.3f')
	nx.draw(H)
	pylab.savefig('./Input/Network_Small_World_Newman.ps')
	fout1.close()
	fout2.close()
d=nNodes-45 #Degree of each node in this random regular
if(Which_network=='RANDOM_REGULAR'):
	G=nx.random_regular_graph(d,nNodes,seed=Seed)
	H=G.to_directed()
	for u,v,d in H.edges(data=True):
		d['weight']=1./H.in_degree(u)
	np.savetxt(fout1,nx.adjacency_matrix(H),'%5.3f')
	np.savetxt(fout2,nx.adjacency_matrix(H),'%5.3f')
	#nx.draw_circular(G)
	nx.draw(H)
	pylab.savefig('./Input/Network_Random_Regular.ps')
	fout1.close()
	fout2.close()
	
