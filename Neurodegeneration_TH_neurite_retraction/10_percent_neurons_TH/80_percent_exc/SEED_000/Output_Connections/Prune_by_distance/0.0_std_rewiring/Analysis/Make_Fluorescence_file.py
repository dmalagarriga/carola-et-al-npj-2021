#!/usr/bin/env python
# -*- coding: utf-8 -*-

from math import floor
from datetime import datetime
import os
import numpy as np
#from scipy.optimize import curve_fit
#import matplotlib.pyplot as plt
#import matplotlib.image as mpimg
#import random
#import matplotlib.cm as cm

####

vers='_v27_'

#use_sub = True
use_sub = False


#use_noise = True
use_noise = True


sample='000'
Dir=os.getenv('P_Dir')
I = int(os.getenv('I'))

## Input file 
File = np.loadtxt('%s/Data_Raster_%04i.dat' %(Dir,I),unpack=True)

##

#
#   Parameters: are constants for each spike record
#

#F = 2*50.0        # frame rate in Herz
#F = 20.0 #25.0 # 25.0 #2000.0 #50.0/2        # frame rate in Herz
F = 20.0
CONV_F = (F/1000)

NMIN=0
nNeurons=int(max(File[1]))

N = nNeurons+2   #     100        # number of neurons
#N = 2100                  # number of neurons
MAX_T = max(File[0]) # A[0][-1]+1000 # 6900
print N, MAX_T
#MAX_T = NMIN*60*1000       #  5*60*1000    # simulation length in ms
max_TJ =  int( floor( MAX_T * CONV_F ) )
N_FRAMES = max_TJ

TAU   = 1000    # in ms
SIGMA = 0.01    # std dev of noise : unit ?

KD = 300    #   rate constant for saturation of calcium in the cell

#noiseON = 1

#####

if use_noise:
    noiseON = 1.0

if use_sub:
    SUB='Subset'
    N=N/10
else:
    SUB=''

#did = "v5circ_6r1307_"+SUB
did = SUB

herz = str(int(F))+'hz_'
times = str(NMIN) + "m_"

#pref = 'calcium_fluo_' + times + herz +did+vers
#pref = 'calcium_fluo_' + times + herz +did
pref = './Analysis/Fluorescence/fluorescence_'+sample+ '_' + times + herz +did

sp_pref = './Analysis/Fluorescence/spikes_'+sample+ '_' + times + herz +did


###


def main():
    
    print '# run this filtering '+vers+' at : ', datetime.now()
    startTime = datetime.now()
    
    in_spikes = '%s/Data_Raster_%04i' %(Dir,I)
    #in_net='NetworkStructure.txt'
    #in_spikes='../../SpikeRecord'+SUB
    #in_spikes='SpikeRecordSubset'
    #in_spikes='SpikeRecordSubset100'
    #in_spikes='SpikeRecord_filtered500_v03_v5circ1307_t600000_v01'
    #in_spikes='SpikeRecord_filtered500_v03_v5circ1307_t'+str(MAX_T)+'_v01'
    #
    m = 0 # number of header comments
    #
    make_fluo(in_spikes,m)
    #
    print '#  it took :', datetime.now() - startTime
#
#
#
def test():
    
    print '# make calcium fluorescence '+vers+' at : ', datetime.now()
    startTime = datetime.now()
    
    #in_net='NetworkStructure.txt'
    in_spikes='SpikeRecordSubset'
    m = 0 # number of header comments

    make_fluo(in_spikes,m)
    
    print '#  it took :', datetime.now() - startTime
#
#    #in_net='NetworkStructure.txt'
#    in_conn='cons2000.txt'
#    m = 9 # number of header comments
#
#    convert_connect(in_conn,m)
#
##

def     make_fluo(in_spikes,m):
    """ Function doc """
    
    startTime = datetime.now()
    print '#  1( load data '
    #
    #   load data: directly into n skip spike trace
    #
    spikes_binned = load_spikes(in_spikes,m)
    
    nextTime=datetime.now()
    print '#    it took :', nextTime - startTime
    startTime=nextTime

    #print len(spikes_binned)
    #print ' in data shape', spikes_binned.shape
    #print ' in data size', spikes_binned.size

    #
    #   show the type and size of spikes_binned
    #
    print '#        type of spikes', type(spikes_binned)
    print '#        size of spikes', len(spikes_binned)
    print '#        size of spikes', spikes_binned.size
    print '#        size of spikes', spikes_binned.shape

    file_output_spikes = sp_pref+'_all.txt'
    file_output_noburst = './Analysis/Fluorescence/noburst_conditioning.txt'
    file_output_bcond = './Analysis/Fluorescence/burst_conditioning.txt'

    file_output_bcond_diff = './Analysis/Fluorescence/burst_conditioning_differential.txt'

    #
    #   convert to 0 and 1
    #
    
    spl = (spikes_binned > 0).astype(int)
    print '#        size of spl', spl.shape
    
    spl=np.transpose(spl) + 1

    np.savetxt( file_output_spikes, spl, fmt='%d', delimiter=' ' )
    #np.savetxt( file_output_spikes, spikes_binned, fmt='%.6g', delimiter=' ' )

    nb = spl[:,0]*0+1

    np.savetxt( file_output_noburst, nb, fmt='%d', delimiter=' ' )





    print '#  2( calculate calcium '
    #
    #   calculate calcium part
    #
    ca = calc_fluo(spikes_binned)

    nextTime=datetime.now()
    print '#    it took :', nextTime - startTime
    startTime=nextTime
    
    print '#  3( add noise '
    #
    #   add noise part
    #   
    fluo = ca/( ca + KD ) + noiseON*noise()

    nextTime=datetime.now()
    print '#    it took :', nextTime - startTime
    startTime=nextTime
    #
    #   TODO: add light scattering
    #
    ## put the call here 
    
    print '#  4( get average fluo signal '
    #
    #   get average fluo signal
    #
    fluo_ave = np.mean(fluo,0)

    nextTime=datetime.now()
    print '#    it took :', nextTime - startTime    
    startTime=nextTime


    #
    #   get conditioning estimate and plot
    #
    #   a) conditioning on the values
    #

    CL = 0.15

    bcond = ( fluo_ave > CL ) + 1

    np.savetxt( file_output_bcond, bcond, fmt='%d', delimiter=' ' )


    #
    #   b) conditioning on the differences
    #
    #
    #   convert to 0 and 1
    #

    CL = CL*1.0/2.0
    
    fa_shift = np.append ( fluo_ave[1:], [0] )
    
    d_fa  = fluo_ave - fa_shift
    
    bcond_d = ( d_fa > CL ) + 1

    np.savetxt( file_output_bcond_diff, bcond_d, fmt='%d', delimiter=' ' )



    print '#  5( write output '
    #
    #   Write output
    #
    mt=''
    file_output_1 = pref+'_avg.txt'
    np.savetxt( file_output_1, fluo_ave, fmt='%.6g', delimiter=' ' )
    #fmt='%.18e', 
    #write_data(file_output_1,  fluo_ave(1,:) )

    file_output_2 = pref+'_all.dat'
    tfluo = np.transpose(fluo)
    #write_data(file_output_2,  fluo, mt )
    np.savetxt( file_output_2, tfluo, fmt='%.7g', delimiter=' ' )
    #numpy.savetxt(fname, X, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ')
    
    nextTime=datetime.now()
    print '#    it took :', nextTime - startTime

    print '# done. '
    
    
def  load_spikes( infile, k ):
    """ Function doc """
    #max_TJ =  floor( MAX_T[0] * CONV_F )
    z = np.zeros( (N,N_FRAMES) )

    with open(infile+'.dat','r') as f:
        #
        for line in f:
            [ j, tj ] = get_spike(line)
            #z += truc 
            z[j,tj] += 1
    return z    
    
    
#
def  get_spike(myline):
    """         get data
            format: index, t
    """
    values = myline.split()
    #print values
    i = int( values[1] )
    ti = int( floor( float(values[0]) * CONV_F ) ) -1 
    return [ i, ti ]
    
    
def  calc_fluo(n):
    """ Function doc 
        calculate the calcium fluorescence signal of spikes
    """       
    #   TODO in fine use this, 
    #   a = (1-dt/tau);
    # at 50 Hz (frames per second) a  = 0.98 . 
    dt = 1.0/F*1000
    a  = (1 - dt/TAU)
    #a  = 0.98
    b  = 50    #    in muM
    
    ca = np.zeros( (N,N_FRAMES) )
    ca[:,0] = b * n[:,0]

    for j in xrange(1,max_TJ):
        k = j-1
        ca[:,j] = a * ca[:,k] + b * n[:,j]
        
    return ca
    
    
def  noise():
    """ Function doc 
        calculate the noise term in fluorescence signal
    """
    w = SIGMA * np.random.randn(N,N_FRAMES);
    return w

####
#
####

####
#
##
#
###



####
#
####

####
#
#
##################################################################

####
#   OUTPUT
#
def write_data(fname, outdata, mt):
    
    f=file(fname, 'w')
    #f=file(fname, 'a')
    #f.write(mt+'\n')
    
    for line in outdata:
        #print line
        #s=to_string(line)
        s="{0: 11.8f}, {1: 11.8f}\n".format(line[1],line[2])
        f.write(s)
        # print s
    f.close()

def to_string(p):
    #s=" "+p[0]
    s=" "
    for x in p[1:]:
        s+="{0: 11.8f}, ".format(x)
    s+=' \n'
    return s
#
##


    
####
#   INPUT: read data
#
#   gets data in format: format: index, x y (z)
#
def  get_line(myline):
    """         get data
            format: index, x y (z)
    """
    values=myline.split()
    s=values[0]
    data_line=[int(s)]+[float(x) for x  in values[1:] ]
    return data_line
#
#
####    


#
#
#
##
#
###





###
#
#test()
#
main()
#
# eof
