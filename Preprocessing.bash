#!/bin/bash

export Percentage_retracted=50
export Percentage_retraction=80
export Number_std=0.0
export Prop_exc_inh=0.8 # Proportion of excitation | 1.0 --> All exc, 0.0 no exc
export Configuration=Output_Connections # Input_Connections # Output_Connections
export SEED=002

#mkdir -p Neurodegeneration_TH_neurite_retraction/`echo $Percentage_retracted`_percent_neurons_TH/`echo $Number_std`_std_rewiring


###### Exporting directories


export Dir_Input=Neurodegeneration_TH_neurite_retraction/`echo $Percentage_retracted`_percent_neurons_TH/80_percent_exc/SEED_`echo $SEED`/`echo $Configuration`/Random_pruning/`echo $Percentage_retraction`_percent_retraction/Input
#export Dir_Input=Neurodegeneration_TH_neurite_retraction/No_retraction_longer_neurite/Input
#export Dir_Input=Neurodegeneration_TH_neurite_retraction/`echo $Percentage_retracted`_percent_neurons_TH/80_percent_exc/`echo $Configuration`/Prune_by_distance/`echo $Number_std`_std_rewiring/Input
#export Dir_Input=Neurodegeneration_TH_neurite_retraction/`echo $Percentage_retracted`_percent_neurons_TH/80_percent_exc/`echo $Configuration`/Rewiring/`echo $Number_std`_std_rewiring/Input


##### Python executables

#python ./Analysis/Reduce_number_of_connections_by_distance.py #<- First
#python ./Analysis/Reduce_number_of_connections_randomly.py #<- First
#python ./Analysis/Rewire_connections.py
python ./Analysis/Convert_from_neurongent.py #<- Second

