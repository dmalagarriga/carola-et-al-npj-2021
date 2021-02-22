#!/bin/bash

export Percentage_retracted=10
export Percentage_retraction=80
export Number_std=0.0
export Prop_exc_inh=0.8 # Proportion of excitation | 1.0 --> All exc, 0.0 no exc
export Configuration=Output_Connections #Input_Connections # Output_Connections
export SEED=002

#mkdir -p Neurodegeneration_TH_neurite_retraction/`echo $Percentage_retracted`_percent_neurons_TH/`echo $Number_std`_std_rewiring


###### Exporting directories

#export P_Dir=Neurodegeneration_TH_neurite_retraction/`echo $Percentage_retracted`_percent_neurons_TH/80_percent_exc/`echo $Configuration`/Prune_by_distance/`echo $Number_std`_std_rewiring/Input
#export Data_Dir=Neurodegeneration_TH_neurite_retraction/`echo $Percentage_retracted`_percent_neurons_TH/80_percent_exc/`echo $Configuration`/Prune_by_distance/`echo $Number_std`_std_rewiring/Analysis/Fluorescence_FunData
#export Fluo_Dir=Neurodegeneration_TH_neurite_retraction/`echo $Percentage_retracted`_percent_neurons_TH/80_percent_exc/`echo $Configuration`/Prune_by_distance/`echo $Number_std`_std_rewiring/Analysis/Fluorescence

export P_Dir=Neurodegeneration_TH_neurite_retraction/`echo $Percentage_retracted`_percent_neurons_TH/80_percent_exc/SEED_`echo $SEED`/`echo $Configuration`/Random_pruning/`echo $Percentage_retraction`_percent_retraction/Input
export Data_Dir=Neurodegeneration_TH_neurite_retraction/`echo $Percentage_retracted`_percent_neurons_TH/80_percent_exc/SEED_`echo $SEED`/`echo $Configuration`/Random_pruning/`echo $Percentage_retraction`_percent_retraction/Analysis/Fluorescence_FunData
export Fluo_Dir=Neurodegeneration_TH_neurite_retraction/`echo $Percentage_retracted`_percent_neurons_TH/80_percent_exc/SEED_`echo $SEED`/`echo $Configuration`/Random_pruning/`echo $Percentage_retraction`_percent_retraction/Analysis/Fluorescence
export Dir_Input=Neurodegeneration_TH_neurite_retraction/`echo $Percentage_retracted`_percent_neurons_TH/80_percent_exc/SEED_`echo $SEED`/`echo $Configuration`/Random_pruning/`echo $Percentage_retraction`_percent_retraction/Input

#export P_Dir=Neurodegeneration_TH_neurite_retraction/`echo $Percentage_retracted`_percent_neurons_TH/80_percent_exc/`echo $Configuration`/Rewiring/`echo $Number_std`_std_rewiring/Input
#export Data_Dir=Neurodegeneration_TH_neurite_retraction/`echo $Percentage_retracted`_percent_neurons_TH/80_percent_exc/`echo $Configuration`/Rewiring/`echo $Number_std`_std_rewiring/Analysis/Fluorescence_FunData
#export Fluo_Dir=Neurodegeneration_TH_neurite_retraction/`echo $Percentage_retracted`_percent_neurons_TH/80_percent_exc/`echo $Configuration`/Rewiring/`echo $Number_std`_std_rewiring/Analysis/Fluorescence

#export P_Dir=Neurodegeneration_TH_neurite_retraction/No_retraction_longer_neurite/Input
#export Data_Dir=Neurodegeneration_TH_neurite_retraction/No_retraction_longer_neurite/Analysis/Fluorescence_FunData
#export Fluo_Dir=Neurodegeneration_TH_neurite_retraction/No_retraction_longer_neurite/Analysis/Fluorescence

##### Python executables

#python ./Analysis/MUA.py
python ./Analysis/Draw_network_positions.py
#python ./Analysis/Analysis_Networks.py
#python ./Analysis/Avalanche_analysis.py
#python ./Analysis/Frequency_ABRUPT_EVENTS.py
#python ./Analysis/Ratio_ABRUPT_EVENTS.py