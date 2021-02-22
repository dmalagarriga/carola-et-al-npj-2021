#!/bin/bash

#source Environment_Variables.bash

#function Crea_Arxiu {
#cat >> ${File_Out} << EOF
#200 
#EOF
#return
#}

#printf "Creating input files ... "
#for i in $(seq 1)
#    do
#    Num=$(echo "${i}" | awk '{printf "%06.3f",$1*0.01}')
#    File_Out="Input/P.dat"
#    Crea_Arxiu
#done
#printf "done\n"



export Output_Map_Dir="../../neurongen-master/Input"
export Dir_Output="Output"
export Dir_Input="Input"
export NumPar_Dat="Input/NumPar.dat"
export VoxPar_Dat="Input/VoxPar.dat"
export Parameters_Dat="Input/Parameters.dat"
export Kex1_Dat="Input/Kex1.dat"
export Kex2_Dat="Input/Kex2.dat"
export Kex_Dat="Input/Kex.dat"
export Kin_Dat="Input/Kin.dat"
export Kth_Dat="Input/Kth.dat"
export Delay_Matrix_Dat="Input/Delay_Matrix.dat"
export Data_Neurons_Dat="Output/Data_Neurons.dat"
export Data_Calcium_Dat="Output/Data_Calcium.dat"
export Data_Raster_Dat="Output/Data_Raster.dat"
export Data_Voxels_Dat="Output/Data_Voxel.dat"
export Data_Exc_Inh_Dat="Input/Exc_Inh_List.dat"
export Data_Dopaminergic_Dat="Input/Dopaminergic_List.dat"
export Data_Number_of_Connections_Dat="Input/Number_of_Connections_all.dat"
export Data_Connections_Dat="Input/Connections_all.dat"
export Prop_exc_inh=0.8
export I=0
export KTH=0
export GSL_RNG_SEE=20
export SEED
LANG=C
export LANG
LC_ALL=C
export LC_ALL


Compile=0
Run=1
Preprocess=0
Analysis=0

printf "Running ...\n"


if [ $Preprocess -eq 1 ]; then
printf "Preprocessing and creating connections...\n"
#python ./Preprocessing/Import_Positions.py
python ./Preprocessing/Convert_from_neurongen.py
printf "..Created!\n"
fi

if [ $Compile -eq 1 ]; then
printf "..Compiling...\n"
make
printf "...Compiled!\n"
fi


for KTH in 1 
	do
	printf -v Num_kth "%04i" $KTH
	mkdir -p "Output/Kth${Num_kth}"
	export P_Dir="Output/Kth${Num_kth}"
	for SEED in 1
		do
		for (( I=0; I <=0; I+=1)) 
			do
			printf -v Num_I "%04i" $I
			Data_Neurons_Dat="Output/Kth${Num_kth}/Data_Neurons_${Num_I}.dat" 
			Data_Calcium_Dat="Output/Kth${Num_kth}/Data_Calcium_${Num_I}.dat"
			Data_Raster_Dat="Output/Kth${Num_kth}/Data_Raster_${Num_I}.dat"
			if [ $Run -eq 1 ]; then
			printf "...Running simulation... %04i\n" $SEED
			./nmm.exe
			fi
			if [ $Analysis -eq 1 ]; then
			printf "...Analysing data..%04i\n"
			./Analysis/Analysis.bash
			printf "...Analysed!\n"
			fi
			done
	done
done
#./figures.bash
printf "..all done!\n"
#rm -f nmm.exe

	


