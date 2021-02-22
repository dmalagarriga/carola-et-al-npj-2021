#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_
#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include "../Headers/neuron.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"


using namespace std;

class Parameters{
//Atributs
	public:
		
	int it; //Contador
	int nt;	// Nombre de passos total
    int ntout; // Nombre de pasos per l'output
	int neq; //Number of equations

	float Prop_exc_inh;

	float dt_neuron; 

	float t_neuron;
 
	float KTH;
    string Dir_Output;
    string Dir_Input;
    string NumPar_Dat;
	string Data_Neurons_Dat;
	string Data_Raster_Dat;
	string Data_Calcium_Dat;
	string Data_Exc_Inh_Dat;
	string Data_Dopaminergic_Dat;
	
	string Delay_Matrix_Dat;
	string Data_Connections_Dat;
	string Data_Number_of_Connections_Dat;
	
    
	        
	int nNeurons;
	int nDopaminergic;
	int nConnections;
   
	
	int ** Connections;
	int * Delay_Matrix;
	int * Exc_Inh_array;
	int * Dopaminergic_array;
	int * Number_of_Connections;
	
	
	neuron * N; //Neuron
	
	
	ofstream fout_data_neurons;
	ofstream fout_data_calcium;
	ofstream fout_raster;
		

    public:
    Parameters();
    ~Parameters();
    
	void Read_Delay_Matrix_Dat(void);
	void Read_NumPar_Dat(void);
	void Read_Data_Exc_Inh_Dat(void);
	void Read_Data_Dopaminergic_Dat(void);
	void Read_Data_Connections_Dat(void);
	void Read_Data_Number_of_Connections_Dat(void);
	void RefreshDelay(void);
    
	void WriteData();
	void openfiles();
	void closefiles();
	

};

#endif
