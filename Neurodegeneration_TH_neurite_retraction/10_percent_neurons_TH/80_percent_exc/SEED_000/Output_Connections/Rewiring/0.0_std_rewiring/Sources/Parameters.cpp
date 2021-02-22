#include "../Headers/Parameters.h"
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <iomanip>

/* Constructor */
Parameters::Parameters(){
cout << "Creating parameters object." << endl;


  Dir_Output=string( getenv("Dir_Output") ); 
  Dir_Input=string( getenv("Dir_Input") ); 
  NumPar_Dat=string( getenv("NumPar_Dat") );
  Data_Neurons_Dat=string( getenv("Data_Neurons_Dat") );
  Data_Calcium_Dat=string( getenv("Data_Calcium_Dat") );
  Data_Exc_Inh_Dat=string( getenv("Data_Exc_Inh_Dat") );
  Data_Raster_Dat=string( getenv ("Data_Raster_Dat") );
  Data_Dopaminergic_Dat=string( getenv("Data_Dopaminergic_Dat") );
  Data_Connections_Dat=string( getenv("Data_Connections_Dat") );
  Data_Number_of_Connections_Dat=string( getenv("Data_Number_of_Connections_Dat") );
 
  Prop_exc_inh=atof(getenv("Prop_exc_inh"));
  const char* e = getenv("Delay_Matrix_Dat");
  Delay_Matrix_Dat=string(e);
 
  t_neuron=0;

  Read_NumPar_Dat();
  Read_Data_Exc_Inh_Dat();
  Read_Data_Connections_Dat();
  Read_Data_Number_of_Connections_Dat();
  Read_Delay_Matrix_Dat();
  Read_Data_Dopaminergic_Dat();
  
  
  

  N=new neuron[nNeurons];
  

  //Create (give dimensions) voxels
  
  openfiles();
  
  
};

/* Destructor */ 
Parameters::~Parameters(){

cout << "Destroying parameters object." << endl;

delete [] Dopaminergic_array;
delete [] Exc_Inh_array;
delete [] Delay_Matrix;
delete [] Connections;
delete [] Number_of_Connections;

delete [] N;


closefiles();
};


/* Read Data to build the Parameter object */

void Parameters::Read_Data_Dopaminergic_Dat(void){

string str;
ifstream fin (Data_Dopaminergic_Dat.c_str(),ifstream::in);
fin >> nDopaminergic;     getline(fin,str);
cout << "nDopaminergic:" << nDopaminergic << endl;

Dopaminergic_array=new int[nDopaminergic];
for (int i=0;i<nDopaminergic;i++){fin>>Dopaminergic_array[i];}

fin.close();

return;
};


void Parameters::Read_Data_Exc_Inh_Dat(void){

string str;
ifstream fin (Data_Exc_Inh_Dat.c_str(),ifstream::in);
fin >> nNeurons;     getline(fin,str);
cout << "nNeurons:" << nNeurons << endl;


Exc_Inh_array=new int[nNeurons-nDopaminergic];
for (int i=0;i<nNeurons-nDopaminergic;i++){fin>>Exc_Inh_array[i];}

fin.close();


return;
};


void Parameters::Read_Data_Connections_Dat(void){

 string str;
 ifstream fin (Data_Connections_Dat.c_str(),ifstream::in);
 fin >> nConnections;     getline(fin,str);
 cout << "nConnections:" << nConnections << endl;

 Connections=new int*[nConnections];
 
 for (int i=0;i<nConnections;i++){Connections[i]=new int[nConnections];}
 
 for (int i=0;i<nConnections;i++){
 		
 		for (int j=0;j<2;j++){fin>>Connections[i][j];}
 		
 		}

fin.close();
return;

};

void Parameters::Read_Data_Number_of_Connections_Dat(void){

 ifstream fin (Data_Number_of_Connections_Dat.c_str(),ifstream::in);
 Number_of_Connections=new int[nNeurons];
 
 
 for (int i=0;i<nNeurons;i++){fin>> Number_of_Connections[i];}

fin.close();
return;

};


void Parameters::Read_Delay_Matrix_Dat(void){
 
 ifstream fin (Delay_Matrix_Dat.c_str(),ifstream::in);
 
 Delay_Matrix=new int[nConnections];
	

	for (int i=0;i<nConnections;i++){fin >> Delay_Matrix[i];}		  
				  
 fin.close();

 return;
};
/* Read Data to build the Parameter object */
void Parameters::Read_NumPar_Dat(void){
 
 string str;
 ifstream fin (NumPar_Dat.c_str(),ifstream::in);
 fin >> nt;  getline(fin,str);
 cout << "nt: " << nt << endl;
 fin >> ntout; getline(fin,str);
 cout << "ntout: " <<ntout<<endl;
 fin >> dt_neuron; getline(fin,str);
 cout << "dt_neuron: " <<dt_neuron<<endl;

 KTH = atof(getenv("KTH"));
 cout << "KTH: " << KTH << endl;
 fin.close();

 return;
};

void Parameters::RefreshDelay(void){

	
    for (int j=0; j<nNeurons; j++){
	
        for (int i=N->nDelays; i>0; i--){	
		
            N[j].V_0[i]=N[j].V_0[i-1];
		
                                        }				    
	N[j].V_0[0]=N[j].V[0];
	
		
    }			
				     				    
				};
//Write data to file
void Parameters::openfiles(void){
fout_data_neurons.open(Data_Neurons_Dat.c_str(),ofstream::out);
fout_data_calcium.open(Data_Calcium_Dat.c_str(),ofstream::out);
fout_raster.open(Data_Raster_Dat.c_str(),ofstream::out);
 return;
};

void Parameters::closefiles(void){
 fout_data_neurons.close();
fout_data_calcium.close();
fout_raster.close();
  return;
};


void Parameters::WriteData(void){

fout_data_neurons.precision(8);
fout_raster.precision(8);

fout_data_neurons << t_neuron ;

/* WRITE RASTER PLOT */
	//for (int i=0;i < nNeurons; i++){if (N[Exc_Inh_array[i]].V[0]>=N[Exc_Inh_array[i]].Vpeak){fout_raster << t_neuron << " " << Exc_Inh_array[i] << endl;}}
	for (int i=0;i < nNeurons; i++){if (N[i].V[0]>=N[i].Vpeak){fout_raster << t_neuron << " " << i << endl;}}
	//fout_raster << endl;

/* WRITE DATA SPIKES OF FEW NEURONS */

	for (int i=0; i<3;i++){fout_data_neurons << " " << N[Exc_Inh_array[i]].V[0];}
	if (Prop_exc_inh!=0){
	for (int i=((nNeurons-nDopaminergic)-Prop_exc_inh*(nNeurons-nDopaminergic)+1); i<((nNeurons-nDopaminergic)-Prop_exc_inh*(nNeurons-nDopaminergic)+3);i++){fout_data_neurons << " " << N[Exc_Inh_array[i]].V[0];}
	}for (int i = 0; i<3;i++){fout_data_neurons << " " << N[Dopaminergic_array[i]].V[0];}
	//for (int i=123; i<123+3;i++){fout_data_neurons << " " << N[Exc_Inh_array[i]].V[0];}
	fout_data_neurons << endl;
	
	

/* WRITE RASTER AS MATRIX OF 0s AND 1s */
	//fout_raster << t_neuron;
	//	for (int i=0;i < nNeurons; i++){if (N[Exc_Inh_array[i]].V[0]>=N[Exc_Inh_array[i]].Vpeak){fout_raster << " " << 1;}else{fout_raster << " " << 0;}}
	//fout_raster << endl;

/* WRITE CALCIUM FLUORESCENCE SIGNAL */
	//fout_data_calcium << t_neuron ;
	//for(int i=0; i< nNeurons;i++){fout_data_neurons << " " << N[i].V[0];/*fout_data_calcium << " " << N[i].Conc_Ca[0]/(N[i].Conc_Ca[0]+300);*/}
	//fout_data_calcium << endl;
	
 return;
};


