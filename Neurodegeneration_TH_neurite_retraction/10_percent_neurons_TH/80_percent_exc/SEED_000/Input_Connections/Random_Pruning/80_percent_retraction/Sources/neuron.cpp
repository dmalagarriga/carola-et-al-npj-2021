//#include <iostream>
//#include <cstdlib>

#include "../Headers/neuron.h"
#include "../Headers/Parameters.h"

using namespace std;

/*Constructor*/
neuron::neuron(){
	
    const char* s =getenv("Parameters_Dat");
    Parameters_Dat =string(s);
    
    
	nDelays = 200;
	Dim();
	set_Initial_Conditions();
    Read_Parameters_Dat();
    
   
};

/*Constructors*/
void neuron::Dim(void){
	V = new double[2];	
	V_old = new double[2];
	V_der = new double[2];	
	V_der_old = new double[2];
	
	
	n_K = new double[1];
	n_K_old =new double[1];
	n_K_i = new double[1];
	n_K_i_old =new double[1];
	n_K_der = new double[1];
	n_K_der_old = new double[1];
	
	h_Na = new double[1];
	h_Na_i= new double[1];
	h_Na_old = new double[1];
	h_Na_i_old= new double[1];
	h_Na_der = new double[1];
	h_Na_der_old= new double[1];
	
	m_Na = new double[1];
	m_Na_old = new double[1];
	m_Na_der = new double[1];
	m_Na_der_old = new double[1];
	
	h_A=new double[1];
	m_KS=new double[1];
	Conc_Ca=new double[1];
	Conc_Na=new double[1];
	s_AMPA=new double[1];
	s_NMDA=new double[1];
	s_GABA=new double[1];
	x_AMPA=new double[1];
	x_NMDA=new double[1];
	s_ext_AMPA=new double[1];
	s_ext_NMDA=new double[1];
	x_ext_AMPA=new double[1];
	x_ext_NMDA=new double[1];
	
	h_A_old=new double[1];
	m_KS_old=new double[1];
	Conc_Ca_old=new double[1];
	Conc_Na_old=new double[1];
	s_AMPA_old=new double[1];
	s_NMDA_old=new double[1];
	s_ext_AMPA_old=new double[1];
	s_ext_NMDA_old=new double[1];
	
	s_GABA_old=new double[1];
	
	x_AMPA_old=new double[1];
	x_NMDA_old=new double[1];
	x_ext_AMPA_old=new double[1];
	x_ext_NMDA_old=new double[1];
	
	h_A_der=new double[1];
	m_KS_der=new double[1];
	Conc_Ca_der=new double[1];
	Conc_Na_der=new double[1];
	s_AMPA_der=new double[1];
	s_NMDA_der=new double[1];
	s_ext_AMPA_der=new double[1];
	s_ext_NMDA_der=new double[1];
	
	s_GABA_der=new double[1];
	
	x_AMPA_der=new double[1];
	x_NMDA_der=new double[1];
	x_ext_AMPA_der=new double[1];
	x_ext_NMDA_der=new double[1];
	
	h_A_der_old=new double[1];
	m_KS_der_old=new double[1];
	Conc_Ca_der_old=new double[1];
	Conc_Na_der_old=new double[1];
	s_AMPA_der_old=new double[1];
	s_NMDA_der_old=new double[1];
	s_ext_AMPA_der_old=new double[1];
	s_ext_NMDA_der_old=new double[1];
	
	s_GABA_der_old=new double[1];
	
	x_AMPA_der_old=new double[1];
	x_NMDA_der_old=new double[1];
	x_ext_AMPA_der_old=new double[1];
	x_ext_NMDA_der_old=new double[1];
	
	Xi = new double[1];
	Xi_old = new double[1];
	Noisy_rate = new double[1];
	Noisy_rate_old = new double[1];
	
	
	V_0 = new double[1];
	
		
//Delay arrays
	V_0 = new double[nDelays];
	for (int i=0; i < nDelays; i++){
	V_0[i]=0;
		
	}
	return;	
};
void neuron::Read_Parameters_Dat(void){
    
    string str;
    ifstream fin (Parameters_Dat.c_str(),ifstream::in);
    fin >> g_syn_soma ; getline(fin,str);
	fin >> g_syn_den; getline(fin,str);
	fin >> C; getline(fin,str);
	fin >> A_soma; getline(fin,str);
	fin >> A_den; getline(fin,str);
	fin >> A_i; getline(fin,str);
	fin >> Temp_fact; getline(fin,str);
	fin >> tau_h_A ; getline(fin,str);
	fin >> alpha_Ca; getline(fin,str);
	fin >> tau_Ca; getline(fin,str);
	fin >> alpha_Na; getline(fin,str);
	fin >> R_pump ; getline(fin,str);
	fin >> D ; getline(fin,str);
	fin >> Conc_Na_eq; getline(fin,str);
	fin >> alpha_x_AMPA; getline(fin,str);
	fin >> alpha_x_NMDA; getline(fin,str);
	fin >> alpha_s_AMPA; getline(fin,str);
	fin >> alpha_s_NMDA; getline(fin,str);
	fin >> alpha_GABA; getline(fin,str);
	fin >> tau_s_AMPA; getline(fin,str);
	fin >> tau_s_NMDA; getline(fin,str);
	fin >> tau_s_GABA; getline(fin,str);
	fin >> tau_x_AMPA; getline(fin,str);
	fin >> tau_x_NMDA; getline(fin,str);
	fin >> g_K; getline(fin,str);
	fin >> g_Na; getline(fin,str);
	fin >> g_L; getline(fin,str);
	fin >> g_A; getline(fin,str);
	fin >> g_KS; getline(fin,str);
	fin >> g_NaP; getline(fin,str);
	fin >> g_AR; getline(fin,str);
	fin >> g_Ca; getline(fin,str);
	fin >> g_KCa; getline(fin,str);
	fin >> g_KNa; getline(fin,str);
	fin >> g_ext_AMPA_exc; getline(fin,str);
	fin >> g_ext_NMDA_exc; getline(fin,str);
	fin >> g_ext_AMPA_inh; getline(fin,str);
	fin >> g_ext_NMDA_inh; getline(fin,str);
	fin >> g_AMPA_exc; getline(fin,str);
	fin >> g_NMDA_exc; getline(fin,str);
	fin >> g_GABA_exc; getline(fin,str);
	fin >> g_AMPA_inh; getline(fin,str);
	fin >> g_NMDA_inh; getline(fin,str);
	fin >> g_GABA_inh; getline(fin,str);
	fin >> g_K_i ; getline(fin,str);
	fin >> g_Na_i; getline(fin,str);
	fin >> g_L_i ; getline(fin,str);
	fin >> K_D; getline(fin,str);
	fin >> V_Ca ; getline(fin,str);
	fin >> V_K ; getline(fin,str);
	fin >> V_Na ; getline(fin,str);
	fin >> V_L ; getline(fin,str);
	fin >> V_K_i; getline(fin,str);
	fin >> V_Na_i; getline(fin,str);
	fin >> V_L_i ; getline(fin,str);
	fin >> V_syn_GABA; getline(fin,str);
	fin >> epsilon; getline(fin,str);
	fin >> Vpeak; getline(fin,str);
	fin >> Conc_Mg; getline(fin,str);
	fin >> tau_ext_rate; getline(fin,str);
	fin >> ext_rate; getline(fin,str);
    fin >> amplitude_ext_rate; getline(fin,str); 
     
    fin.close();
    
    return;
};


/* Destructor */

neuron::~neuron(){

	delete [] V_0;
	delete [] V;	
	delete [] V_old ;
	delete [] n_K;
	delete [] n_K_i;
	delete [] n_K_old;
	delete [] n_K_i_old;
	delete [] n_K_der;
	delete [] n_K_der_old ;
	
	
	delete [] h_Na;
	delete [] h_Na_i;
	delete [] h_Na_old;
	delete [] h_Na_i_old;
	delete [] h_Na_der;
	delete [] h_Na_der_old;
	
	delete [] m_Na;
	delete [] m_Na_old;
	delete [] m_Na_der;
	delete [] m_Na_der_old;
	
	delete [] h_A;
	delete [] m_KS;
	delete [] Conc_Ca;
	delete [] Conc_Na;
	delete [] s_AMPA;
	delete [] s_NMDA;
	delete [] s_ext_AMPA;
	delete [] s_ext_NMDA;
	delete [] s_GABA;
	delete [] x_AMPA;
	delete [] x_NMDA;
	delete [] x_ext_AMPA;
	delete [] x_ext_NMDA;
	
	delete [] h_A_der;
	delete [] m_KS_der;
	delete [] Conc_Ca_der;
	delete [] Conc_Na_der;
	delete [] s_AMPA_der;
	delete [] s_NMDA_der;
	delete [] s_GABA_der;
	delete [] x_AMPA_der;
	delete [] x_NMDA_der;
	delete [] s_ext_AMPA_der;
	delete [] s_ext_NMDA_der;

	delete [] x_ext_AMPA_der;
	delete [] x_ext_NMDA_der;
	
	
	delete [] h_A_der_old;
	delete [] m_KS_der_old;
	delete [] Conc_Ca_der_old;
	delete [] Conc_Na_der_old;
	delete [] s_AMPA_der_old;
	delete [] s_NMDA_der_old;
	delete [] s_GABA_der_old;
	delete [] x_AMPA_der_old;
	delete [] x_NMDA_der_old;
	delete [] s_ext_AMPA_der_old;
	delete [] s_ext_NMDA_der_old;

	delete [] x_ext_AMPA_der_old;
	delete [] x_ext_NMDA_der_old;
	
	delete [] Xi;
	delete [] Xi_old;
	delete [] Noisy_rate;
	delete [] Noisy_rate_old;
	
};

void neuron::set_Initial_Conditions(void){
V[0]=-60.0;
V[1]=-60.0;
V_old[0] = 0;
V_der[0] = 0;	
V_der_old[0] = 0;
V_old[1] = 0;
V_der[1] = 0;	
V_der_old[1] = 0;
	
n_K[0]= 0;
n_K_old[0]= 0;
n_K_i[0]= 0;
n_K_i_old[0]= 0;
n_K_der[0]= 0;
n_K_der_old[0]= 0;

h_Na[0]= 0;
h_Na_i[0]= 0;
h_Na_old[0]= 0;
h_Na_i_old[0]= 0;
h_Na_der[0]= 0;
h_Na_der_old[0]= 0;

h_A[0]=0;
m_KS[0]=0;
Conc_Ca[0]=0;
Conc_Na[0]=0;
s_AMPA[0]=0;
s_NMDA[0]=0;
s_GABA[0]=0;
x_AMPA[0]=0;
x_NMDA[0]=0;
s_ext_AMPA[0]=0;
s_ext_NMDA[0]=0;
x_AMPA[0]=0;
x_NMDA[0]=0;
x_ext_AMPA[0]=0;
x_ext_NMDA[0]=0;

h_A_old[0]=0;
m_KS_old[0]=0;
Conc_Ca_old[0]=0;
Conc_Na_old[0]=0;
s_AMPA_old[0]=0;
s_NMDA_old[0]=0;
s_GABA_old[0]=0;
x_AMPA_old[0]=0;
x_NMDA_old[0]=0;
s_ext_AMPA_old[0]=0;
s_ext_NMDA_old[0]=0;
x_ext_AMPA_old[0]=0;
x_ext_NMDA_old[0]=0;

h_A_der[0]=0;
m_KS_der[0]=0;
Conc_Ca_der[0]=0;
Conc_Na_der[0]=0;
s_AMPA_der[0]=0;
s_NMDA_der[0]=0;
s_GABA_der[0]=0;
x_AMPA_der[0]=0;
x_NMDA_der[0]=0;
s_ext_AMPA_der[0]=0;
s_ext_NMDA_der[0]=0;
x_ext_AMPA_der[0]=0;
x_ext_NMDA_der[0]=0;

h_A_der_old[0]=0;
m_KS_der_old[0]=0;
Conc_Ca_der_old[0]=0;
Conc_Na_der_old[0]=0;
s_AMPA_der_old[0]=0;
s_NMDA_der_old[0]=0;
s_GABA_der_old[0]=0;
x_AMPA_der_old[0]=0;
x_NMDA_der_old[0]=0;
s_ext_AMPA_der_old[0]=0;
s_ext_NMDA_der_old[0]=0;
x_ext_AMPA_der_old[0]=0;
x_ext_NMDA_der_old[0]=0;

Xi[0] = 0;
Xi_old[0] = 0;
Noisy_rate[0] = 0;
Noisy_rate_old[0] = 0;

V_0[0] = 0;
return;
};	
