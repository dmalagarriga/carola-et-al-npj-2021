#ifndef _NEURON_H_
#define _NEURON_H_

//Class Neuron

//#include <cstdio>
#include <string>
#include <iostream>
#include <cstdlib>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

extern gsl_rng * r;

using namespace std;
class neuron
{
	public:
	
	double * V ;	
	double * V_old ;
	double * V_der;
	double * V_der_old;
	
	double * n_K ;
	double * n_K_old ;
	double * n_K_i ;
	double * n_K_i_old ;
	double * n_K_der ;
	double * n_K_der_old;
	
	double * h_Na ;
	double * h_Na_i;
	double * h_Na_old;
	double * h_Na_i_old;
	double * h_Na_der ;
	double * h_Na_der_old;
	
	double * m_Na;
	double * m_Na_old;
	double * m_Na_der;
	double * m_Na_der_old;
	
	
	
	
	double * h_A;
	double * m_KS;
	double * Conc_Ca;
	double * Conc_Na;
	double * s_AMPA;
	double * s_NMDA;
	double * s_GABA;
	double * x_AMPA;
	double * x_NMDA;
	double * s_ext_AMPA;
	double * s_ext_NMDA;
	double * x_ext_AMPA;
	double * x_ext_NMDA;
	
	double * h_A_old;
	double * m_KS_old;
	double * Conc_Ca_old;
	double * Conc_Na_old;
	double * s_AMPA_old;
	double * s_NMDA_old;
	double * s_GABA_old;
	double * x_AMPA_old;
	double * x_NMDA_old;
	double * s_ext_AMPA_old;
	double * s_ext_NMDA_old;
	double * x_ext_AMPA_old;
	double * x_ext_NMDA_old;
	
	double * h_A_der;
	double * m_KS_der;
	double * Conc_Ca_der;
	double * Conc_Na_der;
	double * s_AMPA_der;
	double * s_NMDA_der;
	double * s_GABA_der;
	double * x_AMPA_der;
	double * x_NMDA_der;
	double * s_ext_AMPA_der;
	double * s_ext_NMDA_der;
	double * x_ext_AMPA_der;
	double * x_ext_NMDA_der;
	
	
	double * h_A_der_old;
	double * m_KS_der_old;
	double * Conc_Ca_der_old;
	double * Conc_Na_der_old;
	double * s_AMPA_der_old;
	double * s_NMDA_der_old;
	double * s_GABA_der_old;
	double * x_AMPA_der_old;
	double * x_NMDA_der_old;
	double * s_ext_AMPA_der_old;
	double * s_ext_NMDA_der_old;
	double * x_ext_AMPA_der_old;
	double * x_ext_NMDA_der_old;
	
	double * Xi ;
	double * Xi_old ;
	double * Noisy_rate;
	double * Noisy_rate_old;
	
	
	double * V_0;
	
	
	int nDelays;
	int Delay;
	
	double Rate;
	double ext_rate;
	double tau_ext_rate;
	double amplitude_ext_rate;
	double g_syn_soma;
	double g_syn_den;
	double C;
	double A_soma;
	double A_den;
	double A_i;
	double Temp_fact;
	double tau_h_A;
	double alpha_Ca;
	double tau_Ca;
	double alpha_Na;
	double R_pump;
	double D;
	double Conc_Na_eq;
	double alpha_x_AMPA;
	double alpha_x_NMDA;
	double alpha_s_AMPA;
	double alpha_s_NMDA;
	double alpha_GABA;
	double tau_s_AMPA;
	double tau_s_NMDA;
	double tau_s_GABA;
	double tau_x_AMPA;
	double tau_x_NMDA;
	double g_K;
	double g_Na;
	double g_L;
	double g_A;
	double g_KS;
	double g_NaP;
	double g_AR;
	double g_Ca;
	double g_KCa;
	double g_KNa;
	double g_AMPA; 
	double g_NMDA;
	double g_GABA;
	double g_AMPA_exc;
	double g_NMDA_exc;
	double g_ext_AMPA_exc;
	double g_ext_NMDA_exc;
	double g_ext_NMDA;
	double g_ext_AMPA;
	double g_GABA_exc;
	double g_AMPA_inh;
	double g_NMDA_inh;
	double g_GABA_inh;
	double g_ext_AMPA_inh;
	double g_ext_NMDA_inh;
	double g_K_i ;
	double g_Na_i;
	double g_L_i;
	double K_D;
	double V_Ca;
	double V_K;
	double V_Na;
	double V_L;
	double V_K_i;
	double V_Na_i;
	double V_L_i;
	double V_syn_GABA;
	double epsilon;
	double I_K_i;
	double I_K;
	double I_Na;
	double I_L;
	double I_A;
	double I_KS;
	double I_NaP;
	double I_AR;
	double I_Ca;
	double I_KCa;
	double I_KNa;
	double I_syn_soma;
	double I_syn_den;
	double I_syn_i;
	double I_AMPA;
	double I_NMDA;
	double I_ext_AMPA;
	double I_ext_NMDA;
	double I_GABA;
	double I_AMPA_exc;
	double I_NMDA_exc;
	double I_GABA_exc;
	double I_AMPA_inh;
	double I_NMDA_inh;
	double I_GABA_inh;
	double I_Na_i;
	double I_L_i;
	double m_inf_Na;
	double m_inf_Na_i;
	double Vpeak;
	double Conc_Mg;
	
	string Parameters_Dat;
		
	 
	public:	
	
	
	void set_Initial_Conditions(void);
	void Dim(void);
	void Read_Parameters_Dat(void);
    
		
	neuron();
	~neuron();
};
#endif
