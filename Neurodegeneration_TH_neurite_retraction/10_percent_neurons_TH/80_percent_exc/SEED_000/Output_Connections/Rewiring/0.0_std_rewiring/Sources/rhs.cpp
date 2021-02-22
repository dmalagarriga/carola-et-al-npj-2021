#include "../Headers/rhs.h"
#include "../Headers/Parameters.h"
#include "../Headers/neuron.h"
#include <math.h>
#include <numeric>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

/* Functions for excitatory subpopulation */
double alpha_n_K(neuron *N, double nu){
	
	double value;
	value=0.01*((nu+34.0)/(1-exp(-(nu+34.0)/10.0)));
	    return value;
}
double beta_n_K(neuron *N, double nu){
	
	double value;
	value=0.125*exp(-(nu+44.0)/25.0);
        return value;
}

double alpha_m_Na(neuron *N, double nu){
	
	double value;
	value=0.1*((nu+33.0)/(1-exp(-(nu+33.0)/10.0)));
        return value;
}
double beta_m_Na(neuron *N, double nu){
	
	double value;
	value=4.0*exp(-(nu+53.7)/12.0);
        return value;
}

double alpha_h_Na(neuron *N, double nu){
	
	double value;
	value=0.07*exp(-(nu+50.0)/10.0);
        return value;
}
double beta_h_Na(neuron *N, double nu){
	
	double value;
	value=1.0/(1.0+exp(-(nu+20.0)/10.0));
        return value;
}


double m_inf_A(neuron *N, double nu){
	
	double value;
	value=1.0/(1.0+exp(-(nu+50.0)/20.0));
        return value;
}
double h_inf_A(neuron *N, double nu){
	
	double value;
	value=1.0/(1.0+exp((nu+80.0)/6.0));
        return value;
}

double tau_m_A(neuron *N, double nu){
	
	double value;
	value=1.0/(1.0+exp(-(nu+50.0)/20.0));
        return value;
}

double m_inf_KS(neuron *N, double nu){
	
	double value;
	value=1.0/(1.0+exp(-(nu+34.0)/6.5));
        return value;
}
double tau_m_KS(neuron *N, double nu){
	
	double value;
	value=8.0/(exp(-(nu+55.0)/30.0)+exp(+(nu+55.0)/30.0));
        return value;
}

double m_inf_NaP(neuron *N, double nu){
	
	double value;
	value=1.0/(1.0+exp(-(nu+55.7)/7.7));
        return value;
}

double h_inf_AR(neuron *N, double nu){
	
	double value;
	value=1.0/(1.0+exp((nu+75.0)/4.0));
        return value;
}

double m_inf_Ca(neuron *N, double nu){
	
	double value;
	value=1.0/(1.0+exp(-(nu+20.0)/9.0));
        return value;
}

double w_inf(neuron *N, double nu){
	
	double value;
	value=0.37/(1.0+pow(38.7/nu,3.5));
        return value;
}

/* Functions for inhibitory subpopulation */

double alpha_n_K_i(neuron *N, double nu){
	
	double value;
	value=0.05*((nu+34.0)/(1-exp(-(nu+34.0)/10.0)));
        return value;
}
double beta_n_K_i(neuron *N, double nu){
	
	double value;
	value=0.625*exp(-(nu+44.0)/80.0);
        return value;
}

double alpha_m_Na_i(neuron *N, double nu){
	
	double value;
	value=0.5*((nu+35.0)/(1-exp(-(nu+35.0)/10.0)));
        return value;
}
double beta_m_Na_i(neuron *N, double nu){
	
	double value;
	value=20.0*exp(-(nu+60.0)/18.0);
        return value;
}

double alpha_h_Na_i(neuron *N, double nu){
	
	double value;
	value=0.35*exp(-(nu+58.0)/20.0);
        return value;
}
double beta_h_Na_i(neuron *N, double nu){
	
	double value;
	value=5.0/(1.0+exp(-(nu+28.0)/10.0));
        return value;
}



void rhs (Parameters *P){
/*
  for (int i=0; i<(P->nNeurons-P->nDopaminergic);i++){
   P->N[P->Exc_Inh_array[i]].Rate = P->N[P->Exc_Inh_array[i]].ext_rate + 0.5*P->N[P->Exc_Inh_array[i]].amplitude_ext_rate*P->N[P->Exc_Inh_array[i]].Noisy_rate[0];
   if (P->N[P->Exc_Inh_array[i]].Rate<0.0){P->N[P->Exc_Inh_array[i]].Rate=0.0;}
  }*/
    /* EXCITATORY NEURONS */
    for (int i=0; i<(P->nNeurons-P->nDopaminergic)-(1.0-P->Prop_exc_inh)*(P->nNeurons-P->nDopaminergic);i++){
    	 P->N[P->Exc_Inh_array[i]].Rate = P->N[P->Exc_Inh_array[i]].ext_rate + 0.5*P->N[P->Exc_Inh_array[i]].amplitude_ext_rate*P->N[P->Exc_Inh_array[i]].Noisy_rate[0];
    	 if (P->N[P->Exc_Inh_array[i]].Rate<0.0){P->N[P->Exc_Inh_array[i]].Rate=0.0;}
    	
    	/* DERIVATIVES */
    	
        //Dopaminergic (exc) neurons derivatives and rhs
        P->N[P->Exc_Inh_array[i]].V_der[0] = (-P->N[P->Exc_Inh_array[i]].A_soma*(P->N[P->Exc_Inh_array[i]].I_L + P->N[P->Exc_Inh_array[i]].I_Na + P->N[P->Exc_Inh_array[i]].I_K + P->N[P->Exc_Inh_array[i]].I_A + P->N[P->Exc_Inh_array[i]].I_KS + P->N[P->Exc_Inh_array[i]].I_KNa) - P->N[P->Exc_Inh_array[i]].I_syn_soma - P->N[P->Exc_Inh_array[i]].g_syn_soma*(P->N[P->Exc_Inh_array[i]].V[0]-P->N[P->Exc_Inh_array[i]].V[1]))/(P->N[P->Exc_Inh_array[i]].C*P->N[P->Exc_Inh_array[i]].A_soma);
         
        P->N[P->Exc_Inh_array[i]].V_der[1] = (-P->N[P->Exc_Inh_array[i]].A_den*(P->N[P->Exc_Inh_array[i]].I_Ca + P->N[P->Exc_Inh_array[i]].I_KCa + P->N[P->Exc_Inh_array[i]].I_NaP + P->N[P->Exc_Inh_array[i]].I_AR) - P->N[P->Exc_Inh_array[i]].I_syn_den - P->N[P->Exc_Inh_array[i]].g_syn_den*(P->N[P->Exc_Inh_array[i]].V[1]-P->N[P->Exc_Inh_array[i]].V[0]))/(P->N[P->Exc_Inh_array[i]].C*P->N[P->Exc_Inh_array[i]].A_den);
        
        
    	// Derivatives of gating variables
    	//Exc (soma)
    	P->N[P->Exc_Inh_array[i]].n_K_der[0] = P->N[P->Exc_Inh_array[i]].Temp_fact*(alpha_n_K(&P->N[P->Exc_Inh_array[i]],P->N[P->Exc_Inh_array[i]].V[0])*(1-P->N[P->Exc_Inh_array[i]].n_K[0]) - beta_n_K(&P->N[P->Exc_Inh_array[i]],P->N[P->Exc_Inh_array[i]].V[0])*P->N[P->Exc_Inh_array[i]].n_K[0]);
    	P->N[P->Exc_Inh_array[i]].h_Na_der[0] = P->N[P->Exc_Inh_array[i]].Temp_fact*(alpha_h_Na(&P->N[P->Exc_Inh_array[i]],P->N[P->Exc_Inh_array[i]].V[0])*(1-P->N[P->Exc_Inh_array[i]].h_Na[0]) - beta_h_Na(&P->N[P->Exc_Inh_array[i]],P->N[P->Exc_Inh_array[i]].V[0])*P->N[P->Exc_Inh_array[i]].h_Na[0]);
    	P->N[P->Exc_Inh_array[i]].m_Na_der[0] = P->N[P->Exc_Inh_array[i]].Temp_fact*(alpha_m_Na(&P->N[P->Exc_Inh_array[i]],P->N[P->Exc_Inh_array[i]].V[0])*(1-P->N[P->Exc_Inh_array[i]].m_Na[0]) - beta_m_Na(&P->N[P->Exc_Inh_array[i]],P->N[P->Exc_Inh_array[i]].V[0])*P->N[P->Exc_Inh_array[i]].m_Na[0]);
    	
    	
    	P->N[P->Exc_Inh_array[i]].h_A_der[0] = (1.0/P->N[P->Exc_Inh_array[i]].tau_h_A)*(h_inf_A(&P->N[P->Exc_Inh_array[i]],P->N[P->Exc_Inh_array[i]].V[0])-P->N[P->Exc_Inh_array[i]].h_A[0]);
    	P->N[P->Exc_Inh_array[i]].m_KS_der[0] = (1.0/tau_m_KS(&P->N[P->Exc_Inh_array[i]],P->N[P->Exc_Inh_array[i]].V[0]))*(m_inf_KS(&P->N[P->Exc_Inh_array[i]],P->N[P->Exc_Inh_array[i]].V[0])-P->N[P->Exc_Inh_array[i]].m_KS[0]);
    	
    	
    	// Derivative of the concentration of Ca and Na
    	
    	P->N[P->Exc_Inh_array[i]].Conc_Ca_der[0] = -P->N[P->Exc_Inh_array[i]].alpha_Ca*P->N[P->Exc_Inh_array[i]].A_den*P->N[P->Exc_Inh_array[i]].I_Ca - P->N[P->Exc_Inh_array[i]].Conc_Ca[0]/P->N[P->Exc_Inh_array[i]].tau_Ca;
    	
    	P->N[P->Exc_Inh_array[i]].Conc_Na_der[0] = -P->N[P->Exc_Inh_array[i]].alpha_Na*(P->N[P->Exc_Inh_array[i]].A_soma*P->N[P->Exc_Inh_array[i]].I_Na + P->N[P->Exc_Inh_array[i]].A_den*P->N[P->Exc_Inh_array[i]].I_NaP) - (P->N[P->Exc_Inh_array[i]].R_pump + P->N[P->Exc_Inh_array[i]].D*P->N[P->Exc_Inh_array[i]].Xi[0])*( pow(P->N[P->Exc_Inh_array[i]].Conc_Na[0],3)/(pow(P->N[P->Exc_Inh_array[i]].Conc_Na[0],3)+pow(15.0,3)) - pow(P->N[P->Exc_Inh_array[i]].Conc_Na_eq,3)/(pow(P->N[P->Exc_Inh_array[i]].Conc_Na_eq,3)+pow(15.0,3) ) );
    	
    	// Derivative of time-dependent conductances
    	P->N[P->Exc_Inh_array[i]].s_AMPA_der[0] = P->N[P->Exc_Inh_array[i]].alpha_s_AMPA*P->N[P->Exc_Inh_array[i]].x_AMPA[0]*(1-P->N[P->Exc_Inh_array[i]].s_AMPA[0])-P->N[P->Exc_Inh_array[i]].s_AMPA[0]/P->N[P->Exc_Inh_array[i]].tau_s_AMPA;
    	P->N[P->Exc_Inh_array[i]].s_NMDA_der[0] = P->N[P->Exc_Inh_array[i]].alpha_s_NMDA*P->N[P->Exc_Inh_array[i]].x_NMDA[0]*(1-P->N[P->Exc_Inh_array[i]].s_NMDA[0])-P->N[P->Exc_Inh_array[i]].s_NMDA[0]/P->N[P->Exc_Inh_array[i]].tau_s_NMDA;
    	
    	P->N[P->Exc_Inh_array[i]].s_ext_AMPA_der[0] = P->N[P->Exc_Inh_array[i]].alpha_s_AMPA*P->N[P->Exc_Inh_array[i]].x_ext_AMPA[0]*(1-P->N[P->Exc_Inh_array[i]].s_ext_AMPA[0])-P->N[P->Exc_Inh_array[i]].s_ext_AMPA[0]/P->N[P->Exc_Inh_array[i]].tau_s_AMPA;
    	P->N[P->Exc_Inh_array[i]].s_ext_NMDA_der[0] = P->N[P->Exc_Inh_array[i]].alpha_s_NMDA*P->N[P->Exc_Inh_array[i]].x_ext_NMDA[0]*(1-P->N[P->Exc_Inh_array[i]].s_ext_NMDA[0])-P->N[P->Exc_Inh_array[i]].s_ext_NMDA[0]/P->N[P->Exc_Inh_array[i]].tau_s_NMDA;
    	
    	
    	P->N[P->Exc_Inh_array[i]].s_GABA_der[0] = -P->N[P->Exc_Inh_array[i]].s_GABA[0]/P->N[P->Exc_Inh_array[i]].tau_s_GABA;
    	
    	
    	P->N[P->Exc_Inh_array[i]].x_AMPA_der[0] = -P->N[P->Exc_Inh_array[i]].x_AMPA[0]/P->N[P->Exc_Inh_array[i]].tau_x_AMPA;
    	P->N[P->Exc_Inh_array[i]].x_NMDA_der[0] = -P->N[P->Exc_Inh_array[i]].x_NMDA[0]/P->N[P->Exc_Inh_array[i]].tau_x_NMDA; 
    	
    	P->N[P->Exc_Inh_array[i]].x_ext_AMPA_der[0] = -P->N[P->Exc_Inh_array[i]].x_ext_AMPA[0]/P->N[P->Exc_Inh_array[i]].tau_x_AMPA;
    	P->N[P->Exc_Inh_array[i]].x_ext_NMDA_der[0] = -P->N[P->Exc_Inh_array[i]].x_ext_NMDA[0]/P->N[P->Exc_Inh_array[i]].tau_x_NMDA; 
    	
    	// Gating variable exc 
    	P->N[P->Exc_Inh_array[i]].m_inf_Na = alpha_m_Na(&P->N[P->Exc_Inh_array[i]],P->N[P->Exc_Inh_array[i]].V[0])/(alpha_m_Na(&P->N[P->Exc_Inh_array[i]],P->N[P->Exc_Inh_array[i]].V[0])+beta_m_Na(&P->N[P->Exc_Inh_array[i]],P->N[P->Exc_Inh_array[i]].V[0]));
    	
    	// Currents
    	// Soma
    	P->N[P->Exc_Inh_array[i]].I_K = P->N[P->Exc_Inh_array[i]].g_K*pow(P->N[P->Exc_Inh_array[i]].n_K[0],4)*(P->N[P->Exc_Inh_array[i]].V[0]-P->N[P->Exc_Inh_array[i]].V_K); 
    	
    	//P->N[P->Exc_Inh_array[i]].I_Na = P->N[P->Exc_Inh_array[i]].g_Na*pow(P->N[P->Exc_Inh_array[i]].m_Na[0],3)*P->N[P->Exc_Inh_array[i]].h_Na[0]*(P->N[P->Exc_Inh_array[i]].V[0]-P->N[P->Exc_Inh_array[i]].V_Na);
    	P->N[P->Exc_Inh_array[i]].I_Na = P->N[P->Exc_Inh_array[i]].g_Na*pow(P->N[P->Exc_Inh_array[i]].m_inf_Na,3)*P->N[P->Exc_Inh_array[i]].h_Na[0]*(P->N[P->Exc_Inh_array[i]].V[0]-P->N[P->Exc_Inh_array[i]].V_Na);
    	
    	P->N[P->Exc_Inh_array[i]].I_L = P->N[P->Exc_Inh_array[i]].g_L*(P->N[P->Exc_Inh_array[i]].V[0]-P->N[P->Exc_Inh_array[i]].V_L); 
    	
    	P->N[P->Exc_Inh_array[i]].I_A = P->N[P->Exc_Inh_array[i]].g_A*pow(m_inf_A(&P->N[P->Exc_Inh_array[i]],P->N[P->Exc_Inh_array[i]].V[0]),3)*P->N[P->Exc_Inh_array[i]].h_A[0]*(P->N[P->Exc_Inh_array[i]].V[0]-P->N[P->Exc_Inh_array[i]].V_K);
    	
    	P->N[P->Exc_Inh_array[i]].I_KS = P->N[P->Exc_Inh_array[i]].g_KS*P->N[P->Exc_Inh_array[i]].m_KS[0]*(P->N[P->Exc_Inh_array[i]].V[0]-P->N[P->Exc_Inh_array[i]].V_K);
    	
    	P->N[P->Exc_Inh_array[i]].I_KNa = P->N[P->Exc_Inh_array[i]].g_KNa*w_inf(&P->N[P->Exc_Inh_array[i]],P->N[P->Exc_Inh_array[i]].Conc_Na[0])*(P->N[P->Exc_Inh_array[i]].V[0]-P->N[P->Exc_Inh_array[i]].V_K);
    	
    	 	
    	// Dendrites
    	
    	P->N[P->Exc_Inh_array[i]].I_NaP = P->N[P->Exc_Inh_array[i]].g_NaP*pow(m_inf_NaP(&P->N[P->Exc_Inh_array[i]], P->N[P->Exc_Inh_array[i]].V[1]),3)*(P->N[P->Exc_Inh_array[i]].V[1]-P->N[P->Exc_Inh_array[i]].V_Na);
    	
    	P->N[P->Exc_Inh_array[i]].I_AR = P->N[P->Exc_Inh_array[i]].g_AR*h_inf_AR(&P->N[P->Exc_Inh_array[i]],P->N[P->Exc_Inh_array[i]].V[1])*(P->N[P->Exc_Inh_array[i]].V[1]-P->N[P->Exc_Inh_array[i]].V_K);
    	
    	P->N[P->Exc_Inh_array[i]].I_Ca = P->N[P->Exc_Inh_array[i]].g_Ca*pow(m_inf_Ca(&P->N[P->Exc_Inh_array[i]],P->N[P->Exc_Inh_array[i]].V[1]),2)*(P->N[P->Exc_Inh_array[i]].V[1]-P->N[P->Exc_Inh_array[i]].V_Ca);
    	
    	P->N[P->Exc_Inh_array[i]].I_KCa = P->N[P->Exc_Inh_array[i]].g_KCa*(P->N[P->Exc_Inh_array[i]].Conc_Ca[0]/(P->N[P->Exc_Inh_array[i]].Conc_Ca[0]+P->N[P->Exc_Inh_array[i]].K_D))*(P->N[P->Exc_Inh_array[i]].V[1]-P->N[P->Exc_Inh_array[i]].V_K);
    	
    	
    	
    	//Synaptic
    	P->N[P->Exc_Inh_array[i]].I_syn_soma = P->N[P->Exc_Inh_array[i]].I_GABA;
    	P->N[P->Exc_Inh_array[i]].I_syn_den = P->N[P->Exc_Inh_array[i]].I_AMPA + P->N[P->Exc_Inh_array[i]].I_NMDA + P->N[P->Exc_Inh_array[i]].I_ext_AMPA + P->N[P->Exc_Inh_array[i]].I_ext_NMDA;
    	
    	P->N[P->Exc_Inh_array[i]].I_AMPA = P->N[P->Exc_Inh_array[i]].g_AMPA*P->N[P->Exc_Inh_array[i]].s_AMPA[0]*(P->N[P->Exc_Inh_array[i]].V[1]);
    	P->N[P->Exc_Inh_array[i]].I_NMDA = P->N[P->Exc_Inh_array[i]].g_NMDA*P->N[P->Exc_Inh_array[i]].s_NMDA[0]*(P->N[P->Exc_Inh_array[i]].V[1])*(1.0/(1+P->N[P->Exc_Inh_array[i]].Conc_Mg*exp((-0.062*P->N[P->Exc_Inh_array[i]].V[1]/3.57))));
    	
    	P->N[P->Exc_Inh_array[i]].I_ext_AMPA = P->N[P->Exc_Inh_array[i]].g_ext_AMPA*P->N[P->Exc_Inh_array[i]].s_ext_AMPA[0]*(P->N[P->Exc_Inh_array[i]].V[1]);
    	P->N[P->Exc_Inh_array[i]].I_ext_NMDA = P->N[P->Exc_Inh_array[i]].g_ext_NMDA*P->N[P->Exc_Inh_array[i]].s_ext_NMDA[0]*(P->N[P->Exc_Inh_array[i]].V[1])*(1.0/(1+P->N[P->Exc_Inh_array[i]].Conc_Mg*exp((-0.062*P->N[P->Exc_Inh_array[i]].V[1]/3.57))));
    	
    	P->N[P->Exc_Inh_array[i]].I_GABA = P->N[P->Exc_Inh_array[i]].g_GABA*P->N[P->Exc_Inh_array[i]].s_GABA[0]*(P->N[P->Exc_Inh_array[i]].V[0]-P->N[P->Exc_Inh_array[i]].V_syn_GABA);
    	
    	
    	/* EXTERNAL TRAIN OF POISSON SPIKES */
    	
    	if((P->N[P->Exc_Inh_array[i]].Rate*P->dt_neuron/1000)>=gsl_ran_flat(r,0,1)){
    	
    	P->N[P->Exc_Inh_array[i]].x_ext_NMDA_der[0] += P->N[P->Exc_Inh_array[i]].alpha_x_NMDA;//P->N[P->Exc_Inh_array[i]].x_AMPA_der[0] += P->N[P->Exc_Inh_array[i]].alpha_x_AMPA;P->N[P->Exc_Inh_array[i]].s_GABA_der[0] += P->N[P->Exc_Inh_array[i]].alpha_GABA*(1-P->N[P->Exc_Inh_array[i]].s_GABA[0]);
    	P->N[P->Exc_Inh_array[i]].x_ext_AMPA_der[0] += P->N[P->Exc_Inh_array[i]].alpha_x_AMPA;
    	
    	}
    	
    	
    	/* EXCITATORY NEURONS INPUTS */
    	
    	int sum=0;
    	
    	for (int a=0;a<i;a++){sum+=P->Number_of_Connections[a-0];}
       
        
        for (int j=sum;j<sum+P->Number_of_Connections[i];j++){
        
       				
        
        			if(P->N[P->Connections[j][0]].V_0[P->Delay_Matrix[P->Connections[j][0]]] >= P->N[P->Connections[j][0]].Vpeak){
        			//if(P->N[P->Connections[j][0]].V[0] >= P->N[P->Connections[j][0]].Vpeak){
        			//P->N[P->Exc_Inh_array[i]].x_AMPA_der[0] += P->N[P->Exc_Inh_array[i]].alpha_x_AMPA;
        			//P->N[P->Exc_Inh_array[i]].s_GABA_der[0] += P->N[P->Exc_Inh_array[i]].alpha_GABA*(1-P->N[P->Exc_Inh_array[i]].s_GABA[0]);
        			//P->N[P->Exc_Inh_array[i]].x_NMDA_der[0] += P->N[P->Exc_Inh_array[i]].alpha_x_NMDA;
        			
        			// AMPA input
        			if (P->Connections[j][1] == 1){
        			P->N[P->Exc_Inh_array[i]].x_AMPA_der[0] += P->N[P->Exc_Inh_array[i]].alpha_x_AMPA;
        			
        			}
        			// GABA input
        			if (P->Connections[j][1] == 0){
        			P->N[P->Exc_Inh_array[i]].s_GABA_der[0] += P->N[P->Exc_Inh_array[i]].alpha_GABA*(1-P->N[P->Exc_Inh_array[i]].s_GABA[0]);
        			
        			}
        			// NMDA input
        			if (P->Connections[j][1] == 2){
        			P->N[P->Exc_Inh_array[i]].x_NMDA_der[0] += P->N[P->Exc_Inh_array[i]].alpha_x_NMDA;
        			}
        			
        			}
        
        
        
        }
        
        
    	}
    /* ---------------------------*/
    
    /* INHIBITORY NEURONS */
    
    
    for (int i=(P->nNeurons-P->nDopaminergic)-(1.0-P->Prop_exc_inh)*(P->nNeurons-P->nDopaminergic); i<(P->nNeurons-P->nDopaminergic);i++){
    	
    	P->N[P->Exc_Inh_array[i]].Rate = P->N[P->Exc_Inh_array[i]].ext_rate + 0.5*P->N[P->Exc_Inh_array[i]].amplitude_ext_rate*P->N[P->Exc_Inh_array[i]].Noisy_rate[0] ;
    	if (P->N[P->Exc_Inh_array[i]].Rate<0.0){P->N[P->Exc_Inh_array[i]].Rate=0.0;}
    	
    	/* DERIVATIVES */
    	
    	//Inhibitory neurons derivatives
    	P->N[P->Exc_Inh_array[i]].V_der[0] = (-P->N[P->Exc_Inh_array[i]].A_i*(P->N[P->Exc_Inh_array[i]].I_L_i + P->N[P->Exc_Inh_array[i]].I_Na_i + P->N[P->Exc_Inh_array[i]].I_K_i ) - P->N[P->Exc_Inh_array[i]].I_syn_i)/(P->N[P->Exc_Inh_array[i]].C*P->N[P->Exc_Inh_array[i]].A_i);
    	
    	//Inh gating variables
    	P->N[P->Exc_Inh_array[i]].n_K_der[0] = P->N[P->Exc_Inh_array[i]].Temp_fact*(alpha_n_K_i(&P->N[P->Exc_Inh_array[i]],P->N[P->Exc_Inh_array[i]].V[0])*(1-P->N[P->Exc_Inh_array[i]].n_K[0]) - beta_n_K_i(&P->N[P->Exc_Inh_array[i]],P->N[P->Exc_Inh_array[i]].V[0])*P->N[P->Exc_Inh_array[i]].n_K[0]);
    	P->N[P->Exc_Inh_array[i]].h_Na_der[0] = P->N[P->Exc_Inh_array[i]].Temp_fact*(alpha_h_Na_i(&P->N[P->Exc_Inh_array[i]],P->N[P->Exc_Inh_array[i]].V[0])*(1-P->N[P->Exc_Inh_array[i]].h_Na[0]) - beta_h_Na_i(&P->N[P->Exc_Inh_array[i]],P->N[P->Exc_Inh_array[i]].V[0])*P->N[P->Exc_Inh_array[i]].h_Na[0]);
    	P->N[P->Exc_Inh_array[i]].m_Na_der[0] = P->N[P->Exc_Inh_array[i]].Temp_fact*(alpha_m_Na_i(&P->N[P->Exc_Inh_array[i]],P->N[P->Exc_Inh_array[i]].V[0])*(1-P->N[P->Exc_Inh_array[i]].m_Na[0]) - beta_m_Na_i(&P->N[P->Exc_Inh_array[i]],P->N[P->Exc_Inh_array[i]].V[0])*P->N[P->Exc_Inh_array[i]].m_Na[0]);
    	
    	// Derivative of time-dependent conductances
    	P->N[P->Exc_Inh_array[i]].s_AMPA_der[0] = P->N[P->Exc_Inh_array[i]].alpha_s_AMPA*P->N[P->Exc_Inh_array[i]].x_AMPA[0]*(1-P->N[P->Exc_Inh_array[i]].s_AMPA[0])-P->N[P->Exc_Inh_array[i]].s_AMPA[0]/P->N[P->Exc_Inh_array[i]].tau_s_AMPA;
    	P->N[P->Exc_Inh_array[i]].s_NMDA_der[0] = P->N[P->Exc_Inh_array[i]].alpha_s_NMDA*P->N[P->Exc_Inh_array[i]].x_NMDA[0]*(1-P->N[P->Exc_Inh_array[i]].s_NMDA[0])-P->N[P->Exc_Inh_array[i]].s_NMDA[0]/P->N[P->Exc_Inh_array[i]].tau_s_NMDA;
    	
    	P->N[P->Exc_Inh_array[i]].s_ext_AMPA_der[0] = P->N[P->Exc_Inh_array[i]].alpha_s_AMPA*P->N[P->Exc_Inh_array[i]].x_ext_AMPA[0]*(1-P->N[P->Exc_Inh_array[i]].s_ext_AMPA[0])-P->N[P->Exc_Inh_array[i]].s_ext_AMPA[0]/P->N[P->Exc_Inh_array[i]].tau_s_AMPA;
    	P->N[P->Exc_Inh_array[i]].s_ext_NMDA_der[0] = P->N[P->Exc_Inh_array[i]].alpha_s_NMDA*P->N[P->Exc_Inh_array[i]].x_ext_NMDA[0]*(1-P->N[P->Exc_Inh_array[i]].s_ext_NMDA[0])-P->N[P->Exc_Inh_array[i]].s_ext_NMDA[0]/P->N[P->Exc_Inh_array[i]].tau_s_NMDA;
    	
    	
    	P->N[P->Exc_Inh_array[i]].s_GABA_der[0] = -P->N[P->Exc_Inh_array[i]].s_GABA[0]/P->N[P->Exc_Inh_array[i]].tau_s_GABA;
    	
    	
    	P->N[P->Exc_Inh_array[i]].x_AMPA_der[0] = -P->N[P->Exc_Inh_array[i]].x_AMPA[0]/P->N[P->Exc_Inh_array[i]].tau_x_AMPA;
    	P->N[P->Exc_Inh_array[i]].x_NMDA_der[0] = -P->N[P->Exc_Inh_array[i]].x_NMDA[0]/P->N[P->Exc_Inh_array[i]].tau_x_NMDA;  
    	
    	P->N[P->Exc_Inh_array[i]].x_ext_AMPA_der[0] = -P->N[P->Exc_Inh_array[i]].x_ext_AMPA[0]/P->N[P->Exc_Inh_array[i]].tau_x_AMPA;
    	P->N[P->Exc_Inh_array[i]].x_ext_NMDA_der[0] = -P->N[P->Exc_Inh_array[i]].x_ext_NMDA[0]/P->N[P->Exc_Inh_array[i]].tau_x_NMDA;  
    	
    	//Currents
    	// Inhibitory subpopulation
    	
    	P->N[P->Exc_Inh_array[i]].I_K_i = P->N[P->Exc_Inh_array[i]].g_K_i*pow(P->N[P->Exc_Inh_array[i]].n_K[0],4)*(P->N[P->Exc_Inh_array[i]].V[0]-P->N[P->Exc_Inh_array[i]].V_K_i); 
    	
    	P->N[P->Exc_Inh_array[i]].I_Na_i = P->N[P->Exc_Inh_array[i]].g_Na_i*pow(P->N[P->Exc_Inh_array[i]].m_Na[0],3)*P->N[P->Exc_Inh_array[i]].h_Na[0]*(P->N[P->Exc_Inh_array[i]].V[0]-P->N[P->Exc_Inh_array[i]].V_Na_i);
    	
    	P->N[P->Exc_Inh_array[i]].I_L_i = P->N[P->Exc_Inh_array[i]].g_L_i*(P->N[P->Exc_Inh_array[i]].V[0]-P->N[P->Exc_Inh_array[i]].V_L_i); 
    	
    	// Gating variable inh
    	P->N[P->Exc_Inh_array[i]].m_inf_Na_i = alpha_m_Na_i(&P->N[P->Exc_Inh_array[i]],P->N[P->Exc_Inh_array[i]].V[0])/(alpha_m_Na_i(&P->N[P->Exc_Inh_array[i]],P->N[P->Exc_Inh_array[i]].V[0])+beta_m_Na_i(&P->N[P->Exc_Inh_array[i]],P->N[P->Exc_Inh_array[i]].V[0]));
    	
    	//Synaptic
    	P->N[P->Exc_Inh_array[i]].I_syn_i = P->N[P->Exc_Inh_array[i]].I_GABA + P->N[P->Exc_Inh_array[i]].I_AMPA + P->N[P->Exc_Inh_array[i]].I_NMDA + P->N[P->Exc_Inh_array[i]].I_ext_NMDA + P->N[P->Exc_Inh_array[i]].I_ext_AMPA;
    	
    	
    	P->N[P->Exc_Inh_array[i]].I_AMPA = P->N[P->Exc_Inh_array[i]].g_AMPA*P->N[P->Exc_Inh_array[i]].s_AMPA[0]*(P->N[P->Exc_Inh_array[i]].V[0]);
    	P->N[P->Exc_Inh_array[i]].I_NMDA = P->N[P->Exc_Inh_array[i]].g_NMDA*P->N[P->Exc_Inh_array[i]].s_NMDA[0]*(P->N[P->Exc_Inh_array[i]].V[0])*(1/(1+P->N[P->Exc_Inh_array[i]].Conc_Mg*exp((-0.062*P->N[P->Exc_Inh_array[i]].V[1]/3.57))));
    	
    	P->N[P->Exc_Inh_array[i]].I_ext_AMPA = P->N[P->Exc_Inh_array[i]].g_ext_AMPA*P->N[P->Exc_Inh_array[i]].s_ext_AMPA[0]*(P->N[P->Exc_Inh_array[i]].V[0]);
    	P->N[P->Exc_Inh_array[i]].I_ext_NMDA = P->N[P->Exc_Inh_array[i]].g_ext_NMDA*P->N[P->Exc_Inh_array[i]].s_ext_NMDA[0]*(P->N[P->Exc_Inh_array[i]].V[0])*(1/(1+P->N[P->Exc_Inh_array[i]].Conc_Mg*exp((-0.062*P->N[P->Exc_Inh_array[i]].V[1]/3.57))));
    	
    	
    	P->N[P->Exc_Inh_array[i]].I_GABA = P->N[P->Exc_Inh_array[i]].g_GABA*P->N[P->Exc_Inh_array[i]].s_GABA[0]*(P->N[P->Exc_Inh_array[i]].V[0]-P->N[P->Exc_Inh_array[i]].V_syn_GABA);
    	
    	/* EXTERNAL TRAIN OF POISSON SPIKES */
    	if((P->N[P->Exc_Inh_array[i]].Rate*P->dt_neuron/1000)>=gsl_ran_flat(r,0,1)){
    	
    	P->N[P->Exc_Inh_array[i]].x_ext_NMDA_der[0] += P->N[P->Exc_Inh_array[i]].alpha_x_NMDA;//P->N[P->Exc_Inh_array[i]].x_AMPA_der[0] += P->N[P->Exc_Inh_array[i]].alpha_x_AMPA;P->N[P->Exc_Inh_array[i]].s_GABA_der[0] += P->N[P->Exc_Inh_array[i]].alpha_GABA*(1-P->N[P->Exc_Inh_array[i]].s_GABA[0]);
    	P->N[P->Exc_Inh_array[i]].x_ext_AMPA_der[0] += P->N[P->Exc_Inh_array[i]].alpha_x_AMPA;
    	}
    	
    	/* INHIBITORY NEURONS INPUTS */
    	
    	int sum=0;
    	
    	for (int a=0;a<i;a++){sum+=P->Number_of_Connections[a-0];}
       
        
        for (int j=sum;j<sum+P->Number_of_Connections[i];j++){
        
       				
        
        			if(P->N[P->Connections[j][0]].V_0[P->Delay_Matrix[P->Connections[j][0]]] >= P->N[P->Connections[j][0]].Vpeak){
        			//if(P->N[P->Connections[j][0]].V[0] >= P->N[P->Connections[j][0]].Vpeak){
        			//P->N[P->Exc_Inh_array[i]].x_AMPA_der[0] += P->N[P->Exc_Inh_array[i]].alpha_x_AMPA;
        			//P->N[P->Exc_Inh_array[i]].s_GABA_der[0] += P->N[P->Exc_Inh_array[i]].alpha_GABA*(1-P->N[P->Exc_Inh_array[i]].s_GABA[0]);
        			//P->N[P->Exc_Inh_array[i]].x_NMDA_der[0] += P->N[P->Exc_Inh_array[i]].alpha_x_NMDA;
        			
        			// AMPA input
        			if (P->Connections[j][1] == 1){
        			P->N[P->Exc_Inh_array[i]].x_AMPA_der[0] += P->N[P->Exc_Inh_array[i]].alpha_x_AMPA;
        			
        			}
        			// GABA input
        			if (P->Connections[j][1] == 0){
        			P->N[P->Exc_Inh_array[i]].s_GABA_der[0] += P->N[P->Exc_Inh_array[i]].alpha_GABA*(1-P->N[P->Exc_Inh_array[i]].s_GABA[0]);
        			
        			}
        			// NMDA input
        			if (P->Connections[j][1] == 2){
        			P->N[P->Exc_Inh_array[i]].x_NMDA_der[0] += P->N[P->Exc_Inh_array[i]].alpha_x_NMDA;
        			}
        			
        			}
        
        
        
        }
    	
    	
    	}
    	/*
    	for (int k=0; k<(P->nNeurons-P->nDopaminergic);k++){
    	if((P->N[P->Exc_Inh_array[k]].Rate*P->dt_neuron/1000)>=gsl_ran_flat(r,0,1)){
    	
    	P->N[P->Exc_Inh_array[k]].x_ext_NMDA_der[0] += P->N[P->Exc_Inh_array[k]].alpha_x_NMDA;//P->N[P->Exc_Inh_array[i]].x_AMPA_der[0] += P->N[P->Exc_Inh_array[i]].alpha_x_AMPA;P->N[P->Exc_Inh_array[i]].s_GABA_der[0] += P->N[P->Exc_Inh_array[i]].alpha_GABA*(1-P->N[P->Exc_Inh_array[i]].s_GABA[0]);
    	P->N[P->Exc_Inh_array[k]].x_ext_AMPA_der[0] += P->N[P->Exc_Inh_array[k]].alpha_x_AMPA;
    	
    	}
    	}*/
    /*--------------------------*/
    	
    /* DOPAMINERGIC NEURONS */
    for (int i=0; i<P->nDopaminergic;i++){
    	
    	P->N[P->Dopaminergic_array[i]].Rate = P->N[P->Dopaminergic_array[i]].ext_rate + 0.5*P->N[P->Dopaminergic_array[i]].amplitude_ext_rate*P->N[P->Dopaminergic_array[i]].Noisy_rate[0]  ;
    	if (P->N[P->Dopaminergic_array[i]].Rate<0.0){P->N[P->Dopaminergic_array[i]].Rate=0.0;}
    	/* DERIVATIVES */
    	
         //Dopaminergic (exc) neurons derivatives and rhs
        P->N[P->Dopaminergic_array[i]].V_der[0] = (-P->N[P->Dopaminergic_array[i]].A_soma*(P->N[P->Dopaminergic_array[i]].I_L + P->N[P->Dopaminergic_array[i]].I_Na + P->N[P->Dopaminergic_array[i]].I_K + P->N[P->Dopaminergic_array[i]].I_A + P->N[P->Dopaminergic_array[i]].I_KS + P->N[P->Dopaminergic_array[i]].I_KNa) - P->N[P->Dopaminergic_array[i]].I_syn_soma - P->N[P->Dopaminergic_array[i]].g_syn_soma*(P->N[P->Dopaminergic_array[i]].V[0]-P->N[P->Dopaminergic_array[i]].V[1]))/(P->N[P->Dopaminergic_array[i]].C*P->N[P->Dopaminergic_array[i]].A_soma);
         
        P->N[P->Dopaminergic_array[i]].V_der[1] = (-P->N[P->Dopaminergic_array[i]].A_den*(P->N[P->Dopaminergic_array[i]].I_Ca + P->N[P->Dopaminergic_array[i]].I_KCa + P->N[P->Dopaminergic_array[i]].I_NaP + P->N[P->Dopaminergic_array[i]].I_AR) - P->N[P->Dopaminergic_array[i]].I_syn_den - P->N[P->Dopaminergic_array[i]].g_syn_den*(P->N[P->Dopaminergic_array[i]].V[1]-P->N[P->Dopaminergic_array[i]].V[0]))/(P->N[P->Dopaminergic_array[i]].C*P->N[P->Dopaminergic_array[i]].A_den);
        
        
    	// Derivatives of gating variables
    	//Exc (soma)
    	P->N[P->Dopaminergic_array[i]].n_K_der[0] = P->N[P->Dopaminergic_array[i]].Temp_fact*(alpha_n_K(&P->N[P->Dopaminergic_array[i]],P->N[P->Dopaminergic_array[i]].V[0])*(1-P->N[P->Dopaminergic_array[i]].n_K[0]) - beta_n_K(&P->N[P->Dopaminergic_array[i]],P->N[P->Dopaminergic_array[i]].V[0])*P->N[P->Dopaminergic_array[i]].n_K[0]);
    	P->N[P->Dopaminergic_array[i]].h_Na_der[0] = P->N[P->Dopaminergic_array[i]].Temp_fact*(alpha_h_Na(&P->N[P->Dopaminergic_array[i]],P->N[P->Dopaminergic_array[i]].V[0])*(1-P->N[P->Dopaminergic_array[i]].h_Na[0]) - beta_h_Na(&P->N[P->Dopaminergic_array[i]],P->N[P->Dopaminergic_array[i]].V[0])*P->N[P->Dopaminergic_array[i]].h_Na[0]);
    	P->N[P->Dopaminergic_array[i]].m_Na_der[0] = P->N[P->Dopaminergic_array[i]].Temp_fact*(alpha_m_Na(&P->N[P->Dopaminergic_array[i]],P->N[P->Dopaminergic_array[i]].V[0])*(1-P->N[P->Dopaminergic_array[i]].m_Na[0]) - beta_m_Na(&P->N[P->Dopaminergic_array[i]],P->N[P->Dopaminergic_array[i]].V[0])*P->N[P->Dopaminergic_array[i]].m_Na[0]);
    	
    	
    	P->N[P->Dopaminergic_array[i]].h_A_der[0] = (1.0/P->N[P->Dopaminergic_array[i]].tau_h_A)*(h_inf_A(&P->N[P->Dopaminergic_array[i]],P->N[P->Dopaminergic_array[i]].V[0])-P->N[P->Dopaminergic_array[i]].h_A[0]);
    	P->N[P->Dopaminergic_array[i]].m_KS_der[0] = (1.0/tau_m_KS(&P->N[P->Dopaminergic_array[i]],P->N[P->Dopaminergic_array[i]].V[0]))*(m_inf_KS(&P->N[P->Dopaminergic_array[i]],P->N[P->Dopaminergic_array[i]].V[0])-P->N[P->Dopaminergic_array[i]].m_KS[0]);
    	
    	
    	// Derivative of the concentration of Ca and Na
    	
    	P->N[P->Dopaminergic_array[i]].Conc_Ca_der[0] = -P->N[P->Dopaminergic_array[i]].alpha_Ca*P->N[P->Dopaminergic_array[i]].A_den*P->N[P->Dopaminergic_array[i]].I_Ca - P->N[P->Dopaminergic_array[i]].Conc_Ca[0]/P->N[P->Dopaminergic_array[i]].tau_Ca;
    	
    	P->N[P->Dopaminergic_array[i]].Conc_Na_der[0] = -P->N[P->Dopaminergic_array[i]].alpha_Na*(P->N[P->Dopaminergic_array[i]].A_soma*P->N[P->Dopaminergic_array[i]].I_Na + P->N[P->Dopaminergic_array[i]].A_den*P->N[P->Dopaminergic_array[i]].I_NaP) - (P->N[P->Dopaminergic_array[i]].R_pump + P->N[P->Dopaminergic_array[i]].D*P->N[P->Dopaminergic_array[i]].Xi[0])*( pow(P->N[P->Dopaminergic_array[i]].Conc_Na[0],3)/(pow(P->N[P->Dopaminergic_array[i]].Conc_Na[0],3)+pow(15.0,3)) - pow(P->N[P->Dopaminergic_array[i]].Conc_Na_eq,3)/(pow(P->N[P->Dopaminergic_array[i]].Conc_Na_eq,3)+pow(15.0,3) ) );
    	
    	// Derivative of time-dependent conductances
    	P->N[P->Dopaminergic_array[i]].s_AMPA_der[0] = P->N[P->Dopaminergic_array[i]].alpha_s_AMPA*P->N[P->Dopaminergic_array[i]].x_AMPA[0]*(1-P->N[P->Dopaminergic_array[i]].s_AMPA[0])-P->N[P->Dopaminergic_array[i]].s_AMPA[0]/P->N[P->Dopaminergic_array[i]].tau_s_AMPA;
    	P->N[P->Dopaminergic_array[i]].s_NMDA_der[0] = P->N[P->Dopaminergic_array[i]].alpha_s_NMDA*P->N[P->Dopaminergic_array[i]].x_NMDA[0]*(1-P->N[P->Dopaminergic_array[i]].s_NMDA[0])-P->N[P->Dopaminergic_array[i]].s_NMDA[0]/P->N[P->Dopaminergic_array[i]].tau_s_NMDA;
    	
    	P->N[P->Dopaminergic_array[i]].s_ext_AMPA_der[0] = P->N[P->Dopaminergic_array[i]].alpha_s_AMPA*P->N[P->Dopaminergic_array[i]].x_ext_AMPA[0]*(1-P->N[P->Dopaminergic_array[i]].s_ext_AMPA[0])-P->N[P->Dopaminergic_array[i]].s_ext_AMPA[0]/P->N[P->Dopaminergic_array[i]].tau_s_AMPA;
    	P->N[P->Dopaminergic_array[i]].s_ext_NMDA_der[0] = P->N[P->Dopaminergic_array[i]].alpha_s_NMDA*P->N[P->Dopaminergic_array[i]].x_ext_NMDA[0]*(1-P->N[P->Dopaminergic_array[i]].s_ext_NMDA[0])-P->N[P->Dopaminergic_array[i]].s_ext_NMDA[0]/P->N[P->Dopaminergic_array[i]].tau_s_NMDA;
    	
    	
    	P->N[P->Dopaminergic_array[i]].s_GABA_der[0] = -P->N[P->Dopaminergic_array[i]].s_GABA[0]/P->N[P->Dopaminergic_array[i]].tau_s_GABA;
    	
    	
    	P->N[P->Dopaminergic_array[i]].x_AMPA_der[0] = -P->N[P->Dopaminergic_array[i]].x_AMPA[0]/P->N[P->Dopaminergic_array[i]].tau_x_AMPA;
    	P->N[P->Dopaminergic_array[i]].x_NMDA_der[0] = -P->N[P->Dopaminergic_array[i]].x_NMDA[0]/P->N[P->Dopaminergic_array[i]].tau_x_NMDA; 
    	
    	P->N[P->Dopaminergic_array[i]].x_ext_AMPA_der[0] = -P->N[P->Dopaminergic_array[i]].x_ext_AMPA[0]/P->N[P->Dopaminergic_array[i]].tau_x_AMPA;
    	P->N[P->Dopaminergic_array[i]].x_ext_NMDA_der[0] = -P->N[P->Dopaminergic_array[i]].x_ext_NMDA[0]/P->N[P->Dopaminergic_array[i]].tau_x_NMDA; 
    	
    	// Gating variable exc 
    	P->N[P->Dopaminergic_array[i]].m_inf_Na = alpha_m_Na(&P->N[P->Dopaminergic_array[i]],P->N[P->Dopaminergic_array[i]].V[0])/(alpha_m_Na(&P->N[P->Dopaminergic_array[i]],P->N[P->Dopaminergic_array[i]].V[0])+beta_m_Na(&P->N[P->Dopaminergic_array[i]],P->N[P->Dopaminergic_array[i]].V[0]));
    	
    	// Currents
    	// Soma
    	P->N[P->Dopaminergic_array[i]].I_K = P->N[P->Dopaminergic_array[i]].g_K*pow(P->N[P->Dopaminergic_array[i]].n_K[0],4)*(P->N[P->Dopaminergic_array[i]].V[0]-P->N[P->Dopaminergic_array[i]].V_K); 
    	
    	//P->N[P->Dopaminergic_array[i]].I_Na = P->N[P->Dopaminergic_array[i]].g_Na*pow(P->N[P->Dopaminergic_array[i]].m_Na[0],3)*P->N[P->Dopaminergic_array[i]].h_Na[0]*(P->N[P->Dopaminergic_array[i]].V[0]-P->N[P->Dopaminergic_array[i]].V_Na);
    	P->N[P->Dopaminergic_array[i]].I_Na = P->N[P->Dopaminergic_array[i]].g_Na*pow(P->N[P->Dopaminergic_array[i]].m_inf_Na,3)*P->N[P->Dopaminergic_array[i]].h_Na[0]*(P->N[P->Dopaminergic_array[i]].V[0]-P->N[P->Dopaminergic_array[i]].V_Na);
    	
    	P->N[P->Dopaminergic_array[i]].I_L = P->N[P->Dopaminergic_array[i]].g_L*(P->N[P->Dopaminergic_array[i]].V[0]-P->N[P->Dopaminergic_array[i]].V_L); 
    	
    	P->N[P->Dopaminergic_array[i]].I_A = P->N[P->Dopaminergic_array[i]].g_A*pow(m_inf_A(&P->N[P->Dopaminergic_array[i]],P->N[P->Dopaminergic_array[i]].V[0]),3)*P->N[P->Dopaminergic_array[i]].h_A[0]*(P->N[P->Dopaminergic_array[i]].V[0]-P->N[P->Dopaminergic_array[i]].V_K);
    	
    	P->N[P->Dopaminergic_array[i]].I_KS = P->N[P->Dopaminergic_array[i]].g_KS*P->N[P->Dopaminergic_array[i]].m_KS[0]*(P->N[P->Dopaminergic_array[i]].V[0]-P->N[P->Dopaminergic_array[i]].V_K);
    	
    	P->N[P->Dopaminergic_array[i]].I_KNa = P->N[P->Dopaminergic_array[i]].g_KNa*w_inf(&P->N[P->Dopaminergic_array[i]],P->N[P->Dopaminergic_array[i]].Conc_Na[0])*(P->N[P->Dopaminergic_array[i]].V[0]-P->N[P->Dopaminergic_array[i]].V_K);
    	
    	 	
    	// Dendrites
    	
    	P->N[P->Dopaminergic_array[i]].I_NaP = P->N[P->Dopaminergic_array[i]].g_NaP*pow(m_inf_NaP(&P->N[P->Dopaminergic_array[i]], P->N[P->Dopaminergic_array[i]].V[1]),3)*(P->N[P->Dopaminergic_array[i]].V[1]-P->N[P->Dopaminergic_array[i]].V_Na);
    	
    	P->N[P->Dopaminergic_array[i]].I_AR = P->N[P->Dopaminergic_array[i]].g_AR*h_inf_AR(&P->N[P->Dopaminergic_array[i]],P->N[P->Dopaminergic_array[i]].V[1])*(P->N[P->Dopaminergic_array[i]].V[1]-P->N[P->Dopaminergic_array[i]].V_K);
    	
    	P->N[P->Dopaminergic_array[i]].I_Ca = P->N[P->Dopaminergic_array[i]].g_Ca*pow(m_inf_Ca(&P->N[P->Dopaminergic_array[i]],P->N[P->Dopaminergic_array[i]].V[1]),2)*(P->N[P->Dopaminergic_array[i]].V[1]-P->N[P->Dopaminergic_array[i]].V_Ca);
    	
    	P->N[P->Dopaminergic_array[i]].I_KCa = P->N[P->Dopaminergic_array[i]].g_KCa*(P->N[P->Dopaminergic_array[i]].Conc_Ca[0]/(P->N[P->Dopaminergic_array[i]].Conc_Ca[0]+P->N[P->Dopaminergic_array[i]].K_D))*(P->N[P->Dopaminergic_array[i]].V[1]-P->N[P->Dopaminergic_array[i]].V_K);
    	
    	
    	
    	//Synaptic
    	P->N[P->Dopaminergic_array[i]].I_syn_soma = P->N[P->Dopaminergic_array[i]].I_GABA;
    	P->N[P->Dopaminergic_array[i]].I_syn_den = P->N[P->Dopaminergic_array[i]].I_AMPA + P->N[P->Dopaminergic_array[i]].I_NMDA + P->N[P->Dopaminergic_array[i]].I_ext_AMPA + P->N[P->Dopaminergic_array[i]].I_ext_NMDA;
    	
    	P->N[P->Dopaminergic_array[i]].I_AMPA = P->N[P->Dopaminergic_array[i]].g_AMPA*P->N[P->Dopaminergic_array[i]].s_AMPA[0]*(P->N[P->Dopaminergic_array[i]].V[1]);
    	P->N[P->Dopaminergic_array[i]].I_NMDA = P->N[P->Dopaminergic_array[i]].g_NMDA*P->N[P->Dopaminergic_array[i]].s_NMDA[0]*(P->N[P->Dopaminergic_array[i]].V[1])*(1.0/(1+P->N[P->Dopaminergic_array[i]].Conc_Mg*exp((-0.062*P->N[P->Dopaminergic_array[i]].V[1]/3.57))));
    	
    	P->N[P->Dopaminergic_array[i]].I_ext_AMPA = P->N[P->Dopaminergic_array[i]].g_ext_AMPA*P->N[P->Dopaminergic_array[i]].s_ext_AMPA[0]*(P->N[P->Dopaminergic_array[i]].V[1]);
    	P->N[P->Dopaminergic_array[i]].I_ext_NMDA = P->N[P->Dopaminergic_array[i]].g_ext_NMDA*P->N[P->Dopaminergic_array[i]].s_ext_NMDA[0]*(P->N[P->Dopaminergic_array[i]].V[1])*(1.0/(1+P->N[P->Dopaminergic_array[i]].Conc_Mg*exp((-0.062*P->N[P->Dopaminergic_array[i]].V[1]/3.57))));
    	
    	
    	
    	/* EXTERNAL TRAIN OF POISSON SPIKES */
    	if((P->N[P->Dopaminergic_array[i]].Rate*P->dt_neuron/1000)>=gsl_ran_flat(r,0,1)){
    	
    	P->N[P->Dopaminergic_array[i]].x_ext_NMDA_der[0] += P->N[P->Dopaminergic_array[i]].alpha_x_NMDA;/*P->N[P->Dopaminergic_array[i]].x_AMPA_der[0] += P->N[P->Dopaminergic_array[i]].alpha_x_AMPA;P->N[P->Dopaminergic_array[i]].s_GABA_der[0] += P->N[P->Dopaminergic_array[i]].alpha_GABA*(1-P->N[P->Dopaminergic_array[i]].s_GABA[0]);*/
    	P->N[P->Dopaminergic_array[i]].x_ext_AMPA_der[0] += P->N[P->Dopaminergic_array[i]].alpha_x_AMPA;
    	}
    	
    	
    	/* DOPAMINERGIC NEURONS INPUTS */
    	
    	int sum=0;
    	
    	for (int a=0;a<i;a++){sum+=P->Number_of_Connections[a-0];}
       
        
        for (int j=sum;j<sum+P->Number_of_Connections[i];j++){
        
       				
        
        			if(P->N[P->Connections[j][0]].V_0[P->Delay_Matrix[P->Connections[j][0]]] >= P->N[P->Connections[j][0]].Vpeak){
        			//if(P->N[P->Connections[j][0]].V[0] >= P->N[P->Connections[j][0]].Vpeak){
        			// AMPA input
        			if (P->Connections[j][1] == 1){
        			P->N[P->Dopaminergic_array[i]].x_AMPA_der[0] += P->N[P->Dopaminergic_array[i]].alpha_x_AMPA;
        			
        			}
        			// GABA input
        			if (P->Connections[j][1] == 0){
        			P->N[P->Dopaminergic_array[i]].s_GABA_der[0] += P->N[P->Dopaminergic_array[i]].alpha_GABA*(1-P->N[P->Dopaminergic_array[i]].s_GABA[0]);
        			
        			}
        			// NMDA input
        			if (P->Connections[j][1] == 2){
        			P->N[P->Dopaminergic_array[i]].x_NMDA_der[0] += P->N[P->Dopaminergic_array[i]].alpha_x_NMDA;
        			}
        			
        			}
        
        
        
        }

    	}	
 
};
	

