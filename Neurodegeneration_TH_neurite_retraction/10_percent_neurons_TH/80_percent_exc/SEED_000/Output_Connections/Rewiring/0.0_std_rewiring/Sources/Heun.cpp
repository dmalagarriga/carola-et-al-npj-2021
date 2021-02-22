#include "../Headers/Heun.h"
#include "../Headers/rhs.h"
#include "../Headers/neuron.h"
#include "../Headers/Parameters.h"

#include <math.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"


void Heun (Parameters *P){
	
		for (int i=0; i < P->nNeurons; i++){
		
		P->N[i].Xi_old[0]=P->N[i].Xi[0];
		P->N[i].Noisy_rate_old[0]=P->N[i].Noisy_rate[0];
		
		
			
											 }	
											 
	
	rhs(P); // With this we have all derivatives
	
		for (int iV=0; iV < P->nNeurons; iV++){	
			
			P->N[iV].V_old[0]=P->N[iV].V[0];
			P->N[iV].V_old[1]=P->N[iV].V[1];
			P->N[iV].n_K_old[0] = P->N[iV].n_K[0];
			P->N[iV].h_Na_old[0] = P->N[iV].h_Na[0];
			P->N[iV].m_Na_old[0] = P->N[iV].m_Na[0];
			P->N[iV].h_A_old[0] = P->N[iV].h_A[0];
			P->N[iV].m_KS_old[0] = P->N[iV].m_KS[0];
			P->N[iV].Conc_Ca_old[0] = P->N[iV].Conc_Ca[0];
			P->N[iV].Conc_Na_old[0] = P->N[iV].Conc_Na[0];
			P->N[iV].s_AMPA_old[0] = P->N[iV].s_AMPA[0];
			P->N[iV].s_NMDA_old[0] = P->N[iV].s_NMDA[0];
			P->N[iV].s_ext_AMPA_old[0] = P->N[iV].s_ext_AMPA[0];
			P->N[iV].s_ext_NMDA_old[0] = P->N[iV].s_ext_NMDA[0];
			P->N[iV].s_GABA_old[0] = P->N[iV].s_GABA[0];
			P->N[iV].x_AMPA_old[0] = P->N[iV].x_AMPA[0];
			P->N[iV].x_NMDA_old[0] = P->N[iV].x_NMDA[0];
			P->N[iV].x_ext_AMPA_old[0] = P->N[iV].x_ext_AMPA[0];
			P->N[iV].x_ext_NMDA_old[0] = P->N[iV].x_ext_NMDA[0];
			
			P->N[iV].V_der_old[0]=P->N[iV].V_der[0];
			P->N[iV].V_der_old[1]=P->N[iV].V_der[1];
			P->N[iV].n_K_der_old[0] =P->N[iV].n_K_der[0];
			P->N[iV].m_Na_der_old[0] = P->N[iV].m_Na_der[0];
			P->N[iV].h_Na_der_old[0] = P->N[iV].h_Na_der[0];
			P->N[iV].h_A_der_old[0] =P->N[iV].h_A_der[0];
			P->N[iV].m_KS_der_old[0] =P->N[iV].m_KS_der[0];
			P->N[iV].Conc_Ca_der_old[0] = P->N[iV].Conc_Ca_der[0];
			P->N[iV].Conc_Na_der_old[0] = P->N[iV].Conc_Na_der[0];
			P->N[iV].s_AMPA_der_old[0] = P->N[iV].s_AMPA_der[0];
			P->N[iV].s_NMDA_der_old[0] = P->N[iV].s_NMDA_der[0];
			P->N[iV].s_ext_AMPA_der_old[0] = P->N[iV].s_ext_AMPA_der[0];
			P->N[iV].s_ext_NMDA_der_old[0] = P->N[iV].s_ext_NMDA_der[0];
			P->N[iV].s_GABA_der_old[0] = P->N[iV].s_GABA_der[0];
			P->N[iV].x_AMPA_der_old[0] = P->N[iV].x_AMPA_der[0];
			P->N[iV].x_NMDA_der_old[0] = P->N[iV].x_NMDA_der[0];
			P->N[iV].x_ext_AMPA_der_old[0] = P->N[iV].x_ext_AMPA_der[0];
			P->N[iV].x_ext_NMDA_der_old[0] = P->N[iV].x_ext_NMDA_der[0];
				
			P->N[iV].V[0] += P->N[iV].V_der[0]*P->dt_neuron;
			P->N[iV].V[1] += P->N[iV].V_der[1]*P->dt_neuron;
			P->N[iV].n_K[0] +=P->N[iV].n_K_der[0]*P->dt_neuron;
			P->N[iV].h_Na[0]+=P->N[iV].h_Na_der[0]*P->dt_neuron;
			P->N[iV].m_Na[0]+=P->N[iV].m_Na_der[0]*P->dt_neuron;
			P->N[iV].h_A[0]+=P->N[iV].h_A_der[0]*P->dt_neuron;
			P->N[iV].m_KS[0]+=P->N[iV].m_KS_der[0]*P->dt_neuron;
			P->N[iV].Conc_Ca[0]+=P->N[iV].Conc_Ca_der[0]*P->dt_neuron;
			P->N[iV].Conc_Na[0]+=P->N[iV].Conc_Na_der[0]*P->dt_neuron;
			
			P->N[iV].s_AMPA[0]+=P->N[iV].s_AMPA_der[0]*P->dt_neuron;
			P->N[iV].s_NMDA[0]+=P->N[iV].s_NMDA_der[0]*P->dt_neuron;
			P->N[iV].s_ext_AMPA[0]+=P->N[iV].s_ext_AMPA_der[0]*P->dt_neuron;
			P->N[iV].s_ext_NMDA[0]+=P->N[iV].s_ext_NMDA_der[0]*P->dt_neuron;
			P->N[iV].s_GABA[0]+=P->N[iV].s_GABA_der[0]*P->dt_neuron;
			P->N[iV].x_AMPA[0]+=P->N[iV].x_AMPA_der[0]*P->dt_neuron;
			P->N[iV].x_NMDA[0]+=P->N[iV].x_NMDA_der[0]*P->dt_neuron;
			P->N[iV].x_ext_AMPA[0]+=P->N[iV].x_ext_AMPA_der[0]*P->dt_neuron;
			P->N[iV].x_ext_NMDA[0]+=P->N[iV].x_ext_NMDA_der[0]*P->dt_neuron;
			
			//P->N[iV].Conc_Na[0]+=P->N[iV].D*P->N[iV].Xi[0]*( pow(P->N[iV].Conc_Na[0],3)/(pow(P->N[iV].Conc_Na[0],3)+pow(15,3)) - pow(P->N[iV].Conc_Na_eq,3)/(pow(P->N[iV].Conc_Na_eq,3)+pow(15,3) ) );
			}
	
	rhs(P);	// This is f_tild
	
	for (int iV=0; iV < P->nNeurons; iV++){
	
			P->N[iV].V[0]=P->N[iV].V_old[0]+0.5*(P->N[iV].V_der_old[0]+P->N[iV].V_der[0])*P->dt_neuron;
			P->N[iV].V[1]=P->N[iV].V_old[1]+0.5*(P->N[iV].V_der_old[1]+P->N[iV].V_der[1])*P->dt_neuron;
			P->N[iV].n_K[0]=P->N[iV].n_K_old[0]+0.5*(P->N[iV].n_K_der_old[0]+P->N[iV].n_K_der[0])*P->dt_neuron;
			P->N[iV].h_Na[0]=P->N[iV].h_Na_old[0]+0.5*(P->N[iV].h_Na_der_old[0]+P->N[iV].h_Na_der[0])*P->dt_neuron;
			P->N[iV].m_Na[0]=P->N[iV].m_Na_old[0]+0.5*(P->N[iV].m_Na_der_old[0]+P->N[iV].m_Na_der[0])*P->dt_neuron;
			P->N[iV].h_A[0]=P->N[iV].h_A_old[0] + 0.5*(P->N[iV].h_A_der_old[0]+P->N[iV].h_A_der[0])*P->dt_neuron;
			P->N[iV].m_KS[0]=P->N[iV].m_KS_old[0] + 0.5*(P->N[iV].m_KS_der_old[0]+P->N[iV].m_KS_der[0])*P->dt_neuron;
			P->N[iV].Conc_Ca[0]=P->N[iV].Conc_Ca_old[0] + 0.5*(P->N[iV].Conc_Ca_der_old[0]+P->N[iV].Conc_Ca_der[0])*P->dt_neuron;
			P->N[iV].Conc_Na[0]=P->N[iV].Conc_Na_old[0] + 0.5*(P->N[iV].Conc_Na_der_old[0]+P->N[iV].Conc_Na_der[0])*P->dt_neuron;
			
			P->N[iV].s_AMPA[0]=P->N[iV].s_AMPA_old[0] + 0.5*(P->N[iV].s_AMPA_der_old[0]+P->N[iV].s_AMPA_der[0])*P->dt_neuron;
			P->N[iV].s_NMDA[0]=P->N[iV].s_NMDA_old[0] + 0.5*(P->N[iV].s_NMDA_der_old[0]+P->N[iV].s_NMDA_der[0])*P->dt_neuron;
			P->N[iV].s_ext_AMPA[0]=P->N[iV].s_ext_AMPA_old[0] + 0.5*(P->N[iV].s_ext_AMPA_der_old[0]+P->N[iV].s_ext_AMPA_der[0])*P->dt_neuron;
			P->N[iV].s_ext_NMDA[0]=P->N[iV].s_ext_NMDA_old[0] + 0.5*(P->N[iV].s_ext_NMDA_der_old[0]+P->N[iV].s_ext_NMDA_der[0])*P->dt_neuron;
			P->N[iV].s_GABA[0]=P->N[iV].s_GABA_old[0] + 0.5*(P->N[iV].s_GABA_der_old[0]+P->N[iV].s_GABA_der[0])*P->dt_neuron;
			P->N[iV].x_AMPA[0]=P->N[iV].x_AMPA_old[0] + 0.5*(P->N[iV].x_AMPA_der_old[0]+P->N[iV].x_AMPA_der[0])*P->dt_neuron;
			P->N[iV].x_NMDA[0]=P->N[iV].x_NMDA_old[0] + 0.5*(P->N[iV].x_NMDA_der_old[0]+P->N[iV].x_NMDA_der[0])*P->dt_neuron;
			P->N[iV].x_ext_AMPA[0]=P->N[iV].x_ext_AMPA_old[0] + 0.5*(P->N[iV].x_ext_AMPA_der_old[0]+P->N[iV].x_ext_AMPA_der[0])*P->dt_neuron;
			P->N[iV].x_ext_NMDA[0]=P->N[iV].x_ext_NMDA_old[0] + 0.5*(P->N[iV].x_ext_NMDA_der_old[0]+P->N[iV].x_ext_NMDA_der[0])*P->dt_neuron;
			
			/*P->N[iV].s_AMPA[0]+=P->N[iV].s_AMPA_der[0]*P->dt_neuron;
			P->N[iV].s_NMDA[0]+=P->N[iV].s_NMDA_der[0]*P->dt_neuron;
			P->N[iV].s_GABA[0]+=P->N[iV].s_GABA_der[0]*P->dt_neuron;
			P->N[iV].x_AMPA[0]+=P->N[iV].x_AMPA_der[0]*P->dt_neuron;
			P->N[iV].x_NMDA[0]+=P->N[iV].x_NMDA_der[0]*P->dt_neuron;*/
			//P->N[iV].Conc_Na[0]+=P->N[iV].D*P->N[iV].Xi[0]*( pow(P->N[iV].Conc_Na[0],3)/(pow(P->N[iV].Conc_Na[0],3)+pow(15,3)) - pow(P->N[iV].Conc_Na_eq,3)/(pow(P->N[iV].Conc_Na_eq,3)+pow(15,3) ) );
	
	}

//P->N[0].Noisy_rate[0]=1.0*(P->N[0].Noisy_rate_old[0]*exp(-P->dt_neuron/(P->N[0].tau_ext_rate))+sqrt((1-exp(-(2*P->dt_neuron)/(P->N[0].tau_ext_rate)))/(2*P->N[0].tau_ext_rate))*gsl_ran_gaussian(r,1.0));
//P->N[0].Xi[0]=1.0*(P->N[0].Xi_old[0]*exp(-P->dt_neuron/(P->N[0].epsilon))+sqrt((1-exp(-(2*P->dt_neuron)/(P->N[0].epsilon)))/(2*P->N[0].epsilon))*gsl_ran_gaussian(r,1.0));
for (int i=0;i<P->nNeurons;i++){
		/* Noise in the channel */
		//P->N[i].Xi[0]=P->N[0].Xi[0];
		//P->N[i].Noisy_rate[0]=P->N[0].Noisy_rate[0];
		P->N[i].Xi[0]=1.0*(P->N[i].Xi_old[0]*exp(-P->dt_neuron/(P->N[i].epsilon))+sqrt((1-exp(-(2*P->dt_neuron)/(P->N[i].epsilon)))/(2*P->N[i].epsilon))*gsl_ran_gaussian(r,1.0));
		P->N[i].Noisy_rate[0]=1.0*(P->N[i].Noisy_rate_old[0]*exp(-P->dt_neuron/(P->N[i].tau_ext_rate))+sqrt((1-exp(-(2*P->dt_neuron)/(P->N[i].tau_ext_rate)))/(2*P->N[i].tau_ext_rate))*gsl_ran_gaussian(r,1.0));
		//P->N[i].Xi[0]=2.0*(P->N[i].Xi_old[0]*exp(-P->dt_neuron/(P->N[i].epsilon))+sqrt((1-exp(-(2*P->dt_neuron)/(P->N[i].epsilon)))/(2*P->N[i].epsilon))*gsl_ran_flat(r,0.0,1.0));
		}
	
};
