#include <iostream>
#include <cstdlib>
#include <numeric>
#include <algorithm>
#include <vector>
#include "../Headers/Parameters.h"
#include "../Headers/Heun.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

using namespace std;

Parameters *P;
const gsl_rng_type * T;
gsl_rng * r;
unsigned long int seed=atof(getenv("SEED"));


int main (int argc, char **argv){

	P=new Parameters;
	gsl_rng_env_setup();
	T = gsl_rng_default;
    r = gsl_rng_alloc (T);
   
// Initial conditions
    
	gsl_rng_set(r,seed); 
	
	
	
	//Initial conditions for all the neurons


	for (int i=0; i<P->nNeurons; i++){
		
		
		
			
			P->N[i].V[1] =  -60.0+(gsl_ran_flat(r,0,1)-0.5)*5.0;
			P->N[i].Conc_Na[0]	= P->N[i].Conc_Na_eq+(gsl_ran_flat(r,0,1)-0.5)*2.0;
			P->N[i].Conc_Ca[0]	= 0.5;		
 			P->N[i].n_K[0] = 0.0;
 			P->N[i].m_Na[0] = 0.0;
 			P->N[i].h_Na[0] = 0.0;
 			P->N[i].Xi[0] = sqrt(1.0/P->N[i].epsilon)*gsl_ran_gaussian(r,1.0);
 			P->N[i].Noisy_rate[0] = sqrt(1.0/P->N[i].tau_ext_rate);
 			
 			P->N[i].alpha_Na = P->N[i].alpha_Na + 2.0*gsl_ran_gaussian(r,1);
			P->N[i].R_pump = P->N[i].R_pump + 1.8e-3*gsl_ran_gaussian(r,1);
			P->N[i].g_syn_den = P->N[i].g_syn_den + gsl_ran_gaussian(r,1)*0.1e-3;
			P->N[i].g_syn_soma = P->N[i].g_syn_den + gsl_ran_gaussian(r,1)*0.1e-3;
			
		}

    for (int i=0; i<(P->nNeurons-P->nDopaminergic)-(1.0-P->Prop_exc_inh)*(P->nNeurons-P->nDopaminergic);i++){
			
			P->N[P->Exc_Inh_array[i]].V[0] = -86.0+(gsl_ran_flat(r,0,1)-0.5)*5.0;
			
			P->N[P->Exc_Inh_array[i]].g_AMPA = P->N[P->Exc_Inh_array[i]].g_AMPA_exc;
			P->N[P->Exc_Inh_array[i]].g_NMDA = P->N[P->Exc_Inh_array[i]].g_NMDA_exc;
			P->N[P->Exc_Inh_array[i]].g_ext_AMPA = P->N[P->Exc_Inh_array[i]].g_ext_AMPA_exc;
			P->N[P->Exc_Inh_array[i]].g_ext_NMDA = P->N[P->Exc_Inh_array[i]].g_ext_NMDA_exc;
			P->N[P->Exc_Inh_array[i]].g_GABA = P->N[P->Exc_Inh_array[i]].g_GABA_inh;
			P->N[P->Exc_Inh_array[i]].V_L = P->N[P->Exc_Inh_array[i]].V_L + 0.3*gsl_ran_gaussian(r,1);
			P->N[P->Exc_Inh_array[i]].g_L = P->N[P->Exc_Inh_array[i]].g_L + 6.7e-3*gsl_ran_gaussian(r,1);
			P->N[P->Exc_Inh_array[i]].V_K = (8.314*(34+273.15)/96.845)*log(2.45/150);
			P->N[P->Exc_Inh_array[i]].Vpeak = 0.0;
			P->N[P->Exc_Inh_array[i]].Temp_fact = 4.0;
			
    		
    }

    for (int i=(P->nNeurons-P->nDopaminergic)-(1.0-P->Prop_exc_inh)*(P->nNeurons-P->nDopaminergic); i<=(P->nNeurons-P->nDopaminergic);i++){
    		
			P->N[P->Exc_Inh_array[i]].V[0] = -60.0+(gsl_ran_flat(r,0,1)-0.5)*5.0;
	
			P->N[P->Exc_Inh_array[i]].g_AMPA = P->N[P->Exc_Inh_array[i]].g_AMPA_inh;
			P->N[P->Exc_Inh_array[i]].g_NMDA = P->N[P->Exc_Inh_array[i]].g_NMDA_inh;
			P->N[P->Exc_Inh_array[i]].g_ext_AMPA = P->N[P->Exc_Inh_array[i]].g_ext_AMPA_inh;
			P->N[P->Exc_Inh_array[i]].g_ext_NMDA = P->N[P->Exc_Inh_array[i]].g_ext_NMDA_inh;
			P->N[P->Exc_Inh_array[i]].g_GABA = P->N[P->Exc_Inh_array[i]].g_GABA_inh;
			P->N[P->Exc_Inh_array[i]].V_L_i = P->N[P->Exc_Inh_array[i]].V_L_i + 0.15*gsl_ran_gaussian(r,1);
			P->N[P->Exc_Inh_array[i]].g_L_i = P->N[P->Exc_Inh_array[i]].g_L_i + 0.0025*gsl_ran_gaussian(r,1);
			P->N[P->Exc_Inh_array[i]].V_K_i = (8.314*(34+273.15)/96.485)*log(2.45/150);
			P->N[P->Exc_Inh_array[i]].Vpeak = -20.0;
			
    }
    
    for (int i=0; i<P->nDopaminergic;i++){
    		
			P->N[P->Dopaminergic_array[i]].V[0] = -86.0+(gsl_ran_flat(r,0,1)-0.5)*5.0;
			
			P->N[P->Dopaminergic_array[i]].g_AMPA = P->N[P->Dopaminergic_array[i]].g_AMPA_exc;
			P->N[P->Dopaminergic_array[i]].g_NMDA = P->N[P->Dopaminergic_array[i]].g_NMDA_exc;
			P->N[P->Dopaminergic_array[i]].g_ext_AMPA = P->N[P->Dopaminergic_array[i]].g_ext_AMPA_exc;
			P->N[P->Dopaminergic_array[i]].g_ext_NMDA = P->N[P->Dopaminergic_array[i]].g_ext_NMDA_exc;
			P->N[P->Dopaminergic_array[i]].g_GABA = P->N[P->Dopaminergic_array[i]].g_GABA_inh;
			P->N[P->Dopaminergic_array[i]].V_L = P->N[P->Dopaminergic_array[i]].V_L + 0.3*gsl_ran_gaussian(r,1);
			P->N[P->Dopaminergic_array[i]].g_L = P->N[P->Dopaminergic_array[i]].g_L + 6.7e-3*gsl_ran_gaussian(r,1);
			P->N[P->Dopaminergic_array[i]].V_K = (8.314*(34+273.15)/96.485)*log(2.45/150);
			P->N[P->Dopaminergic_array[i]].Vpeak = 0.0;
			P->N[P->Dopaminergic_array[i]].Temp_fact = 4.0;
    }
	
	// The integration loop is here 
 	
      for (P->it=0; (P->it) < (P->nt); P->it++){		
	  
			
 			
            Heun(P);
          	
	 		 		
	    	P->RefreshDelay();
	   		P->t_neuron+=P->dt_neuron;
	   	
		
			if (P->it % P->ntout ==0)
				{
		 		P->WriteData(); 	
				}
											
					
        
	}
	
	
	gsl_rng_free (r);
	delete P;

return EXIT_SUCCESS;
}
