/*****************************************************************************
 * photon_data.cpp:                                                          *
 *    class member functions for:                                            *
 *       photon_data:    photon transport data (function of energy)          *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 00/02/22        *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include "definitions.h"
#include "global.h"
#include "photon_data.h"

/*****************************************************************************
 * member functions for photon_data:                                         *
 *****************************************************************************/

// use (copy) transport data from input array to calculate the array
// of transport parameters for MC simulation (better resolution)

extern compton_XS_total tot_comp; // total compton cross section for water
extern pair_XS_total tot_pair;    // total pair cross section data for water
void photon_data_init(photon_data &A, unsigned int number, photon_data_inp *inp, real energy_max)
{
   real ee,delta,p,q;
   int  lower,upper;
   if (energy_max < p_cut)
     xvmc_error("photon_data::photon_data","energy_max < p_cut",8);

   // allocate memory for MC simulation photon data
   A.num       = number;
   A.delta_e   = (1.01*energy_max - p_cut)/double(A.num);
   A.inverse_delta_e = ONE/A.delta_e;
   bool error = false;
   if ( (A.energy    = new real[A.num]) == NULL ) error = true;
   if ( (A.mu_comp   = new real[A.num]) == NULL ) error = true;
   if ( (A.mu_pair   = new real[A.num]) == NULL ) error = true;
   if ( (A.mu_phot   = new real[A.num]) == NULL ) error = true;
   if ( (A.mu_tot    = new real[A.num]) == NULL ) error = true;
   if ( (A.mu_en     = new real[A.num]) == NULL ) error = true;

   if (error)
   {
      xvmc_error("photon_data::photon_data",
                 "cannot allocate memory for photon transport data",8);
   }

   // calculate interpolated parameters
   ee      = p_cut - A.delta_e/TWO;
   for (register unsigned int i=0; i<A.num; ++i)
   {
      ee = ee + A.delta_e;
      get_index(inp->num,inp->energy,ee,lower,upper);
      delta = inp->energy[upper] - inp->energy[lower];
      if (delta > ZERO)
      {
         p = (ee - inp->energy[lower])/delta;
         q = ONE-p;
      }
      else
      {
         p = ONE;
         q = ZERO;
      }
      A.energy[i]  = ee;
      // compton cross section is calculated according to Klein and Nishina
      A.mu_comp[i] = tot_comp.get(ee);
      // pair cross section is calculated according to Bethe and Heitler
      A.mu_pair[i] = tot_pair.get(ee);
      // photo cross section is taken from ICRU 46
      A.mu_phot[i] = q*inp->mu_phot[lower] + p*inp->mu_phot[upper];
      // total cross section is the sum of 3 contributions
      A.mu_tot[i]  = A.mu_comp[i] + A.mu_pair[i] + A.mu_phot[i];
      // energy absorption coefficient is taken from ICRU 46
      A.mu_en[i]   = q*inp->mu_en[lower]   + p*inp->mu_en[upper];
   }
}

photon_data::photon_data(int nTableSize,
			 photon_data_inp& inp,
			 real energy_max,
			 real p_cut,
			 const compton_XS_total& tot_comp, // total compton cross section for water
			 const pair_XS_total& tot_pair) :
  num(nTableSize)
{
  real ee,delta,p,q;
  int  lower,upper;
  if (energy_max < p_cut)
    xvmc_error("photon_data::photon_data","energy_max < p_cut",8);

  // allocate memory for MC simulation photon data
  delta_e   = (1.01*energy_max - p_cut)/double(num);
  inverse_delta_e = ONE/delta_e;
  bool error = false;
  if ( (energy    = new real[num]) == NULL ) error = true;
  if ( (mu_comp   = new real[num]) == NULL ) error = true;
  if ( (mu_pair   = new real[num]) == NULL ) error = true;
  if ( (mu_phot   = new real[num]) == NULL ) error = true;
  if ( (mu_tot    = new real[num]) == NULL ) error = true;
  if ( (mu_en     = new real[num]) == NULL ) error = true;

  if (error)
     {
       xvmc_error("photon_data::photon_data",
		  "cannot allocate memory for photon transport data",8);
     }

  // calculate interpolated parameters
  ee      = p_cut - delta_e/TWO;
  for (register unsigned int i=0; i<num; ++i)
    {
      ee = ee + delta_e;
      get_index(inp.num,inp.energy,ee,lower,upper);
      delta = inp.energy[upper] - inp.energy[lower];
      if (delta > ZERO)
	{
	  p = (ee - inp.energy[lower])/delta;
	  q = ONE-p;
	}
      else
	{
	  p = ONE;
	  q = ZERO;
      }
      energy[i]  = ee;
      // compton cross section is calculated according to Klein and Nishina
      mu_comp[i] = tot_comp.get(ee);
      // pair cross section is calculated according to Bethe and Heitler
      mu_pair[i] = tot_pair.get(ee);
      // photo cross section is taken from ICRU 46
      mu_phot[i] = q*inp.mu_phot[lower] + p*inp.mu_phot[upper];
      // total cross section is the sum of 3 contributions
      mu_tot[i]  = mu_comp[i] + mu_pair[i] + mu_phot[i];
      // energy absorption coefficient is taken from ICRU 46
      mu_en[i]   = q*inp.mu_en[lower]   + p*inp.mu_en[upper];
    }
}

// delete data
photon_data::~photon_data()
{
   num = 0;
   delta_e = ZERO;
   inverse_delta_e = ZERO;
   if(energy!=NULL)  delete [] energy;
   if(mu_comp!=NULL) delete [] mu_comp;
   if(mu_pair!=NULL) delete [] mu_pair;
   if(mu_phot!=NULL) delete [] mu_phot;
   if(mu_tot!=NULL)  delete [] mu_tot;
   if(mu_en!=NULL)   delete [] mu_en;
}
