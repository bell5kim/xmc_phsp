/*****************************************************************************
 * photo_XS_total.cpp:                                                       *
 *    class member functions for:                                            *
 *       photo_XS_total: total photo-absorption cross section                *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 14.12.1999      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <fstream>
using namespace std;

#include "global.h"
#include "photo_XS_total.h"

/*****************************************************************************
 * member functions for class photo_XS_total:                                *
 *    total photoelectric absorption cross section                           *
 *****************************************************************************/

// delete total photo data
photo_XS_total::~photo_XS_total(void)
{
   delete [] log_sigma; log_sigma = NULL;
   n_energy = 0;
}

// input total photo data from NIST (XCOM) file
void photo_XS_total::nist(char *file_name)
{
   // delete old cross section data
   if (log_sigma != NULL)
   {
      delete [] log_sigma; log_sigma = NULL;
   }

   // open data file
   ifstream file(file_name, ios::in);
   if (!file)
   {
      xvmc_error("photo_XS_total::nist","cannot open file",8);
   }

   // number of energy entries in file
   unsigned num_inp = 0;
   file >> num_inp;

   // mass density
   real density = ZERO;
   file >> density;

   // electron density
   real electron_density = ZERO;
   file >> electron_density;

   // allocate memory for energy input data
   real *log_energy_inp = NULL;
   if ( (log_energy_inp  = new real[num_inp]) == NULL )
   {
      xvmc_error("photo_XS_total::nist",
                 "cannot allocate memory for photon energy input",8);
   }

   // allocate memory for total photo cross section input data
   real *log_mu_phot_inp = NULL;
   if ( (log_mu_phot_inp  = new real[num_inp]) == NULL )
   {
      xvmc_error("photo_XS_total::nist",
                 "cannot allocate memory for photo cross section input",8);
   }

   // read data from file
   real energy_inp  = ZERO;
   real mu_phot_inp = ZERO;
   real dummy       = ZERO;
   for (register unsigned int i=0; i<num_inp; ++i)
   {
      file >> energy_inp
           >> dummy            // coherent scattering
           >> dummy            // Compton  scattering
           >> mu_phot_inp      // photoelectric absorption
           >> dummy            // pair production in nuclear field
           >> dummy            // pair production in electron field
           >> dummy            // total attenuation with coherent scatter
           >> dummy;           // total attenuation without coherent scatter

      // log energy scale
      log_energy_inp[i]  = log(energy_inp);

      // mass -> linear attenuation (log scale)
      log_mu_phot_inp[i] = log(mu_phot_inp*density);
   }

   // close file
   file.close();

   // set number of photon energy bins
   n_energy = 2000;

   // minimum energy
   log_energy_min = log_energy_inp[0];
   energy_min     = exp(log_energy_min);

   // maximum energy
   log_energy_max = log_energy_inp[num_inp-1];
   energy_max     = exp(log_energy_max);
   if (energy_max <= energy_min)
   {
      xvmc_error("photo_XS_total::nist",
                 "maximum energy in input file too small",8);
   }

   // calculate inverse bin size for photon energy (log scale)
   inv_delta = float(n_energy-1)/(log_energy_max - log_energy_min);

   // bin size
   real delta = ONE/inv_delta;

   // allocate memory for photo cross section data
   if ( (log_sigma = new real[n_energy]) == NULL )
   {
      xvmc_error("photo_XS_total::nist",
                 "cannot allocate memory for photo cross section data",8);
   }

   // interpolate total photo cross section data
   real log_energy = log_energy_min;
   for (register int j=0; j<n_energy; ++j)
   {
      // input array indices
      int lower,upper;
      get_index(num_inp,log_energy_inp,log_energy,lower,upper);

      // interpolate
      real delta_inp = log_energy_inp[upper] - log_energy_inp[lower];
      real p,q;
      if (delta_inp > ZERO)
      {
         p = (log_energy - log_energy_inp[lower])/delta_inp;
         q = ONE-p;
      }
      else
      {
         p = ONE;
         q = ZERO;
      }

      // calculate logarithm and store
      log_sigma[j] = q*log_mu_phot_inp[lower] + p*log_mu_phot_inp[upper];

      // next energy
      log_energy += delta;
   }

   // delete input arrays
   delete [] log_energy_inp;  log_energy_inp  = NULL;
   delete [] log_mu_phot_inp; log_mu_phot_inp = NULL;
}

// get total photo cross section for the specified energy
real photo_XS_total::get(real energy)
{
   real  log_energy;     // ln(energy)
   int   i_bin;          // bin number, integer value
   real  r_bin;          // bin number, real value
   real  sigma;          // cross section

   if (energy <= energy_min) return(exp(log_sigma[0]));

   if (energy >= energy_max) return(exp(log_sigma[n_energy-1]));

   // energy_min < energy < energy_max
   log_energy = log(energy);
   r_bin      = (log_energy-log_energy_min)*inv_delta;
   i_bin      = int(r_bin);
   r_bin      = r_bin - float(i_bin);
   sigma      = exp( (ONE-r_bin)*log_sigma[i_bin] + r_bin*log_sigma[i_bin+1] );
   return(sigma);
}
