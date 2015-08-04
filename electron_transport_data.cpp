/*****************************************************************************
 * electron_transport_data.cpp:                                              *
 *    class member functions for:                                            *
 *       electron_transport_data: electron stopping powers and ranges        *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 06.12.2001      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <fstream>
using namespace std;

#include "global.h"
#include "electron_transport_data.h"

/*****************************************************************************
 * member functions for class electron_transport_data:                       *
 *    electron stopping powers and ranges                                    *
 *****************************************************************************/

// delete electron transport data data
electron_transport_data::~electron_transport_data(void)
{
   delete [] log_dedx; log_dedx = NULL;
   delete [] log_csda; log_csda = NULL;
   n_energy = 0;
}

// input electron transport data from ESTAR (NIST, ICRU) file
void electron_transport_data::estar(char *file_name)
{
   // delete old electron transport data
   if (log_dedx != NULL)
   {
      delete [] log_dedx; log_dedx = NULL;
   }
   if (log_csda != NULL)
   {
      delete [] log_csda; log_csda = NULL;
   }

   // open data file
   ifstream file(file_name, ios::in);
   if (!file)
   {
      xvmc_error("electron_transport_data::estar","cannot open file",8);
   }

   // number of energy entries in file
   unsigned num_inp = 0;
   file >> num_inp;

   // mass density
   real density = ZERO;
   file >> density;

   // allocate memory for energy input data
   real *log_energy_inp = NULL;
   if ( (log_energy_inp  = new real[num_inp]) == NULL )
   {
      xvmc_error("electron_transport_data::estar",
                 "cannot allocate memory for electron energy input",8);
   }

   // allocate memory for total stopping power input data
   real *log_dedx_inp = NULL;
   if ( (log_dedx_inp  = new real[num_inp]) == NULL )
   {
      xvmc_error("electron_transport_data::estar",
                 "cannot allocate memory for total stopping power input",8);
   }

   // allocate memory for CSDA range input data
   real *log_csda_inp = NULL;
   if ( (log_csda_inp  = new real[num_inp]) == NULL )
   {
      xvmc_error("electron_transport_data::estar",
                 "cannot allocate memory for CSDA range input",8);
   }

   // read data from file
   real energy_inp  = ZERO;
   real dedx_inp    = ZERO;
   real csda_inp    = ZERO;
   real dummy       = ZERO;
   for (register unsigned int i=0; i<num_inp; ++i)
   {
      file >> energy_inp
           >> dummy            // collision stopping power (MeV cm2/g)
           >> dummy            // radiative stopping power (MeV cm2/g)
           >> dedx_inp         // total     stopping power (MeV cm2/g)
           >> csda_inp         // CSDA range (g/cm2)
           >> dummy            // radiation yield
           >> dummy;           // density effect parameter

      // log scale
      log_energy_inp[i] = log(energy_inp);

      // mass stopping power (MeV cm2/g) -> linear stopping power (MeV/cm)
      // -> log scale
      log_dedx_inp[i] = log(dedx_inp*density);

      // density*range (g/cm2) -> range (cm) -> log scale
      log_csda_inp[i] = log(csda_inp/density);
   }

   // close file
   file.close();

   // set number of electron energy bins
   n_energy = 2000;

   // minimum energy
   log_energy_min = log_energy_inp[0];
   energy_min     = exp(log_energy_min);

   // maximum energy
   log_energy_max = log_energy_inp[num_inp-1];
   energy_max     = exp(log_energy_max);
   if (energy_max <= energy_min)
   {
      xvmc_error("electron_transport_data::estar",
                 "maximum energy in input file too small",8);
   }

   // calculate inverse bin size for electron energy (log scale)
   inv_delta = float(n_energy-1)/(log_energy_max - log_energy_min);

   // bin size
   real delta = ONE/inv_delta;

   // allocate memory for total stopping power data
   if ( (log_dedx = new real[n_energy]) == NULL )
   {
      xvmc_error("electron_transport_data::estar",
                 "cannot allocate memory for total stopping power data",8);
   }

   // allocate memory for CSDA range data
   if ( (log_csda = new real[n_energy]) == NULL )
   {
      xvmc_error("electron_transport_data::estar",
                 "cannot allocate memory for CSDA range data",8);
   }

   // interpolate electron transport data
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
      log_dedx[j] = q*log_dedx_inp[lower] + p*log_dedx_inp[upper];
      log_csda[j] = q*log_csda_inp[lower] + p*log_csda_inp[upper];

      // next energy
      log_energy += delta;
   }

   // delete input arrays
   delete [] log_energy_inp;  log_energy_inp  = NULL;
   delete [] log_dedx_inp;    log_dedx_inp    = NULL;
   delete [] log_csda_inp;    log_csda_inp    = NULL;

   return;
}

// get total stopping power for the specified energy
real electron_transport_data::get_dedx(real energy)
{
   real  log_energy;     // ln(energy)
   int   i_bin;          // bin number, integer value
   real  r_bin;          // bin number, real value
   real  dedx;           // total stopping power

   if (energy <= energy_min) return(exp(log_dedx[0]));

   if (energy >= energy_max) return(exp(log_dedx[n_energy-1]));

   // energy_min < energy < energy_max
   log_energy = log(energy);
   r_bin      = (log_energy-log_energy_min)*inv_delta;
   i_bin      = int(r_bin);
   r_bin      = r_bin - float(i_bin);
   dedx       = exp( (ONE-r_bin)*log_dedx[i_bin] + r_bin*log_dedx[i_bin+1] );
   return(dedx);
}

// get continuous slowing down (CSDA) range for the specified energy
real electron_transport_data::get_csda(real energy)
{
   real  log_energy;     // ln(energy)
   int   i_bin;          // bin number, integer value
   real  r_bin;          // bin number, real value
   real  csda;           // CSDA range

   if (energy <= energy_min) return(exp(log_csda[0]));

   if (energy >= energy_max) return(exp(log_csda[n_energy-1]));

   // energy_min < energy < energy_max
   log_energy = log(energy);
   r_bin      = (log_energy-log_energy_min)*inv_delta;
   i_bin      = int(r_bin);
   r_bin      = r_bin - float(i_bin);
   csda       = exp( (ONE-r_bin)*log_csda[i_bin] + r_bin*log_csda[i_bin+1] );
   return(csda);
}
