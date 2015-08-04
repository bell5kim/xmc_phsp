/*****************************************************************************
 * pair_XS_total.cpp:                                                        *
 *    class member functions for:                                            *
 *       pair_XS_total:  total pair cross section                            *
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
#include <string>
using namespace std;

#include "global.h"
#include "pair_XS_total.h"

/*****************************************************************************
 * member functions for class pair_XS_total:                                 *
 *    total pair cross section data                                          *
 *****************************************************************************/


PairTotData::PairTotData(int nEnergy, real energyMin, real energyMax) :
  n_energy(nEnergy),       // number of bins for photon energy
  energy_min(energyMin),   // minimum energy
  energy_max(energyMax),   // maximum energy
  sigma(NULL)          // pointer to total cross section data
{
  sigma = new real[n_energy];
}

//##################################################################################

PairTotData::~PairTotData()
{
  delete[] sigma;
}

//##################################################################################

PairTotData* PairTotFactory(const std::string& file_name)

{
  int n_energy;
  real energy_min, energy_max;

  // open data file
  ifstream file(file_name.c_str(), ios::in);
  if (!file)
    {
      xvmc_error("pair_XS_total::read","cannot open file",8);
    }

  // read the number of bins
  file >> n_energy;
  // read minimum and maximum energy and calculate logarithms
  file >> energy_min >> energy_max;

  PairTotData* ppxst = new PairTotData(n_energy, energy_min, energy_max);

  // read data from file
  int i;
  for (i=0; i<n_energy; ++i) {
    file >> ppxst->sigma[i];
  }

   // close file
   file.close();

   if(ppxst->invariant() == false) {
    delete ppxst;
    return NULL;
  }

  return ppxst;
}

//##################################################################################

pair_XS_total::pair_XS_total(const PairTotData& A) :
  n_energy(A.n_energy),       // number of bins for photon energy
  energy_min(A.energy_min),   // minimum energy
  energy_max(A.energy_max),   // maximum energy
  sigma(NULL)                 // pointer to total cross section data
{
  log_energy_min = log(energy_min);
  log_energy_max = log(energy_max);

  // calculate inverse bin size for photon energy (log scale)
  inv_delta = float(n_energy-1)/(log_energy_max - log_energy_min);

  sigma = new real[n_energy];

  int i;
  for (i=0; i<n_energy; ++i) {
    sigma[i] = A.sigma[i];
  }
}

//##################################################################################

pair_XS_total::~pair_XS_total(void)
{
   delete [] sigma;
}

//##################################################################################

// input total pair data from NIST (XCOM) file
void pair_XS_total::nist(char *file_name)
{
   // delete old cross section data
   if (sigma != NULL)
   {
      delete [] sigma; sigma = NULL;
   }

   // open data file
   ifstream file(file_name, ios::in);
   if (!file)
   {
      xvmc_error("pair_XS_total::nist","cannot open file",8);
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
      xvmc_error("pair_XS_total::nist",
                 "cannot allocate memory for photon energy input",8);
   }

   // allocate memory for total pair cross section input data
   real *mu_pair_inp = NULL;
   if ( (mu_pair_inp  = new real[num_inp]) == NULL )
   {
      xvmc_error("pair_XS_total::nist",
                 "cannot allocate memory for pair cross section input",8);
   }

   // read data from file
   real energy_inp  = ZERO;
   real dummy       = ZERO;
   real mu_pair_nuc = ZERO;
   real mu_pair_ele = ZERO;
   for (register unsigned int i=0; i<num_inp; ++i)
   {
      file >> energy_inp
           >> dummy            // coherent scattering
           >> dummy            // Compton  scattering
           >> dummy            // photoelectric absorption
           >> mu_pair_nuc      // pair production in nuclear field
           >> mu_pair_ele      // pair production in electron field
           >> dummy            // total attenuation with coherent scatter
           >> dummy;           // total attenuation without coherent scatter

      // log scale
      log_energy_inp[i] = log(energy_inp);

      // mass -> linear attenuation
      mu_pair_inp[i] = (mu_pair_nuc + mu_pair_ele)*density;
   }

   // close file
   file.close();

   // set number of photon energy bins
   n_energy = 2000;

   // minimum energy
   log_energy_min = log_energy_inp[0];
   energy_min     = exp(log_energy_min);
   if (energy_min <= TWOxEMASS)
   {
      energy_min     = TWOxEMASS;
      log_energy_min = log(energy_min);
   }

   // maximum energy
   log_energy_max = log_energy_inp[num_inp-1];
   energy_max     = exp(log_energy_max);
   if (energy_max <= energy_min)
   {
      xvmc_error("pair_XS_total::nist",
                 "maximum energy in input file too small",8);
   }

   // calculate inverse bin size for photon energy (log scale)
   inv_delta = float(n_energy-1)/(log_energy_max - log_energy_min);

   // bin size
   real delta = ONE/inv_delta;

   // allocate memory for pair data
   if ( (sigma = new real[n_energy]) == NULL )
   {
      xvmc_error("pair_XS_total::nist",
                 "cannot allocate memory for pair cross section data",8);
   }

   // interpolate total pair cross section data
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
      real mu_pair = q*mu_pair_inp[lower] + p*mu_pair_inp[upper];

      // divide mu_pair by factor (1-2m/k)^3 and store
      real sig2 = ONE - TWOxEMASS/exp(log_energy);
      sigma[j] = mu_pair/(sig2*sig2*sig2);

      // next energy
      log_energy += delta;
   }

   // delete input arrays
   delete [] log_energy_inp;  log_energy_inp  = NULL;
   delete [] mu_pair_inp;     mu_pair_inp     = NULL;
}

// input total pair data from file
void pair_XS_total::read(char *file_name)
{
   // delete old cross section data
   if (sigma != NULL)
   {
      delete [] sigma; sigma = NULL;
   }

   // open data file
   ifstream file(file_name, ios::in);
   if (!file)
   {
      xvmc_error("pair_XS_total::read","cannot open file",8);
   }

   // read the number of bins
   file >> n_energy;

   // read minimum and maximum energy and calculate logarithms
   file >> energy_min >> energy_max;
   log_energy_min = log(energy_min);
   log_energy_max = log(energy_max);

   // calculate inverse bin size for photon energy (log scale)
   inv_delta = float(n_energy-1)/(log_energy_max - log_energy_min);

   // allocate memory for pair data
   if ( (sigma = new real[n_energy]) == NULL )
   {
      xvmc_error("pair_XS_total::read",
                 "cannot allocate memory for pair cross section data",8);
   }

   // read data from file
   for (register int i=0; i<n_energy; ++i) {
      file >> sigma[i]; }

   // close file
   file.close();
}

// get total pair cross section for the specified energy
real pair_XS_total::get(real energy) const
{
   real  log_energy;     // ln(energy)
   int   i_bin;          // bin number, integer value
   real  r_bin;          // bin number, real value
   real  sig1,sig2;

   if (energy <= energy_min)
   {
      return(ZERO);
   }

   if (energy >= energy_max)
   {
      sig1 = sigma[n_energy-1];
   }
   else
   {
      // energy_min < energy < energy_max
      log_energy = log(energy);
      r_bin      = (log_energy-log_energy_min)*inv_delta;
      i_bin      = int(r_bin);
      r_bin      = r_bin - float(i_bin);
      sig1       = (ONE-r_bin)*sigma[i_bin] + r_bin*sigma[i_bin+1];
   }
   sig2 = ONE - TWOxEMASS/energy;
   return(sig1*sig2*sig2*sig2);
}
