/*****************************************************************************
 * photon_datai_inp.cpp:                                                     *
 *    class member functions for:                                            *
 *       photon_data_inp: input photon transport data                        *
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

#include <fstream>
using namespace std;

//#include "definitions.h"
#include "global.h"
#include "photon_data_inp.h"

/*****************************************************************************
 * member functions for photon_data_inp:                                     *
 *****************************************************************************/

// input photon transport data from file
photon_data_inp::photon_data_inp(char *file_name)
{
   char comma;     // the numbers are separated by commas
   real dummy;     // dummy variable

   ifstream file(file_name, ios::in);  // open data file
   if (!file)
   {
      xvmc_error("photon_data_inp::photon_data_inp","cannot open file",8);
   }

   // read the number of file entries
   file >> num;

   // allocate memory for input photon data
   bool error = false;
   if ( (energy  = new real[num]) == NULL ) error = true;
   if ( (mu_comp = new real[num]) == NULL ) error = true;
   if ( (mu_pair = new real[num]) == NULL ) error = true;
   if ( (mu_phot = new real[num]) == NULL ) error = true;
   if ( (mu_tot  = new real[num]) == NULL ) error = true;
   if ( (mu_en   = new real[num]) == NULL ) error = true;

   if (error)
   {
      xvmc_error("photon_data_inp::photon_data_inp",
                 "cannot allocate memory for photon input data",8);
   }

   for (register unsigned int i=0; i<num; ++i)
   {
      file >> energy[i]   >> comma
           >> mu_phot[i]  >> comma
           >> dummy       >> comma  // reseved for Rayleigh cross section
           >> mu_comp[i]  >> comma
           >> mu_pair[i]  >> comma
           >> mu_tot[i]   >> comma
           >> mu_en[i]    >> comma
           >> dummy       >> comma
           >> dummy       >> comma
           >> dummy       >> comma
           >> dummy       >> comma
           >> dummy       >> comma
           >> dummy;

      // change units: 1 m^2/kg = 10 cm^2/g
      mu_phot[i] = 10.0*mu_phot[i];
      mu_comp[i] = 10.0*mu_comp[i];
      mu_pair[i] = 10.0*mu_pair[i];
      // don't use mu_tot from file
      mu_tot[i]  = mu_phot[i] + mu_comp[i] + mu_pair[i];
      mu_en[i]   = 10.0*mu_en[i]*energy[i];  // mu_en --> KERMA
   }

   file.close();
}

//################################################################################

photon_data_inp* PhotonDataFactory(const char* filename)
{
   char comma;     // the numbers are separated by commas
   real dummy;     // dummy variable
   int num;

   ifstream file(filename, ios::in);  // open data file
   if (!file)
   {
      xvmc_error("photon_data_inp::photon_data_inp","cannot open file",8);
   }

   // read the number of file entries
   file >> num;

   photon_data_inp* ppdi = new photon_data_inp(num);

   int i;
   for (i=0; i<num; ++i) {
     file >> ppdi->energy[i]   >> comma
	  >> ppdi->mu_phot[i]  >> comma
	  >> dummy       >> comma  // reseved for Rayleigh cross section
	  >> ppdi->mu_comp[i]  >> comma
	  >> ppdi->mu_pair[i]  >> comma
	  >> ppdi->mu_tot[i]   >> comma
	  >> ppdi->mu_en[i]    >> comma
	  >> dummy       >> comma
	  >> dummy       >> comma
	  >> dummy       >> comma
	  >> dummy       >> comma
	  >> dummy       >> comma
	  >> dummy;

      // change units: 1 m^2/kg = 10 cm^2/g
      ppdi->mu_phot[i] *= 10.0;
      ppdi->mu_comp[i] *= 10.0;
      ppdi->mu_pair[i] *= 10.0;
      // don't use mu_tot from file
      ppdi->mu_tot[i]  = ppdi->mu_phot[i] + ppdi->mu_comp[i] + ppdi->mu_pair[i];
      ppdi->mu_en[i]   = 10.0*ppdi->mu_en[i]*ppdi->energy[i];  // mu_en --> KERMA
   }

   file.close();

   return ppdi;
}

//################################################################################

photon_data_inp::photon_data_inp(int nTableSize) :
  num(nTableSize)
{
  // allocate memory for input photon data
  bool error = false;
  if ( (energy  = new real[num]) == NULL ) error = true;
  if ( (mu_comp = new real[num]) == NULL ) error = true;
  if ( (mu_pair = new real[num]) == NULL ) error = true;
  if ( (mu_phot = new real[num]) == NULL ) error = true;
  if ( (mu_tot  = new real[num]) == NULL ) error = true;
  if ( (mu_en   = new real[num]) == NULL ) error = true;

  if (error) {
    xvmc_error("photon_data_inp::photon_data_inp",
	       "cannot allocate memory for photon input data",8);
  }
}

//################################################################################

photon_data_inp::~photon_data_inp()
{
   num = 0;
   delete [] energy;
   delete [] mu_comp;
   delete [] mu_pair;
   delete [] mu_phot;
   delete [] mu_tot;
   delete [] mu_en;
}
