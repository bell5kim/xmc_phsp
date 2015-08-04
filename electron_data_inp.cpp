/*****************************************************************************
 * electron_data_inp.cpp:                                                    *
 *    class member functions for:                                            *
 *       electron_data_inp: input electron transport data                    *
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

#include "definitions.h"
#include "global.h"
#include "electron_data_inp.h"

/*****************************************************************************
 * member functions for electron_data_inp:                                   *
 *****************************************************************************/

//########################NEWNEWNEWNEWNEWNEW##################################

electron_data_inp::electron_data_inp(int nTableSize) :
  num(static_cast<unsigned int>(nTableSize))
{
  // allocate memory for input electron data
  bool error = false;
  if ( (energy  = new real[num]) == NULL ) error = true;
  if ( (s_col   = new real[num]) == NULL ) error = true;
  if ( (s_rad   = new real[num]) == NULL ) error = true;
  if ( (s_tot   = new real[num]) == NULL ) error = true;
  if ( (s_air   = new real[num]) == NULL ) error = true;
  if ( (s_photo = new real[num]) == NULL ) error = true;
  if ( (s_scat  = new real[num]) == NULL ) error = true;

  if (error)
    {
      xvmc_error("electron_data_inp::electron_data_inp",
                 "cannot allocate memory for electron input data",8);
    }
}

//########################NEWNEWNEWNEWNEWNEW##################################

electron_data_inp* ElectronDataFactory(const char* filename)
{
  int num;
  char comma;     // the numbers are separated by commas
  real dummy;     // dummy variable
  real ee;        // kinetic energy
  real p2;        // momentum squared
  real xx;        // inverse screening angle squared
  real beta2;     // beta squared, (v/c)^2
  real ts;        // scattering power

  ifstream file(filename, ios::in);   // open data file
  if (!file)
    {
      xvmc_error("electron_data_inp::electron_data_inp","cannot open file",8);
    }

  // read the number of file entries
  file >> num;

  electron_data_inp* pedi = new electron_data_inp(num);

  int i=0;
  for (i=0; i<num; ++i)
    {
      file >> pedi->energy[i]   >> comma
           >> dummy             >> comma
           >> dummy             >> comma
           >> dummy             >> comma
           >> dummy             >> comma
           >> dummy             >> comma
           >> dummy             >> comma
           >> pedi->s_col[i]    >> comma
           >> pedi->s_rad[i]    >> comma
           >> pedi->s_tot[i]    >> comma
           >> pedi->s_scat[i]   >> comma
           >> pedi->s_air[i]    >> comma
           >> pedi->s_photo[i];

      // change units: 1 MeV m^2/kg = 10 MeV cm^2/g
      pedi->s_col[i] *= 10.0;
      pedi->s_rad[i] *= 10.0;
      pedi->s_tot[i]  = pedi->s_col[i] + pedi->s_rad[i]; // don't use s_tot from file
      pedi->s_scat[i] *= 10.0;
      // ratio s_air/s_h2o to calculate ionization instead of dose
      // s_air from ICRU report 37 page 120 (air, dry, near sea level)
      pedi->s_air[i]  = 10.0*pedi->s_air[i]/pedi->s_col[i];
      // ratio s_photo/s_h2o to calculate film dose
      // s_photo from ICRU report 37 page 172 (photographic emulsion)
      pedi->s_photo[i]  = 10.0*pedi->s_photo[i]/pedi->s_col[i];

      // re-calculate electron scattering power for linear interpolation
      ee    = pedi->energy[i];            // kinetic energy
      p2    = ee*(ee+TWOxEMASS);    // momentum squared
      xx    = 4.0*p2/AMU2;          // inverse screening angle squared
      beta2 = p2/(p2+EMASSxEMASS);  // beta squared, (v/c)^2
      // estimated scattering power
      ts    = 4.0*BC*((ONE+ONE/xx)*log(ONE+xx)-ONE)/xx/beta2;
      pedi->s_scat[i]=pedi->s_scat[i]/ts;
   }

   file.close();

   return pedi;
}

// input electron transport data from file
electron_data_inp::electron_data_inp(char *file_name)
{
   char comma;     // the numbers are separated by commas
   real dummy;     // dummy variable
   real ee;        // kinetic energy
   real p2;        // momentum squared
   real xx;        // inverse screening angle squared
   real beta2;     // beta squared, (v/c)^2
   real ts;        // scattering power

   ifstream file(file_name, ios::in);   // open data file
   if (!file)
   {
      xvmc_error("electron_data_inp::electron_data_inp","cannot open file",8);
   }

   // read the number of file entries
   file >> num;

   // allocate memory for input electron data
   bool error = false;
   if ( (energy  = new real[num]) == NULL ) error = true;
   if ( (s_col   = new real[num]) == NULL ) error = true;
   if ( (s_rad   = new real[num]) == NULL ) error = true;
   if ( (s_tot   = new real[num]) == NULL ) error = true;
   if ( (s_air   = new real[num]) == NULL ) error = true;
   if ( (s_photo = new real[num]) == NULL ) error = true;
   if ( (s_scat  = new real[num]) == NULL ) error = true;

   if (error)
   {
      xvmc_error("electron_data_inp::electron_data_inp",
                 "cannot allocate memory for electron input data",8);
   }

   for (register unsigned int i=0; i<num; ++i)
   {
      file >> energy[i]   >> comma
           >> dummy       >> comma
           >> dummy       >> comma
           >> dummy       >> comma
           >> dummy       >> comma
           >> dummy       >> comma
           >> dummy       >> comma
           >> s_col[i]    >> comma
           >> s_rad[i]    >> comma
           >> s_tot[i]    >> comma
           >> s_scat[i]   >> comma
           >> s_air[i]    >> comma
           >> s_photo[i];

      // change units: 1 MeV m^2/kg = 10 MeV cm^2/g
      s_col[i]  = 10.0*s_col[i];
      s_rad[i]  = 10.0*s_rad[i];
      s_tot[i]  = s_col[i] + s_rad[i]; // don't use s_tot from file
      s_scat[i] = 10.0*s_scat[i];
      // ratio s_air/s_h2o to calculate ionization instead of dose
      // s_air from ICRU report 37 page 120 (air, dry, near sea level)
      s_air[i]  = 10.0*s_air[i]/s_col[i];
      // ratio s_photo/s_h2o to calculate film dose
      // s_photo from ICRU report 37 page 172 (photographic emulsion)
      s_photo[i]  = 10.0*s_photo[i]/s_col[i];

      // re-calculate electron scattering power for linear interpolation
      ee    = energy[i];            // kinetic energy
      p2    = ee*(ee+TWOxEMASS);    // momentum squared
      xx    = 4.0*p2/AMU2;          // inverse screening angle squared
      beta2 = p2/(p2+EMASSxEMASS);  // beta squared, (v/c)^2
      // estimated scattering power
      ts    = 4.0*BC*((ONE+ONE/xx)*log(ONE+xx)-ONE)/xx/beta2;
      s_scat[i]=s_scat[i]/ts;

   }

   file.close();

}

// delete data
electron_data_inp::~electron_data_inp()
{
   num = 0;
   delete [] energy;
   delete [] s_col;
   delete [] s_rad;
   delete [] s_tot;
   delete [] s_air;
   delete [] s_photo;
   delete [] s_scat;
}
