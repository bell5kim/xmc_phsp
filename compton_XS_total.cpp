/*****************************************************************************
 * compton_XS_total.cpp:                                                     *
 *    class member functions for:                                            *
 *       compton_XS_total: total compton cross section                       *
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
#include "compton_XS_total.h"

/*****************************************************************************
 * member functions for class compton_XS_total:                              *
 *    total compton cross section                                            *
 *****************************************************************************/

// construct the Compton cross section by reading the electron density
// from the cross section file
// attention: we do not use the Compton data from this file!
compton_XS_total::compton_XS_total(char *file_name)
{
   // open data file
   ifstream file(file_name, ios::in);
   if (!file)
   {
      xvmc_error("compton_XS_total::compton_XS_total",
                 "cannot open file",8);
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

   // close file
   file.close();

   // initialize cross section
   init(electron_density);
}

// get total compton cross section for the specified energy
real compton_XS_total::get(real energy) const
{
   real kappa; // electron rest mass divided by photon energy (EMASS/energy)
   real sigma; // total cross section

   kappa = EMASS/energy;

   sigma = 4.0*kappa + 2.0*(1.0+kappa)/(2.0+kappa)/(2.0+kappa)
         + (1.0 - 2.0*kappa*(1.0+kappa))*log((2.0+kappa)/kappa);
   return(sigma*kappa*factorc);
}
