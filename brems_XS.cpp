/*****************************************************************************
 * brems_XS.cpp:                                                             *
 *    class member functions for:                                            *
 *       brems_XS:       bremsstrahlung cross section                        *
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

#include "brems_XS.h"

/*****************************************************************************
 * member functions for class brems_XS:                                      *
 *    bremsstrahlung cross section                                           *
 *****************************************************************************/

// construct the bremsstrahlung cross section by the minimum
// bremsstrahlung photon energy (in EGS: AP) and the
// radiation length (water: 36.0863)
brems_XS::brems_XS(real ini_k_min, real rad_length) :
   factorb(4.0/(3.0*rad_length)),
   k_min(ini_k_min)
{}

// calculate total bremsstrahlung cross section for the specified energy
real brems_XS::total(real &energy) const
{
   if (energy <= k_min) return(ZERO);

   // total cross section
   real sigma= factorb;

#ifdef BREMS_EXACT
   sigma *= (log(energy/k_min) - 0.625 + k_min/energy
                               - 0.375*k_min*k_min/(energy*energy));
#else
   sigma *= (log(energy/k_min) - 1.000 + k_min/energy);
#endif

   return(sigma);
}

// calculate first moment for the given energy
real brems_XS::first(real &energy) const
{
   if (energy <= k_min) return(ZERO);

   // first moment
   real moment = factorb;

#ifdef BREMS_EXACT
   moment *= (0.75*energy - k_min + 0.5*k_min*k_min/energy
                          - 0.25*k_min*k_min*k_min/(energy*energy));
#else
   moment *= (0.50*energy - k_min + 0.5*k_min*k_min/energy);
#endif

   return(moment);
}
