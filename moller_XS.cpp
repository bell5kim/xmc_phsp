/*****************************************************************************
 * moller_XS.cpp:                                                            *
 *    class member functions for:                                            *
 *       moller_XS:      Moller cross section                                *
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

#include "moller_XS.h"

/*****************************************************************************
 * member functions for class moller_XS:                                     *
 *    Moller cross section                                                   *
 *****************************************************************************/

// construct the Moller cross section by the minimum kinetic delta
// electron energy (in EGS: AE-0.511) and the electron density in
// units of 10^23 cm^-3 (water: 3.343)
moller_XS::moller_XS(real ini_e_min, real edens) :
   factor1(edens * ONE_PIxR0xR0),
   e_min(ini_e_min),
   TWOe_min(TWO*ini_e_min)
{}

// calculate total Moller cross section for the specified energy
real moller_XS::total(real &energy) const
{
   if (energy <= TWOe_min) return(ZERO);

   // total electron energy, kinetic + rest (EMASS)
   real t = energy + EMASS;

   // additional factor: 2m/(beta**2), beta = v/c
   real factor2 = TWOxEMASS*t*t/(energy*(energy+TWOxEMASS));

   // total cross section
   real sigma = factor1*factor2/energy;

#ifdef MOLLER_EXACT
   // exact Moller cross section (see SLAC-265 page 58)
   real xi1   = e_min/energy;
   real xi1p  = ONE - xi1;
   real c1    = energy*energy/(t*t);
   real c2    = EMASS*(energy+t)/(t*t);

   sigma *= (c1*(ONE_HALF-xi1) + ONE/xi1 - ONE/xi1p - c2*log(xi1p/xi1));
#else
   // Moller cross section approximated by the leading (1/E^2) term
   sigma *= (energy/e_min - TWO);
#endif /* MOLLER_EXACT */

   return(sigma);
}

// calculate first moment for the given energy
real moller_XS::first(real &energy) const
{
   if (energy <= TWOe_min) return(ZERO);

   // total electron energy, kinetic + rest (EMASS)
   real t = energy + EMASS;

   // additional factor: 2m/(beta**2), beta = v/c
   real factor2 = TWOxEMASS*t*t/(energy*(energy+TWOxEMASS));

   // first moment
   real moment = factor1*factor2;

#ifdef MOLLER_EXACT
   // exact Moller cross section (see SLAC-265 page 58)
   real xi1   = e_min/energy;
   real xi1p  = ONE - xi1;
   real c1    = energy*energy/(t*t);
   real c2    = EMASS*(energy+t)/(t*t);

   moment *= (TWO - ONE/xi1p - log(4.0*xi1*xi1p) +
              0.5*c1*(0.25-xi1*xi1) - c2*log(2.0*xi1p));
#else
   // Moller cross section approximated by the leading (1/E^2) term
   moment *= log(energy/TWO/e_min);
#endif /* MOLLER_EXACT */

   return(moment);
}
