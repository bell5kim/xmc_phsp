/*****************************************************************************
 * bhabha_XS.cpp:                                                            *
 *    class member functions for:                                            *
 *       bhabha_XS:      Bhabha cross section                                *
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

#include "bhabha_XS.h"

/*****************************************************************************
 * member functions for class bhabha_XS:                                     *
 *    Bhabha cross section                                                   *
 *****************************************************************************/

// construct the Bhabha cross section by the minimum kinetic delta
// electron energy (in EGS: AE-0.511) and the electron density in
// units of 10^23 cm^-3 (water: 3.343)
bhabha_XS::bhabha_XS(real ini_e_min, real edens) :
   factor0(edens * TWOxEMASS * ONE_PIxR0xR0),
   e_min(ini_e_min),
   TWOe_min(TWO*ini_e_min)
{}

// calculate total bhabha cross section for the specified positron energy
// (see SLAC-265 pages 59-61, PIRS-701 page 66)
real bhabha_XS::total(real &energy) const
{
   if (energy <= TWOe_min) return(ZERO);

   real tau       = energy/EMASS;
   real yy        = ONE/(tau+TWO);
   real inv_beta2 = (tau+ONE)*(tau+ONE)/(tau+TWO)/tau;
   real b0        = ONE-TWO*yy;
   real b1        = TWO-yy*yy;
   real b2        = b0*(3.0+yy*yy);
   real b3        = TWO*(ONE-yy)*b0*b0;
   real b4        = b0*b0*b0;
   real xi        = e_min/energy;
   real inv_xi    = ONE/xi;

   // total cross section
   real sigma  = inv_beta2*(inv_xi-ONE) - b1*log(inv_xi) + b2*(ONE-xi)
               + b4/3.0-0.5*b3 - xi*xi*(xi*b4/3.0-0.5*b3);
   sigma      *= factor0/energy;

   return(sigma);
}
