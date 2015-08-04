#ifndef _MOLLER_XS_H_
#define _MOLLER_XS_H_

/*****************************************************************************
 * moller_XS.h:                                                              *
 *    class declarations and inline member functions for:                    *
 *     class moller_XS:  Moller cross section                                *
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

#include "definitions.h"
#include "ranmar.h"

#define MOLLER_EXACT

/*****************************************************************************
 * class moller_XS:                                                          *
 *    Moller cross section                                                   *
 *****************************************************************************/

class moller_XS
{
   public:
      // construct the Moller cross section by the minimum kinetic delta
      // electron energy (in EGS: AE-0.511) and the electron density in
      // units of 10^23 cm^-3 (water: 3.343)
      moller_XS(real ini_e_min, real edens=3.343);

      // calculate total Moller cross section for the specified energy
      real total(real &energy) const;

      // calculate first moment for the given energy
      real first(real &energy) const;

      // sample Moller interaction parameters
      inline bool interaction(real energy, ranmar &rndm,
                     real &energy_n, real &cos_t_n, real &sin_t_n,
                     real &energy_d, real &cos_t_d, real &sin_t_d) const;

   private:
      moller_XS();        // default constructor not implemented
      real     factor1;   // material dependent factor (n_e x pi x r_0 x r_0)
      real     e_min;     // minimum kinetic delta electron energy
      real     TWOe_min;  // 2 x e_min
};

// sample Moller interaction parameters
inline bool moller_XS::interaction(real energy, ranmar  &rndm,
                          real &energy_n, real &cos_t_n, real &sin_t_n,
                          real &energy_d, real &cos_t_d, real &sin_t_d) const
{
   // if the energy is too small we don't have a Moller interaction
   if (energy <= TWOe_min)
   {
      energy_n = energy;
      cos_t_n  = ONE;
      sin_t_n  = ZERO;
      energy_d = ZERO;
      cos_t_d  = ZERO;
      sin_t_d  = ONE;
      return(false);
   }

#ifdef MOLLER_EXACT
   // sample secondary electron energy from the exact Moller cross section
   // using a method described in the EGSnrc manual (PIRS-701)
   real t   = energy + EMASS; // total electron energy, kinetic + rest (EMASS)
   real c2e = EMASS*(energy+t)/(energy*t*t); // c2/energy
   real c3  = ONE/(TWO*t*t);                 // 1/(2*t*t)
   // maximum to normalize the rejection function
   real gmax = 3.0*energy*energy/(2.0*t*t);

   real rej_weight; // rejection weight
   do
   {
      // pick a random number
      real eta = rndm.number();

      // sample delta energy from (1/energy_d)^2 distribution
      energy_d = (energy-e_min)*e_min/((energy-e_min)*(ONE-eta)+e_min*eta);

      // calculate rejection weight (function)
      rej_weight = (ONE+(c3*energy_d-c2e)*energy_d)/gmax;
   }
   while (rndm.number() > rej_weight);

   // the secondary electron is the electron with smaller energy
   energy_d = min_of(energy_d,energy-energy_d);
#else
   // sample energy from approximated Moller cross section (1/E^2)
   real eta = rndm.number();
   energy_d = energy*e_min/(energy*(ONE-eta)+TWOe_min*eta);
#endif /* MOLLER_EXACT */

   // new energy of the primary electron
   energy_n = energy - energy_d;

   // calculate polar scattering angle of the delta electron
   cos_t_d = energy_d*(energy+TWOxEMASS)/energy/(energy_d+TWOxEMASS);
   if (cos_t_d < ONE)
   {
      sin_t_d = sqrt(ONE-cos_t_d);
      cos_t_d = sqrt(cos_t_d);
   }
   else
   {
      sin_t_d = ZERO;
      cos_t_d = ONE;
   }

   // calculate polar scattering angle of the primary electron
   cos_t_n = energy_n*(energy+TWOxEMASS)/energy/(energy_n+TWOxEMASS);
   if (cos_t_n < ONE)
   {
      sin_t_n = sqrt(ONE-cos_t_n);
      cos_t_n = sqrt(cos_t_n);
   }
   else
   {
      sin_t_n = ZERO;
      cos_t_n = ONE;
   }

   return(true);
}

#endif /* _MOLLER_XS_H_ */
