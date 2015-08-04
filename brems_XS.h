#ifndef _BREMS_XS_H_
#define _BREMS_XS_H_

/*****************************************************************************
 * brems_XS.h:                                                               *
 *    class declarations and inline member functions for:                    *
 *     class brems_XS:   bremsstrahlung cross section                        *
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

/*****************************************************************************
 * class brems_XS:                                                           *
 *    bremsstrahlung cross section                                           *
 *****************************************************************************/

class brems_XS
{
   public:
      // construct the bremsstrahlung cross section by the minimum
      // bremsstrahlung photon energy (in EGS: AP) and the
      // radiation length (water: 36.0863)
      brems_XS(real ini_k_min, real rad_length=36.0863);

      // calculate total bremsstrahlung cross section for the specified energy
      real total(real &energy) const;

      // calculate first moment for the given energy
      real first(real &energy) const;

      // sample bremsstrahlung photon energy
      inline real sample(real energy, ranmar &rndm) const;

      // sample bremsstrahlung interaction parameters
      inline bool interaction(real energy, ranmar &rndm,
                     real &energy_n, real &cos_t_n, real &sin_t_n,
                     real &energy_x, real &cos_t_x, real &sin_t_x) const;

   private:
      brems_XS();         // default constructor not implemented
      real     factorb;   // material dependent factor
      real     k_min;     // minimum bremsstrahlung photon energy
};

// sample bremsstrahlung photon energy
inline real brems_XS::sample(real energy, ranmar &rndm) const
{
   if (energy <= k_min) return(ZERO);

#ifdef BREMS_EXACT
   real energy2    = energy*energy;     // energy squared
#endif
   real e_brems    = ZERO;    // photon energy
   real rej_weight = ONE;     // rejection weight

   do
   {
      e_brems    = k_min*exp(rndm.number()*log(energy/k_min));
#ifdef BREMS_EXACT
      rej_weight = ONE - e_brems/energy + 0.75*e_brems*e_brems/energy2;
#else
      rej_weight = ONE - e_brems/energy;
#endif
   }
   while (rndm.number() > rej_weight);

   return(e_brems);
}

// sample bremsstrahlung interaction parameters
inline bool brems_XS::interaction(real energy, ranmar &rndm,
                         real &energy_n, real &cos_t_n, real &sin_t_n,
                         real &energy_x, real &cos_t_x, real &sin_t_x) const
{
   if (energy <= k_min)
   {
      energy_n = energy;
      cos_t_n  = ONE;
      sin_t_n  = ZERO;
      energy_x = ZERO;
      cos_t_x  = ONE;
      sin_t_x  = ZERO;
      return(false);
   }

#ifdef BREMS_EXACT
   real energy2    = energy*energy;     // energy squared
#endif
   real rej_weight = ONE;               // rejection weight
   do
   {
      energy_x   = k_min*exp(rndm.number()*log(energy/k_min));
#ifdef BREMS_EXACT
      rej_weight = ONE - energy_x/energy + 0.75*energy_x*energy_x/energy2;
#else
      rej_weight = ONE - energy_x/energy;
#endif
   }
   while (rndm.number() > rej_weight);

   // new energy of the primary electron
   energy_n = energy - energy_x;

   // the direction of the electron is unchanged (as in EGS4)
   cos_t_n  = ONE;
   sin_t_n  = ZERO;

   // fixed photon angle (as in EGS4)
   sin_t_x = sin(EMASS/energy);
   cos_t_x = sqrt(max_of(ZERO,(ONE-sin_t_x)*(ONE+sin_t_x)));

   return(true);
}

#endif /* _BREMS_XS_H_ */
