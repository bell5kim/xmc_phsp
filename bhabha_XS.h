#ifndef _BHABHA_XS_H_
#define _BHABHA_XS_H_

/*****************************************************************************
 * bhabha_XS.h:                                                              *
 *    class declarations and inline member functions for:                    *
 *     class bhabha_XS:  Bhabha cross section                                *
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
 * class bhabha_XS:                                                          *
 *    Bhabha cross section                                                   *
 *****************************************************************************/

class bhabha_XS
{
   public:
      // construct the Bhabha cross section by the minimum kinetic delta
      // electron energy (in EGS: AE-0.511) and the electron density in
      // units of 10^23 cm^-3 (water: 3.343)
      bhabha_XS(real ini_e_min, real edens=3.343);

      // calculate total Bhabha cross section for the specified energy
      real total(real &energy) const;

      // sample Bhabha interaction parameters
      inline bool interaction(real energy, ranmar &rndm,
                     real &energy_n, real &cos_t_n, real &sin_t_n,
                     real &energy_d, real &cos_t_d, real &sin_t_d) const;

   private:
      bhabha_XS();  // default constructor not implemented
      real factor0; // material dependent factor (n_e x 2 x pi x r_0 x r_0 x m)
      real e_min;    // minimum kinetic delta electron energy
      real TWOe_min; // 2 x e_min
};

// sample Bhabha interaction parameters
inline bool bhabha_XS::interaction(real energy, ranmar  &rndm,
                          real &energy_n, real &cos_t_n, real &sin_t_n,
                          real &energy_d, real &cos_t_d, real &sin_t_d) const
{
   // if the energy is too small we don't have a Bhabha interaction
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

   // sample energy from approximated Bhabha cross section (1/E^2)
   real eta = rndm.number();
   energy_d = energy*e_min/(energy*(ONE-eta)+TWOe_min*eta);

   // new energy of the positron
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

   // calculate polar scattering angle of the positron
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

#endif /* _BHABHA_XS_H_ */
