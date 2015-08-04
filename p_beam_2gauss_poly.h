#ifndef _P_BEAM_2GAUSS_POLY_H_
#define _P_BEAM_2GAUSS_POLY_H_

/*****************************************************************************
 * p_beam_2gauss_poly.h:                                                     *
 *    class declarations and inline member functions for:                    *
 *       p_beam_2gauss_poly: two Gaussian sources, poly-energetic photons    *
 *                                                                           *
 * Copyright (C) 2002    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 31.01.2002      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include "p_beam_2gauss_mono.h"

// ###########################################################################
// HYPERION PHOTON DOSE MODELS
// ###########################################################################

// two Gaussian sources, poly-energetic photons
class p_beam_2gauss_poly : public p_beam_2gauss_mono
{
   public:

      p_beam_2gauss_poly(const beam_core *this_beam)
         : p_beam_2gauss_mono(this_beam)
      {
         init_energy();
      }

      p_beam_2gauss_poly(const p_beam_2gauss_poly &);
      real average_energy(void) const { return(energy_mean); }
      real maximum_energy(void) const { return(energy_max); }

   protected:
      // parameters for energy spectrum
      real      energy_min;    // minimum energy for spectrum (MeV)
      real      energy_max;    // maximum energy for spectrum (MeV)
      real      energy_prob;   // most probable energy (MeV)
      real      energy_mean;   // mean energy (MeV)
      real      l,b,s;         // auxiliary spectrum parameters

      // initialize energy spectrum, calculate b, l and s
      // (for an explanation see comments in file "p_beam_model.cpp")
      void      init_energy(void);

      // sample energy from spectrum
      real      sample_energy(ranmar &);
};

// sample energy from spectrum
// adapted from function "gamdev(ia,idum)" (see Numerical Recipes in C)
inline real p_beam_2gauss_poly::sample_energy(ranmar &rndm)
{
   real v1,v2,x,y,energy;

   do {
      do {
         do {
            do {
                v1 = TWO*rndm.number() - ONE;
                v2 = TWO*rndm.number() - ONE;
            } while (((v1*v1+v2*v2) > ONE) || (v1 == ZERO));
            y = v2/v1;
            x = s*y + l;
         } while (x <= ZERO);
      } while (rndm.number() > (ONE+y*y)*exp(l*log(x/l)-s*y));
      energy = x/b;
   } while ((energy < energy_min) || (energy > energy_max));

   return(energy);
}

#endif  /* _P_BEAM_2GAUSS_POLY_H_ */
