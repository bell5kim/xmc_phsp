#ifndef _E_BEAM_TRIPLE_POLY_H_
#define _E_BEAM_TRIPLE_POLY_H_

/*****************************************************************************
 * e_beam_triple_poly.h:                                                     *
 *    class declarations and inline member functions for:                    *
 *       e_beam_triple_poly: triple electron source model, poly-energetic    *
 *                                                                           *
 * Copyright (C) 2002    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 08.07.2002      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************
 
#include "e_beam_triple_mono.h"

// triple electron source model, poly-energetic
// (two Gaussian primary sources with different sigma_theta_x values
//  plus one area source to model scatter from the applicator,
//  poly-energetic electrons)
class e_beam_triple_poly : public e_beam_triple_mono
{
   public:
      e_beam_triple_poly(const beam_core *this_beam)
         : e_beam_triple_mono(this_beam)
      {
         init_energy();
         init_photons();
      }

      real average_energy(void) const { return(energy_mean); }
      real maximum_energy(void) const { return(energy_max); }

      // add photon background to the dose distribution
      bool add_photons(real &, int);

   protected:
      // parameters for the energy spectrum
      real      energy_X;       // energy spectrum parameter X (MeV)
      real      spectrum_A;     // energy spectrum parameter A
      real      spectrum_alpha; // energy spectrum parameter alpha
      real      spectrum_beta;  // energy spectrum parameter beta
      real      photo_a;        // parameter a for photon background
      real      photo_b;        // parameter b for photon background

      // calculated spectrum parameters
      real      energy_mean;    // mean energy (MeV)
      real      w1,w2,w3;       // weights of the different contributions
      real      b,c1,c2,d1,d2,e3,r1,r2; // sampling parameters

      // initialize energy spectrum
      void      init_energy(void);

      // initialize photon background
      void      init_photons(void);

      // sample energy from spectrum
      real      sample_energy(ranmar &);
};

// sample energy from spectrum
// (equivalent to "e_beam_1point_poly::sample_energy")
inline real e_beam_triple_poly::sample_energy(ranmar &rndm)
{
   // energy
   real energy = ZERO;

   // we need two random numbers
   real rndm1 = rndm.number();
   real rndm2 = rndm.number();

   if (rndm1 < w1)
   {
      // main contribution
      energy = energy_X - b/sqrt(rndm2*c2+(ONE-rndm2)*c1);
   }
   else
   {
      if (rndm1 < w2)
      {
         // linear contribution
         energy = sqrt(rndm2*d1+(ONE-rndm2)*d2);
      }
      else
      {
         if (rndm1 < w3)
         {
            // constant contribution
            energy = energy_min + rndm2*e3;
         }
         else
         {
            // low energy contribution
            energy = -log(rndm2*r2+(ONE-rndm2)*r1)/spectrum_beta;
            while (energy > energy_max)
            {
               rndm2  =  rndm.number();
               energy = -log(rndm2*r2+(ONE-rndm2)*r1)/spectrum_beta;
            }
            if (energy < energy_min) energy = energy_min;
         }
      }
   }

   return(energy);
}

#endif  /* _E_BEAM_TRIPLE_POLY_H_ */
