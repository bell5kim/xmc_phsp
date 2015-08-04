#ifndef _E_BEAM_1POINT_POLY_H_
#define _E_BEAM_1POINT_POLY_H_

/*****************************************************************************
 * e_beam_1point_poly.h:                                                     *
 *    class declarations and inline member functions for:                    *
 *       e_beam_1point_poly: point source beam, poly-energetic electrons     *
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
 
#include "e_beam_model.h"

// point source beam, poly-energetic electrons
class e_beam_1point_poly : public e_beam_model
{
   public:
      e_beam_1point_poly(const beam_core *this_beam)
         : e_beam_model(this_beam)
      {
         find_base_data(this_beam);
         get_base_data();
         close_base_file();
      }

      real average_energy(void) const { return(energy_mean); }
      real maximum_energy(void) const { return(energy_max); }

      // get particle parameters from the beam
      bool emit(particle_parameters &,
#ifdef USE_SOBOL
                sobseq &, int &,
#endif
                ranmar &);

      // get particle type, particle energy and set flag if the
      // particle (electron) is from the primary source
      void emit(particle_type &type,    real   &energy,
                bool &primary_particle, ranmar &rndm);

      // get particle weight, starting position, direction and voxel index
      bool emit(real &, real_3 &, real_3 &, int_3 &, bool &,
#ifdef USE_SOBOL
                sobseq &, int &,
#endif
                ranmar &);

      // add photon background to the dose distribution
      bool add_photons(real &, int);

   protected:
      // electron point source model parameters
      real      sigma_theta_x;  // angular spread (radians)
      real      energy_min;     // minimum energy for spectrum (MeV)
      real      energy_max;     // maximum and most probable energy (MeV)
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

      // get base data from file, initialize energy spectrum
      void      get_base_data(void);

      // sample energy from spectrum
      real      sample_energy(ranmar &) const;
};

// sample energy from spectrum
inline real e_beam_1point_poly::sample_energy(ranmar &rndm) const
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

// get particle type, particle energy and set flag if the
// particle (electron) is from the primary source
inline void e_beam_1point_poly::emit(particle_type &type,    real   &energy,
                                     bool &primary_particle, ranmar &rndm)
{
   type   = ELECTRON;
   energy = sample_energy(rndm);

   // here all electrons are primary particles
   primary_particle = true;

   return;
}

#endif  /* _E_BEAM_1POINT_POLY_H_ */
