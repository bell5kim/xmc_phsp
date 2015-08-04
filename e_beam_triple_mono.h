#ifndef _E_BEAM_TRIPLE_MONO_H_
#define _E_BEAM_TRIPLE_MONO_H_

/*****************************************************************************
 * e_beam_triple_mono.h:                                                     *
 *    class declarations and inline member functions for:                    *
 *       e_beam_triple_mono: triple electron source model, mono-energetic    *
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
 
#include "e_beam_model.h"

// triple electron source model, mono-energetic
// (two Gaussian primary sources with different sigma_theta_x values
//  plus one area source to model scatter from the applicator,
//  mono-energetic electrons)
class e_beam_triple_mono : public e_beam_model
{
   public:
      e_beam_triple_mono(const beam_core *this_beam)
         : e_beam_model(this_beam)
      {
         // get beam model parameters from base data file
         find_base_data(this_beam);
         get_base_data();
         close_base_file();
      }

      virtual real average_energy(void) const { return(nominal_energy); }
      virtual real maximum_energy(void) const { return(nominal_energy); }

      // get particle parameters from the beam
      bool emit(particle_parameters &,
#ifdef USE_SOBOL
                sobseq &, int &,
#endif
                ranmar &);

      // get particle type, particle energy and set flag if the
      // particle (electron) is from the primary source
      void emit(particle_type &, real &, bool &, ranmar &);

      // get particle weight, starting position, direction and voxel index
      bool emit(real &, real_3 &, real_3 &, int_3 &, bool &,
#ifdef USE_SOBOL
                sobseq &, int &,
#endif
                ranmar &);

   protected:
      // electron triple source model parameters
      real      pri_weight1;    // weight of primary source 1 relative to
                                // the sum of both primary sources
      real      pri_xo2;        // Gaussian width squared of position fluct.
      real      sigma_theta_x1; // angular spread of primary source 1 (radians)
      real      sigma_theta_x2; // angular spread of primary source 2 (radians)
      real      sct_weight;     // weight of scatter source relative to the
                                // contribution of all sources (prim. & scat.)
      real      sct_app_dist;   // scatter source to applicator distance

      // minimum and maximum energies of the spectrum, in the mono-energetic
      // case both variables are taken to calculate the contribution of
      // the scatter source which is not!!! mono-energetic, the energy
      // spectrum of the scatter source is always uniform between
      // "energy_min" and "energy_max" also for "e_beam_triple_mono"!!!
      // "energy_max" is also the maximum and the most probable energy
      // of the spectrum in derived class "e_beam_triple_poly"!
      real      energy_min;     // minimum energy for spectrum (MeV)
      real      energy_max;     // maximum energy for spectrum (MeV)

      // further parameters from base data file needed by the derived class
      real      energy_X_aux,spectrum_A_aux,
                spectrum_alpha_aux,spectrum_beta_aux,
                photo_a_aux,photo_b_aux;

      // get base data from file
      void      get_base_data();

      // here "sample_energy" just returns the nominal energy,
      // however, derived beam models will overload this function
      virtual real sample_energy(ranmar &rndm) { return(nominal_energy); }
};

#endif  /* _E_BEAM_TRIPLE_MONO_H_ */
