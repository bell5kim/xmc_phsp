#ifndef _E_BEAM_1POINT_MONO_H_
#define _E_BEAM_1POINT_MONO_H_

/*****************************************************************************
 * e_beam_1point_mono.h:                                                     *
 *    class declarations and inline member functions for:                    *
 *       e_beam_1point_mono: point source beam, mono-energetic electrons     *
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

// point source beam, mono-energetic electrons
class e_beam_1point_mono : public e_beam_model
{
   public:
      e_beam_1point_mono(const beam_core *this_beam)
         : e_beam_model(this_beam) {}

      real average_energy(void) const { return(nominal_energy); }
      real maximum_energy(void) const { return(nominal_energy); }

      // get particle parameters from the beam
      bool emit(particle_parameters &,
#ifdef USE_SOBOL
                sobseq &, int &,
#endif
                ranmar &);

      // get particle type, particle energy and set flag if the
      // particle (electron) is from the primary source
      void emit(particle_type &type,    real   &energy,
                bool &primary_particle, ranmar &rndm)
      {
         type             = ELECTRON;
         energy           = nominal_energy;

         // here all electrons are primary particles
         primary_particle = true;

         return;
      }

      // get particle weight, starting position, direction and voxel index
      bool emit(real &, real_3 &, real_3 &, int_3 &, bool &,
#ifdef USE_SOBOL
                sobseq &, int &,
#endif
                ranmar &);
};

#endif  /* _E_BEAM_1POINT_MONO_H_ */
