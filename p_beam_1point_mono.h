#ifndef _P_BEAM_1POINT_MONO_H_
#define _P_BEAM_1POINT_MONO_H_

/*****************************************************************************
 * p_beam_1point_mono.h:                                                     *
 *    class declarations and inline member functions for:                    *
 *       p_beam_1point_mono: point source beam, mono-energetic photons       *
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
 
#include "p_beam_model.h"

// point source beam, mono-energetic photons
class p_beam_1point_mono : public p_beam_model
{
   public:
      p_beam_1point_mono(const beam_core *this_beam)
         : p_beam_model(this_beam) {}

      real average_energy(void) const { return(nominal_energy); }
      real maximum_energy(void) const { return(nominal_energy); }

      bool emit(particle_parameters &,
#ifdef USE_SOBOL
                sobseq &, int &,
#endif
                ranmar &);
   protected:
      // here "sample_energy" just returns the nominal energy,
      // however, derived beam models will overload this function
      virtual real sample_energy(ranmar &rndm) { return(nominal_energy); }
};

#endif  /* _P_BEAM_1POINT_MONO_H_ */
