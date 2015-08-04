#ifndef _E_BEAM_MODEL_H_
#define _E_BEAM_MODEL_H_

/*****************************************************************************
 * e_beam_model.h:                                                           *
 *    class declarations and inline member functions for:                    *
 *       e_beam_model:   abstract base class to model electron beams         *
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
 
#include "beam_model.h"

#include "beam_modifier.h"

// abstract base class to model electron beams
class e_beam_model : public beam_model
{
   public:
      e_beam_model(const beam_core *this_beam)
         : beam_model(this_beam)
      {
         app_distance = 95.0;   // default: 95 cm
         app_x1       = open_x1*app_distance/iso_distance;
         app_x2       = open_x2*app_distance/iso_distance;
         app_y1       = open_y1*app_distance/iso_distance;
         app_y2       = open_y2*app_distance/iso_distance;
         app_width_x  = app_x2 - app_x1;
         app_width_y  = app_y2 - app_y1;
         if ( (modifier =
               new beam_modifier(app_distance-ONE,app_distance)) == NULL)
         {
            xvmc_error("e_beam_model::e_beam_model",
                       "cannot create empty beam modifier",8);
         }
      }

      // pointer to the beam modifier is initialized by an empty modifier
      // later we can replace this by a contour, block etc.
      beam_modifier *modifier;

      // get source to applicator distance
      real get_app_distance(void) const { return(app_distance); }

      // get particle parameters from the beam
      virtual bool emit(particle_parameters &,
#ifdef USE_SOBOL
                        sobseq &, int &,
#endif
                        ranmar &) = 0;

      // get particle type, particle energy and set flag if the
      // particle (electron) is from the primary source
      virtual void emit(particle_type &, real &, bool &, ranmar &) = 0;

      // get particle weight, starting position, direction and voxel index
      virtual bool emit(real &, real_3 &, real_3 &, int_3 &, bool &,
#ifdef USE_SOBOL
                        sobseq &, int &,
#endif
                        ranmar &) = 0;

      // add photon background to the dose distribution
      virtual bool add_photons(real &cpu_time, int n_batch)
      {
         return(false);
      }

   protected:
      real      app_distance;   // origin to applicator distance
      real      app_x1;         // applicator X1 position
      real      app_x2;         // applicator X2 position
      real      app_y1;         // applicator Y1 position
      real      app_y2;         // applicator Y2 position
      real      app_width_x;    // app_x2 - app_x1
      real      app_width_y;    // app_y2 - app_y1
};

#endif  /* _E_BEAM_MODEL_H_ */
