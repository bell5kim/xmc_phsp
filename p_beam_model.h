#ifndef _P_BEAM_MODEL_H_
#define _P_BEAM_MODEL_H_

/*****************************************************************************
 * p_beam_model.h:                                                           *
 *    class declarations and inline member functions for:                    *
 *       p_beam_model:   abstract base class to model photon beams           *
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

// ######################################################################
// ONLY PHOTON BEAMS OF TYPE p_beam_2gauss_mono/poly ARE USED IN HYPERION
// everything else are abstract classes
// ######################################################################

// abstract base class to model photon beams
class p_beam_model : public beam_model
{
   public:
      p_beam_model() {}

      p_beam_model(const beam_core *this_beam)
         : beam_model(this_beam)
      {
         col_mdistance = 50.0;   // default: 50 cm
         col_cdistance = 50.0;   // default: 50 cm
         col_xdistance = 50.0;   // default: 50 cm
         col_ydistance = 50.0;   // default: 50 cm
         col_x1       = open_x1*col_xdistance/iso_distance;
         col_x2       = open_x2*col_xdistance/iso_distance;
         col_y1       = open_y1*col_ydistance/iso_distance;
         col_y2       = open_y2*col_ydistance/iso_distance;
         col_width_x  = col_x2 - col_x1;
         col_width_y  = col_y2 - col_y1;
         if ( (modifier = new beam_modifier()) == NULL)
         {
            xvmc_error("p_beam_model::p_beam_model",
            "cannot create empty beam modifier",8);
         }
      }

      virtual ~p_beam_model(void) {;}

      // pointer to the beam modifier is initialized by an empty modifier
      // later we can replace this by an MLC, compensator, etc.
      beam_modifier *modifier;
#ifdef PHSP_READ
      real get_col_C_distance() {return col_cdistance;}
#endif

   protected:
      // abstract base class constructor for Hyperion - does very little...
      // real constructor defined in p_beam_2gauss_mono
      real  col_mdistance;  // distance from origin to plane above the MLC or
                            // other modifier (starting plane of MC transport)
      real  col_cdistance;  // distance from origin to cookie cutter plane
      real  col_xdistance;  // distance from origin to X-collimator lower plane
      real  col_ydistance;  // distance from origin to Y-collimator lower plane
      real  col_x1;         // collimator X1 position
      real  col_x2;         // collimator X2 position
      real  col_y1;         // collimator Y1 position
      real  col_y2;         // collimator Y2 position
      real  col_width_x;    // col_x2 - col_x1
      real  col_width_y;    // col_y2 - col_y1
};

#endif  /* _P_BEAM_MODEL_H_ */
