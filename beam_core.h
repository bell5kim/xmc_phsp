#ifndef _BEAM_CORE_H_
#define _BEAM_CORE_H_

/*****************************************************************************
 * beam_core.h:                                                              *
 *    class declarations and inline member functions for:                    *
 *       beam_core:      parameters for one electron or photon beam          *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 99/12/14        *
 *    added pointer to portal dose image                  MF 03/03/13        *
 *    added pointer to jaws                               JK Feb 21, 2008    *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include "definitions.h"
#include "global.h"
#include "contour.h"
#include "multi_leaf.h"
#include "MC_jaws.h"
#include "portal_dose.h"

class beam_core
{
   public:
      int            id;                 // present beam id
      particle_type  type;               // PHOTON or ELECTRON beam
      // not implemented yet             // check latch bits (phase space file)
      // not implemented yet             // number of latch bits
      int            model_id;           // beam model id
      char          *base_key;           // base data file name (without ext.)
      char          *applicator;         // applicator id
      real           weight;             // beam weight
      long           n_history;          // number of histories to simulate
      int            n_repeat;           // number of history repetitions
      int            n_rotate;           // number of additional repetitions in
                                         // the case of gantry rotation
      int            n_batch;            // number of simulation batches
      real           energy;             // nominal beam energy (MeV)
      real_3         iso_center;         // iso-center position (cm)
      bool           gantry_rotation;    // true for gantry rotation
      real           start_gantry_angle; // start gantry angle (deg)
      real           stop_gantry_angle;  // stop gantry angle for rotation (deg)
      real           table_angle;        // table angle (deg)
      real           collimator_angle;   // collimator angle (deg)
      real           open_x1;            // beam limit X1 (cm, iso-center plane)
      real           open_x2;            // beam limit X2 (cm, iso-center plane)
      real           open_y1;            // beam limit Y1 (cm, iso-center plane)
      real           open_y2;            // beam limit Y2 (cm, iso-center plane)
      // irregular beam contour (defined at the iso-center plane, unit: cm)
      contour       *irregular;
      // multi-leaf collimator (leaf positions defined at the iso-center plane)
      multi_leaf    *mlc;
      // physical jaws (jaw positions defined at the iso-center plane) Added by JK Feb 21, 2008
      MC_jaw        *jaw;
      // pointer to the portal dose image
      portal_dose   *portal;
      // not implemented yet             // fluence profile
      real           cpu_time;           // total CPU time for this beam

  
      // return a pointer to the next beam
      beam_core     *next(void) { return(next_beam_core); }
  
   private:
      beam_core     *next_beam_core;     // MA pointer to the next beam
      beam_core     *previous_beam_core; // MA pointer to the previous beam
  
      beam_core(void);                   // create empty beam
      ~beam_core(void);                  // delete beam, free memory
  
   friend class treatment_plan;
};

#endif  /* _BEAM_CORE_H_ */
