/*****************************************************************************
 * beam_core.cpp:                                                            *
 *    class member functions for:                                            *
 *       beam_core:      parameters for one electron or photon beam          *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 99/12/14        *
 *    added pointer to portal dose image                  MF 03/03/13        *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <stdlib.h>
#include <string.h>
#include "definitions.h"
#include "global.h"
#include "beam_core.h"

// ****************************************
// member functions of class beam_core
// ****************************************

// create empty beam
beam_core::beam_core(void)
{
   id                    = 0;
   type                  = PHOTON;        // photon beam
   model_id              = 0;             // mono-energetic point source
   base_key              = NULL;          // no base data file
   applicator            = NULL;          // no applicator
   weight                = 0.0;
   n_history             = 0;
   n_repeat              = 0;
   n_rotate              = 0;
   n_batch               = 0;
   energy                = 0.0;
   iso_center.x          = 0.0;
   iso_center.y          = 0.0;
   iso_center.z          = 0.0;
   gantry_rotation       = false;
   start_gantry_angle    = 0.0;
   stop_gantry_angle     = 0.0;
   table_angle           = 0.0;
   collimator_angle      = 0.0;
   open_x1               = 0.0;
   open_x2               = 0.0;
   open_y1               = 0.0;
   open_y2               = 0.0;
   irregular             = NULL;
   mlc                   = NULL;
   portal                = NULL;
   jaw                   = NULL;
   cpu_time              = 0.0;

   next_beam_core          = NULL;
   previous_beam_core      = NULL;
}


// delete beam, free memory
beam_core::~beam_core(void)
{
   if (base_key != NULL)   { delete [] base_key;   base_key   = NULL; }
   if (applicator != NULL) { delete [] applicator; applicator = NULL; }
   if (irregular != NULL)  { delete irregular;     irregular  = NULL; }
   if (mlc != NULL)        { delete mlc;           mlc        = NULL; }
   if (jaw != NULL)        { delete jaw;           jaw        = NULL; }

   // we don't free the memory for the portal dose image
   // because it belongs to the treatment_plan
   portal = NULL;
}
