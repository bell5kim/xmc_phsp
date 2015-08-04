/*****************************************************************************
 * treatment_plan.cpp:                                                       *
 *    class member functions for:                                            *
 *       treatment_plan: treatment plan parameters                           *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 99/12/14        *
 *    added portal dose image                             MF 03/03/13        *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <stdlib.h>
#include <string.h>
#include "definitions.h"
#include "global.h"
#include "electron_data_inp.h"
#include "photon_data_inp.h"
#include "electron_data.h"
#include "photon_data.h"
#include "treatment_plan.h"

// ****************************************
// member functions of class treatment_plan
// ****************************************

// construct empty treatment plan
treatment_plan::treatment_plan(void)
{
   patient_id       = NULL;
   plan_id          = NULL;
   calculate        = false;
   out3d_beam       = false;
   out3d_plan       = false;
   first_plane      = NULL;
   last_plane       = NULL;
   first_profile    = NULL;
   last_profile     = NULL;
   total_number     = 0;
   sum_weight       = 0.0;
   sum_monunits     = 0.0;
   num_fractions    = 1;
   result           = DOSE;
   film_factor      = 1.0;
   norm             = NONE;
   norm_point.x     = 0.0;
   norm_point.y     = 0.0;
   norm_point.z     = 0.0;
   cpu_time         = 0.0;
   dose_max         = 0.0;
   d_max.x          = 0;
   d_max.y          = 0;
   d_max.z          = 0;
   ini_rndm1        = 12;
   ini_rndm2        = 34;
   ini_rndm3        = 56;
   ini_rndm4        = 78;
   portal           = NULL;
   first_beam       = NULL;
   last_beam        = NULL;
}

// free memory
treatment_plan::~treatment_plan()
{
   extern electron_data_inp *e_inp; // input electron transport data
   extern photon_data_inp   *p_inp; // input photon transport data

   delete e_inp; e_inp = NULL;
   delete p_inp; p_inp = NULL;

   if (portal != NULL) { delete portal; portal = NULL; }
}

// initialize treatment plan without beams by patient name and plan name
void treatment_plan::init(const char *patient_name, const char *plan_name)
{
   patient_id       = strdup(patient_name);
   plan_id          = strdup(plan_name);
   calculate        = false;
   out3d_beam       = false;
   out3d_plan       = false;
   first_plane      = NULL;
   last_plane       = NULL;
   first_profile    = NULL;
   last_profile     = NULL;
   total_number     = 0;
   sum_weight       = 0.0;
   sum_monunits     = 0.0;
   num_fractions    = 1;
   result           = DOSE;
   film_factor      = 1.0;
   norm             = NONE;
   norm_point.x     = 0.0;
   norm_point.y     = 0.0;
   norm_point.z     = 0.0;
   cpu_time         = 0.0;
   dose_max         = 0.0;
   d_max.x          = 0;
   d_max.y          = 0;
   d_max.z          = 0;
   ini_rndm1        = 12;
   ini_rndm2        = 34;
   ini_rndm3        = 56;
   ini_rndm4        = 78;
   portal           = NULL;
   first_beam       = NULL;
   last_beam        = NULL;
}

// add new beam on bottom of the list (FIFO)
void treatment_plan::add_beam(real weight)
{
   beam_core    *new_beam;                   // pointer to the new beam

   new_beam                  = new beam_core; // create new beam
   new_beam->id              = total_number;  // beam id
   new_beam->weight          = weight;        // set beam weight
   // first beam is on top of the list
   if (first_beam == NULL) first_beam = new_beam;
   // the old last beam becomes the last but one
   if (last_beam != NULL) last_beam->next_beam_core = new_beam;
   // the new beam is placed on bottom of the list
   new_beam->previous_beam_core  = last_beam;
   // the new beam becomes the last beam
   last_beam                = new_beam;

   ++total_number;
   sum_weight += new_beam->weight;
}

// delete first beam from list (FIFO)
void treatment_plan::erase_first(void)
{
   beam_core            *old_beam; // pointer to the old first beam
   extern electron_data *e_h2o;    // electron transport data for this beam
   extern photon_data   *p_h2o;    // photon transport data for this beam

   if (first_beam == NULL)  // there isn't any beam to delete
   {
      xvmc_error("treatment_plan::erase_first",
                 "there isn't any beam to erase",8);
   }
   else
   {
      old_beam   = first_beam;             // the beam to delete
      first_beam = old_beam->next_beam_core;    // the new first beam
      if (first_beam == NULL)              // there isn't any further beam
      {
         last_beam = NULL;
      }
      else
      {
         first_beam->previous_beam_core = NULL;
      }

      sum_weight -= old_beam->weight;
      --total_number;
      // free memory
      delete old_beam;
      delete e_h2o; e_h2o = NULL;
      delete p_h2o; p_h2o = NULL;
   }
}

// normalize beam weights
void treatment_plan::normalize_weights(void)
{
   beam_core *this_beam = first_beam;
   while (this_beam != NULL)
   {
      this_beam->weight /= sum_weight;
      this_beam = this_beam->next();
   }
#ifdef MONITORBSF   
   if (result == ABS_DOSE || result == ABS_DOSE_MBSF) sum_monunits = sum_weight;
#else
   if (result == ABS_DOSE) sum_monunits = sum_weight;
#endif      
   sum_weight = ONE;
}
