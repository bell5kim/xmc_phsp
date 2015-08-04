#ifndef _TREATMENT_PLAN_H_
#define _TREATMENT_PLAN_H_

/*****************************************************************************
 * treatment_plan.h:                                                         *
 *    class declarations and inline member functions for:                    *
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

#include "definitions.h"
#include "global.h"
#include "beam_core.h"
#include "portal_dose.h"

class treatment_plan
{
   public:
      char         *patient_id;      // patient name
      char         *plan_id;         // treatment plan name
      bool          calculate;       // calculate or read dose distribution

      bool          out3d_beam;      // write 3D dose matrix for each beam
      bool          out3d_plan;      // write 3D dose matrix for the whole plan

      // linked list of parameters to plot xy, xz and yz planes
      plane_parameters  *first_plane;
      plane_parameters  *last_plane;

      // linked list of parameters to plot x, y and z profiles
      profile_parameters  *first_profile;
      profile_parameters  *last_profile;

      // not implemented yet         // change electron density
      // not implemented yet         // scale X-section
      // not implemented yet         // analyse region

      int           total_number;    // total number of beams
      real          sum_weight;      // sum of beam weights

      // number of monitor units per fraction for all beams
      real          sum_monunits;
      //
      // the total number of monitor units is: "sum_monunits * num_fractions"
      //

      //  number of dose fractions (default: 1)
      int           num_fractions;

      // simulation result type:
      // DOSE, DOSE_TO_H2O, IONIZATION, FILM, ABS_DOSE, ABD_DOSE_MBSF
      result_type   result;

      real          film_factor;         // factor for film dose
      norm_type     norm;                // normalization type
                                         // 100% -> NONE, MAXIMUM, VALUE, POINT
      real_3        norm_point;          // normalization point
      real          cpu_time;            // total CPU time
      real          dose_max;            // total dose maximum
      int_3         d_max;               // indices of the dose maximum voxel
      int           ini_rndm1,ini_rndm2, // integers to initialize the
                    ini_rndm3,ini_rndm4; // random number generator

      // portal dose image
      portal_dose  *portal;

      beam_core    *first_beam;          // pointer to the first beam
      beam_core    *last_beam;           // pointer to the last beam

      // construct empty treatment plan
      treatment_plan(void);
      // free memory
      ~treatment_plan();
      // for patient name initialize treatment plan without beams
      void init(const char *, const char *);
      // add new beam on bottom of the list (FIFO)
      void add_beam(real);
      // delete first beam from list (FIFO)
      void erase_first(void);
      // normalize beam weights
      void normalize_weights(void);
};

#endif  /* _TREATMENT_PLAN_H_ */
