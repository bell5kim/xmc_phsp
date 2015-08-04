#ifndef _BEAM_MODIFIER_H_
#define _BEAM_MODIFIER_H_

/*****************************************************************************
 * beam_modifier.h:                                                          *
 *    class declarations and inline member functions for:                    *
 *       beam_modifier:       beam modifier base class (empty modifier)      *
 *                                                                           *
 * Copyright (C) 2002    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 31.01.2002      *
 *    dynamic modifier mode implemented                   MF 09.10.2003      *
 *    Added MODIFIER_REAL_JAW                             JK Feb 22, 2008    *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include "definitions.h"
#include "ranmar.h"

#include "contour.h"

enum modifier_type { MODIFIER_EMPTY,      MODIFIER_IRREGULAR,
                     MODIFIER_SIMPLE_MLC, MODIFIER_REAL_MLC,
                     MODIFIER_SIMPLE_JAW, MODIFIER_REAL_JAW   };
enum modifier_mode { MODIFIER_STATIC,     MODIFIER_DYNAMIC    };

// beam modifier base class (empty modifier)
class beam_modifier
{
   public:
      // construct (empty) beam modifier, default input and output planes
      beam_modifier(void)
      {
         type           = MODIFIER_EMPTY;
         mode           = MODIFIER_STATIC;
         upper_distance = 30.0;
         lower_distance = 40.0;
      }

      // construct (empty) beam modifier, define input and output planes
      beam_modifier(const real z_upper, const real z_lower)
      {
         type           = MODIFIER_EMPTY;
         mode           = MODIFIER_STATIC;
         upper_distance = z_upper;
         lower_distance = z_lower;
      }

      // delete beam modifier
      virtual ~beam_modifier(void) {;}

      // get beam modifier type
      modifier_type get_modifier_type(void) { return(type); }

      // get beam modifier mode
      modifier_mode get_modifier_mode(void) { return(mode); }

      // adjust dynamic beam modifier, in static mode do nothing
      virtual void adjust(long i_history, long n_history) { return; }

      // transport particle through the (empty) modifier, return true
      virtual bool transport(particle_parameters &p, ranmar &rndm);

      // transport particle position and direction through the
      // (empty) modifier, return true
      virtual bool transport(real_3 &pos, real_3 &dir, ranmar &rndm);

   protected:
      // type of beam modifier
      modifier_type type;

      // mode of beam modifier (static or dynamic)
      modifier_mode mode;

      // origin to beam modifier upper plane distance
      real          upper_distance;

      // origin to beam modifier lower plane distance
      real          lower_distance;
};

#endif  /* _BEAM_MODIFIER_H_ */
