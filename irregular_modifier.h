#ifndef _IRREGULAR_MODIFIER_H_
#define _IRREGULAR_MODIFIER_H_

/*****************************************************************************
 * irregular_modifier.h:                                                     *
 *    class declarations and inline member functions for:                    *
 *       irregular_modifier:  modifier based on an irregular beam contour    *
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

#include "definitions.h"
#include "ranmar.h"

#include "contour.h"
#include "beam_modifier.h"

// modifier based on an irregular beam contour
class irregular_modifier : public beam_modifier
{
   public:
      // create irregular beam modifier based on the irregular beam
      // contour at the iso-center plane
      irregular_modifier(contour *, real);

      // create irregular beam modifier based on the irregular beam
      // contour at the iso-center plane and the position of this contour
      irregular_modifier(contour *, real, real);

      // delete irregular beam modifier
      ~irregular_modifier(void) { delete lower_contour; }

      // transport particle through the irregular beam modifier,
      // if the particle is deleted return false, otherwise return true
      bool transport(particle_parameters &, ranmar &);

      // transport position and direction through the irregular beam modifier,
      // if the particle is deleted return false, otherwise return true
      bool transport(real_3 &, real_3 &, ranmar &);

   protected:
      // irregular beam contour at the lower modifier plane
      contour *lower_contour;
};

#endif  /* _IRREGULAR_MODIFIER_H_ */
