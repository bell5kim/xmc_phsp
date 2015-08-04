#ifndef _REAL_MLC_H_
#define _REAL_MLC_H_

/*****************************************************************************
 * real_jaw.h:                                                               *
 *    class declarations and inline member functions for:                    *
 *       real_jaw:            real 3D Monte Carlo JAW model                  *
 *                                                                           *
 * Copyright (C) 2002    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 31.01.2002      *
 *    real dynamic MLC                                    MF 23.10.2003      *
 *    Adoped from real_mlc.h                              JK Feb 22, 2008    *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include "definitions.h"
#include "ranmar.h"

#include "beam_modifier.h"
#include "MC_jaws.h"

// real 3D Monte Carlo JAW model
class real_jaw : public beam_modifier
{
   public:
      // create JAW based on the jaw positions at the iso-center plane
      real_jaw(MC_jaw *, real);

      // delete MLC
      ~real_jaw(void) { delete jaw; }

      // set beam modifier and JAW modes using a new beam modifier mode
      void set_mode(modifier_mode new_mode);

      // set beam modifier and JAW modes using a new JAW mode
      void set_mode(jaw_mode new_mode);

      // change position of one jaw pair,
      // the new jaw opening position is defined at the iso-center plane
      void change_opening(const real &new_left, const real &new_right) {
         jaw->change_opening(new_left, new_right); }

      // change left position of jaw,
      // the new jaw position is defined at the iso-center plane
      void change_left(const real &new_left) {
         jaw->change_left(new_left); }

      // change right position of jaw,
      // the new jaw position is defined at the iso-center plane
      void change_right(const real &new_right) {
         jaw->change_right(new_right); }

      // adjust real dynamic JAW, in static mode do nothing
      void adjust(long i_history, long n_history);

      // transport particle through jaws,
      // if the particle is deleted return false, otherwise return true
      bool transport(particle_parameters &, ranmar &);

      // transport particle position and direction through the JAW,
      // always return false because this function is not implemented
      bool transport(real_3 &pos, real_3 &dir, ranmar &rndm) { return(false); }

   protected:
      // the JAW object
      MC_jaw *jaw;
};

#endif  /* _REAL_MLC_H_ */
