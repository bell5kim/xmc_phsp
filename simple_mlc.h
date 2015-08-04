#ifndef _SIMPLE_MLC_H_
#define _SIMPLE_MLC_H_

/*****************************************************************************
 * simple_mlc.h:                                                             *
 *    class declarations and inline member functions for:                    *
 *       simple_mlc:     simple multi-leaf collimator model                  *
 *                                                                           *
 * Copyright (C) 2002    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 31.01.2002      *
 *    simple dynamic MLC                                  MF 09.10.2003      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include "definitions.h"
#include "ranmar.h"

#include "multi_leaf.h"
#include "beam_modifier.h"

// simple multi-leaf collimator model
class simple_mlc : public beam_modifier
{
   public:
      // create MLC based on the leaf positions at the iso-center plane,
      // the source to iso-center distance, the upper and lower MLC limits
      simple_mlc(multi_leaf *iso_mlc, const real ini_iso_distance);

      // delete MLC
      ~simple_mlc(void) { delete lower_mlc; }

      // set beam modifier and MLC modes using a new beam modifier mode
      void set_mode(modifier_mode new_mode);

      // set beam modifier and MLC modes using a new MLC mode
      void set_mode(mlc_mode new_mode);

      // change position of one leaf pair,
      // the new leaf positions are defined at the iso-center plane
      void change_pair(const int  &pair_index,
                       const real &new_left,
                       const real &new_right);

      // change starting position of one leaf pair,
      // the new leaf positions are defined at the iso-center plane
      void change_start_pair(const int  &pair_index,
                             const real &new_start_left,
                             const real &new_start_right);

      // change stopping position of one leaf pair,
      // the new leaf positions are defined at the iso-center plane
      void change_stop_pair(const int  &pair_index,
                            const real &new_stop_left,
                            const real &new_stop_right);

      // change leaf positions
      void change_leaf_positions(multi_leaf *);

      // adjust simple dynamic MLC, in static mode do nothing
      void adjust(long i_history, long n_history);

      // transport particle through multi-leaf collimator,
      // if the particle is deleted return false, otherwise return true
      bool transport(particle_parameters &, ranmar &);

      // transport particle position and direction through the MLC,
      // if the particle is deleted return false, otherwise return true
      bool transport(real_3 &pos, real_3 &dir, ranmar &rndm);

   protected:
      // MLC opening at the lower modifier plane
      multi_leaf *lower_mlc;

      // source to iso-center distance
      real        iso_distance;
};

#endif  /* _SIMPLE_MLC_H_ */
