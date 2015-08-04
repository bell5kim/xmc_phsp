#ifndef _REAL_MLC_H_
#define _REAL_MLC_H_

/*****************************************************************************
 * real_mlc.h:                                                               *
 *    class declarations and inline member functions for:                    *
 *       real_mlc:            real 3D Monte Carlo MLC model                  *
 *                                                                           *
 * Copyright (C) 2002    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 31.01.2002      *
 *    real dynamic MLC                                    MF 23.10.2003      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include "definitions.h"
#include "ranmar.h"

#include "multi_leaf.h"
#include "beam_modifier.h"
#include "MC_mlc.h"
#include "MC_mlc_2focus.h"
#include "MC_mlc_rfocus.h"
#include "MC_mlc_elekta.h"
#include "MC_mlc_varian.h"

// real 3D Monte Carlo MLC model
class real_mlc : public beam_modifier
{
   public:
      // create MLC based on the leaf positions at the iso-center plane
      real_mlc(multi_leaf *, real);

      // delete MLC
      ~real_mlc(void) { delete mlc; }

      // set beam modifier and MLC modes using a new beam modifier mode
      void set_mode(modifier_mode new_mode);

      // set beam modifier and MLC modes using a new MLC mode
      void set_mode(mlc_mode new_mode);

      // change position of one leaf pair,
      // the new leaf positions are defined at the iso-center plane
      void change_pair(const int  &pair_index,
                       const real &new_left,
                       const real &new_right) {
         mlc->change_pair(pair_index, new_left, new_right); }

      // change starting position of one leaf pair,
      // the new leaf positions are defined at the iso-center plane
      void change_start_pair(const int  &pair_index,
                             const real &new_start_left,
                             const real &new_start_right) {
         mlc->change_start_pair(pair_index, new_start_left, new_start_right); }

      // change stopping position of one leaf pair,
      // the new leaf positions are defined at the iso-center plane
      void change_stop_pair(const int  &pair_index,
                            const real &new_stop_left,
                            const real &new_stop_right) {
         mlc->change_stop_pair(pair_index, new_stop_left, new_stop_right); }

      // change leaf positions
      void change_leaf_positions(multi_leaf *);

      // adjust real dynamic MLC, in static mode do nothing
      void adjust(long i_history, long n_history);

      // transport particle through multi-leaf collimator,
      // if the particle is deleted return false, otherwise return true
      bool transport(particle_parameters &, ranmar &);

      // transport particle position and direction through the MLC,
      // always return false because this function is not implemented
      bool transport(real_3 &pos, real_3 &dir, ranmar &rndm) { return(false); }

   protected:
      // the MLC object
      MC_mlc *mlc;
};

#endif  /* _REAL_MLC_H_ */
