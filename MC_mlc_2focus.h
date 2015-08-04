#ifndef _MC_MLC_2FOCUS_H_
#define _MC_MLC_2FOCUS_H_

/*****************************************************************************
 * MC_mlc_2focus.h:                                                          *
 *    class declarations and inline member functions for:                    *
 *       MC_mlc_2focus:  double focussing MLC                                *
 *                                                                           *
 * Copyright (C) 2001    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 19.11.2001      *
 *                                                                           *
 *****************************************************************************/

#include "definitions.h"
#include "MC_volume_6p.h"
#include "MC_mlc.h"

// ****************************************************************
// class MC_mlc_2focus: double focussing MLC
// ****************************************************************

class MC_mlc_2focus : public MC_mlc
{
   private:
      // initialize MLC
      void init(void);

   public:
      // define MLC using the nominal MLC, the opening material name
      // and the source to iso-center distance
      MC_mlc_2focus(multi_leaf *ini_mlc, const char *open_material,
                    const real  ini_iso_dist)
         : MC_mlc(ini_mlc, open_material, ini_iso_dist,
                  3*ini_mlc->get_num()+7, 3*ini_mlc->get_num()+2) { init(); }

      // delete MLC
      ~MC_mlc_2focus(void) {;}

      // estimate the region index of a point
      int estimate_region_index(const real_3 &);

      // change position of one leaf pair,
      // the new leaf positions are defined at the iso-center plane
      void change_pair(const int &, const real &, const real &);
};

// estimate the region index of a point
inline int MC_mlc_2focus::estimate_region_index(const real_3 &p0)
{
   // scale x and y to the iso-center plane
   real x0 = p0.x*iso_distance/p0.z;
   real y0 = p0.y*iso_distance/p0.z;

   return( 3*nominal_mlc->get_index(x0,y0)+3 );
}

#endif /* _MC_MLC_2FOCUS_H_ */
