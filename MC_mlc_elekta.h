#ifndef _MC_MLC_ELEKTA_H_
#define _MC_MLC_ELEKTA_H_

/*****************************************************************************
 * MC_mlc_elekta.h:                                                          *
 *    class declaration and inline member functions for:                     *
 *       MC_mlc_elekta:  MLC for ELEKTA medical linear accelerators          *
 *                                                                           *
 * Copyright (C) 2001    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 29.11.2001      *
 *    implementation of interleaf leakage                 MF 03.03.2003      *
 *                                                                           *
 *****************************************************************************/

#include "definitions.h"
#include "MC_halfcyl.h"
#include "MC_volume_6p.h"
#include "MC_mlc.h"

// ****************************************************************
// class MC_mlc_elekta: MLC for ELEKTA medical linear accelerators
// ****************************************************************

class MC_mlc_elekta : public MC_mlc
{
   private:
      // initialize ELEKTA MLC
      void init(void);

   public:
      // define ELEKTA MLC using the nominal MLC, the opening material name
      // and the source to iso-center distance
      MC_mlc_elekta(multi_leaf *ini_mlc, const char *open_material,
                    const real  ini_iso_dist)
         : MC_mlc(ini_mlc, open_material, ini_iso_dist,
                  6*ini_mlc->get_num()+12, 12*ini_mlc->get_num()+7) { init(); }

      // delete ELEKTA MLC
      ~MC_mlc_elekta(void) {;}

      // estimate the region index of a point
      int estimate_region_index(const real_3 &);

      // change position of one leaf pair,
      // the new leaf positions are defined at the iso-center plane
      void change_pair(const int &, const real &, const real &);

   private:
      // center z-position and radius of leaf curvature
      real   z_center_curve, radius_curve;

      // z-positions for tongue and groove
      real   z_tongue, z_groove;

      // tonge and groove width at the iso-center
      real   width_tg;

      // width of the interleaf leakage gap
      real   width_gap;

      // distance of the perpendicular leaf shift
      real   d_shift;

      // shift distance of the leaf side focus point
      real   f_shift;

      // corresponding shift at the iso-center (the shift of the leaf
      // side focus point is caused by tilting the whole MLC)
      real   i_shift;

      // tongue and groove planes
      MC_plane *plane_z_tongue, *plane_z_groove;
};

// estimate the region index of a point
inline int MC_mlc_elekta::estimate_region_index(const real_3 &p0)
{
   // scale x and y to the iso-center plane
   real x0 = p0.x*iso_distance/p0.z;
   real y0 = p0.y*iso_distance/p0.z;

   return( 12*nominal_mlc->get_index(x0,y0)+11 );
}

#endif /* _MC_MLC_ELEKTA_H_ */
