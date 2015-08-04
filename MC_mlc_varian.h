#ifndef _MC_MLC_VARIAN_H_
#define _MC_MLC_VARIAN_H_

/*****************************************************************************
 * MC_mlc_varian.h:                                                          *
 *    class declaration and inline member functions for:                     *
 *       MC_mlc_varian:  MLC for VARIAN medical linear accelerators          *
 *                                                                           *
 * Copyright (C) 2001    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding (adapted from ELEKTA MLC)            MF 10.07.2003      *
 *                                                                           *
 *****************************************************************************/

#include "definitions.h"
#include "MC_halfcyl.h"
#include "MC_volume_6p.h"
#include "MC_mlc.h"

// ****************************************************************
// class MC_mlc_varian: MLC for VARIAN medical linear accelerators
// ****************************************************************

class MC_mlc_varian : public MC_mlc
{
   private:
      // initialize VARIAN MLC
      void init(void);

   public:
      // define VARIAN MLC using the nominal MLC, the opening material name
      // and the source to iso-center distance
      MC_mlc_varian(multi_leaf *ini_mlc, const char *open_material,
                    const real  ini_iso_dist)
         : MC_mlc(ini_mlc, open_material, ini_iso_dist,
                  6*ini_mlc->get_num()+14, 20*ini_mlc->get_num()+11) { init(); }

      // delete VARIAN MLC
      ~MC_mlc_varian(void) {;}

      // estimate the region index of a point
      int estimate_region_index(const real_3 &);

      // change position of one leaf pair,
      // the new leaf positions are defined at the iso-center plane
      void change_pair(const int &, const real &, const real &);

   private:
      // center z-position and radius of leaf curvature
      real   z_center_curve, radius_curve;

      // upper and lower z-positions for tongue and groove
      real   z_tongue_upper, z_tongue_lower;
      real   z_groove_upper, z_groove_lower;

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
      MC_plane *plane_z_tongue_upper, *plane_z_tongue_lower;
      MC_plane *plane_z_groove_upper, *plane_z_groove_lower;
};

// estimate the region index of a point
inline int MC_mlc_varian::estimate_region_index(const real_3 &p0)
{
   // scale x and y to the iso-center plane
   real x0 = p0.x*iso_distance/p0.z;
   real y0 = p0.y*iso_distance/p0.z;

   return( 20*nominal_mlc->get_index(x0,y0)+18 );
}

#endif /* _MC_MLC_VARIAN_H_ */
