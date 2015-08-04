#ifndef _MC_JAWS_FOCUS_H_
#define _MC_JAWS_FOCUS_H_

/*****************************************************************************
 * MC_jaws_focus.h:                                                          *
 *    class declarations and inline member functions for:                    *
 *       MC_jaws_focus:  a pair of focussing jaws                            *
 *                                                                           *
 * Copyright (C) 2001    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 29.08.2001      *
 *                                                                           *
 *****************************************************************************/

#include "definitions.h"
#include "MC_volume_6p.h"
#include "MC_object.h"

// ****************************************************************
// class MC_jaws_focus: a pair of focussing jaws
// ****************************************************************

class MC_jaws_focus : public MC_object
{
   private:
      // initialize jaws pair
      void init(const axis &,
                const real &, const real &,
                const real &, const real &,
                const real &, const real &,
                const real &, const real &,
                const real &, const real &,
                const char *, const char *);

   public:
      // define x or y jaws pair with given opening, z boundaries,
      // jaws bar material and opening material names
      MC_jaws_focus(const axis &, const real &, const real &,
                                  const real &, const real &,
                                  const char *, const char *);

      // delete focussing jaws
      ~MC_jaws_focus(void);

      // estimate the region index of a point
      int estimate_region_index(const real_3 &);

      // change left jaw position
      void change_left(const real &);

      // change right jaw position
      void change_right(const real &);

      // change left and right jaw positions
      void change_opening(const real &new_left, const real &new_right)
         { change_left(new_left); change_right(new_right); }

   private:
      // X or Y jaws
      axis   type;

      // ESTAR material electron transport data file paths
      char *bars_material_estar;
      char *open_material_estar;

      // NIST material cross section file paths
      char *bars_material_nist;
      char *open_material_nist;

      // differential Compton and pair cross section file paths
      char *compton_file;
      char *pair_file;

      // the outer object dimensions
      real   x_min,x_max,y_min,y_max,z_min,z_max;

      // origin to iso-center distance (cm)
      real   iso_distance;

      // z position of the focus
      real   z_focus;

      // left and right jaws positions projected to the iso center plane (cm)
      real   open_left,open_right;

      // plane pointers
      MC_plane *plane_x_min, *plane_x_max, *plane_y_min, *plane_y_max,
               *plane_z_min, *plane_z_max, *plane_left,  *plane_right;
};

// estimate the region index of a point
inline int MC_jaws_focus::estimate_region_index(const real_3 &p0)
{
   if (type == X)
   {
      // scale x to the iso-center plane
      real x0 = p0.x*iso_distance/p0.z;

      // determine region
      if (x0 < open_left)  return(0);
      if (x0 > open_right) return(1);
   }
   else
   {
      // scale y to the iso-center plane
      real y0 = p0.y*iso_distance/p0.z;

      // determine region
      if (y0 < open_left)  return(0);
      if (y0 > open_right) return(1);
   }

   // return intermediate region
   return(2);
}

#endif /* _MC_JAWS_FOCUS_H_ */
