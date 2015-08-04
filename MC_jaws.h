#ifndef _MC_JAWS_H_
#define _MC_JAWS_H_

/*****************************************************************************
 * MC_jaws.h:                                                                *
 *    class declarations and inline member functions for:                    *
 *       MC_jaws:        physical collimator (JAW) base class                *
 *                                                                           *
 * Copyright (C) 2001    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 03.12.2001      *
 *    dynamic MLC                                         MF 23.10.2003      *
 *    Adopted from MC_mlc.h and MC_jaws_focus.h           JK Feb 22, 2008    *
 *                                                                           *
 *****************************************************************************/

#include "definitions.h"
#include "MC_volume_6p.h"
#include "MC_object.h"
#include "MC_material.h"

enum jaw_type {SIMPLE_JAW, DBLFOCUS_JAW, REAL_JAW};
enum jaw_mode {STATIC_JAW, DYNAMIC_JAW};

// ****************************************************************
// class MC_jaw: multi-leaf collimator (MLC) base class
// ****************************************************************

class MC_jaw : public MC_object
{
   protected:
      // initialize JAW
      void init(const axis &ini_type, 
                const jaw_type   &int_jawtype,  const jaw_mode   &ini_jawmode,
                const real &ini_left,           const real &ini_right,
                const real &ini_x_min,          const real &ini_x_max,
                const real &ini_y_min,          const real &ini_y_max,
                const real &ini_z_min,          const real &ini_z_max,
                const real &ini_iso_dist,       const real &ini_z_focus,
                const char *bars_material_name, 
                const char *open_material_name);

   public:
      // define JAW, the opening material name,
      // the source to iso-center distance as well as
      // the numbers of object planes and regions
      MC_jaw(const axis &ini_type, 
             const jaw_type   &int_jawtype,  const jaw_mode   &ini_jawmode,
             const real &ini_left,  const real &ini_right,
             const real &ini_z_min, const real &ini_z_max,
             const char *bars_material,
             const char *open_material);
      // delete MLC
      virtual ~MC_jaw(void);

      // get JAW type
      jaw_type get_type(void) { return(jawtype); }

      // get JAW mode
      jaw_mode get_mode(void) { return(jawmode); }

      // set JAW mode
      void set_mode(jaw_mode new_mode) {
         jawmode = new_mode; set_mode(new_mode); }

      // get jaw moving direction X or Y
      char get_xy(void) {
         return( (xytype == X) ? 'X' : ( (xytype == Y) ? 'Y' : 'Z') ); }

      // get starting position of left jaw
      real get_start_left() { return(get_start_left()); }

      // get starting position of right jaw
      real get_start_right() { return(get_start_right()); }

      // get stopping position of left jaw
      real get_stop_left() { return(get_stop_left()); }

      // get stopping position of right jaw
      real get_stop_right() { return(get_stop_right()); }

      // estimate the region index of a point
      int estimate_region_index(const real_3 &);

      // change left jaw position
      bool change_left(const real &new_left);

      // change right jaw position
      bool change_right(const real &new_right);

      // change position of jaw,
      // the new jaw position is defined at the iso-center plane,
      // this base class function only performs checks and changes the
      // nominal JAW data, the real jaw positions must be changed by
      // the corresponding derived class member function
      virtual bool change_opening(const real &new_left,
                                  const real &new_right);

      // change starting position of jaw pair (iso-center plane),
      // only the nominal JAW data will be changed
      bool change_start_pair(const real &new_start_left,
                             const real &new_start_right);

      // change stopping position of one leaf pair (iso-center plane),
      // only the nominal JAW data will be changed
      bool change_stop_pair(const real &new_stop_left,
                            const real &new_stop_right);

      // get upper boundary 
      real get_upper(void) { return(z_min); }

      // get lower boundary
      real get_lower(void) { return(z_max); }

      // get open left 
      real get_open_left(void) { return(open_left); }

      // get open right 
      real get_open_right(void) { return(open_right); }

   protected:
      // JAW type (SIMPLE_JAW, DBLFOCUS_JAW, REAL_JAW)
      jaw_type  jawtype;

      // JAW mode (STATIC_JAW, DYNAMIC_JAW)
      jaw_mode  jawmode;

      // X or Y JAW (jaw moving direction)
      axis      xytype;

      // cross section data bases for bars and opening materials
      MC_material *bars_material;
      MC_material *open_material;

      // the outer object dimensions
      real   x_min,x_max,y_min,y_max,z_min,z_max;

      // origin to iso-center distance (cm)
      real   iso_distance;

      // z positions of the focus
      real   z_focus;

      // left and right jaw positions projected to the iso center plane (cm)
      real   open_left,open_right;

      // plane pointers (outer object planes)
      MC_plane *plane_x_min, *plane_x_max, *plane_y_min, *plane_y_max,
               *plane_z_min, *plane_z_max;

      // array of pointers to the left and right leaf edge planes
      MC_plane *plane_left, *plane_right;

};

// estimate the region index of a point
inline int MC_jaw::estimate_region_index(const real_3 &p0)
{
   if (xytype == X)
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
#endif /* _MC_JAWS_H_ */
