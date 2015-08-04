#ifndef _MC_MLC_H_
#define _MC_MLC_H_

/*****************************************************************************
 * MC_mlc.h:                                                                 *
 *    class declarations and inline member functions for:                    *
 *       MC_mlc:         multi-leaf collimator (MLC) base class              *
 *                                                                           *
 * Copyright (C) 2001    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 03.12.2001      *
 *    dynamic MLC                                         MF 23.10.2003      *
 *                                                                           *
 *****************************************************************************/

#include "definitions.h"
#include "multi_leaf.h"
#include "MC_object.h"
#include "MC_material.h"

// ****************************************************************
// class MC_mlc: multi-leaf collimator (MLC) base class
// ****************************************************************

class MC_mlc : public MC_object
{
   protected:
      // initialize MLC
      void init(multi_leaf *ini_mlc,
                const real &ini_x_min, const real &ini_x_max,
                const real &ini_y_min, const real &ini_y_max,
                const real &ini_iso_dist,
                const real &ini_z_focus_wall,
                const real &ini_z_focus_edge,
                const char *open_material_name);

   public:
      // define MLC using the nominal MLC, the opening material name,
      // the source to iso-center distance as well as
      // the numbers of object planes and regions
      MC_mlc(multi_leaf *ini_mlc,      const char *open_material_name,
             const real  ini_iso_dist,
             unsigned    n_planes,     unsigned    n_regions);

      // delete MLC
      virtual ~MC_mlc(void);

      // get MLC type
      mlc_type get_type(void) { return(mlctype); }

      // get MLC mode
      mlc_mode get_mode(void) { return(mlcmode); }

      // set MLC mode
      void set_mode(mlc_mode new_mode) {
         mlcmode = new_mode; nominal_mlc->set_mode(new_mode); }

      // get leaf moving direction X or Y
      char get_xy(void) {
         return( (xytype == X) ? 'X' : ( (xytype == Y) ? 'Y' : 'Z') ); }

      // get number of leaf pairs
      int get_num(void) {  return(num_pairs); }

      // get starting position of left leaf i
      real get_start_left(int i) { return(nominal_mlc->get_start_left(i)); }

      // get starting position of right leaf i
      real get_start_right(int i) { return(nominal_mlc->get_start_right(i)); }

      // get stopping position of left leaf i
      real get_stop_left(int i) { return(nominal_mlc->get_stop_left(i)); }

      // get stopping position of right leaf i
      real get_stop_right(int i) { return(nominal_mlc->get_stop_right(i)); }

      // estimate the region index of a point
      virtual int estimate_region_index(const real_3 &)=0;

      // change position of one leaf pair,
      // the new leaf positions are defined at the iso-center plane,
      // this base class function only performs checks and changes the
      // nominal MLC data, the real leaf positions must be changed by
      // the corresponding derived class member function
      virtual void change_pair(const int  &pair_index,
                               const real &new_left,
                               const real &new_right);

      // change starting position of one leaf pair (iso-center plane),
      // only the nominal MLC data will be changed
      void change_start_pair(const int  &pair_index,
                             const real &new_start_left,
                             const real &new_start_right);

      // change stopping position of one leaf pair (iso-center plane),
      // only the nominal MLC data will be changed
      void change_stop_pair(const int  &pair_index,
                            const real &new_stop_left,
                            const real &new_stop_right);

   protected:
      // MLC type (SIMPLE_MLC, DBLFOCUS_MLC, RNDFOCUS_MLC, ELEKTA_MLC,
      // VARIAN_MLC, etc.)
      mlc_type  mlctype;

      // MLC mode (STATIC_MLC, DYNAMIC_MLC)
      mlc_mode  mlcmode;

      // X or Y MLC (leaf moving direction)
      axis      xytype;

      // cross section data bases for leaf and opening materials
      MC_material *leaf_material;
      MC_material *open_material;

      // the outer object dimensions
      real   x_min,x_max,y_min,y_max,z_min,z_max;

      // origin to iso-center distance (cm)
      real   iso_distance;

      // z positions of the two focus points
      real   z_focus_wall, z_focus_edge;

      // number of leaf pairs
      int    num_pairs;

      // pointer to the nominal MLC data
      multi_leaf *nominal_mlc;

      // plane pointers (outer object planes)
      MC_plane *plane_x_min, *plane_x_max, *plane_y_min, *plane_y_max,
               *plane_z_min, *plane_z_max;

      // array of pointers to the left and right leaf edge planes
      MC_plane **left_edge_planes, **right_edge_planes;
};

#endif /* _MC_MLC_H_ */
