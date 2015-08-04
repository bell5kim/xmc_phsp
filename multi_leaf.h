#ifndef _MULTI_LEAF_H_
#define _MULTI_LEAF_H_

/*****************************************************************************
 * multi_leaf.h:                                                             *
 *    class declarations and inline member functions for:                    *
 *       class multi_leaf: nominal multi-leaf collimator (iso-center plane)  *
 *                                                                           *
 * Copyright (C) 2001    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 15.06.2001      *
 *    dynamic MLC mode implemented                        MF 09.10.2003      *
 *                                                                           *
 *****************************************************************************/

#include "definitions.h"
#include "global.h"

enum mlc_type {SIMPLE_MLC, DBLFOCUS_MLC, RNDFOCUS_MLC, ELEKTA_MLC, VARIAN_MLC};
enum mlc_mode {STATIC_MLC, DYNAMIC_MLC};

// nominal multi-leaf collimator (iso-center plane)
class multi_leaf
{
   public:
      // define MLC by type, mode, leaf moving direction (X or Y),
      // number of leaf pairs, MLC material,
      // upper and lower limits of the MLC,
      // center z-position and radius of leaf curvature (if there is any)
      multi_leaf(const mlc_type  inp_type, const mlc_mode  inp_mode,
                 const char      xy,
                 const int       num,      const char     *inp_material,
                 const real     *z_upper,  const real     *z_lower,
                 const real     *z_curve,  const real     *r_curve) {
         init(inp_type,inp_mode,xy,num,inp_material,
              z_upper,z_lower,z_curve,r_curve); }

      // delete MLC
      ~multi_leaf(void);

      // reset leaf positions (zero beam opening)
      void reset(void);

      // change width of one leaf pair,
      // if this is impossible -> return false
      bool change_width(int, real);

      // change position of one leaf pair,
      // if this is impossible -> return false
      bool change_pair(int, real, real);

      // change starting position of one leaf pair,
      // if this is impossible -> return false
      bool change_start_pair(int, real, real);

      // change stopping position of one leaf pair,
      // if this is impossible -> return false
      bool change_stop_pair(int, real, real);

      // get MLC type
      mlc_type get_type(void) { return(type); }

      // get MLC mode
      mlc_mode get_mode(void) { return(mode); }

      // set MLC mode
      void set_mode(mlc_mode new_mode) { mode = new_mode; }

      // get leaf moving direction X or Y
      char get_xy(void) { return(x_leafs ? 'X' : 'Y'); }

      // get MLC material
      char *get_material(void) { return(material); }

      // get upper MLC limit
      real *get_upper(void) { return(upper_limit); }

      // get lower MLC limit
      real *get_lower(void) { return(lower_limit); }

      // get center z-position of leaf curvature
      real *get_z_center_curve(void) { return(z_center_curve); }

      // get radius of leaf curvature
      real *get_radius_curve(void) { return(radius_curve); }

      // get number of leaf pairs
      int get_num(void) {  return(num_pairs); }

      real get_width(int);       // get width of leaf i
      real get_lperp(int);       // get perpendicular leaf limit i
      real get_left(int);        // get position of left leaf i
      real get_right(int);       // get position of right leaf i
      real get_start_left(int);  // get starting position of left leaf i
      real get_start_right(int); // get starting position of right leaf i
      real get_stop_left(int);   // get stopping position of left leaf i
      real get_stop_right(int);  // get stopping position of right leaf i

      // for a point (x0,y0) get the corresponding leaf pair index,
      // return -1 if the point is outside the MLC
      int get_index(real, real);

      // return true if a point (x0,y0) is located inside the beam limits
      bool inside(real, real);

   private:
      // MLC type (SIMPLE_MLC, DBLFOCUS_MLC, RNDFOCUS_MLC, ELEKTA_MLC,
      // VARIAN_MLC, etc.)
      mlc_type  type;

      // MLC modality (static or dynamic)
      mlc_mode  mode;

      // leaf moving direction: true  if the leafs are parallel to the X axis
      //                        false if the leafs are parallel to the Y axis
      bool      x_leafs;

      // MLC material
      char *material;

      // upper and lower MLC limits
      real *upper_limit, *lower_limit;

      // center z-position and radius of leaf curvature
      real *z_center_curve, *radius_curve;

      int       num_pairs;       // number of leaf pairs
      real     *width;           // array of leaf widths
      real     *lperp;           // array of perpendicular leaf limits
      real     *left;            // array of left leaf positions
      real     *right;           // array of right leaf positions
      real      min_left;        // minimum left  leaf position (max. opening)
      real      max_right;       // maximum right leaf position (max. opening)
      real     *start_left;      // array of left  leaf starting positions
      real     *start_right;     // array of right leaf starting positions
      real      min_start_left;  // minimum left  leaf starting position
      real      max_start_right; // maximum right leaf starting position
      real     *stop_left;       // array of left  leaf stopping positions
      real     *stop_right;      // array of right leaf stopping positions
      real      min_stop_left;   // minimum left  leaf stopping position
      real      max_stop_right;  // maximum right leaf stopping position

      // minimum and maximum beam openings perpendicular to the leaf direction
      real      min_perp,max_perp;

      // initialize MLC by type, mode, leaf moving direction (X or Y),
      // number of leaf pairs, MLC material,
      // upper and lower limits of the MLC,
      // center z-position and radius of leaf curvature (if there is any)
      void init(const mlc_type, const mlc_mode,
                const char,
                const int,      const char *,
                const real *,   const real *,
                const real *,   const real *);
};

// for a point (x0,y0) get the corresponding leaf pair index,
// return -1 if the point is outside the MLC
inline int multi_leaf::get_index(real x0, real y0)
{
   // initialize index
   int lower = 0;
   int upper = num_pairs;

   if (x_leafs) // X leafs
   {
      // test, whether (x0,y0) is outside the MLC
      if (y0 < min_perp)  return(-1);
      if (y0 > max_perp)  return(-1);

      // calculate leaf pair number
      while ( upper > lower+1 )
      {
         int m = (lower+upper)/2;
         if (y0 < lperp[m]) upper = m;
         else               lower = m;
      }
   }
   else         // Y leafs
   {
      // test, whether (x0,y0) is outside the MLC
      if (x0 < min_perp)  return(-1);
      if (x0 > max_perp)  return(-1);

      // calculate leaf pair number
      while ( upper > lower+1 )
      {
         int m = (lower+upper)/2;
         if (x0 < lperp[m]) upper = m;
         else               lower = m;
      }
   }

   return(upper-1);
}

// return true if a point (x0,y0) is located inside the beam limits
inline bool multi_leaf::inside(real x0, real y0)
{
   if (x_leafs) // X leafs
   {
      // test, whether (x0,y0) is outside the rectangle defined by
      // min_left,max_right,min_perp,max_perp (to save compuation time)
      if (x0 < min_left)  return(false);
      if (x0 > max_right) return(false);
      if (y0 < min_perp)  return(false);
      if (y0 > max_perp)  return(false);

      // calculate leaf pair number
      int lower = 0;
      int upper = num_pairs;
      while ( upper > lower+1 )
      {
         int m = (lower+upper)/2;
         if (y0 < lperp[m]) upper = m;
         else               lower = m;
      }
      int i_pair = upper - 1;

      // check whether left[i] < x0 < right[i]
      if ( (left[i_pair] < x0) && (x0 < right[i_pair]) ) return(true);
   }
   else         // Y leafs
   {
      // test, whether (x0,y0) is outside the rectangle defined by
      // min_left,max_right,min_perp,max_perp (to save compuation time)
      if (x0 < min_perp)  return(false);
      if (x0 > max_perp)  return(false);
      if (y0 < min_left)  return(false);
      if (y0 > max_right) return(false);

      // calculate leaf pair number
      int lower = 0;
      int upper = num_pairs;
      while ( upper > lower+1 )
      {
         int m = (lower+upper)/2;
         if (x0 < lperp[m]) upper = m;
         else               lower = m;
      }
      int i_pair = upper - 1;

      // check whether left[i] < y0 < right[i]
      if ( (left[i_pair] < y0) && (y0 < right[i_pair]) ) return(true);
   }

   // the point (x0,y0) is outside
   return(false);
}

#endif /* _MULTI_LEAF_H_ */
