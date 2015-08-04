/*****************************************************************************
 * multi_leaf.cpp:                                                           *
 *    class member functions for:                                            *
 *       multi_leaf: nominal multi-leaf collimator (iso-center plane)        *
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

#include <string.h>
#include <new>
using namespace std;

#include "definitions.h"
#include "global.h"
#include "multi_leaf.h"

// ****************************************
// member functions of class multi_leaf
// ****************************************

// initialize MLC by type, mode, leaf moving direction (X or Y),
// number of leaf pairs, MLC material,
// upper and lower limits of the MLC,
// center z-position and radius of leaf curvature (if there is any)
void multi_leaf::init(const mlc_type  inp_type, const mlc_mode  inp_mode,
                      const char      xy,
                      const int       num,      const char     *inp_material,
                      const real     *z_upper,  const real     *z_lower,
                      const real     *z_curve,  const real     *r_curve)
{
   // set MLC type
   type = inp_type;

   // set MLC mode
   mode = inp_mode;

   // check leaf moving direction (X or Y)
   switch (xy)
   {
   case 'x':    // X leafs
   case 'X':    // X leafs
      x_leafs = true;
      break;
   case 'y':    // Y leafs
   case 'Y':    // Y leafs
      x_leafs = false;
      break;
   default:
      xvmc_error("multi_leaf::init","there are only X or Y MLCs",8);
      break;
   }

   // set MLC material
   material = NULL;
   if (inp_material != NULL)
   {
      if ( (material = new (nothrow) char[strlen(inp_material)+1]) == NULL )
         xvmc_error("multi_leaf::init",
                    "cannot allocate memory for MLC material",8);
      strcpy(material,inp_material);
   }

   // set upper and lower MLC limits
   upper_limit = NULL;
   lower_limit = NULL;
   if (z_upper != NULL)
   {
      if ( (upper_limit = new (nothrow) real) == NULL )
         xvmc_error("multi_leaf::init",
                    "cannot allocate memory for upper limit",8);
      *upper_limit = *z_upper;
   }
   if (z_lower != NULL)
   {
      if ( (lower_limit = new (nothrow) real) == NULL )
         xvmc_error("multi_leaf::init",
                    "cannot allocate memory for lower limit",8);
      *lower_limit = *z_lower;
   }

   // center z-position of leaf curvature
   z_center_curve = NULL;
   if ( (type == RNDFOCUS_MLC) ||
        (type == ELEKTA_MLC)   ||
        (type == VARIAN_MLC)      )
   {
      if (z_curve != NULL)
      {
         if ( (z_center_curve = new (nothrow) real) == NULL )
            xvmc_error("multi_leaf::init",
            "cannot allocate memory for center z-position of leaf curvature",8);
         *z_center_curve = *z_curve;
      }
   }

   // radius of leaf curvature
   radius_curve = NULL;
   if ( (type == RNDFOCUS_MLC) ||
        (type == ELEKTA_MLC)   ||
        (type == VARIAN_MLC)      )
   {
      if (r_curve != NULL)
      {
         if ( (radius_curve = new (nothrow) real) == NULL )
            xvmc_error("multi_leaf::init",
               "cannot allocate memory for radius of leaf curvature",8);
         *radius_curve = *r_curve;
      }
   }

   // number of leaf pairs
   num_pairs = num;

   // allocate memory for the leaf width array
   width = NULL;
   if ( (width = new (nothrow) real[num_pairs]) == NULL)
   {
      xvmc_error("multi_leaf::init",
                 "cannot allocate memory for leaf widths",8);
   }

   // allocate memory for the perpendicular leaf limits
   lperp = NULL;
   if ( (lperp = new (nothrow) real[num_pairs+1]) == NULL)
   {
      xvmc_error("multi_leaf::init",
                 "cannot allocate memory for perpendicular leaf limits",8);
   }

   // initialize MLC width for default leaf width of 1 cm
   max_perp  =  real(num_pairs)*ONE/TWO;
   min_perp  = -max_perp;

   // default leaf width is 1 cm at the iso-center plane
   lperp[0] = min_perp;
   for (register int i_pair=0; i_pair<num_pairs; ++i_pair)
   {
      width[i_pair]   = ONE;
      lperp[i_pair+1] = width[i_pair] + lperp[i_pair];
   }

   // allocate memory for leaf position arrays
   left = NULL;
   if ( (left = new (nothrow) real[num_pairs]) == NULL)
   {
      xvmc_error("multi_leaf::init",
                 "cannot allocate memory for left leaf positions",8);
   }
   right = NULL;
   if ( (right = new (nothrow) real[num_pairs]) == NULL)
   {
      xvmc_error("multi_leaf::init",
                 "cannot allocate memory for rigth leaf positions",8);
   }

   // allocate memory for leaf starting position arrays
   start_left = NULL;
   if ( (start_left = new (nothrow) real[num_pairs]) == NULL)
   {
      xvmc_error("multi_leaf::init",
                 "cannot allocate memory for left leaf starting positions",8);
   }
   start_right = NULL;
   if ( (start_right = new (nothrow) real[num_pairs]) == NULL)
   {
      xvmc_error("multi_leaf::init",
                 "cannot allocate memory for rigth leaf starting positions",8);
   }

   // allocate memory for leaf stopping position arrays
   stop_left = NULL;
   if ( (stop_left = new (nothrow) real[num_pairs]) == NULL)
   {
      xvmc_error("multi_leaf::init",
                 "cannot allocate memory for left leaf stopping positions",8);
   }
   stop_right = NULL;
   if ( (stop_right = new (nothrow) real[num_pairs]) == NULL)
   {
      xvmc_error("multi_leaf::init",
                 "cannot allocate memory for rigth leaf stopping positions",8);
   }

   // reset MLC, initialize leaf positions by zero beam opening
   reset();
}

// delete MLC
multi_leaf::~multi_leaf(void)
{
   if (material != NULL) {
      delete [] material; material = NULL; }
   if (upper_limit != NULL) {
      delete upper_limit; upper_limit = NULL; }
   if (lower_limit != NULL) {
      delete lower_limit; lower_limit = NULL; }
   if (z_center_curve != NULL) {
      delete z_center_curve; z_center_curve = NULL; }
   if (radius_curve != NULL) {
      delete radius_curve; radius_curve = NULL; }
   delete [] width; width = NULL;
   delete [] lperp; lperp = NULL;
   delete [] left;  left  = NULL;
   delete [] right; right = NULL;
   delete [] start_left;  start_left  = NULL;
   delete [] start_right; start_right = NULL;
   delete [] stop_left;   stop_left   = NULL;
   delete [] stop_right;  stop_right  = NULL;
}

// reset leaf positions (zero beam opening)
void multi_leaf::reset(void)
{
   // initialize leaf positions by zero openings
   for (register int i=0; i<num_pairs; ++i)
   {
      left[i]        = ZERO;
      right[i]       = ZERO;
      start_left[i]  = ZERO;
      start_right[i] = ZERO;
      stop_left[i]   = ZERO;
      stop_right[i]  = ZERO;
   }

   // initialize minimum and maximum leaf openings
   max_right       = ZERO;
   min_left        = ZERO;
   max_start_right = ZERO;
   min_start_left  = ZERO;
   max_stop_right  = ZERO;
   min_stop_left   = ZERO;
}

// change width of one leaf pair, if this is impossible -> return false
bool multi_leaf::change_width(int i_new, real w_new)
{
   // check leaf pair index
   if (i_new < 0)          return(false);
   if (i_new >= num_pairs) return(false);

   // new leaf width
   width[i_new]  = w_new;

   // adjust minimum and maximum beam openings
   // perpendicular to the leaf direction
   real total_width = ZERO;
   for (register int i_pair=0; i_pair<num_pairs; ++i_pair)
   {
      total_width += width[i_pair];
   }
   max_perp = total_width/TWO;
   min_perp = -max_perp;

   // calculate new perpendicular leaf limits
   lperp[0] = min_perp;
   for (register int k_pair=0; k_pair<num_pairs; ++k_pair)
   {
      lperp[k_pair+1] = width[k_pair] + lperp[k_pair];
   }

   // one leaf width is changed
   return(true);
}

// change position of one leaf pair, if this is impossible -> return false
bool multi_leaf::change_pair(int i_new, real l_new, real r_new)
{
   // check leaf pair index
   if (i_new < 0)          return(false);
   if (i_new >= num_pairs) return(false);

   // check leaf positions
   if (l_new > r_new)
   {
      xvmc_warning("multi_leaf::change_pair",
         "cannot change leaf positions, left position > right position",1);
      return(false);
   }

   // new leaf positions
   left[i_new]  = l_new;
   right[i_new] = r_new;

   // adjust minimum and maximum beam openings
   min_left  = min_of(min_left, l_new);
   max_right = max_of(max_right,r_new);

   // one leaf pair is changed
   return(true);
}

// change starting position of one leaf pair,
// if this is impossible -> return false
bool multi_leaf::change_start_pair(int i_new, real l_new, real r_new)
{
   // check leaf pair index
   if (i_new < 0)          return(false);
   if (i_new >= num_pairs) return(false);

   // check leaf positions
   if (l_new > r_new)
   {
      xvmc_warning("multi_leaf::change_start_pair",
         "cannot change leaf positions, left position > right position",1);
      return(false);
   }

   // new leaf positions
   start_left[i_new]  = l_new;
   start_right[i_new] = r_new;

   // adjust minimum and maximum beam openings
   min_start_left  = min_of(min_start_left, l_new);
   max_start_right = max_of(max_start_right,r_new);

   // one leaf pair is changed
   return(true);
}

// change stopping position of one leaf pair,
// if this is impossible -> return false
bool multi_leaf::change_stop_pair(int i_new, real l_new, real r_new)
{
   // check leaf pair index
   if (i_new < 0)          return(false);
   if (i_new >= num_pairs) return(false);

   // check leaf positions
   if (l_new > r_new)
   {
      xvmc_warning("multi_leaf::change_stop_pair",
         "cannot change leaf positions, left position > right position",1);
      return(false);
   }

   // new leaf positions
   stop_left[i_new]  = l_new;
   stop_right[i_new] = r_new;

   // adjust minimum and maximum beam openings
   min_stop_left  = min_of(min_stop_left, l_new);
   max_stop_right = max_of(max_stop_right,r_new);

   // one leaf pair is changed
   return(true);
}

// get width of leaf i
real multi_leaf::get_width(int i)
{
   if (i < 0)
   {
      xvmc_warning("multi_leaf::get_width",
                   "i < 0, return zero",1);
      return(ZERO);
   }

   if (i >= num_pairs)
   {
      xvmc_warning("multi_leaf::get_width",
                   "i >= num_pairs, return zero",1);
      return(ZERO);
   }

   return(width[i]);
}

// get perpendicular leaf limit i
real multi_leaf::get_lperp(int i)
{
   if (i < 0)
   {
      xvmc_warning("multi_leaf::get_lperp",
                   "i < 0, return zero",1);
      return(ZERO);
   }

   if (i > num_pairs)
   {
      xvmc_warning("multi_leaf::get_lperp",
                   "i > num_pairs, return zero",1);
      return(ZERO);
   }

   return(lperp[i]);
}

// get position of left leaf i
real multi_leaf::get_left(int i)
{
   if (i < 0)
   {
      xvmc_warning("multi_leaf::get_left",
                   "i < 0, return zero",1);
      return(ZERO);
   }

   if (i >= num_pairs)
   {
      xvmc_warning("multi_leaf::get_left",
                   "i >= num_pairs, return zero",1);
      return(ZERO);
   }

   return(left[i]);
}

// get position of right leaf i
real multi_leaf::get_right(int i)
{
   if (i < 0)
   {
      xvmc_warning("multi_leaf::get_right",
                   "i < 0, return zero",1);
      return(ZERO);
   }

   if (i >= num_pairs)
   {
      xvmc_warning("multi_leaf::get_right",
                   "i >= num_pairs, return zero",1);
      return(ZERO);
   }

   return(right[i]);
}

// get starting position of left leaf i
real multi_leaf::get_start_left(int i)
{
   if (i < 0)
   {
      xvmc_warning("multi_leaf::get_start_left",
                   "i < 0, return zero",1);
      return(ZERO);
   }

   if (i >= num_pairs)
   {
      xvmc_warning("multi_leaf::get_start_left",
                   "i >= num_pairs, return zero",1);
      return(ZERO);
   }

   return(start_left[i]);
}

// get starting position of right leaf i
real multi_leaf::get_start_right(int i)
{
   if (i < 0)
   {
      xvmc_warning("multi_leaf::get_start_right",
                   "i < 0, return zero",1);
      return(ZERO);
   }

   if (i >= num_pairs)
   {
      xvmc_warning("multi_leaf::get_start_right",
                   "i >= num_pairs, return zero",1);
      return(ZERO);
   }

   return(start_right[i]);
}

// get stopping position of left leaf i
real multi_leaf::get_stop_left(int i)
{
   if (i < 0)
   {
      xvmc_warning("multi_leaf::get_stop_left",
                   "i < 0, return zero",1);
      return(ZERO);
   }

   if (i >= num_pairs)
   {
      xvmc_warning("multi_leaf::get_stop_left",
                   "i >= num_pairs, return zero",1);
      return(ZERO);
   }

   return(stop_left[i]);
}

// get stopping position of right leaf i
real multi_leaf::get_stop_right(int i)
{
   if (i < 0)
   {
      xvmc_warning("multi_leaf::get_stop_right",
                   "i < 0, return zero",1);
      return(ZERO);
   }

   if (i >= num_pairs)
   {
      xvmc_warning("multi_leaf::get_stop_right",
                   "i >= num_pairs, return zero",1);
      return(ZERO);
   }

   return(stop_right[i]);
}
