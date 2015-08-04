/*****************************************************************************
 * jaw.cpp:                                                                  *
 *    class member functions for:                                            *
 *       jaws: nominal jaws (iso-center plane)                               *
 *                                                                           *
 * Copyright (C) 2001    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 15.06.2001      *
 *    dynamic MLC mode implemented                        MF 09.10.2003      *
 *    modified from multi_leaf.cpp                        JK 01.30.2008      *
 *                                                                           *
 *****************************************************************************/

#include <string.h>
#include <new>
using namespace std;

#include "definitions.h"
#include "global.h"
#include "jaws.h"

// ****************************************
// member functions of class jaws
// ****************************************

// initialize JAW by type, mode, jaw moving direction (X or Y),
// number of jaw pairs (must be 1), JAW material,
// upper and lower limits of the JAW
void jaws::init(const jaw_type  inp_type, const jaw_mode  inp_mode,
                const char      xy,
                const int       num,      const char     *inp_material,
                const real      z_upper,  const real     z_lower)
{
   // set JAW type
   type = inp_type;

   // set JAW mode
   mode = inp_mode;

   // check jaw moving direction (X or Y)
   switch (xy)
   {
   case 'x':    // X jaw
   case 'X':    // X jaw
      x_jaw = true;
      break;
   case 'y':    // Y jaw
   case 'Y':    // Y jaw
      x_jaw = false;
      break;
   default:
      xvmc_error("jaws::init","there are only X or Y JAWs",8);
      break;
   }

   // set JAW material
   material = NULL;
   if (inp_material != NULL)
   {
      if ( (material = new (nothrow) char[strlen(inp_material)+1]) == NULL )
         xvmc_error("jaws::init",
                    "cannot allocate memory for MLC material",8);
      strcpy(material,inp_material);
   }

   // set upper and lower JAW limits
   upper_limit = z_upper;;
   lower_limit = z_lower;

   // number of jaw pairs (it must be ONE)
   num_pairs = 1;
   if (num_pairs != ONE) {
      xvmc_error("jaws::init",
                 "number of jaw pairs must be ONE",8);
   }

   // initialize JAW width for default jaw width of 30 cm
   width = 30.0;

   // initialize Jaw positions
   left = ZERO;
   right = ZERO;

   // jaw starting positions
   start_left = ZERO;
   start_right = ZERO;

   // jaw stopping position arrays
   stop_left = -20.0;
   stop_right = 20.0;

   // reset MLC, initialize jaw positions by zero beam opening
   reset();
}

// delete MLC
jaws::~jaws(void)
{
   material = NULL;
}

// reset jaw positions (zero beam opening)
void jaws::reset(void)
{
   // initialize jaw positions by zero openings
   left        = ZERO;
   right       = ZERO;
   start_left  = ZERO;
   start_right = ZERO;
   stop_left   = ZERO;
   stop_right  = ZERO;

   // initialize minimum and maximum jaw openings
   max_width       = ZERO;
   min_width       = ZERO;
   max_start_right = ZERO;
   min_start_left  = ZERO;
   max_stop_right  = ZERO;
   min_stop_left   = ZERO;
}

// change position of jaw pair, if this is impossible -> return false
bool jaws::change_jaws(real l_new, real r_new)
{
   // check jaw positions
   if (l_new > r_new)
   {
      xvmc_warning("jaws::change_jaws",
         "cannot change jaw positions, left position > right position",1);
      return(false);
   }

   // new jaw positions
   left  = l_new;
   right = r_new;

   // adjust minimum and maximum beam openings
   min_width = min_of(min_width,l_new);
   max_width = max_of(max_width,r_new);

   // jaw pair is changed
   return(true);
}

// change starting position of jaw pair,
// if this is impossible -> return false
bool jaws::change_start_jaws(real l_new, real r_new)
{
   // check jaw positions
   if (l_new > r_new)
   {
      xvmc_warning("jaws::change_start_jaws",
         "cannot change jaw positions, left position > right position",1);
      return(false);
   }

   // new jaw positions
   start_left  = l_new;
   start_right = r_new;

   // adjust minimum and maximum beam openings
   min_start_left  = min_of(min_start_left, l_new);
   max_start_right = max_of(max_start_right,r_new);

   // jaw pair is changed
   return(true);
}

// change stopping position of one jaw pair,
// if this is impossible -> return false
bool jaws::change_stop_jaws(real l_new, real r_new)
{
   // check jaw positions
   if (l_new > r_new)
   {
      xvmc_warning("jaws::change_stop_pair",
         "cannot change jaw positions, left position > right position",1);
      return(false);
   }

   // new jaw positions
   stop_left  = l_new;
   stop_right = r_new;

   // adjust minimum and maximum beam openings
   min_stop_left  = min_of(min_stop_left, l_new);
   max_stop_right = max_of(max_stop_right,r_new);

   // jaw pair is changed
   return(true);
}
