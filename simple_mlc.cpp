/*****************************************************************************
 * simple_mlc.cpp:                                                           *
 *    class member functions for:                                            *
 *       simple_mlc:         simple multi-leaf collimator model              *
 *                                                                           *
 * Copyright (C) 2002    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 31.01.2002      *
 *    simple dynamic MLC                                  MF 09.10.2003      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <new>
using namespace std;

#include "simple_mlc.h"

// ****************************************
// declare functions and global variables
// ****************************************

// ********************************************
// member functions of class simple_mlc
// ********************************************

// create MLC based on the leaf positions at the iso-center plane,
// the source to iso-center distance, the upper and lower MLC limits
simple_mlc::simple_mlc(multi_leaf *iso_mlc, const real ini_iso_distance)
 : beam_modifier()
{
   // set beam modifier type
   type = MODIFIER_SIMPLE_MLC;

   // set beam modifier mode
   switch (iso_mlc->get_mode())
   {
   case STATIC_MLC:
      mode = MODIFIER_STATIC;
      break;
   case DYNAMIC_MLC:
      mode = MODIFIER_DYNAMIC;
      break;
   default:
      xvmc_error("simple_mlc::simple_mlc","unknown MLC mode",8);
      break;
   }

   // check MLC type
   if (iso_mlc->get_type() != SIMPLE_MLC)
   {
      xvmc_error("simple_mlc::simple_mlc",
                 "this is not a simple MLC",8);
   }

   // set source to iso-center distance
   iso_distance = ini_iso_distance;

   // change upper and lower beam modifier distances
   if (iso_mlc->get_upper() != NULL) upper_distance = *iso_mlc->get_upper();
   if (iso_mlc->get_lower() != NULL) lower_distance = *iso_mlc->get_lower();

   // scaling factor (iso-center plane to lower modifier plane)
   real scale = lower_distance/iso_distance;

   // get leaf moving direction (X or Y)
   char mlc_xy = iso_mlc->get_xy();

   // get number of leaf pairs
   int num_pairs = iso_mlc->get_num();

   // create MLC
   lower_mlc = NULL;
   if ( (lower_mlc = new (nothrow) multi_leaf(SIMPLE_MLC,iso_mlc->get_mode(),
                         mlc_xy,num_pairs,
                         NULL,NULL,NULL,NULL,NULL)) == NULL)
   {
      xvmc_error("simple_mlc::simple_mlc",
                 "cannot create MLC",8);
   }

   // scale leaf widths and positions (iso-center plane to lower modifier plane)
   for (register int i=0; i<num_pairs; ++i)
   {
      real width = scale*iso_mlc->get_width(i);
      if (!lower_mlc->change_width(i,width))
      {
         xvmc_error("simple_mlc::simple_mlc",
                    "cannot change leaf width",8);
      }

      real left  = scale*iso_mlc->get_left(i);
      real right = scale*iso_mlc->get_right(i);
      if (!lower_mlc->change_pair(i,left,right))
      {
         xvmc_error("simple_mlc::simple_mlc",
                    "cannot change leaf position",8);
      }

      real start_left  = scale*iso_mlc->get_start_left(i);
      real start_right = scale*iso_mlc->get_start_right(i);
      if (!lower_mlc->change_start_pair(i,start_left,start_right))
      {
         xvmc_error("simple_mlc::simple_mlc",
                    "cannot change leaf starting position",8);
      }

      real stop_left  = scale*iso_mlc->get_stop_left(i);
      real stop_right = scale*iso_mlc->get_stop_right(i);
      if (!lower_mlc->change_stop_pair(i,stop_left,stop_right))
      {
         xvmc_error("simple_mlc::simple_mlc",
                    "cannot change leaf stopping position",8);
      }
   }
}

// set beam modifier and MLC modes using a new beam modifier mode
void simple_mlc::set_mode(modifier_mode new_mode)
{
   // set beam modifier mode
   mode = new_mode;

   // set MLC mode
   switch (mode)
   {
   case MODIFIER_STATIC:
      lower_mlc->set_mode(STATIC_MLC);
      break;
   case MODIFIER_DYNAMIC:
      lower_mlc->set_mode(DYNAMIC_MLC);
      break;
   default:
      xvmc_error("simple_mlc::set_mode","unknown beam modifier mode",8);
      break;
   }

   return;
}

// set beam modifier and MLC modes using a new MLC mode
void simple_mlc::set_mode(mlc_mode new_mode)
{
   // set beam modifier mode
   switch (new_mode)
   {
   case STATIC_MLC:
      mode = MODIFIER_STATIC;
      break;
   case DYNAMIC_MLC:
      mode = MODIFIER_DYNAMIC;
      break;
   default:
      xvmc_error("simple_mlc::set_mode","unknown MLC mode",8);
      break;
   }

   // set MLC mode
   lower_mlc->set_mode(new_mode);

   return;
}

// change position of one leaf pair,
// the new leaf positions are defined at the iso-center plane
void simple_mlc::change_pair(const int  &pair_index,
                             const real &new_left,
                             const real &new_right)
{
   // scaling factor (iso-center plane to lower modifier plane)
   real scale = lower_distance/iso_distance;

   // scale new leaf positions
   real left  = scale*new_left;
   real right = scale*new_right;
   if (!lower_mlc->change_pair(pair_index,left,right))
   {
      xvmc_error("simple_mlc::change_pair",
                 "cannot change leaf positions",8);
   }

   return;
}

// change starting position of one leaf pair,
// the new leaf positions are defined at the iso-center plane
void simple_mlc::change_start_pair(const int  &pair_index,
                                   const real &new_start_left,
                                   const real &new_start_right)
{
   // scaling factor (iso-center plane to lower modifier plane)
   real scale = lower_distance/iso_distance;

   // scale new leaf positions
   real start_left  = scale*new_start_left;
   real start_right = scale*new_start_right;
   if (!lower_mlc->change_start_pair(pair_index,start_left,start_right))
   {
      xvmc_error("simple_mlc::change_start_pair",
                 "cannot change leaf starting positions",8);
   }

   return;
}

// change stopping position of one leaf pair,
// the new leaf positions are defined at the iso-center plane
void simple_mlc::change_stop_pair(const int  &pair_index,
                             const real &new_stop_left,
                             const real &new_stop_right)
{
   // scaling factor (iso-center plane to lower modifier plane)
   real scale = lower_distance/iso_distance;

   // scale new leaf positions
   real stop_left  = scale*new_stop_left;
   real stop_right = scale*new_stop_right;
   if (!lower_mlc->change_stop_pair(pair_index,stop_left,stop_right))
   {
      xvmc_error("simple_mlc::change_stop_pair",
                 "cannot change leaf stopping positions",8);
   }

   return;
}

// change leaf positions
void simple_mlc::change_leaf_positions(multi_leaf *iso_mlc)
{
   // check MLC type
   if (iso_mlc->get_type() != SIMPLE_MLC)
   {
      xvmc_error("simple_mlc::change_leaf_positions",
                 "inconsistent MLC types",8);
   }

   // check MLC mode
   if (lower_mlc->get_mode() != iso_mlc->get_mode())
   {
      xvmc_error("simple_mlc::change_leaf_positions",
                 "inconsistent MLC modes",8);
   }

   // check upper and lower beam modifier distances
   if (iso_mlc->get_upper() != NULL)
   {
      if (upper_distance != *iso_mlc->get_upper())
      {
         xvmc_error("simple_mlc::change_leaf_positions",
                    "upper MLC distances are inconsistent",8);
      }
   }
   if (iso_mlc->get_lower() != NULL)
   {
      if (lower_distance != *iso_mlc->get_lower())
      {
         xvmc_error("simple_mlc::change_leaf_positions",
                    "lower MLC distances are inconsistent",8);
      }
   }

   // check leaf moving direction (X or Y)
   if (lower_mlc->get_xy() != iso_mlc->get_xy())
   {
      xvmc_error("simple_mlc::change_leaf_positions",
                 "leaf moving directions are inconsistent",8);
   }

   // check number of leaf pairs
   int num_pairs = lower_mlc->get_num();
   if (num_pairs != iso_mlc->get_num())
   {
      xvmc_error("simple_mlc::change_leaf_positions",
                 "number of leaf pairs are inconsistent",8);
   }

   // change leaf positions
   for (register int i=0; i<num_pairs; ++i)
   {
      change_pair(i,iso_mlc->get_left(i),
                    iso_mlc->get_right(i));
      change_start_pair(i,iso_mlc->get_start_left(i),
                          iso_mlc->get_start_right(i));
      change_stop_pair(i,iso_mlc->get_stop_left(i),
                         iso_mlc->get_stop_right(i));
   }

   return;
}

// adjust simple dynamic MLC, in static mode do nothing
void simple_mlc::adjust(long i_history, long n_history)
{
   // in static mode do nothing
   if (mode == MODIFIER_STATIC) return;

   // number of steps to move the leafes
   long n_steps = 1000;

   // number of histories per step
   long n_history_per_step = n_history/n_steps;
   if (n_history_per_step < 1) n_history_per_step = 1;

   if ( (i_history%n_history_per_step) != 0 ) return;

   // calculate scaling factor to adjust the leafes
   real scale = real(i_history)/real(n_history);

   // change leaf positions
   int num_pairs = lower_mlc->get_num();
   for (register int i=0; i<num_pairs; ++i)
   {
      // get starting positions of left and right leaf
      real start_left  = lower_mlc->get_start_left(i);
      real start_right = lower_mlc->get_start_right(i);

      // get moving distance of left and right leaf
      real dist_left  = lower_mlc->get_stop_left(i)  - start_left;
      real dist_right = lower_mlc->get_stop_right(i) - start_right;

      // new left and right leaf positions
      real left  = start_left  + scale*dist_left;
      real right = start_right + scale*dist_right;
      lower_mlc->change_pair(i,left,right);
   }

   return;
}

// transport particle through multi-leaf collimator,
// if the particle is deleted return false, otherwise return true
bool simple_mlc::transport(particle_parameters &p, ranmar &rndm)
{
   // length of line from actual particle position to lower plane
   real length = (lower_distance-p.pos.z)/p.dir.z;

   // intersection point of particle path with lower plane
   p.pos.x += length*p.dir.x;
   p.pos.y += length*p.dir.y;
   p.pos.z  = lower_distance;

   // test whether the particle is inside or outside the MLC opening
   if (lower_mlc->inside(p.pos.x,p.pos.y)) return(true);

   // the particle is outside
   return(false);
}

// transport particle position and direction through the MLC,
// if the particle is deleted return false, otherwise return true
bool simple_mlc::transport(real_3 &pos, real_3 &dir, ranmar &rndm)
{
   // length of line from actual particle position to lower plane
   real length = (lower_distance-pos.z)/dir.z;

   // intersection point of particle path with lower plane
   pos.x += length*dir.x;
   pos.y += length*dir.y;
   pos.z  = lower_distance;

   // test whether the particle is inside or outside the MLC opening
   if (lower_mlc->inside(pos.x,pos.y)) return(true);

   // the particle is outside
   return(false);
}
