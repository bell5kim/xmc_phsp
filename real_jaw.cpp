/*****************************************************************************
 * real_jaw.cpp:                                                             *
 *    class member functions for:                                            *
 *       real_jaw:           real 3D Monte Carlo JAW model                   *
 *                                                                           *
 * Copyright (C) 2002    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 31.01.2002      *
 *    real dynamic MLC                                    MF 23.10.2003      *
 *    fixed bug in transport (bug if distance <= 0)       MF 22.01.2004      *
 *    Adopted from real_mlc.cpp                           JK Feb 22, 2008    *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <new>
using namespace std;

#include "real_jaw.h"

// ****************************************
// declare functions and global variables
// ****************************************

// ********************************************
// member functions of class real_jaw
// ********************************************

// create JAW based on the jaw positions at the iso-center plane
real_jaw::real_jaw(MC_jaw *iso_jaw, const real iso_distance)
 : beam_modifier()
{
   // set beam modifier type
   type = MODIFIER_REAL_JAW;

   // set beam modifier mode
   switch (iso_jaw->get_mode())
   {
   case STATIC_JAW:
      mode = MODIFIER_STATIC;
      break;
   case DYNAMIC_JAW:
      mode = MODIFIER_DYNAMIC;
      break;
   default:
      xvmc_error("real_jaw::real_jaw","unknown JAW mode",8);
      break;
   }

   // change upper and lower beam modifier distance
   if (iso_jaw->get_upper() != ZERO) upper_distance = iso_jaw->get_upper();
   if (iso_jaw->get_lower() != ZERO) lower_distance = iso_jaw->get_lower();

   // change jaw opening
   real open_left  = iso_jaw->get_open_left();
   real open_right = iso_jaw->get_open_right();

   // check type and create JAW
   jaw = NULL;
   if ( (jaw = new (nothrow) MC_jaw(X, REAL_JAW, STATIC_JAW,
                                    open_left, open_right,   // open_left and open_right
                                    upper_distance, lower_distance,   // z_min and z_max
                                    "tungsten", "air")) == NULL)
   {
      xvmc_error("real_jaw::real_jaw",
                 "cannot create double focussing MLC",8);
   }

}

// set beam modifier and JAW modes using a new beam modifier mode
void real_jaw::set_mode(modifier_mode new_mode)
{
   // set beam modifier mode
   mode = new_mode;

   // set JAW mode
   switch (mode)
   {
   case MODIFIER_STATIC:
      jaw->set_mode(STATIC_JAW);
      break;
   case MODIFIER_DYNAMIC:
      jaw->set_mode(DYNAMIC_JAW);
      break;
   default:
      xvmc_error("real_jaw::set_mode","unknown beam modifier mode",8);
      break;
   }

   return;
}

// set beam modifier and JAW modes using a new JAW mode
void real_jaw::set_mode(jaw_mode new_mode)
{
   // set beam modifier mode
   switch (new_mode)
   {
   case STATIC_JAW:
      mode = MODIFIER_STATIC;
      break;
   case DYNAMIC_JAW:
      mode = MODIFIER_DYNAMIC;
      break;
   default:
      xvmc_error("real_jaw::set_mode","unknown JAW mode",8);
      break;
   }

   // set JAW mode
   jaw->set_mode(new_mode);

   return;
}

// adjust real dynamic JAW, in static mode do nothing
void real_jaw::adjust(long i_history, long n_history)
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

   // change jaw positions
 
   // -------------------------------------------------------- 
   //   NOT Implemented YET!
   //   Must be changed to statistical sampling for dynamic
   // -------------------------------------------------------- 

   return;
}

// transport particle through jaws,
// if the particle is deleted return false, otherwise return true
bool real_jaw::transport(particle_parameters &p, ranmar &rndm)
{
   // get the starting plane
   MC_plane *starting_plane = jaw->get_starting_plane();

   // get plane normal vector
   real_3 normal = starting_plane->get_normal();

   // get distance of plane to origin
   real dist2zero = starting_plane->get_dist2zero();

   // scalar products
   real product_pos = p.pos.x*normal.x + p.pos.y*normal.y + p.pos.z*normal.z;
   real product_dir = p.dir.x*normal.x + p.dir.y*normal.y + p.dir.z*normal.z;

   // if particle direction is parallel to the plane (perpendicular to the
   // normal vector) --> delete particle, i.e. return false
   if (product_dir == ZERO) return(false);

   // get distance of the photon to the starting plane
   real distance = (dist2zero-product_pos)/product_dir;

   // move photon to this plane
   p.pos.x += distance*p.dir.x;
   p.pos.y += distance*p.dir.y;
   p.pos.z += distance*p.dir.z;

   // transport charged particle ...
   if (p.type != PHOTON) return(jaw->primary_electron(p,starting_plane,rndm));

   // ... or photon through the JAW
   return(jaw->primary_photon(p,starting_plane,rndm));
}
