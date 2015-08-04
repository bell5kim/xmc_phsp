/*****************************************************************************
 * real_mlc.cpp:                                                             *
 *    class member functions for:                                            *
 *       real_mlc:           real 3D Monte Carlo MLC model                   *
 *                                                                           *
 * Copyright (C) 2002    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 31.01.2002      *
 *    real dynamic MLC                                    MF 23.10.2003      *
 *    fixed bug in transport (bug if distance <= 0)       MF 22.01.2004      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <new>
using namespace std;

#include "real_mlc.h"

// ****************************************
// declare functions and global variables
// ****************************************

// ********************************************
// member functions of class real_mlc
// ********************************************

// create MLC based on the leaf positions at the iso-center plane
real_mlc::real_mlc(multi_leaf *iso_mlc, const real iso_distance)
 : beam_modifier()
{
   // set beam modifier type
   type = MODIFIER_REAL_MLC;

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
      xvmc_error("real_mlc::real_mlc","unknown MLC mode",8);
      break;
   }

   // change upper and lower beam modifier distance
   if (iso_mlc->get_upper() != NULL) upper_distance = *iso_mlc->get_upper();
   if (iso_mlc->get_lower() != NULL) lower_distance = *iso_mlc->get_lower();

   // check type and create MLC
   mlc = NULL;
   switch (iso_mlc->get_type())
   {
   case DBLFOCUS_MLC:
      if ( (mlc = new (nothrow) MC_mlc_2focus(iso_mlc,"air",iso_distance))
                == NULL)
      {
         xvmc_error("real_mlc::real_mlc",
                    "cannot create double focussing MLC",8);
      }
      break;
   case RNDFOCUS_MLC:
      if ( (mlc = new (nothrow) MC_mlc_rfocus(iso_mlc,"air",iso_distance))
                == NULL)
      {
         xvmc_error("real_mlc::real_mlc",
                    "cannot create MLC with curved leaf ends",8);
      }
      break;
   case ELEKTA_MLC:
      if ( (mlc = new (nothrow) MC_mlc_elekta(iso_mlc,"air",iso_distance))
                == NULL)
      {
         xvmc_error("real_mlc::real_mlc",
                    "cannot create ELEKTA MLC",8);
      }
      break;
   case VARIAN_MLC:
      if ( (mlc = new (nothrow) MC_mlc_varian(iso_mlc,"air",iso_distance))
                == NULL)
      {
         xvmc_error("real_mlc::real_mlc",
                    "cannot create VARIAN MLC",8);
      }
      break;
   default:
      xvmc_error("real_mlc::real_mlc",
                 "this is no real 3D Monte Carlo MLC model",8);
      break;
   }
}

// set beam modifier and MLC modes using a new beam modifier mode
void real_mlc::set_mode(modifier_mode new_mode)
{
   // set beam modifier mode
   mode = new_mode;

   // set MLC mode
   switch (mode)
   {
   case MODIFIER_STATIC:
      mlc->set_mode(STATIC_MLC);
      break;
   case MODIFIER_DYNAMIC:
      mlc->set_mode(DYNAMIC_MLC);
      break;
   default:
      xvmc_error("real_mlc::set_mode","unknown beam modifier mode",8);
      break;
   }

   return;
}

// set beam modifier and MLC modes using a new MLC mode
void real_mlc::set_mode(mlc_mode new_mode)
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
      xvmc_error("real_mlc::set_mode","unknown MLC mode",8);
      break;
   }

   // set MLC mode
   mlc->set_mode(new_mode);

   return;
}

// change leaf positions
void real_mlc::change_leaf_positions(multi_leaf *iso_mlc)
{
   // check MLC type
   if (mlc->get_type() != iso_mlc->get_type())
   {
      xvmc_error("real_mlc::change_leaf_positions",
                 "MLC types are inconsistent",8);
   }

   // check MLC mode
   if (mlc->get_mode() != iso_mlc->get_mode())
   {
      xvmc_error("real_mlc::change_leaf_positions",
                 "MLC modes are inconsistent",8);
   }

   // check upper and lower beam modifier distances
   if (iso_mlc->get_upper() != NULL)
   {
      if (upper_distance != *iso_mlc->get_upper())
      {
         xvmc_error("real_mlc::change_leaf_positions",
                    "upper MLC distances are inconsistent",8);
      }
   }
   if (iso_mlc->get_lower() != NULL)
   {
      if (lower_distance != *iso_mlc->get_lower())
      {
         xvmc_error("real_mlc::change_leaf_positions",
                    "lower MLC distances are inconsistent",8);
      }
   }

   // check leaf moving direction (X or Y)
   if (mlc->get_xy() != iso_mlc->get_xy())
   {
      xvmc_error("real_mlc::change_leaf_positions",
                 "leaf moving directions are inconsistent",8);
   }

   // check number of leaf pairs
   int num_pairs = mlc->get_num();
   if (num_pairs != iso_mlc->get_num())
   {
      xvmc_error("real_mlc::change_leaf_positions",
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

// adjust real dynamic MLC, in static mode do nothing
void real_mlc::adjust(long i_history, long n_history)
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
   int num_pairs = mlc->get_num();
   for (register int i=0; i<num_pairs; ++i)
   {
      // get starting positions of left and right leaf (iso-center plane)
      real start_left  = mlc->get_start_left(i);
      real start_right = mlc->get_start_right(i);

      // get moving distance of left and right leaf (iso-center plane)
      real dist_left  = mlc->get_stop_left(i)  - start_left;
      real dist_right = mlc->get_stop_right(i) - start_right;

      // new left and right leaf positions (iso-center plane)
      real left  = start_left  + scale*dist_left;
      real right = start_right + scale*dist_right;
      change_pair(i,left,right);
   }

   return;
}

// transport particle through multi-leaf collimator,
// if the particle is deleted return false, otherwise return true
bool real_mlc::transport(particle_parameters &p, ranmar &rndm)
{
   // get the starting plane
   MC_plane *starting_plane = mlc->get_starting_plane();

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
   if (p.type != PHOTON) return(mlc->primary_electron(p,starting_plane,rndm));

   // ... or photon through the MLC
   return(mlc->primary_photon(p,starting_plane,rndm));
}
