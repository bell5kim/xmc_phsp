/*****************************************************************************
 * MC_jaws_focus.cpp:                                                        *
 *    class member functions for:                                            *
 *       MC_jaws_focus:  a pair of focussing jaws                            *
 *                                                                           *
 * Copyright (C) 2001    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 31.08.2001      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <new>
using namespace std;

#include "MC_jaws_focus.h"

// ****************************************
// member functions of class MC_jaws_focus
// ****************************************

// initialize jaws pair
void MC_jaws_focus::init(const axis &ini_type,
                         const real &ini_left,      const real &ini_right,
                         const real &ini_x_min,     const real &ini_x_max,
                         const real &ini_y_min,     const real &ini_y_max,
                         const real &ini_z_min,     const real &ini_z_max,
                         const real &ini_iso_dist,  const real &ini_z_focus,
                         const char *bars_material, const char *open_material)
{
   // check and set jaws type (X or Y)
   if (ini_type == Z)
   {
      xvmc_error("MC_jaws_focus::init",
                 "there are only X or Y jaws",8);
   }
   type = ini_type;

   // determine (ESTAR) material electron data file paths
   bars_material_estar = NULL;
   if ( (bars_material_estar = get_file_path(bars_material,"estar")) == NULL )
   {
      xvmc_error("MC_jaws_focus::init",
                 "cannot determine bar material electron data file path",8);
   }

   open_material_estar = NULL;
   if ( (open_material_estar = get_file_path(open_material,"estar")) == NULL )
   {
      xvmc_error("MC_jaws_focus::init",
         "cannot determine opening material electron data file path",8);
   }

   // determine (NIST) material cross section file paths
   bars_material_nist = NULL;
   if ( (bars_material_nist = get_file_path(bars_material,"nist")) == NULL )
   {
      xvmc_error("MC_jaws_focus::init",
                 "cannot determine bar material cross section file path",8);
   }

   open_material_nist = NULL;
   if ( (open_material_nist = get_file_path(open_material,"nist")) == NULL )
   {
      xvmc_error("MC_jaws_focus::init",
         "cannot determine opening material cross section file path",8);
   }

   // determine differential Compton cross section file path
   compton_file = NULL;
   if ( (compton_file = get_file_path("compton","data")) == NULL )
   {
      xvmc_error("MC_jaws_focus::init",
         "cannot determine differential Compton cross section file path",8);
   }

   // determine differential pair cross section file path
   pair_file = NULL;
   if ( (pair_file = get_file_path("pair","data")) == NULL )
   {
      xvmc_error("MC_jaws_focus::init",
         "cannot determine differential pair cross section file path",8);
   }

   // default origin to iso-center distance and focus position
   iso_distance = ini_iso_dist;
   z_focus      = ini_z_focus;

   // check and set outer object dimensions in x direction
   if (ini_x_max > ini_x_min)
   {
      x_min = ini_x_min;
      x_max = ini_x_max;
   }
   else
   {
      xvmc_error("MC_jaws_focus::init",
                 "the maximum x boundary must be larger than the minimum",8);
   }

   // check and set outer object dimensions in y direction
   if (ini_y_max > ini_y_min)
   {
      y_min = ini_y_min;
      y_max = ini_y_max;
   }
   else
   {
      xvmc_error("MC_jaws_focus::init",
                 "the maximum y boundary must be larger than the minimum",8);
   }

   // check and set outer object dimensions in z direction
   if (ini_z_max > ini_z_min)
   {
      z_min = ini_z_min;
      z_max = ini_z_max;
   }
   else
   {
      xvmc_error("MC_jaws_focus::init",
                 "the maximum z boundary must be larger than the minimum",8);
   }

   // check and set jaws opening
   if (type == X)
   {
      if (ini_left*z_max/iso_distance <= x_min)
      {
         xvmc_error("MC_jaws_focus::init",
                    "the left jaws opening is too large",8);
      }
      if (ini_right*z_max/iso_distance >= x_max)
      {
         xvmc_error("MC_jaws_focus::init",
                    "the right jaws position is too large",8);
      }
   }
   else
   {
      if (ini_left*z_max/iso_distance <= y_min)
      {
         xvmc_error("MC_jaws_focus::init",
                    "the left jaws opening is too large",8);
      }
      if (ini_right*z_max/iso_distance >= y_max)
      {
         xvmc_error("MC_jaws_focus::init",
                    "the right jaws position is too large",8);
      }
   }

   if (ini_right > ini_left)
   {
      open_left  = ini_left;
      open_right = ini_right;
   }
   else
   {
      xvmc_error("MC_jaws_focus::init",
                 "the right jaw position must be larger than the left",8);
   }

   // create the 6 outer planes
   plane_x_min = NULL;
   if ( (plane_x_min = new (nothrow) MC_plane(X,x_min,-1)) == NULL )
   {
      xvmc_error("MC_jaws_focus::init",
                 "cannot create minimum x plane",8);
   }
   plane_x_max = NULL;
   if ( (plane_x_max = new (nothrow) MC_plane(X,x_max,-1)) == NULL )
   {
      xvmc_error("MC_jaws_focus::init",
                 "cannot create maximum x plane",8);
   }
   plane_y_min = NULL;
   if ( (plane_y_min = new (nothrow) MC_plane(Y,y_min,-1)) == NULL )
   {
      xvmc_error("MC_jaws_focus::init",
                 "cannot create minimum y plane",8);
   }
   plane_y_max = NULL;
   if ( (plane_y_max = new (nothrow) MC_plane(Y,y_max,-1)) == NULL )
   {
      xvmc_error("MC_jaws_focus::init",
                 "cannot create maximum y plane",8);
   }

   // typically the starting plane of the particle transport
   plane_z_min = NULL;
   if ( (plane_z_min = new (nothrow) MC_plane(Z,z_min,-1)) == NULL )
   {
      xvmc_error("MC_jaws_focus::init",
                 "cannot create minimum z plane",8);
   }
   starting_plane = plane_z_min;

   // typically the final plane of the particle transport
   plane_z_max = NULL;
   if ( (plane_z_max = new (nothrow) MC_plane(Z,z_max,-1)) == NULL )
   {
      xvmc_error("MC_jaws_focus::init",
                 "cannot create maximum z plane",8);
   }
   final_plane = plane_z_max;

   // the focus point
   real_3 p_focus;
   p_focus.x = ZERO; p_focus.y = ZERO; p_focus.z = z_focus;

   // we need 2 further points for the left jaws plane
   real_3 p1_left,p2_left;
   if (type == X)
   {
      p1_left.x = open_left; p1_left.y = -10.0; p1_left.z = iso_distance;
      p2_left.x = open_left; p2_left.y =  10.0; p2_left.z = iso_distance;
   }
   else
   {
      p1_left.x = -10.0; p1_left.y = open_left; p1_left.z = iso_distance;
      p2_left.x =  10.0; p2_left.y = open_left; p2_left.z = iso_distance;
   }

   // create left jaws plane
   plane_left  = NULL;
   if ( (plane_left = new (nothrow)
         MC_plane(p_focus,p1_left,p2_left,-1)) == NULL )
   {
      xvmc_error("MC_jaws_focus::init",
                 "cannot create left jaws plane",8);
   }

   // we need 2 further points for the right jaws plane
   real_3 p1_right,p2_right;
   if (type == X)
   {
      p1_right.x = open_right; p1_right.y = -10.0; p1_right.z = iso_distance;
      p2_right.x = open_right; p2_right.y =  10.0; p2_right.z = iso_distance;
   }
   else
   {
      p1_right.x = -10.0; p1_right.y = open_right; p1_right.z = iso_distance;
      p2_right.x =  10.0; p2_right.y = open_right; p2_right.z = iso_distance;
   }

   // create right jaws plane
   plane_right = NULL;
   if ( (plane_right = new (nothrow)
         MC_plane(p_focus,p1_right,p2_right,-1)) == NULL )
   {
      xvmc_error("MC_jaws_focus::init",
                 "cannot create right jaws plane",8);
   }

   // assign object plane pointers
   separator[0] = plane_x_min;
   separator[1] = plane_x_max;
   separator[2] = plane_y_min;
   separator[3] = plane_y_max;
   separator[4] = plane_z_min;
   separator[5] = plane_z_max;
   separator[6] = plane_left;
   separator[7] = plane_right;

   // define the reference points to identify the regions
   real r_left  = open_left*(z_min+z_max)/TWO/iso_distance;
   real r_right = open_right*(z_min+z_max)/TWO/iso_distance;
   real_3 p_ref;
   p_ref.x = ZERO;
   p_ref.y = ZERO;
   p_ref.z = ZERO;

   // create regions
   if (type == X)   // X jaws
   {
      // create left jaws bar, a region bounded by 6 planes
      p_ref.x = (x_min+r_left)/TWO;
      p_ref.y = (y_min+y_max)/TWO;
      p_ref.z = (z_min+z_max)/TWO;
      if ( (piece[0] = new (nothrow) MC_volume_6p(
                                        plane_x_min, plane_left,
                                        plane_y_min, plane_y_max,
                                        plane_z_min, plane_z_max,
                                        p_ref,
                                        bars_material_estar,
                                        bars_material_nist,
                                        compton_file, pair_file)  ) == NULL )
      {
         xvmc_error("MC_jaws_focus::init",
                    "cannot create left jaws bar",8);
      }

      // create right jaws bar, a region bounded by 6 planes
      p_ref.x = (r_right+x_max)/TWO;
      p_ref.y = (y_min+y_max)/TWO;
      p_ref.z = (z_min+z_max)/TWO;
      if ( (piece[1] = new (nothrow) MC_volume_6p(
                                        plane_right, plane_x_max,
                                        plane_y_min, plane_y_max,
                                        plane_z_min, plane_z_max,
                                        p_ref,
                                        bars_material_estar,
                                        bars_material_nist,
                                        compton_file, pair_file)  ) == NULL )
      {
         xvmc_error("MC_jaws_focus::init",
                    "cannot create right jaws bar",8);
      }

      // create intermediate region between the two jaw bars,
      // a region bounded by 6 planes
      p_ref.x = (r_left+r_right)/TWO;
      p_ref.y = (y_min+y_max)/TWO;
      p_ref.z = (z_min+z_max)/TWO;
      if ( (piece[2] = new (nothrow) MC_volume_6p(
                                        plane_left,  plane_right,
                                        plane_y_min, plane_y_max,
                                        plane_z_min, plane_z_max,
                                        p_ref,
                                        open_material_estar,
                                        open_material_nist,
                                        compton_file, pair_file)  ) == NULL )
      {
         xvmc_error("MC_jaws_focus::init",
                    "cannot create intermediate region",8);
      }
   }
   else             // Y jaws
   {
      // create left jaws bar, a region bounded by 6 planes
      p_ref.x = (x_min+x_max)/TWO;
      p_ref.y = (y_min+r_left)/TWO;
      p_ref.z = (z_min+z_max)/TWO;
      if ( (piece[0] = new (nothrow) MC_volume_6p(
                                        plane_y_min, plane_left,
                                        plane_x_min, plane_x_max,
                                        plane_z_min, plane_z_max,
                                        p_ref,
                                        bars_material_estar,
                                        bars_material_nist,
                                        compton_file, pair_file)  ) == NULL )
      {
         xvmc_error("MC_jaws_focus::init",
                    "cannot create left jaws bar",8);
      }

      // create right jaws bar, a region bounded by 6 planes
      p_ref.x = (x_min+x_max)/TWO;
      p_ref.y = (r_right+y_max)/TWO;
      p_ref.z = (z_min+z_max)/TWO;
      if ( (piece[1] = new (nothrow) MC_volume_6p(
                                        plane_right, plane_y_max,
                                        plane_x_min, plane_x_max,
                                        plane_z_min, plane_z_max,
                                        p_ref,
                                        bars_material_estar,
                                        bars_material_nist,
                                        compton_file, pair_file)  ) == NULL )
      {
         xvmc_error("MC_jaws_focus::init",
                    "cannot create right jaws bar",8);
      }

      // create intermediate region between the two jaw bars,
      // a region bounded by 6 planes
      p_ref.x = (x_min+x_max)/TWO;
      p_ref.y = (r_left+r_right)/TWO;
      p_ref.z = (z_min+z_max)/TWO;
      if ( (piece[2] = new (nothrow) MC_volume_6p(
                                        plane_left,  plane_right,
                                        plane_x_min, plane_x_max,
                                        plane_z_min, plane_z_max,
                                        p_ref,
                                        open_material_estar,
                                        open_material_nist,
                                        compton_file, pair_file)  ) == NULL )
      {
         xvmc_error("MC_jaws_focus::init",
                    "cannot create intermediate region",8);
      }
   }

   // set bit masks and bit patterns of all regions using the reference points
   set_bits();
}

// define x or y jaws pair with given opening, z boundaries,
// jaws bar material and opening material names
MC_jaws_focus::MC_jaws_focus(const axis &ini_type,
                             const real &ini_left,  const real &ini_right,
                             const real &ini_z_min, const real &ini_z_max,
                             const char *bars_material,
                             const char *open_material)
             : MC_object(8,3)
{
   // default origin to iso-center distance and focus position
   real ini_iso_dist = 100.0;
   real ini_z_focus  = ZERO;

   // default outer (x and y) object dimensions
   real ini_x_min = -15.0;
   real ini_x_max =  15.0;
   real ini_y_min = -15.0;
   real ini_y_max =  15.0;

   // initialize
   init(ini_type,  ini_left,  ini_right,
        ini_x_min, ini_x_max, ini_y_min, ini_y_max, ini_z_min, ini_z_max,
        ini_iso_dist, ini_z_focus, bars_material, open_material);
}

// delete focussing jaws
MC_jaws_focus::~MC_jaws_focus(void)
{
   // free memory for file names
   delete [] bars_material_estar; bars_material_estar = NULL;
   delete [] open_material_estar; open_material_estar = NULL;
   delete [] bars_material_nist;  bars_material_nist  = NULL;
   delete [] open_material_nist;  open_material_nist  = NULL;
   delete [] compton_file;        compton_file        = NULL;
   delete [] pair_file;           pair_file           = NULL;
}

// change left jaw position
void MC_jaws_focus::change_left(const real &new_left)
{
   // check outer object dimensions
   if (type == X)
   {
      if (new_left*z_max/iso_distance <= x_min)
      {
         xvmc_error("MC_jaws_focus::change_left",
                    "the left jaws opening is too large",8);
      }
   }
   else
   {
      if (new_left*z_max/iso_distance <= y_min)
      {
         xvmc_error("MC_jaws_focus::change_left",
                    "the left jaws opening is too large",8);
      }
   }

   // check and set jaw position
   if (new_left < open_right) open_left = new_left;
   else
   {
      xvmc_error("MC_jaws_focus::change_left",
                 "the left jaw position must be smaller than the right",8);
   }

   // the focus point
   real_3 p_focus;
   p_focus.x = ZERO; p_focus.y = ZERO; p_focus.z = z_focus;

   // we need 2 further points for the left jaws plane
   real_3 p1_left,p2_left;
   if (type == X)
   {
      p1_left.x = open_left; p1_left.y = -10.0; p1_left.z = iso_distance;
      p2_left.x = open_left; p2_left.y =  10.0; p2_left.z = iso_distance;
   }
   else
   {
      p1_left.x = -10.0; p1_left.y = open_left; p1_left.z = iso_distance;
      p2_left.x =  10.0; p2_left.y = open_left; p2_left.z = iso_distance;
   }

   // now move the left plane to the new position
   plane_left->set(p_focus,p1_left,p2_left);
}

// change right jaw position
void MC_jaws_focus::change_right(const real &new_right)
{
   // check outer object dimensions
   if (type == X)
   {
      if (new_right*z_max/iso_distance >= x_max)
      {
         xvmc_error("MC_jaws_focus::change_right",
                    "the right jaws opening is too large",8);
      }
   }
   else
   {
      if (new_right*z_max/iso_distance >= y_max)
      {
         xvmc_error("MC_jaws_focus::change_right",
                    "the right jaws opening is too large",8);
      }
   }

   // check and set jaw position
   if (open_left < new_right) open_right = new_right;
   else
   {
      xvmc_error("MC_jaws_focus::change_right",
                 "the right jaw position must be larger than the left",8);
   }

   // the focus point
   real_3 p_focus;
   p_focus.x = ZERO; p_focus.y = ZERO; p_focus.z = z_focus;

   // we need 2 further points for the right jaws plane
   real_3 p1_right,p2_right;
   if (type == X)
   {
      p1_right.x = open_right; p1_right.y = -10.0; p1_right.z = iso_distance;
      p2_right.x = open_right; p2_right.y =  10.0; p2_right.z = iso_distance;
   }
   else
   {
      p1_right.x = -10.0; p1_right.y = open_right; p1_right.z = iso_distance;
      p2_right.x =  10.0; p2_right.y = open_right; p2_right.z = iso_distance;
   }

   // now move the right plane to the new position
   plane_right->set(p_focus,p1_right,p2_right);
}
