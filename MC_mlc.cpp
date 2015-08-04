/*****************************************************************************
 * MC_mlc.cpp:                                                               *
 *    class member functions for:                                            *
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

// ****************************************
// includes
// ****************************************

#include <new>
using namespace std;

#include "MC_mlc.h"

// ****************************************
// member functions of class MC_mlc
// ****************************************

// initialize MLC
void MC_mlc::init(multi_leaf *ini_mlc,
                  const real &ini_x_min,     const real &ini_x_max,
                  const real &ini_y_min,     const real &ini_y_max,
                  const real &ini_iso_dist,
                  const real &ini_z_focus_wall,
                  const real &ini_z_focus_edge,
                  const char *open_material_name)
{
   // set MLC type
   mlctype = ini_mlc->get_type();

   // set MLC mode
   mlcmode = ini_mlc->get_mode();

   // check and set MLC leaf moving direction (X or Y)
   switch (ini_mlc->get_xy())
   {
   case 'x':  // leafes in X direction
   case 'X':
      xytype = X;
      break;
   case 'y':  // leafes in Y direction
   case 'Y':
      xytype = Y;
      break;
   default:
      xvmc_error("MC_mlc::init","there are only X or Y MLCs",8);
      break;
   }

   // determine leaf material name
   char *leaf_material_name = ini_mlc->get_material();

   // determine (ESTAR) material electron data file paths
   char *leaf_material_estar = NULL;
   if ( (leaf_material_estar=get_file_path(leaf_material_name,"estar"))==NULL )
   {
      xvmc_error("MC_mlc::init",
                 "cannot determine leaf material electron data file path",8);
   }

   char *open_material_estar = NULL;
   if ( (open_material_estar=get_file_path(open_material_name,"estar"))==NULL )
   {
      xvmc_error("MC_mlc::init",
         "cannot determine opening material electron data file path",8);
   }

   // determine (NIST) material cross section file paths
   char *leaf_material_nist = NULL;
   if ( (leaf_material_nist=get_file_path(leaf_material_name,"nist"))==NULL )
   {
      xvmc_error("MC_mlc::init",
                 "cannot determine leaf material cross section file path",8);
   }

   char *open_material_nist = NULL;
   if ( (open_material_nist=get_file_path(open_material_name,"nist"))==NULL )
   {
      xvmc_error("MC_mlc::init",
         "cannot determine opening material cross section file path",8);
   }

   // determine differential Compton cross section file path
   char *compton_file = NULL;
   if ( (compton_file = get_file_path("compton","data")) == NULL )
   {
      xvmc_error("MC_mlc::init",
         "cannot determine differential Compton cross section file path",8);
   }

   // determine differential pair cross section file path
   char *pair_file = NULL;
   if ( (pair_file = get_file_path("pair","data")) == NULL )
   {
      xvmc_error("MC_mlc::init",
         "cannot determine differential pair cross section file path",8);
   }

   // construct cross section data base for the leaf material
   leaf_material = NULL;
   if ( (leaf_material = new (nothrow)
         MC_material(leaf_material_estar, leaf_material_nist,
                     compton_file,        pair_file))          == NULL )
   {
      xvmc_error("MC_mlc::init",
         "cannot construct cross section data base for the leaf material",8);
   }

   // construct cross section data base for the opening material
   open_material = NULL;
   if ( (open_material = new (nothrow)
         MC_material(open_material_estar, open_material_nist,
                     compton_file,        pair_file))          == NULL )
   {
      xvmc_error("MC_mlc::init",
         "cannot construct cross section data base for the opening material",8);
   }

   // free memory for file names
   delete [] leaf_material_estar; leaf_material_estar = NULL;
   delete [] open_material_estar; open_material_estar = NULL;
   delete [] leaf_material_nist;  leaf_material_nist  = NULL;
   delete [] open_material_nist;  open_material_nist  = NULL;
   delete [] compton_file;        compton_file        = NULL;
   delete [] pair_file;           pair_file           = NULL;

   // check and set outer object dimensions in x direction
   if (ini_x_max > ini_x_min)
   {
      x_min = ini_x_min;
      x_max = ini_x_max;
   }
   else
   {
      xvmc_error("MC_mlc::init",
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
      xvmc_error("MC_mlc::init",
                 "the maximum y boundary must be larger than the minimum",8);
   }

   // check and set outer object dimensions in z direction
   if (ini_mlc->get_upper() == NULL)
   {
      xvmc_error("MC_mlc::init",
                 "upper MLC limit not set",8);
   }
   if (ini_mlc->get_lower() == NULL)
   {
      xvmc_error("MC_mlc::init",
                 "lower MLC limit not set",8);
   }
   real ini_z_min = *ini_mlc->get_upper();
   real ini_z_max = *ini_mlc->get_lower();
   if (ini_z_max > ini_z_min)
   {
      z_min = ini_z_min;
      z_max = ini_z_max;
   }
   else
   {
      xvmc_error("MC_mlc::init",
                 "the maximum z boundary must be larger than the minimum",8);
   }

   // origin to iso-center distance and focus positions
   iso_distance = ini_iso_dist;
   z_focus_wall = ini_z_focus_wall;
   z_focus_edge = ini_z_focus_edge;

   // number of leaf pairs
   num_pairs = ini_mlc->get_num();

   // check leaf positions
   real scale = z_max/iso_distance;
   if (xytype == X)
   {
      if ( scale*ini_mlc->get_lperp(0) <= y_min )
      {
         xvmc_error("MC_mlc::init",
                    "first leaf position exceeds object dimensions",8);
      }
      if ( scale*ini_mlc->get_lperp(num_pairs) >= y_max )
      {
         xvmc_error("MC_mlc::init",
                    "last leaf position exceeds object dimensions",8);
      }
      for (int i_pair=0; i_pair<num_pairs; ++i_pair)
      {
         real ini_left = ini_mlc->get_left(i_pair);
         if (ini_left*scale <= x_min)
         {
            xvmc_error("MC_mlc::init",
                       "the left leaf opening is too large",8);
         }
         real ini_right = ini_mlc->get_right(i_pair);
         if (ini_right*scale >= x_max)
         {
            xvmc_error("MC_mlc::init",
                       "the right leaf opening is too large",8);
         }
         if (ini_left >= ini_right)
         {
            xvmc_error("MC_mlc::init",
               "the right leaf position must be larger than the left",8);
         }
      }
   }
   else
   {
      if ( scale*ini_mlc->get_lperp(0) <= x_min )
      {
         xvmc_error("MC_mlc::init",
                    "first leaf position exceeds object dimensions",8);
      }
      if ( scale*ini_mlc->get_lperp(num_pairs) >= x_max )
      {
         xvmc_error("MC_mlc::init",
                    "last leaf position exceeds object dimensions",8);
      }
      for (int i_pair=0; i_pair<num_pairs; ++i_pair)
      {
         real ini_left = ini_mlc->get_left(i_pair);
         if (ini_left*scale <= y_min)
         {
            xvmc_error("MC_mlc::init",
                       "the left leaf opening is too large",8);
         }
         real ini_right = ini_mlc->get_right(i_pair);
         if (ini_right*scale >= y_max)
         {
            xvmc_error("MC_mlc::init",
                       "the right leaf opening is too large",8);
         }
         if (ini_left >= ini_right)
         {
            xvmc_error("MC_mlc::init",
               "the right leaf position must be larger than the left",8);
         }
      }
   }

   // construct nominal MLC data
   nominal_mlc = NULL;
   if ( (nominal_mlc = new (nothrow) multi_leaf(
                                      ini_mlc->get_type(),
                                      ini_mlc->get_mode(),
                                      ini_mlc->get_xy(),
                                      num_pairs,
                                      ini_mlc->get_material(),
                                      ini_mlc->get_upper(),
                                      ini_mlc->get_lower(),
                                      ini_mlc->get_z_center_curve(),
                                      ini_mlc->get_radius_curve())) == NULL)
   {
      xvmc_error("MC_mlc::init",
                 "cannot create nominal MLC data",8);
   }

   // set nominal leaf widths and positions
   for (register int ip=0; ip<num_pairs; ++ip)
   {
      real width = ini_mlc->get_width(ip);
      if (!nominal_mlc->change_width(ip,width))
      {
         xvmc_error("MC_mlc::init",
                    "cannot change nominal leaf width",8);
      }

      // present position
      real left  = ini_mlc->get_left(ip);
      real right = ini_mlc->get_right(ip);
      if (!nominal_mlc->change_pair(ip,left,right))
      {
         xvmc_error("MC_mlc::init",
                    "cannot change nominal leaf position",8);
      }

      // starting position
      real start_left  = ini_mlc->get_start_left(ip);
      real start_right = ini_mlc->get_start_right(ip);
      if (!nominal_mlc->change_start_pair(ip,start_left,start_right))
      {
         xvmc_error("MC_mlc::init",
                    "cannot change nominal leaf starting position",8);
      }

      // stopping position
      real stop_left  = ini_mlc->get_stop_left(ip);
      real stop_right = ini_mlc->get_stop_right(ip);
      if (!nominal_mlc->change_stop_pair(ip,stop_left,stop_right))
      {
         xvmc_error("MC_mlc::init",
                    "cannot change nominal leaf stopping position",8);
      }
   }

   // initialize pointers for the outer object planes
   plane_x_min = NULL;   plane_x_max = NULL;
   plane_y_min = NULL;   plane_y_max = NULL;
   plane_z_min = NULL;   plane_z_max = NULL;

   // allocate memory for the leaf edge plane pointer arrays
   left_edge_planes = NULL;
   if ( (left_edge_planes = new (nothrow) MC_plane*[num_pairs]) == NULL )
   {
      xvmc_error("MC_mlc::init",
      "cannot allocate memory for the left leaf edge plane pointer array",8);
   }
   right_edge_planes = NULL;
   if ( (right_edge_planes = new (nothrow) MC_plane*[num_pairs]) == NULL )
   {
      xvmc_error("MC_mlc::init",
      "cannot allocate memory for the right leaf edge plane pointer array",8);
   }

   // initialize leaf edge plane pointers by the NULL pointer
   for (register int i_pair=0; i_pair<num_pairs; ++i_pair)
   {
      left_edge_planes[i_pair]  = NULL;
      right_edge_planes[i_pair] = NULL;
   }

   return;
}

// define MLC using the nominal MLC, the opening material name,
// the source to iso-center distance as well as
// the numbers of object planes and regions
MC_mlc::MC_mlc(multi_leaf *ini_mlc,      const char *open_material_name,
               const real  ini_iso_dist,
               unsigned    n_planes,     unsigned    n_regions)
      : MC_object(n_planes,n_regions)
{
   // default focus positions
   real ini_z_focus_wall = ZERO;
   real ini_z_focus_edge = ZERO;

   // default outer (x and y) object dimensions
   real ini_x_min = -15.0;
   real ini_x_max =  15.0;
   real ini_y_min = -15.0;
   real ini_y_max =  15.0;

   // initialize
   init(ini_mlc, ini_x_min,     ini_x_max,
                 ini_y_min,     ini_y_max,
                 ini_iso_dist,  ini_z_focus_wall, ini_z_focus_edge,
                 open_material_name);
}

// delete MLC
MC_mlc::~MC_mlc(void)
{
   // memory for pointers to leaf edge planes
   delete [] left_edge_planes;    left_edge_planes    = NULL;
   delete [] right_edge_planes;   right_edge_planes   = NULL;

   // free memory for cross section data bases
   delete leaf_material; leaf_material = NULL;
   delete open_material; open_material = NULL;

   // free memory for nominal MLC data
   delete nominal_mlc; nominal_mlc = NULL;
}

// change position of one leaf pair,
// the new leaf positions are defined at the iso-center plane,
// this base class function only performs checks and changes the
// nominal MLC data, the real leaf positions must be changed by
// the corresponding derived class member function
void MC_mlc::change_pair(const int  &pair_index,
                         const real &new_left,
                         const real &new_right)
{
   // check leaf pair index
   if (pair_index < 0)
   {
      xvmc_error("MC_mlc::change_pair",
                 "leaf pair index less than zero",8);
   }
   if (pair_index >= num_pairs)
   {
      xvmc_error("MC_mlc::change_pair",
                 "leaf pair index to large",8);
   }

   // check leaf positions
   if (new_left >= new_right)
   {
      xvmc_error("MC_mlc::change_pair",
         "cannot change leaf positions, left position >= right position",8);
   }

   // check outer object dimensions
   if (xytype == X)
   {
      if (new_left*z_max/iso_distance <= x_min)
      {
         xvmc_error("MC_mlc::change_pair",
                    "the new left leaf opening is too large",8);
      }
      if (new_right*z_max/iso_distance >= x_max)
      {
         xvmc_error("MC_mlc::change_pair",
                    "the new right leaf opening is too large",8);
      }
   }
   else
   {
      if (new_left*z_max/iso_distance <= y_min)
      {
         xvmc_error("MC_mlc::change_pair",
                    "the new left leaf opening is too large",8);
      }
      if (new_right*z_max/iso_distance >= y_max)
      {
         xvmc_error("MC_mlc::change_pair",
                    "the new right leaf opening is too large",8);
      }
   }

   // set new nominal leaf positions
   if ( !nominal_mlc->change_pair(pair_index,new_left,new_right) )
   {
      xvmc_error("MC_mlc::change_pair",
                 "cannot set new nominal leaf positions",8);
   }

   return;
}

// change starting position of one leaf pair (iso-center plane),
// only the nominal MLC data will be changed
void MC_mlc::change_start_pair(const int  &pair_index,
                               const real &new_start_left,
                               const real &new_start_right)
{
   // check leaf pair index
   if (pair_index < 0)
   {
      xvmc_error("MC_mlc::change_start_pair",
                 "leaf pair index less than zero",8);
   }
   if (pair_index >= num_pairs)
   {
      xvmc_error("MC_mlc::change_start_pair",
                 "leaf pair index to large",8);
   }

   // check leaf positions
   if (new_start_left >= new_start_right)
   {
      xvmc_error("MC_mlc::change_start_pair",
         "cannot change leaf positions, left position >= right position",8);
   }

   // check outer object dimensions
   if (xytype == X)
   {
      if (new_start_left*z_max/iso_distance <= x_min)
      {
         xvmc_error("MC_mlc::change_start_pair",
                    "the new left leaf opening is too large",8);
      }
      if (new_start_right*z_max/iso_distance >= x_max)
      {
         xvmc_error("MC_mlc::change_start_pair",
                    "the new right leaf opening is too large",8);
      }
   }
   else
   {
      if (new_start_left*z_max/iso_distance <= y_min)
      {
         xvmc_error("MC_mlc::change_start_pair",
                    "the new left leaf opening is too large",8);
      }
      if (new_start_right*z_max/iso_distance >= y_max)
      {
         xvmc_error("MC_mlc::change_start_pair",
                    "the new right leaf opening is too large",8);
      }
   }

   // set new nominal leaf positions
   if ( !nominal_mlc->change_start_pair(pair_index,
                                        new_start_left,
                                        new_start_right) )
   {
      xvmc_error("MC_mlc::change_start_pair",
                 "cannot set new nominal leaf positions",8);
   }

   return;
}

// change stopping position of one leaf pair (iso-center plane),
// only the nominal MLC data will be changed
void MC_mlc::change_stop_pair(const int  &pair_index,
                              const real &new_stop_left,
                              const real &new_stop_right)
{
   // check leaf pair index
   if (pair_index < 0)
   {
      xvmc_error("MC_mlc::change_stop_pair",
                 "leaf pair index less than zero",8);
   }
   if (pair_index >= num_pairs)
   {
      xvmc_error("MC_mlc::change_stop_pair",
                 "leaf pair index to large",8);
   }

   // check leaf positions
   if (new_stop_left >= new_stop_right)
   {
      xvmc_error("MC_mlc::change_stop_pair",
         "cannot change leaf positions, left position >= right position",8);
   }

   // check outer object dimensions
   if (xytype == X)
   {
      if (new_stop_left*z_max/iso_distance <= x_min)
      {
         xvmc_error("MC_mlc::change_stop_pair",
                    "the new left leaf opening is too large",8);
      }
      if (new_stop_right*z_max/iso_distance >= x_max)
      {
         xvmc_error("MC_mlc::change_stop_pair",
                    "the new right leaf opening is too large",8);
      }
   }
   else
   {
      if (new_stop_left*z_max/iso_distance <= y_min)
      {
         xvmc_error("MC_mlc::change_stop_pair",
                    "the new left leaf opening is too large",8);
      }
      if (new_stop_right*z_max/iso_distance >= y_max)
      {
         xvmc_error("MC_mlc::change_stop_pair",
                    "the new right leaf opening is too large",8);
      }
   }

   // set new nominal leaf positions
   if ( !nominal_mlc->change_stop_pair(pair_index,
                                       new_stop_left,
                                       new_stop_right) )
   {
      xvmc_error("MC_mlc::change_stop_pair",
                 "cannot set new nominal leaf positions",8);
   }

   return;
}
