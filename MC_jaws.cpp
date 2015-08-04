/*****************************************************************************
 * MC_jaws.cpp:                                                              *
 *    class member functions for:                                            *
 *       MC_jaw:         physical collimator (JAW) base class                *
 *                                                                           *
 * Copyright (C) 2001    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 03.12.2001      *
 *    dynamic MLC                                         MF 23.10.2003      *
 *    Adopted from MC_mlc.cpp and MC_jaws_focus.cpp       JK Feb 22, 2008    *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <new>
using namespace std;

#include "MC_jaws.h"

// ****************************************
// member functions of class MC_jaw
// ****************************************

// initialize JAW
void MC_jaw::init(const axis &ini_type, 
                  const jaw_type    &int_jawtype, const jaw_mode    &ini_jawmode,
                  const real &ini_left,           const real &ini_right,
                  const real &ini_x_min,          const real &ini_x_max,
                  const real &ini_y_min,          const real &ini_y_max,
                  const real &ini_z_min,          const real &ini_z_max,
                  const real &ini_iso_dist,       const real &ini_z_focus,
                  const char *bars_material_name, 
                  const char *open_material_name)
{
   // set JAW type
   jawtype = get_type();

   // set JAW mode
   jawmode = get_mode();

   // check and set JAW moving direction (X or Y)
   switch (get_xy())
   {
   case 'x':  // jaw in X direction
   case 'X':
      xytype = X;
      break;
   case 'y':  // jaw in Y direction
   case 'Y':
      xytype = Y;
      break;
   default:
      xvmc_error("MC_jaw::init","there are only X or Y JAWs",8);
      break;
   }

   // determine (ESTAR) material electron data file paths
   char *bars_material_estar = NULL;
   if ( (bars_material_estar=get_file_path(bars_material_name,"estar"))==NULL )
   {
      xvmc_error("MC_jaw::init",
                 "cannot determine bars material electron data file path",8);
   }

   char *open_material_estar = NULL;
   if ( (open_material_estar=get_file_path(open_material_name,"estar"))==NULL )
   {
      xvmc_error("MC_jaw::init",
         "cannot determine opening material electron data file path",8);
   }

   // determine (NIST) material cross section file paths
   char *bars_material_nist = NULL;
   if ( (bars_material_nist=get_file_path(bars_material_name,"nist"))==NULL )
   {
      xvmc_error("MC_jaw::init",
                 "cannot determine bars material cross section file path",8);
   }

   char *open_material_nist = NULL;
   if ( (open_material_nist=get_file_path(open_material_name,"nist"))==NULL )
   {
      xvmc_error("MC_jaw::init",
         "cannot determine opening material cross section file path",8);
   }

   // determine differential Compton cross section file path
   char *compton_file = NULL;
   if ( (compton_file = get_file_path("compton","data")) == NULL )
   {
      xvmc_error("MC_jaw::init",
         "cannot determine differential Compton cross section file path",8);
   }

   // determine differential pair cross section file path
   char *pair_file = NULL;
   if ( (pair_file = get_file_path("pair","data")) == NULL )
   {
      xvmc_error("MC_jaw::init",
         "cannot determine differential pair cross section file path",8);
   }

   // construct cross section data base for the bars material
   bars_material = NULL;
   if ( (bars_material = new (nothrow)
         MC_material(bars_material_estar, bars_material_nist,
                     compton_file,        pair_file))          == NULL )
   {
      xvmc_error("MC_jaw::init",
         "cannot construct cross section data base for the bars material",8);
   }

   // construct cross section data base for the opening material
   open_material = NULL;
   if ( (open_material = new (nothrow)
         MC_material(open_material_estar, open_material_nist,
                     compton_file,        pair_file))          == NULL )
   {
      xvmc_error("MC_jaw::init",
         "cannot construct cross section data base for the opening material",8);
   }

   // free memory for file names
   delete [] bars_material_estar; bars_material_estar = NULL;
   delete [] open_material_estar; open_material_estar = NULL;
   delete [] bars_material_nist;  bars_material_nist  = NULL;
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
      xvmc_error("MC_jaw::init",
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
      xvmc_error("MC_jaw::init",
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
      xvmc_error("MC_jaw::init",
                 "the maximum z boundary must be larger than the minimum",8);
   }

   // origin to iso-center distance and focus positions
   iso_distance = ini_iso_dist;
   z_focus = ini_z_focus;

   // check and set jaw opening positions
   real scale = z_max/iso_distance;
   if (xytype == X)
   {
      if ( scale*ini_left <= y_min )
      {
         xvmc_error("MC_jaw::init",
                    "left jaw position in y direction exceeds object dimensions",8);
      }
      if ( scale*ini_right >= y_max )
      {
         xvmc_error("MC_jaw::init",
                    "right jaw position in y direction exceeds object dimensions",8);
      }
   }
   else
   {
      if ( scale*ini_left <= x_min )
      {
         xvmc_error("MC_jaw::init",
                    "left jaw position in x direction exceeds object dimensions",8);
      }
      if ( scale*ini_right >= x_max )
      {
         xvmc_error("MC_jaw::init",
                    "right jaw position in x direction exceeds object dimensions",8);
      }
   }

   if (ini_right > ini_left)
   {
      open_left  = ini_left;
      open_right = ini_right;
   }
   else
   {
      xvmc_error("MC_jaw::init",
                 "the right jaw position must be larger than the left",8);
   }

   // create the 6 outer planes
   plane_x_min = NULL;
   if ( (plane_x_min = new (nothrow) MC_plane(X,x_min,-1)) == NULL )
   {
      xvmc_error("MC_jaw::init",
                 "cannot create minimum x plane",8);
   }
   plane_x_max = NULL;
   if ( (plane_x_max = new (nothrow) MC_plane(X,x_max,-1)) == NULL )
   {
      xvmc_error("MC_jaw::init",
                 "cannot create maximum x plane",8);
   }
   plane_y_min = NULL;
   if ( (plane_y_min = new (nothrow) MC_plane(Y,y_min,-1)) == NULL )
   {
      xvmc_error("MC_jaw::init",
                 "cannot create minimum y plane",8);
   }
   plane_y_max = NULL;
   if ( (plane_y_max = new (nothrow) MC_plane(Y,y_max,-1)) == NULL )
   {
      xvmc_error("MC_jaw::init",
                 "cannot create maximum y plane",8);
   }

   // typically the starting plane of the particle transport
   plane_z_min = NULL;
   if ( (plane_z_min = new (nothrow) MC_plane(Z,z_min,-1)) == NULL )
   {
      xvmc_error("MC_jaw::init",
                 "cannot create minimum z plane",8);
   }
   starting_plane = plane_z_min;

   // typically the final plane of the particle transport
   plane_z_max = NULL;
   if ( (plane_z_max = new (nothrow) MC_plane(Z,z_max,-1)) == NULL )
   {
      xvmc_error("MC_jaw::init",
                 "cannot create maximum z plane",8);
   }
   final_plane = plane_z_max;

   // the focus point
   real_3 p_focus;
   p_focus.x = ZERO; p_focus.y = ZERO; p_focus.z = z_focus;

   // we need 2 further points for the left jaws plane
   real_3 p1_left,p2_left;
   if (xytype == X)
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
      xvmc_error("MC_jaw::init",
                 "cannot create left jaws plane",8);
   }

   // we need 2 further points for the right jaws plane
   real_3 p1_right,p2_right;
   if (xytype == X)
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
      xvmc_error("MC_jaw::init",
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
   if (xytype == X)   // X jaws
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
         xvmc_error("MC_jaw::init",
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
         xvmc_error("MC_jaw::init",
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
         xvmc_error("MC_jaw::init",
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
         xvmc_error("MC_jaw::init",
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
         xvmc_error("MC_jaw::init",
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
         xvmc_error("MC_jaw::init",
                    "cannot create intermediate region",8);
      }
   }

   // set bit masks and bit patterns of all regions using the reference points
   set_bits();

   // set previous jaw positions
   // current position
   real left  = ini_left;
   real right = ini_right;
   if (!change_opening(left,right))
   {
      xvmc_error("MC_jaw::init",
                 "cannot change nominal jaw opening position",8);
   }

   // starting position
   real start_left  = get_start_left();
   real start_right = get_start_right();
   if (!change_start_pair(start_left,start_right))
   {
      xvmc_error("MC_jaw::init",
                 "cannot change nominal jaw starting positions",8);
   }

   // stopping position
   real stop_left  = get_stop_left();
   real stop_right = get_stop_right();
   if (!change_stop_pair(stop_left,stop_right))
   {
      xvmc_error("MC_jaw::init",
                 "cannot change nominal jaw stopping positions",8);
   }

   // initialize jaw edge plane pointers by the NULL pointer
   plane_left  = NULL;
   plane_right = NULL;

   return;
}

// define JAW using the nominal JAW, the opening material name,
// the source to iso-center distance as well as
// the numbers of object planes and regions
MC_jaw::MC_jaw(const axis &ini_type,
               const jaw_type   &ini_jawtype,  const jaw_mode   &ini_jawmode,
               const real &ini_left,           const real &ini_right,
               const real &ini_z_min,          const real &ini_z_max,
               const char *bars_material,      const char *open_material)
             : MC_object(8,3)
{
   // default isocenter distance
   real ini_iso_dist = 100.0;
   // default focus positions
   real ini_z_focus = ZERO;

   // default outer (x and y) object dimensions
   real ini_x_min = -20.0;
   real ini_x_max =  20.0;
   real ini_y_min = -20.0;
   real ini_y_max =  20.0;

   // initialize
   init(ini_type, ini_jawtype, ini_jawmode, ini_left, ini_right,
        ini_x_min, ini_x_max, ini_y_min, ini_y_max, ini_z_min, ini_z_max,
        ini_iso_dist, ini_z_focus, bars_material, open_material);
}

// delete MLC
MC_jaw::~MC_jaw(void)
{
   // memory for pointers to jaw edge planes
   delete plane_left;    plane_left    = NULL;
   delete plane_right;   plane_right   = NULL;

   // free memory for cross section data bases
   delete bars_material; bars_material = NULL;
   delete open_material; open_material = NULL;

}

// change position of jaw pair,
// the new jaw positions are defined at the iso-center plane,
// this base class function only performs checks and changes the
// nominal JAW data, the real jaw positions must be changed by
// the corresponding derived class member function
bool MC_jaw::change_opening(const real &new_left,
                            const real &new_right)
{
   // check jaw positions
   if (new_left >= new_right)
   {
      xvmc_error("MC_jaw::change_opening",
         "cannot change jaw positions, left position >= right position",8);
   }

   // check outer object dimensions
   if (xytype == X)
   {
      if (new_left*z_max/iso_distance <= x_min)
      {
         xvmc_error("MC_jaw::change_opening",
                    "the new left jaw opening is too large",8);
      }
      if (new_right*z_max/iso_distance >= x_max)
      {
         xvmc_error("MC_jaw::change_opening",
                    "the new right jaw opening is too large",8);
      }
   }
   else
   {
      if (new_left*z_max/iso_distance <= y_min)
      {
         xvmc_error("MC_jaw::change_opening",
                    "the new left jaw opening is too large",8);
      }
      if (new_right*z_max/iso_distance >= y_max)
      {
         xvmc_error("MC_jaw::change_opening",
                    "the new right jaw opening is too large",8);
      }
   }

   // set new current jaw positions
   if ( !change_opening(new_left,new_right) )
   {
      xvmc_error("MC_jaw::change_opening",
                 "cannot set new nominal jaw positions",8);
   }

   // set current opeining
   open_left  = new_left;
   open_right = new_right;

   return (true);
}

// change left jaw position
bool MC_jaw::change_left(const real &new_left)
{
   // check outer object dimensions
   if (xytype == X)
   {
      if (new_left*z_max/iso_distance <= x_min)
      {
         xvmc_error("MC_jaw::change_left",
                    "the left jaws opening is too large",8);
      }
   }
   else
   {
      if (new_left*z_max/iso_distance <= y_min)
      {
         xvmc_error("MC_jaw::change_left",
                    "the left jaws opening is too large",8);
      }
   }

   // check and set jaw position
   if (new_left < open_right) open_left = new_left;
   else
   {
      xvmc_error("MC_jaw::change_left",
                 "the left jaw position must be smaller than the right",8);
   }

   // the focus point
   real_3 p_focus;
   p_focus.x = ZERO; p_focus.y = ZERO; p_focus.z = z_focus;

   // we need 2 further points for the left jaws plane
   real_3 p1_left,p2_left;
   if (xytype == X)
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

   // set current opeining
   open_left  = new_left;

   return (true);
}

// change right jaw position
bool MC_jaw::change_right(const real &new_right)
{
   // check outer object dimensions
   if (xytype == X)
   {
      if (new_right*z_max/iso_distance >= x_max)
      {
         xvmc_error("MC_jaw::change_right",
                    "the right jaws opening is too large",8);
      }
   }
   else
   {
      if (new_right*z_max/iso_distance >= y_max)
      {
         xvmc_error("MC_jaw::init::change_right",
                    "the right jaws opening is too large",8);
      }
   }

   // check and set jaw position
   if (open_left < new_right) open_right = new_right;
   else
   {
      xvmc_error("MC_jaw::init::change_right",
                 "the right jaw position must be larger than the left",8);
   }

   // the focus point
   real_3 p_focus;
   p_focus.x = ZERO; p_focus.y = ZERO; p_focus.z = z_focus;

   // we need 2 further points for the right jaws plane
   real_3 p1_right,p2_right;
   if (xytype == X)
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

   // set current opeining
   open_right = new_right;

   return(true);
}

// change starting position of jaw pair (iso-center plane),
// only the nominal JAW data will be changed
bool MC_jaw::change_start_pair(const real &new_start_left,
                               const real &new_start_right)
{
   // check jaw positions
   if (new_start_left >= new_start_right)
   {
      xvmc_error("MC_jaw::change_start_pair",
         "cannot change jaw positions, left position >= right position",8);
   }

   // set new jaw start positions
   if (!change_left(new_start_left))
   {
      xvmc_error("MC_jaw::change_start_pair",
                 "cannot set new left jaw start positions",8);
   }

   if (!change_right(new_start_left))
   {
      xvmc_error("MC_jaw::change_start_pair",
                 "cannot set new right jaw start positions",8);
   }

   return (true);
}

// change stopping position of jaw pair (iso-center plane),
// only the nominal MLC data will be changed
bool MC_jaw::change_stop_pair(const real &new_stop_left,
                              const real &new_stop_right)
{
   // check jaw positions
   if (new_stop_left >= new_stop_right)
   {
      xvmc_error("MC_jaw::change_stop_pair",
         "cannot change jaw positions, jaw position >= right position",8);
   }

   // set new jaw stop positions
   if (!change_left(new_stop_left))
   {
      xvmc_error("MC_jaw::change_stop_pair",
                 "cannot set new left jaw stop positions",8);
   }

   if (!change_right(new_stop_left))
   {
      xvmc_error("MC_jaw::change_stop_pair",
                 "cannot set new right jaw stop positions",8);
   }

   return (true);
}
