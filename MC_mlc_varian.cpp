/*****************************************************************************
 * MC_mlc_varian.cpp:                                                        *
 *    class member functions for:                                            *
 *       MC_mlc_varian:  MLC for VARIAN medical linear accelerators          *
 *                                                                           *
 * Copyright (C) 2001    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding (adapted from ELEKTA MLC)            MF 10.07.2003      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <new>
using namespace std;

#include "MC_mlc_varian.h"

// ****************************************
// member functions of class MC_mlc_varian
// ****************************************

// initialize VARIAN MLC
void MC_mlc_varian::init(void)
{
   // set positions for tongue and groove
   z_groove_upper = z_min + 1.41;
   z_tongue_upper = z_min + 1.61;
   z_tongue_lower = z_min + 4.18;
   z_groove_lower = z_min + 4.38;
   if (z_groove_lower >= z_max)
   {
     xvmc_error("MC_mlc_varian::init",
     "the maximum z boundary must be larger than the lower groove position",8);
   }

   // set tongue and groove width (1.25 mm at the iso-center plane)
   width_tg = 0.125;

   // set width of the interleaf leakage gap (0.5 mm at the iso-center plane)
   width_gap = 0.05;

   // set distance of the perpendicular leaf shift (unknown)
   d_shift = 0.0;

   // set shift distance of the leaf side focus point (unknown)
   f_shift = 0.0;

   // calculate corresponding shift at the iso-center (the shift of
   // the leaf side focus point is caused by tilting the whole MLC)
   i_shift = f_shift*(ONE - TWO*iso_distance/(z_min+z_max));

   // check and set parameters for the leaf end curvature
   if (nominal_mlc->get_z_center_curve() == NULL)
   {
      xvmc_error("MC_mlc_varian::init",
                 "center z-position of leaf curvature not set",8);
   }
   if (nominal_mlc->get_radius_curve() == NULL)
   {
      xvmc_error("MC_mlc_varian::init",
                 "radius of leaf curvature not set",8);
   }
   z_center_curve = *nominal_mlc->get_z_center_curve();
   radius_curve   = *nominal_mlc->get_radius_curve();
   if ( (z_center_curve <= z_min) || (z_center_curve >= z_max) )
   {
      xvmc_error("MC_mlc_varian::init",
                 "the center z-position of leaf curvature is incorrect",8);
   }

   // counter for the object separator planes and object pieces
   unsigned n_separator = 0, n_piece = 0;

   // create the 6 outer planes
   plane_x_min = NULL;
   if ( (plane_x_min = new (nothrow) MC_plane(X,x_min,-1)) == NULL )
   {
      xvmc_error("MC_mlc_varian::init",
                 "cannot create minimum x plane",8);
   }
   separator[n_separator] = plane_x_min; ++n_separator;

   plane_x_max = NULL;
   if ( (plane_x_max = new (nothrow) MC_plane(X,x_max,-1)) == NULL )
   {
      xvmc_error("MC_mlc_varian::init",
                 "cannot create maximum x plane",8);
   }
   separator[n_separator] = plane_x_max; ++n_separator;

   plane_y_min = NULL;
   if ( (plane_y_min = new (nothrow) MC_plane(Y,y_min,-1)) == NULL )
   {
      xvmc_error("MC_mlc_varian::init",
                 "cannot create minimum y plane",8);
   }
   separator[n_separator] = plane_y_min; ++n_separator;

   plane_y_max = NULL;
   if ( (plane_y_max = new (nothrow) MC_plane(Y,y_max,-1)) == NULL )
   {
      xvmc_error("MC_mlc_varian::init",
                 "cannot create maximum y plane",8);
   }
   separator[n_separator] = plane_y_max; ++n_separator;

   // typically the starting plane of the particle transport
   plane_z_min = NULL;
   if ( (plane_z_min = new (nothrow) MC_plane(Z,z_min,-1)) == NULL )
   {
      xvmc_error("MC_mlc_varian::init",
                 "cannot create minimum z plane",8);
   }
   separator[n_separator] = plane_z_min; ++n_separator;
   starting_plane = plane_z_min;

   // typically the final plane of the particle transport
   plane_z_max = NULL;
   if ( (plane_z_max = new (nothrow) MC_plane(Z,z_max,-1)) == NULL )
   {
      xvmc_error("MC_mlc_varian::init",
                 "cannot create maximum z plane",8);
   }
   separator[n_separator] = plane_z_max; ++n_separator;
   final_plane = plane_z_max;

   // upper groove plane
   plane_z_groove_upper = NULL;
   if ( (plane_z_groove_upper = new (nothrow)
      MC_plane(Z,z_groove_upper,-1)) == NULL )
   {
      xvmc_error("MC_mlc_varian::init",
                 "cannot create upper groove z plane",8);
   }
   separator[n_separator] = plane_z_groove_upper; ++n_separator;

   // upper tongue plane
   plane_z_tongue_upper = NULL;
   if ( (plane_z_tongue_upper = new (nothrow)
      MC_plane(Z,z_tongue_upper,-1)) == NULL )
   {
      xvmc_error("MC_mlc_varian::init",
                 "cannot create upper tongue z plane",8);
   }
   separator[n_separator] = plane_z_tongue_upper; ++n_separator;

   // lower tongue plane
   plane_z_tongue_lower = NULL;
   if ( (plane_z_tongue_lower = new (nothrow)
      MC_plane(Z,z_tongue_lower,-1)) == NULL )
   {
      xvmc_error("MC_mlc_varian::init",
                 "cannot create lower tongue z plane",8);
   }
   separator[n_separator] = plane_z_tongue_lower; ++n_separator;

   // lower groove plane
   plane_z_groove_lower = NULL;
   if ( (plane_z_groove_lower = new (nothrow)
      MC_plane(Z,z_groove_lower,-1)) == NULL )
   {
      xvmc_error("MC_mlc_varian::init",
                 "cannot create lower groove z plane",8);
   }
   separator[n_separator] = plane_z_groove_lower; ++n_separator;

   // define pointers to the previous and present (upper/lower and central)
   // wall planes
   MC_plane  *previous_plane_uplower = NULL;
   MC_plane  *previous_plane_central = NULL;
   MC_plane  *present_plane_uplower  = NULL;
   MC_plane  *present_plane_central  = NULL;

   // positions of the previous and present (upper/lower and central)
   // wall planes
   real previous_pos_uplower = ZERO;
   real previous_pos_central = ZERO;
   real present_pos_uplower  = ZERO;
   real present_pos_central  = ZERO;

   // define pointers to the left and right leaf edge planes
   MC_plane  *plane_left      = NULL;
   MC_plane  *plane_right     = NULL;

   // pointer to the present region
   MC_region *present_region  = NULL;

   // the focus points
   real_3 p_focus_wall,p_focus_edge;
   if (xytype == X)
   {
      p_focus_wall.x = ZERO;
      p_focus_wall.y = f_shift;
      p_focus_wall.z = z_focus_wall;
   }
   else
   {
      p_focus_wall.x = f_shift;
      p_focus_wall.y = ZERO;
      p_focus_wall.z = z_focus_wall;
   }
   p_focus_edge.x = ZERO;
   p_focus_edge.y = ZERO;
   p_focus_edge.z = z_focus_edge;

   // we need 4 further points for the upper/lower and central leaf walls
   real_3 p1_wall_uplower,p2_wall_uplower;
   real_3 p1_wall_central,p2_wall_central;
   p1_wall_uplower.x = ZERO;           p2_wall_uplower.x = ZERO;
   p1_wall_uplower.y = ZERO;           p2_wall_uplower.y = ZERO;
   p1_wall_uplower.z = iso_distance;   p2_wall_uplower.z = iso_distance;
   p1_wall_central.x = ZERO;           p2_wall_central.x = ZERO;
   p1_wall_central.y = ZERO;           p2_wall_central.y = ZERO;
   p1_wall_central.z = iso_distance;   p2_wall_central.z = iso_distance;

   // 6 further points for the left and right leaf ends
   real_3 p1_left,p2_left,p3_left,p1_right,p2_right,p3_right;
   p1_left.x  = ZERO;           p1_right.x = ZERO;
   p1_left.y  = ZERO;           p1_right.y = ZERO;
   p1_left.z  = z_center_curve; p1_right.z = z_center_curve;
   p2_left.x  = ZERO;           p2_right.x = ZERO;
   p2_left.y  = ZERO;           p2_right.y = ZERO;
   p2_left.z  = z_center_curve; p2_right.z = z_center_curve;
   p3_left.x  = ZERO;           p3_right.x = ZERO;
   p3_left.y  = ZERO;           p3_right.y = ZERO;
   p3_left.z  = z_center_curve; p3_right.z = z_center_curve;

   // define 5 reference points to identify the regions
   real_3 p_ref_upper;
   real_3 p_ref_joint_upper;
   real_3 p_ref_central;
   real_3 p_ref_joint_lower;
   real_3 p_ref_lower;

   p_ref_upper.x       = ZERO;
   p_ref_upper.y       = ZERO;
   // the upper reference plane
   p_ref_upper.z       = (z_min+z_groove_upper)/TWO;

   p_ref_joint_upper.x = ZERO;
   p_ref_joint_upper.y = ZERO;
   // the upper leaf joint reference plane
   p_ref_joint_upper.z = (z_groove_upper+z_tongue_upper)/TWO;

   p_ref_central.x = ZERO;
   p_ref_central.y = ZERO;
   // the central reference plane
   p_ref_central.z = (z_tongue_upper+z_tongue_lower)/TWO;

   p_ref_joint_lower.x = ZERO;
   p_ref_joint_lower.y = ZERO;
   // the lower leaf joint reference plane
   p_ref_joint_lower.z = (z_tongue_lower+z_groove_lower)/TWO;

   p_ref_lower.x = ZERO;
   p_ref_lower.y = ZERO;
   // the lower reference plane
   p_ref_lower.z = (z_groove_lower+z_max)/TWO;

   // create regions
   if (xytype == X)
   {                // ************ X MLC ***************** //

      // set upper/lower MLC wall points
      present_pos_uplower = nominal_mlc->get_lperp(0)
                          + i_shift + d_shift - width_gap;
      p1_wall_uplower.x = -10.0;
      p1_wall_uplower.y = present_pos_uplower;
      p1_wall_uplower.z = iso_distance;
      p2_wall_uplower.x = 10.0;
      p2_wall_uplower.y = present_pos_uplower;
      p2_wall_uplower.z = iso_distance;

      // create upper/lower MLC wall plane
      present_plane_uplower = NULL;
      if ( (present_plane_uplower = new (nothrow)
         MC_plane(p_focus_wall,p1_wall_uplower,p2_wall_uplower,-1)) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create upper/lower MLC wall plane",8);
      }
      separator[n_separator] = present_plane_uplower; ++n_separator;

      // set central MLC wall points
      present_pos_central = nominal_mlc->get_lperp(0)
                          + i_shift + d_shift + width_tg;
      p1_wall_central.x = -10.0;
      p1_wall_central.y = present_pos_central;
      p1_wall_central.z = iso_distance;
      p2_wall_central.x = 10.0;
      p2_wall_central.y = present_pos_central;
      p2_wall_central.z = iso_distance;

      // create central MLC wall plane
      present_plane_central = NULL;
      if ( (present_plane_central = new (nothrow)
         MC_plane(p_focus_wall,p1_wall_central,p2_wall_central,-1)) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create central MLC wall plane",8);
      }
      separator[n_separator] = present_plane_central; ++n_separator;

      // set reference point for the first MLC wall (upper part)
      p_ref_upper.x = ZERO;
      p_ref_upper.y = (y_min + present_pos_uplower)/TWO - f_shift;
      p_ref_upper.y = p_ref_upper.y*p_ref_upper.z/iso_distance + f_shift;

      // create upper part of first MLC wall, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_x_min, plane_x_max,
                                    plane_y_min, present_plane_uplower,
                                    plane_z_min, plane_z_tongue_upper,
                                    p_ref_upper,
                                    leaf_material) ) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create upper part of the first MLC wall",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set reference point for the first MLC wall (central part)
      p_ref_central.x = ZERO;
      p_ref_central.y = (y_min + present_pos_central)/TWO - f_shift;
      p_ref_central.y = p_ref_central.y*p_ref_central.z/iso_distance + f_shift;

      // create central part of first MLC wall, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_x_min, plane_x_max,
                                    plane_y_min, present_plane_central,
                                    plane_z_tongue_upper, plane_z_tongue_lower,
                                    p_ref_central,
                                    leaf_material) ) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create central part of the first MLC wall",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set reference point for the first MLC wall (lower part)
      p_ref_lower.x = ZERO;
      p_ref_lower.y = (y_min + present_pos_uplower)/TWO - f_shift;
      p_ref_lower.y = p_ref_lower.y*p_ref_lower.z/iso_distance + f_shift;

      // create lower part of first MLC wall, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_x_min, plane_x_max,
                                    plane_y_min, present_plane_uplower,
                                    plane_z_tongue_lower, plane_z_max,
                                    p_ref_lower,
                                    leaf_material) ) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create lower part of the first MLC wall",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set first upper/lower leaf wall points
      previous_pos_uplower = present_pos_uplower;
      present_pos_uplower  = nominal_mlc->get_lperp(0) + i_shift + d_shift;
      p1_wall_uplower.x = -10.0;
      p1_wall_uplower.y = present_pos_uplower;
      p1_wall_uplower.z = iso_distance;
      p2_wall_uplower.x = 10.0;
      p2_wall_uplower.y = present_pos_uplower;
      p2_wall_uplower.z = iso_distance;

      // create first upper/lower leaf wall plane
      previous_plane_uplower = present_plane_uplower;
      present_plane_uplower  = NULL;
      if ( (present_plane_uplower = new (nothrow)
         MC_plane(p_focus_wall,p1_wall_uplower,p2_wall_uplower,-1)) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create first upper/lower leaf wall plane",8);
      }
      separator[n_separator] = present_plane_uplower; ++n_separator;

      // set first central leaf wall points
      previous_pos_central = present_pos_central;
      present_pos_central  = nominal_mlc->get_lperp(0)
                           + i_shift + d_shift + width_tg + width_gap;
      p1_wall_central.x = -10.0;
      p1_wall_central.y = present_pos_central;
      p1_wall_central.z = iso_distance;
      p2_wall_central.x = 10.0;
      p2_wall_central.y = present_pos_central;
      p2_wall_central.z = iso_distance;

      // create first central leaf wall plane
      previous_plane_central = present_plane_central;
      present_plane_central  = NULL;
      if ( (present_plane_central = new (nothrow)
         MC_plane(p_focus_wall,p1_wall_central,p2_wall_central,-1))==NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create first central leaf wall plane",8);
      }
      separator[n_separator] = present_plane_central; ++n_separator;

      // set reference point for the first upper air gap
      p_ref_upper.x = ZERO;
      p_ref_upper.y = (previous_pos_uplower + present_pos_uplower)/TWO
                    - f_shift;
      p_ref_upper.y = p_ref_upper.y*p_ref_upper.z/iso_distance + f_shift;

      // create first upper air gap, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                 plane_x_min, plane_x_max,
                                 previous_plane_uplower, present_plane_uplower,
                                 plane_z_min, plane_z_groove_upper,
                                 p_ref_upper,
                                 open_material) ) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create first upper air gap",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set reference point for the first upper leaf joint air gap
      p_ref_joint_upper.x = ZERO;
      p_ref_joint_upper.y = (previous_pos_uplower + present_pos_central)/TWO
                          - f_shift;
      p_ref_joint_upper.y = p_ref_joint_upper.y*p_ref_joint_upper.z/iso_distance                          + f_shift;

      // create first upper leaf joint air gap, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                 plane_x_min, plane_x_max,
                                 previous_plane_uplower, present_plane_central,
                                 plane_z_groove_upper, plane_z_tongue_upper,
                                 p_ref_joint_upper,
                                 open_material) ) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create first upper leaf joint air gap",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set reference point for the first central air gap
      p_ref_central.x = ZERO;
      p_ref_central.y = (previous_pos_central + present_pos_central)/TWO
                      - f_shift;
      p_ref_central.y = p_ref_central.y*p_ref_central.z/iso_distance + f_shift;

      // create first central air gap, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                 plane_x_min, plane_x_max,
                                 previous_plane_central, present_plane_central,
                                 plane_z_tongue_upper, plane_z_tongue_lower,
                                 p_ref_central,
                                 open_material) ) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create first central air gap",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set reference point for the first lower leaf joint air gap
      p_ref_joint_lower.x = ZERO;
      p_ref_joint_lower.y = (previous_pos_uplower + present_pos_central)/TWO
                          - f_shift;
      p_ref_joint_lower.y = p_ref_joint_lower.y*p_ref_joint_lower.z/iso_distance                          + f_shift;

      // create first lower leaf joint air gap, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                 plane_x_min, plane_x_max,
                                 previous_plane_uplower, present_plane_central,
                                 plane_z_tongue_lower, plane_z_groove_lower,
                                 p_ref_joint_lower,
                                 open_material) ) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create first lower leaf joint air gap",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set reference point for the first lower air gap
      p_ref_lower.x = ZERO;
      p_ref_lower.y = (previous_pos_uplower + present_pos_uplower)/TWO
                    - f_shift;
      p_ref_lower.y = p_ref_lower.y*p_ref_lower.z/iso_distance + f_shift;

      // create first lower air gap, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                 plane_x_min, plane_x_max,
                                 previous_plane_uplower, present_plane_uplower,
                                 plane_z_groove_lower, plane_z_max,
                                 p_ref_lower,
                                 open_material) ) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create first lower air gap",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // start loop for all leaf pairs
      for (int i_pair=0; i_pair<num_pairs; ++i_pair)
      {
         // set next upper/lower leaf wall points
         previous_pos_uplower = present_pos_uplower;
         present_pos_uplower  = nominal_mlc->get_lperp(i_pair+1)
                              + i_shift + d_shift - width_gap;
         p1_wall_uplower.x = -10.0;
         p1_wall_uplower.y = present_pos_uplower;
         p1_wall_uplower.z = iso_distance;
         p2_wall_uplower.x = 10.0;
         p2_wall_uplower.y = present_pos_uplower;
         p2_wall_uplower.z = iso_distance;

         // create next upper/lower leaf wall plane
         previous_plane_uplower = present_plane_uplower;
         present_plane_uplower  = NULL;
         if ( (present_plane_uplower = new (nothrow)
            MC_plane(p_focus_wall,p1_wall_uplower,p2_wall_uplower,-1))==NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create next upper/lower leaf wall plane",8);
         }
         separator[n_separator] = present_plane_uplower; ++n_separator;

         // set next central leaf wall points
         previous_pos_central = present_pos_central;
         present_pos_central  = nominal_mlc->get_lperp(i_pair+1)
                              + i_shift + d_shift + width_tg;
         p1_wall_central.x = -10.0;
         p1_wall_central.y = present_pos_central;
         p1_wall_central.z = iso_distance;
         p2_wall_central.x = 10.0;
         p2_wall_central.y = present_pos_central;
         p2_wall_central.z = iso_distance;

         // create next central leaf wall plane
         previous_plane_central = present_plane_central;
         present_plane_central  = NULL;
         if ( (present_plane_central = new (nothrow)
            MC_plane(p_focus_wall,p1_wall_central,p2_wall_central,-1))==NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create next central leaf wall plane",8);
         }
         separator[n_separator] = present_plane_central; ++n_separator;

         // left leaf position at the iso-center
         real open_left = nominal_mlc->get_left(i_pair);
         real tan_left  = open_left/iso_distance;
         real c_left    = z_center_curve*tan_left
                        - radius_curve*sqrt(ONE+tan_left*tan_left);
         real t_left    = c_left + radius_curve;
         p1_left.x = c_left; p1_left.y =  0.0;
         p2_left.x = c_left; p2_left.y = 10.0;
         p3_left.x = t_left; p3_left.y =  0.0;

         // create left leaf edge plane
         plane_left = NULL;
         if ( (plane_left = new (nothrow)
            MC_halfcyl(p_focus_edge,p1_left,p2_left,p3_left,-1))==NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create left leaf edge plane",8);
         }
         separator[n_separator]   = plane_left; ++n_separator;
         left_edge_planes[i_pair] = plane_left;

         // right leaf position at the iso-center
         real open_right = nominal_mlc->get_right(i_pair);
         real tan_right  = open_right/iso_distance;
         real c_right    = z_center_curve*tan_right
                         + radius_curve*sqrt(ONE+tan_right*tan_right);
         real t_right    = c_right - radius_curve;
         p1_right.x = c_right; p1_right.y =  0.0;
         p2_right.x = c_right; p2_right.y = 10.0;
         p3_right.x = t_right; p3_right.y =  0.0;

         // create right leaf edge plane
         plane_right = NULL;
         if ( (plane_right = new (nothrow)
            MC_halfcyl(p_focus_edge,p1_right,p2_right,p3_right,-1))==NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create right leaf edge plane",8);
         }
         separator[n_separator]    = plane_right; ++n_separator;
         right_edge_planes[i_pair] = plane_right;

         // the left and right leaf positions
         // in the upper reference plane
         real rl_upper = open_left*p_ref_upper.z/iso_distance;
         real rr_upper = open_right*p_ref_upper.z/iso_distance;

         // in the upper leaf joint reference plane
         real rl_joint_upper = open_left*p_ref_joint_upper.z/iso_distance;
         real rr_joint_upper = open_right*p_ref_joint_upper.z/iso_distance;

         // in the central reference plane
         real rl_central = open_left*p_ref_central.z/iso_distance;
         real rr_central = open_right*p_ref_central.z/iso_distance;

         // in the lower leaf joint reference plane
         real rl_joint_lower = open_left*p_ref_joint_lower.z/iso_distance;
         real rr_joint_lower = open_right*p_ref_joint_lower.z/iso_distance;

         // in the lower reference plane
         real rl_lower = open_left*p_ref_lower.z/iso_distance;
         real rr_lower = open_right*p_ref_lower.z/iso_distance;

         // the two upper/lower leaf wall positions
         // in the upper reference plane
         real r1_upper = f_shift
            + (previous_pos_uplower-f_shift)*p_ref_upper.z/iso_distance;
         real r2_upper = f_shift
            + (present_pos_uplower-f_shift)*p_ref_upper.z/iso_distance;

         // the two upper leaf joint wall positions
         // in the upper joint reference plane
         real r1_joint_upper = f_shift
            + (previous_pos_central-f_shift)*p_ref_joint_upper.z/iso_distance;
         real r2_joint_upper = f_shift
            + (present_pos_uplower-f_shift)*p_ref_joint_upper.z/iso_distance;

         // the two central leaf wall positions in the central reference plane
         real r1_central = f_shift
            + (previous_pos_central-f_shift)*p_ref_central.z/iso_distance;
         real r2_central = f_shift
            + (present_pos_central-f_shift)*p_ref_central.z/iso_distance;

         // the two lower leaf joint wall positions
         // in the lower joint reference plane
         real r1_joint_lower = f_shift
            + (previous_pos_central-f_shift)*p_ref_joint_lower.z/iso_distance;
         real r2_joint_lower = f_shift
            + (present_pos_uplower-f_shift)*p_ref_joint_lower.z/iso_distance;

         // the two upper/lower leaf wall positions
         // in the lower reference plane
         real r1_lower = f_shift
            + (previous_pos_uplower-f_shift)*p_ref_lower.z/iso_distance;
         real r2_lower = f_shift
            + (present_pos_uplower-f_shift)*p_ref_lower.z/iso_distance;

         // set reference point for the left leaf (upper part)
         p_ref_upper.x = (x_min + rl_upper)/TWO;
         p_ref_upper.y = (r1_upper + r2_upper)/TWO;

         // create upper part of left leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_x_min,    plane_left,
                                       previous_plane_uplower,
                                       present_plane_uplower,
                                       plane_z_min,    plane_z_groove_upper,
                                       p_ref_upper,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create upper part of left leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the left upper leaf joint
         p_ref_joint_upper.x = (x_min + rl_joint_upper)/TWO;
         p_ref_joint_upper.y = (r1_joint_upper + r2_joint_upper)/TWO;

         // create upper leaf joint of left leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_x_min,    plane_left,
                                       previous_plane_central,
                                       present_plane_uplower,
                                       plane_z_groove_upper,
                                       plane_z_tongue_upper,
                                       p_ref_joint_upper,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create upper leaf joint of left leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the left leaf (central part)
         p_ref_central.x = (x_min + rl_central)/TWO;
         p_ref_central.y = (r1_central + r2_central)/TWO;

         // create central part of left leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_x_min,    plane_left,
                                       previous_plane_central,
                                       present_plane_central,
                                       plane_z_tongue_upper,
                                       plane_z_tongue_lower,
                                       p_ref_central,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create central part of left leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the left lower leaf joint
         p_ref_joint_lower.x = (x_min + rl_joint_lower)/TWO;
         p_ref_joint_lower.y = (r1_joint_lower + r2_joint_lower)/TWO;

         // create lower leaf joint of left leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_x_min,    plane_left,
                                       previous_plane_central,
                                       present_plane_uplower,
                                       plane_z_tongue_lower,
                                       plane_z_groove_lower,
                                       p_ref_joint_lower,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create lower leaf joint of left leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the left leaf (lower part)
         p_ref_lower.x = (x_min + rl_lower)/TWO;
         p_ref_lower.y = (r1_lower + r2_lower)/TWO;

         // create lower part of left leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_x_min,    plane_left,
                                       previous_plane_uplower,
                                       present_plane_uplower,
                                       plane_z_groove_lower, plane_z_max,
                                       p_ref_lower,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create lower part of left leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the right leaf (upper part)
         p_ref_upper.x = (rr_upper + x_max)/TWO;
         p_ref_upper.y = (r1_upper + r2_upper)/TWO;

         // create upper part of right leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_right,    plane_x_max,
                                       previous_plane_uplower,
                                       present_plane_uplower,
                                       plane_z_min,    plane_z_groove_upper,
                                       p_ref_upper,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create upper part of right leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the right upper leaf joint
         p_ref_joint_upper.x = (rr_joint_upper + x_max)/TWO;
         p_ref_joint_upper.y = (r1_joint_upper + r2_joint_upper)/TWO;

         // create upper leaf joint of right leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_right,    plane_x_max,
                                       previous_plane_central,
                                       present_plane_uplower,
                                       plane_z_groove_upper,
                                       plane_z_tongue_upper,
                                       p_ref_joint_upper,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create upper leaf joint of right leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the right leaf (central part)
         p_ref_central.x = (rr_central + x_max)/TWO;
         p_ref_central.y = (r1_central + r2_central)/TWO;

         // create lower part of right leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_right,    plane_x_max,
                                       previous_plane_central,
                                       present_plane_central,
                                       plane_z_tongue_upper,
                                       plane_z_tongue_lower,
                                       p_ref_central,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create central part of right leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the right lower leaf joint
         p_ref_joint_lower.x = (rr_joint_lower + x_max)/TWO;
         p_ref_joint_lower.y = (r1_joint_lower + r2_joint_lower)/TWO;

         // create lower leaf joint of right leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_right,    plane_x_max,
                                       previous_plane_central,
                                       present_plane_uplower,
                                       plane_z_tongue_lower,
                                       plane_z_groove_lower,
                                       p_ref_joint_lower,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create lower leaf joint of right leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the right leaf (lower part)
         p_ref_lower.x = (rr_lower + x_max)/TWO;
         p_ref_lower.y = (r1_lower + r2_lower)/TWO;

         // create lower part of right leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_right,    plane_x_max,
                                       previous_plane_uplower,
                                       present_plane_uplower,
                                       plane_z_groove_lower,
                                       plane_z_max,
                                       p_ref_lower,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create lower part of right leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the intermediate leaf region (upper part)
         p_ref_upper.x = (rl_upper + rr_upper)/TWO;
         p_ref_upper.y = (r1_upper + r2_upper)/TWO;

         // create upper part of intermediate leaf region,
         // a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_left,     plane_right,
                                       previous_plane_uplower,
                                       present_plane_uplower,
                                       plane_z_min,    plane_z_groove_upper,
                                       p_ref_upper,
                                       open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                     "cannot create upper part of intermediate leaf region",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the intermediate leaf region
         // (upper leaf joint)
         p_ref_joint_upper.x = (rl_joint_upper + rr_joint_upper)/TWO;
         p_ref_joint_upper.y = (r1_joint_upper + r2_joint_upper)/TWO;

         // create upper leaf joint of intermediate leaf region,
         // a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_left,     plane_right,
                                       previous_plane_central,
                                       present_plane_uplower,
                                       plane_z_groove_upper,
                                       plane_z_tongue_upper,
                                       p_ref_joint_upper,
                                       open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
               "cannot create upper leaf joint of intermediate leaf region",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the intermediate leaf region (central part)
         p_ref_central.x = (rl_central + rr_central)/TWO;
         p_ref_central.y = (r1_central + r2_central)/TWO;

         // create central part of intermediate leaf region,
         // a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_left,     plane_right,
                                       previous_plane_central,
                                       present_plane_central,
                                       plane_z_tongue_upper,
                                       plane_z_tongue_lower,
                                       p_ref_central,
                                       open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                 "cannot create central part of intermediate leaf region",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the intermediate leaf region
         // (lower leaf joint)
         p_ref_joint_lower.x = (rl_joint_lower + rr_joint_lower)/TWO;
         p_ref_joint_lower.y = (r1_joint_lower + r2_joint_lower)/TWO;

         // create lower leaf joint of intermediate leaf region,
         // a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_left,     plane_right,
                                       previous_plane_central,
                                       present_plane_uplower,
                                       plane_z_tongue_lower,
                                       plane_z_groove_lower,
                                       p_ref_joint_lower,
                                       open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
               "cannot create lower leaf joint of intermediate leaf region",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the intermediate leaf region (lower part)
         p_ref_lower.x = (rl_lower + rr_lower)/TWO;
         p_ref_lower.y = (r1_lower + r2_lower)/TWO;

         // create lower part of intermediate leaf region,
         // a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_left,     plane_right,
                                       previous_plane_uplower,
                                       present_plane_uplower,
                                       plane_z_groove_lower, plane_z_max,
                                       p_ref_lower,
                                       open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                     "cannot create lower part of intermediate leaf region",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set next upper/lower leaf wall points
         previous_pos_uplower = present_pos_uplower;
         present_pos_uplower  = nominal_mlc->get_lperp(i_pair+1)
                              + i_shift + d_shift;
         p1_wall_uplower.x = -10.0;
         p1_wall_uplower.y = present_pos_uplower;
         p1_wall_uplower.z = iso_distance;
         p2_wall_uplower.x = 10.0;
         p2_wall_uplower.y = present_pos_uplower;
         p2_wall_uplower.z = iso_distance;

         // create next upper/lower leaf wall plane
         previous_plane_uplower = present_plane_uplower;
         present_plane_uplower  = NULL;
         if ( (present_plane_uplower = new (nothrow)
            MC_plane(p_focus_wall,p1_wall_uplower,p2_wall_uplower,-1))==NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create next upper/lower leaf wall plane",8);
         }
         separator[n_separator] = present_plane_uplower; ++n_separator;

         // set next central leaf wall points
         previous_pos_central = present_pos_central;
         present_pos_central  = nominal_mlc->get_lperp(i_pair+1)
                              + i_shift + d_shift + width_tg + width_gap;
         p1_wall_central.x = -10.0;
         p1_wall_central.y = present_pos_central;
         p1_wall_central.z = iso_distance;
         p2_wall_central.x = 10.0;
         p2_wall_central.y = present_pos_central;
         p2_wall_central.z = iso_distance;

         // create next central leaf wall plane
         previous_plane_central = present_plane_central;
         present_plane_central  = NULL;
         if ( (present_plane_central = new (nothrow)
            MC_plane(p_focus_wall,p1_wall_central,p2_wall_central,-1))==NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create first central leaf wall plane",8);
         }
         separator[n_separator] = present_plane_central; ++n_separator;

         // set reference point for the next upper air gap
         p_ref_upper.x = ZERO;
         p_ref_upper.y = (previous_pos_uplower + present_pos_uplower)/TWO
                       - f_shift;
         p_ref_upper.y = p_ref_upper.y*p_ref_upper.z/iso_distance + f_shift;

         // create next upper air gap, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_x_min, plane_x_max,
                                    previous_plane_uplower,
                                    present_plane_uplower,
                                    plane_z_min, plane_z_groove_upper,
                                    p_ref_upper,
                                    open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create next upper air gap",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the next upper leaf joint air gap
         p_ref_joint_upper.x = ZERO;
         p_ref_joint_upper.y = (previous_pos_uplower + present_pos_central)/TWO
                             - f_shift;
         p_ref_joint_upper.y = f_shift
            + p_ref_joint_upper.y*p_ref_joint_upper.z/iso_distance;

         // create next upper leaf joint air gap, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_x_min,    plane_x_max,
                                       previous_plane_uplower,
                                       present_plane_central,
                                       plane_z_groove_upper,
                                       plane_z_tongue_upper,
                                       p_ref_joint_upper,
                                       open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create next upper leaf joint air gap",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the next central air gap
         p_ref_central.x = ZERO;
         p_ref_central.y = (previous_pos_central + present_pos_central)/TWO
                         - f_shift;
         p_ref_central.y = p_ref_central.y*p_ref_central.z/iso_distance
                         + f_shift;

         // create next central air gap, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_x_min, plane_x_max,
                                    previous_plane_central,
                                    present_plane_central,
                                    plane_z_tongue_upper,
                                    plane_z_tongue_lower,
                                    p_ref_central,
                                    open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create next central air gap",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the next lower leaf joint air gap
         p_ref_joint_lower.x = ZERO;
         p_ref_joint_lower.y = (previous_pos_uplower + present_pos_central)/TWO
                             - f_shift;
         p_ref_joint_lower.y = f_shift
            + p_ref_joint_lower.y*p_ref_joint_lower.z/iso_distance;

         // create next lower leaf joint air gap, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_x_min,    plane_x_max,
                                       previous_plane_uplower,
                                       present_plane_central,
                                       plane_z_tongue_lower,
                                       plane_z_groove_lower,
                                       p_ref_joint_lower,
                                       open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create next lower leaf joint air gap",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the next lower air gap
         p_ref_lower.x = ZERO;
         p_ref_lower.y = (previous_pos_uplower + present_pos_uplower)/TWO
                       - f_shift;
         p_ref_lower.y = p_ref_lower.y*p_ref_lower.z/iso_distance + f_shift;

         // create next lower air gap, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_x_min, plane_x_max,
                                    previous_plane_uplower,
                                    present_plane_uplower,
                                    plane_z_groove_lower,
                                    plane_z_max,
                                    p_ref_lower,
                                    open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create next lower air gap",8);
         }
         piece[n_piece] = present_region; ++n_piece;

      } // end of leaf pair loop

      // set reference point for the last MLC wall (upper part)
      p_ref_upper.x = ZERO;
      p_ref_upper.y = (y_max + present_pos_uplower)/TWO - f_shift;
      p_ref_upper.y = p_ref_upper.y*p_ref_upper.z/iso_distance + f_shift;

      // create upper part of last MLC wall, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_x_min, plane_x_max,
                                    present_plane_uplower, plane_y_max,
                                    plane_z_min, plane_z_groove_upper,
                                    p_ref_upper,
                                    leaf_material) ) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create upper part of the last MLC wall",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set reference point for the last MLC wall (central part)
      p_ref_central.x = ZERO;
      p_ref_central.y = (y_max + present_pos_central)/TWO - f_shift;
      p_ref_central.y = p_ref_central.y*p_ref_central.z/iso_distance + f_shift;

      // create central part of last MLC wall, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_x_min, plane_x_max,
                                    present_plane_central, plane_y_max,
                                    plane_z_groove_upper,
                                    plane_z_groove_lower,
                                    p_ref_central,
                                    leaf_material) ) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create central part of the last MLC wall",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set reference point for the last MLC wall (lower part)
      p_ref_lower.x = ZERO;
      p_ref_lower.y = (y_max + present_pos_uplower)/TWO - f_shift;
      p_ref_lower.y = p_ref_lower.y*p_ref_lower.z/iso_distance + f_shift;

      // create lower part of last MLC wall, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_x_min, plane_x_max,
                                    present_plane_uplower, plane_y_max,
                                    plane_z_groove_lower,  plane_z_max,
                                    p_ref_lower,
                                    leaf_material) ) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create lower part of the last MLC wall",8);
      }
      piece[n_piece] = present_region; ++n_piece;

   }                // ******** end of X MLC ************** //
   else
   {                // ************ Y MLC ***************** //

      // set upper/lower MLC wall points
      present_pos_uplower = nominal_mlc->get_lperp(0)
                          + i_shift + d_shift - width_gap;
      p1_wall_uplower.x = present_pos_uplower;
      p1_wall_uplower.y = -10.0;
      p1_wall_uplower.z = iso_distance;
      p2_wall_uplower.x = present_pos_uplower;
      p2_wall_uplower.y = 10.0;
      p2_wall_uplower.z = iso_distance;

      // create upper/lower MLC wall plane
      present_plane_uplower = NULL;
      if ( (present_plane_uplower = new (nothrow)
         MC_plane(p_focus_wall,p1_wall_uplower,p2_wall_uplower,-1)) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create upper/lower MLC wall plane",8);
      }
      separator[n_separator] = present_plane_uplower; ++n_separator;

      // set central MLC wall points
      present_pos_central = nominal_mlc->get_lperp(0)
                          + i_shift + d_shift + width_tg;
      p1_wall_central.x = present_pos_central;
      p1_wall_central.y = -10.0;
      p1_wall_central.z = iso_distance;
      p2_wall_central.x = present_pos_central;
      p2_wall_central.y = 10.0;
      p2_wall_central.z = iso_distance;

      // create central MLC wall plane
      present_plane_central = NULL;
      if ( (present_plane_central = new (nothrow)
         MC_plane(p_focus_wall,p1_wall_central,p2_wall_central,-1)) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create central MLC wall plane",8);
      }
      separator[n_separator] = present_plane_central; ++n_separator;

      // set reference point for the first MLC wall (upper part)
      p_ref_upper.x = (x_min + present_pos_uplower)/TWO - f_shift;
      p_ref_upper.x = p_ref_upper.x*p_ref_upper.z/iso_distance + f_shift;
      p_ref_upper.y = ZERO;

      // create upper part of first MLC wall, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_y_min, plane_y_max,
                                    plane_x_min, present_plane_uplower,
                                    plane_z_min, plane_z_tongue_upper,
                                    p_ref_upper,
                                    leaf_material) ) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create upper part of the first MLC wall",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set reference point for the first MLC wall (central part)
      p_ref_central.x = (x_min + present_pos_central)/TWO - f_shift;
      p_ref_central.x = p_ref_central.x*p_ref_central.z/iso_distance + f_shift;
      p_ref_central.y = ZERO;

      // create central part of first MLC wall, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_y_min, plane_y_max,
                                    plane_x_min, present_plane_central,
                                    plane_z_tongue_upper, plane_z_tongue_lower,
                                    p_ref_central,
                                    leaf_material) ) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create central part of the first MLC wall",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set reference point for the first MLC wall (lower part)
      p_ref_lower.x = (x_min + present_pos_uplower)/TWO - f_shift;
      p_ref_lower.x = p_ref_lower.x*p_ref_lower.z/iso_distance + f_shift;
      p_ref_lower.y = ZERO;

      // create lower part of first MLC wall, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_y_min, plane_y_max,
                                    plane_x_min, present_plane_uplower,
                                    plane_z_tongue_lower, plane_z_max,
                                    p_ref_lower,
                                    leaf_material) ) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create lower part of the first MLC wall",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set first upper/lower leaf wall points
      previous_pos_uplower = present_pos_uplower;
      present_pos_uplower  = nominal_mlc->get_lperp(0) + i_shift + d_shift;
      p1_wall_uplower.x = present_pos_uplower;
      p1_wall_uplower.y = -10.0;
      p1_wall_uplower.z = iso_distance;
      p2_wall_uplower.x = present_pos_uplower;
      p2_wall_uplower.y = 10.0;
      p2_wall_uplower.z = iso_distance;

      // create first upper/lower leaf wall plane
      previous_plane_uplower = present_plane_uplower;
      present_plane_uplower  = NULL;
      if ( (present_plane_uplower = new (nothrow)
         MC_plane(p_focus_wall,p1_wall_uplower,p2_wall_uplower,-1)) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create first upper/lower leaf wall plane",8);
      }
      separator[n_separator] = present_plane_uplower; ++n_separator;

      // set first central leaf wall points
      previous_pos_central = present_pos_central;
      present_pos_central  = nominal_mlc->get_lperp(0)
                           + i_shift + d_shift + width_tg + width_gap;
      p1_wall_central.x = present_pos_central;
      p1_wall_central.y = -10.0;
      p1_wall_central.z = iso_distance;
      p2_wall_central.x = present_pos_central;
      p2_wall_central.y = 10.0;
      p2_wall_central.z = iso_distance;

      // create first central leaf wall plane
      previous_plane_central = present_plane_central;
      present_plane_central  = NULL;
      if ( (present_plane_central = new (nothrow)
         MC_plane(p_focus_wall,p1_wall_central,p2_wall_central,-1)) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create first central leaf wall plane",8);
      }
      separator[n_separator] = present_plane_central; ++n_separator;

      // set reference point for the first upper air gap
      p_ref_upper.x = (previous_pos_uplower + present_pos_uplower)/TWO
                    - f_shift;
      p_ref_upper.x = p_ref_upper.x*p_ref_upper.z/iso_distance + f_shift;
      p_ref_upper.y = ZERO;

      // create first upper air gap, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                 plane_y_min, plane_y_max,
                                 previous_plane_uplower, present_plane_uplower,
                                 plane_z_min, plane_z_groove_upper,
                                 p_ref_upper,
                                 open_material) ) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create first upper air gap",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set reference point for the first upper leaf joint air gap
      p_ref_joint_upper.x = (previous_pos_uplower + present_pos_central)/TWO
                          - f_shift;
      p_ref_joint_upper.x = p_ref_joint_upper.x*p_ref_joint_upper.z/iso_distance                          + f_shift;
      p_ref_joint_upper.y = ZERO;

      // create first upper leaf joint air gap, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                 plane_y_min, plane_y_max,
                                 previous_plane_uplower, present_plane_central,
                                 plane_z_groove_upper, plane_z_tongue_upper,
                                 p_ref_joint_upper,
                                 open_material) ) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create first upper leaf joint air gap",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set reference point for the first central air gap
      p_ref_central.x = (previous_pos_central + present_pos_central)/TWO
                      - f_shift;
      p_ref_central.x = p_ref_central.x*p_ref_central.z/iso_distance + f_shift;
      p_ref_central.y = ZERO;

      // create first central air gap, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                 plane_y_min, plane_y_max,
                                 previous_plane_central, present_plane_central,
                                 plane_z_tongue_upper, plane_z_tongue_lower,
                                 p_ref_central,
                                 open_material) ) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create first central air gap",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set reference point for the first lower leaf joint air gap
      p_ref_joint_lower.x = (previous_pos_uplower + present_pos_central)/TWO
                          - f_shift;
      p_ref_joint_lower.x = p_ref_joint_lower.x*p_ref_joint_lower.z/iso_distance                          + f_shift;
      p_ref_joint_lower.y = ZERO;

      // create first lower leaf joint air gap, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                 plane_y_min, plane_y_max,
                                 previous_plane_uplower, present_plane_central,
                                 plane_z_tongue_lower, plane_z_groove_lower,
                                 p_ref_joint_lower,
                                 open_material) ) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create first lower leaf joint air gap",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set reference point for the first lower air gap
      p_ref_lower.x = (previous_pos_uplower + present_pos_uplower)/TWO
                    - f_shift;
      p_ref_lower.x = p_ref_lower.x*p_ref_lower.z/iso_distance + f_shift;
      p_ref_lower.y = ZERO;

      // create first lower air gap, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                 plane_y_min, plane_y_max,
                                 previous_plane_uplower, present_plane_uplower,
                                 plane_z_groove_lower, plane_z_max,
                                 p_ref_lower,
                                 open_material) ) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create first lower air gap",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // start loop for all leaf pairs
      for (int i_pair=0; i_pair<num_pairs; ++i_pair)
      {
         // set next upper/lower leaf wall points
         previous_pos_uplower = present_pos_uplower;
         present_pos_uplower  = nominal_mlc->get_lperp(i_pair+1)
                              + i_shift + d_shift - width_gap;
         p1_wall_uplower.x = present_pos_uplower;
         p1_wall_uplower.y = -10.0;
         p1_wall_uplower.z = iso_distance;
         p2_wall_uplower.x = present_pos_uplower;
         p2_wall_uplower.y = 10.0;
         p2_wall_uplower.z = iso_distance;

         // create next upper/lower leaf wall plane
         previous_plane_uplower = present_plane_uplower;
         present_plane_uplower  = NULL;
         if ( (present_plane_uplower = new (nothrow)
            MC_plane(p_focus_wall,p1_wall_uplower,p2_wall_uplower,-1))==NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create next upper/lower leaf wall plane",8);
         }
         separator[n_separator] = present_plane_uplower; ++n_separator;

         // set next central wall points
         previous_pos_central = present_pos_central;
         present_pos_central  = nominal_mlc->get_lperp(i_pair+1)
                              + i_shift + d_shift + width_tg;
         p1_wall_central.x = present_pos_central;
         p1_wall_central.y = -10.0;
         p1_wall_central.z = iso_distance;
         p2_wall_central.x = present_pos_central;
         p2_wall_central.y = 10.0;
         p2_wall_central.z = iso_distance;

         // create next central leaf wall plane
         previous_plane_central = present_plane_central;
         present_plane_central  = NULL;
         if ( (present_plane_central = new (nothrow)
            MC_plane(p_focus_wall,p1_wall_central,p2_wall_central,-1))==NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create next central leaf wall plane",8);
         }
         separator[n_separator] = present_plane_central; ++n_separator;

         // left leaf position at the iso-center
         real open_left = nominal_mlc->get_left(i_pair);
         real tan_left  = open_left/iso_distance;
         real c_left    = z_center_curve*tan_left
                        - radius_curve*sqrt(ONE+tan_left*tan_left);
         real t_left    = c_left + radius_curve;
         p1_left.x =  0.0; p1_left.y = c_left;
         p2_left.x = 10.0; p2_left.y = c_left;
         p3_left.x =  0.0; p3_left.y = t_left;

         // create left leaf edge plane
         plane_left = NULL;
         if ( (plane_left = new (nothrow)
            MC_halfcyl(p_focus_edge,p1_left,p2_left,p3_left,-1))==NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create left leaf edge plane",8);
         }
         separator[n_separator]   = plane_left; ++n_separator;
         left_edge_planes[i_pair] = plane_left;

         // right leaf position at the iso-center
         real open_right = nominal_mlc->get_right(i_pair);
         real tan_right  = open_right/iso_distance;
         real c_right    = z_center_curve*tan_right
                         + radius_curve*sqrt(ONE+tan_right*tan_right);
         real t_right    = c_right - radius_curve;
         p1_right.x =  0.0; p1_right.y = c_right;
         p2_right.x = 10.0; p2_right.y = c_right;
         p3_right.x =  0.0; p3_right.y = t_right;

         // create right leaf edge plane
         plane_right = NULL;
         if ( (plane_right = new (nothrow)
            MC_halfcyl(p_focus_edge,p1_right,p2_right,p3_right,-1))==NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create right leaf edge plane",8);
         }
         separator[n_separator]    = plane_right; ++n_separator;
         right_edge_planes[i_pair] = plane_right;

         // the left and right leaf positions
         // in the upper reference plane
         real rl_upper = open_left*p_ref_upper.z/iso_distance;
         real rr_upper = open_right*p_ref_upper.z/iso_distance;

         // in the upper leaf joint reference plane
         real rl_joint_upper = open_left*p_ref_joint_upper.z/iso_distance;
         real rr_joint_upper = open_right*p_ref_joint_upper.z/iso_distance;

         // in the central reference plane
         real rl_central = open_left*p_ref_central.z/iso_distance;
         real rr_central = open_right*p_ref_central.z/iso_distance;

         // in the lower leaf joint reference plane
         real rl_joint_lower = open_left*p_ref_joint_lower.z/iso_distance;
         real rr_joint_lower = open_right*p_ref_joint_lower.z/iso_distance;

         // in the lower reference plane
         real rl_lower = open_left*p_ref_lower.z/iso_distance;
         real rr_lower = open_right*p_ref_lower.z/iso_distance;

         // the two upper/lower leaf wall positions
         // in the upper reference plane
         real r1_upper = f_shift
            + (previous_pos_uplower-f_shift)*p_ref_upper.z/iso_distance;
         real r2_upper = f_shift
            + (present_pos_uplower-f_shift)*p_ref_upper.z/iso_distance;

         // the two upper leaf joint wall positions
         // in the upper joint reference plane
         real r1_joint_upper = f_shift
            + (previous_pos_central-f_shift)*p_ref_joint_upper.z/iso_distance;
         real r2_joint_upper = f_shift
            + (present_pos_uplower-f_shift)*p_ref_joint_upper.z/iso_distance;

         // the two central leaf wall positions in the central reference plane
         real r1_central = f_shift
            + (previous_pos_central-f_shift)*p_ref_central.z/iso_distance;
         real r2_central = f_shift
            + (present_pos_central-f_shift)*p_ref_central.z/iso_distance;

         // the two lower leaf joint wall positions
         // in the lower joint reference plane
         real r1_joint_lower = f_shift
            + (previous_pos_central-f_shift)*p_ref_joint_lower.z/iso_distance;
         real r2_joint_lower = f_shift
            + (present_pos_uplower-f_shift)*p_ref_joint_lower.z/iso_distance;

         // the two upper/lower leaf wall positions
         // in the lower reference plane
         real r1_lower = f_shift
            + (previous_pos_uplower-f_shift)*p_ref_lower.z/iso_distance;
         real r2_lower = f_shift
            + (present_pos_uplower-f_shift)*p_ref_lower.z/iso_distance;

         // set reference point for the left leaf (upper part)
         p_ref_upper.x = (r1_upper + r2_upper)/TWO;
         p_ref_upper.y = (y_min + rl_upper)/TWO;

         // create upper part of left leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_y_min,    plane_left,
                                       previous_plane_uplower,
                                       present_plane_uplower,
                                       plane_z_min,    plane_z_groove_upper,
                                       p_ref_upper,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create upper part of left leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the left upper leaf joint
         p_ref_joint_upper.x = (r1_joint_upper + r2_joint_upper)/TWO;
         p_ref_joint_upper.y = (y_min + rl_joint_upper)/TWO;

         // create upper leaf joint of left leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_y_min,    plane_left,
                                       previous_plane_central,
                                       present_plane_uplower,
                                       plane_z_groove_upper,
                                       plane_z_tongue_upper,
                                       p_ref_joint_upper,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create upper leaf joint of left leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the left leaf (central part)
         p_ref_central.x = (r1_central + r2_central)/TWO;
         p_ref_central.y = (y_min + rl_central)/TWO;

         // create central part of left leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_y_min,    plane_left,
                                       previous_plane_central,
                                       present_plane_central,
                                       plane_z_tongue_upper,
                                       plane_z_tongue_lower,
                                       p_ref_central,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create central part of left leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the left lower leaf joint
         p_ref_joint_lower.x = (r1_joint_lower + r2_joint_lower)/TWO;
         p_ref_joint_lower.y = (y_min + rl_joint_lower)/TWO;

         // create lower leaf joint of left leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_y_min,    plane_left,
                                       previous_plane_central,
                                       present_plane_uplower,
                                       plane_z_tongue_lower,
                                       plane_z_groove_lower,
                                       p_ref_joint_lower,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create lower leaf joint of left leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the left leaf (lower part)
         p_ref_lower.x = (r1_lower + r2_lower)/TWO;
         p_ref_lower.y = (y_min + rl_lower)/TWO;

         // create lower part of left leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_y_min,    plane_left,
                                       previous_plane_uplower,
                                       present_plane_uplower,
                                       plane_z_groove_lower, plane_z_max,
                                       p_ref_lower,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create lower part of left leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the right leaf (upper part)
         p_ref_upper.x = (r1_upper + r2_upper)/TWO;
         p_ref_upper.y = (rr_upper + y_max)/TWO;

         // create upper part of right leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_right,    plane_y_max,
                                       previous_plane_uplower,
                                       present_plane_uplower,
                                       plane_z_min,    plane_z_groove_upper,
                                       p_ref_upper,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create upper part of right leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the right upper leaf joint
         p_ref_joint_upper.x = (r1_joint_upper + r2_joint_upper)/TWO;
         p_ref_joint_upper.y = (rr_joint_upper + y_max)/TWO;

         // create upper leaf joint of right leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_right,    plane_y_max,
                                       previous_plane_central,
                                       present_plane_uplower,
                                       plane_z_groove_upper,
                                       plane_z_tongue_upper,
                                       p_ref_joint_upper,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create upper leaf joint of right leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the right leaf (central part)
         p_ref_central.x = (r1_central + r2_central)/TWO;
         p_ref_central.y = (rr_central + y_max)/TWO;

         // create central part of right leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_right,    plane_y_max,
                                       previous_plane_central,
                                       present_plane_central,
                                       plane_z_tongue_upper,
                                       plane_z_tongue_lower,
                                       p_ref_central,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create central part of right leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the right lower leaf joint
         p_ref_joint_lower.x = (r1_joint_lower + r2_joint_lower)/TWO;
         p_ref_joint_lower.y = (rr_joint_lower + y_max)/TWO;

         // create lower leaf joint of right leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_right,    plane_y_max,
                                       previous_plane_central,
                                       present_plane_uplower,
                                       plane_z_tongue_lower,
                                       plane_z_groove_lower,
                                       p_ref_joint_lower,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create lower leaf joint of right leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the right leaf (lower part)
         p_ref_lower.x = (r1_lower + r2_lower)/TWO;
         p_ref_lower.y = (rr_lower + y_max)/TWO;

         // create lower part of right leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_right,    plane_y_max,
                                       previous_plane_uplower,
                                       present_plane_uplower,
                                       plane_z_groove_lower,
                                       plane_z_max,
                                       p_ref_lower,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create lower part of right leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the intermediate leaf region (upper part)
         p_ref_upper.x = (r1_upper + r2_upper)/TWO;
         p_ref_upper.y = (rl_upper + rr_upper)/TWO;

         // create upper part of intermediate leaf region,
         // a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_left,     plane_right,
                                       previous_plane_uplower,
                                       present_plane_uplower,
                                       plane_z_min,    plane_z_groove_upper,
                                       p_ref_upper,
                                       open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                     "cannot create upper part of intermediate leaf region",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the intermediate leaf region
         // (upper leaf joint)
         p_ref_joint_upper.x = (r1_joint_upper + r2_joint_upper)/TWO;
         p_ref_joint_upper.y = (rl_joint_upper + rr_joint_upper)/TWO;

         // create upper leaf joint of intermediate leaf region,
         // a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_left,     plane_right,
                                       previous_plane_central,
                                       present_plane_uplower,
                                       plane_z_groove_upper,
                                       plane_z_tongue_upper,
                                       p_ref_joint_upper,
                                       open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
               "cannot create upper leaf joint of intermediate leaf region",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the intermediate leaf region (central part)
         p_ref_central.x = (r1_central + r2_central)/TWO;
         p_ref_central.y = (rl_central + rr_central)/TWO;

         // create central part of intermediate leaf region,
         // a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_left,     plane_right,
                                       previous_plane_central,
                                       present_plane_central,
                                       plane_z_tongue_upper,
                                       plane_z_tongue_lower,
                                       p_ref_central,
                                       open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                 "cannot create central part of intermediate leaf region",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the intermediate leaf region
         // (lower leaf joint)
         p_ref_joint_lower.x = (r1_joint_lower + r2_joint_lower)/TWO;
         p_ref_joint_lower.y = (rl_joint_lower + rr_joint_lower)/TWO;

         // create lower leaf joint of intermediate leaf region,
         // a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_left,     plane_right,
                                       previous_plane_central,
                                       present_plane_uplower,
                                       plane_z_tongue_lower,
                                       plane_z_groove_lower,
                                       p_ref_joint_lower,
                                       open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
               "cannot create lower leaf joint of intermediate leaf region",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the intermediate leaf region (lower part)
         p_ref_lower.x = (r1_lower + r2_lower)/TWO;
         p_ref_lower.y = (rl_lower + rr_lower)/TWO;

         // create lower part of intermediate leaf region,
         // a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_left,     plane_right,
                                       previous_plane_uplower,
                                       present_plane_uplower,
                                       plane_z_groove_lower, plane_z_max,
                                       p_ref_lower,
                                       open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                     "cannot create lower part of intermediate leaf region",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set next upper/lower leaf wall points
         previous_pos_uplower = present_pos_uplower;
         present_pos_uplower  = nominal_mlc->get_lperp(i_pair+1)
                              + i_shift + d_shift;
         p1_wall_uplower.x = present_pos_uplower;
         p1_wall_uplower.y = -10.0;
         p1_wall_uplower.z = iso_distance;
         p2_wall_uplower.x = present_pos_uplower;
         p2_wall_uplower.y = 10.0;
         p2_wall_uplower.z = iso_distance;

         // create next upper/lower leaf wall plane
         previous_plane_uplower = present_plane_uplower;
         present_plane_uplower  = NULL;
         if ( (present_plane_uplower = new (nothrow)
            MC_plane(p_focus_wall,p1_wall_uplower,p2_wall_uplower,-1))==NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create next upper/lower leaf wall plane",8);
         }
         separator[n_separator] = present_plane_uplower; ++n_separator;

         // set next central leaf wall points
         previous_pos_central = present_pos_central;
         present_pos_central  = nominal_mlc->get_lperp(i_pair+1)
                              + i_shift + d_shift + width_tg + width_gap;
         p1_wall_central.x = present_pos_central;
         p1_wall_central.y = -10.0;
         p1_wall_central.z = iso_distance;
         p2_wall_central.x = present_pos_central;
         p2_wall_central.y = 10.0;
         p2_wall_central.z = iso_distance;

         // create next central leaf wall plane
         previous_plane_central = present_plane_central;
         present_plane_central  = NULL;
         if ( (present_plane_central = new (nothrow)
            MC_plane(p_focus_wall,p1_wall_central,p2_wall_central,-1))==NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create first central leaf wall plane",8);
         }
         separator[n_separator] = present_plane_central; ++n_separator;

         // set reference point for the next upper air gap
         p_ref_upper.x = (previous_pos_uplower + present_pos_uplower)/TWO
                       - f_shift;
         p_ref_upper.x = p_ref_upper.x*p_ref_upper.z/iso_distance + f_shift;
         p_ref_upper.y = ZERO;

         // create next upper air gap, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_y_min, plane_y_max,
                                    previous_plane_uplower,
                                    present_plane_uplower,
                                    plane_z_min, plane_z_groove_upper,
                                    p_ref_upper,
                                    open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create next upper air gap",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the next upper leaf joint air gap
         p_ref_joint_upper.x = (previous_pos_uplower + present_pos_central)/TWO
                             - f_shift;
         p_ref_joint_upper.x = f_shift
            + p_ref_joint_upper.x*p_ref_joint_upper.z/iso_distance;
         p_ref_joint_upper.y = ZERO;

         // create next upper leaf joint air gap, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_y_min,    plane_y_max,
                                       previous_plane_uplower,
                                       present_plane_central,
                                       plane_z_groove_upper,
                                       plane_z_tongue_upper,
                                       p_ref_joint_upper,
                                       open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create next upper leaf joint air gap",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the next central air gap
         p_ref_central.x = (previous_pos_central + present_pos_central)/TWO
                         - f_shift;
         p_ref_central.x = p_ref_central.x*p_ref_central.z/iso_distance
                         + f_shift;
         p_ref_central.y = ZERO;

         // create next central air gap, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_y_min, plane_y_max,
                                    previous_plane_central,
                                    present_plane_central,
                                    plane_z_tongue_upper,
                                    plane_z_tongue_lower,
                                    p_ref_central,
                                    open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create next central air gap",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the next lower leaf joint air gap
         p_ref_joint_lower.x = (previous_pos_uplower + present_pos_central)/TWO
                             - f_shift;
         p_ref_joint_lower.x = f_shift
            + p_ref_joint_lower.x*p_ref_joint_lower.z/iso_distance;
         p_ref_joint_lower.y = ZERO;

         // create next lower leaf joint air gap, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_y_min,    plane_y_max,
                                       previous_plane_uplower,
                                       present_plane_central,
                                       plane_z_tongue_lower,
                                       plane_z_groove_lower,
                                       p_ref_joint_lower,
                                       open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create next lower leaf joint air gap",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the next lower air gap
         p_ref_lower.x = (previous_pos_uplower + present_pos_uplower)/TWO
                       - f_shift;
         p_ref_lower.x = p_ref_lower.x*p_ref_lower.z/iso_distance + f_shift;
         p_ref_lower.y = ZERO;

         // create next lower air gap, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_y_min, plane_y_max,
                                    previous_plane_uplower,
                                    present_plane_uplower,
                                    plane_z_groove_lower,
                                    plane_z_max,
                                    p_ref_lower,
                                    open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_varian::init",
                       "cannot create next lower air gap",8);
         }
         piece[n_piece] = present_region; ++n_piece;

      } // end of leaf pair loop

      // set reference point for the last MLC wall (upper part)
      p_ref_upper.x = (x_max + present_pos_uplower)/TWO - f_shift;
      p_ref_upper.x = p_ref_upper.x*p_ref_upper.z/iso_distance + f_shift;
      p_ref_upper.y = ZERO;

      // create upper part of last MLC wall, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_y_min, plane_y_max,
                                    present_plane_uplower, plane_x_max,
                                    plane_z_min, plane_z_groove_upper,
                                    p_ref_upper,
                                    leaf_material) ) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create upper part of the last MLC wall",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set reference point for the last MLC wall (central part)
      p_ref_central.x = (x_max + present_pos_central)/TWO - f_shift;
      p_ref_central.x = p_ref_central.x*p_ref_central.z/iso_distance + f_shift;
      p_ref_central.y = ZERO;

      // create central part of last MLC wall, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_y_min, plane_y_max,
                                    present_plane_central, plane_x_max,
                                    plane_z_groove_upper,
                                    plane_z_groove_lower,
                                    p_ref_central,
                                    leaf_material) ) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create central part of the last MLC wall",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set reference point for the last MLC wall (lower part)
      p_ref_lower.x = (x_max + present_pos_uplower)/TWO - f_shift;
      p_ref_lower.x = p_ref_lower.x*p_ref_lower.z/iso_distance + f_shift;
      p_ref_lower.y = ZERO;

      // create lower part of last MLC wall, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_y_min, plane_y_max,
                                    present_plane_uplower, plane_x_max,
                                    plane_z_groove_lower,  plane_z_max,
                                    p_ref_lower,
                                    leaf_material) ) == NULL )
      {
         xvmc_error("MC_mlc_varian::init",
                    "cannot create lower part of the last MLC wall",8);
      }
      piece[n_piece] = present_region; ++n_piece;

   }                // ******** end of Y MLC ************** //

   // check number of planes and regions
   if (num_planes != n_separator)
   {
      xvmc_error("MC_mlc_varian::init",
                 "the number of object separator planes is incorrect",8);
   }
   if (num_regions != n_piece)
   {
      xvmc_error("MC_mlc_varian::init",
                 "the number of regions is incorrect",8);
   }

   // set bit masks and bit patterns of all regions using the reference points
   set_bits();

   return;
}

// change position of one leaf pair,
// the new leaf positions are defined at the iso-center plane
void MC_mlc_varian::change_pair(const int  &pair_index,
                                const real &new_left,   const real &new_right)
{
   // perform checks and set nominal MLC data
   MC_mlc::change_pair(pair_index, new_left, new_right);

   // the leaf edge focus point
   real_3 p_focus;
   p_focus.x = ZERO; p_focus.y = ZERO; p_focus.z = z_focus_edge;

   // four further points for the left and right leaf edges
   real_3 p1_left,p2_left,p1_right,p2_right;
   real tan_left  = new_left/iso_distance;
   real c_left    = z_center_curve*tan_left
                  - radius_curve*sqrt(ONE+tan_left*tan_left);
   real tan_right = new_right/iso_distance;
   real c_right   = z_center_curve*tan_right
                  + radius_curve*sqrt(ONE+tan_right*tan_right);
   if (xytype == X)
   {
      p1_left.x  = c_left;          p2_left.x  = c_left;
      p1_left.y  = 0.0;             p2_left.y  = 10.0;
      p1_left.z  = z_center_curve;  p2_left.z  = z_center_curve;
      p1_right.x = c_right;         p2_right.x = c_right;
      p1_right.y = 0.0;             p2_right.y = 10.0;
      p1_right.z = z_center_curve;  p2_right.z = z_center_curve;
   }
   else
   {
      p1_left.x  = 0.0;             p2_left.x  = 10.0;
      p1_left.y  = c_left;          p2_left.y  = c_left;
      p1_left.z  = z_center_curve;  p2_left.z  = z_center_curve;
      p1_right.x = 0.0;             p2_right.x = 10.0;
      p1_right.y = c_right;         p2_right.y = c_right;
      p1_right.z = z_center_curve;  p2_right.z = z_center_curve;
   }

   // now move the left and right leaf edge planes to the new positions
   left_edge_planes[pair_index]->set(p_focus,p1_left,p2_left);
   right_edge_planes[pair_index]->set(p_focus,p1_right,p2_right);

   return;
}
