/*****************************************************************************
 * MC_mlc_elekta.cpp:                                                        *
 *    class member functions for:                                            *
 *       MC_mlc_elekta:  MLC for ELEKTA medical linear accelerators          *
 *                                                                           *
 * Copyright (C) 2001    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 29.11.2001      *
 *    implementation of interleaf leakage                 MF 03.03.2003      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <new>
using namespace std;

#include "MC_mlc_elekta.h"

// ****************************************
// member functions of class MC_mlc_elekta
// ****************************************

// initialize ELEKTA MLC
void MC_mlc_elekta::init(void)
{
   // set positions for tongue and groove
   z_tongue = z_min + 3.94;
   z_groove = z_min + 4.09;
   if (z_groove >= z_max)
   {
      xvmc_error("MC_mlc_elekta::init",
         "the maximum z boundary must be larger than the groove position",8);
   }

   // set tongue and groove width (1 mm at the iso-center plane)
   width_tg = 0.1;

   // set width of the interleaf leakage gap (0.775 mm at the iso-center plane)
   width_gap = 0.0775;

   // set distance of the perpendicular leaf shift (0.5 mm to the left)
   d_shift = -0.05;

   // set shift distance of the leaf side focus point (2 mm to the right)
   f_shift = 0.156;

   // calculate corresponding shift at the iso-center (the shift of
   // the leaf side focus point is caused by tilting the whole MLC)
   i_shift = f_shift*(ONE - TWO*iso_distance/(z_min+z_max));

   // check and set parameters for the leaf end curvature
   if (nominal_mlc->get_z_center_curve() == NULL)
   {
      xvmc_error("MC_mlc_elekta::init",
                 "center z-position of leaf curvature not set",8);
   }
   if (nominal_mlc->get_radius_curve() == NULL)
   {
      xvmc_error("MC_mlc_elekta::init",
                 "radius of leaf curvature not set",8);
   }
   z_center_curve = *nominal_mlc->get_z_center_curve();
   radius_curve   = *nominal_mlc->get_radius_curve();
   if ( (z_center_curve <= z_min) || (z_center_curve >= z_max) )
   {
      xvmc_error("MC_mlc_elekta::init",
                 "the center z-position of leaf curvature is incorrect",8);
   }

   // counter for the object separator planes and object pieces
   unsigned n_separator = 0, n_piece = 0;

   // create the 6 outer planes
   plane_x_min = NULL;
   if ( (plane_x_min = new (nothrow) MC_plane(X,x_min,-1)) == NULL )
   {
      xvmc_error("MC_mlc_elekta::init",
                 "cannot create minimum x plane",8);
   }
   separator[n_separator] = plane_x_min; ++n_separator;

   plane_x_max = NULL;
   if ( (plane_x_max = new (nothrow) MC_plane(X,x_max,-1)) == NULL )
   {
      xvmc_error("MC_mlc_elekta::init",
                 "cannot create maximum x plane",8);
   }
   separator[n_separator] = plane_x_max; ++n_separator;

   plane_y_min = NULL;
   if ( (plane_y_min = new (nothrow) MC_plane(Y,y_min,-1)) == NULL )
   {
      xvmc_error("MC_mlc_elekta::init",
                 "cannot create minimum y plane",8);
   }
   separator[n_separator] = plane_y_min; ++n_separator;

   plane_y_max = NULL;
   if ( (plane_y_max = new (nothrow) MC_plane(Y,y_max,-1)) == NULL )
   {
      xvmc_error("MC_mlc_elekta::init",
                 "cannot create maximum y plane",8);
   }
   separator[n_separator] = plane_y_max; ++n_separator;

   // typically the starting plane of the particle transport
   plane_z_min = NULL;
   if ( (plane_z_min = new (nothrow) MC_plane(Z,z_min,-1)) == NULL )
   {
      xvmc_error("MC_mlc_elekta::init",
                 "cannot create minimum z plane",8);
   }
   separator[n_separator] = plane_z_min; ++n_separator;
   starting_plane = plane_z_min;

   // typically the final plane of the particle transport
   plane_z_max = NULL;
   if ( (plane_z_max = new (nothrow) MC_plane(Z,z_max,-1)) == NULL )
   {
      xvmc_error("MC_mlc_elekta::init",
                 "cannot create maximum z plane",8);
   }
   separator[n_separator] = plane_z_max; ++n_separator;
   final_plane = plane_z_max;

   // tongue plane
   plane_z_tongue = NULL;
   if ( (plane_z_tongue = new (nothrow) MC_plane(Z,z_tongue,-1)) == NULL )
   {
      xvmc_error("MC_mlc_elekta::init",
                 "cannot create tongue z plane",8);
   }
   separator[n_separator] = plane_z_tongue; ++n_separator;

   // groove plane
   plane_z_groove = NULL;
   if ( (plane_z_groove = new (nothrow) MC_plane(Z,z_groove,-1)) == NULL )
   {
      xvmc_error("MC_mlc_elekta::init",
                 "cannot create groove z plane",8);
   }
   separator[n_separator] = plane_z_groove; ++n_separator;

   // define pointers to the previous and present (upper and lower) wall planes
   MC_plane  *previous_plane_upper  = NULL;
   MC_plane  *previous_plane_lower  = NULL;
   MC_plane  *present_plane_upper   = NULL;
   MC_plane  *present_plane_lower   = NULL;

   // positions of the previous and present (upper and lower) wall planes
   real previous_pos_upper = ZERO;
   real previous_pos_lower = ZERO;
   real present_pos_upper  = ZERO;
   real present_pos_lower  = ZERO;

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

   // we need 4 further points for the upper and lower leaf walls
   real_3 p1_wall_upper,p2_wall_upper,p1_wall_lower,p2_wall_lower;
   p1_wall_upper.x = ZERO;           p2_wall_upper.x = ZERO;
   p1_wall_upper.y = ZERO;           p2_wall_upper.y = ZERO;
   p1_wall_upper.z = iso_distance;   p2_wall_upper.z = iso_distance;
   p1_wall_lower.x = ZERO;           p2_wall_lower.x = ZERO;
   p1_wall_lower.y = ZERO;           p2_wall_lower.y = ZERO;
   p1_wall_lower.z = iso_distance;   p2_wall_lower.z = iso_distance;

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

   // define 3 reference points to identify the regions
   real_3 p_ref_upper,p_ref_joint,p_ref_lower;
   p_ref_upper.x = ZERO;
   p_ref_upper.y = ZERO;
   p_ref_upper.z = (z_min+z_tongue)/TWO;    // the upper reference plane
   p_ref_joint.x = ZERO;
   p_ref_joint.y = ZERO;
   p_ref_joint.z = (z_tongue+z_groove)/TWO; // the leaf joint reference plane
   p_ref_lower.x = ZERO;
   p_ref_lower.y = ZERO;
   p_ref_lower.z = (z_groove+z_max)/TWO;    // the lower reference plane

   // create regions
   if (xytype == X)
   {                // ************ X MLC ***************** //

      // set upper MLC wall points
      present_pos_upper = nominal_mlc->get_lperp(0)
                        + i_shift + d_shift - width_gap;
      p1_wall_upper.x = -10.0;             p2_wall_upper.x = 10.0;
      p1_wall_upper.y = present_pos_upper; p2_wall_upper.y = present_pos_upper;
      p1_wall_upper.z = iso_distance;      p2_wall_upper.z = iso_distance;

      // create upper MLC wall plane
      present_plane_upper = NULL;
      if ( (present_plane_upper = new (nothrow)
         MC_plane(p_focus_wall,p1_wall_upper,p2_wall_upper,-1)) == NULL )
      {
         xvmc_error("MC_mlc_elekta::init",
                    "cannot create upper MLC wall plane",8);
      }
      separator[n_separator] = present_plane_upper; ++n_separator;

      // set lower MLC wall points
      present_pos_lower = nominal_mlc->get_lperp(0)
                        + i_shift + d_shift + width_tg;
      p1_wall_lower.x = -10.0;             p2_wall_lower.x = 10.0;
      p1_wall_lower.y = present_pos_lower; p2_wall_lower.y = present_pos_lower;
      p1_wall_lower.z = iso_distance;      p2_wall_lower.z = iso_distance;

      // create lower MLC wall plane
      present_plane_lower = NULL;
      if ( (present_plane_lower = new (nothrow)
         MC_plane(p_focus_wall,p1_wall_lower,p2_wall_lower,-1)) == NULL )
      {
         xvmc_error("MC_mlc_elekta::init",
                    "cannot create lower MLC wall plane",8);
      }
      separator[n_separator] = present_plane_lower; ++n_separator;

      // set reference point for the first MLC wall (upper part)
      p_ref_upper.x = ZERO;
      p_ref_upper.y = (y_min + present_pos_upper)/TWO - f_shift;
      p_ref_upper.y = p_ref_upper.y*p_ref_upper.z/iso_distance + f_shift;

      // create upper part of first MLC wall, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_x_min, plane_x_max,
                                    plane_y_min, present_plane_upper,
                                    plane_z_min, plane_z_groove,
                                    p_ref_upper,
                                    leaf_material) ) == NULL )
      {
         xvmc_error("MC_mlc_elekta::init",
                    "cannot create upper part of the first MLC wall",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set reference point for the first MLC wall (lower part)
      p_ref_lower.x = ZERO;
      p_ref_lower.y = (y_min + present_pos_lower)/TWO - f_shift;
      p_ref_lower.y = p_ref_lower.y*p_ref_lower.z/iso_distance + f_shift;

      // create lower part of first MLC wall, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_x_min, plane_x_max,
                                    plane_y_min, present_plane_lower,
                                    plane_z_groove, plane_z_max,
                                    p_ref_lower,
                                    leaf_material) ) == NULL )
      {
         xvmc_error("MC_mlc_elekta::init",
                    "cannot create lower part of the first MLC wall",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set first upper leaf wall points
      previous_pos_upper = present_pos_upper;
      present_pos_upper  = nominal_mlc->get_lperp(0) + i_shift + d_shift;
      p1_wall_upper.x = -10.0;           p2_wall_upper.x = 10.0;
      p1_wall_upper.y=present_pos_upper; p2_wall_upper.y=present_pos_upper;
      p1_wall_upper.z = iso_distance;    p2_wall_upper.z = iso_distance;

      // create first upper leaf wall plane
      previous_plane_upper = present_plane_upper;
      present_plane_upper  = NULL;
      if ( (present_plane_upper = new (nothrow)
         MC_plane(p_focus_wall,p1_wall_upper,p2_wall_upper,-1))==NULL )
      {
         xvmc_error("MC_mlc_elekta::init",
                    "cannot create first upper leaf wall plane",8);
      }
      separator[n_separator] = present_plane_upper; ++n_separator;

      // set first lower leaf wall points
      previous_pos_lower = present_pos_lower;
      present_pos_lower  = nominal_mlc->get_lperp(0)
                         + i_shift + d_shift + width_tg + width_gap;
      p1_wall_lower.x = -10.0;           p2_wall_lower.x = 10.0;
      p1_wall_lower.y=present_pos_lower; p2_wall_lower.y=present_pos_lower;
      p1_wall_lower.z = iso_distance;    p2_wall_lower.z = iso_distance;

      // create first lower leaf wall plane
      previous_plane_lower = present_plane_lower;
      present_plane_lower  = NULL;
      if ( (present_plane_lower = new (nothrow)
         MC_plane(p_focus_wall,p1_wall_lower,p2_wall_lower,-1))==NULL )
      {
         xvmc_error("MC_mlc_elekta::init",
                    "cannot create first lower leaf wall plane",8);
      }
      separator[n_separator] = present_plane_lower; ++n_separator;

      // set reference point for the first upper air gap
      p_ref_upper.x = ZERO;
      p_ref_upper.y = (previous_pos_upper + present_pos_upper)/TWO - f_shift;
      p_ref_upper.y = p_ref_upper.y*p_ref_upper.z/iso_distance + f_shift;

      // create first upper air gap, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_x_min, plane_x_max,
                                    previous_plane_upper, present_plane_upper,
                                    plane_z_min, plane_z_tongue,
                                    p_ref_upper,
                                    open_material) ) == NULL )
      {
         xvmc_error("MC_mlc_elekta::init",
                    "cannot create first upper air gap",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set reference point for the first central air gap (at the leaf joint)
      p_ref_joint.x = ZERO;
      p_ref_joint.y = (previous_pos_upper + present_pos_lower)/TWO - f_shift;
      p_ref_joint.y = p_ref_joint.y*p_ref_joint.z/iso_distance + f_shift;

      // create first central air gap, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_x_min, plane_x_max,
                                    previous_plane_upper, present_plane_lower,
                                    plane_z_tongue, plane_z_groove,
                                    p_ref_joint,
                                    open_material) ) == NULL )
      {
         xvmc_error("MC_mlc_elekta::init",
                    "cannot create first central air gap",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set reference point for the first lower air gap
      p_ref_lower.x = ZERO;
      p_ref_lower.y = (previous_pos_lower + present_pos_lower)/TWO - f_shift;
      p_ref_lower.y = p_ref_lower.y*p_ref_lower.z/iso_distance + f_shift;

      // create first lower air gap, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_x_min, plane_x_max,
                                    previous_plane_lower, present_plane_lower,
                                    plane_z_groove, plane_z_max,
                                    p_ref_lower,
                                    open_material) ) == NULL )
      {
         xvmc_error("MC_mlc_elekta::init",
                    "cannot create first lower air gap",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // start loop for all leaf pairs
      for (int i_pair=0; i_pair<num_pairs; ++i_pair)
      {
         // set next upper leaf wall points
         previous_pos_upper = present_pos_upper;
         present_pos_upper  = nominal_mlc->get_lperp(i_pair+1)
                            + i_shift + d_shift - width_gap;
         p1_wall_upper.x = -10.0;           p2_wall_upper.x = 10.0;
         p1_wall_upper.y=present_pos_upper; p2_wall_upper.y=present_pos_upper;
         p1_wall_upper.z = iso_distance;    p2_wall_upper.z = iso_distance;

         // create next upper leaf wall plane
         previous_plane_upper = present_plane_upper;
         present_plane_upper  = NULL;
         if ( (present_plane_upper = new (nothrow)
            MC_plane(p_focus_wall,p1_wall_upper,p2_wall_upper,-1))==NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                       "cannot create next upper leaf wall plane",8);
         }
         separator[n_separator] = present_plane_upper; ++n_separator;

         // set next lower leaf wall points
         previous_pos_lower = present_pos_lower;
         present_pos_lower  = nominal_mlc->get_lperp(i_pair+1)
                            + i_shift + d_shift + width_tg;
         p1_wall_lower.x = -10.0;           p2_wall_lower.x = 10.0;
         p1_wall_lower.y=present_pos_lower; p2_wall_lower.y=present_pos_lower;
         p1_wall_lower.z = iso_distance;    p2_wall_lower.z = iso_distance;

         // create next lower leaf wall plane
         previous_plane_lower = present_plane_lower;
         present_plane_lower  = NULL;
         if ( (present_plane_lower = new (nothrow)
            MC_plane(p_focus_wall,p1_wall_lower,p2_wall_lower,-1))==NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                       "cannot create next lower leaf wall plane",8);
         }
         separator[n_separator] = present_plane_lower; ++n_separator;

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
            xvmc_error("MC_mlc_elekta::init",
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
            xvmc_error("MC_mlc_elekta::init",
                       "cannot create right leaf edge plane",8);
         }
         separator[n_separator]    = plane_right; ++n_separator;
         right_edge_planes[i_pair] = plane_right;

         // the left and right leaf positions
         // in the upper reference plane
         real rl_upper = open_left*p_ref_upper.z/iso_distance;
         real rr_upper = open_right*p_ref_upper.z/iso_distance;

         // in the leaf joint reference plane
         real rl_joint = open_left*p_ref_joint.z/iso_distance;
         real rr_joint = open_right*p_ref_joint.z/iso_distance;

         // in the lower reference plane
         real rl_lower = open_left*p_ref_lower.z/iso_distance;
         real rr_lower = open_right*p_ref_lower.z/iso_distance;

         // the two upper leaf wall positions in the upper reference plane
         real r1_upper =
            (previous_pos_upper-f_shift)*p_ref_upper.z/iso_distance + f_shift;
         real r2_upper =
            (present_pos_upper-f_shift)*p_ref_upper.z/iso_distance + f_shift;

         // the two leaf joint wall positions in the joint reference plane
         real r1_joint =
            (previous_pos_lower-f_shift)*p_ref_joint.z/iso_distance + f_shift;
         real r2_joint =
            (present_pos_upper-f_shift)*p_ref_joint.z/iso_distance + f_shift;

         // the two lower leaf wall positions in the lower reference plane
         real r1_lower =
            (previous_pos_lower-f_shift)*p_ref_lower.z/iso_distance + f_shift;
         real r2_lower =
            (present_pos_lower-f_shift)*p_ref_lower.z/iso_distance + f_shift;

         // set reference point for the left leaf (upper part)
         p_ref_upper.x = (x_min + rl_upper)/TWO;
         p_ref_upper.y = (r1_upper + r2_upper)/TWO;

         // create upper part of left leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_x_min,    plane_left,
                                       previous_plane_upper,
                                       present_plane_upper,
                                       plane_z_min,    plane_z_tongue,
                                       p_ref_upper,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                       "cannot create upper part of left leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the left leaf (leaf joint)
         p_ref_joint.x = (x_min + rl_joint)/TWO;
         p_ref_joint.y = (r1_joint + r2_joint)/TWO;

         // create leaf joint of left leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_x_min,    plane_left,
                                       previous_plane_lower,
                                       present_plane_upper,
                                       plane_z_tongue, plane_z_groove,
                                       p_ref_joint,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                       "cannot create leaf joint of left leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the left leaf (lower part)
         p_ref_lower.x = (x_min + rl_lower)/TWO;
         p_ref_lower.y = (r1_lower + r2_lower)/TWO;

         // create lower part of left leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_x_min,    plane_left,
                                       previous_plane_lower,
                                       present_plane_lower,
                                       plane_z_groove, plane_z_max,
                                       p_ref_lower,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
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
                                       previous_plane_upper,
                                       present_plane_upper,
                                       plane_z_min,    plane_z_tongue,
                                       p_ref_upper,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                       "cannot create upper part of right leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the right leaf (leaf joint)
         p_ref_joint.x = (rr_joint + x_max)/TWO;
         p_ref_joint.y = (r1_joint + r2_joint)/TWO;

         // create leaf joint of right leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_right,    plane_x_max,
                                       previous_plane_lower,
                                       present_plane_upper,
                                       plane_z_tongue, plane_z_groove,
                                       p_ref_joint,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                       "cannot create leaf joint of right leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the right leaf (lower part)
         p_ref_lower.x = (rr_lower + x_max)/TWO;
         p_ref_lower.y = (r1_lower + r2_lower)/TWO;

         // create lower part of right leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_right,    plane_x_max,
                                       previous_plane_lower,
                                       present_plane_lower,
                                       plane_z_groove, plane_z_max,
                                       p_ref_lower,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
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
                                       previous_plane_upper,
                                       present_plane_upper,
                                       plane_z_min,    plane_z_tongue,
                                       p_ref_upper,
                                       open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                     "cannot create upper part of intermediate leaf region",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the intermediate leaf region (leaf joint)
         p_ref_joint.x = (rl_joint + rr_joint)/TWO;
         p_ref_joint.y = (r1_joint + r2_joint)/TWO;

         // create leaf joint of intermediate leaf region,
         // a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_left,     plane_right,
                                       previous_plane_lower,
                                       present_plane_upper,
                                       plane_z_tongue, plane_z_groove,
                                       p_ref_joint,
                                       open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
               "cannot create leaf joint of intermediate leaf region",8);
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
                                       previous_plane_lower,
                                       present_plane_lower,
                                       plane_z_groove, plane_z_max,
                                       p_ref_lower,
                                       open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                     "cannot create lower part of intermediate leaf region",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set next upper leaf wall points
         previous_pos_upper = present_pos_upper;
         present_pos_upper  = nominal_mlc->get_lperp(i_pair+1)
                            + i_shift + d_shift;
         p1_wall_upper.x = -10.0;           p2_wall_upper.x = 10.0;
         p1_wall_upper.y=present_pos_upper; p2_wall_upper.y=present_pos_upper;
         p1_wall_upper.z = iso_distance;    p2_wall_upper.z = iso_distance;

         // create next upper leaf wall plane
         previous_plane_upper = present_plane_upper;
         present_plane_upper  = NULL;
         if ( (present_plane_upper = new (nothrow)
            MC_plane(p_focus_wall,p1_wall_upper,p2_wall_upper,-1))==NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                       "cannot create next upper leaf wall plane",8);
         }
         separator[n_separator] = present_plane_upper; ++n_separator;

         // set next lower leaf wall points
         previous_pos_lower = present_pos_lower;
         present_pos_lower  = nominal_mlc->get_lperp(i_pair+1)
                            + i_shift + d_shift + width_tg + width_gap;
         p1_wall_lower.x = -10.0;           p2_wall_lower.x = 10.0;
         p1_wall_lower.y=present_pos_lower; p2_wall_lower.y=present_pos_lower;
         p1_wall_lower.z = iso_distance;    p2_wall_lower.z = iso_distance;

         // create next lower leaf wall plane
         previous_plane_lower = present_plane_lower;
         present_plane_lower  = NULL;
         if ( (present_plane_lower = new (nothrow)
            MC_plane(p_focus_wall,p1_wall_lower,p2_wall_lower,-1))==NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                       "cannot create first lower leaf wall plane",8);
         }
         separator[n_separator] = present_plane_lower; ++n_separator;

         // set reference point for the next upper air gap
         p_ref_upper.x = ZERO;
         p_ref_upper.y = (previous_pos_upper + present_pos_upper)/TWO - f_shift;
         p_ref_upper.y = p_ref_upper.y*p_ref_upper.z/iso_distance + f_shift;

         // create next upper air gap, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_x_min, plane_x_max,
                                    previous_plane_upper, present_plane_upper,
                                    plane_z_min, plane_z_tongue,
                                    p_ref_upper,
                                    open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                       "cannot create next upper air gap",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the next central air gap
         p_ref_joint.x = ZERO;
         p_ref_joint.y = (previous_pos_upper + present_pos_lower)/TWO - f_shift;
         p_ref_joint.y = p_ref_joint.y*p_ref_joint.z/iso_distance + f_shift;

         // create next central air gap, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_x_min,    plane_x_max,
                                       previous_plane_upper,
                                       present_plane_lower,
                                       plane_z_tongue, plane_z_groove,
                                       p_ref_joint,
                                       open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                       "cannot create next central air gap",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the next lower air gap
         p_ref_lower.x = ZERO;
         p_ref_lower.y = (previous_pos_lower + present_pos_lower)/TWO - f_shift;
         p_ref_lower.y = p_ref_lower.y*p_ref_lower.z/iso_distance + f_shift;

         // create next lower air gap, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_x_min, plane_x_max,
                                    previous_plane_lower, present_plane_lower,
                                    plane_z_groove, plane_z_max,
                                    p_ref_lower,
                                    open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                       "cannot create next lower air gap",8);
         }
         piece[n_piece] = present_region; ++n_piece;

      } // end of leaf pair loop

      // set reference point for the last MLC wall (upper part)
      p_ref_upper.x = ZERO;
      p_ref_upper.y = (y_max + present_pos_upper)/TWO - f_shift;
      p_ref_upper.y = p_ref_upper.y*p_ref_upper.z/iso_distance + f_shift;

      // create upper part of last MLC wall, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_x_min, plane_x_max,
                                    present_plane_upper, plane_y_max,
                                    plane_z_min, plane_z_tongue,
                                    p_ref_upper,
                                    leaf_material) ) == NULL )
      {
         xvmc_error("MC_mlc_elekta::init",
                    "cannot create upper part of the last MLC wall",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set reference point for the last MLC wall (lower part)
      p_ref_lower.x = ZERO;
      p_ref_lower.y = (y_max + present_pos_lower)/TWO - f_shift;
      p_ref_lower.y = p_ref_lower.y*p_ref_lower.z/iso_distance + f_shift;

      // create lower part of last MLC wall, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_x_min, plane_x_max,
                                    present_plane_lower, plane_y_max,
                                    plane_z_tongue, plane_z_max,
                                    p_ref_lower,
                                    leaf_material) ) == NULL )
      {
         xvmc_error("MC_mlc_elekta::init",
                    "cannot create lower part of the last MLC wall",8);
      }
      piece[n_piece] = present_region; ++n_piece;

   }                // ******** end of X MLC ************** //
   else
   {                // ************ Y MLC ***************** //

      // set upper MLC wall points
      present_pos_upper = nominal_mlc->get_lperp(0)
                        + i_shift + d_shift - width_gap;
      p1_wall_upper.x = present_pos_upper; p2_wall_upper.x = present_pos_upper;
      p1_wall_upper.y = -10.0;             p2_wall_upper.y = 10.0;
      p1_wall_upper.z = iso_distance;      p2_wall_upper.z = iso_distance;

      // create upper MLC wall plane
      present_plane_upper = NULL;
      if ( (present_plane_upper = new (nothrow)
         MC_plane(p_focus_wall,p1_wall_upper,p2_wall_upper,-1)) == NULL )
      {
         xvmc_error("MC_mlc_elekta::init",
                    "cannot create upper MLC wall plane",8);
      }
      separator[n_separator] = present_plane_upper; ++n_separator;

      // set lower MLC wall points
      present_pos_lower = nominal_mlc->get_lperp(0)
                        + i_shift + d_shift + width_tg;
      p1_wall_lower.x = present_pos_lower; p2_wall_lower.x = present_pos_lower;
      p1_wall_lower.y = -10.0;             p2_wall_lower.y = 10.0;
      p1_wall_lower.z = iso_distance;      p2_wall_lower.z = iso_distance;

      // create lower MLC wall plane
      present_plane_lower = NULL;
      if ( (present_plane_lower = new (nothrow)
         MC_plane(p_focus_wall,p1_wall_lower,p2_wall_lower,-1)) == NULL )
      {
         xvmc_error("MC_mlc_elekta::init",
                    "cannot create lower MLC wall plane",8);
      }
      separator[n_separator] = present_plane_lower; ++n_separator;

      // set reference point for the first MLC wall (upper part)
      p_ref_upper.x = (x_min + present_pos_upper)/TWO - f_shift;
      p_ref_upper.x = p_ref_upper.x*p_ref_upper.z/iso_distance + f_shift;
      p_ref_upper.y = ZERO;

      // create upper part of first MLC wall, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_y_min, plane_y_max,
                                    plane_x_min, present_plane_upper,
                                    plane_z_min, plane_z_groove,
                                    p_ref_upper,
                                    leaf_material) ) == NULL )
      {
         xvmc_error("MC_mlc_elekta::init",
                    "cannot create upper part of the first MLC wall",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set reference point for the first MLC wall (lower part)
      p_ref_lower.x = (x_min + present_pos_lower)/TWO - f_shift;
      p_ref_lower.x = p_ref_lower.x*p_ref_lower.z/iso_distance + f_shift;
      p_ref_lower.y = ZERO;

      // create lower part of first MLC wall, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_y_min, plane_y_max,
                                    plane_x_min, present_plane_lower,
                                    plane_z_groove, plane_z_max,
                                    p_ref_lower,
                                    leaf_material) ) == NULL )
      {
         xvmc_error("MC_mlc_elekta::init",
                    "cannot create lower part of the first MLC wall",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set first upper leaf wall points
      previous_pos_upper = present_pos_upper;
      present_pos_upper  = nominal_mlc->get_lperp(0) + i_shift + d_shift;
      p1_wall_upper.x=present_pos_upper; p2_wall_upper.x=present_pos_upper;
      p1_wall_upper.y = -10.0;           p2_wall_upper.y = 10.0;
      p1_wall_upper.z = iso_distance;    p2_wall_upper.z = iso_distance;

      // create first upper leaf wall plane
      previous_plane_upper = present_plane_upper;
      present_plane_upper  = NULL;
      if ( (present_plane_upper = new (nothrow)
         MC_plane(p_focus_wall,p1_wall_upper,p2_wall_upper,-1))==NULL )
      {
         xvmc_error("MC_mlc_elekta::init",
                    "cannot create first upper leaf wall plane",8);
      }
      separator[n_separator] = present_plane_upper; ++n_separator;

      // set first lower leaf wall points
      previous_pos_lower = present_pos_lower;
      present_pos_lower  = nominal_mlc->get_lperp(0)
                         + i_shift + d_shift + width_tg + width_gap;
      p1_wall_lower.x=present_pos_lower; p2_wall_lower.x=present_pos_lower;
      p1_wall_lower.y = -10.0;           p2_wall_lower.y = 10.0;
      p1_wall_lower.z = iso_distance;    p2_wall_lower.z = iso_distance;

      // create first lower leaf wall plane
      previous_plane_lower = present_plane_lower;
      present_plane_lower  = NULL;
      if ( (present_plane_lower = new (nothrow)
         MC_plane(p_focus_wall,p1_wall_lower,p2_wall_lower,-1))==NULL )
      {
         xvmc_error("MC_mlc_elekta::init",
                    "cannot create first lower leaf wall plane",8);
      }
      separator[n_separator] = present_plane_lower; ++n_separator;

      // set reference point for the first upper air gap
      p_ref_upper.x = (previous_pos_upper + present_pos_upper)/TWO - f_shift;
      p_ref_upper.x = p_ref_upper.x*p_ref_upper.z/iso_distance + f_shift;
      p_ref_upper.y = ZERO;

      // create first upper air gap, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_y_min, plane_y_max,
                                    previous_plane_upper, present_plane_upper,
                                    plane_z_min, plane_z_tongue,
                                    p_ref_upper,
                                    open_material) ) == NULL )
      {
         xvmc_error("MC_mlc_elekta::init",
                    "cannot create first upper air gap",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set reference point for the first central air gap (at the leaf joint)
      p_ref_joint.x = (previous_pos_upper + present_pos_lower)/TWO - f_shift;
      p_ref_joint.x = p_ref_joint.x*p_ref_joint.z/iso_distance + f_shift;
      p_ref_joint.y = ZERO;

      // create first central air gap, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_y_min, plane_y_max,
                                    previous_plane_upper, present_plane_lower,
                                    plane_z_tongue, plane_z_groove,
                                    p_ref_joint,
                                    open_material) ) == NULL )
      {
         xvmc_error("MC_mlc_elekta::init",
                    "cannot create first central air gap",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set reference point for the first lower air gap
      p_ref_lower.x = (previous_pos_lower + present_pos_lower)/TWO - f_shift;
      p_ref_lower.x = p_ref_lower.x*p_ref_lower.z/iso_distance + f_shift;
      p_ref_lower.y = ZERO;

      // create first lower air gap, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_y_min, plane_y_max,
                                    previous_plane_lower, present_plane_lower,
                                    plane_z_groove, plane_z_max,
                                    p_ref_lower,
                                    open_material) ) == NULL )
      {
         xvmc_error("MC_mlc_elekta::init",
                    "cannot create first lower air gap",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // start loop for all leaf pairs
      for (int i_pair=0; i_pair<num_pairs; ++i_pair)
      {
         // set next upper leaf wall points
         previous_pos_upper = present_pos_upper;
         present_pos_upper  = nominal_mlc->get_lperp(i_pair+1)
                            + i_shift + d_shift - width_gap;
         p1_wall_upper.x=present_pos_upper; p2_wall_upper.x=present_pos_upper;
         p1_wall_upper.y = -10.0;           p2_wall_upper.y = 10.0;
         p1_wall_upper.z = iso_distance;    p2_wall_upper.z = iso_distance;

         // create next upper leaf wall plane
         previous_plane_upper = present_plane_upper;
         present_plane_upper  = NULL;
         if ( (present_plane_upper = new (nothrow)
            MC_plane(p_focus_wall,p1_wall_upper,p2_wall_upper,-1))==NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                       "cannot create next upper leaf wall plane",8);
         }
         separator[n_separator] = present_plane_upper; ++n_separator;

         // set next lower wall points
         previous_pos_lower = present_pos_lower;
         present_pos_lower  = nominal_mlc->get_lperp(i_pair+1)
                            + i_shift + d_shift + width_tg;
         p1_wall_lower.x=present_pos_lower; p2_wall_lower.x=present_pos_lower;
         p1_wall_lower.y = -10.0;           p2_wall_lower.y = 10.0;
         p1_wall_lower.z = iso_distance;    p2_wall_lower.z = iso_distance;

         // create next lower leaf wall plane
         previous_plane_lower = present_plane_lower;
         present_plane_lower  = NULL;
         if ( (present_plane_lower = new (nothrow)
            MC_plane(p_focus_wall,p1_wall_lower,p2_wall_lower,-1))==NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                       "cannot create next lower leaf wall plane",8);
         }
         separator[n_separator] = present_plane_lower; ++n_separator;

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
            xvmc_error("MC_mlc_elekta::init",
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
            xvmc_error("MC_mlc_elekta::init",
                       "cannot create right leaf edge plane",8);
         }
         separator[n_separator]    = plane_right; ++n_separator;
         right_edge_planes[i_pair] = plane_right;

         // the left and right leaf positions
         // in the upper reference plane
         real rl_upper = open_left*p_ref_upper.z/iso_distance;
         real rr_upper = open_right*p_ref_upper.z/iso_distance;

         // in the leaf joint reference plane
         real rl_joint = open_left*p_ref_joint.z/iso_distance;
         real rr_joint = open_right*p_ref_joint.z/iso_distance;

         // in the lower reference plane
         real rl_lower = open_left*p_ref_lower.z/iso_distance;
         real rr_lower = open_right*p_ref_lower.z/iso_distance;

         // the two upper leaf wall positions in the upper reference plane
         real r1_upper =
            (previous_pos_upper-f_shift)*p_ref_upper.z/iso_distance + f_shift;
         real r2_upper =
            (present_pos_upper-f_shift)*p_ref_upper.z/iso_distance + f_shift;

         // the two leaf joint wall positions in the joint reference plane
         real r1_joint =
            (previous_pos_lower-f_shift)*p_ref_joint.z/iso_distance + f_shift;
         real r2_joint =
            (present_pos_upper-f_shift)*p_ref_joint.z/iso_distance + f_shift;

         // the two lower leaf wall positions in the lower reference plane
         real r1_lower =
            (previous_pos_lower-f_shift)*p_ref_lower.z/iso_distance + f_shift;
         real r2_lower =
            (present_pos_lower-f_shift)*p_ref_lower.z/iso_distance + f_shift;

         // set reference point for the left leaf (upper part)
         p_ref_upper.x = (r1_upper + r2_upper)/TWO;
         p_ref_upper.y = (y_min + rl_upper)/TWO;

         // create upper part of left leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_y_min,    plane_left,
                                       previous_plane_upper,
                                       present_plane_upper,
                                       plane_z_min,    plane_z_tongue,
                                       p_ref_upper,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                       "cannot create upper part of left leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the left leaf (leaf joint)
         p_ref_joint.x = (r1_joint + r2_joint)/TWO;
         p_ref_joint.y = (y_min + rl_joint)/TWO;

         // create leaf joint of left leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_y_min,    plane_left,
                                       previous_plane_lower,
                                       present_plane_upper,
                                       plane_z_tongue, plane_z_groove,
                                       p_ref_joint,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                       "cannot create leaf joint of left leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the left leaf (lower part)
         p_ref_lower.x = (r1_lower + r2_lower)/TWO;
         p_ref_lower.y = (y_min + rl_lower)/TWO;

         // create lower part of left leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_y_min,    plane_left,
                                       previous_plane_lower,
                                       present_plane_lower,
                                       plane_z_groove, plane_z_max,
                                       p_ref_lower,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
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
                                       previous_plane_upper,
                                       present_plane_upper,
                                       plane_z_min,    plane_z_tongue,
                                       p_ref_upper,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                       "cannot create upper part of right leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the right leaf (leaf joint)
         p_ref_joint.x = (r1_joint + r2_joint)/TWO;
         p_ref_joint.y = (rr_joint + y_max)/TWO;

         // create leaf joint of right leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_right,    plane_y_max,
                                       previous_plane_lower,
                                       present_plane_upper,
                                       plane_z_tongue, plane_z_groove,
                                       p_ref_joint,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                       "cannot create leaf joint of right leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the right leaf (lower part)
         p_ref_lower.x = (r1_lower + r2_lower)/TWO;
         p_ref_lower.y = (rr_lower + y_max)/TWO;

         // create lower part of right leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_right,    plane_y_max,
                                       previous_plane_lower,
                                       present_plane_lower,
                                       plane_z_groove, plane_z_max,
                                       p_ref_lower,
                                       leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
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
                                       previous_plane_upper,
                                       present_plane_upper,
                                       plane_z_min,    plane_z_tongue,
                                       p_ref_upper,
                                       open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                     "cannot create upper part of intermediate leaf region",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the intermediate leaf region (leaf joint)
         p_ref_joint.x = (r1_joint + r2_joint)/TWO;
         p_ref_joint.y = (rl_joint + rr_joint)/TWO;

         // create leaf joint of intermediate leaf region,
         // a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_left,     plane_right,
                                       previous_plane_lower,
                                       present_plane_upper,
                                       plane_z_tongue, plane_z_groove,
                                       p_ref_joint,
                                       open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
               "cannot create leaf joint of intermediate leaf region",8);
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
                                       previous_plane_lower,
                                       present_plane_lower,
                                       plane_z_groove, plane_z_max,
                                       p_ref_lower,
                                       open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                     "cannot create lower part of intermediate leaf region",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set next upper leaf wall points
         previous_pos_upper = present_pos_upper;
         present_pos_upper  = nominal_mlc->get_lperp(i_pair+1)
                            + i_shift + d_shift;
         p1_wall_upper.x=present_pos_upper; p2_wall_upper.x=present_pos_upper;
         p1_wall_upper.y = -10.0;           p2_wall_upper.y = 10.0;
         p1_wall_upper.z = iso_distance;    p2_wall_upper.z = iso_distance;

         // create next upper leaf wall plane
         previous_plane_upper = present_plane_upper;
         present_plane_upper  = NULL;
         if ( (present_plane_upper = new (nothrow)
            MC_plane(p_focus_wall,p1_wall_upper,p2_wall_upper,-1))==NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                       "cannot create next upper leaf wall plane",8);
         }
         separator[n_separator] = present_plane_upper; ++n_separator;

         // set next lower leaf wall points
         previous_pos_lower = present_pos_lower;
         present_pos_lower  = nominal_mlc->get_lperp(i_pair+1)
                            + i_shift + d_shift + width_tg + width_gap;
         p1_wall_lower.x=present_pos_lower; p2_wall_lower.x=present_pos_lower;
         p1_wall_lower.y = -10.0;           p2_wall_lower.y = 10.0;
         p1_wall_lower.z = iso_distance;    p2_wall_lower.z = iso_distance;

         // create next lower leaf wall plane
         previous_plane_lower = present_plane_lower;
         present_plane_lower  = NULL;
         if ( (present_plane_lower = new (nothrow)
            MC_plane(p_focus_wall,p1_wall_lower,p2_wall_lower,-1))==NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                       "cannot create first lower leaf wall plane",8);
         }
         separator[n_separator] = present_plane_lower; ++n_separator;

         // set reference point for the next upper air gap
         p_ref_upper.x = (previous_pos_upper + present_pos_upper)/TWO - f_shift;
         p_ref_upper.x = p_ref_upper.x*p_ref_upper.z/iso_distance + f_shift;
         p_ref_upper.y = ZERO;

         // create next upper air gap, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_y_min, plane_y_max,
                                    previous_plane_upper, present_plane_upper,
                                    plane_z_min, plane_z_tongue,
                                    p_ref_upper,
                                    open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                       "cannot create next upper air gap",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the next central air gap
         p_ref_joint.x = (previous_pos_upper + present_pos_lower)/TWO - f_shift;
         p_ref_joint.x = p_ref_joint.x*p_ref_joint.z/iso_distance + f_shift;
         p_ref_joint.y = ZERO;

         // create next central air gap, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                       plane_y_min,    plane_y_max,
                                       previous_plane_upper,
                                       present_plane_lower,
                                       plane_z_tongue, plane_z_groove,
                                       p_ref_joint,
                                       open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                       "cannot create next central air gap",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the next lower air gap
         p_ref_lower.x = (previous_pos_lower + present_pos_lower)/TWO - f_shift;
         p_ref_lower.x = p_ref_lower.x*p_ref_lower.z/iso_distance + f_shift;
         p_ref_lower.y = ZERO;

         // create next lower air gap, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_y_min, plane_y_max,
                                    previous_plane_lower, present_plane_lower,
                                    plane_z_groove, plane_z_max,
                                    p_ref_lower,
                                    open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_elekta::init",
                       "cannot create next lower air gap",8);
         }
         piece[n_piece] = present_region; ++n_piece;

      } // end of leaf pair loop

      // set reference point for the last MLC wall (upper part)
      p_ref_upper.x = (x_max + present_pos_upper)/TWO - f_shift;
      p_ref_upper.x = p_ref_upper.x*p_ref_upper.z/iso_distance + f_shift;
      p_ref_upper.y = ZERO;

      // create upper part of last MLC wall, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_y_min, plane_y_max,
                                    present_plane_upper, plane_x_max,
                                    plane_z_min, plane_z_tongue,
                                    p_ref_upper,
                                    leaf_material) ) == NULL )
      {
         xvmc_error("MC_mlc_elekta::init",
                    "cannot create upper part of the last MLC wall",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // set reference point for the last MLC wall (lower part)
      p_ref_lower.x = (x_max + present_pos_lower)/TWO - f_shift;
      p_ref_lower.x = p_ref_lower.x*p_ref_lower.z/iso_distance + f_shift;
      p_ref_lower.y = ZERO;

      // create lower part of last MLC wall, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                    plane_y_min, plane_y_max,
                                    present_plane_lower, plane_x_max,
                                    plane_z_tongue, plane_z_max,
                                    p_ref_lower,
                                    leaf_material) ) == NULL )
      {
         xvmc_error("MC_mlc_elekta::init",
                    "cannot create lower part of the last MLC wall",8);
      }
      piece[n_piece] = present_region; ++n_piece;

   }                // ******** end of Y MLC ************** //

   // check number of planes and regions
   if (num_planes != n_separator)
   {
      xvmc_error("MC_mlc_elekta::init",
                 "the number of object separator planes is incorrect",8);
   }
   if (num_regions != n_piece)
   {
      xvmc_error("MC_mlc_elekta::init",
                 "the number of regions is incorrect",8);
   }

   // set bit masks and bit patterns of all regions using the reference points
   set_bits();

   return;
}

// change position of one leaf pair,
// the new leaf positions are defined at the iso-center plane
void MC_mlc_elekta::change_pair(const int  &pair_index,
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
