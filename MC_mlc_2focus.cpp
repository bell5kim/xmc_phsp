/*****************************************************************************
 * MC_mlc_2focus.cpp:                                                        *
 *    class member functions for:                                            *
 *       MC_mlc_2focus:  double focussing MLC                                *
 *                                                                           *
 * Copyright (C) 2001    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 19.11.2001      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <new>
using namespace std;

#include "MC_mlc_2focus.h"

// ****************************************
// member functions of class MC_mlc_2focus
// ****************************************

// initialize MLC
void MC_mlc_2focus::init(void)
{
   // counter for the object separator planes and object pieces
   unsigned n_separator = 0, n_piece = 0;

   // create the 6 outer planes
   plane_x_min = NULL;
   if ( (plane_x_min = new (nothrow) MC_plane(X,x_min,-1)) == NULL )
   {
      xvmc_error("MC_mlc_2focus::init",
                 "cannot create minimum x plane",8);
   }
   separator[n_separator] = plane_x_min; ++n_separator;

   plane_x_max = NULL;
   if ( (plane_x_max = new (nothrow) MC_plane(X,x_max,-1)) == NULL )
   {
      xvmc_error("MC_mlc_2focus::init",
                 "cannot create maximum x plane",8);
   }
   separator[n_separator] = plane_x_max; ++n_separator;

   plane_y_min = NULL;
   if ( (plane_y_min = new (nothrow) MC_plane(Y,y_min,-1)) == NULL )
   {
      xvmc_error("MC_mlc_2focus::init",
                 "cannot create minimum y plane",8);
   }
   separator[n_separator] = plane_y_min; ++n_separator;

   plane_y_max = NULL;
   if ( (plane_y_max = new (nothrow) MC_plane(Y,y_max,-1)) == NULL )
   {
      xvmc_error("MC_mlc_2focus::init",
                 "cannot create maximum y plane",8);
   }
   separator[n_separator] = plane_y_max; ++n_separator;

   // typically the starting plane of the particle transport
   plane_z_min = NULL;
   if ( (plane_z_min = new (nothrow) MC_plane(Z,z_min,-1)) == NULL )
   {
      xvmc_error("MC_mlc_2focus::init",
                 "cannot create minimum z plane",8);
   }
   separator[n_separator] = plane_z_min; ++n_separator;
   starting_plane = plane_z_min;

   // typically the final plane of the particle transport
   plane_z_max = NULL;
   if ( (plane_z_max = new (nothrow) MC_plane(Z,z_max,-1)) == NULL )
   {
      xvmc_error("MC_mlc_2focus::init",
                 "cannot create maximum z plane",8);
   }
   separator[n_separator] = plane_z_max; ++n_separator;
   final_plane = plane_z_max;

   // define pointers to the previous and present wall planes
   MC_plane  *previous_plane  = NULL;
   MC_plane  *present_plane   = NULL;

   // positions of the previous and present wall planes
   real previous_pos = ZERO;
   real present_pos  = ZERO;

   // define pointers to the left and right leaf edge planes
   MC_plane  *plane_left      = NULL;
   MC_plane  *plane_right     = NULL;

   // pointer to the present region
   MC_region *present_region  = NULL;

   // the focus points
   real_3 p_focus_wall,p_focus_edge;
   p_focus_wall.x = ZERO; p_focus_wall.y = ZERO; p_focus_wall.z = z_focus_wall;
   p_focus_edge.x = ZERO; p_focus_edge.y = ZERO; p_focus_edge.z = z_focus_edge;

   // we need 2 further points for the leaf walls
   real_3 p1_wall,p2_wall;
   p1_wall.x = ZERO;           p2_wall.x = ZERO;
   p1_wall.y = ZERO;           p2_wall.y = ZERO;
   p1_wall.z = iso_distance;   p2_wall.z = iso_distance;

   // four further points for the left and right leaf edges
   real_3 p1_left,p2_left,p1_right,p2_right;
   p1_left.x  = ZERO;          p2_left.x  = ZERO;
   p1_left.y  = ZERO;          p2_left.y  = ZERO;
   p1_left.z  = iso_distance;  p2_left.z  = iso_distance;
   p1_right.x = ZERO;          p2_right.x = ZERO;
   p1_right.y = ZERO;          p2_right.y = ZERO;
   p1_right.z = iso_distance;  p2_right.z = iso_distance;

   // define reference point to identify the regions
   real_3 p_ref;
   p_ref.x = ZERO;
   p_ref.y = ZERO;
   p_ref.z = (z_min+z_max)/TWO; // the reference plane

   // create regions
   if (xytype == X)
   {                // ************ X MLC ***************** //

      // set wall points
      present_pos = nominal_mlc->get_lperp(0);
      p1_wall.x = -10.0;          p2_wall.x = 10.0;
      p1_wall.y = present_pos;    p2_wall.y = present_pos;
      p1_wall.z = iso_distance;   p2_wall.z = iso_distance;

      // create first leaf wall plane
      present_plane = NULL;
      if ( (present_plane = new (nothrow)
         MC_plane(p_focus_wall,p1_wall,p2_wall,-1)) == NULL )
      {
         xvmc_error("MC_mlc_2focus::init",
                    "cannot create first leaf wall plane",8);
      }
      separator[n_separator] = present_plane; ++n_separator;

      // set reference point for the upper MLC wall
      p_ref.x = ZERO;
      p_ref.y = (y_min + present_pos*p_ref.z/iso_distance)/TWO;

      // create upper MLC wall, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                              plane_x_min, plane_x_max,
                                              plane_y_min, present_plane,
                                              plane_z_min, plane_z_max,
                                              p_ref,
                                              leaf_material) ) == NULL )
      {
         xvmc_error("MC_mlc_2focus::init",
                    "cannot create upper MLC wall",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // start loop for all leaf pairs
      for (int i_pair=0; i_pair<num_pairs; ++i_pair)
      {
         // set next wall points
         previous_pos = present_pos;
         present_pos  = nominal_mlc->get_lperp(i_pair+1);
         p1_wall.x = -10.0;          p2_wall.x = 10.0;
         p1_wall.y = present_pos;    p2_wall.y = present_pos;
         p1_wall.z = iso_distance;   p2_wall.z = iso_distance;

         // create next leaf wall plane
         previous_plane = present_plane;
         present_plane  = NULL;
         if ( (present_plane = new (nothrow)
            MC_plane(p_focus_wall,p1_wall,p2_wall,-1)) == NULL )
         {
            xvmc_error("MC_mlc_2focus::init",
                       "cannot create next leaf wall plane",8);
         }
         separator[n_separator] = present_plane; ++n_separator;

         // left leaf position at the iso-center
         real open_left = nominal_mlc->get_left(i_pair);
         p1_left.x = open_left; p1_left.y = -10.0;
         p2_left.x = open_left; p2_left.y =  10.0;

         // create left leaf edge plane
         plane_left = NULL;
         if ( (plane_left = new (nothrow)
            MC_plane(p_focus_edge,p1_left,p2_left,-1)) == NULL )
         {
            xvmc_error("MC_mlc_2focus::init",
                       "cannot create left leaf edge plane",8);
         }
         separator[n_separator]   = plane_left; ++n_separator;
         left_edge_planes[i_pair] = plane_left;

         // right leaf position at the iso-center
         real open_right = nominal_mlc->get_right(i_pair);
         p1_right.x = open_right; p1_right.y = -10.0;
         p2_right.x = open_right; p2_right.y =  10.0;

         // create right leaf edge plane
         plane_right = NULL;
         if ( (plane_right = new (nothrow)
            MC_plane(p_focus_edge,p1_right,p2_right,-1)) == NULL )
         {
            xvmc_error("MC_mlc_2focus::init",
                       "cannot create right leaf edge plane",8);
         }
         separator[n_separator]    = plane_right; ++n_separator;
         right_edge_planes[i_pair] = plane_right;

         // the left and right leaf positions in the reference plane
         real r_left  = open_left*p_ref.z/iso_distance;
         real r_right = open_right*p_ref.z/iso_distance;

         // the upper and lower leaf walls in the reference plane
         real r_upper = previous_pos*p_ref.z/iso_distance;
         real r_lower = present_pos*p_ref.z/iso_distance;

         // set reference point for the left leaf
         p_ref.x = (x_min+r_left)/TWO;
         p_ref.y = (r_upper+r_lower)/TWO;

         // create left leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                                 plane_x_min,    plane_left,
                                                 previous_plane, present_plane,
                                                 plane_z_min,    plane_z_max,
                                                 p_ref,
                                                 leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_2focus::init",
                       "cannot create left leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the right leaf
         p_ref.x = (r_right+x_max)/TWO;
         p_ref.y = (r_upper+r_lower)/TWO;

         // create right leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                                 plane_right,    plane_x_max,
                                                 previous_plane, present_plane,
                                                 plane_z_min,    plane_z_max,
                                                 p_ref,
                                                 leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_2focus::init",
                       "cannot create right leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the intermediate leaf region
         p_ref.x = (r_left+r_right)/TWO;
         p_ref.y = (r_upper+r_lower)/TWO;

         // create intermediate leaf region, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                                 plane_left,     plane_right,
                                                 previous_plane, present_plane,
                                                 plane_z_min,    plane_z_max,
                                                 p_ref,
                                                 open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_2focus::init",
                       "cannot create intermediate leaf region",8);
         }
         piece[n_piece] = present_region; ++n_piece;

      } // end of leaf pair loop

      // set reference point for the lower MLC wall
      p_ref.x = ZERO;
      p_ref.y = (y_max + present_pos*p_ref.z/iso_distance)/TWO;

      // create lower MLC wall, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                              plane_x_min,   plane_x_max,
                                              present_plane, plane_y_max,
                                              plane_z_min,   plane_z_max,
                                              p_ref,
                                              leaf_material) ) == NULL )
      {
         xvmc_error("MC_mlc_2focus::init",
                    "cannot create lower MLC wall",8);
      }
      piece[n_piece] = present_region; ++n_piece;

   }                // ******** end of X MLC ************** //
   else
   {                // ************ Y MLC ***************** //

      // set wall points
      present_pos = nominal_mlc->get_lperp(0);
      p1_wall.x = present_pos;    p2_wall.x = present_pos;
      p1_wall.y = -10.0;          p2_wall.y = 10.0;
      p1_wall.z = iso_distance;   p2_wall.z = iso_distance;

      // create first leaf wall plane
      present_plane = NULL;
      if ( (present_plane = new (nothrow)
         MC_plane(p_focus_wall,p1_wall,p2_wall,-1)) == NULL )
      {
         xvmc_error("MC_mlc_2focus::init",
                    "cannot create first leaf wall plane",8);
      }
      separator[n_separator] = present_plane; ++n_separator;

      // set reference point for the upper MLC wall
      p_ref.x = (x_min + present_pos*p_ref.z/iso_distance)/TWO;
      p_ref.y = ZERO;

      // create upper MLC wall, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                              plane_y_min, plane_y_max,
                                              plane_x_min, present_plane,
                                              plane_z_min, plane_z_max,
                                              p_ref,
                                              leaf_material) ) == NULL )
      {
         xvmc_error("MC_mlc_2focus::init",
                    "cannot create upper MLC wall",8);
      }
      piece[n_piece] = present_region; ++n_piece;

      // start loop for all leaf pairs
      for (int i_pair=0; i_pair<num_pairs; ++i_pair)
      {
         // set next wall points
         previous_pos = present_pos;
         present_pos  = nominal_mlc->get_lperp(i_pair+1);
         p1_wall.x = present_pos;    p2_wall.x = present_pos;
         p1_wall.y = -10.0;          p2_wall.y = 10.0;
         p1_wall.z = iso_distance;   p2_wall.z = iso_distance;

         // create next leaf wall plane
         previous_plane = present_plane;
         present_plane  = NULL;
         if ( (present_plane = new (nothrow)
            MC_plane(p_focus_wall,p1_wall,p2_wall,-1)) == NULL )
         {
            xvmc_error("MC_mlc_2focus::init",
                       "cannot create next leaf wall plane",8);
         }
         separator[n_separator] = present_plane; ++n_separator;

         // left leaf position at the iso-center
         real open_left = nominal_mlc->get_left(i_pair);
         p1_left.x = -10.0; p1_left.y = open_left;
         p2_left.x =  10.0; p2_left.y = open_left;

         // create left leaf edge plane
         plane_left = NULL;
         if ( (plane_left = new (nothrow)
            MC_plane(p_focus_edge,p1_left,p2_left,-1)) == NULL )
         {
            xvmc_error("MC_mlc_2focus::init",
                       "cannot create left leaf edge plane",8);
         }
         separator[n_separator]   = plane_left; ++n_separator;
         left_edge_planes[i_pair] = plane_left;

         // right leaf position at the iso-center
         real open_right = nominal_mlc->get_right(i_pair);
         p1_right.x = -10.0; p1_right.y = open_right;
         p2_right.x =  10.0; p2_right.y = open_right;

         // create right leaf edge plane
         plane_right = NULL;
         if ( (plane_right = new (nothrow)
            MC_plane(p_focus_edge,p1_right,p2_right,-1)) == NULL )
         {
            xvmc_error("MC_mlc_2focus::init",
                       "cannot create right leaf edge plane",8);
         }
         separator[n_separator]    = plane_right; ++n_separator;
         right_edge_planes[i_pair] = plane_right;

         // the left and right leaf positions in the reference plane
         real r_left  = open_left*p_ref.z/iso_distance;
         real r_right = open_right*p_ref.z/iso_distance;

         // the upper and lower leaf walls in the reference plane
         real r_upper = previous_pos*p_ref.z/iso_distance;
         real r_lower = present_pos*p_ref.z/iso_distance;

         // set reference point for the left leaf
         p_ref.x = (r_upper+r_lower)/TWO;
         p_ref.y = (y_min+r_left)/TWO;

         // create left leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                                 plane_y_min,    plane_left,
                                                 previous_plane, present_plane,
                                                 plane_z_min,    plane_z_max,
                                                 p_ref,
                                                 leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_2focus::init",
                       "cannot create left leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the right leaf
         p_ref.x = (r_upper+r_lower)/TWO;
         p_ref.y = (r_right+y_max)/TWO;

         // create right leaf, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                                 plane_right,    plane_y_max,
                                                 previous_plane, present_plane,
                                                 plane_z_min,    plane_z_max,
                                                 p_ref,
                                                 leaf_material) ) == NULL )
         {
            xvmc_error("MC_mlc_2focus::init",
                       "cannot create right leaf",8);
         }
         piece[n_piece] = present_region; ++n_piece;

         // set reference point for the intermediate leaf region
         p_ref.x = (r_upper+r_lower)/TWO;
         p_ref.y = (r_left+r_right)/TWO;

         // create intermediate leaf region, a region bounded by 6 planes
         present_region = NULL;
         if ( (present_region = new (nothrow) MC_volume_6p(
                                                 plane_left,     plane_right,
                                                 previous_plane, present_plane,
                                                 plane_z_min,    plane_z_max,
                                                 p_ref,
                                                 open_material) ) == NULL )
         {
            xvmc_error("MC_mlc_2focus::init",
                       "cannot create intermediate leaf region",8);
         }
         piece[n_piece] = present_region; ++n_piece;

      } // end of leaf pair loop

      // set reference point for the lower MLC wall
      p_ref.x = (x_max + present_pos*p_ref.z/iso_distance)/TWO;
      p_ref.y = ZERO;

      // create lower MLC wall, a region bounded by 6 planes
      present_region = NULL;
      if ( (present_region = new (nothrow) MC_volume_6p(
                                              plane_y_min,   plane_y_max,
                                              present_plane, plane_x_max,
                                              plane_z_min,   plane_z_max,
                                              p_ref,
                                              leaf_material) ) == NULL )
      {
         xvmc_error("MC_mlc_2focus::init",
                    "cannot create lower MLC wall",8);
      }
      piece[n_piece] = present_region; ++n_piece;

   }                // ******** end of Y MLC ************** //

   // check number of planes and regions
   if (num_planes != n_separator)
   {
      xvmc_error("MC_mlc_2focus::init",
                 "the number of object separator planes is incorrect",8);
   }
   if (num_regions != n_piece)
   {
      xvmc_error("MC_mlc_2focus::init",
                 "the number of regions is incorrect",8);
   }

   // set bit masks and bit patterns of all regions using the reference points
   set_bits();

   return;
}

// change position of one leaf pair,
// the new leaf positions are defined at the iso-center plane
void MC_mlc_2focus::change_pair(const int  &pair_index,
                                const real &new_left,   const real &new_right)
{
   // perform checks and set nominal MLC data
   MC_mlc::change_pair(pair_index, new_left, new_right);

   // the leaf edge focus point
   real_3 p_focus;
   p_focus.x = ZERO; p_focus.y = ZERO; p_focus.z = z_focus_edge;

   // four further points for the left and right leaf edges
   real_3 p1_left,p2_left,p1_right,p2_right;
   if (xytype == X)
   {
      p1_left.x  = new_left;      p2_left.x  = new_left;
      p1_left.y  = -10.0;         p2_left.y  = 10.0;
      p1_left.z  = iso_distance;  p2_left.z  = iso_distance;
      p1_right.x = new_right;     p2_right.x = new_right;
      p1_right.y = -10.0;         p2_right.y = 10.0;
      p1_right.z = iso_distance;  p2_right.z = iso_distance;
   }
   else
   {
      p1_left.x  = -10.0;         p2_left.x  = 10.0;
      p1_left.y  = new_left;      p2_left.y  = new_left;
      p1_left.z  = iso_distance;  p2_left.z  = iso_distance;
      p1_right.x = -10.0;         p2_right.x = 10.0;
      p1_right.y = new_right;     p2_right.y = new_right;
      p1_right.z = iso_distance;  p2_right.z = iso_distance;
   }

   // now move the left and right leaf edge planes to the new positions
   left_edge_planes[pair_index]->set(p_focus,p1_left,p2_left);
   right_edge_planes[pair_index]->set(p_focus,p1_right,p2_right);

   return;
}
