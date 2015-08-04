/*****************************************************************************
 * MC_halfcyl.cpp:                                                           *
 *    class member functions for:                                            *
 *       MC_halfcyl:     plane with half cylinder                            *
 *                                                                           *
 * Copyright (C) 2001    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 14.11.2001      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include "MC_halfcyl.h"

// ****************************************
// member functions of class MC_halfcyl
// ****************************************

// define plane with half cylinder by 4 points in 3D space and index:
// the first 3 points (p0, p1, p2) are for the plane, the last
// 3 points (p1, p2, p3) are for the cylinder axis, radius and to
// determine the cylinder type (half cylinder before or behind the plane)
MC_halfcyl::MC_halfcyl(const real_3 &p0, const real_3 &p1, const real_3 &p2,
                       const real_3 &p3, const int &ind)
          : MC_plane(p0,p1,p2,ind)
{
   // set plane parameters
   set(p0,p1,p2);

   // to determine the half cylinder type we check the position
   // of point p3 relative to the plane
   int p3_position = MC_plane::relationship(p3);
   if (p3_position > 0)
   {
      before_plane = true;
   }
   else
   {
      if (p3_position < 0)
      {
         before_plane = false;
      }
      else
      {
         xvmc_error("MC_halfcyl::MC_halfcyl",
                    "point 3 must not be part of the plane",8);
      }
   }

   // projection of point p3 to the cylinder axis
   real proj = p3.x*axis_dir.x + p3.y*axis_dir.y + p3.z*axis_dir.z;

   // distance squared of point p3 to the axis point
   real dist2 = (axis_point.x-p3.x)*(axis_point.x-p3.x)
              + (axis_point.y-p3.y)*(axis_point.y-p3.y)
              + (axis_point.z-p3.z)*(axis_point.z-p3.z);

   // calculate radius using the Theorem of Pythagoras
   radius = dist2 - proj*proj;
   if (radius > ZERO)
   {
      radius = sqrt(radius);
   }
   else
   {
      xvmc_error("MC_halfcyl::MC_halfcyl","radius squared <= 0",8);
   }
}

// set new plane parameters by 3 points in 3D space,
// the radius and half cylinder type remain unchanged
void MC_halfcyl::set(const real_3 &p0, const real_3 &p1, const real_3 &p2)
{
   // set plane parameters
   MC_plane::set(p0,p1,p2);

   // we use the last two points to define the cylinder axis
   axis_dir.x = p2.x - p1.x;
   axis_dir.y = p2.y - p1.y;
   axis_dir.z = p2.z - p1.z;

   // normalize axis direction vector
   real norm =
      sqrt(axis_dir.x*axis_dir.x+axis_dir.y*axis_dir.y+axis_dir.z*axis_dir.z);

   // the normalization factor must be larger than zero
   const real EPSILON  = 1.0e-05;
   if (norm > EPSILON)
   {
      axis_dir.x /= norm;
      axis_dir.y /= norm;
      axis_dir.z /= norm;
   }
   else
   {
      xvmc_error("MC_halfcyl::set",
                 "the 3 points are improperly defined",8);
   }

   // determine the axis point with the smallest distance to zero
   real product = p1.x*axis_dir.x + p1.y*axis_dir.y + p1.z*axis_dir.z;
   axis_point.x = p1.x - product*axis_dir.x;
   axis_point.y = p1.y - product*axis_dir.y;
   axis_point.z = p1.z - product*axis_dir.z;
}
