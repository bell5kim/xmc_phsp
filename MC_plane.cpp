/*****************************************************************************
 * MC_plane.cpp:                                                             *
 *    class member functions for:                                            *
 *       MC_plane:       plane in 3D space                                   *
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

#include "MC_plane.h"

// ****************************************
// member functions of class MC_plane
// ****************************************

// define plane perpendicular to the x, y or z axis
MC_plane::MC_plane(const axis &xyz, const real &d, const int &ind)
{
   // set index
   index = ind;

   // initialize normal vector by zeros
   normal.x = ZERO;
   normal.y = ZERO;
   normal.z = ZERO;

   // if d < 0 the plane is located at the negative axis part
   dist2zero = fabs(d);
   if (d < ZERO)
   {
      if (xyz == X) normal.x = -ONE;
      if (xyz == Y) normal.y = -ONE;
      if (xyz == Z) normal.z = -ONE;
   }
   else
   {
      if (xyz == X) normal.x = ONE;
      if (xyz == Y) normal.y = ONE;
      if (xyz == Z) normal.z = ONE;
   }
}

// define plane by 3 points in 3D space
MC_plane::MC_plane(const real_3 &p0, const real_3 &p1, const real_3 &p2,
                   const int &ind)
{
   // set index
   index = ind;

   // set plane parameters using the 3 points
   set(p0,p1,p2);
}

// define plane parallel to an existing plane and one point
MC_plane::MC_plane(const MC_plane *old_plane, const real_3 &point,
                   const int &ind)
{
   // set index
   index = ind;

   // the normal vectors are identical
   normal = old_plane->normal;

   // distance to the origin is given by the scalar product
   dist2zero = point.x*normal.x + point.y*normal.y + point.z*normal.z;

   // by definition dist2zero must be larger than or equal to zero
   if (dist2zero < ZERO)
   {
      normal.x  *= -ONE;
      normal.y  *= -ONE;
      normal.z  *= -ONE;
      dist2zero *= -ONE;
   }
}

// set plane parameters by 3 new points in 3D space
void MC_plane::set(const real_3 &p0, const real_3 &p1, const real_3 &p2)
{
   // calculate vector perpendicular to the plane,
   // after normalization this will become the normal vector
   normal.x = (p1.y-p0.y)*(p2.z-p0.z) - (p2.y-p0.y)*(p1.z-p0.z);
   normal.y = (p1.z-p0.z)*(p2.x-p0.x) - (p2.z-p0.z)*(p1.x-p0.x);
   normal.z = (p1.x-p0.x)*(p2.y-p0.y) - (p2.x-p0.x)*(p1.y-p0.y);

   // distance to the origin point before normalization
   dist2zero = normal.x*p0.x + normal.y*p0.y + normal.z*p0.z;

   // calculate normalization factor
   real norm = sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);

   // the normalization factor must be larger than zero
   const real EPSILON  = 1.0e-05;
   if (norm > EPSILON)
   {
      normal.x  /= norm;
      normal.y  /= norm;
      normal.z  /= norm;
      dist2zero /= norm;
   }
   else
   {
      xvmc_error("MC_plane::set",
                 "the 3 points are improperly defined",8);
   }

   // by definition dist2zero must be larger than or equal to zero
   if (dist2zero < ZERO)
   {
      normal.x  *= -ONE;
      normal.y  *= -ONE;
      normal.z  *= -ONE;
      dist2zero *= -ONE;
   }
}
