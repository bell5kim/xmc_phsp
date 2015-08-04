#ifndef _MC_HALFCYL_H_
#define _MC_HALFCYL_H_

/*****************************************************************************
 * MC_halfcyl.h:                                                             *
 *    class declarations and inline member functions for:                    *
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

#include "MC_plane.h"

// ****************************************************************
// class MC_halfcyl: plane with half cylinder
// ****************************************************************

class MC_halfcyl : public MC_plane
{
   public:
      // define plane with half cylinder by 4 points in 3D space and index:
      // the first 3 points (p0, p1, p2) are for the plane, the last
      // 3 points (p1, p2, p3) are for the cylinder axis, radius and to
      // determine the cylinder type (half cylinder before or behind the plane)
      MC_halfcyl(const real_3 &, const real_3 &, const real_3 &,
                 const real_3 &, const int &);

      // calculate ralationship (<0, =0 or >0) of a point to the plane
      // with half cylinder
      inline int relationship(const real_3 &);

      // calculate linear distance of a particle to the plane with half cylinder
      inline real distance(const particle_parameters &,
                           const MC_plane *);

      // set new plane parameters by 3 points in 3D space,
      // the radius and half cylinder type remain unchanged
      void set(const real_3 &, const real_3 &, const real_3 &);

   private:
      // type of the half cylinder:
      // true  if the half cylinder points to the plane normal direction
      // false if the half cylinder is behind the plane
      bool    before_plane;

      // radius of the cylinder
      real    radius;

      // direction vector of the cylinder axis
      real_3  axis_dir;

      // the point at the cylinder axis with the smallest distance to the
      // coordinate origin point, i.e. (axis_dir * axis_point) = 0
      real_3  axis_point;
};

// *******************************************
// inline member functions of class MC_halfcyl
// *******************************************

// calculate ralationship (<0, =0 or >0) of a point to the plane
// with half cylinder
// +2: the point is located before the half cylinder surface
//     but behind the plane (only for before_plane=false)
// +1: the point is located before the half cylinder surface
//     and before the plane
//  0: the point is located at the object surface
// -1: the point is located behind the half cylinder surface
//     and behind the plane
// -2: the point is located behind the half cylinder surface
//     but before the plane (only for before_plane=true)
inline int MC_halfcyl::relationship(const real_3 &point)
{
   // scalar product
   real product = point.x*normal.x + point.y*normal.y + point.z*normal.z;

   if (before_plane)
   {
      // the half cylinder is before the plane, check relationship to the plane
      if (product < dist2zero) return(-1);

      // the point is before or at the plane,
      // calculate and check distance to the cylinder axis
      real proj = point.x*axis_dir.x + point.y*axis_dir.y + point.z*axis_dir.z;
      real dist = (axis_point.x-point.x)*(axis_point.x-point.x)
                + (axis_point.y-point.y)*(axis_point.y-point.y)
                + (axis_point.z-point.z)*(axis_point.z-point.z);
      dist -= proj*proj;
#ifdef DEBUG
      if (dist < ZERO) xvmc_error("MC_halfcyl::MC_relationship",
                                  "distance squared <= 0",8);
#endif
      dist = sqrt(dist);
      if (dist <  radius) return(-2); // the point is behind the half cylinder
      if (dist == radius) return(0);  // the point is at the surface

      // the point is outside of the cylinder,
      // check again relationship to the plane
      if (product > dist2zero) return(+1);
      else                     return(0);
   }
   else
   {
      // the half cylinder is behind the plane, check relationship to the plane
      if (product > dist2zero) return(+1);

      // the point is behind or at the plane,
      // calculate and check distance to the cylinder axis
      real proj = point.x*axis_dir.x + point.y*axis_dir.y + point.z*axis_dir.z;
      real dist = (axis_point.x-point.x)*(axis_point.x-point.x)
                + (axis_point.y-point.y)*(axis_point.y-point.y)
                + (axis_point.z-point.z)*(axis_point.z-point.z);
      dist -= proj*proj;
#ifdef DEBUG
      if (dist < ZERO) xvmc_error("MC_halfcyl::MC_relationship",
                                  "distance squared <= 0",8);
#endif
      dist = sqrt(dist);
      if (dist <  radius) return(+2); // the point is before the half cylinder
      if (dist == radius) return(0);  // the point is at the surface

      // the point is outside of the cylinder,
      // check again relationship to the plane
      if (product < dist2zero) return(-1);
      else                     return(0);
   }
}

// calculate linear distance of a particle to the plane with half cylinder
inline real MC_halfcyl::distance(const particle_parameters &p,
                                 const MC_plane            *p_plane)
{
   // constants
   const real MC_INFINITY = 1.0e+10;

   // at first, we calculate the distance to the plane,
   // i.e. we need the scalar products
   real product_pos = p.pos.x*normal.x + p.pos.y*normal.y + p.pos.z*normal.z;
   real product_dir = p.dir.x*normal.x + p.dir.y*normal.y + p.dir.z*normal.z;

   // if particle direction is parallel to the plane (perpendicular to the
   // normal vector) --> set distance to infinity
   real dist2plane = MC_INFINITY;

   // if particle is not at this plane and if particle direction
   // is not parallel to this plane
   if ( (product_pos != dist2zero) && (product_dir != ZERO) )
   {
      // linear distance of particle to the plane
      real result = (dist2zero-product_pos)/product_dir;

      // we leave the distance infinity if the plane is behind the particle
      // (i.e. result <= 0)

      // set distance if the plane is before the particle
      if (result > ZERO) dist2plane = result;
   }

   // now we calculate the distance(s) to the cylinder surface, we need
   // three parameters (see paper by A.F. Bielajew, HOWFAR and HOWNEAR:
   // Geometry Modeling for Monte Carlo Particle Transport, PIRS-341)
   real pdir = p.dir.x*axis_dir.x + p.dir.y*axis_dir.y + p.dir.z*axis_dir.z;
   real proj = p.pos.x*axis_dir.x + p.pos.y*axis_dir.y + p.pos.z*axis_dir.z;
   real_3 diff; diff.x = p.pos.x - axis_point.x;
                diff.y = p.pos.y - axis_point.y;
                diff.z = p.pos.z - axis_point.z;

   // A >= 0 always,
   // A == 0 means the particle moves parallel to the cylinder axis
   real parA = ONE - pdir*pdir;

   real parB = p.dir.x*diff.x + p.dir.y*diff.y + p.dir.z*diff.z - pdir*proj;

   // C <  0: the particle is inside the cylinder
   // C == 0: the particle is at the cylinder surface
   // C >  0: the particle is outside
   real parC = diff.x*diff.x  + diff.y*diff.y  + diff.z*diff.z
             - proj*proj - radius*radius;

   // adjust C if we know that we are at this plane
   const real EPSILON  = 1.0e-03;
   if (p_plane == this)
   {
      // C must be larger than or equal to zero,
      // C must be zero if the particle is not at the plane
      if (fabs(product_pos-dist2zero) > EPSILON)
      {
         parC = ZERO;
      }
      else
      {
         product_pos = dist2zero;
         dist2plane  = MC_INFINITY;
      }
   }
   else
   {
      if (fabs(parC) < EPSILON) parC = ZERO;
   }

   // if B2AC <  0: the particle does not intersect the cylinder surface
   // if B2AC == 0: the particle touches the cylinder surface
   // in both cases we return the distance to the plane,
   // i.e. no interaction with the cylinder
   real B2AC = parB*parB - parA*parC;
   if (B2AC <= ZERO) return(dist2plane);

   // now we are sure that B2AC > 0

   // the particle is outside of the cylinder
   if (parC > ZERO)
   {
      // A < 0 is impossible, A = 0 means the particle moves parallel to
      // the cylinder axis, i.e. no intersection, return distance to plane
      if (parA <= ZERO) return(dist2plane);

      // if B >= 0 there is no intersection with the cylinder,
      // therefore we return the distance to the plane
      if (parB >= ZERO) return(dist2plane);

      // A > 0, B < 0, C > 0

      // we have two intersections with the cylinder surface
      // depending on the half cylinder type we have to find the right one
      if (before_plane)
      {
         // the half cylinder is before the plane

         // the particle is at the plane
         if (product_pos == dist2zero)
         {
            // we return the closest intersection point,
            // if the particle points before the plane
            real dist2cyl = (-parB-sqrt(B2AC))/parA;
            return( product_dir > ZERO ? dist2cyl : MC_INFINITY );
         }

         // the particle is before the plane (like the half cylinder)
         if (product_pos > dist2zero)
         {
            // i.e. we take the closest intersection point
            real dist2cyl = (-parB-sqrt(B2AC))/parA;

            // return the minimum of the distances to the plane and the
            // cylinder surface
            return( dist2plane <= dist2cyl ? dist2plane : dist2cyl );
         }

         // the particle is behind the plane

         // the closest intersection point ...
         real dist2cyl = (-parB-sqrt(B2AC))/parA;

         // ... is at the other side of the plane,
         // return distance to plane
         if (dist2cyl >= dist2plane) return(dist2plane);

         // the furthermost intersection point ...
         dist2cyl = (-parB+sqrt(B2AC))/parA;

         // ... is at this side of the plane, return distance to plane
         if (dist2cyl <= dist2plane) return(dist2plane);

         // now we can return the furthermost intersection point
         return(dist2cyl);
      }

      // the half cylinder is behind the plane

      // the particle is at the plane
      if (product_pos == dist2zero)
      {
         // we return the closest intersection point,
         // if the particle points behind the plane
         real dist2cyl = (-parB-sqrt(B2AC))/parA;
         return( product_dir < ZERO ? dist2cyl : MC_INFINITY );
      }

      // the particle is behind the plane (like the half cylinder)
      if (product_pos < dist2zero)
      {
         // i.e. we take the closest intersection point
         real dist2cyl = (-parB-sqrt(B2AC))/parA;

         // return the minimum of the distances to the plane and the
         // cylinder surface
         return( dist2plane <= dist2cyl ? dist2plane : dist2cyl );
      }

      // the particle is before the plane

      // the closest intersection point ...
      real dist2cyl = (-parB-sqrt(B2AC))/parA;

      // ... is at the other side of the plane,
      // return distance to plane
      if (dist2cyl >= dist2plane) return(dist2plane);

      // the furthermost intersection point ...
      dist2cyl = (-parB+sqrt(B2AC))/parA;

      // ... is at this side of the plane, return distance to plane
      if (dist2cyl <= dist2plane) return(dist2plane);

      // now we can return the furthermost intersection point
      return(dist2cyl);
   }

   // C <= 0

   // the particle is inside of the cylinder (not at the surface)
   if (parC < ZERO)
   {

      // A < 0 is impossible, A = 0 means the particle moves parallel to
      // the cylinder axis, i.e. no intersection, return distance to plane
      if (parA <= ZERO) return(dist2plane);

      // A > 0, C < 0: there is one positive solution
      real dist2cyl = (-parB+sqrt(B2AC))/parA;

      if (before_plane)
      {
         // the half cylinder is before the plane

         // the particle is also before the plane
         if (product_pos > dist2zero)
         {
            return( dist2cyl <= dist2plane ? dist2cyl : MC_INFINITY );
         }

         // the particle is behind the plane
         if (product_pos < dist2zero)
         {
            return( dist2cyl >= dist2plane ? dist2cyl : dist2plane );
         }

         // the particle is at the plane
         return( product_dir >= ZERO ? dist2cyl : MC_INFINITY );
      }

      // the half cylinder is behind the plane

      // the particle is also behind the plane
      if (product_pos < dist2zero)
      {
         return( dist2cyl <= dist2plane ? dist2cyl : MC_INFINITY );
      }

      // the particle is before the plane
      if (product_pos > dist2zero)
      {
         return( dist2cyl >= dist2plane ? dist2cyl : dist2plane );
      }

      // the particle is at the plane
      return( product_dir <= ZERO ? dist2cyl : MC_INFINITY );
   }

   // C = 0, the particle is at the cylinder surface

   // A < 0 is impossible, A = 0 means the particle moves parallel to
   // the cylinder axis, i.e. no intersection, return distance to plane
   if (parA <= ZERO) return(dist2plane);

   // in principle we have at least one solution (dist2cyl = 0)
   // but we want to know the distance to the next plane,
   // therfore we return distance to plane if this is the only solution
   if (parB >= ZERO) return(dist2plane);

   // A > 0, B < 0, C = 0: in this case we have a second solution
   real dist2cyl = -TWO*parB/parA;
   if (dist2cyl <= ZERO) return(dist2plane);

   if (before_plane)
   {
      // the half cylinder is before the plane

      // the particle is also before the plane
      if (product_pos > dist2zero)
      {
         return( dist2cyl <= dist2plane ? dist2cyl : MC_INFINITY );
      }

      // the particle is behind the plane
      if (product_pos < dist2zero)
      {
         return( dist2cyl >= dist2plane ? dist2cyl : dist2plane );
      }

      // the particle is at the plane
      return( product_dir >= ZERO ? dist2cyl : MC_INFINITY );
   }

   // the half cylinder is behind the plane

   // the particle is also behind the plane
   if (product_pos < dist2zero)
   {
      return( dist2cyl <= dist2plane ? dist2cyl : MC_INFINITY );
   }

   // the particle is before the plane
   if (product_pos > dist2zero)
   {
      return( dist2cyl >= dist2plane ? dist2cyl : dist2plane );
   }

   // the particle is at the plane
   return( product_dir <= ZERO ? dist2cyl : MC_INFINITY );
}

#endif /* _MC_HALFCYL_H_ */
