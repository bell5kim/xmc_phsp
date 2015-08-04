#ifndef _MC_PLANE_H_
#define _MC_PLANE_H_

/*****************************************************************************
 * MC_plane.h:                                                               *
 *    class declarations and inline member functions for:                    *
 *       MC_plane:       plane in 3D space                                   *
 *                                                                           *
 * Copyright (C) 2001    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 29.08.2001      *
 *                                                                           *
 *****************************************************************************/

#include "definitions.h"
#include "global.h"

enum axis {X, Y, Z};

// ****************************************************************
// class MC_plane: plane in 3D space
// ****************************************************************

class MC_plane
{
   public:
      // define plane perpendicular to the x, y or z axis
      MC_plane(const axis &, const real &, const int &);

      // define plane by 3 points in 3D space
      MC_plane(const real_3 &, const real_3 &, const real_3 &, const int &);

      // define plane parallel to an existing plane and one point
      MC_plane(const MC_plane *, const real_3 &, const int &);

      // get normal vector
      real_3 get_normal(void) {return(normal);}

      // get distance to origin
      real get_dist2zero(void) {return(dist2zero);}

      // calculate ralationship (<0, =0 or >0) of a point to the plane
      virtual inline int relationship(const real_3 &);

      // calculate linear distance of a particle to the plane, if we know that
      // the particle is at this plane (p_plane == this) --> return infinity
      virtual inline real distance(const particle_parameters &,
                                   const MC_plane *);

      // set plane parameters by 3 new points in 3D space
      virtual void set(const real_3 &, const real_3 &, const real_3 &);

      // plane index (index < 0 means: the index is not set)
      int     index;

   protected:
      real_3  normal;          // normal vector of the plane
      real    dist2zero;       // distance of the plane to origin
      // by definition dist2zero must be larger than or equal to ZERO,
      // the normal vector points away from origin, i.e. the origin is
      // located "behind" the plane
};

// *******************************************
// inline member functions of class MC_plane
// *******************************************

// calculate ralationship (<0, =0 or >0) of a point to the plane
// +1: the point is located before the plane, where the normal vector points to
//  0: the point is located at the plane
// -1: the point is located behind the plane
// by definition the coordinate origin (zero point) is always behind
// the plane if the zero point is not part of the plane itself
inline int MC_plane::relationship(const real_3 &point)
{
   // scalar product
   real product = point.x*normal.x + point.y*normal.y + point.z*normal.z;

   // compare scalar product and plane distance to zero
   if      (product > dist2zero) return(+1);
   else if (product < dist2zero) return(-1);
   else                          return(0);
}

// calculate linear distance of a particle to the plane, if we know that
// the particle is at this plane (p_plane == this) --> return infinity
inline real MC_plane::distance(const particle_parameters &p,
                               const MC_plane            *p_plane)
{
   // constants
   const real MC_INFINITY = 1.0e+10;

   // if the particle is located directly at the plane, the distance is
   // either zero or infinity -- we return infinity
   if (p_plane == this) return(MC_INFINITY);

   // scalar products
   real product_pos = p.pos.x*normal.x + p.pos.y*normal.y + p.pos.z*normal.z;
   real product_dir = p.dir.x*normal.x + p.dir.y*normal.y + p.dir.z*normal.z;

   // if particle direction is parallel to the plane (perpendicular to the
   // normal vector) --> set distance to infinity
   if (product_dir == ZERO) return(MC_INFINITY);

   // linear distance of particle to the plane
   real result = (dist2zero-product_pos)/product_dir;

   // we set the distance to infinity if the plane is behind the particle
   if (result <= ZERO) return(MC_INFINITY);
   else                return(result);
}

#endif /* _MC_PLANE_ */
