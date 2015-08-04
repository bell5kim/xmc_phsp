#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

/*****************************************************************************
 * geometry.h:                                                               *
 *    class declarations and inline member functions for:                    *
 *       MC_plane:      plane in 3D space                                    *
 *       MC_halfcyl:    plane with half cylinder                             *
 *       MC_region:     geometrical region defined by an array of planes     *
 *       MC_slab:       slab defined by two parallel planes                  *
 *       MC_volume_6p:  region (volume) defined by 6 planes                  *
 *       MC_object:     geometrical object defined by regions and planes     *
 *       MC_doubleslab: two slabs defined by 3 parallel planes               *
 *       MC_jaws_focus: a pair of focussing jaws                             *
 *       MC_mlc:        multi-leaf collimator (MLC) base class               *
 *                                                                           *
 * Copyright (C) 2001    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 29.08.2001      *
 *    class MC_halfcyl                                    MF 14.11.2001      *
 *    class MC_mlc                                        MF 03.12.2001      *
 *                                                                           *
 *****************************************************************************/

#define SECONDARY_MODIFIER_ELECTRONS
//#define CPU_WATCH

#include <math.h>
#include "definitions.h"
#include "global.h"

#ifdef CPU_WATCH
#include "cpu_watch.h"
#endif

#include "bitset.h"
#include "xvmc_util.h"
#include "rndm.h"
#include "interactions.h"
#include "contour.h"

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

// *******************************************************************
// class MC_region: geometrical region defined by an array of planes
// *******************************************************************

class MC_region
{
   public:
      // allocate new region with a given number of surface planes,
      // a material specified by the ESTAR (NIST, ICRU) electron
      // stopping power file, the NIST photon cross section file,
      // the differential Compton and pair cross section files
      MC_region(unsigned, char *, char *, char *, char *);

      // delete region
      ~MC_region(void);

      // get number of planes defining the region
      unsigned  get_num_planes(void) { return(num_planes); }

      // get the reference point
      real_3 get_p_ref(void) { return(p_ref); }

      // calculate ralationship (<0, =0 or >0) of a point to the region
      // +1: the point is outside the region
      //  0: the point is located at one of the surface planes
      // -1: the point is inside
      // inline int relationship(const real_3 &);

      // calculate distance of a particle to the region surface by taking
      // into account every surface plane and the present particle plane given
      // by the plane pointer, this pointer also returns the target plane
      inline real distance(const particle_parameters &, MC_plane *&);

      // array of pointers to plane defining the surface of the region
      MC_plane    **surface;

      // Compton cross section
      compton_XS_total *tot_comp;
      compton_XS_diff  *compton;

      // pair cross section
      pair_XS_total *tot_pair;
      pair_XS_diff  *pair;

      // total photo cross section
      photo_XS_total *tot_phot;

      // electron transport data (stopping powers and ranges)
      electron_transport_data *e_data;

   protected:
      // number of planes
      unsigned   num_planes;

      // to identify the inside of the region we need a reference point within
      // the region to calculate the reference point to surface plane
      // relationships
      real_3   p_ref;
};

// *******************************************
// inline member functions of class MC_region
// *******************************************

// calculate distance of a particle to the region surface by taking
// into account every surface plane and the present particle plane given
// by the plane pointer, this pointer also returns the target plane
inline real MC_region::distance(const particle_parameters &p,
                                MC_plane *&present_plane)
{
   // constants
   const real MC_INFINITY = 1.0e+10;

   // if the number of planes is 0, i.e. if the region is an exterior region,
   // return infinity distance and plane pointer to NULL
   if (num_planes == 0)
   {
      present_plane = NULL;
      return(MC_INFINITY);
   }

   // the number of planes is larger than 0
   real       dist_min  = MC_INFINITY;
   MC_plane  *new_plane = NULL;
   for (register unsigned i=0; i<num_planes; ++i)
   {
      // old version: if (surface[i] != present_plane)
      // old version: {
      // calculate linear distance of the particle to plane i
      real distance = surface[i]->distance(p,present_plane);
      if (distance < dist_min)
      {
         dist_min  = distance;
         new_plane = surface[i];
      }
      // old version: }
   }
   present_plane = new_plane;
   return(dist_min);
}

// ****************************************************************
// class MC_slab: slab defined by two parallel planes
// ****************************************************************

class MC_slab : public MC_region
{
   private:
      // initialize slab by two existing parallel planes
      void init(MC_plane *, MC_plane *);

   public:
      // define slab by two existing parallel planes
      MC_slab(MC_plane *plane0, MC_plane *plane1,
              char *material_estar, char *material_nist,
              char *compton_file,   char *pair_file)
         : MC_region(2,material_estar,material_nist,compton_file,pair_file)
              { init(plane0,plane1); }

      // define slab perpendicular to the x, y or z axis
      MC_slab(const axis &, const real &, const real &,
              char *, char *, char *, char *);

      // define slab by four points in 3D space
      // the first three points define one plane, the fourth point is to
      // define the second plane parallel to plane one
      MC_slab(const real_3 &, const real_3 &,
              const real_3 &, const real_3 &,
              char *, char *, char *, char *);
};

// ****************************************************************
// class MC_volume_6p: region (volume) defined by 6 planes
// ****************************************************************

class MC_volume_6p : public MC_region
{
   private:
      // initialize volume by 6 planes
      void init(MC_plane *, MC_plane *, MC_plane *,
                MC_plane *, MC_plane *, MC_plane *);

   public:
      // define volume by 6 planes
      MC_volume_6p(MC_plane *plane0, MC_plane *plane1,  MC_plane *plane2,
                   MC_plane *plane3, MC_plane *plane4,  MC_plane *plane5,
                   const real_3 &p_inp,
                   char *material_estar, char *material_nist,
                   char *compton_file,   char *pair_file)
         : MC_region(6,material_estar,material_nist,compton_file,pair_file) {
              p_ref.x = p_inp.x; p_ref.y = p_inp.y; p_ref.z = p_inp.z;
              init(plane0,plane1,plane2,plane3,plane4,plane5); }
};

// *******************************************************************
// class MC_object: geometrical object defined by regions and planes
// *******************************************************************

class MC_object
{
   public:
      MC_object(unsigned, unsigned);  // allocate new object
      virtual ~MC_object(void);       // delete object

      // get number of planes defining the object
      unsigned  get_num_planes(void) { return(num_planes); }

      // get bit mask of one region
      MCBitSet  get_bit_mask(const unsigned i_region)
                   { return(*bit_mask[i_region]); }

      // get bit pattern of one region
      MCBitSet  get_bit_pattern(const unsigned i_region)
                   { return(*bit_pattern[i_region]); }

      // get number of regions defining the object
      unsigned  get_num_regions(void) { return(num_regions); }

      // get particle starting plane of the object
      MC_plane *get_starting_plane(void) { return(starting_plane); }

      // get final particle transport plane of the object
      MC_plane *get_final_plane(void) { return(final_plane); }

      // estimate the region index of a point, default: -1 == no estimation
      virtual int estimate_region_index(const real_3 &p0) { return(-1); }

      // calculate particle bit pattern by determining the plane relationships,
      // only the bit associated with the specified plane will be set explicitly
      inline MCBitSet calc_bit_pattern(const particle_parameters &,
                                       const MC_plane *, const bool &);

      // determine the region index of the particle by comparing the
      // the masked particle bit pattern with the region bit patterns,
      // return -1 if the particle is outside the object,
      // the integer input variable can be used for an index estimation
      inline int get_region_index(const MCBitSet &, const int &);

      // transport photon through the object, the plane pointer points to
      // the starting plane of the photon (it must not be the NULL pointer),
      // if the photon survives (or one daughter particle) -> return true,
      // if no particle survives  -> return false,
      // the particle_parameters struct returns at most one particle
      inline bool primary_photon(particle_parameters &,
                                 MC_plane *, ranmar &);

      // transport electron through the object, the plane pointer points to
      // the starting plane of the electron (it must not be the NULL pointer),
      // if the electron survives (and/or daughter particles) -> return true,
      // if no particle survives  -> return false,
      // this is just a Continuous Slowing Down Approximation (CSDA),
      // the random number generator is reserved for later use
      inline bool primary_electron(particle_parameters &,
                                   MC_plane *, ranmar &);

      // transport secondary electron through the object,
      // the bit pattern and the region index of the electron must be known,
      // if the electron survives (and/or daughter particles) -> return true,
      // if no particle survives  -> return false,
      // this is just a Continuous Slowing Down Approximation (CSDA),
      // the random number generator is reserved for later use
      inline bool secondary_electron(particle_parameters &,
                                     MCBitSet &, int &, ranmar &);

      // array of pointers to plane defining the object
      MC_plane    **separator;

      // array of pointers to region
      MC_region   **piece;

   protected:
      // number of planes
      unsigned      num_planes;

      // number of regions, including the outer region with index 0
      unsigned      num_regions;

      // pointer to the starting plane
      // (generally the particle transport begins at this plane)
      MC_plane     *starting_plane;

      // pointer to the final plane
      // (the transport will be stopped if the particle crosses this plane)
      MC_plane     *final_plane;

      // set plane indices for all planes of the object,
      // set bit masks of all regions by checking the surface plane indices and
      // set bit patterns of all regions by checking the relationships of the
      // region reference points to the corresponding surface planes
      void set_bits(void);

   private:
      // for each region of this object we define a bit mask,
      // set bit if the corresponding plane is a surface plane of the region,
      // clear bit if the corresponding plane has nothing to do with the plane
      MCBitSet **bit_mask;

      // for each region of this object we define a bit pattern
      // providing the region to plane relationships,
      // set bit to 1 (or true) if the corresponding plane is a surface plane
      // of the region (see bit_mask) and the region is behind this plane,
      // clear bit (set to 0 or to false) if the corresponding plane is a
      // surface plane of the region (see bit_mask) and the region is before
      // this plane, clear also bit if the corresponding plane has nothing
      // to do with the region (see bit_mask)
      MCBitSet **bit_pattern;
};

// *******************************************
// inline member functions of class MC_object
// *******************************************

// calculate particle bit pattern by determining the plane relationships,
// only the bit associated with the specified plane will be set explicitly
inline MCBitSet MC_object::calc_bit_pattern(const particle_parameters &p,
                                            const MC_plane *p_plane,
                                            const bool &p_plane_bit)
{
#ifdef CPU_WATCH
   extern cpu_watch t_cpu_pbp; t_cpu_pbp.start();
#endif

   MCBitSet p_bit_pattern(num_planes);

   for (register unsigned i_plane=0; i_plane<num_planes; ++i_plane)
   {
      if (separator[i_plane] == p_plane)
      {
         // set bit explicitly for plane p_plane
         if (p_plane_bit) p_bit_pattern.set(i_plane);
         else             p_bit_pattern.clear(i_plane);
      }
      else
      {
         // calculate bits for all planes with the exception of plane p_plane
         if (separator[i_plane]->relationship(p.pos) < 0)
         {
            // the particle is behind the plane, set bit
            p_bit_pattern.set(i_plane);
         }
         else
         {
            // the particle is either before or at the plane
            if (separator[i_plane]->relationship(p.pos) > 0)
            {
               // the particle is before the plane, clear bit
               p_bit_pattern.clear(i_plane);
            }
            else
            {
               // the particle is at the plane, therefore we take the particle
               // direction into account, compare with the plane normal vector
               real_3 normal = separator[i_plane]->get_normal();
            
               // scalar product of particle direction and plane normal vectors
               real product =
                  p.dir.x*normal.x + p.dir.y*normal.y + p.dir.z*normal.z;

               if (product > ZERO)
               {
                  // the particle points to the region before the plane,
                  // clear bit
                  p_bit_pattern.clear(i_plane);
               }
               else
               {
                  // the particle points to the region behind the plane,
                  // set bit,
                  // or it will stay at the plane, set also bit because
                  // in this case we define the particle to be behind the plane
                  p_bit_pattern.set(i_plane);
               }
            }
         }
      }
   }


#ifdef CPU_WATCH
   t_cpu_pbp.stop();
#endif

   // return new particle bit pattern
   return(p_bit_pattern);
}

// determine the region index of the particle by comparing the
// the masked particle bit pattern with the region bit patterns,
// return -1 if the particle is outside the object,
// the integer input variable can be used for an index estimation
inline int MC_object::get_region_index(const MCBitSet &p_bit_pattern,
                                       const int      &estimate_i_region)
{
#ifdef CPU_WATCH
   extern cpu_watch t_cpu_ind; t_cpu_ind.start();
#endif

   MCBitSet test_bit_set(p_bit_pattern);

   // test whether we can use the estimated index
   if ( (estimate_i_region > 0) && (estimate_i_region < int(num_regions)) )
   {
      // use estimated index

      // the new region index should be close to the estimated region index,
      // therefore we start the index search from this estimated index
      for (register int i1=estimate_i_region-1, i2=estimate_i_region;
                        ( (0 <= i1) || (i2<int(num_regions)) ); --i1, ++i2)
      {
         if (i1 >= 0)
         {
            if ( ( test_bit_set.bit_and(p_bit_pattern,*bit_mask[i1]) )
                                                        == *bit_pattern[i1] )
            {
#ifdef CPU_WATCH
               t_cpu_ind.stop();
#endif
               return(i1);
            }
         }

         if (i2 < int(num_regions))
         {
            if ( ( test_bit_set.bit_and(p_bit_pattern,*bit_mask[i2]) )
                                                     == *bit_pattern[i2] )
            {
#ifdef CPU_WATCH
               t_cpu_ind.stop();
#endif
               return(i2);
            }
         }
      }
   }
   else
   {
      // don't use estimated index
      for (register unsigned i0=0; i0<num_regions; ++i0)
      {
         if ( ( test_bit_set.bit_and(p_bit_pattern,*bit_mask[i0]) )
                                                     == *bit_pattern[i0] )
         {
#ifdef CPU_WATCH
            t_cpu_ind.stop();
#endif
            return(i0);
         }
      }
   }
 
#ifdef CPU_WATCH
   t_cpu_ind.stop();
#endif

   return(-1);
}

// transport photon through the object, the plane pointer points to
// the starting plane of the photon (it must not be the NULL pointer),
// if the photon survives (or one daughter particle) -> return true,
// if no particle survives  -> return false,
// the particle_parameters struct returns at most one particle
inline bool MC_object::primary_photon(particle_parameters &p,
                                      MC_plane *p_plane, ranmar &rndm)
{
   // at the beginning of this routine the primary photon is still alive
   bool primary_photon_alive = true;

#ifdef SECONDARY_MODIFIER_ELECTRONS
   // we create a linked list of secondary particles,
   // the pointers point to the first, last and new entries in this list
   particle_parameters *p_first = NULL;
   particle_parameters *p_last  = NULL;
   particle_parameters *p_new   = NULL;

   // number of secondary particles in linked list
   int num_secondaries = 0;

   // Compton electron
   particle_parameters e_comp;
   e_comp.type     = ELECTRON;
   e_comp.weight   = p.weight;
   e_comp.i.x      = 0;
   e_comp.i.y      = 0;
   e_comp.i.z      = 0;
   e_comp.next     = NULL;

   // pair electron
   particle_parameters e_pair;
   e_pair.type     = ELECTRON;
   e_pair.weight   = p.weight;
   e_pair.i.x      = 0;
   e_pair.i.y      = 0;
   e_pair.i.z      = 0;
   e_pair.next     = NULL;

   // pair positron
   particle_parameters p_pair;
   p_pair.type     = POSITRON;
   p_pair.weight   = p.weight;
   p_pair.i.x      = 0;
   p_pair.i.y      = 0;
   p_pair.i.z      = 0;
   p_pair.next     = NULL;

   // photo electron
   particle_parameters e_phot;
   e_phot.type     = ELECTRON;
   e_phot.weight   = p.weight;
   e_phot.i.x      = 0;
   e_phot.i.y      = 0;
   e_phot.i.z      = 0;
   e_phot.next     = NULL;
#endif

   // there must be an error if p_plane points to NULL
   if (p_plane == NULL)
   {
      xvmc_error("MC_object::primary_photon",
                 "the starting plane is undefined",8);
   }

   // determine the relationship of the photon to p_plane by taking only the
   // photon direction into account, compare with the plane normal vector
   real_3 normal = p_plane->get_normal();
            
   // scalar product of photon direction and plane normal vector
   real product = p.dir.x*normal.x + p.dir.y*normal.y + p.dir.z*normal.z;

   bool p_plane_bit;
   if (product > ZERO)
   {
      // the photon points to the region before the plane, clear bit
      p_plane_bit = false;
   }
   else
   {
      // the photon points to the region behind the plane, set bit
      // or it will stay at the plane, set also bit because
      // in this case we define the photon to be behind the plane
      p_plane_bit = true;
   }

   // determine the photon bit pattern for all planes
   // including the starting plane (set by bool p_plane_bit)
   MCBitSet p_bit_pattern = calc_bit_pattern(p,p_plane,p_plane_bit);

   // estimate and get region index of the photon
   int old_i_region = estimate_region_index(p.pos);
   int p_i_region   = get_region_index(p_bit_pattern,old_i_region);

   // there must be an error if the photon is outside of the object
   if (p_i_region < 0)
   {
      xvmc_error("MC_object::primary_photon",
                 "the photon is outside of the object",8);
   }

   // number of mean free paths to the next interaction point
   real num_mfp = -log(ONE-rndm.number());

   // now trace photon through the object
   while (p_i_region >= 0)
   {
      // distance to the next plane and next plane
      real distance = piece[p_i_region]->distance(p,p_plane);

      if (p_plane != NULL)
      {
         // total Compton cross section
         real mu_comp = piece[p_i_region]->tot_comp->get(p.energy);

         // total pair cross section
         real mu_pair = piece[p_i_region]->tot_pair->get(p.energy);

         // total photo cross section
         real mu_phot = piece[p_i_region]->tot_phot->get(p.energy);

         // total attenuation coefficient
         real mu_tot = mu_comp + mu_pair + mu_phot;

         // Compton interaction probability
         real prob_comp = mu_comp/mu_tot;

         // pair interaction probability
         real prob_pair = mu_pair/mu_tot;

         // first interval:  Compton interaction
         // second interval: pair production
         // third interval:  photo-electric absorption
         prob_pair += prob_comp;

         // path length to the next interaction point
         real path_length = num_mfp/mu_tot;

         // compare photon path length and distance
         if (path_length < distance)
         {
            // move photon to the interaction site
            p.pos.x += path_length*p.dir.x;
            p.pos.y += path_length*p.dir.y;
            p.pos.z += path_length*p.dir.z;

            // the photon is within the region (not at a surface plane)
            p_plane = NULL;

            // determine interaction type
            real eta_itype = rndm.number();
            if (eta_itype <= prob_comp)
            {
               // Compton scattering
               real rnno = rndm.number();  // just a random number
               real energy_cx;             // Compton photon energy
               real energy_ce;             // Compton electron energy
               real cos_t_cx,sin_t_cx;     // Compton photon scattering angle
               real cos_t_ce,sin_t_ce;     // Compton electron scattering angle
               piece[p_i_region]->compton->interaction(p.energy, rnno, rndm,
                                              energy_cx, cos_t_cx, sin_t_cx,
                                              energy_ce, cos_t_ce, sin_t_ce);

               // azimuthal scattering angle
               real phi = TWO_PI * rndm.number();
               real sin_phi = sin(phi);
               real cos_phi = cos(phi);

#ifdef SECONDARY_MODIFIER_ELECTRONS
               if (energy_ce > e_cut)
               {
                  // simulate Compton electron
                  e_comp.energy = energy_ce;
                  e_comp.pos.x  = p.pos.x;
                  e_comp.pos.y  = p.pos.y;
                  e_comp.pos.z  = p.pos.z;
                  e_comp.dir.x  = p.dir.x;
                  e_comp.dir.y  = p.dir.y;
                  e_comp.dir.z  = p.dir.z;

                  // change electron direction
                  rotate(cos_t_ce, sin_t_ce, cos_phi, sin_phi, e_comp.dir);

                  // set electron region identifiers
                  MCBitSet ec_bit_pattern = p_bit_pattern;
                  int      ec_i_region    = p_i_region;

                  // transport electron
                  if ( secondary_electron(e_comp,
                                          ec_bit_pattern, ec_i_region,
                                          rndm) )
                  {
                     // the electron survived, create new list entry
                     ++num_secondaries;
                     p_new  = NULL;
                     if ( (p_new = new particle_parameters) == NULL )
                     {
                        xvmc_error("MC_object::primary_photon",
                           "cannot allocate memory for Compton electron",8);
                     }

                     // store parameters
                     p_new->type     = e_comp.type;
                     p_new->energy   = e_comp.energy;
                     p_new->weight   = e_comp.weight;
                     p_new->pos.x    = e_comp.pos.x;
                     p_new->pos.y    = e_comp.pos.y;
                     p_new->pos.z    = e_comp.pos.z;
                     p_new->dir.x    = e_comp.dir.x;
                     p_new->dir.y    = e_comp.dir.y;
                     p_new->dir.z    = e_comp.dir.z;
                     p_new->i.x      = e_comp.i.x;
                     p_new->i.y      = e_comp.i.y;
                     p_new->i.z      = e_comp.i.z;
                     p_new->next     = e_comp.next;

                     // assign to p_first if this is the first secondary,
                     // assign to the next list entry otherwise
                     if (p_first == NULL) p_first      = p_new;
                     else                 p_last->next = p_new;

                     // the electron is now the last particle in the list
                     p_last = p_new;
                  }
               }
#endif
               // stop transport if the Compton photon energy is below p_cut
               if (energy_cx <= p_cut)
               {
                  // stop transport
                  p_i_region = -1;

                  // and kill photon
                  primary_photon_alive = false;
               }
               else
               {
                  // rotate direction of the photon
                  rotate(cos_t_cx, sin_t_cx, cos_phi, sin_phi, p.dir);

                  // change energy
                  p.energy = energy_cx;

                  // sample new number of mean free paths
                  num_mfp = -log(ONE-rndm.number());
               }
            }
            else
            {
#ifdef SECONDARY_MODIFIER_ELECTRONS
               // pair production or photoelectric absorption,
               if (eta_itype < prob_pair)
               {
                  // pair production
                  real rnno = rndm.number();  // just a random number
                  real energy_pe;             // pair electron energy
                  real energy_pp;             // pair positron energy
                  real cos_t_pe,sin_t_pe;     // pair electron scattering angle
                  real cos_t_pp,sin_t_pp;     // pair positron scattering angle
                  piece[p_i_region]->pair->interaction(p.energy, rnno, rndm,
                                              energy_pe, cos_t_pe, sin_t_pe,
                                              energy_pp, cos_t_pp, sin_t_pp);

                  // azimuthal scattering angle
                  real phi = TWO_PI * rndm.number();
                  real sin_phi = sin(phi);
                  real cos_phi = cos(phi);

                  if (energy_pe > e_cut)
                  {
                     // simulate pair electron
                     e_pair.energy = energy_pe;
                     e_pair.pos.x  = p.pos.x;
                     e_pair.pos.y  = p.pos.y;
                     e_pair.pos.z  = p.pos.z;
                     e_pair.dir.x  = p.dir.x;
                     e_pair.dir.y  = p.dir.y;
                     e_pair.dir.z  = p.dir.z;

                     // change electron direction
                     rotate(cos_t_pe, sin_t_pe, cos_phi, sin_phi, e_pair.dir);

                     // set pair electron region identifiers
                     MCBitSet ep_bit_pattern = p_bit_pattern;
                     int      ep_i_region    = p_i_region;

                     // transport electron
                     if ( secondary_electron(e_pair,
                                             ep_bit_pattern, ep_i_region,
                                             rndm) )
                     {
                        // the electron survived, create new list entry
                        ++num_secondaries;
                        p_new  = NULL;
                        if ( (p_new = new particle_parameters) == NULL )
                        {
                           xvmc_error("MC_object::primary_photon",
                              "cannot allocate memory for pair electron",8);
                        }

                        // store parameters
                        p_new->type     = e_pair.type;
                        p_new->energy   = e_pair.energy;
                        p_new->weight   = e_pair.weight;
                        p_new->pos.x    = e_pair.pos.x;
                        p_new->pos.y    = e_pair.pos.y;
                        p_new->pos.z    = e_pair.pos.z;
                        p_new->dir.x    = e_pair.dir.x;
                        p_new->dir.y    = e_pair.dir.y;
                        p_new->dir.z    = e_pair.dir.z;
                        p_new->i.x      = e_pair.i.x;
                        p_new->i.y      = e_pair.i.y;
                        p_new->i.z      = e_pair.i.z;
                        p_new->next     = e_pair.next;

                        // assign to p_first if this is the first secondary,
                        // assign to the next list entry otherwise
                        if (p_first == NULL) p_first      = p_new;
                        else                 p_last->next = p_new;

                        // the electron is now the last particle in the list
                        p_last = p_new;
                     }
                  }

                  if (energy_pp > e_cut)
                  {
                     // simulate pair positron
                     p_pair.energy = energy_pp;
                     p_pair.pos.x  = p.pos.x;
                     p_pair.pos.y  = p.pos.y;
                     p_pair.pos.z  = p.pos.z;
                     p_pair.dir.x  = p.dir.x;
                     p_pair.dir.y  = p.dir.y;
                     p_pair.dir.z  = p.dir.z;

                     // change positron direction
                     sin_phi = -sin_phi;
                     cos_phi = -cos_phi;
                     rotate(cos_t_pp, sin_t_pp, cos_phi, sin_phi, p_pair.dir);

                     // set pair positron region identifiers
                     MCBitSet pp_bit_pattern = p_bit_pattern;
                     int      pp_i_region    = p_i_region;

                     // transport positron
                     if ( secondary_electron(p_pair,
                                             pp_bit_pattern, pp_i_region,
                                             rndm) )
                     {
                        // the positron survived, create new list entry
                        ++num_secondaries;
                        p_new  = NULL;
                        if ( (p_new = new particle_parameters) == NULL )
                        {
                           xvmc_error("MC_object::primary_photon",
                              "cannot allocate memory for pair positron",8);
                        }

                        // store parameters
                        p_new->type     = p_pair.type;
                        p_new->energy   = p_pair.energy;
                        p_new->weight   = p_pair.weight;
                        p_new->pos.x    = p_pair.pos.x;
                        p_new->pos.y    = p_pair.pos.y;
                        p_new->pos.z    = p_pair.pos.z;
                        p_new->dir.x    = p_pair.dir.x;
                        p_new->dir.y    = p_pair.dir.y;
                        p_new->dir.z    = p_pair.dir.z;
                        p_new->i.x      = p_pair.i.x;
                        p_new->i.y      = p_pair.i.y;
                        p_new->i.z      = p_pair.i.z;
                        p_new->next     = p_pair.next;

                        // assign to p_first if this is the first secondary,
                        // assign to the next list entry otherwise
                        if (p_first == NULL) p_first      = p_new;
                        else                 p_last->next = p_new;

                        // the positron is now the last particle in the list
                        p_last = p_new;
                     }
                  }

                  // stop transport
                  p_i_region = -1;

                  // because the photon has been killed
                  primary_photon_alive = false;
               }
               else
               {
                  // photoelectric absorption, transform photon into electron
                  if (p.energy > e_cut)
                  {
                     // simulate photo electron
                     e_phot.energy = p.energy;
                     e_phot.pos.x  = p.pos.x;
                     e_phot.pos.y  = p.pos.y;
                     e_phot.pos.z  = p.pos.z;
                     e_phot.dir.x  = p.dir.x;
                     e_phot.dir.y  = p.dir.y;
                     e_phot.dir.z  = p.dir.z;

                     // set photo electron region identifiers
                     MCBitSet ea_bit_pattern = p_bit_pattern;
                     int      ea_i_region    = p_i_region;

                     // transport electron
                     if ( secondary_electron(e_phot,
                                             ea_bit_pattern, ea_i_region,
                                             rndm) )
                     {
                        // the electron survived, create new list entry
                        ++num_secondaries;
                        p_new  = NULL;
                        if ( (p_new = new particle_parameters) == NULL )
                        {
                           xvmc_error("MC_object::primary_photon",
                              "cannot allocate memory for photo electron",8);
                        }

                        // store parameters
                        p_new->type     = e_phot.type;
                        p_new->energy   = e_phot.energy;
                        p_new->weight   = e_phot.weight;
                        p_new->pos.x    = e_phot.pos.x;
                        p_new->pos.y    = e_phot.pos.y;
                        p_new->pos.z    = e_phot.pos.z;
                        p_new->dir.x    = e_phot.dir.x;
                        p_new->dir.y    = e_phot.dir.y;
                        p_new->dir.z    = e_phot.dir.z;
                        p_new->i.x      = e_phot.i.x;
                        p_new->i.y      = e_phot.i.y;
                        p_new->i.z      = e_phot.i.z;
                        p_new->next     = e_phot.next;

                        // assign to p_first if this is the first secondary,
                        // assign to the next list entry otherwise
                        if (p_first == NULL) p_first      = p_new;
                        else                 p_last->next = p_new;

                        // the electron is now the last particle in the list
                        p_last = p_new;
                     }
                  }

                  // stop transport
                  p_i_region = -1;

                  // because the photon has been killed
                  primary_photon_alive = false;
               }
#else
               // pair production or photoelectric absorption,
               // stop transport because secondary electrons are neglected
               p_i_region = -1;

               // and kill photon
               primary_photon_alive = false;
#endif
            }
         }
         else
         {
            // move the photon to this plane
            p.pos.x += distance*p.dir.x;
            p.pos.y += distance*p.dir.y;
            p.pos.z += distance*p.dir.z;

            if (p_plane == final_plane)
            {
               // the photon is at the final plane, stop transport
               p_i_region = -1;
            }
            else
            {
               // reduce the number of mean free paths
               num_mfp -= distance*mu_tot;

               // the photon crosses the plane, check the corresponding bit
               // if the bit is set the new bit will be cleared (false)
               // if the bit is not set the new bit will be set (true)
               // that is, we invert the corresponding bit
               if (p_bit_pattern.test(p_plane->index)) p_plane_bit = false;
               else                                    p_plane_bit = true;

               // set this bit and calculate all other bits
               p_bit_pattern = calc_bit_pattern(p,p_plane,p_plane_bit);

               // determine the new region index of the photon, it should
               // be close to the present region, therefore we start the
               // region search at the present region (index "old_i_region")
               old_i_region = p_i_region;
               p_i_region = get_region_index(p_bit_pattern,old_i_region);
            }
         }
      }
      else
      {
         // the distance to the next plane is infinity
         p_i_region = -1;
      }
   }

#ifdef SECONDARY_MODIFIER_ELECTRONS
   // we sample only one of all primary and secondary particles
   if (primary_photon_alive)
   {
      // the primary photon is still alive
      if (p_first != NULL)
      {
         // there are also secondary particles
         int   num_particles = num_secondaries + 1;
         float sum_particles = float(num_particles);
         int p_index = int(rndm.number()*sum_particles); // int random number
         p_index     = max_of(0,p_index);
         p_index     = min_of(num_particles-1,p_index);
         p_new       = &p;           // initialize loop pointer
         register int i_entry = 0;   // index of list entry
         while (p_new != NULL)
         {
            if (p_index == i_entry)
            {
               // move parameters only for secondary particles,
               // change only the weight if we return the primary photon
               if (p_new == &p)
               {
                  p.weight *= sum_particles;
               }
               else
               {
                  p.type     = p_new->type;
                  p.energy   = p_new->energy;
                  p.weight   = p_new->weight*sum_particles;
                  p.pos.x    = p_new->pos.x;
                  p.pos.y    = p_new->pos.y;
                  p.pos.z    = p_new->pos.z;
                  p.dir.x    = p_new->dir.x;
                  p.dir.y    = p_new->dir.y;
                  p.dir.z    = p_new->dir.z;
               }
            }

            // old particle in linked list
            p_last = p_new;

            // next active particle
            if (p_last == &p)
            {
               p_new = p_first;
            }
            else
            {
               p_new = p_last->next;
               delete  p_last;  // delete old particle
            }

            // increase list entry index
            ++i_entry;
         }
      }

      // there is at least one surviving particle
      return(true);
   }

   if (p_first != NULL)
   {
      // there are only secondary patricles

      // we need an interger random number
      float sum_particles = float(num_secondaries);
      int p_index = int(rndm.number()*sum_particles);
      p_index     = max_of(0,p_index);
      p_index     = min_of(num_secondaries-1,p_index);
      p_new       = p_first;  // initialize loop pointer
      register int i_entry = 0;   // index of list entry
      while (p_new != NULL)
      {
         if (p_index == i_entry)
         {
            // store parameters
            p.type     = p_new->type;
            p.energy   = p_new->energy;
            p.weight   = p_new->weight*sum_particles;
            p.pos.x    = p_new->pos.x;
            p.pos.y    = p_new->pos.y;
            p.pos.z    = p_new->pos.z;
            p.dir.x    = p_new->dir.x;
            p.dir.y    = p_new->dir.y;
            p.dir.z    = p_new->dir.z;
         }

         // old particle in linked list
         p_last = p_new;

         // next active particle
         p_new  = p_last->next;

         // delete old particle
         delete p_last;

         // increase list entry index
         ++i_entry;
      }

      // there is a surviving secondary particle
      return(true);
   }

   // there is neither a primary nor a secondary surviving particle
   p.energy = ZERO;
   p.weight = ZERO;
   return(false);
#else
   if (primary_photon_alive) return(true);

   // there is no surviving particle, set energy and weight to zero
   p.energy = ZERO;
   p.weight = ZERO;
   return(false);
#endif
}

// transport primary electron through the object, the plane pointer points to
// the starting plane of the electron (it must not be the NULL pointer),
// if the electron survives (and/or daughter particles) -> return true,
// if no particle survives  -> return false,
// this is just a Continuous Slowing Down Approximation (CSDA),
// the random number generator is reserved for later use
inline bool MC_object::primary_electron(particle_parameters &e,
                                        MC_plane *e_plane, ranmar &rndm)
{
   // there must be an error if e_plane points to NULL
   if (e_plane == NULL)
   {
      xvmc_error("MC_object::primary_electron",
                 "the starting plane is undefined",8);
   }

   // determine the relationship of the electron to e_plane by taking only the
   // electron direction into account, compare with the plane normal vector
   real_3 normal = e_plane->get_normal();
            
   // scalar product of electron direction and plane normal vector
   real product = e.dir.x*normal.x + e.dir.y*normal.y + e.dir.z*normal.z;

   bool e_plane_bit;
   if (product > ZERO)
   {
      // the electron points to the region before the plane, clear bit
      e_plane_bit = false;
   }
   else
   {
      // the electron points to the region behind the plane, set bit
      // or it will stay at the plane, set also bit because
      // in this case we define the electron to be behind the plane
      e_plane_bit = true;
   }

   // determine the electron bit pattern for all planes
   // including the starting plane (set by bool e_plane_bit)
   MCBitSet e_bit_pattern = calc_bit_pattern(e,e_plane,e_plane_bit);

   // estimate and get region index of the electron
   int old_i_region = estimate_region_index(e.pos);
   int e_i_region   = get_region_index(e_bit_pattern,old_i_region);

   // there must be an error if the electron is outside of the object
   if (e_i_region < 0)
   {
      xvmc_error("MC_object::primary_electron",
                 "the electron is outside of the object",8);
   }

   // now trace electron through the object
   while (e_i_region >= 0)
   {
      // distance to the next plane and next plane
      real distance = piece[e_i_region]->distance(e,e_plane);

      if (e_plane != NULL)
      {
         // CSDA range of the electron for the material of the present region
         real csda_range = piece[e_i_region]->e_data->get_csda(e.energy);

         // compare electron CSDA range and distance
         if (csda_range < distance)
         {
            // kill electron
            return(false);
         }
         else
         {
            // move electron to the next plane
            e.pos.x += distance*e.dir.x;
            e.pos.y += distance*e.dir.y;
            e.pos.z += distance*e.dir.z;

            // estimate the new energy of the electron
            real new_energy = e.energy -
                 distance*piece[e_i_region]->e_data->get_dedx(e.energy);
            if (new_energy <= e_cut) return(false);
            e.energy = new_energy;

            if (e_plane == final_plane)
            {
               // the electron is at the final plane, stop transport
               e_i_region = -1;
            }
            else
            {
               // the electron crosses the plane, check the corresponding bit
               // if the bit is set the new bit will be cleared (false)
               // if the bit is not set the new bit will be set (true)
               // that is, we invert the corresponding bit
               if (e_bit_pattern.test(e_plane->index)) e_plane_bit = false;
               else                                    e_plane_bit = true;

               // set this bit and calculate all other bits
               e_bit_pattern = calc_bit_pattern(e,e_plane,e_plane_bit);

               // determine the new region index of the electron, it should
               // be close to the present region, therefore we start the
               // region search at the present region (index "old_i_region")
               old_i_region = e_i_region;
               e_i_region = get_region_index(e_bit_pattern,old_i_region);
            }
         }
      }
      else
      {
         // the distance to the next plane is infinity
         e_i_region = -1;
      }
   }

   // the electron survived, return true
   return(true);
}

// transport secondary electron through the object,
// the bit pattern and the region index of the electron must be known,
// if the electron survives (and/or daughter particles) -> return true,
// if no particle survives  -> return false,
// this is just a Continuous Slowing Down Approximation (CSDA),
// the random number generator is reserved for later use
inline bool MC_object::secondary_electron(particle_parameters &e,
                                          MCBitSet &e_bit_pattern,
                                          int      &e_i_region,
                                          ranmar   &rndm)
{
   // there must be an error if the electron is outside of the object
   if (e_i_region < 0)
   {
      xvmc_error("MC_object::secondary_electron",
                 "the electron is outside of the object",8);
   }

   MC_plane *e_plane = NULL;
   bool      e_plane_bit;
   int       old_i_region;

   // now trace electron through the object
   while (e_i_region >= 0)
   {
      // distance to the next plane and next plane
      real distance = piece[e_i_region]->distance(e,e_plane);

      if (e_plane != NULL)
      {
         // CSDA range of the electron for the material of the present region
         real csda_range = piece[e_i_region]->e_data->get_csda(e.energy);

         // compare electron CSDA range and distance
         if (csda_range < distance)
         {
            // kill electron
            return(false);
         }
         else
         {
            // move electron to the next plane
            e.pos.x += distance*e.dir.x;
            e.pos.y += distance*e.dir.y;
            e.pos.z += distance*e.dir.z;

            // estimate the new energy of the electron
            real new_energy = e.energy -
                 distance*piece[e_i_region]->e_data->get_dedx(e.energy);
            if (new_energy <= e_cut) return(false);
            e.energy = new_energy;

            if (e_plane == final_plane)
            {
               // the electron is at the final plane, stop transport
               e_i_region = -1;
            }
            else
            {
               // the electron crosses the plane, check the corresponding bit
               // if the bit is set the new bit will be cleared (false)
               // if the bit is not set the new bit will be set (true)
               // that is, we invert the corresponding bit
               if (e_bit_pattern.test(e_plane->index)) e_plane_bit = false;
               else                                    e_plane_bit = true;

               // set this bit and calculate all other bits
               e_bit_pattern = calc_bit_pattern(e,e_plane,e_plane_bit);

               // determine the new region index of the electron, it should
               // be close to the present region, therefore we start the
               // region search at the present region (index "old_i_region")
               old_i_region = e_i_region;
               e_i_region = get_region_index(e_bit_pattern,old_i_region);
            }
         }
      }
      else
      {
         // the distance to the next plane is infinity
         e_i_region = -1;
      }
   }

   // the electron survived, return true
   return(true);
}

// ****************************************************************
// class MC_doubleslab: two slabs defined by 3 parallel planes
// ****************************************************************

class MC_doubleslab : public MC_object
{
   private:
      // initialize doubleslab by three existing parallel planes
      // and two slab materials
      void init(MC_plane *, MC_plane *, MC_plane *, char *, char *);

   public:
      // define doubleslab by three existing parallel planes and two materials
      MC_doubleslab(MC_plane *plane0, MC_plane *plane1, MC_plane *plane2,
                    char *material0,  char *material1)
         : MC_object(3,2) { init(plane0,plane1,plane2,material0,material1); }

      // define double slab perpendicular to the x, y or z axis
      // and two materials
      MC_doubleslab(const axis &, const real &, const real &, const real &,
                    char *,  char *);

      // define doubleslab by five points in 3D space and two materials,
      // the first three points define one plane, the two other points are to
      // define two further planes parallel to plane one
      MC_doubleslab(const real_3 &, const real_3 &, const real_3 &,
                    const real_3 &, const real_3 &,
                    char *,  char *);
};

// ****************************************************************
// class MC_jaws_focus: a pair of focussing jaws
// ****************************************************************

class MC_jaws_focus : public MC_object
{
   private:
      // initialize jaws pair
      void init(const axis &,
                const real &, const real &,
                const real &, const real &,
                const real &, const real &,
                const real &, const real &,
                const real &, const real &,
                const char *, const char *);

   public:
      // define x or y jaws pair with given opening, z boundaries,
      // jaws bar material and opening material names
      MC_jaws_focus(const axis &, const real &, const real &,
                                  const real &, const real &,
                                  const char *, const char *);

      // delete focussing jaws
      ~MC_jaws_focus(void);

      // estimate the region index of a point
      int estimate_region_index(const real_3 &);

      // change left jaw position
      void change_left(const real &);

      // change right jaw position
      void change_right(const real &);

      // change left and right jaw positions
      void change_opening(const real &new_left, const real &new_right)
         { change_left(new_left); change_right(new_right); }

   private:
      // X or Y jaws
      axis   type;

      // ESTAR material electron transport data file paths
      char *bars_material_estar;
      char *open_material_estar;

      // NIST material cross section file paths
      char *bars_material_nist;
      char *open_material_nist;

      // differential Compton and pair cross section file paths
      char *compton_file;
      char *pair_file;

      // the outer object dimensions
      real   x_min,x_max,y_min,y_max,z_min,z_max;

      // origin to iso-center distance (cm)
      real   iso_distance;

      // z position of the focus
      real   z_focus;

      // left and right jaws positions projected to the iso center plane (cm)
      real   open_left,open_right;

      // plane pointers
      MC_plane *plane_x_min, *plane_x_max, *plane_y_min, *plane_y_max,
               *plane_z_min, *plane_z_max, *plane_left,  *plane_right;
};

// estimate the region index of a point
inline int MC_jaws_focus::estimate_region_index(const real_3 &p0)
{
   if (type == X)
   {
      // scale x to the iso-center plane
      real x0 = p0.x*iso_distance/p0.z;

      // determine region
      if (x0 < open_left)  return(0);
      if (x0 > open_right) return(1);
   }
   else
   {
      // scale y to the iso-center plane
      real y0 = p0.y*iso_distance/p0.z;

      // determine region
      if (y0 < open_left)  return(0);
      if (y0 > open_right) return(1);
   }

   // return intermediate region
   return(2);
}

// ****************************************************************
// class MC_mlc: multi-leaf collimator (MLC) base class
// ****************************************************************

class MC_mlc : public MC_object
{
   protected:
      // initialize MLC
      void init(multi_leaf *,
                const real &, const real &,
                const real &, const real &,
                const real &, const real &, const real &,
                const char *);

   public:
      // define MLC using the nominal MLC, the opening material name,
      // the source to iso-center distance as well as
      // the numbers of object planes and regions
      MC_mlc(multi_leaf *, const char *, const real, unsigned, unsigned);

      // delete MLC
      virtual ~MC_mlc(void);

      // get MLC type
      mlc_type get_type(void) { return(mlctype); }

      // get leaf moving direction X or Y
      char get_xy(void) {
         return( (xytype == X) ? 'X' : ( (xytype == Y) ? 'Y' : 'Z') ); }

      // get number of leaf pairs
      int get_num(void) {  return(num_pairs); }

      // estimate the region index of a point
      virtual int estimate_region_index(const real_3 &)=0;

      // change position of one leaf pair,
      // the new leaf positions are defined at the iso-center plane,
      // this base class function only performs checks and changes the
      // nominal MLC data, the real leaf positions must be changed by
      // the corresponding derived class member function
      virtual void change_pair(const int &, const real &, const real &);

   protected:
      // MLC type (SIMPLE_MLC, DBLFOCUS_MLC, RNDFOCUS_MLC, ELEKTA_MLC, etc.)
      mlc_type  mlctype;

      // X or Y MLC (leaf moving direction)
      axis      xytype;

      // ESTAR material electron transport data file paths
      char *leaf_material_estar;
      char *open_material_estar;

      // NIST material cross section file paths
      char *leaf_material_nist;
      char *open_material_nist;

      // differential Compton and pair cross section file paths
      char *compton_file;
      char *pair_file;

      // the outer object dimensions
      real   x_min,x_max,y_min,y_max,z_min,z_max;

      // origin to iso-center distance (cm)
      real   iso_distance;

      // z positions of the two focus points
      real   z_focus_wall, z_focus_edge;

      // number of leaf pairs
      int    num_pairs;

      // pointer to the nominal MLC data
      multi_leaf *nominal_mlc;

      // plane pointers (outer object planes)
      MC_plane *plane_x_min, *plane_x_max, *plane_y_min, *plane_y_max,
               *plane_z_min, *plane_z_max;

      // array of pointers to the left and right leaf edge planes
      MC_plane **left_edge_planes, **right_edge_planes;
};


// ****************************************************************
// class MC_jaw: jaws (JAW) base class Added by JK Feb 21, 2008
// ****************************************************************

class MC_jaw : public MC_object
{
   protected:
      // initialize JAW
      void init(jaws *,
                const real &, const real &,
                const real &, const real &,
                const real &, const real &, const real &,
                const char *);

   public:
      // define JAW using the nominal JAW, the opening material name,
      // the source to iso-center distance as well as
      // the numbers of object planes and regions
      MC_jaw(jaws *, const char *, const real, unsigned, unsigned);

      // delete JAW
      virtual ~MC_jaw(void);

      // get JAW type
      jaw_type get_type(void) { return(jawtype); }

      // get jaw moving direction X or Y
      char get_xy(void) {
         return( (xytype == X) ? 'X' : ( (xytype == Y) ? 'Y' : 'Z') ); }

      // get number of jaw pairs (Must be 1)
      int get_num(void) {  return(num_pairs); }

      // estimate the region index of a point
      virtual int estimate_region_index(const real_3 &)=0;

      // change position of jaw pair,
      // the new jaw positions are defined at the iso-center plane,
      // this base class function only performs checks and changes the
      // nominal JAW data, the real jaw positions must be changed by
      // the corresponding derived class member function
      virtual void change_jaws(const int &, const real &, const real &);

   protected:
      // JAW type (SIMPLE_JAW, DBLFOCUS_JAW, REAL_JAW)
      jaw_type  jawtype;

      // X or Y MLC (leaf moving direction)
      axis      xytype;

      // ESTAR material electron transport data file paths
      char *jaw_material_estar;
      char *open_material_estar;

      // NIST material cross section file paths
      char *jaw_material_nist;
      char *open_material_nist;

      // differential Compton and pair cross section file paths
      char *compton_file;
      char *pair_file;

      // the outer object dimensions
      real   x_min,x_max,y_min,y_max,z_min,z_max;

      // origin to iso-center distance (cm)
      real   iso_distance;

      // z positions of the two focus points
      real   z_focus_wall, z_focus_edge;

      // number of leaf pairs
      int    num_pairs;

      // pointer to the nominal MLC data
      jaws *nominal_jaw;

      // plane pointers (outer object planes)
      MC_plane *plane_x_min, *plane_x_max, *plane_y_min, *plane_y_max,
               *plane_z_min, *plane_z_max;

      // pointers to the left and right jaw edge planes
      MC_plane *left_edge_planes, *right_edge_planes;
};

#endif /* _GEOMETRY_H_ */
