/*****************************************************************************
 * MC_slab.cpp:                                                              *
 *    class member functions for:                                            *
 *       MC_slab:        slab defined by two parallel planes                 *
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

#include <new>
using namespace std;

#include "MC_slab.h"

// ****************************************
// member functions of class MC_slab
// ****************************************

// initialize slab by two existing parallel planes
void MC_slab::init(MC_plane *plane0, MC_plane *plane1)
{
   // check that the two planes are parallel
   const real EPSILON = 1.0e-05;
   real_3 nor0    = plane0->get_normal();
   real_3 nor1    = plane1->get_normal();
   real   product = nor0.x*nor1.x + nor0.y*nor1.y + nor0.z*nor1.z;
   if (fabs(fabs(product)-ONE) > EPSILON)
   {
      xvmc_error("MC_slab::init","the planes are not parallel",8);
   }

   // check distances to the zero point
   real dist0 = plane0->get_dist2zero();
   real dist1 = plane1->get_dist2zero();
   if (product < ZERO)
   {
      if (fabs(dist1+dist0) < EPSILON)
      {
         xvmc_error("MC_slab::init","the planes are identical",8);
      }
   }
   else
   {
      if (fabs(dist1-dist0) < EPSILON)
      {
         xvmc_error("MC_slab::init","the planes are identical",8);
      }
   }  

   // now assign plane pointers
   surface[0] = plane0;
   surface[1] = plane1;

   // calculate reference point
   p_ref.x = (dist0*nor0.x+dist1*nor1.x)/TWO;
   p_ref.y = (dist0*nor0.y+dist1*nor1.y)/TWO;
   p_ref.z = (dist0*nor0.z+dist1*nor1.z)/TWO;
}

// define slab perpendicular to the x, y or z axis and
// the cross section file names
MC_slab::MC_slab(const axis &xyz,
                 const real &d0, const real &d1,
                 char *material_estar, char *material_nist,
                 char *compton_file,   char *pair_file)
       : MC_region(2,material_estar,material_nist,compton_file,pair_file)
{
   // the two plane positions must not be identical
   if (d0 == d1)
   {
      xvmc_error("MC_slab::MC_slab",
                 "cannot create slab with two identical planes",8);
   }

   // create two parallel planes
   MC_plane *plane0 = NULL;
   if ( (plane0 = new (nothrow) MC_plane(xyz,d0,-1)) == NULL )
   {
      xvmc_error("MC_slab::MC_slab","cannot create plane 0",8);
   }
   MC_plane *plane1 = NULL;
   if ( (plane1 = new (nothrow) MC_plane(xyz,d1,-1)) == NULL )
   {
      xvmc_error("MC_slab::MC_slab","cannot create plane 1",8);
   }

   // initialize slab by the two planes
   init(plane0,plane1);
}

// define slab perpendicular to the x, y or z axis and
// a cross section data base for the slab material
MC_slab::MC_slab(const axis &xyz, const real &d0, const real &d1,
                 MC_material *material) : MC_region(2,material)
{
   // the two plane positions must not be identical
   if (d0 == d1)
   {
      xvmc_error("MC_slab::MC_slab",
                 "cannot create slab with two identical planes",8);
   }

   // create two parallel planes
   MC_plane *plane0 = NULL;
   if ( (plane0 = new (nothrow) MC_plane(xyz,d0,-1)) == NULL )
   {
      xvmc_error("MC_slab::MC_slab","cannot create plane 0",8);
   }
   MC_plane *plane1 = NULL;
   if ( (plane1 = new (nothrow) MC_plane(xyz,d1,-1)) == NULL )
   {
      xvmc_error("MC_slab::MC_slab","cannot create plane 1",8);
   }

   // initialize slab by the two planes
   init(plane0,plane1);
}

// define slab by four points in 3D space
// the first three points define one plane, the fourth point is to
// define the second plane parallel to plane one,
// use cross section file names
MC_slab::MC_slab(const real_3 &p0, const real_3 &p1,
                 const real_3 &p2, const real_3 &p3,
                 char *material_estar, char *material_nist,
                 char *compton_file,   char *pair_file)
       : MC_region(2,material_estar,material_nist,compton_file,pair_file)
{
   // create first plane by three points
   MC_plane *plane0 = NULL;
   if ( (plane0 = new (nothrow) MC_plane(p0,p1,p2,-1)) == NULL )
   {
      xvmc_error("MC_slab::MC_slab",
                 "cannot create plane 0 by three points",8);
   }

   // p3 must not be part of plane0
   if (plane0->relationship(p3) == 0)
   {
      xvmc_error("MC_slab::MC_slab","point 3 is part of plane 0",8);
   }

   // create plane1 by plane0 and p3
   MC_plane *plane1 = NULL;
   if ( (plane1 = new (nothrow) MC_plane(plane0,p3,-1)) == NULL )
   {
      xvmc_error("MC_slab::MC_slab",
                 "cannot create plane 1 by plane 0 and point 3",8);
   }

   // initialize slab by the two planes
   init(plane0,plane1);
}

// define slab by four points in 3D space
// the first three points define one plane, the fourth point is to
// define the second plane parallel to plane one,
// use cross section data base for the slab material
MC_slab::MC_slab(const real_3 &p0, const real_3 &p1,
                 const real_3 &p2, const real_3 &p3,
                 MC_material *material) : MC_region(2,material)
{
   // create first plane by three points
   MC_plane *plane0 = NULL;
   if ( (plane0 = new (nothrow) MC_plane(p0,p1,p2,-1)) == NULL )
   {
      xvmc_error("MC_slab::MC_slab",
                 "cannot create plane 0 by three points",8);
   }

   // p3 must not be part of plane0
   if (plane0->relationship(p3) == 0)
   {
      xvmc_error("MC_slab::MC_slab","point 3 is part of plane 0",8);
   }

   // create plane1 by plane0 and p3
   MC_plane *plane1 = NULL;
   if ( (plane1 = new (nothrow) MC_plane(plane0,p3,-1)) == NULL )
   {
      xvmc_error("MC_slab::MC_slab",
                 "cannot create plane 1 by plane 0 and point 3",8);
   }

   // initialize slab by the two planes
   init(plane0,plane1);
}
