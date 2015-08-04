/*****************************************************************************
 * MC_doubleslab.cpp:                                                        *
 *    class member functions for:                                            *
 *       MC_doubleslab:  two slabs defined by 3 parallel planes              *
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

#include "MC_doubleslab.h"

// ****************************************
// member functions of class MC_doubleslab
// ****************************************

// initialize doubleslab by three existing parallel planes
// and two slab materials
void MC_doubleslab::init(MC_plane *plane0, MC_plane *plane1, MC_plane *plane2,
                         char *material0,  char *material1)
{
   // check that the three planes are parallel, plane0 is the reference plane
   const real EPSILON = 1.0e-05;
   real_3 nor0      = plane0->get_normal();
   real_3 nor1      = plane1->get_normal();
   real_3 nor2      = plane2->get_normal();
   real   product01 = nor0.x*nor1.x + nor0.y*nor1.y + nor0.z*nor1.z;
   if (fabs(fabs(product01)-ONE) > EPSILON)
   {
      xvmc_error("MC_doubleslab::init","planes 0 and 1 are not parallel",8);
   }
   real   product02 = nor0.x*nor2.x + nor0.y*nor2.y + nor0.z*nor2.z;
   if (fabs(fabs(product02)-ONE) > EPSILON)
   {
      xvmc_error("MC_doubleslab::init","planes 0 and 2 are not parallel",8);
   }
   real   product12 = nor1.x*nor2.x + nor1.y*nor2.y + nor1.z*nor2.z;
   if (fabs(fabs(product12)-ONE) > EPSILON)
   {
      xvmc_error("MC_doubleslab::init","planes 1 and 2 are not parallel",8);
   }

   // check distances to the zero point
   real dist0 = plane0->get_dist2zero();
   real dist1 = plane1->get_dist2zero();
   real dist2 = plane2->get_dist2zero();

   // sort planes
   MC_plane *new_plane0 = NULL;
   MC_plane *new_plane1 = NULL;
   MC_plane *new_plane2 = NULL;
   if (product01 < ZERO)
   {
      // ensure that the planes are really different
      if (fabs(dist1+dist0) < EPSILON)
      {
         xvmc_error("MC_doubleslab::init","planes 0 and 1 are identical",8);
      }

      // planes 0 and 1 are at different sides of the origin
      if (product02 < ZERO)
      {
         // ensure that the planes are really different
         if (fabs(dist2+dist0) < EPSILON)
         {
            xvmc_error("MC_doubleslab::init","planes 0 and 2 are identical",8);
         }
         if (fabs(dist2-dist1) < EPSILON)
         {
            xvmc_error("MC_doubleslab::init","planes 1 and 2 are identical",8);
         }

         // planes 1 and 2 are at the same side of the origin
         // planes 0 is at the other side
         // choose plane0 to be the new plane0
         new_plane0 = plane0;

         // sort planes 1 and 2
         if (fabs(dist1) < fabs(dist2))
         {
            new_plane1 = plane1;
            new_plane2 = plane2;
         }
         else
         {
            new_plane1 = plane2;
            new_plane2 = plane1;
         }
      }
      else
      {
         // ensure that the planes are really different
         if (fabs(dist2-dist0) < EPSILON)
         {
            xvmc_error("MC_doubleslab::init","planes 0 and 2 are identical",8);
         }
         if (fabs(dist2+dist1) < EPSILON)
         {
            xvmc_error("MC_doubleslab::init","planes 1 and 2 are identical",8);
         }

         // planes 0 and 2 are at the same side of the origin
         // planes 1 is at the other side
         // choose plane1 to be the new plane0
         new_plane0 = plane1;

         // sort planes 0 and 2
         if (fabs(dist0) < fabs(dist2))
         {
            new_plane1 = plane0;
            new_plane2 = plane2;
         }
         else
         {
            new_plane1 = plane2;
            new_plane2 = plane0;
         }
      }
   }
   else
   {
      // ensure that the planes are really different
      if (fabs(dist1-dist0) < EPSILON)
      {
         xvmc_error("MC_doubleslab::init","planes 0 and 1 are identical",8);
      }

      // planes 0 and 1 are at the same side of the origin
      if (product02 < ZERO)
      {
         // ensure that the planes are really different
         if (fabs(dist2+dist0) < EPSILON)
         {
            xvmc_error("MC_doubleslab::init","planes 0 and 2 are identical",8);
         }
         if (fabs(dist2+dist1) < EPSILON)
         {
            xvmc_error("MC_doubleslab::init","planes 1 and 2 are identical",8);
         }

         // planes 0 and 1 are at the same side of the origin
         // planes 2 is at the other side
         // choose plane2 to be the new plane0
         new_plane0 = plane2;

         // sort planes 0 and 1
         if (fabs(dist0) < fabs(dist1))
         {
            new_plane1 = plane0;
            new_plane2 = plane1;
         }
         else
         {
            new_plane1 = plane1;
            new_plane2 = plane0;
         }
      }
      else
      {
         // ensure that the planes are really different
         if (fabs(dist2-dist0) < EPSILON)
         {
            xvmc_error("MC_doubleslab::init","planes 0 and 2 are identical",8);
         }
         if (fabs(dist2-dist1) < EPSILON)
         {
            xvmc_error("MC_doubleslab::init","planes 1 and 2 are identical",8);
         }

         // planes 0, 1 and 2 are at the same side of the origin
         // sort planes
         if (fabs(dist0) < fabs(dist1))
         {
            if (fabs(dist0) < fabs(dist2))
            {
               new_plane0 = plane0;
               if (fabs(dist1) < fabs(dist2))
               {
                  new_plane1 = plane1;
                  new_plane2 = plane2;
               }
               else
               {
                  new_plane1 = plane2;
                  new_plane2 = plane1;
               }
            }
            else
            {
               new_plane0 = plane2;
               new_plane1 = plane0;
               new_plane2 = plane1;
            }
         }
         else
         {
            if (fabs(dist0) < fabs(dist2))
            {
               new_plane0 = plane1;
               new_plane1 = plane0;
               new_plane2 = plane2;
            }
            else
            {
               if (fabs(dist1) < fabs(dist2))
               {
                  new_plane0 = plane1;
                  new_plane1 = plane2;
               }
               else
               {
                  new_plane0 = plane2;
                  new_plane1 = plane1;
               }
               new_plane2 = plane0;
            }
         }
      }
   }  

   // new normal vectors
   nor0 = new_plane0->get_normal();
   nor1 = new_plane1->get_normal();
   nor2 = new_plane2->get_normal();

   // new distances to the origin
   dist0 = new_plane0->get_dist2zero();
   dist1 = new_plane1->get_dist2zero();
   dist2 = new_plane2->get_dist2zero();

   // assign plane pointers
   separator[0] = new_plane0;
   separator[1] = new_plane1;
   separator[2] = new_plane2;

   // determine (ESTAR) material electron data file paths
   char *material0_estar = NULL;
   if ( (material0_estar = get_file_path(material0,"estar")) == NULL )
   {
      xvmc_error("MC_doubleslab::init",
                 "cannot determine electron data file path for material 0",8);
   }
   char *material1_estar = NULL;
   if ( (material1_estar = get_file_path(material1,"estar")) == NULL )
   {
      xvmc_error("MC_doubleslab::init",
                 "cannot determine electron data file path for material 1",8);
   }

   // determine (NIST) material cross section file paths
   char *material0_nist = NULL;
   if ( (material0_nist = get_file_path(material0,"nist")) == NULL )
   {
      xvmc_error("MC_doubleslab::init",
                 "cannot determine cross section file path for material 0",8);
   }
   char *material1_nist = NULL;
   if ( (material1_nist = get_file_path(material1,"nist")) == NULL )
   {
      xvmc_error("MC_doubleslab::init",
                 "cannot determine cross section file path for material 1",8);
   }

   // determine differential Compton cross section file path
   char *compton_file = NULL;
   if ( (compton_file = get_file_path("compton","data")) == NULL )
   {
      xvmc_error("MC_doubleslab::init",
         "cannot determine differential Compton cross section file path",8);
   }

   // determine differential pair cross section file path
   char *pair_file = NULL;
   if ( (pair_file = get_file_path("pair","data")) == NULL )
   {
      xvmc_error("MC_doubleslab::init",
         "cannot determine differential pair cross section file path",8);
   }

   // create first slab using planes 0 and 1
   piece[0] = NULL;
   if ( (piece[0] = new (nothrow) MC_slab(
                                separator[0],    separator[1],
                                material0_estar, material0_nist,
                                compton_file,    pair_file)) == NULL )
   {
      xvmc_error("MC_doubleslab::init","cannot create slab 0",8);
   }

   // create second slab using planes 1 and 2
   piece[1] = NULL;
   if ( (piece[1] = new (nothrow) MC_slab(
                                separator[1],    separator[2],
                                material1_estar, material1_nist,
                                compton_file,    pair_file)) == NULL )
   {
      xvmc_error("MC_doubleslab::init","cannot create slab 1",8);
   }

   // delete cross section and electron transport data file names
   delete [] material0_estar; material0_estar = NULL;
   delete [] material1_estar; material1_estar = NULL;
   delete [] material0_nist; material0_nist = NULL;
   delete [] material1_nist; material1_nist = NULL;
   delete [] compton_file;   compton_file   = NULL;
   delete [] pair_file;      pair_file      = NULL;

   // set bit mask and bit patterns of all regions using the reference points
   set_bits();
}

// define double slab perpendicular to the x, y or z axis
MC_doubleslab::MC_doubleslab(const axis &xyz, const real &d0,
                             const real &d1,  const real &d2,
                             char *material0, char *material1)
             : MC_object(3,2)
{
   // the plane positions must not be identical
   if (d0 == d1)
   {
      xvmc_error("MC_doubleslab::MC_doubleslab",
                 "cannot create double slab with d0 == d1",8);
   }
   if (d0 == d2)
   {
      xvmc_error("MC_doubleslab::MC_doubleslab",
                 "cannot create double slab with d0 == d2",8);
   }
   if (d1 == d2)
   {
      xvmc_error("MC_doubleslab::MC_doubleslab",
                 "cannot create double slab with d1 == d2",8);
   }

   // create 3 parallel planes
   MC_plane *plane0 = NULL;
   if ( (plane0 = new (nothrow) MC_plane(xyz,d0,-1)) == NULL )
   {
      xvmc_error("MC_doubleslab::MC_doubleslab","cannot create plane 0",8);
   }
   MC_plane *plane1 = NULL;
   if ( (plane1 = new (nothrow) MC_plane(xyz,d1,-1)) == NULL )
   {
      xvmc_error("MC_doubleslab::MC_doubleslab","cannot create plane 1",8);
   }
   MC_plane *plane2 = NULL;
   if ( (plane2 = new (nothrow) MC_plane(xyz,d2,-1)) == NULL )
   {
      xvmc_error("MC_doubleslab::MC_doubleslab","cannot create plane 2",8);
   }

   // initialize doubleslab by the three planes
   init(plane0,plane1,plane2,material0,material1);
}

// define doubleslab by five points in 3D space
// the first three points define one plane, the two other points are to
// define two further planes parallel to plane one
MC_doubleslab::MC_doubleslab(const real_3 &p0, const real_3 &p1,
                             const real_3 &p2, const real_3 &p3,
                             const real_3 &p4,
                             char *material0,  char *material1)
             : MC_object(3,2)
{
   // create first plane by three points
   MC_plane *plane0 = NULL;
   if ( (plane0 = new (nothrow) MC_plane(p0,p1,p2,-1)) == NULL )
   {
      xvmc_error("MC_doubleslab::MC_doubleslab",
                 "cannot create plane 0 by three points",8);
   }

   // p3 must not be part of plane0
   if (plane0->relationship(p3) == 0)
   {
      xvmc_error("MC_doubleslab::MC_doubleslab",
                 "point 3 is part of plane 0",8);
   }

   // create plane1 by plane0 and p3
   MC_plane *plane1 = NULL;
   if ( (plane1 = new (nothrow) MC_plane(plane0,p3,-1)) == NULL )
   {
      xvmc_error("MC_doubleslab::MC_doubleslab",
                 "cannot create plane 1 by plane 0 and point 3",8);
   }

   // p4 must not be part of plane0
   if (plane0->relationship(p4) == 0)
   {
      xvmc_error("MC_doubleslab::MC_doubleslab",
                 "point 4 is part of plane 0",8);
   }

   // p4 must not be part of plane1
   if (plane1->relationship(p4) == 0)
   {
      xvmc_error("MC_doubleslab::MC_doubleslab",
                 "point 4 is part of plane 1",8);
   }

   // create plane2 by plane0 and p4
   MC_plane *plane2 = NULL;
   if ( (plane2 = new (nothrow) MC_plane(plane0,p4,-1)) == NULL )
   {
      xvmc_error("MC_doubleslab::MC_doubleslab",
                 "cannot create plane 2 by plane 0 and point 4",8);
   }

   // initialize doubleslab by the three planes
   init(plane0,plane1,plane2,material0,material1);
}
