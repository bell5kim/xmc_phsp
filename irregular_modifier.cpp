/*****************************************************************************
 * irregular_modifier.cpp:                                                   *
 *    class member functions for:                                            *
 *       irregular_modifier: modifier based on an irregular beam contour     *
 *                                                                           *
 * Copyright (C) 2002    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 31.01.2002      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include "irregular_modifier.h"

// ****************************************
// declare functions and global variables
// ****************************************

// ********************************************
// member functions of class irregular_modifier
// ********************************************

// create irregular beam modifier based on the irregular beam
// contour at the iso-center plane
irregular_modifier::irregular_modifier(contour *iso_contour, real iso_distance)
 : beam_modifier()
{
   // set beam modifier type
   type = MODIFIER_IRREGULAR;

   // get number of contour points
   int num_points = iso_contour->get_num();

   // create contour
   if ( (lower_contour = new contour(num_points)) == NULL)
   {
      xvmc_error("irregular_modifier::irregular_modifier",
                 "cannot create contour",8);
   }

   // scale contour points (iso-center plane to lower modifier plane)
   real scale = lower_distance/iso_distance;
   for (register int i=0; i<num_points; ++i)
   {
      real x = scale*iso_contour->get_x(i);
      real y = scale*iso_contour->get_y(i);
      if (!lower_contour->change_point(i,x,y))
      {
         xvmc_error("irregular_modifier::irregular_modifier",
                    "cannot change irregular beam contour point",8);
      }
   }

   // close contour
   if (!lower_contour->close())
   {
      xvmc_error("irregular_modifier::irregular_modifier",
                 "cannot close irregular beam contour",8);
   }
}

// create irregular beam modifier based on the irregular beam
// contour at the iso-center plane and the position of this contour
irregular_modifier::irregular_modifier(contour *iso_contour,
                                       real iso_distance,
                                       real distance)
 : beam_modifier(distance-ONE,distance)
{
   // set beam modifier type
   type = MODIFIER_IRREGULAR;

   // get number of contour points
   int num_points = iso_contour->get_num();

   // create contour
   if ( (lower_contour = new contour(num_points)) == NULL)
   {
      xvmc_error("irregular_modifier::irregular_modifier",
                 "cannot create contour",8);
   }

   // scale contour points (iso-center plane to lower modifier plane)
   real scale = lower_distance/iso_distance;
   for (register int i=0; i<num_points; ++i)
   {
      real x = scale*iso_contour->get_x(i);
      real y = scale*iso_contour->get_y(i);
      if (!lower_contour->change_point(i,x,y))
      {
         xvmc_error("irregular_modifier::irregular_modifier",
                    "cannot change irregular beam contour point",8);
      }
   }

   // close contour
   if (!lower_contour->close())
   {
      xvmc_error("irregular_modifier::irregular_modifier",
                 "cannot close irregular beam contour",8);
   }
}

// transport particle through the irregular beam modifier,
// if the particle is deleted return false, otherwise return true
bool irregular_modifier::transport(particle_parameters &p, ranmar &rndm)
{
   // length of line from actual particle position to lower plane
   real length = (lower_distance-p.pos.z)/p.dir.z;

   // intersection point of particle path with lower plane
   p.pos.x += length*p.dir.x;
   p.pos.y += length*p.dir.y;
   p.pos.z  = lower_distance;

   // test whether the particle is inside or outside the irregular contour
   if (lower_contour->inside(p.pos.x,p.pos.y)) return(true);

   // the particle is outside
   return(false);
}

// transport position and direction through the irregular beam modifier,
// if the particle is deleted return false, otherwise return true
bool irregular_modifier::transport(real_3 &pos, real_3 &dir, ranmar &rndm)
{
   // length of line from actual particle position to lower plane
   real length = (lower_distance-pos.z)/dir.z;

   // intersection point of particle path with lower plane
   pos.x += length*dir.x;
   pos.y += length*dir.y;
   pos.z  = lower_distance;

   // test whether the particle is inside or outside the irregular contour
   if (lower_contour->inside(pos.x,pos.y)) return(true);

   // the particle is outside
   return(false);
}
