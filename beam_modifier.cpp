/*****************************************************************************
 * beam_modifier.cpp:                                                        *
 *    class member functions for:                                            *
 *       beam_modifier:      beam modifier base class (empty modifier)       *
 *                                                                           *
 * Copyright (C) 2002    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 31.01.2002      *
 *    dynamic modifier mode implemented                   MF 09.10.2003      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include "beam_modifier.h"

// ****************************************
// declare functions and global variables
// ****************************************

// ****************************************
// member functions of class beam_modifier
// ****************************************

// transport particle through the (empty) modifier, return true
bool beam_modifier::transport(particle_parameters &p, ranmar &rndm)
{
   // length of line from actual particle position to lower plane
   real length = (lower_distance-p.pos.z)/p.dir.z;

   // intersection point of particle path with lower plane
   p.pos.x += length*p.dir.x;
   p.pos.y += length*p.dir.y;
   p.pos.z  = lower_distance;

   return(true);
}

// transport particle position and direction through the
// (empty) modifier, return true
bool beam_modifier::transport(real_3 &pos, real_3 &dir, ranmar &rndm)
{
   // length of line from actual particle position to lower plane
   real length = (lower_distance-pos.z)/dir.z;

   // intersection point of particle path with lower plane
   pos.x += length*dir.x;
   pos.y += length*dir.y;
   pos.z  = lower_distance;

   return(true);
}
