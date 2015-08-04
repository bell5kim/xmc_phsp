/*****************************************************************************
 * e_beam_1point_mono.cpp:                                                   *
 *    class member functions for:                                            *
 *       e_beam_1point_mono: point source beam, mono-energetic electrons     *
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

#include "e_beam_1point_mono.h"

// ***********************************************
// member functions of class e_beam_1point_mono
// ***********************************************

// get particle parameters from mono-energetic electron beam
bool e_beam_1point_mono::emit(particle_parameters &p,
#ifdef USE_SOBOL
                              sobseq &sobol, int &sobol_count,
#endif
                              ranmar &rndm)
{
   real origin_point_dist; // origin to point distance

   // random (Sobol) numbers to sample initial position
#ifdef USE_SOBOL
   real rx = sobol.number(sobol_count); ++sobol_count;
   real ry = sobol.number(sobol_count); ++sobol_count;
#else
   real rx = rndm.number();
   real ry = rndm.number();
#endif 

   // set particle parameters
   p.type   = ELECTRON;
   p.energy = nominal_energy;
   p.weight = ONE;

   // position
   p.pos.x  = app_x1 + app_width_x*rx;
   p.pos.y  = app_y1 + app_width_y*ry;
   p.pos.z  = app_distance;

   // origin to point distance
   origin_point_dist =
      sqrt(p.pos.x*p.pos.x + p.pos.y*p.pos.y + p.pos.z*p.pos.z);

   // direction
   p.dir.x  = p.pos.x/origin_point_dist;
   p.dir.y  = p.pos.y/origin_point_dist;
   p.dir.z  = p.pos.z/origin_point_dist;

   if ( modifier->transport(p,rndm) )
   {
      // trace particle to the simulation grid (take beam angles into account)
      // return true if the particle hits the calculation cube
      return( trace2cube(p.pos, p.dir, p.i, origin_point_dist, rndm) );
   }

   // there isn't any particle to emit, return false
   return(false);
}

// get particle weight, starting position, direction and voxel index
bool e_beam_1point_mono::emit(real   &weight,
                              real_3 &pos, real_3 &dir, int_3 &i,
                              bool   &primary_particle,
#ifdef USE_SOBOL
                              sobseq &sobol, int &sobol_count,
#endif
                              ranmar &rndm)
{
   real origin_point_dist; // origin to point distance

   // random (Sobol) numbers to sample initial position
#ifdef USE_SOBOL
   real rx = sobol.number(sobol_count); ++sobol_count;
   real ry = sobol.number(sobol_count); ++sobol_count;
#else
   real rx = rndm.number();
   real ry = rndm.number();
#endif 

   // set particle parameters
   weight = ONE;

   // position
   pos.x  = app_x1 + app_width_x*rx;
   pos.y  = app_y1 + app_width_y*ry;
   pos.z  = app_distance;

   // origin to point distance
   origin_point_dist = sqrt(pos.x*pos.x + pos.y*pos.y + pos.z*pos.z);

   // direction
   dir.x  = pos.x/origin_point_dist;
   dir.y  = pos.y/origin_point_dist;
   dir.z  = pos.z/origin_point_dist;

   if ( modifier->transport(pos,dir,rndm) )
   {
      // trace particle to the simulation grid (take beam angles into account)
      // return true if the particle hits the calculation cube
      return( trace2cube(pos, dir, i, origin_point_dist, rndm) );
   }

   // there isn't any particle to emit, return false
   return(false);
}
