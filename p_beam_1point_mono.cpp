/*****************************************************************************
 * p_beam_1point_mono.cpp:                                                   *
 *    class member functions for:                                            *
 *       p_beam_1point_mono: point source beam, mono-energetic photons       *
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

#include "p_beam_1point_mono.h"

// ***********************************************
// member functions of class p_beam_1point_mono
// ***********************************************

// get particle parameters from mono-energetic photon beam
bool p_beam_1point_mono::emit(particle_parameters &p,
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
   p.type   = type;
   p.energy = sample_energy(rndm);
   p.weight = ONE;

   // sample photon position in upper collimator plane

   // Z-position in the MC starting plane (above the collimators)
   p.pos.z = col_mdistance;

   // X-position in the X-collimator plane
   p.pos.x = col_x1 + col_width_x*rx;
   // X-position in the MC starting plane (above the collimators)
   p.pos.x = p.pos.x*col_mdistance/col_xdistance;

   // Y-position in the Y-collimator plane
   p.pos.y = col_y1 + col_width_y*ry;
   // Y-position in the MC starting plane (above the collimators)
   p.pos.y = p.pos.y*col_mdistance/col_ydistance;

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
