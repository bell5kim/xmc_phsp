#ifndef _DEFINITIONS_H_
#define _DEFINITIONS_H_

/*****************************************************************************
 * definitions.h: definitions and simple type declarations:                  *
 *       real:                float or double                                *
 *       int_3:               3D integer vector                              *
 *       real_3:              3D real vector                                 *
 *       particle_type:       ELECTRON or PHOTON                             *
 *       result_type:         DOSE, DOSE_TO_H2O, IONIZATION, FILM, ABS_DOSE  *
 *       norm_type:           NONE, MAXIMUM, VALUE or POINT                  *
 *       plane_type:          XY_PLANE, XZ_PLANE, YZ_PLANE                   *
 *       plane_parameters:    output plane type, position etc.               *
 *       profile_type:        X_PROFILE, Y_PROFILE, Z_PROFILE                *
 *       profile_parameters:  output profile type, position etc.             *
 *       particle_parameters: particle energy, position, direction etc.      *
 *       sum_energy_type:     test energy conservation                       *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 00/02/22        *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <stdlib.h>

// in debug mode test energy conservation
#ifdef DEBUG
#define CHECK_ENERGY
#define ONE_PROCESS
#endif

#ifdef ALL_INLINE
#define RANDOM_INLINE
#define ROTATE_INLINE
#define INTERPOLATE_INLINE
#endif

// real is an alias name for double or float (default)
#ifdef DOUBLE_REAL
typedef double real;
#else
typedef float real;
#endif /* DOUBLE_REAL */

// 3d integer vector
struct int_3
{
   int x;
   int y;
   int z;
};

// 3d real (float) vector
struct real_3
{
   real x;
   real y;
   real z;
};

// particle type (photon or electron)
enum particle_type {ELECTRON = -1, PHOTON = 0, POSITRON = 1, NOP = -32767};

// result type
#define MONITORBSF
#ifndef MONITORBSF
enum result_type {DOSE, DOSE_TO_H2O, IONIZATION, FILM, ABS_DOSE};
#else
// Option for Monitor Back Scatter Factor
enum result_type {DOSE, DOSE_TO_H2O, IONIZATION, FILM, ABS_DOSE, ABS_DOSE_MBSF};
#endif
// normalization type
enum norm_type {NONE, MAXIMUM, VALUE, POINT};

// plane type
enum plane_type {XY_PLANE, XZ_PLANE, YZ_PLANE};

// output plane parameters
struct plane_parameters
{
   plane_type           type;    // xy, xz or yz plane
   real                 pos;     // z, y or x position of the plane (cm)
   plane_parameters    *next;    // next plane in linked list
};

// profile type
enum profile_type {X_PROFILE, Y_PROFILE, Z_PROFILE};

// output profile parameters
struct profile_parameters
{
   profile_type         type;    // x, y or z profile
   real_3               pos;     // x,y,z position (cm)
   profile_parameters  *next;    // next profile in linked list
};

// particle parameters
struct particle_parameters
{
   particle_type type;        // PHOTON or ELECTRON
   real          energy;      // kinetic energy (MeV)
   real          weight;      // particle weight
   real_3        pos;         // x,y,z position (cm)
   real_3        dir;         // direction cosine
   int_3         i;           // voxel indices
   particle_parameters *next; // MA next particle in linked list
};

// test energy conservation
#ifdef CHECK_ENERGY
class sum_energy_type
{
   public:
      double start;  // the initial particle energy (photons and electrons)
      double edepo;  // deposited electron energy
      double pdepo;  // directly deposited photon energy
      double eloss;  // energy loss if the electron leaves the simulation cube
      double ploss;  // energy loss if the photon leaves the simulation cube

      void init(void) { start = 0.0;
         edepo = 0.0; pdepo = 0.0; eloss = 0.0; ploss = 0.0; }
      sum_energy_type(void) { init(); }
};
#endif // CHECK_ENERGY

#endif /* _DEFINITIONS_H_ */
