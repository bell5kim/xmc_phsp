/*****************************************************************************
 * portal_dose.cpp:                                                          *
 *    class member functions for:                                            *
 *       portal_dose:    portal dose image (dummy)                           *
 *                                                                           *
 * Copyright (C) 2003    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 03/03/13        *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <iostream>
// Added by JOKim 25NOV2009 -------------
#include <fstream>
#include <sstream>
#include <iomanip>
using namespace std;
/*
#include <xvmc_util.h>
#include <photon_data.h>
#include <math.h>
#include "definitions.h"
#include "global.h"
#include <compton_XS_diff.h>
#include <pair_XS_diff.h>
#include <electron_data.h>
#include <moller_XS.h>
#include <bhabha_XS.h>
#include <brems_XS.h>
*/

#define ELECTRON_TRACK_FULL_STOP 1
// End of Added -------------------------
#include "portal_dose.h"

// ****************************************
// declare global variables
// ****************************************

// Added by JOKim 25NOV2009 -------------

// simulate photon histories by KERMA approximation
void kerma_photon_portal(particle_parameters &p, ranmar &rndm,
#ifdef CHECK_ENERGY
                  sum_energy_type &sum_energy,
#endif // CHECK_ENERGY
                  array_3d<double> *batch_dose_portal);

// simulate photon histories by multiple photon transport
void multi_photon_portal(particle_parameters &, int, ranmar &,
#ifdef USE_SOBOL
                  real *,
#endif // USE_SOBOL
#ifdef CHECK_ENERGY
                  sum_energy_type &,
#endif // CHECK_ENERGY
                  array_3d<double> *);

// simulate electron history
void one_electron_portal(particle_parameters &, ranmar &,
#ifdef CHECK_ENERGY
                  sum_energy_type &,
#endif // CHECK_ENERGY
                  array_3d<double> *);


void   swap_bytes( float &); // swap bytes
// End of Added -------------------------
// ****************************************
// member functions of class portal_dose
// ****************************************

// define portal dose image plane by the distance to the target
// the image matrix dimension, the image matrix resolution and
// the BMP image size
portal_dose::portal_dose(double       new_distance,
                         unsigned int new_dim_x,
                         unsigned int new_dim_y,
                         unsigned int new_dim_z,
                         double       new_res_x,
                         double       new_res_y,
                         double       new_res_z,
                         int          new_bmp_size)
{
// Added by JOKim  25NOV2009  -----------------------------------------------------
   distance = new_distance;
   dim_portal.x    = new_dim_x;
   dim_portal.y    = new_dim_y;
   dim_portal.z    = new_dim_z;          // Set z dimension to 1
   voxel_size_portal.x    = new_res_x;
   voxel_size_portal.y    = new_res_y;
   voxel_size_portal.z    = new_res_z;  // Set z Resolution to voxel_size_portal.x
   bmp_size = new_bmp_size;

   batch_dose_portal = new array_3d<double>(dim_portal, ZERO);
   beam_dose_portal  = new array_3d<float>(dim_portal, ZERO);
   beam_error_portal = new array_3d<float>(dim_portal, ZERO);
//   density_portal    = new array_3d<float>(dim_portal, ONE);

// End of Added ------------------------------------------------------------------
//   portal_exit();
}

// destructor
portal_dose::~portal_dose(void)
{
   portal_exit();
}

// clear batch dose matrix
void portal_dose::reset_batch_dose_portal(int n_repeat)
{
// Added by JOKim 25NOV2009 --------------------------------------------------
   // portal_exit();
   n_repeat_portal = n_repeat;
   // Reset Portal Imager after Recording
   memset((void *) batch_dose_portal->data, 0, batch_dose_portal->size*sizeof(double));
// End of Added --------------------------------------------------------------
}

// move the portal dose plane to the correct place using the
// coordinates of the target and the beam direction angles
void portal_dose::move(real_3 origin, real gantry_angle, real table_angle)
{

   // Beam origin after rotated from beam_modifier
   r_origin.x = origin.x;
   r_origin.y = origin.y;
   r_origin.z = origin.z;

   // Gantry and Table angles
   portal_gantry_angle = gantry_angle;
   portal_table_angle  = table_angle;

   // Portal Imager Distance in Beam Coordinetes
   portal_distance.x = 0.0;  // No Shift
   portal_distance.y = 0.0;  // No Shift
   portal_distance.z = 100.0 + distance;

   // Portal Imager Size
   cube_size_portal.x = voxel_size_portal.x*dim_portal.x;
   cube_size_portal.y = voxel_size_portal.y*dim_portal.y;
   cube_size_portal.z = voxel_size_portal.z*dim_portal.z;

   // Portal Imager Origin Position in Beam Coordinates (Corner)
   portal_center.x = cube_size_portal.x * 0.5;
   portal_center.y = cube_size_portal.y * 0.5;
   portal_center.z = 0.0;

   // origin changed to imager coordinates
   beam_origin.x = portal_center.x - portal_distance.x;
   beam_origin.y = portal_center.y - portal_distance.y;
   beam_origin.z = portal_center.z - portal_distance.z;

   cout << "XVMC> Portal Dose Imager Setup: " << endl;
   cout << "      Distance = " << distance << " (cm) from Isocenter; "
        << portal_distance.z << " (cm) from Origin" << endl;
   cout << "      Diemnsion = " << dim_portal.x << " (pixels) x "
				<< dim_portal.y << " (pixels) x "
				<< dim_portal.z << " (pixels)"<< endl;
   cout << "                  " << cube_size_portal.x << " (cm) x "
                                << cube_size_portal.y << " (cm) x "
                                << cube_size_portal.z << " (cm)"<< endl;
   cout << "      resolution = " << voxel_size_portal.x << " (cm) x "
                                 << voxel_size_portal.y << " (cm) x "
                                 << voxel_size_portal.z << " (cm)"<< endl;
   cout << "      Bitmap Size = " << bmp_size << endl;
   cout << "      Beam Origin = " << beam_origin.x << ", "
                                  << beam_origin.y << ", "
                                  << beam_origin.z << " (in Imager Coordinates)" << endl;
   // portal_exit();
}

// update dose and dose error for one batch
void portal_dose::update_batch(double dose_factor, double batch_weight)
{
// Added by JOKim 25NOV2009  -----------------------------------------------------------
   // portal_exit();
   real dose = 0.0;
   real voxel_vol_portal = voxel_size_portal.x * voxel_size_portal.y * voxel_size_portal.z;
   dose_factor /= voxel_vol_portal;
   for (register int k=0; k<dim_portal.z; ++k) {
   	for (register int j=0; j<dim_portal.y; ++j) {
      	for (register int i=0; i<dim_portal.x; ++i) {
         	dose = batch_dose_portal->matrix[i][j][k]*dose_factor/batch_weight;
         	beam_dose_portal->matrix[i][j][k]  += dose;
         	beam_error_portal->matrix[i][j][k] += dose*dose;
      	}
   	}
   }

// End of Added ------------------------------------------------------------------------
}

// Added by JOKim 25NOV2009  -----------------------------------------------------------
// add portal dose contribution of one photon
void portal_dose::add(particle_parameters p, ranmar rndm)
{
   // parameters of particle
   particle_parameters pt;     // Particle in Portal Imager Coordinates
//   pt.type   = p.type;
//   pt.weight = p.weight;
//   pt.energy = p.energy;
//   pt.dir    = p.dir;
   pt = p;

  // std::cout << pt.pos.x << " " << pt.pos.y << " " << pt.pos.z << std::endl;

   // Patient Coordinates change back to Beam coordinates pt.pos
   // vectors, p.pos and origin, are in the patient coordinates yet
   pt.pos.x = p.pos.x - r_origin.x;
   pt.pos.y = p.pos.y - r_origin.y;
   pt.pos.z = p.pos.z - r_origin.z;

   // distance from particle to rotated beam origin
   real origin_point_dist
          = sqrt(pt.pos.x*pt.pos.x+pt.pos.y*pt.pos.y+pt.pos.z*pt.pos.z);

   // Normalized Vector of Particle wrt Current Origin
   pt.pos.x = pt.pos.x/origin_point_dist;
   pt.pos.y = pt.pos.y/origin_point_dist;
   pt.pos.z = pt.pos.z/origin_point_dist;

   // Rotate back (reverse rotation)
   cos_alpha    = cos(-portal_gantry_angle);
   sin_alpha    = sin(-portal_gantry_angle);
   cos_beta     = cos(-portal_table_angle);
   sin_beta     = sin(-portal_table_angle);

   real temp = 0.0; // For temporary variable

   // rotate position vector by the gantry angle
   temp     = pt.pos.x*sin_alpha;
   pt.pos.x = pt.pos.x*cos_alpha - pt.pos.z*sin_alpha;
   pt.pos.z = temp               + pt.pos.z*cos_alpha;

   // rotate direction vector by the gantry angle
   temp     = pt.dir.x*sin_alpha;
   pt.dir.x = pt.dir.x*cos_alpha - pt.dir.z*sin_alpha;
   pt.dir.z = temp               + pt.dir.z*cos_alpha;

   // rotate position vector by the table angle
   temp     = pt.pos.x*sin_beta;
   pt.pos.x = pt.pos.x*cos_beta + pt.pos.y*sin_beta;
   pt.pos.y = -temp             + pt.pos.y*cos_beta;

   // rotate direction vector by the table angle
   temp     = pt.dir.x*sin_beta;
   pt.dir.x = pt.dir.x*cos_beta + pt.dir.y*sin_beta;
   pt.dir.y = -temp             + pt.dir.y*cos_beta;


   // particle position shifted to imager coordinates
   pt.pos.x += beam_origin.x;
   pt.pos.y += beam_origin.y;
   pt.pos.z += beam_origin.z;

   real_3 escape;

   if ( (pt.pos.x > ZERO) && (pt.pos.x < cube_size_portal.x) &&
        (pt.pos.y > ZERO) && (pt.pos.y < cube_size_portal.y) &&
        (pt.pos.z > ZERO) && (pt.pos.z < cube_size_portal.z) )
   {
      // the particle is within the cube, calculate voxel indices and return
      pt.i.x = int(pt.pos.x/voxel_size_portal.x);
      pt.i.y = int(pt.pos.y/voxel_size_portal.y);
      pt.i.z = int(pt.pos.z/voxel_size_portal.z);
   }
   else
   {
   	// the particle is outside the calculation cube, therefore
   	// find the surface point where this particle hits the cube

   	escape.x = pt.pos.x;
   	escape.y = pt.pos.y;
   	escape.z = pt.pos.z;

   	if (pt.dir.z > ZERO) // check in +z direction
   	{
      	temp  = escape.z/pt.dir.z;
      	pt.pos.x = escape.x - temp*pt.dir.x;
      	pt.pos.y = escape.y - temp*pt.dir.y;
      	if ( (pt.pos.x > ZERO) && (pt.pos.x < cube_size_portal.x) &&
         	  (pt.pos.y > ZERO) && (pt.pos.y < cube_size_portal.y) )
      	{
         	pt.pos.z = ZERO;
         	pt.i.x = int(pt.pos.x/voxel_size_portal.x);
         	pt.i.y = int(pt.pos.y/voxel_size_portal.y);
         	pt.i.z = 0;
      	}
   	}

      else if (pt.dir.z < ZERO) // check in -z direction
   	{
      	temp  = (cube_size_portal.z-escape.z)/pt.dir.z;
      	pt.pos.x = escape.x + temp*pt.dir.x;
      	pt.pos.y = escape.y + temp*pt.dir.y;
      	if ( (pt.pos.x > ZERO) && (pt.pos.x < cube_size_portal.x) &&
         	  (pt.pos.y > ZERO) && (pt.pos.y < cube_size_portal.y) )
      	{
         	pt.pos.z = cube_size_portal.z;
         	pt.i.x = int(pt.pos.x/voxel_size_portal.x);
         	pt.i.y = int(pt.pos.y/voxel_size_portal.y);
         	pt.i.z = dim_portal.z-1;
      	}
   	}

   	// the particle does not go trough the upper/lower side
   	else if (pt.dir.x > ZERO) // check in +x direction
   	{
      	temp  = escape.x/pt.dir.x;
      	pt.pos.y = escape.y - temp*pt.dir.y;
      	pt.pos.z = escape.z - temp*pt.dir.z;
      	if ( (pt.pos.y > ZERO) && (pt.pos.y < cube_size_portal.y) &&
         	  (pt.pos.z > ZERO) && (pt.pos.z < cube_size_portal.z) )
      	{
         	pt.pos.x = ZERO;
         	pt.i.x = 0;
         	pt.i.y = int(pt.pos.y/voxel_size_portal.y);
         	pt.i.z = int(pt.pos.z/voxel_size_portal.z);
      	}
   	}

      else if (pt.dir.x < ZERO)  // check in -x direction
   	{
      	temp  = (cube_size_portal.x-escape.x)/pt.dir.x;
      	pt.pos.y = escape.y + temp*pt.dir.y;
      	pt.pos.z = escape.z + temp*pt.dir.z;
      	if ( (pt.pos.y > ZERO) && (pt.pos.y < cube_size_portal.y) &&
         	  (pt.pos.z > ZERO) && (pt.pos.z < cube_size_portal.z) )
      	{
         	pt.pos.x = cube_size_portal.x;
         	pt.i.x = dim_portal.x-1;
         	pt.i.y = int(pt.pos.y/voxel_size_portal.y);
         	pt.i.z = int(pt.pos.z/voxel_size_portal.z);
      	}
   	}

   	// the particle does not go trough the left/right side
   	else if (pt.dir.y > ZERO)  // check in +y direction
   	{
      	temp  = escape.y/pt.dir.y;
      	pt.pos.x = escape.x - temp*pt.dir.x;
      	pt.pos.z = escape.z - temp*pt.dir.z;
      	if ( (pt.pos.x > ZERO) && (pt.pos.x < cube_size_portal.x) &&
         	  (pt.pos.z > ZERO) && (pt.pos.z < cube_size_portal.z) )
      	{
         	pt.pos.y = ZERO;
         	pt.i.x = int(pt.pos.x/voxel_size_portal.x);
         	pt.i.y = 0;
         	pt.i.z = int(pt.pos.z/voxel_size_portal.z);
      	}
   	}

      else if (pt.dir.y < ZERO)  // check in -y direction
   	{
      	temp  = (cube_size_portal.y-escape.y)/pt.dir.y;
      	pt.pos.x = escape.x + temp*pt.dir.x;
      	pt.pos.z = escape.z + temp*pt.dir.z;
      	if ( (pt.pos.x > ZERO) && (pt.pos.x < cube_size_portal.x) &&
         	  (pt.pos.z > ZERO) && (pt.pos.z < cube_size_portal.z) )
      	{
         	pt.pos.y = cube_size_portal.y;
         	pt.i.x = int(pt.pos.x/voxel_size_portal.x);
         	pt.i.y = dim_portal.y-1;
         	pt.i.z = int(pt.pos.z/voxel_size_portal.z);
      	}
   	}
      else
      {
         return;
      }
   } // End of else

   //if (pt.pos.z < 0) return;
   //std::cout << pt.pos.x << " " << pt.pos.y << " " << pt.pos.z << "   "
   //     << pt.i.x   << " " << pt.i.y   << " " << pt.i.z   << std::endl;

#ifdef USE_SOBOL
   // store Sobol random numbers to call "multi_photon" recursively
   real   zeta[SOBOL_DIM];
#else
   real   zeta[3];        // array to store random numbers
   for (register int i=0; i<3; ++i) zeta[i] = rndm.number();
#endif // USE_SOBOL

   if (pt.type == PHOTON) {
      if (pt.energy <= k0_cut) {
         kerma_photon_portal(pt, rndm,
#ifdef CHECK_ENERGY
                   sum_energy,
#endif // CHECK_ENERGY
                   batch_dose_portal);
      }
      else {
         int n_repeat = n_repeat_portal;
         multi_photon_portal(pt, n_repeat, rndm,
#ifdef USE_SOBOL
                               zeta,
#endif // USE_SOBOL
#ifdef CHECK_ENERGY
                               sum_energy,
#endif // CHECK_ENERGY
                               batch_dose_portal);
      }
	}
	else {
      // this is a charged particle
      one_electron_portal(pt, rndm,
#ifdef CHECK_ENERGY
                   sum_energy,
#endif // CHECK_ENERGY
                   batch_dose_portal);
	}

   return;
   // portal_exit();
   // return(false);
}
// End of Addes  --------------------------------------------------------------------------


