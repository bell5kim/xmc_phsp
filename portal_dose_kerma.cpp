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
#include "xvmc_util.h"
#include "photon_data.h"
#include <math.h>
#include "definitions.h"
#include "global.h"
#include "compton_XS_diff.h"
#include "pair_XS_diff.h"
#include "electron_data.h"
#include "moller_XS.h"
#include "bhabha_XS.h"
#include "brems_XS.h"

#define ELECTRON_TRACK_FULL_STOP 1
// End of Added -------------------------
#include "portal_dose.h"

// ****************************************
// declare global variables
// ****************************************
// Added by JOKim 25NOV2009 -------------

extern photon_data   *p_h2o;    // water photon transport data for this beam
extern compton_XS_diff  compton; // differential Compton cross section data
extern pair_XS_diff     dif_pair; // differential pair cross section data

// water electron transport data for this beam
extern electron_data *e_h2o;
// Moller cross section for water
extern moller_XS h2o_moller;
// Bhabha cross section for water
extern bhabha_XS h2o_bhabha;
// bremsstrahlung cross section for water
extern brems_XS h2o_brems;

// sample multiple scattering angle
bool mscat(const real &, const real &, real &, real &, ranmar &);

real rho = 1.0;
real DENS_SCOL = 1.0;
real DENS_SRAD = 1.0/DENS_SCOL;
real DENS_CCOL = 0.0;
real DENS_FCHI = 1.0;
real DENS_COMP = 0.99*rho+0.01*rho*rho;
real DENS_PAIR = 1.0;
real DENS_PHOT = 1.0;

void   swap_bytes( float &); // swap bytes
// End of Added -------------------------
// ****************************************
// member functions of class portal_dose
// ****************************************

// define portal dose image plane by the distance to the target
// the image matrix dimension, the image matrix resolution and
// the BMP image size
portal_dose::portal_dose(double       new_distance,
                         unsigned int new_dim_x, unsigned int new_dim_y,
                         double       new_res_x, double       new_res_y,
                         int          new_bmp_size)
{
// Added by JOKim  25NOV2009  -----------------------------------------------------
   distance = new_distance;
   p_dim.x    = new_dim_x;
   p_dim.y    = new_dim_y;
   p_dim.z    = 1;          // Set z dimension to 1
   p_size.x    = new_res_x;
   p_size.y    = new_res_y;
   p_size.z    = new_res_x;  // Set z Resolution to p_size.x
   bmp_size = new_bmp_size;

   batch_portal_dose = new array_2d<float>(p_dim.x, p_dim.y, ZERO);
   beam_portal_dose  = new array_2d<float>(p_dim.x, p_dim.y, ZERO);
   beam_portal_error = new array_2d<float>(p_dim.x, p_dim.y, ZERO);

   // Portal Imager Center Position in Beam Coordinates
   p_pos.x = 0.0;
   p_pos.y = 0.0;
   p_pos.z = distance + 100.0;

   // Define 4 Points of Portal Imager in Beam Coordinates
   x1.x = p_pos.x - p_dim.x * p_size.x * 0.5;
   x1.y = p_pos.y;
   x1.z = p_pos.z;

   x2.x = p_pos.x + p_dim.x * p_size.x * 0.5;
   x2.y = p_pos.y;
   x2.z = p_pos.z;

   y1.x = p_pos.x;
   y1.y = p_pos.y - p_dim.y * p_size.y * 0.5;
   y1.z = p_pos.z;

   y2.x = p_pos.x;
   y2.y = p_pos.y + p_dim.y * p_size.y * 0.5;
   y2.z = p_pos.z;

   cout << "XVMC> Portal Dose Imager Setup: Distance = "
        << distance << " dim = " << p_dim.x << " x " << p_dim.y
	<< "  resolution = " << p_size.x << " x " << p_size.y
	<< "  Bitmap Size = " << bmp_size << endl;

// End of Added ------------------------------------------------------------------
//   portal_exit();
}

// destructor
portal_dose::~portal_dose(void)
{
   portal_exit();
}

// clear batch dose matrix
void portal_dose::reset_batch_dose_portal(void)
{
// Added by JOKim 25NOV2009 --------------------------------------------------
   // portal_exit();
	// Reset Portal Imager after Recording
   memset((void *) batch_portal_dose->data, 0, batch_portal_dose->size*sizeof(float));
// End of Added --------------------------------------------------------------
}

// move the portal dose plane to the correct place using the
// coordinates of the target and the beam direction angles
void portal_dose::move(real_3 origin, real gantry_angle, real table_angle)
{
   cos_alpha    = cos(gantry_angle);
   sin_alpha    = sin(gantry_angle);
   cos_beta     = cos(table_angle);
   sin_beta     = sin(table_angle);
   if(cos_alpha < 1.745329252e-6) cos_alpha = 1.745329252e-6;
   if(sin_alpha < 1.745329252e-6) sin_alpha = 1.745329252e-6;

   // normalize position vector
   real origin_point_dist = distance + 100.0;
   p_dir.x  = p_pos.x/origin_point_dist;
   p_dir.y  = p_pos.y/origin_point_dist;
   p_dir.z  = p_pos.z/origin_point_dist;

   real temp = 0.0; // Temporary Variable
   // rotate position vector by the gantry angle
   temp    = p_dir.x*sin_alpha;
   p_dir.x = p_dir.x*cos_alpha - p_dir.z*sin_alpha;
   p_dir.z = temp              + p_dir.z*cos_alpha;

   // rotate position vector by the table angle
   temp    = p_dir.x*sin_beta;
   p_dir.x = p_dir.x*cos_beta + p_dir.y*sin_beta;
   p_dir.y = -temp            + p_dir.y*cos_beta;

   // calculate the accelerator escaping position of this particle
   p_pos.x = origin.x + origin_point_dist*p_dir.x;
   p_pos.y = origin.y + origin_point_dist*p_dir.y;
   p_pos.z = origin.z + origin_point_dist*p_dir.z;

   real_3 dummy;

   // X1 in Patient Model Coordinates ----------------------------------
   dummy.x = x1.x;  dummy.y = x1.y; dummy.z = x1.z;

   real dist = sqrt(dummy.x*dummy.x + dummy.y*dummy.y + dummy.z*dummy.z);

   p_dir.x  = dummy.x/dist;
   p_dir.y  = dummy.y/dist;
   p_dir.z  = dummy.z/dist;

   // rotate position vector by the gantry angle
   temp    = p_dir.x*sin_alpha;
   p_dir.x = p_dir.x*cos_alpha - p_dir.z*sin_alpha;
   p_dir.z = temp              + p_dir.z*cos_alpha;

   // rotate position vector by the table angle
   temp    = p_dir.x*sin_beta;
   p_dir.x = p_dir.x*cos_beta + p_dir.y*sin_beta;
   p_dir.y = -temp            + p_dir.y*cos_beta;

   // calculate the accelerator escaping position of this particle
   x1.x = origin.x + dist*p_dir.x;
   x1.y = origin.y + dist*p_dir.y;
   x1.z = origin.z + dist*p_dir.z;

   // X2  in Patient Model Coordinates ----------------------------------
   dummy.x = x2.x;  dummy.y = x2.y; dummy.z = x2.z;

   dist = sqrt(dummy.x*dummy.x + dummy.y*dummy.y + dummy.z*dummy.z);

   p_dir.x  = dummy.x/dist;
   p_dir.y  = dummy.y/dist;
   p_dir.z  = dummy.z/dist;

   // rotate position vector by the gantry angle
   temp    = p_dir.x*sin_alpha;
   p_dir.x = p_dir.x*cos_alpha - p_dir.z*sin_alpha;
   p_dir.z = temp              + p_dir.z*cos_alpha;

   // rotate position vector by the table angle
   temp    = p_dir.x*sin_beta;
   p_dir.x = p_dir.x*cos_beta + p_dir.y*sin_beta;
   p_dir.y = -temp            + p_dir.y*cos_beta;

   // calculate the accelerator escaping position of this particle
   x2.x = origin.x + dist*p_dir.x;
   x2.y = origin.y + dist*p_dir.y;
   x2.z = origin.z + dist*p_dir.z;

   // Y1  in Patient Model Coordinates ----------------------------------
   dummy.x = y1.x;  dummy.y = y1.y; dummy.z = y1.z;

   dist = sqrt(dummy.x*dummy.x + dummy.y*dummy.y + dummy.z*dummy.z);

   p_dir.x  = dummy.x/dist;
   p_dir.y  = dummy.y/dist;
   p_dir.z  = dummy.z/dist;

   // rotate position vector by the gantry angle
   temp    = p_dir.x*sin_alpha;
   p_dir.x = p_dir.x*cos_alpha - p_dir.z*sin_alpha;
   p_dir.z = temp              + p_dir.z*cos_alpha;

   // rotate position vector by the table angle
   temp    = p_dir.x*sin_beta;
   p_dir.x = p_dir.x*cos_beta + p_dir.y*sin_beta;
   p_dir.y = -temp            + p_dir.y*cos_beta;

   // calculate the accelerator escaping position of this particle
   y1.x = origin.x + dist*p_dir.x;
   y1.y = origin.y + dist*p_dir.y;
   y1.z = origin.z + dist*p_dir.z;

   // Y2  in Patient Model Coordinates ----------------------------------
   dummy.x = y2.x;  dummy.y = y2.y; dummy.z = y2.z;

   dist = sqrt(dummy.x*dummy.x + dummy.y*dummy.y + dummy.z*dummy.z);

   p_dir.x  = dummy.x/dist;
   p_dir.y  = dummy.y/dist;
   p_dir.z  = dummy.z/dist;

   // rotate position vector by the gantry angle
   temp    = p_dir.x*sin_alpha;
   p_dir.x = p_dir.x*cos_alpha - p_dir.z*sin_alpha;
   p_dir.z = temp              + p_dir.z*cos_alpha;

   // rotate position vector by the table angle
   temp    = p_dir.x*sin_beta;
   p_dir.x = p_dir.x*cos_beta + p_dir.y*sin_beta;
   p_dir.y = -temp            + p_dir.y*cos_beta;

   // calculate the accelerator escaping position of this particle
   y2.x = origin.x + dist*p_dir.x;
   y2.y = origin.y + dist*p_dir.y;
   y2.z = origin.z + dist*p_dir.z;

	cout << "XVMC> Position of Portal Imager Edges after Rotation" << endl;
   cout << "      x1 = " << x1.x << "  " << x1.y << "  " << x1.z << endl;
   cout << "      x2 = " << x2.x << "  " << x2.y << "  " << x2.z << endl;
   cout << "      y1 = " << y1.x << "  " << y1.y << "  " << y1.z << endl;
   cout << "      y2 = " << y2.x << "  " << y2.y << "  " << y2.z << endl;

   // Portal Imager Limits in Patient Coordinates
   min_limit.x  = x1.x;
   max_limit.x  = x2.x;
   min_limit.y  = y1.y;
   max_limit.y  = y2.y;
   min_limit.z  = x1.z;
   max_limit.z  = x2.z;

   // New Resolution of Pixel After Rotation
   p_size_rot.x = fabs(max_limit.x - min_limit.x)/p_dim.x;
   p_size_rot.y = fabs(max_limit.y - min_limit.y)/p_dim.y;
   p_size_rot.z = fabs(max_limit.z - min_limit.z)/p_dim.x;

   // Normal Vector of Portal Imager in Patient Coordinates
   normal.x = (x2.y - x1.y)*(y1.z - x1.z) - (y1.y - x1.y)*(x2.z - x1.z);
   normal.y = (x2.z - x1.z)*(y1.x - x1.x) - (y1.z - x1.z)*(x2.x - x1.x);
   normal.z = (x2.x - x1.x)*(y1.y - x1.y) - (y1.x - x1.x)*(x2.y - x1.y);

   // distance to the origin point before normalization
   dist2zero = normal.x*x1.x + normal.y*x1.y + normal.z*x1.z;

   // calculate normalization factor
   real norm = sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);

   if (norm > 1.0e-5) {
      normal.x  /= norm;
      normal.y  /= norm;
      normal.z  /= norm;
      dist2zero  /= norm;
   }
   else{
      xvmc_error("portal_dose::move",
                 "the 3 points are improperly defined",8);
   }

	// by definition dist2zero must be larger than or equal to zero
   if (dist2zero < ZERO) {
      normal.x  *= -ONE;
      normal.y  *= -ONE;
      normal.z  *= -ONE;
      dist2zero *= -ONE;
   }

   cout << "XVMC> Source Position in Patient Model Coordinates:" << origin.x << ", " << origin.y << ", " << origin.z << endl;
   cout << "XVMC> Portal Imager Limits in Parient Model Coordinates" << endl;
   cout << "XVMC>      x = " << p_pos.x << " (" << min_limit.x << ":" << max_limit.x << ") by " << p_size_rot.x << endl;
   cout << "XVMC       y = " << p_pos.y << " (" << min_limit.y << ":" << max_limit.y << ") by " << p_size_rot.y << endl;
   cout << "XVMC       z = " << p_pos.z << " (" << min_limit.z << ":" << max_limit.z << ") by " << p_size_rot.z << endl;
   cout << "XVMC> Portal Imager Normal Vector = " << normal.x << ", " << normal.y << ", " << normal.z << endl;
   // cout << "XVMC> Beam Vector = " << beam_dir.x << ", " << beam_dir.y << ", " << beam_dir.z << endl;
   cout << "XVMC> Distance of Plane to Patient Model Origin: " << dist2zero << endl;
   // portal_exit();
}

// update dose and dose error for one batch
void portal_dose::update_batch(double batch_weight)
{
// Added by JOKim 25NOV2009  -----------------------------------------------------------
   // portal_exit();
   for (register int j=0; j<p_dim.y; ++j) {
      for (register int i=0; i<p_dim.x; ++i) {
         real dose = batch_portal_dose->matrix[i][j]/batch_weight;
         beam_portal_dose->matrix[i][j]  += dose;
         beam_portal_error->matrix[i][j] += dose*dose;
      }
   }

// End of Added ------------------------------------------------------------------------
}

// write portal dose files
void portal_dose::write_files(const char *hed_name,
                              const char *dat_name,
                              const char *err_name,
                              const char *bmp_name, int n_batch)
{
   // portal_exit();
// Added by JOKim 25NOV2009 ---------------------------------------------------------------
   cout << hed_name << "  " << dat_name << "  "  << err_name << "  " << bmp_name << endl;
   cout << "XVMC> Portal Dose Write: " << dat_name << endl;
   // Write Header File
   ofstream hed_file(hed_name,ios::out);
   if (!hed_file){
      xvmc_error("read_inp_file","cannot open density header file",8);
   }
   hed_file.setf( ios::fixed, ios::floatfield );

   hed_file << setprecision(6);
   hed_file << "PIXELSIZE      |" << setw(11) << p_size.x
                                  << setw(11) << p_size.y << endl;
   hed_file << "DIMENSION      |" << setw(11) << p_dim.x
                                  << setw(11) << p_dim.y << endl;
   hed_file << "DISTANCE       |" << setw(11) << distance << endl;

   // Write Portal Dose Matrix
   float *out_data = NULL;     // pointer to portal dose output data
   float *err_data = NULL;     // pointer to portal dose error output data
   float maxDose = -1.0;
   int   iMax = 0;
   int   jMax = 0;

   // allocate memory for dose output data
   if ( (out_data = new float[p_dim.x*p_dim.y]) == NULL ){
      xvmc_error("write_dose_file",
                 "cannot allocate memory for portal dose output data",8);
   }

   // allocate memory for dose error output data
   if ( (err_data = new float[p_dim.x*p_dim.y]) == NULL ){
      xvmc_error("write_dose_file",
                 "cannot allocate memory for dose error output data",8);
   }

   // write portal dose and error data to output memory
   float dose,error;
   float *out_count = out_data;
   float *err_count = err_data;

   real dose_factor = 1.602 / (p_size.x * p_size.y) / n_batch;
   dose_factor = 1.0;
   for (register int j=0; j<p_dim.y; ++j) {
      for (register int i=0; i<p_dim.x; ++i) {
         dose  = beam_portal_dose->matrix[i][j];
         // convert average dose squared to standard deviation (error)
         error = float(n_batch)*beam_portal_error->matrix[i][j];
         error = sqrt(fabs(error-dose*dose)/float(n_batch-1));
         dose  *= dose_factor;
	 error *= dose_factor;

	 if (dose > maxDose) {
	    maxDose = dose;
	    iMax = i;
	    jMax = j;
	 }
         swap_bytes(dose);      swap_bytes(error);
         *out_count = dose;     *err_count = error;
         ++out_count;           ++err_count;
      }
   }

   hed_file << "DOSE-MAXIMUM   |" << setw(11) << maxDose << endl;
   hed_file << "POS-MAXIMUM    |" << setw(11) << iMax
                                  << setw(11) << jMax << endl;
   hed_file << "END-INPUT      |" << endl;
   hed_file.close();

   ofstream dat_file(dat_name,ios::out|ios::binary);
   char *out_dat_memory = (char *) out_data;
   dat_file.write( out_dat_memory, p_dim.x*p_dim.y*sizeof(float));
   dat_file.close();

   ofstream err_file(err_name,ios::out|ios::binary);
   char *out_err_memory = (char *) err_data;
   err_file.write( out_err_memory, p_dim.x*p_dim.y*sizeof(float));
   err_file.close();

   // Write BMP
   // write_bmp(bmp_name, p_dim.x, p_dim.y, p_dim.x, p_dim.y, out_dat_memory);
   // free memory
   delete [] out_data; out_data = NULL;
   delete [] err_data; err_data = NULL;

// End of Added ---------------------------------------------------------------------------
}

// Added by JOKim 25NOV2009  -----------------------------------------------------------
// add portal dose contribution of one photon
void portal_dose::add(particle_parameters p, ranmar rndm)
{
   // parameters of particle
   particle_parameters pt;     // Particle in Portal Imager Coordinates
   pt.type = p.type;
   pt.weight = p.weight;
   pt.energy = p.energy;

   // scalar products
   real product_pos = p.pos.x*normal.x + p.pos.y*normal.y + p.pos.z*normal.z;
   real product_dir = p.dir.x*normal.x + p.dir.y*normal.y + p.dir.z*normal.z;

   if (product_dir != ZERO) {
       real R = (dist2zero-product_pos)/product_dir;
       // if (R > 0) {
       // New Particle Position in Portal Image Coordinates
       pt.pos.x = p.pos.x + R*p.dir.x - min_limit.x;
       pt.pos.y = p.pos.y + R*p.dir.y - min_limit.y;
       pt.pos.z = p.pos.z + R*p.dir.z - min_limit.z;

       pt.pos.x /= cos_alpha;
       if (fabs(p_size_rot.z) > fabs(p_size_rot.x)){
	  pt.pos.x = pt.pos.z / sin_alpha;
       }
       pt.pos.z = 0.0;  // Reset Z Position at Portal Imager Surface
       //if ((pt.pos.x >= 0.0 && pt.pos.x <= p_dim.x * p_size.x) &&
       //   (pt.pos.y >= 0.0 && pt.pos.y <= p_dim.y * p_size.y) &&
       //   (pt.pos.z >= 0.0 && pt.pos.z <= 1.0     * p_size.z)) {
       pt.i.x = (int)(pt.pos.x/p_size.x);
       pt.i.y = (int)(pt.pos.y/p_size.y);
       pt.i.z = 0;

       if (pt.i.x >= 0 && pt.i.x < p_dim.x &&
           pt.i.y >= 0 && pt.i.y < p_dim.y) {
	  // Direction Vector of Particle in Portal Imager Coodinates
	  real norm = sqrt(p.dir.x*p.dir.x+p.dir.y*p.dir.y+p.dir.z*p.dir.z);
	  if(norm < 1.745329252e-6) norm = 1.745329252e-6;
	  pt.dir.x = (p.dir.x*normal.x)/norm;
	  pt.dir.y = (p.dir.y*normal.y)/norm;
	  pt.dir.z = (p.dir.z*normal.z)/norm;
          if (pt.dir.z < 1.745329252e-6) return;

	  if (pt.type == PHOTON) {
              if (pt.energy <= k0_cut) {
		 portal_kerma(pt, rndm);
              }
              else {
                 portal_photon(pt, rndm);
              }
	  }
	  else {
            pt.type == ELECTRON;
            one_electron(pt, rndm);
	    // real e_energy = pt.energy - 0.511034;  // Kinetic Energy
	    // batch_portal_dose->matrix[pt.i.x][pt.i.y]  += (e_energy * pt.weight);
	  }
       }
	//}
	     // }
   }

   return;
   // portal_exit();
   // return(false);
}
// End of Addes  --------------------------------------------------------------------------

// Added by JOK 12-25-2009 ---------------------------------------------------------------
// Portal Kerma Transport
void portal_dose::portal_kerma(particle_parameters &p, ranmar &rndm) {

   int    i_e;            // energy index

   real_3 step;           // step sizes to the next x-, y-, z-voxel boundaries
   int_3  istep;          // step direction to the next voxel boundary (-1,0,1)
   int_3  inew;           // new voxel indices

   real   mu_tot_h2o;     // total cross section (attenuation coeff.) for water
   real   mu_en_h2o;      // energy absorption coefficient for water
   real   p_phot_h2o;     // probability for photo electric absorption in water
   real   p_pair_h2o;     // pair production probability in water
   real   p_comp_h2o;     // Compton interaction probability in water
   real   mu_tot;         // total cross section in present voxel
   real   mu_en;          // energy absorption coefficient in present voxel
   real   p_phot;         // probability for photo absorption in present voxel
   real   p_pair;         // pair production probability in present voxel
   real   p_comp;         // Compton interaction probability in present voxel
   real   p_tot;          // sum of p_phot, p_pair, p_comp for normalization
   real   dens_corr_phot; // density correction for photo cross section
   real   dens_corr_pair; // density correction for pair cross section
   real   dens_corr_comp; // density correction for Compton cross section
   real   zeta;           // random number in [0,1]
   real   num_mfp;        // number of mean free paths
   real   path;           // photon path to the next interaction
   real   path_step;      // photon path length to the present position
   real   voxel_step;     // step size to the next voxel boundary
   real   cos_t_x;        // Compton photon scattering angle
   real   sin_t_x;        // Compton photon scattering angle
   real   energy_e;       // energy of the Compton electron
   real   cos_t_e;        // Compton electron scattering angle
   real   sin_t_e;        // Compton electron scattering angle
   real   phi;            // azimuthal angle
   real   cos_phi;        // cos(phi)
   real   sin_phi;        // sin(phi)
   real   rnno;           // auxiliary random number

   // parameters of the Compton photon
   particle_parameters x; x.type = PHOTON;

   // deposited energy times particle weight
   real   energy_depo, tot_energy_depo;

   bool   repeat_step;
   // determine total cross sections and interaction probabilities for water
   i_e=int((p.energy-p_cut)*p_h2o->inverse_delta_e);
   mu_tot_h2o = p_h2o->mu_tot[i_e];
   mu_en_h2o  = p_h2o->mu_en[i_e];
   p_phot_h2o = p_h2o->mu_phot[i_e]/mu_tot_h2o;
   p_pair_h2o = p_h2o->mu_pair[i_e]/mu_tot_h2o;
   if (p.energy <= TWOxEMASS) p_pair_h2o = ZERO;
   p_comp_h2o = ONE - p_pair_h2o - p_phot_h2o;

   // determine number of mean free paths to the next interaction
   // avoid NaN bug
   zeta                   = rndm.number();
   if (zeta == ZERO) zeta = rndm.number();
   num_mfp = -log(ONE-zeta);

   // find the X-, Y-, Z-distance to the next voxel boundary
   // check X-direction ---------------------
   if (p.dir.x > ZERO) {
      istep.x = 1;
      step.x  = (p_size.x*double(p.i.x+1)-p.pos.x)/p.dir.x;
   }
   else {
      if (p.dir.x < ZERO) {
         istep.x = -1;
         step.x  = (p_size.x*double(p.i.x)-p.pos.x)/p.dir.x;
      }
      else {
         istep.x = 0;
         step.x  = HUGE_STEP;
      }
   }

   // check Y-direction ---------------------
   if (p.dir.y > ZERO) {
      istep.y = 1;
      step.y  = (p_size.y*double(p.i.y+1)-p.pos.y)/p.dir.y;
   }
   else {
      if (p.dir.y < ZERO) {
         istep.y = -1;
         step.y  = (p_size.y*double(p.i.y)-p.pos.y)/p.dir.y;
      }
      else {
         istep.y = 0;
         step.y  = HUGE_STEP;
      }
   }

   // check Z-direction ---------------------
   if (p.dir.z > ZERO) {
      istep.z = 1;
      step.z  = (p_size.z*double(p.i.z+1)-p.pos.z)/p.dir.z;
   }
   else {
      if (p.dir.z < ZERO) {
         istep.z = -1;
         step.z  = (p_size.z*double(p.i.z)-p.pos.z)/p.dir.z;
      }
      else {
         istep.z = 0;
         step.z  = HUGE_STEP;
      }
   }

   // update voxel indices
   inew = p.i;

   // start transport
   repeat_step = true;
   path_step   = ZERO;

   // store deposited energy
   tot_energy_depo = ZERO;

   while (repeat_step) {
   // the photon will be traced to the next voxel boundary, if the
   // photon reaches the boundary we repeat this step in the next
   // voxel, the bool variable "repeat_step" is "true" in this case

      // determine the interaction probabilities in the present voxel
      dens_corr_phot = 1.0;  // No correction for water
      dens_corr_pair = 1.0;  // No correction for water
      dens_corr_comp = 1.0;  // No correction for water
      p_phot = p_phot_h2o*dens_corr_phot;
      p_pair = p_pair_h2o*dens_corr_pair;
      p_comp = p_comp_h2o*dens_corr_comp;
      p_tot  = p_phot + p_pair + p_comp;

      // normalize (only "p_comp", we don't need "p_phot" and "p_pair")
      p_comp /= p_tot;

      // total cross section in this voxel
      mu_tot = mu_tot_h2o*p_tot;

      // path to the next interaction
      path = num_mfp/mu_tot;

      // find the next voxel
      if ( (step.z < step.x) && (step.z < step.y) ) {
         voxel_step = step.z;
         if (voxel_step > path) {
            voxel_step  = path;
            step.z     -= voxel_step;
            repeat_step = false;
         }
         else {
            step.z  = p_size.z/fabs(p.dir.z);
            inew.z += istep.z;
         }
         step.x -= voxel_step;
         step.y -= voxel_step;
      }
      else {
         if (step.x < step.y) {
            voxel_step = step.x;
            if (voxel_step > path) {
               voxel_step  = path;
               step.x     -= voxel_step;
               repeat_step = false;
            }
            else {
               step.x  = p_size.x/fabs(p.dir.x);
               inew.x += istep.x;
            }
            step.y -= voxel_step;
            step.z -= voxel_step;
         }
         else {
            voxel_step = step.y;
            if (voxel_step > path) {
               voxel_step  = path;
               step.y     -= voxel_step;
               repeat_step = false;
            }
            else {
               step.y  = p_size.y/fabs(p.dir.y);
               inew.y += istep.y;
            }
            step.x -= voxel_step;
            step.z -= voxel_step;
         }
      }

      // update number of mean free paths and photon path length
      num_mfp   -= voxel_step*mu_tot;
      path_step += voxel_step;

      // calculate deposited energy
      mu_en = mu_en_h2o*dens_corr_comp
            + mu_tot_h2o*p_phot_h2o*p.energy*(dens_corr_phot-dens_corr_comp);
      if (p.energy > TWOxEMASS) {
         mu_en += mu_tot_h2o*p_pair_h2o*(p.energy-TWOxEMASS)*
                     (dens_corr_pair-dens_corr_comp);
      }
      energy_depo = p.weight*mu_en*voxel_step;
      tot_energy_depo                         += energy_depo;
      batch_portal_dose->matrix[p.i.x][p.i.y] += energy_depo;

      if ( (inew.x < 0) || (inew.x >= p_dim.x) ||
           (inew.y < 0) || (inew.y >= p_dim.y) ||
           (inew.z < 0) || (inew.z >= p_dim.z) ){
	 // Out of Portal Imager
         return;
      }

      // update voxel indices
      p.i = inew;

   // if "repeat_step=true" we continue the step in the next voxel
   } // End of while

   if (rndm.number() < p_comp) {
      // Compton scattering
      rnno = rndm.number();
      compton.interaction(p.energy, rnno, rndm,
                          x.energy, cos_t_x, sin_t_x,
                          energy_e, cos_t_e, sin_t_e);

      // azimuthal scattering angle
      phi = TWO_PI * rndm.number();
      sin_phi = sin(phi);
      cos_phi = cos(phi);

      // the weight of the new photon
      x.weight = p.weight;

      // starting position of the new Compton photon
      x.pos.x = p.pos.x + path_step*p.dir.x;
      x.pos.y = p.pos.y + path_step*p.dir.y;
      x.pos.z = p.pos.z + path_step*p.dir.z;

      // rotate direction of the new Compton photon
      x.dir = p.dir;
      rotate(cos_t_x, sin_t_x, cos_phi, sin_phi, x.dir);

      // calculate voxel indices
      x.i.x = int(x.pos.x/p_size.x);
      x.i.y = int(x.pos.y/p_size.y);
      x.i.z = int(x.pos.z/p_size.z);

      // simulate the new photon
      portal_kerma(x, rndm);
   }

   return;
}

void portal_dose::portal_photon (particle_parameters &p,ranmar &rndm) {

   int n_repeat  = 70;
   real   zeta[3];        // array to store random numbers
   for (register int i=0; i<3; ++i) zeta[i] = rndm.number();

   int    i_e;            // energy index

   real_3 step;           // step sizes to the next x-, y-, z-voxel boundaries
   int_3  istep;          // step direction to the next voxel boundary (-1,0,1)
   int_3  inew;           // new voxel indices

   real   mu_tot_h2o;     // total cross section (attenuation coeff.) for water
   real   prob_phot_h2o;  // probability for photo electric absorption in water
   real   prob_pair_h2o;  // pair production probability in water
   real   prob_comp_h2o;  // Compton interaction probability in water
   real   mu_tot;         // total cross section in present voxel
   real   prob_phot;      // probability for photo absorption in present voxel
   real   prob_pair;      // pair production probability in present voxel
   real   prob_comp;      // Compton interaction probability in present voxel
   real   prob_tot;       // sum of prob_phot, prob_pair, prob_comp
   real   eta;            // random number to sample number of mean free paths
   real   delta_eta;      //
   real   eta_prime;      //
   real   num_mfp;        // number of mean free paths
   real   save_mfp;       // save number of mean free paths
   real   path;           // photon path to the next interaction
   real   path_step;      // photon path length to the present position
   real   voxel_step;     // step size to the next voxel boundary
   real   cos_t_cx;       // Compton photon scattering angle
   real   sin_t_cx;       // Compton photon scattering angle
   real   cos_t_ce;       // Compton electron scattering angle
   real   sin_t_ce;       // Compton electron scattering angle
   real   cos_t_pe;       // pair electron scattering angle
   real   sin_t_pe;       // pair electron scattering angle
   real   cos_t_pp;       // pair positron scattering angle
   real   sin_t_pp;       // pair positron scattering angle
   real   phi;            // azimuthal angle
   real   cos_phi;        // cos(phi)
   real   sin_phi;        // sin(phi)

   particle_parameters x_comp; x_comp.type = PHOTON;   // Compton photon
   particle_parameters e_comp; e_comp.type = ELECTRON; // Compton electron
   particle_parameters e_pair; e_pair.type = ELECTRON; // pair electron
   particle_parameters p_pair; p_pair.type = POSITRON; // pair positron
   particle_parameters e_phot; e_phot.type = ELECTRON; // photo electron

   real   x_comp_energy = ZERO;  // save Compton photon energy
   real   e_comp_energy = ZERO;  // save Compton electron energy
   real   e_pair_energy = ZERO;  // save pair electron energy
   real   p_pair_energy = ZERO;  // save pair positron energy
   real   e_phot_energy = ZERO;  // save photo electron energy

   // create and trace electron histories (history repetition)
   MULTI_ELECTRON electron_comp; // Compton electron
   MULTI_ELECTRON electron_pair; // pair electron
   MULTI_ELECTRON positron_pair; // pair positron
   MULTI_ELECTRON electron_phot; // photo electron

   bool next_interaction,repeat_step;

// ****************************************
// begin photon transport
// ****************************************

   // the photon must be within the calculation cube
   if ( (p.i.x < 0) || (p.i.x >= p_dim.x) ||
        (p.i.y < 0) || (p.i.y >= p_dim.y) ||
        (p.i.z < 0) || (p.i.z >= p_dim.z) ) {
      return;
   }

   // kill photon if the energy is below p_cut
   if (p.energy <= p_cut) {
      return;
   }

   // determine total cross sections and interaction probabilities for water
   i_e=int((p.energy-p_cut)*p_h2o->inverse_delta_e);
   mu_tot_h2o = p_h2o->mu_tot[i_e];
   prob_phot_h2o = p_h2o->mu_phot[i_e]/mu_tot_h2o;
   prob_pair_h2o = p_h2o->mu_pair[i_e]/mu_tot_h2o;
   if (p.energy <= TWOxEMASS) prob_pair_h2o = ZERO;
   prob_comp_h2o = ONE - prob_pair_h2o - prob_phot_h2o;

   // random number to sample the number of mean free paths
   const real EPSILON = 1.0e-10;
   eta       = max_of(zeta[0],EPSILON);
   delta_eta = ONE/double(n_repeat);
   eta_prime = (double(n_repeat+1) - eta)*delta_eta;

   // create electron histories in water if the following variables are true
   bool first_comp = true;
   bool first_pair = true;
   bool first_phot = true;

   // count number of history repetitions for portal dose and
   // for checking the energy conservation
   int i_repeat = 0;

   // find the X-, Y-, Z-distance to the next voxel boundary
   // check X-direction
   if (p.dir.x > ZERO) {
      istep.x = 1;
      step.x  = (p_size.x*double(p.i.x+1)-p.pos.x)/p.dir.x;
   }
   else {
      if (p.dir.x < ZERO) {
         istep.x = -1;
         step.x  = (p_size.x*double(p.i.x)-p.pos.x)/p.dir.x;
      }
      else {
         istep.x = 0;
         step.x  = HUGE_STEP;
      }
   }

   // check Y-direction
   if (p.dir.y > ZERO) {
      istep.y = 1;
      step.y  = (p_size.y*double(p.i.y+1)-p.pos.y)/p.dir.y;
   }
   else {
      if (p.dir.y < ZERO) {
         istep.y = -1;
         step.y  = (p_size.y*double(p.i.y)-p.pos.y)/p.dir.y;
      }
      else {
         istep.y = 0;
         step.y  = HUGE_STEP;
      }
   }

   // check Z-direction
   if (p.dir.z > ZERO) {
      istep.z = 1;
      step.z  = (p_size.z*double(p.i.z+1)-p.pos.z)/p.dir.z;
   }
   else
   {
      if (p.dir.z < ZERO) {
         istep.z = -1;
         step.z  = (p_size.z*double(p.i.z)-p.pos.z)/p.dir.z;
      }
      else {
         istep.z = 0;
         step.z  = HUGE_STEP;
      }
   }

   next_interaction = true;
   save_mfp         = ZERO;
   while (next_interaction) {
      // random number to sample path to the next interaction site
      eta_prime -= delta_eta;
      if (eta_prime <= ZERO) break;

      // determine number of mean free paths to the next interaction
      num_mfp   = -log(eta_prime) - save_mfp;
      save_mfp += num_mfp;

      // update voxel indices
      inew = p.i;

      // start transport
      repeat_step = true;
      path_step   = ZERO;

      while (repeat_step) {
      // the photon will be traced to the next voxel boundary, if the
      // photon reaches the boundary we repeat this step in the next
      // voxel, the bool variable "repeat_step" is "true" in this case

         // determine the interaction probabilities in the present voxel
         prob_phot = prob_phot_h2o*1.0;
         prob_pair = prob_pair_h2o*1.0;
         prob_comp = prob_comp_h2o*1.0;
         prob_tot  = prob_phot + prob_pair + prob_comp;

         // total cross section in this voxel
         mu_tot = mu_tot_h2o*prob_tot;

         // path to the next interaction
         path = num_mfp/mu_tot;

         // find the next voxel
         if ( (step.z < step.x) && (step.z < step.y) ) {
            voxel_step = step.z;
            if (voxel_step > path) {
               voxel_step  = path;
               step.z     -= voxel_step;
               repeat_step = false;
            }
            else {
               step.z  = p_size.z/fabs(p.dir.z);
               inew.z += istep.z;
            }
            step.x -= voxel_step;
            step.y -= voxel_step;
         }
         else {
            if (step.x < step.y) {
               voxel_step = step.x;
               if (voxel_step > path) {
                  voxel_step  = path;
                  step.x     -= voxel_step;
                  repeat_step = false;
               }
               else {
                  step.x  = p_size.x/fabs(p.dir.x);
                  inew.x += istep.x;
               }
               step.y -= voxel_step;
               step.z -= voxel_step;
            }
            else {
               voxel_step = step.y;
               if (voxel_step > path) {
                  voxel_step  = path;
                  step.y     -= voxel_step;
                  repeat_step = false;
               }
               else {
                  step.y  = p_size.y/fabs(p.dir.y);
                  inew.y += istep.y;
               }
               step.x -= voxel_step;
               step.z -= voxel_step;
            }
         }

         // update number of mean free paths and photon path length
         num_mfp   -= voxel_step*mu_tot;
         path_step += voxel_step;

         // check voxel index
         if ( (inew.x < 0) || (inew.x >= p_dim.x) ||
              (inew.y < 0) || (inew.y >= p_dim.y) ||
              (inew.z < 0) || (inew.z >= p_dim.z) ) {
            // the photon is leaving the calculation cube, the history
            // ends here, change weight and count portal image dose
            p.weight *= double(n_repeat-i_repeat)/double(n_repeat);
            return;
         }

      // if "repeat_step=true" we continue the step in the next voxel
      }

      // transport the photon
      p.pos.x += path_step*p.dir.x;
      p.pos.y += path_step*p.dir.y;
      p.pos.z += path_step*p.dir.z;

      // the following is in principle not necessary, however, due to
      // the accumulation of truncation errors the actual voxel may be
      // different from the one pointed to by "inew"
      p.i.x = int(p.pos.x/p_size.x);
      p.i.y = int(p.pos.y/p_size.y);
      p.i.z = int(p.pos.z/p_size.z);

      // check voxel index
      if ( (p.i.x < 0) || (p.i.x >= p_dim.x) ||
           (p.i.y < 0) || (p.i.y >= p_dim.y) ||
           (p.i.z < 0) || (p.i.z >= p_dim.z) ) {
         // the photon is leaving the calculation cube, the history
         // ends here, change weight and count portal image dose
         p.weight *= double(n_repeat-i_repeat)/double(n_repeat);

         return;
      }

      // find the X-, Y-, Z-distance to the next voxel boundary
      // check X-direction
      if (p.i.x != inew.x) {
         if (istep.x == 1) {
            step.x  = (p_size.x*double(p.i.x+1)-p.pos.x)/p.dir.x;
         }
         else {
            if (istep.x == -1) {
               step.x  = (p_size.x*double(p.i.x)-p.pos.x)/p.dir.x;
            }
            else {
               step.x  = HUGE_STEP;
            }
         }
      }

      // check Y-direction
      if (p.i.y != inew.y) {
         if (istep.y == 1) {
            step.y  = (p_size.y*double(p.i.y+1)-p.pos.y)/p.dir.y;
         }
         else {
            if (istep.y == -1) {
               step.y  = (p_size.y*double(p.i.y)-p.pos.y)/p.dir.y;
            }
            else {
               step.y  = HUGE_STEP;
            }
         }
      }

      // check Z-direction
      if (p.i.z != inew.z) {
         if (istep.z == 1) {
            step.z  = (p_size.z*double(p.i.z+1)-p.pos.z)/p.dir.z;
         }
         else {
            if (istep.z == -1) {
               step.z  = (p_size.z*double(p.i.z)-p.pos.z)/p.dir.z;
            }
            else {
               step.z  = HUGE_STEP;
            }
         }
      }

      // the photon has reached the interaction site

      // count number of interaction repetitions
      ++i_repeat;

      // normalize "prob_comp" and "prob_pair" (we don't need "prob_phot")
      prob_comp /= prob_tot;
      prob_pair /= prob_tot;

      // first interval:  Compton interaction
      // second interval: pair production
      // third interval:  photo-electric absorption
      prob_pair += prob_comp;

      if (zeta[1] <= prob_comp) {
         // Compton scattering

         if (first_comp) {
            compton.interaction(p.energy, zeta[2], rndm,
                                x_comp_energy, cos_t_cx, sin_t_cx,
                                e_comp_energy, cos_t_ce, sin_t_ce);

            // create the Compton electron history in water
            multi_electron_create(electron_comp, e_comp_energy, rndm);

            first_comp = false;
         }

         // azimuthal scattering angle
         phi = TWO_PI * rndm.number();
         sin_phi = sin(phi);
         cos_phi = cos(phi);

         // electron direction
         e_comp.dir = p.dir;
         rotate(cos_t_ce, sin_t_ce, cos_phi, sin_phi, e_comp.dir);

         // trace the Compton electron
         e_comp.energy = e_comp_energy;
         e_comp.weight = p.weight*delta_eta;
         e_comp.pos    = p.pos;
         e_comp.i      = p.i;
         multi_electron_trace(electron_comp, e_comp, rndm);

         // with the photon play Russian Roulette
         if (rndm.number() < delta_eta) {
            // photon direction
            sin_phi = -sin_phi;
            cos_phi = -cos_phi;
            x_comp.dir = p.dir;
            rotate(cos_t_cx, sin_t_cx, cos_phi, sin_phi, x_comp.dir);

            x_comp.energy = x_comp_energy;
            x_comp.weight = p.weight;
            x_comp.pos    = p.pos;
            x_comp.i      = p.i;
            // simulate Compton photon
            if (x_comp_energy <= k1_cut) {
               // KERMA approximation for energies less than "k1_cut"
               portal_kerma(x_comp, rndm);
            }
            else {
               portal_photon(x_comp, rndm);
            }
         } // end of Russian Roulette

      } // end of Compton interaction
      else {
         // pair production or photo effect
         if (zeta[1] < prob_pair) {
            // pair production
            if (first_pair) {
               dif_pair.interaction(p.energy, zeta[2], rndm,
                                    e_pair_energy, cos_t_pe, sin_t_pe,
                                    p_pair_energy, cos_t_pp, sin_t_pp);

               // create the pair electron history in water
               multi_electron_create(electron_pair, e_pair_energy, rndm);

               // create the pair positron history in water
               // multi_positron_create(positron_pair, p_pair_energy, rndm);
               multi_electron_create(electron_pair, e_pair_energy, rndm);

               first_pair = false;
            }

            // azimuthal scattering angle
            phi = TWO_PI * rndm.number();
            sin_phi = sin(phi);
            cos_phi = cos(phi);

            // electron direction
            e_pair.dir = p.dir;
            rotate(cos_t_pe, sin_t_pe, cos_phi, sin_phi, e_pair.dir);

            // trace the pair electron
            e_pair.energy = e_pair_energy;
            e_pair.weight = p.weight*delta_eta;
            e_pair.pos    = p.pos;
            e_pair.i      = p.i;
            multi_electron_trace(electron_pair, e_pair, rndm);

            // positron direction
            sin_phi = -sin_phi;
            cos_phi = -cos_phi;
            p_pair.dir = p.dir;
            rotate(cos_t_pp, sin_t_pp, cos_phi, sin_phi, p_pair.dir);

            // trace the pair positron
            p_pair.energy = p_pair_energy;
            p_pair.weight = p.weight*delta_eta;
            p_pair.pos    = p.pos;
            p_pair.i      = p.i;
            // multi_positron_trace(positron_pair, p_pair, rndm);
            multi_electron_trace(positron_pair, p_pair, rndm);
         } // end of pair production
         else {
            // photo electric absorption
            if (first_phot) {
               // the photon energy is transfered to the electron
               e_phot_energy = p.energy;

               // create the photo electron history in water
               multi_electron_create(electron_phot, e_phot_energy, rndm);

               first_phot = false;
            }

            // trace the photo electron
            e_phot.energy = e_phot_energy;
            e_phot.weight = p.weight*delta_eta;
            // the photo electron has the direction of th incoming photon
            e_phot.dir    = p.dir;
            e_phot.pos    = p.pos;
            e_phot.i      = p.i;
            multi_electron_trace(electron_phot, e_phot, rndm);

         } // end of photo electric absorption

      }

   } // while(next_interaction)

   return;
}


void portal_dose::multi_electron_create(MULTI_ELECTRON &A,
                  real energy, ranmar &rndm) {
   int    i_e;          // energy index
   int    new_i_e;      // energy index for new_energy

   real   eta;          // random number in [0,1]
   real   zeta;         // random number in [0,1]
   real   num_mfp;      // number of mean free paths
   real   eloss;        // electron energy loss per step
   real   eloss_max;    // maximum electron energy loss per step
   real   new_energy;   // electron energy after one step
   real   sub_energy;   // electron energy after the first sub-step
   real   path;         // electron path for one whole step
   real   p2;           // momentum squared
   real   beta2;        // beta squared, (v/c)^2
   real   omega_0;      // number of elastic scatterings
   real   xr;
   real   cos_z,sin_z;  // multiple scattering angle
   real   sin_t_n;      // sin of discrete scattering angle (prim. electron)
   real   cos_t_n;      // cos of discrete scattering angle (prim. electron)
   real   energy_d;     // delta (Moller) electron energy
   real   sin_t_d;      // sin of delta electron scattering angle
   real   cos_t_d;      // cos of delta electron scattering angle
   real   rad_loss;     // bremsstrahlung energy (radiation loss)
   real   sin_t_p;      // sin of bremsstrahlung photon angle
   real   cos_t_p;      // cos of bremsstrahlung photon angle

   // program control
   bool   next_step,discrete_interaction;

// ****************************************
// begin electron transport
// ****************************************

   new_energy =  ZERO;
   next_step  =  true;
   A.n_step     =  0;
   A.n_delta    =  0;
   num_mfp    = -1.0;

   while (next_step) {
   // generate one step

      if (A.n_step >= A.MAX_STEP-1) {
         xvmc_warning("multi_electron_create",
                      "too many electron steps (n_step >= MAX_STEP-1)",8);
      }

      if ( (energy <= e_cut) || (A.n_step >= A.MAX_STEP-1) ) {
         next_step                  = false;
         A.energy_loss[A.n_step]    = energy;
         A.radiation_loss[A.n_step] = ZERO;
         A.sin_brems[A.n_step]      = ZERO;
         A.cos_brems[A.n_step]      = ZERO;
         A.dose_cor[A.n_step]       = e_h2o->dose_cor[0];
         eta                        = rndm.number();
         if (eta == ZERO)      eta  = rndm.number();
         A.random_number[A.n_step]  = eta;
         path                       = energy/e_h2o->s_res[0];
         A.step_size[A.n_step]      = path;
         A.s_cor[A.n_step]          = e_h2o->s_cor[0];
         A.alpha_r[A.n_step]        = e_h2o->alpha_r[0];
         discrete_interaction       = false;
         sub_energy                 = e_cut;
      }
      else {
      // begin: energy > e_cut

         // determine energy index
         i_e=int((energy-e_cut)*e_h2o->inverse_delta_e);

         // determine number of mean free paths to the next discrete interaction
         // using the fictitious interaction method (see SLAC-265, p. 16)
         if (num_mfp < ZERO) {
            num_mfp = ZERO;
            while (true) {
               zeta                   = rndm.number();
               if (zeta == ZERO) zeta = rndm.number();
               num_mfp   -= log(ONE-zeta);
               eloss      = num_mfp/e_h2o->sigma_max;
               new_energy = energy - eloss;
               if (new_energy > e_cut) {
                  new_i_e = int((new_energy-e_cut)*e_h2o->inverse_delta_e);
                  if (rndm.number() <= e_h2o->sigma_tot[new_i_e]) break;
               }
               else {
                  // the factor TWO is for numerical reasons only,
                  // it avoids any discrete interaction at the end of the step
                  num_mfp = TWO*energy*e_h2o->sigma_max;
                  break;
               }
            }
         }

         // determine energy loss depending on "num_mfp" and "e_h2o->max_loss"
         eloss     = num_mfp/e_h2o->sigma_max;
         eloss_max = e_h2o->max_loss[i_e]*energy;
         if (eloss > eloss_max) {
            // the maximum step size is reached
            eloss                = eloss_max;
            discrete_interaction = false;
         }
         else {
            // there is a discrete interaction at the end of the step
            discrete_interaction = true;
         }

         // determine and check "new_energy"
         new_energy = energy - eloss;
         if (new_energy <= e_cut) {
            eloss                = energy;
            new_energy           = ZERO;
            next_step            = false;
            discrete_interaction = false;
         }

         // the new number of mean free paths
         num_mfp -= eloss*e_h2o->sigma_max;

         // use a random number to divide the step into two sub-steps
         eta = rndm.number();
         if (eta == ZERO) eta = rndm.number();
         sub_energy = energy - eloss*eta; // energy after the first sub-step
         if (sub_energy < e_cut) sub_energy = e_cut;
         new_i_e = int((sub_energy-e_cut)*e_h2o->inverse_delta_e);

         // the electron path for the whole step
         path = eloss/e_h2o->s_res[new_i_e];

         // put the variables in the history scoring arrays
         A.energy_loss[A.n_step]    = eloss;
         A.radiation_loss[A.n_step] = ZERO;
         A.sin_brems[A.n_step]      = ZERO;
         A.cos_brems[A.n_step]      = ZERO;
         A.dose_cor[A.n_step]       = e_h2o->dose_cor[new_i_e];
         A.random_number[A.n_step]  = eta;
         A.step_size[A.n_step]      = path;
         A.s_cor[A.n_step]          = e_h2o->s_cor[new_i_e];
         A.alpha_r[A.n_step]        = e_h2o->alpha_r[new_i_e];

      // end: energy > e_cut
      }

      // at this point "eloss", "path" and "sub_energy" are fixed,
      // "new_energy" is the energy at the end of the step,
      // now perform a multiple scattering between the two sub-steps
      p2      = sub_energy*(sub_energy + TWOxEMASS); // momentum squared
      beta2   = p2/(p2 + EMASSxEMASS);               // beta squared, (v/c)^2
      omega_0 = BC*path/beta2;               // number of elastic scatterings
      xr      = CHI_CC2*path/p2/beta2/TWO;

      // sample multiple scattering angle (A.reduced_angle)
      mscat(omega_0, xr, cos_z, sin_z, rndm);

      // MS angle
      A.reduced_angle[A.n_step] = ONE - cos_z;

      // "energy" is now the energy at the end of the MS step
      energy = new_energy;

      // default: no Moller interaction
      A.moller[A.n_step] = false;

      // perform a discrete interaction
      if (discrete_interaction) {
         // at first, reset number of mean free paths for the next step
         num_mfp = -1.0;

         // determine the type of discrete interaction (Moller or brems. prod.)
         i_e = int((energy-e_cut)*e_h2o->inverse_delta_e);
         if (rndm.number() <= e_h2o->p_moller[i_e]){
            // no Moller interaction if the energy is too small
            // (in this case we jump to the next MS step)
            if ( h2o_moller.interaction(energy,     rndm,
                                        new_energy, cos_t_n, sin_t_n,
                                        energy_d,   cos_t_d, sin_t_d) ) {
               if (A.n_delta >= A.MAX_DELTA) {
                  xvmc_error("multi_electron_create",
                    "too many delta electrons (A.n_delta >= MAX_DELTA)",8);
               }

               // count Moller interaction (delta electron production)
               A.moller[A.n_step]        = true;
               A.i_delta[A.n_step]       = A.n_delta;
               A.energy_delta[A.n_delta] = energy_d;

               // store the polar scattering angles
               // (primary and delta electron)
               A.sin_theta[A.n_step]     = sin_t_n;
               A.cos_theta[A.n_step]     = cos_t_n;
               A.sin_theta_d[A.n_delta]  = sin_t_d;
               A.cos_theta_d[A.n_delta]  = cos_t_d;

               // "energy" is now the energy after the Moller interaction
               energy = new_energy;

               ++A.n_delta;
            }
         } // end of Moller interaction
         else { // bremsstrahlung production
           // if the energy is larger than TC (corresponds to AP in EGS4)
            if ( h2o_brems.interaction(energy,     rndm,
                                       new_energy, cos_t_n, sin_t_n,
                                       rad_loss,   cos_t_p, sin_t_p) ) {
               A.radiation_loss[A.n_step] = rad_loss;
               A.sin_brems[A.n_step]      = sin_t_p;
               A.cos_brems[A.n_step]      = cos_t_p;
               energy = new_energy;
            }
         } // end of bremsstrahlung production
      }
      // end of discrete interaction

      ++A.n_step;
      // end of this step
   }

   return;

}

// trace pre-calculated electron history through the calculation cube

void portal_dose::multi_electron_trace(MULTI_ELECTRON &A,
                          particle_parameters &e, ranmar &rndm) {
   real   eta;          // random number in [0,1]
   real   total_sh2o;   // water step size of the first electron sub-step
   real   total_loss;   // energy loss during the first electron sub-step
   real   sub_sh2o;     // water step size of the second electron sub-step
   real   sub_loss;     // energy loss during the second electron sub-step
   real   alpha;        // A.alpha_r for this step
   real   dose_factor;  // correction for dose, ionization etc.
   real   f_col;        // function f_c(rho) to calculate the collision
                        // stopping power (see eq. (19) in VMC paper I)
   real   f_tot;        // total stopping power factor of the present voxel
   real   factor;       // ratio f_tot/f_col
   real   s_correction; // low energy correction for collision stopping power
   real   path_step;    // electron path length from the beginning of the step
   real   voxel_step;   // step size to the next voxel boundary
   real   voxel_sh2o;   // water step size to the next voxel boundary
   real   voxel_loss;   // energy loss in the present voxel
   real   voxel_depo;   // deposited energy in the present voxel

   // Russian Roulette probability for bremsstrahlung photons
   const real P_ROULETTE = 0.025;

   particle_parameters p; // bremsstrahlung photon
   particle_parameters d; // delta electron

   real_3 step;         // step sizes to the next x-, y-, z-voxel boundaries
   int_3  istep;        // step direction to the next voxel boundary (-1,0,1)
   int_3  inew;         // new voxel indices

   // auxiliary variables
   real   cos_t, sin_t, cos_phi, sin_phi;
   real   temp_phi, temp_x, temp_y, temp_x1, temp_y1;

   // double
   const double DTWO   = 1.99999999;
   const double DSMALL = 1.0e-10;
   double save_eloss;   // saves "total_loss" for later use
   double sum_chi_cc2;  // MS
   double temp, sin_z;

   // program control
   bool repeat_step, go_back, multiple_scattering;

// ****************************************
// begin electron transport
// ****************************************

   // the electron must be within the calculation cube
   if ( (e.i.x < 0) || (e.i.x >= p_dim.x) ||
        (e.i.y < 0) || (e.i.y >= p_dim.y) ||
        (e.i.z < 0) || (e.i.z >= p_dim.z) ) {
      return;
   }

   sub_sh2o = ZERO;
   sub_loss = ZERO;

   for (register int i_step=0; i_step<A.n_step; ++i_step) {
   // simulate electron step

#ifdef ELECTRON_TRACK_FULL_STOP
      // deposit the remaining energy in the present voxel
      if (i_step == A.n_step-1) {
         total_loss = sub_loss + A.energy_loss[i_step];
         batch_portal_dose->matrix[e.i.x][e.i.y] += total_loss*e.weight;
         return;
      }
#endif // ELECTRON_TRACK_FULL_STOP

      // transport electron until zero energy by CSDA

      // get parameters from stack
      eta          = A.random_number[i_step];
      total_sh2o   = sub_sh2o + eta*A.step_size[i_step];
      total_loss   = sub_loss + eta*A.energy_loss[i_step];
      sub_sh2o     = (ONE-eta)*A.step_size[i_step];
      sub_loss     = (ONE-eta)*A.energy_loss[i_step];
      alpha        = A.alpha_r[i_step];
      dose_factor  = A.dose_cor[i_step];
      s_correction = A.s_cor[i_step];

      // sample multiple scattering angle after the first sub-step
      multiple_scattering = true;

      do {
      // we go back to this point if there is a discrete interaction
      // at the and of the step or the present step is the last step,
      // i.e. "i_step=A.n_step-1", the bool variable "go_back" must be
      // "true" in both cases

         // find the X-, Y-, Z-distance to the next voxel boundary
         // check X-direction
         if (e.dir.x > ZERO) {
            istep.x = 1;
            step.x  = (p_size.x*double(e.i.x+1)-e.pos.x)/e.dir.x;
         }
         else {
            if (e.dir.x < ZERO) {
               istep.x = -1;
               step.x  = (p_size.x*double(e.i.x)-e.pos.x)/e.dir.x;
            }
            else {
               istep.x = 0;
               step.x  = HUGE_STEP;
            }
         }

         // check Y-direction
         if (e.dir.y > ZERO) {
            istep.y = 1;
            step.y  = (p_size.y*double(e.i.y+1)-e.pos.y)/e.dir.y;
         }
         else {
            if (e.dir.y < ZERO) {
               istep.y = -1;
               step.y  = (p_size.y*double(e.i.y)-e.pos.y)/e.dir.y;
            }
            else
            {
               istep.y = 0;
               step.y  = HUGE_STEP;
            }
         }

         // check Z-direction
         if (e.dir.z > ZERO) {
            istep.z = 1;
            step.z  = (p_size.z*double(e.i.z+1)-e.pos.z)/e.dir.z;
         }
         else {
            if (e.dir.z < ZERO) {
               istep.z = -1;
               step.z  = (p_size.z*double(e.i.z)-e.pos.z)/e.dir.z;
            }
            else {
               istep.z = 0;
               step.z  = HUGE_STEP;
            }
         }

         // update voxel indices
         inew = e.i;

         // start transport
         repeat_step = true;
         path_step   = ZERO;
         save_eloss  = total_loss;
         sum_chi_cc2 = ZERO;
         while (repeat_step) {
         // the electron will be traced to the next voxel boundary, if the
         // electron reaches the boundary we repeat this step in the next
         // voxel, the bool variable "repeat_step" is "true" in this case

            // find the next voxel
            if ( (step.z < step.x) && (step.z < step.y) ) {
               voxel_step = step.z;
               step.x -= voxel_step;
               step.y -= voxel_step;
               step.z  = p_size.z/fabs(e.dir.z);
               inew.z = e.i.z + istep.z;
            }
            else {
               if (step.x < step.y) {
                  voxel_step = step.x;
                  step.x  = p_size.x/fabs(e.dir.x);
                  step.y -= voxel_step;
                  step.z -= voxel_step;
                  inew.x = e.i.x + istep.x;
               }
               else {
                  voxel_step = step.y;
                  step.x -= voxel_step;
                  step.y  = p_size.y/fabs(e.dir.y);
                  step.z -= voxel_step;
                  inew.y = e.i.y + istep.y;
               }
            }

            // collision stopping power factor of the present voxel
            f_col = DENS_CCOL - DENS_CCOL*s_correction;

            // total stopping power factor of the present voxel
            factor =
            (ONE + alpha*DENS_SRAD)/(ONE + alpha);
            f_tot = f_col*factor;

            // corresponding step size in water
            voxel_sh2o = voxel_step*f_tot;

            if (voxel_sh2o >= total_sh2o) {
               // the step ends in this voxel
               repeat_step = false;
               path_step  += total_sh2o/f_tot;
               voxel_depo  = total_loss/factor;

               // put the energy in the present voxel
               batch_portal_dose->matrix[e.i.x][e.i.y] +=
                  voxel_depo*e.weight*dose_factor;

               // MS
               sum_chi_cc2 +=
                  total_loss*DENS_FCHI;
            }
            else {
               // the step continues in the next voxel
               voxel_loss  = total_loss*voxel_sh2o/total_sh2o;
               voxel_depo  = voxel_loss/factor;
               total_loss -= voxel_loss;
               total_sh2o -= voxel_sh2o;
               path_step  += voxel_step;

               // put the energy in the present voxel
               batch_portal_dose->matrix[e.i.x][e.i.y] +=
                  voxel_depo*e.weight*dose_factor;

               // MS
               sum_chi_cc2 +=
                  voxel_loss*DENS_FCHI;
               if ( (inew.x < 0) || (inew.x >= p_dim.x) ||
                    (inew.y < 0) || (inew.y >= p_dim.y) ||
                    (inew.z < 0) || (inew.z >= p_dim.z) ) {
                  // the electron leaves the calculation cube, the history
                  // ends here, test energy conservation
                  return;
               }

               // update voxel indices
               e.i = inew;

            } // if (voxel_sh2o >= total_sh2o)

         // if "repeat_step=true" we continue the step in the next voxel
         }

         // move the particle
         e.pos.x += path_step*e.dir.x;
         e.pos.y += path_step*e.dir.y;
         e.pos.z += path_step*e.dir.z;

         go_back = false;
         if (multiple_scattering) {
            // sample multiple scattering angle

            if (save_eloss > DSMALL)
                 temp = A.reduced_angle[i_step]*sum_chi_cc2/save_eloss;
            else temp = ZERO;
            if (temp < DTWO) {
               // scattering angles
               cos_t    = ONE - temp;
               sin_t    = sqrt(temp*(ONE+cos_t));
               real phi = TWO_PI * rndm.number(); // azimuthal angle
               sin_phi  = sin(phi);
               cos_phi  = cos(phi);
               sin_z    = (ONE-e.dir.z)*(ONE+e.dir.z);

               // rotate electron
               if (sin_z > DSMALL) {
                  sin_z    = sqrt(sin_z);
                  temp     = sin_t/sin_z;
                  temp_phi = e.dir.z*cos_phi;
                  temp_x   = e.dir.x*cos_t;
                  temp_y   = e.dir.y*cos_t;
                  temp_x1  = temp_phi*e.dir.x - e.dir.y*sin_phi;
                  temp_y1  = temp_phi*e.dir.y + e.dir.x*sin_phi;
                  e.dir.x  = temp*temp_x1 + temp_x;
                  e.dir.y  = temp*temp_y1 + temp_y;
                  e.dir.z  = e.dir.z*cos_t - sin_z*sin_t*cos_phi;
               }
               else {
                  e.dir.x = sin_t*cos_phi;
                  e.dir.y = sin_t*sin_phi;
                  e.dir.z = cos_t;
               } // if (sin_z > DSMALL)
            }
            else {
               temp     = DTWO*rndm.number();
               e.dir.z  = ONE - temp;
               sin_z    = sqrt(temp*(TWO-temp));
               real phi = TWO_PI * rndm.number(); // azimuthal angle
               e.dir.x  = sin_z*cos(phi);
               e.dir.y  = sin_z*sin(phi);
            } // if (temp < DTWO)

            if (A.moller[i_step]) {
               go_back = true;
               multiple_scattering = false;
               total_sh2o = sub_sh2o;
               total_loss = sub_loss;
               sub_sh2o   = ZERO;
               sub_loss   = ZERO;
            }

            if (i_step == A.n_step-1) {
               go_back = true;
               multiple_scattering = false;
               total_sh2o = sub_sh2o;
               total_loss = sub_loss;
               sub_sh2o   = ZERO;
               sub_loss   = ZERO;
            }

         }  // end of multiple scattering

      // if "go_back=true" we perform the second sub-step separately,
      // otherwise the second part of this step is added to the first
      // part of the next step
      }
      while (go_back);

      // follow brems photons with KERMA approximation
      p.energy = A.radiation_loss[i_step];

      if (p.energy > ZERO) {
         // play Russian Roulette
         if (rndm.number() < P_ROULETTE) {
            // get polar photon angle
            sin_t = A.sin_brems[i_step];
            cos_t = A.cos_brems[i_step];

            // sample azimuthal photon angle
            real phi = TWO_PI * rndm.number();
            sin_phi  = sin(phi);
            cos_phi  = cos(phi);

            // rotate bremsstrahlung photon
            p.dir = e.dir;
            rotate(cos_t, sin_t, cos_phi, sin_phi, p.dir);

            // set photon parameters
            p.type   = PHOTON;
            p.weight = e.weight/P_ROULETTE;
            p.pos    = e.pos;
            p.i      = e.i;
            // simulate photon
            portal_kerma(p,rndm);
         }
      }

      if ( A.moller[i_step] ) {
         // add delta electron contribution

         // get array index
         int i_d = A.i_delta[i_step];

         // calculate moving directions of the primary and delta electrons
         real phi = TWO_PI * rndm.number();
         sin_phi  = sin(phi);
         cos_phi  = cos(phi);
         sin_z    = (ONE-e.dir.z)*(ONE+e.dir.z);
         if ( sin_z > DSMALL ) {
            sin_z = sqrt(sin_z);

            // rotate delta electron
            cos_t    = A.cos_theta_d[i_d];
            sin_t    = A.sin_theta_d[i_d];
            temp     = sin_t/sin_z;
            temp_phi = e.dir.z*cos_phi;
            d.dir.x  = e.dir.x*cos_t + temp*(temp_phi*e.dir.x-e.dir.y*sin_phi);
            d.dir.y  = e.dir.y*cos_t + temp*(temp_phi*e.dir.y+e.dir.x*sin_phi);
            d.dir.z  = e.dir.z*cos_t - sin_z*sin_t*cos_phi;

            // rotate primary electron
            cos_t    = A.cos_theta[i_step];
            sin_t    = A.sin_theta[i_step];
            temp     = sin_t/sin_z;
            temp_phi = -e.dir.z*cos_phi;
            temp_x   = e.dir.x*cos_t;
            temp_y   = e.dir.y*cos_t;
            temp_x1  = temp_phi*e.dir.x + e.dir.y*sin_phi;
            temp_y1  = temp_phi*e.dir.y - e.dir.x*sin_phi;
            e.dir.x  = temp*temp_x1 + temp_x;
            e.dir.y  = temp*temp_y1 + temp_y;
            e.dir.z  = sin_z*sin_t*cos_phi + e.dir.z*cos_t;
         }
         else {
            // rotate delta electron
            cos_t    = A.cos_theta_d[i_d];
            sin_t    = A.sin_theta_d[i_d];
            d.dir.x  = sin_t*cos_phi;
            d.dir.y  = sin_t*sin_phi;
            d.dir.z  = cos_t;

            // rotate primary electron
            cos_t   = A.cos_theta[i_step];
            sin_t   = A.sin_theta[i_step];
            e.dir.x = -sin_t*cos_phi;
            e.dir.y = -sin_t*sin_phi;
            e.dir.z = cos_t;

         } // if ( sin_z > DSMALL )

         // set delta electron parameters
         d.type   = ELECTRON;
         d.energy = A.energy_delta[i_d];
         d.weight = e.weight;
         d.pos    = e.pos;
         d.i      = e.i;

         // transport delta electron
         one_electron(d,rndm);

      } // end of delta electron contribution

   // end of this step
   }

   return;
}

// One Electron
void portal_dose::one_electron(particle_parameters &e, ranmar &rndm) {
   int    i_e;          // energy index
   int    new_i_e;      // energy index for new_energy

   real_3 step;         // step sizes to the next x-, y-, z-voxel boundaries
   int_3  istep;        // step direction to the next voxel boundary (-1,0,1)
   int_3  inew;         // new voxel indices

   real   phi;          // azimuthal angle
   real   eta;          // random number in [0,1]
   real   zeta;         // random number in [0,1]
   real   num_mfp;      // number of mean free paths
   real   eloss;        // electron energy loss per step
   real   eloss_max;    // maximum electron energy loss per step
   real   total_sh2o;   // water step size of the first electron sub-step
   real   total_loss;   // energy loss during the first electron sub-step
   real   sub_sh2o;     // water step size of the second electron sub-step
   real   sub_loss;     // energy loss during the second electron sub-step
   real   new_energy;   // electron energy after one step
   real   sub_energy;   // electron energy after the first sub-step
   real   path;         // electron path for one whole step
   real   alpha;        // alpha_r for this step
   real   dose_factor;  // correction for dose, ionization etc.
   real   f_col;        // function f_c(rho) to calculate the collision
                        // stopping power (see eq. (19) in VMC paper I)
   real   f_tot;        // total stopping power factor of the present voxel
   real   factor;       // ratio f_tot/f_col
   real   s_correction; // low energy correction for collision stopping power
   real   path_step;    // electron path length from the beginning of the step
   real   voxel_step;   // step size to the next voxel boundary
   real   voxel_sh2o;   // water step size to the next voxel boundary
   real   voxel_loss;   // energy loss in the present voxel
   real   voxel_depo;   // deposited energy in the present voxel
   real   p2;           // momentum squared
   real   beta2;        // beta squared, (v/c)^2
   real   omega_0;      // number of elastic scatterings
   real   xr;
   real   cos_t,sin_t;     // multiple scattering angle
   real   cos_phi,sin_phi; // multiple scattering angle
   real   sin_t_n;      // sin of discrete scattering angle (prim. electron)
   real   cos_t_n;      // cos of discrete scattering angle (prim. electron)
   real   sin_t_d;      // sin of delta electron scattering angle
   real   cos_t_d;      // cos of delta electron scattering angle
   real   sin_t_p;      // sin of bremsstrahlung photon angle
   real   cos_t_p;      // cos of bremsstrahlung photon angle

   // double
   const double DSMALL = 1.0e-10;
   double save_eloss;   // saves "total_loss" for later use
   double sum_chi_cc2;  // MS

   particle_parameters p; // bremsstrahlung photon
   particle_parameters d; // delta electron

   // Russian Roulette probability for bremsstrahlung photons
   const real P_ROULETTE = 0.025;

   // program control
   bool next_step, repeat_step, go_back;
   bool multiple_scattering, discrete_interaction;

// ****************************************
// begin electron transport
// ****************************************
	//e.energy = 2.0;
	//e.i.x = 50; e.i.y = 50; e.i.z = 0;
	//e.dir.x = 0; e.dir.y = 0; e.dir.z = 1;
	//e.pos.x = 20.0; e.pos.y = 20.0; e.pos.z = 0.0;
	// cout << "---- e.energy = " << e.energy << " at " << e.i.x << " " << e.i.y << " " << e.i.z << endl;
   // the electron must be within the calculation cube
   if ( (e.i.x < 0) || (e.i.x >= p_dim.x) ||
        (e.i.y < 0) || (e.i.y >= p_dim.y) ||
        (e.i.z < 0) || (e.i.z >= p_dim.z) ) {
      return;
   }

   next_step =  true;
   num_mfp   = -1.0;
   sub_sh2o  = ZERO;
   sub_loss  = ZERO;

   while (next_step) {
   // generate one step

      if (e.energy <= e_cut) {
#ifdef ELECTRON_TRACK_FULL_STOP
         // deposit the remaining energy in the present voxel
         batch_portal_dose->matrix[e.i.x][e.i.y] += e.energy*e.weight;
         return;
#else // transport electron until zero energy by CSDA
         next_step            = false;
         discrete_interaction = false;
         path                 = e.energy/e_h2o->s_res[0];
         s_correction         = e_h2o->s_cor[0];
         eloss                = e.energy;
         sub_energy           = e_cut;
         new_energy           = ZERO;
         alpha                = e_h2o->alpha_r[0];
         dose_factor          = e_h2o->dose_cor[0];
         eta                  = rndm.number();
         if (eta == ZERO) eta = rndm.number();
#endif // ELECTRON_TRACK_FULL_STOP
      }
      else {
      // begin: e.energy > e_cut
      //   cout << e.energy << endl;
         // determine energy index
         i_e=int((e.energy-e_cut)*e_h2o->inverse_delta_e);
			// cout << "e.energy = " << e.energy << " " << e_cut << " " <<e_h2o->inverse_delta_e << endl;
         // cout << "i_e = " << i_e << endl;;
         // determine number of mean free paths to the next discrete interaction
         // using the fictitious interaction method (see SLAC-265, p. 16)
         if (num_mfp < ZERO) {
            num_mfp = ZERO;
            while (true) {
               zeta                   = rndm.number();
               if (zeta == ZERO) zeta = rndm.number();
               num_mfp   -= log(ONE-zeta);
               eloss      = num_mfp/e_h2o->sigma_max;
               new_energy = e.energy - eloss;
               if (new_energy > e_cut) {
                  new_i_e = int((new_energy-e_cut)*e_h2o->inverse_delta_e);
                  if (rndm.number() <= e_h2o->sigma_tot[new_i_e]) break;
               }
               else {
                  // the factor TWO is for numerical reasons only,
                  // it avoids any discrete interaction at the end of the step
                  num_mfp = TWO*e.energy*e_h2o->sigma_max;
                  break;
               }
            }
         }

         // determine energy loss depending on "num_mfp" and "e_h2o->max_loss"
         eloss     = num_mfp/e_h2o->sigma_max;
         eloss_max = e_h2o->max_loss[i_e]*e.energy;
         if (eloss > eloss_max) {
            // the maximum step size is reached
            eloss                = eloss_max;
            discrete_interaction = false;
         }
         else {
            // there is a discrete interaction at the end of the step
            discrete_interaction = true;
         }

         // determine and check "new_energy"
         new_energy = e.energy - eloss;
			// cout << "eloss = " << eloss << "  " << eloss_max << endl;
			// cout << "new_energy = " << new_energy << endl;
         if (new_energy <= e_cut) {
            eloss                = e.energy;
            new_energy           = ZERO;
            next_step            = false;
            discrete_interaction = false;
         }

         // the new number of mean free paths
         num_mfp -= eloss*e_h2o->sigma_max;

         // use a random number to divide the step into two sub-steps
         eta = rndm.number();
         if (eta == ZERO) eta = rndm.number();
         sub_energy = e.energy - eloss*eta; // energy after the first sub-step
         if (sub_energy < e_cut) sub_energy = e_cut;
         new_i_e = int((sub_energy-e_cut)*e_h2o->inverse_delta_e);

         // the electron path length for the whole step in water
         path = eloss/e_h2o->s_res[new_i_e];

         // put the variables in the history scoring arrays
         dose_factor    = e_h2o->dose_cor[new_i_e];
         s_correction   = e_h2o->s_cor[new_i_e];
         alpha          = e_h2o->alpha_r[new_i_e];
			// cout <<"path = " << path << " " << dose_factor << " " << s_correction << " " << alpha << endl;
      // end: e.energy > e_cut
      }

      // at this point "eloss", "path" and "sub_energy" are fixed,
      // "sub_energy" is the energy at the end of the first sub-step,
      // "new_energy" is the energy at the end of the whole step

      // now apply this step to the patient
      total_sh2o = sub_sh2o + eta*path;
      total_loss = sub_loss + eta*eloss;
      sub_sh2o   = (ONE-eta)*path;
      sub_loss   = (ONE-eta)*eloss;

      // sample multiple scattering angle after the first sub-step
      multiple_scattering = true;
// cout << "voxel_size = " << p_size.x << " " << p_size.y << " " << p_size.z << endl;
// cout << "pos = " << e.pos.x << " " << e.pos.y << " " << e.pos.z << endl;
      do {
      // we go back to this point if there is a discrete interaction
      // at the and of the step or the present step is the last step,
      // i.e. "next_step=false", the bool variable "go_back" must be
      // "true" in both cases

         // find the X-, Y-, Z-distance to the next voxel boundary
         // check X-direction
         if (e.dir.x > ZERO) {
            istep.x = 1;
            step.x  = (p_size.x*double(e.i.x+1)-e.pos.x)/e.dir.x;
         }
         else {
            if (e.dir.x < ZERO) {
               istep.x = -1;
               step.x  = (p_size.x*double(e.i.x)-e.pos.x)/e.dir.x;
            }
            else {
               istep.x = 0;
               step.x  = HUGE_STEP;
            }
         }

         // check Y-direction
         if (e.dir.y > ZERO) {
            istep.y = 1;
            step.y  = (p_size.y*double(e.i.y+1)-e.pos.y)/e.dir.y;
         }
         else {
            if (e.dir.y < ZERO) {
               istep.y = -1;
               step.y  = (p_size.y*double(e.i.y)-e.pos.y)/e.dir.y;
            }
            else {
               istep.y = 0;
               step.y  = HUGE_STEP;
            }
         }

         // check Z-direction
         if (e.dir.z > ZERO) {
            istep.z = 1;
            step.z  = (p_size.z*double(e.i.z+1)-e.pos.z)/e.dir.z;
         }
         else {
            if (e.dir.z < ZERO) {
               istep.z = -1;
               step.z  = (p_size.z*double(e.i.z)-e.pos.z)/e.dir.z;
            }
            else {
               istep.z = 0;
               step.z  = HUGE_STEP;
            }
         }

         // update voxel indices
         inew = e.i;

         // start transport
         repeat_step = true;
         path_step   = ZERO;
         save_eloss  = total_loss;
         sum_chi_cc2 = ZERO;
         while (repeat_step) {
         // the electron will be traced to the next voxel boundary, if the
         // electron reaches the boundary we repeat this step in the next
         // voxel, the bool variable "repeat_step" is "true" in this case
         // cout << "step = " << step.x << " " << step.y << " " << step.z << endl;
            // find the next voxel
            if ( (step.z < step.x) && (step.z < step.y) ) {
               voxel_step = step.z;
               step.x -= voxel_step;
               step.y -= voxel_step;
               step.z  = p_size.z/fabs(e.dir.z);
               inew.z = e.i.z + istep.z;
            }
            else {
               if (step.x < step.y) {
                  voxel_step = step.x;
                  step.x  = p_size.x/fabs(e.dir.x);
                  step.y -= voxel_step;
                  step.z -= voxel_step;
                  inew.x = e.i.x + istep.x;
               }
               else {
                  voxel_step = step.y;
                  step.x -= voxel_step;
                  step.y  = p_size.y/fabs(e.dir.y);
                  step.z -= voxel_step;
                  inew.y = e.i.y + istep.y;
               }
            }

            // collision stopping power factor of the present voxel
            f_col = DENS_CCOL - DENS_CCOL*s_correction;

            // total stopping power factor of the present voxel
            factor = (ONE + alpha*DENS_SRAD)/(ONE + alpha);
            f_tot = f_col*factor;

            // corresponding step size in water
            voxel_sh2o = voxel_step*f_tot;
// cout << "f_col = " << f_col << " " << factor << " " << voxel_sh2o << " " << total_sh2o << " " << voxel_step << endl;
            if (voxel_sh2o >= total_sh2o) {
               // the step ends in this voxel
               repeat_step = false;
               path_step  += total_sh2o/f_tot;
               voxel_depo  = total_loss/factor;

               // put the energy in the present voxel
               batch_portal_dose->matrix[e.i.x][e.i.y] +=
                  voxel_depo*e.weight*dose_factor;

               // MS
               sum_chi_cc2 += total_loss*DENS_FCHI;
            }
            else {
               // the step continues in the next voxel
               voxel_loss  = total_loss*voxel_sh2o/total_sh2o;
               voxel_depo  = voxel_loss/factor;
               total_loss -= voxel_loss;
               total_sh2o -= voxel_sh2o;
               path_step  += voxel_step;

               // put the energy in the present voxel
               batch_portal_dose->matrix[e.i.x][e.i.y] +=
                  voxel_depo*e.weight*dose_factor;

               // MS
               sum_chi_cc2 +=
                  voxel_loss*DENS_FCHI;

               if ( (inew.x < 0) || (inew.x >= p_dim.x) ||
                    (inew.y < 0) || (inew.y >= p_dim.y) ||
                    (inew.z < 0) || (inew.z >= p_dim.z) ) {
                  // the electron leaves the calculation cube, the history
                  // ends here, test energy conservation
                  return;
               }

               // update voxel indices
               e.i = inew;

            } // if (voxel_sh2o >= total_sh2o)

         // if "repeat_step=true" we continue the step in the next voxel
         }

         // move the particle
         e.pos.x += path_step*e.dir.x;
         e.pos.y += path_step*e.dir.y;
         e.pos.z += path_step*e.dir.z;
//cout << "path_step = " << path_step << " " << e.pos.x << " " << e.pos.y << " " << e.pos.z << endl;

         go_back = false;
         if (multiple_scattering) {
            // sample multiple scattering angle between the two sub-steps

            // momentum squared
            p2      = sub_energy*(sub_energy + TWOxEMASS);
            // beta squared, (v/c)^2
            beta2   = p2/(p2 + EMASSxEMASS);
            // number of elastic scatterings
            omega_0 = BC*path/beta2;
            xr      = CHI_CC2*path/p2/beta2/TWO;
            if (save_eloss > DSMALL) xr *= sum_chi_cc2/save_eloss;

            // sample multiple scattering angle (reduced_angle)

            mscat(omega_0, xr, cos_t, sin_t, rndm);

            phi = TWO_PI * rndm.number();
            sin_phi = sin(phi);
            cos_phi = cos(phi);

            // rotate electron direction by the multiple scattering angle
            rotate(cos_t, sin_t, cos_phi, sin_phi, e.dir);

            if (discrete_interaction) {
               go_back = true;
               multiple_scattering = false;
               total_sh2o = sub_sh2o;
               total_loss = sub_loss;
               sub_sh2o   = ZERO;
               sub_loss   = ZERO;
            }

            if (!next_step) {
               go_back = true;
               multiple_scattering = false;
               total_sh2o = sub_sh2o;
               total_loss = sub_loss;
               sub_sh2o   = ZERO;
               sub_loss   = ZERO;
            }

         }  // end of multiple scattering

      // if "go_back=true" we perform the second sub-step separately,
      // otherwise the second part of this step is added to the first
      // part of the next step

// cout << "go_back = " << go_back << endl;
      }
      while (go_back);
// cout << "after while(go_back) " << endl;
      // "e.energy" is now the energy at the end of the MS step
      e.energy = new_energy;
// cout << "discrete_interaction = " << discrete_interaction << " " << e.energy << endl;
      // perform a discrete interaction
      if (discrete_interaction) {
         // at first, reset number of mean free paths for the next step
         num_mfp = -1.0;

         // determine the type of discrete interaction (Moller or brems. prod.)
         i_e = int((e.energy-e_cut)*e_h2o->inverse_delta_e);
         if (rndm.number() <= e_h2o->p_moller[i_e]) {
            // no Moller interaction if the energy is too small
            // (in this case we jump to the next MS step)
            if ( h2o_moller.interaction(e.energy,   rndm,
                                        new_energy, cos_t_n, sin_t_n,
                                        d.energy,   cos_t_d, sin_t_d) ) {
               // sample azimuthal scattering angle
               phi = TWO_PI * rndm.number();
               sin_phi = sin(phi);
               cos_phi = cos(phi);

               // rotate delta electron
               d.dir = e.dir;
               rotate(cos_t_d, sin_t_d, cos_phi, sin_phi, d.dir);

               // set delta electron parameters
               d.type   = ELECTRON;
               d.weight = e.weight;
               d.pos    = e.pos;
               d.i      = e.i;

               // transport delta electron
               one_electron(d,rndm);

               // direction change of the primary electron
               cos_phi = -cos_phi;
               sin_phi = -sin_phi;

               // rotate primary electron
               rotate(cos_t_n, sin_t_n, cos_phi, sin_phi, e.dir);

               // "e.energy" is now the energy after the Moller interaction
               e.energy = new_energy;
            }
         } // end of Moller interaction
         else { // bremsstrahlung production
           // if the energy is larger than TC (corresponds to AP in EGS4)
            if ( h2o_brems.interaction(e.energy,   rndm,
                                       new_energy, cos_t_n, sin_t_n,
                                       p.energy,   cos_t_p, sin_t_p) ) {
               // "e.energy" is now the energy after bremsstrahlung production
               e.energy = new_energy;

               // play Russian Roulette
               if (rndm.number() < P_ROULETTE) {
                  // sample azimuthal photon angle
                  phi = TWO_PI * rndm.number();
                  sin_phi = sin(phi);
                  cos_phi = cos(phi);

                  // rotate bremsstrahlung photon
                  p.dir = e.dir;
                  rotate(cos_t_p, sin_t_p, cos_phi, sin_phi, p.dir);

                  // set photon parameters
                  p.type   = PHOTON;
                  p.weight = e.weight/P_ROULETTE;
                  p.pos    = e.pos;
                  p.i      = e.i;
                  // simulate photon
                  portal_kerma(p,rndm);
               }
            }

         } // end of bremsstrahlung production

      } // end of discrete interaction

   } // end of this step
   return;
}

