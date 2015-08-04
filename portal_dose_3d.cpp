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
// End of Added -------------------------
#include "portal_dose.h"

// ****************************************
// declare global variables
// ****************************************
// Added by JOKim 25NOV2009 -------------

extern   photon_data   *p_h2o;    // water photon transport data for this beam

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
   p_dim.x  = new_dim_x;
   p_dim.y  = new_dim_y;
   p_dim.z  = 1;
   p_size.x = new_res_x;
   p_size.y = new_res_y;
   p_size.z = new_res_x;
   bmp_size = new_bmp_size;

   batch_portal_dose = new array_3d<float>(p_dim, ZERO);
   beam_portal_dose  = new array_3d<float>(p_dim, ZERO);
   beam_portal_error = new array_3d<float>(p_dim, ZERO);

	// Portal Imager Center Position in Beam Coordinates
   p_pos.x = 0.0;
   p_pos.y = 0.0;
   p_pos.z = distance + 100.0;

	// Portal Imager Size
   real   width = p_dim.x * p_size.x;
   real   length= p_dim.y * p_size.y;

	// Define 4 Points of Portal Imager in Beam Coordinates
   x1.x = p_pos.x - width * 0.5;
   x1.y = p_pos.y;
   x1.z = p_pos.z;

   x2.x = p_pos.x + width * 0.5;
   x2.y = p_pos.y;
   x2.z = p_pos.z;

   y1.x = p_pos.x;
   y1.y = p_pos.y - length* 0.5;
   y1.z = p_pos.z;

   y2.x = p_pos.x;
   y2.y = p_pos.y + length* 0.5;
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

   // Beam Direction Vector in Patient Model Coordinates
   // beam_dir.x   = -sin_alpha*cos_beta;
   // beam_dir.y   =  sin_alpha*sin_beta;
   // beam_dir.z   =  cos_alpha;

   // normalize position vector
   real origin_point_dist = distance + 100.0;
   p_dir.x  = p_pos.x/origin_point_dist;
   p_dir.y  = p_pos.y/origin_point_dist;
   p_dir.z  = p_pos.z/origin_point_dist;

   real temp = 0.0;
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

	// Portal Imager Limits
   min_limit.x  = x1.x;
   max_limit.x  = x2.x;
   min_limit.y  = y1.y;
   max_limit.y  = y2.y;
   min_limit.z  = x1.z;
   max_limit.z  = x2.z;

   // New Resolution of Pixel After Rotation
   p_size_rot.x = (max_limit.x - min_limit.x)/p_dim.x;
   p_size_rot.y = (max_limit.y - min_limit.y)/p_dim.y;
   p_size_rot.z = (max_limit.z - min_limit.z)/p_dim.x;

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
         real dose = batch_portal_dose->matrix[i][j][0]/batch_weight;
         beam_portal_dose->matrix[i][j][0]  += dose;
         beam_portal_error->matrix[i][j][0] += dose*dose;
      }
   }

// End of Added ------------------------------------------------------------------------
}

// Added by JOKim 25NOV2009  -----------------------------------------------------------
// add portal dose contribution of one photon
void portal_dose::add(particle_parameters p)
{

   int    i_e;            // energy index
   real   mu_tot_h2o;     // total cross section (attenuation coeff.) for water
   real   mu_en_h2o;      // energy absorption coefficient for water
   real   p_phot_h2o;     // probability for photo electric absorption in water
   real   p_pair_h2o;     // pair production probability in water
   real   p_comp_h2o;     // Compton interaction probability in water

   // determine total cross sections and interaction probabilities for water
   i_e=int((p.energy-p_cut)*p_h2o->inverse_delta_e);
   mu_tot_h2o = p_h2o->mu_tot[i_e];
   mu_en_h2o  = p_h2o->mu_en[i_e];
   p_phot_h2o = p_h2o->mu_phot[i_e]/mu_tot_h2o;
   p_pair_h2o = p_h2o->mu_pair[i_e]/mu_tot_h2o;
   if (p.energy <= TWOxEMASS) p_pair_h2o = ZERO;
   p_comp_h2o = ONE - p_pair_h2o - p_phot_h2o;



   // scalar products
   real product_pos = p.pos.x*normal.x + p.pos.y*normal.y + p.pos.z*normal.z;
   real product_dir = p.dir.x*normal.x + p.dir.y*normal.y + p.dir.z*normal.z;

   if (product_dir != ZERO) {
       real R = (dist2zero-product_pos)/product_dir;
       // if (R > 0) {
       real x = p.pos.x + R*p.dir.x;
       real y = p.pos.y + R*p.dir.y;
       real z = p.pos.z + R*p.dir.z;
       if (((min_limit.x - x) * (max_limit.x - x) < 0
	      || (min_limit.z - z) * (max_limit.z - z) < 0)
	      && (min_limit.y - y) * (max_limit.y - y) < 0) {
            int i = (int)((x - min_limit.x)/p_size_rot.x);
            int j = (int)((y - min_limit.y)/p_size_rot.y);
            int k = (int)((z - min_limit.z)/p_size_rot.z);
	         if (fabs(p_size_rot.z) > fabs(p_size_rot.z)){
	            i = (int)((z - min_limit.z)/p_size_rot.z);
	         }
            if (i < p_dim.x && j < p_dim.y && k < p_dim.z) {
               real wt = 1.0;
	            real pn = sqrt(p.dir.x*p.dir.x+p.dir.y*p.dir.y+p.dir.z*p.dir.z)
	                    * sqrt(normal.x*normal.x+normal.y*normal.y+normal.z*normal.z);
               if (pn > 1e-5) wt = product_dir/pn;
	            if (p.type == PHOTON) {
		            // determine total cross section
		            int  i_e = int((p.energy - p_cut)*p_h2o->inverse_delta_e);
						real mu_en_h2o = p_h2o->mu_en[i_e];
						// if (bmp_size == 0) mu_en_h2o = 1;
	         		batch_portal_dose->matrix[i][j][0]  += (mu_en_h2o * p.weight);
	      		}
	      		else {
						real e_energy = (p.energy - 0.511034);
						batch_portal_dose->matrix[i][j][0]  += (e_energy * p.weight);
	      		}
           }
	     }
	     // }
   }
   // portal_exit();
   // return(false);
}
// End of Addes  --------------------------------------------------------------------------

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
   for (register int j=0; j<p_dim.y; ++j) {
      for (register int i=0; i<p_dim.x; ++i) {
         dose  = beam_portal_dose->matrix[i][j][0];
         // convert average dose squared to standard deviation (error)
         error = float(n_batch)*beam_portal_error->matrix[i][j][0];
         error = sqrt(fabs(error-dose*dose)/float(n_batch-1));
         dose *= dose_factor;
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
