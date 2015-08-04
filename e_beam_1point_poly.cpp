/*****************************************************************************
 * e_beam_1point_poly.cpp:                                                   *
 *    class member functions for:                                            *
 *       e_beam_1point_poly: point source beam, poly-energetic electrons     *
 *                                                                           *
 * Copyright (C) 2002    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 31.01.2002      *
 *    triple electron source model                        MF 08.07.2002      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <sstream>
using namespace std;

#include "xvmc_util.h"
#include "e_beam_1point_poly.h"

// ****************************************
// declare functions and global variables
// ****************************************

real etime(void);
real trace_line(const real_3 &, const real_3 &, real &);

// ***********************************************
// member functions of class e_beam_1point_poly
// ***********************************************

// get base data from file, initialize energy spectrum
void e_beam_1point_poly::get_base_data()
{
   // define bool variables
   bool app_distance_found   = false;
   bool sigma_theta_x_found  = false;
   bool norm_value_found     = false;
   bool gray_mu_dmax_found   = false;
   bool energy_min_found     = false;
   bool energy_max_found     = false;
   bool energy_X_found       = false;
   bool spectrum_A_found     = false;
   bool spectrum_alpha_found = false;
   bool spectrum_beta_found  = false;
   bool photo_a_found        = false;
   bool photo_b_found        = false;

   // read lines
   char line[81]  = "";    // lines to read from base data file
   bool read_line = true;
   while (read_line)
   {
      if (base_file.eof())
      {
         xvmc_error("e_beam_1point_poly::get_base_data",
                    "end of file reached, base data incomplete",8);
      }
      base_file.getline(line,sizeof(line));
      istringstream line_stream(line);
      char keyword[81] = "";
      line_stream >> keyword;
      if (!strcmp(keyword,"END-PARAMETERS"))
      {
         xvmc_error("e_beam_1point_poly::get_base_data",
                    "end of parameter entry, base data incomplete",8);
      }
      else
      {
         if (!strcmp(keyword,"SAD:"))
         {
            line_stream >> app_distance;
            app_distance_found = true;
         }
         if (!strcmp(keyword,"SIGMA-THETA-X:"))
         {
            // sigma_theta_x in radians
            line_stream >> sigma_theta_x;
            sigma_theta_x_found = true;
         }
         if (!strcmp(keyword,"NORM-VALUE:"))
         {
            line_stream >> norm_value;
            norm_value_found = true;
         }
         if (!strcmp(keyword,"GY/MU-DMAX:"))
         {
            line_stream >> gray_mu_dmax;
            gray_mu_dmax_found = true;
         }
         if (!strcmp(keyword,"ENERGY-MIN:"))
         {
            line_stream >> energy_min;
            energy_min_found = true;
         }
         if (!strcmp(keyword,"ENERGY-MAX:"))
         {
            // maximum and most probable energy (MeV)
            line_stream >> energy_max;
            energy_max_found = true;
         }
         if (!strcmp(keyword,"ENERGY-X:"))
         {
            line_stream >> energy_X;
            energy_X_found = true;
         }
         if (!strcmp(keyword,"SPECTRUM-A:"))
         {
            line_stream >> spectrum_A;
            spectrum_A_found = true;
         }
         if (!strcmp(keyword,"SPECTRUM-ALPHA:"))
         {
            line_stream >> spectrum_alpha;
            spectrum_alpha_found = true;
         }
         if (!strcmp(keyword,"SPECTRUM-BETA:"))
         {
            line_stream >> spectrum_beta;
            spectrum_beta_found = true;
         }
         if (!strcmp(keyword,"GAMMA-A:"))
         {
            line_stream >> photo_a;
            photo_a_found = true;
         }
         if (!strcmp(keyword,"GAMMA-B:"))
         {
            line_stream >> photo_b;
            photo_b_found = true;
         }
      } // if (!strcmp(keyword,"END-PARAMETERS"))

      if ( app_distance_found   &&
           sigma_theta_x_found  &&
           norm_value_found     &&
           gray_mu_dmax_found   &&
           energy_min_found     &&
           energy_max_found     &&
           energy_X_found       &&
           spectrum_A_found     &&
           spectrum_alpha_found &&
           spectrum_beta_found  &&
           photo_a_found        &&
           photo_b_found           ) read_line = false;

   } // while (read_line)

   // parameter input was successful, now reset applicator openings
   app_x1       = open_x1*app_distance/iso_distance;
   app_x2       = open_x2*app_distance/iso_distance;
   app_y1       = open_y1*app_distance/iso_distance;
   app_y2       = open_y2*app_distance/iso_distance;
   app_width_x  = app_x2 - app_x1;
   app_width_y  = app_y2 - app_y1;

   // print beam model parameters
   xvmc_message("app_distance:   ",app_distance,"cm",1);
   xvmc_message("sigma_theta_x:  ",sigma_theta_x,"rad",0);
   xvmc_message("norm_value:     ",norm_value,"",0);
   xvmc_message("gray_mu_dmax:   ",gray_mu_dmax,"Gy",0);
   xvmc_message("energy_min:     ",energy_min,"MeV",0);
   xvmc_message("energy_max:     ",energy_max,"MeV",0);
   xvmc_message("energy_X:       ",energy_X,"MeV",0);
   xvmc_message("spectrum_A:     ",spectrum_A,"",0);
   xvmc_message("spectrum_alpha: ",spectrum_alpha,"",0);
   xvmc_message("spectrum_beta:  ",spectrum_beta,"",0);
   xvmc_message("photo_a:        ",photo_a,"",0);
   xvmc_message("photo_b:        ",photo_b,"",0);

   // initialize energy spectrum
   // sampling parameters
   real e1 = energy_X - energy_min;
   real e2 = energy_X - energy_max;
   e3      = energy_max - energy_min;
   real e4 = energy_max + energy_min;
   b  = e1*e2,
   c1 = e1*e1,
   c2 = e2*e2;
   d1 = energy_min*energy_min;
   d2 = energy_max*energy_max;
   r1 = ZERO;
   r2 = 1.0e-35;;
   
   // weights
   w1 = (c1-c2)/TWO/(c1*c2);
   w2 = spectrum_A*e3*e4/TWO;
   w3 = spectrum_A*TWOxEMASS*e3;
   real w4 = ZERO;
   if (spectrum_alpha > ZERO)
   {
      if( spectrum_beta > ZERO)
      {
         real FIFTY = 50.0;
         real temp1 = spectrum_beta*energy_min;
         if (temp1 < FIFTY)
         {
            r1 = exp(-temp1);
            w4 = r1*spectrum_alpha/spectrum_beta;
            real temp2 = spectrum_beta*e3;
            if (temp2 < FIFTY) w4 *= ONE-exp(-temp2);
            real EIGHTY = 80.0;
            real temp3 = spectrum_beta*energy_max;
            if (temp3 < EIGHTY) r2 = exp(-temp3);
         }
      }
   }

   // normalize and print weights
   real wtot = w1 + w2 + w3 + w4;
   w1 /= wtot;
   w2 /= wtot;
   w3 /= wtot;
   w4 /= wtot;
   xvmc_message("Initialize energy spectrum, contributions:",1);
   xvmc_message("  high energy (main): ",100.0*w1,"%",0);
   xvmc_message("  linear:             ",100.0*w2,"%",0);
   xvmc_message("  constant:           ",100.0*w3,"%",0);
   xvmc_message("  low energy:         ",100.0*w4,"%",0);

   // mean energy
   real e_mean1 = ZERO;
   if (w1 > ZERO)
   {
      e_mean1 = ((e1-energy_min)/c1 - (e2-energy_max)/c2)/TWO/wtot;
   }
   real e_mean2 = ZERO;
   if (w2 > ZERO)
   {
      e_mean2 =
        spectrum_A*(pow(energy_max,THREE)-pow(energy_min,THREE))/THREE/wtot;
   }
   real e_mean3 = ZERO;
   if (w3 > ZERO)
   {
      e_mean3 = spectrum_A*EMASS*e3*e4/wtot;
   }
   real e_mean4 = ZERO;
   if (w4 > ZERO)
   {
      e_mean4 =
        exp(-spectrum_beta*energy_min)*(ONE+spectrum_beta*energy_min);
      e_mean4 -=
        exp(-spectrum_beta*energy_max)*(ONE+spectrum_beta*energy_max);
      e_mean4 *= spectrum_alpha/pow(spectrum_beta,TWO)/wtot;
   }
   energy_mean = e_mean1 + e_mean2 + e_mean3 + e_mean4;
   xvmc_message("  average energy:     ",energy_mean,"MeV",0);
   xvmc_message("  maximum energy:     ",energy_max,"MeV",0);

   // cumulative weigths
   w2 += w1;
   w3 += w2;
}

// get particle parameters from poly-energetic electron point source
bool e_beam_1point_poly::emit(particle_parameters &p,
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
   p.energy = sample_energy(rndm);
   p.weight = ONE;

   // position
   p.pos.x  = app_x1 + app_width_x*rx;
   p.pos.y  = app_y1 + app_width_y*ry;
   p.pos.z  = app_distance;

   // origin to point distance
   origin_point_dist =
      sqrt(p.pos.x*p.pos.x + p.pos.y*p.pos.y + p.pos.z*p.pos.z);

   // preferential moving direction
   p.dir.x  = p.pos.x/origin_point_dist;
   p.dir.y  = p.pos.y/origin_point_dist;
   p.dir.z  = p.pos.z/origin_point_dist;

   // sample fluctuations (d_theta, d_phi) of the moving direction
   const real EPSI = 1.0e-20;
   real d_theta =  sigma_theta_x*SQRT2
                  *sqrt(max_of(ZERO,real(-log(max_of(EPSI,rndm.number())))));
   real d_phi = TWO_PI*rndm.number();

   // rotate electron direction by d_theta and d_phi
   rotate(d_theta, d_phi, p.dir);

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
bool e_beam_1point_poly::emit(real   &weight,
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

   // preferential moving direction
   dir.x  = pos.x/origin_point_dist;
   dir.y  = pos.y/origin_point_dist;
   dir.z  = pos.z/origin_point_dist;

   // sample fluctuations (d_theta, d_phi) of the moving direction
   const real EPSI = 1.0e-20;
   real d_theta = sigma_theta_x*SQRT2
                  *sqrt(max_of(ZERO,real(-log(max_of(EPSI,rndm.number())))));
   real d_phi = TWO_PI*rndm.number();

   // rotate electron direction by d_theta and d_phi
   rotate(d_theta, d_phi, dir);

   if ( modifier->transport(pos,dir,rndm) )
   {
      // trace particle to the simulation grid (take beam angles into account)
      // return true if the particle hits the calculation cube
      return( trace2cube(pos, dir, i, origin_point_dist, rndm) );
   }

   // there isn't any particle to emit, return false
   return(false);
}

// add photon background to the dose distribution
bool e_beam_1point_poly::add_photons(real &cpu_time, int n_batch)
{
   // measure CPU time
   real cpu_start = etime();

   // determine dose maximum
   float dose,dose_max = ZERO;
   for (register int k=0; k<dim.z; ++k)
   {
      for (register int j=0; j<dim.y; ++j)
      {
         for (register int i=0; i<dim.x; ++i)
         {
            dose = beam_dose->matrix[i][j][k];
            if (dose > dose_max) dose_max = dose;
         }
      }
   }
   xvmc_message("Maximum dose before photon correction:",
                 dose_max,"10^-10 Gy cm^2",1);

   // calculate scaled photon background depth parameters
   real p_a = photo_a*dose_max/100.0;
   real p_b = photo_b*dose_max/100.0;

   // calculate minimum and maximum beam angles for photon background
   real tan_x_min=(app_x1 - 10.0)/app_distance;
   real tan_x_max=(app_x2 + 10.0)/app_distance;
   real tan_y_min=(app_y1 - 10.0)/app_distance;
   real tan_y_max=(app_y2 + 10.0)/app_distance;

   // calculate photon background profile parameters
   real rx2 = 1.0 - (app_width_x-4.0)/42.0;
   rx2      = 0.66*(rx2*app_width_x)*(rx2*app_width_x);
   real ry2 = 1.0 - (app_width_y-4.0)/42.0;
   ry2      = 0.66*(ry2*app_width_y)*(ry2*app_width_y);

   // calculate photon background for all voxels
   // (fan lines from origin point to each voxel center)
   real_3 pos; // position of voxel center
   real_3 dop; // difference vector from origin to voxel center
   real   temp;

   pos.z = -voxel_size.z*ONE_HALF;
   for (register int k=0; k<dim.z; ++k)
   {
      pos.z +=  voxel_size.z;
      dop.z  =  pos.z - origin.z;

      pos.y  = -voxel_size.y*ONE_HALF;
      for (register int j=0; j<dim.y; ++j)
      {
         pos.y += voxel_size.y;
         dop.y  = pos.y - origin.y;

         pos.x  = -voxel_size.x*ONE_HALF;
         for (register int i=0; i<dim.x; ++i)
         {
            pos.x += voxel_size.x;
            dop.x  = pos.x - origin.x;

            // origin to voxel center distance
            real dist = sqrt(dop.x*dop.x + dop.y*dop.y + dop.z*dop.z);

            // fan line direction vector
            real_3 dir;
            dir.x = dop.x/dist;
            dir.y = dop.y/dist;
            dir.z = dop.z/dist;

            // rotate back by the table angle
            temp  = dir.x*sin_beta;
            dir.x = dir.x*cos_beta - dir.y*sin_beta;
            dir.y = temp           + dir.y*cos_beta;

            // rotate back by the gantry angle
            temp  = dir.x*sin_alpha;
            dir.x = dir.x*cos_alpha + dir.z*sin_alpha;
            dir.z = -temp           + dir.z*cos_alpha;

            // rotate back by the collimator angle
            temp  = dir.x*sin_gamma;
            dir.x = dir.x*cos_gamma - dir.y*sin_gamma;
            dir.y = temp            + dir.y*cos_gamma;

            // investigate fan lines
            if (fabs(dir.z) > ZERO)
            {
               real tan_x = dir.x/fabs(dir.z);
               real tan_y = dir.y/fabs(dir.z);

               // the fan line must be within the beam cone
               if ((tan_x > tan_x_min) && (tan_x < tan_x_max))
               {
                  if ((tan_y > tan_y_min) && (tan_y < tan_y_max))
                  {
                     // geometrical length within the cube (dummy variable)
                     real length = ZERO;

                     // effective length in water
                     real eff_length = trace_line(origin, pos, length);

                     // profile parameters
                     temp = app_distance*tan_x;
                     real bpx = exp(-temp*temp/rx2);
                     temp = app_distance*tan_y;
                     real bpy = exp(-temp*temp/ry2);

                     // depth parameters
                     real dose = p_b + p_a*eff_length;

                     // total photon dose (must be larger than 0)
                     dose *= bpx*bpy;
                     if (dose < ZERO) dose = ZERO;

                     // dose in vacuum or air is 0
                     real rho = density->matrix[i][j][k];
                     const real rho_min = 0.01;
                     dose = rho > rho_min ? dose : ZERO;

                     // the photon background should not contribute to the
                     // statistical uncertainty, the following update to the
                     // "beam_error->matrix" provides this functionality
                     real error =
                      (TWO*dose*beam_dose->matrix[i][j][k]+dose*dose)/n_batch;

                     // add photon background to beam dose and error matrices
                     beam_dose->matrix[i][j][k]  += dose;
                     beam_error->matrix[i][j][k] += error;
                  }
               }
            }
         }
      }
   }

   // measure CPU time
   cpu_time = etime() - cpu_start;

   return(true);
}
