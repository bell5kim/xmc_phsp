/*****************************************************************************
 * e_beam_bumpy_mono.cpp:                                                    *
 *    class member functions for:                                            *
 *       e_beam_bumpy_mono: triple electron source model, mono-energetic     *
 *                          with bump correction                             *
 *                                                                           *
 * Copyright (C) 2005    Matthias Fippel                                     *
 *                       University of Tuebingen, Germany                    *
 *                                                                           *
 * revisions:                                                                *
 *    coded by modifying e_beam_triple_mono.cpp           MF 27.02.2005      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <sstream>
using namespace std;

#include "xvmc_util.h"
#include "e_beam_bumpy_mono.h"

// ***********************************************
// member functions of class e_beam_bumpy_mono
// ***********************************************

// base data file input
void e_beam_bumpy_mono::get_base_data()
{
   // define bool variables
   bool pri_weight1_found    = false;
   bool app_distance_found   = false;
   bool pri_xo2_found        = false;
   bool sigma_theta_x1_found = false;
   bool sigma_theta_x2_found = false;
   bool sct_weight_found     = false;
   bool sct_app_dist_found   = false;
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
   bool n_bump_found         = false;

   // read lines
   char line[81]  = "";    // lines to read from base data file
   bool read_line = true;
   while (read_line)
   {
      if (base_file.eof())
      {
         xvmc_error("e_beam_triple_mono::get_base_data",
                    "end of file reached, base data incomplete",8);
      }
      base_file.getline(line,sizeof(line));
      istringstream line_stream(line);
      char keyword[81] = "";
      line_stream >> keyword;
      if (!strcmp(keyword,"END-PARAMETERS"))
      {
         xvmc_error("e_beam_triple_mono::get_base_data",
                    "end of parameter entry, base data incomplete",8);
      }
      else
      {
         if (!strcmp(keyword,"PRIMARY-WEIGHT1:"))
         {
            // weight of primary source 1
            line_stream >> pri_weight1;
            pri_weight1_found = true;
         }
         if (!strcmp(keyword,"PRIMARY-SAD:"))
         {
            // source to applicator distance
            line_stream >> app_distance;
            app_distance_found = true;
         }
         if (!strcmp(keyword,"PRIMARY-XO2:"))
         {
            // Gaussian width squared of position fluctuation
            line_stream >> pri_xo2;
            pri_xo2_found = true;
         }
         if (!strcmp(keyword,"PRIMARY-SIGMA1:"))
         {
            // angular spread of primary source 1 (radians)
            line_stream >> sigma_theta_x1;
            sigma_theta_x1_found = true;
         }
         if (!strcmp(keyword,"PRIMARY-SIGMA2:"))
         {
            // angular spread of primary source 2 (radians)
            line_stream >> sigma_theta_x2;
            sigma_theta_x2_found = true;
         }
         if (!strcmp(keyword,"SCATTER-WEIGHT:"))
         {
            // weight of scatter source
            line_stream >> sct_weight;
            sct_weight_found = true;
         }
         if (!strcmp(keyword,"SCATTER-SAD:"))
         {
            // scatter source to applicator distance
            line_stream >> sct_app_dist;
            sct_app_dist_found = true;
         }
         if (!strcmp(keyword,"NORM-VALUE:"))
         {
            // normalization value
            line_stream >> norm_value;
            norm_value_found = true;
         }
         if (!strcmp(keyword,"GY/MU-DMAX:"))
         {
            // Gray per monitor unit at Dmax
            line_stream >> gray_mu_dmax;
            gray_mu_dmax_found = true;
         }
         if (!strcmp(keyword,"ENERGY-MIN:"))
         {
            // minimum energy for spectrum (MeV)
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
            // energy spectrum parameter X (MeV)
            line_stream >> energy_X_aux;
            energy_X_found = true;
         }
         if (!strcmp(keyword,"SPECTRUM-A:"))
         {
            // energy spectrum parameter A
            line_stream >> spectrum_A_aux;
            spectrum_A_found = true;
         }
         if (!strcmp(keyword,"SPECTRUM-ALPHA:"))
         {
            // energy spectrum parameter alpha
            line_stream >> spectrum_alpha_aux;
            spectrum_alpha_found = true;
         }
         if (!strcmp(keyword,"SPECTRUM-BETA:"))
         {
            // energy spectrum parameter beta
            line_stream >> spectrum_beta_aux;
            spectrum_beta_found = true;
         }
         if (!strcmp(keyword,"GAMMA-A:"))
         {
            // parameter a for photon background
            line_stream >> photo_a_aux;
            photo_a_found = true;
         }
         if (!strcmp(keyword,"GAMMA-B:"))
         {
            // parameter b for photon background
            line_stream >> photo_b_aux;
            photo_b_found = true;
         }
         if (!strcmp(keyword,"BUMP-CORRECTION:"))
         {
            // bump correction array size
            line_stream >> n_bump;
            n_bump_found = true;

            // create bump correction arrays
            if (n_bump > 0)
            {
               xy_bump = new real[n_bump];
               wx_bump = new real[n_bump];
               wy_bump = new real[n_bump];

               // input coordinates and bump correction weights
               for (register int i=0; i<n_bump; ++i)
               {
                  char bump_line[200] = "";
                  base_file.getline(bump_line,sizeof(bump_line));
                  istringstream bumps(bump_line);
                  bumps >> xy_bump[i] >> wx_bump[i] >> wy_bump[i];
               }
            }
         }
      } // if (!strcmp(keyword,"END-PARAMETERS"))

      if ( pri_weight1_found    &&
           app_distance_found   &&
           pri_xo2_found        &&
           sigma_theta_x1_found &&
           sigma_theta_x2_found &&
           sct_weight_found     &&
           sct_app_dist_found   &&
           norm_value_found     &&
           gray_mu_dmax_found   &&
           energy_min_found     &&
           energy_max_found     &&
           energy_X_found       &&
           spectrum_A_found     &&
           spectrum_alpha_found &&
           spectrum_beta_found  &&
           photo_a_found        &&
           photo_b_found        &&
           n_bump_found            ) read_line = false;

   } // while (read_line)

   // parameter input was successful, now reset applicator openings
   app_x1       = open_x1*app_distance/iso_distance;
   app_x2       = open_x2*app_distance/iso_distance;
   app_y1       = open_y1*app_distance/iso_distance;
   app_y2       = open_y2*app_distance/iso_distance;
   app_width_x  = app_x2 - app_x1;
   app_width_y  = app_y2 - app_y1;

   // print beam model parameters
   xvmc_message("pri_weight1:    ",pri_weight1*100.0,"%",1);
   xvmc_message("app_distance:   ",app_distance,"cm",0);
   xvmc_message("pri_xo2:        ",pri_xo2,"cm^2",0);
   xvmc_message("sigma_theta_x1: ",sigma_theta_x1,"rad",0);
   xvmc_message("sigma_theta_x2: ",sigma_theta_x2,"rad",0);
   xvmc_message("sct_weight:     ",sct_weight*100.0,"%",0);
   xvmc_message("sct_app_dist:   ",sct_app_dist,"cm",0);
   xvmc_message("norm_value:     ",norm_value,"",0);
   xvmc_message("gray_mu_dmax:   ",gray_mu_dmax,"Gy",0);
   xvmc_message("energy_min:     ",energy_min,"MeV",0);
   xvmc_message("energy_max:     ",energy_max,"MeV",0);
   xvmc_message("energy_X:       ",energy_X_aux,"MeV",0);
   xvmc_message("spectrum_A:     ",spectrum_A_aux,"",0);
   xvmc_message("spectrum_alpha: ",spectrum_alpha_aux,"",0);
   xvmc_message("spectrum_beta:  ",spectrum_beta_aux,"",0);
   xvmc_message("photo_a:        ",photo_a_aux,"",0);
   xvmc_message("photo_b:        ",photo_b_aux,"",0);
   xvmc_message("n_bump:         ",n_bump,"",0);
   for (int i=0; i<n_bump; ++i)
   {
      xvmc_message("x, y: ",xy_bump[i],"cm, weight_x:",wx_bump[i],
                                           "weight_y:",wy_bump[i],"",0);
   }

   if (energy_min > energy_max)
   {
      xvmc_error("e_beam_triple_mono::get_base_data()",
                 "minimum energy > maximum energy",8);
   }
}

// get particle parameters from the beam
bool e_beam_bumpy_mono::emit(particle_parameters &p,
#ifdef USE_SOBOL
                             sobseq &sobol, int &sobol_count,
#endif
                             ranmar &rndm)
{
   const real EPSI = 1.0e-20;    // a small parameter
   real origin_point_dist;       // origin to point distance

   // set particle type and initial weight
   p.type   = ELECTRON;
   p.weight = ONE;

   // random (Sobol) number to determine the source particle type
#ifdef USE_SOBOL
   real rp = sobol.number(sobol_count); ++sobol_count;
#else
   real rp = rndm.number();
#endif

   if (rp > sct_weight)
   {
      // this is an electron from the primary sources
      // choose electron energy
      p.energy = sample_energy(rndm);

      // random (Sobol) numbers to sample initial position
#ifdef USE_SOBOL
      real rx = sobol.number(sobol_count); ++sobol_count;
      real ry = sobol.number(sobol_count); ++sobol_count;
#else
      real rx = rndm.number();
      real ry = rndm.number();
#endif 

      // electron position in the applicator plane
      p.pos.x  = app_x1 + app_width_x*rx;
      p.pos.y  = app_y1 + app_width_y*ry;
      p.pos.z  = app_distance;

      // sample fluctuation of the electron starting position
      // in the applicator plane from Gaussian distribution
#ifdef USE_SOBOL
      real rr   = sobol.number(sobol_count); ++sobol_count;
      real rphi = sobol.number(sobol_count); ++sobol_count;
#else
      real rr   = rndm.number();
      real rphi = rndm.number();
#endif
      real rad_src =
        sqrt(TWO*pri_xo2*max_of(ZERO,real(-log(max_of(EPSI,rr)))));
      real phi_src = TWO_PI*rphi;
      p.pos.x += rad_src*cos(phi_src);
      p.pos.y += rad_src*sin(phi_src);

      // apply bump correction
		real radius = sqrt(p.pos.x*p.pos.x+p.pos.y*p.pos.y);
		real x_radius = radius;
		real y_radius = radius;
		if (p.pos.x < 0) x_radius = -x_radius;
		if (p.pos.y < 0) y_radius = -y_radius;	
		
		if (radius == 0.0) {
         p.weight *= interpolate(n_bump,xy_bump,wx_bump,p.pos.x);
         p.weight *= interpolate(n_bump,xy_bump,wy_bump,p.pos.y);
		} else {
			real thetax = acos(fabs(p.pos.x)/radius);
			real thetay = acos(fabs(p.pos.y)/radius);
			real weightx = interpolate(n_bump,xy_bump,wx_bump,x_radius);
			real weighty = interpolate(n_bump,xy_bump,wy_bump,y_radius);
			p.weight *= ((thetay*weightx+thetax*weighty)/(thetax+thetay));
		}

      // origin to point distance
      origin_point_dist =
         sqrt(p.pos.x*p.pos.x + p.pos.y*p.pos.y + p.pos.z*p.pos.z);

      // calculate preferential moving direction and normalize
      p.dir.x = p.pos.x/origin_point_dist;
      p.dir.y = p.pos.y/origin_point_dist;
      p.dir.z = p.pos.z/origin_point_dist;

      // sample fluctuations (d_theta, d_phi) of the moving direction
      real d_theta = ZERO;
      if (rndm.number() <= pri_weight1)
      {
         d_theta =  sigma_theta_x1*SQRT2
                   *sqrt(max_of(ZERO,real(-log(max_of(EPSI,rndm.number())))));
      }
      else
      {
         d_theta =  sigma_theta_x2*SQRT2
                   *sqrt(max_of(ZERO,real(-log(max_of(EPSI,rndm.number())))));
      }
      real d_phi = TWO_PI*rndm.number();

      // rotate electron direction by d_theta and d_phi
      rotate(d_theta, d_phi, p.dir);
   }
   else
   {
      // this is an electron from the scatter source
      // choose electron energy
      p.energy = energy_min + (energy_max-energy_min)*rndm.number();

      // random (Sobol) numbers to sample initial position
#ifdef USE_SOBOL
      real rx = sobol.number(sobol_count); ++sobol_count;
      real ry = sobol.number(sobol_count); ++sobol_count;
#else
      real rx = rndm.number();
      real ry = rndm.number();
#endif 

      // electron position in the applicator plane
      p.pos.x  = app_x1 + app_width_x*rx;
      p.pos.y  = app_y1 + app_width_y*ry;
      p.pos.z  = app_distance;

      // origin to point distance
      origin_point_dist =
         sqrt(p.pos.x*p.pos.x + p.pos.y*p.pos.y + p.pos.z*p.pos.z);

      // sample x position in the applicator scatter plane
      real x_src = ZERO;
      if (rndm.number() < rx) {
         x_src = app_x1 + (p.pos.x-app_x1)*sqrt(rndm.number()); }
      else {
         x_src = app_x2 + (p.pos.x-app_x2)*sqrt(rndm.number()); }

      // sample y position in the applicator scatter plane
      real y_src = ZERO;
      if (rndm.number() < ry) {
         y_src = app_y1 + (p.pos.y-app_y1)*sqrt(rndm.number()); }
      else {
         y_src = app_y2 + (p.pos.y-app_y2)*sqrt(rndm.number()); }

      // z coordinate is given by the position of scatter source plane
      real z_src = app_distance - sct_app_dist;

      // calculate moving direction and normalize
      p.dir.x = p.pos.x - x_src;
      p.dir.y = p.pos.y - y_src;
      p.dir.z = p.pos.z - z_src;
      real src_point_dist =
           sqrt(p.dir.x*p.dir.x + p.dir.y*p.dir.y + p.dir.z*p.dir.z);
      p.dir.x /= src_point_dist;
      p.dir.y /= src_point_dist;
      p.dir.z /= src_point_dist;
   }

   if ( modifier->transport(p,rndm) )
   {
      // trace particle to the simulation grid (take beam angles into account)
      // return true if the particle hits the calculation cube
      return( trace2cube(p.pos, p.dir, p.i, origin_point_dist, rndm) );
   }

   // there isn't any particle to emit, return false
   return(false);
}

// get particle type, particle energy and set flag if the
// particle (electron) is from the primary source
void e_beam_bumpy_mono::emit(particle_type &type,    real   &energy,
                             bool &primary_particle, ranmar &rndm)
{
   // set particle type
   type   = ELECTRON;

   // determine the source particle type
   if (rndm.number() > sct_weight)
   {
      // this is an electron from the primary sources
      primary_particle = true;

      // choose electron energy
      energy = sample_energy(rndm);
   }
   else
   {
      // this is an electron from the scatter source
      primary_particle = false;

      // choose electron energy
      energy = energy_min + (energy_max-energy_min)*rndm.number();
   }
   return;
}

// get particle weight, starting position, direction and voxel index
bool e_beam_bumpy_mono::emit(real   &weight,
                             real_3 &pos, real_3 &dir, int_3 &i,
                             bool   &primary_particle,
#ifdef USE_SOBOL
                             sobseq &sobol, int &sobol_count,
#endif
                             ranmar &rndm)
{
   const real EPSI = 1.0e-20;    // a small parameter
   real origin_point_dist;       // origin to point distance

   // set particle weight
   weight = ONE;

   if (primary_particle)
   {
      // this is an electron from the primary sources

      // random (Sobol) numbers to sample initial position
#ifdef USE_SOBOL
      real rx = sobol.number(sobol_count); ++sobol_count;
      real ry = sobol.number(sobol_count); ++sobol_count;
#else
      real rx = rndm.number();
      real ry = rndm.number();
#endif 

      // electron position in the applicator plane
      pos.x  = app_x1 + app_width_x*rx;
      pos.y  = app_y1 + app_width_y*ry;
      pos.z  = app_distance;

      // sample fluctuation of the electron starting position
      // in the applicator plane from Gaussian distribution
#ifdef USE_SOBOL
      real rr   = sobol.number(sobol_count); ++sobol_count;
      real rphi = sobol.number(sobol_count); ++sobol_count;
#else
      real rr   = rndm.number();
      real rphi = rndm.number();
#endif
      real rad_src =
        sqrt(TWO*pri_xo2*max_of(ZERO,real(-log(max_of(EPSI,rr)))));
      real phi_src = TWO_PI*rphi;
      pos.x += rad_src*cos(phi_src);
      pos.y += rad_src*sin(phi_src);

      // apply bump correction
		real radius = sqrt(pos.x*pos.x+pos.y*pos.y);
		real x_radius = radius;
		real y_radius = radius;
		if (pos.x < 0) x_radius = -x_radius;
		if (pos.y < 0) y_radius = -y_radius;	
		
		if (radius == 0.0) {
         weight *= interpolate(n_bump,xy_bump,wx_bump,pos.x);
         weight *= interpolate(n_bump,xy_bump,wy_bump,pos.y);
		} else {
			real thetax = acos(fabs(pos.x)/radius);
			real thetay = acos(fabs(pos.y)/radius);
			real weightx = interpolate(n_bump,xy_bump,wx_bump,x_radius);
			real weighty = interpolate(n_bump,xy_bump,wy_bump,y_radius);
			weight *= ((thetay*weightx+thetax*weighty)/(thetax+thetay));
		}
		
      // origin to point distance
      origin_point_dist = sqrt(pos.x*pos.x + pos.y*pos.y + pos.z*pos.z);

      // calculate preferential moving direction and normalize
      dir.x = pos.x/origin_point_dist;
      dir.y = pos.y/origin_point_dist;
      dir.z = pos.z/origin_point_dist;

      // sample fluctuations (d_theta, d_phi) of the moving direction
      real d_theta = ZERO;
      if (rndm.number() <= pri_weight1)
      {
         d_theta =  sigma_theta_x1*SQRT2
                   *sqrt(max_of(ZERO,real(-log(max_of(EPSI,rndm.number())))));
      }
      else
      {
         d_theta =  sigma_theta_x2*SQRT2
                   *sqrt(max_of(ZERO,real(-log(max_of(EPSI,rndm.number())))));
      }
      real d_phi = TWO_PI*rndm.number();

      // rotate electron direction by d_theta and d_phi
      rotate(d_theta, d_phi, dir);
   }
   else
   {
      // this is an electron from the scatter source

      // random (Sobol) numbers to sample initial position
#ifdef USE_SOBOL
      real rx = sobol.number(sobol_count); ++sobol_count;
      real ry = sobol.number(sobol_count); ++sobol_count;
#else
      real rx = rndm.number();
      real ry = rndm.number();
#endif 

      // electron position in the applicator plane
      pos.x  = app_x1 + app_width_x*rx;
      pos.y  = app_y1 + app_width_y*ry;
      pos.z  = app_distance;

      // origin to point distance
      origin_point_dist = sqrt(pos.x*pos.x + pos.y*pos.y + pos.z*pos.z);

      // sample x position in the applicator scatter plane
      real x_src = ZERO;
      if (rndm.number() < rx) {
         x_src = app_x1 + (pos.x-app_x1)*sqrt(rndm.number()); }
      else {
         x_src = app_x2 + (pos.x-app_x2)*sqrt(rndm.number()); }

      // sample y position in the applicator scatter plane
      real y_src = ZERO;
      if (rndm.number() < ry) {
         y_src = app_y1 + (pos.y-app_y1)*sqrt(rndm.number()); }
      else {
         y_src = app_y2 + (pos.y-app_y2)*sqrt(rndm.number()); }

      // z coordinate is given by the position of scatter source plane
      real z_src = app_distance - sct_app_dist;

      // calculate moving direction and normalize
      dir.x = pos.x - x_src;
      dir.y = pos.y - y_src;
      dir.z = pos.z - z_src;
      real src_point_dist = sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
      dir.x /= src_point_dist;
      dir.y /= src_point_dist;
      dir.z /= src_point_dist;
   }

   if ( modifier->transport(pos,dir,rndm) )
   {
      // trace particle to the simulation grid (take beam angles into account)
      // return true if the particle hits the calculation cube
      return( trace2cube(pos, dir, i, origin_point_dist, rndm) );
   }

   // there isn't any particle to emit, return false
   return(false);
}
