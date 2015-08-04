/*****************************************************************************
 * p_beam_2gauss_mono.cpp:                                                   *
 *    class member functions for:                                            *
 *       p_beam_2gauss_mono: two Gaussian sources, mono-energetic photons    *
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
#include <iostream>

#include <sstream>
#include "p_beam_2gauss_mono.h"

// ***********************************************
// member functions of class p_beam_2gauss_mono
// ***********************************************

// base data file input
void p_beam_2gauss_mono::get_base_data()
{
   // define bool variables
   bool p_pri_found         = false;
   bool pri_distance_found  = false;
   bool sigma_pri_found     = false;
   bool horn_0_found        = false;
   bool horn_1_found        = false;
   bool horn_2_found        = false;
   bool horn_3_found        = false;
   bool horn_4_found        = false;
   bool sct_distance_found  = false;
   bool sigma_sct_found     = false;
   bool col_mdistance_found = false;
   bool col_cdistance_found = false;
   bool col_xdistance_found = false;
   bool col_ydistance_found = false;
   bool norm_value_found    = false;
   bool gray_mu_dmax_found  = false;
   bool energy_min_found    = false;
   bool energy_max_found    = false;
   bool l_value_found       = false;
   bool b_value_found       = false;
   bool p_con_found         = false;
   bool distance_con_found  = false;
   bool radius_con_found    = false;
   bool e_mean_con_found    = false;
#ifdef USE_NU
   bool nu_value_found      = false;
	nu_value = 0.45;         // Defualt Value
#endif

   // read lines
   char line[81]  = "";    // lines to read from base data file
   bool read_line = true;
   while (read_line)
   {
      if (base_file.eof())
      {
         xvmc_error("p_beam_2gauss_mono::get_base_data",
                    "end of file reached, base data incomplete",8);
      }
      base_file.getline(line,sizeof(line));
      istringstream line_stream(line);
      char keyword[81] = "";
      line_stream >> keyword;
      if (!strcmp(keyword,"END-PARAMETERS"))
      {
         xvmc_error("p_beam_2gauss_mono::get_base_data",
                    "end of parameter entry, base data incomplete",8);
      }
      else
      {
         if (!strcmp(keyword,"PRIMARY-PHOTONS:"))
         {
            line_stream >> p_pri;
            p_sct = ONE-p_pri;
            p_pri_found = true;
         }
         if (!strcmp(keyword,"PRIMARY-DIST:"))
         {
            line_stream >> pri_distance;
            pri_distance_found = true;
         }
         if (!strcmp(keyword,"PRIMARY-SIGMA:"))
         {
            line_stream >> sigma_pri;
            sigma_pri_found = true;
         }
         if (!strcmp(keyword,"PRIMARY-HORN0:"))
         {
            line_stream >> horn_0;
            horn_0_found = true;
         }
         if (!strcmp(keyword,"PRIMARY-HORN1:"))
         {
            line_stream >> horn_1;
            horn_1_found = true;
         }
         if (!strcmp(keyword,"PRIMARY-HORN2:"))
         {
            line_stream >> horn_2;
            horn_2_found = true;
         }
         if (!strcmp(keyword,"PRIMARY-HORN3:"))
         {
            line_stream >> horn_3;
            horn_3_found = true;
         }
         if (!strcmp(keyword,"PRIMARY-HORN4:"))
         {
            line_stream >> horn_4;
            horn_4_found = true;
         }
         if (!strcmp(keyword,"SCATTER-DIST:"))
         {
            line_stream >> sct_distance;
            sct_distance_found = true;
         }
         if (!strcmp(keyword,"SCATTER-SIGMA:"))
         {
            line_stream >> sigma_sct;
            sigma_sct_found = true;
         }
         if (!strcmp(keyword,"COLM-DIST:"))
         {
            line_stream >> col_mdistance;
            col_mdistance_found = true;
         }
         if (!strcmp(keyword,"COLC-DIST:"))
         {
            line_stream >> col_cdistance;
            col_cdistance_found = true;
         }
         if (!strcmp(keyword,"COLX-DIST:"))
         {
            line_stream >> col_xdistance;
            col_xdistance_found = true;
         }
         if (!strcmp(keyword,"COLY-DIST:"))
         {
            line_stream >> col_ydistance;
            col_ydistance_found = true;
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
            line_stream >> energy_min_aux;
            energy_min_found = true;
         }
         if (!strcmp(keyword,"ENERGY-MAX:"))
         {
            line_stream >> energy_max_aux;
            energy_max_found = true;
         }
         if (!strcmp(keyword,"L-VALUE:"))
         {
            line_stream >> l_aux;
            l_value_found = true;
         }
         if (!strcmp(keyword,"B-VALUE:"))
         {
            line_stream >> b_aux;
            b_value_found = true;
         }
         if (!strcmp(keyword,"CHARGED-PARTICLES:"))
         {
            line_stream >> p_con;
            p_con_found = true;
         }
         if (!strcmp(keyword,"CHARGED-DIST:"))
         {
            line_stream >> distance_con;
            distance_con_found = true;
         }
         if (!strcmp(keyword,"CHARGED-RADIUS:"))
         {
            line_stream >> radius_con;
            radius_con_found = true;
         }
         if (!strcmp(keyword,"CHARGED-E-MEAN:"))
         {
            line_stream >> e_mean_con;
            e_mean_con_found = true;
         }
#ifdef USE_NU
         if (!strcmp(keyword,"NU-VALUE:"))
         {
            line_stream >> nu_value;
            nu_value_found = true;
         }
#endif
      } // if (!strcmp(keyword,"END-PARAMETERS"))

      if ( p_pri_found         &&
           pri_distance_found  &&
           sigma_pri_found     &&
           horn_0_found        &&
           horn_1_found        &&
           horn_2_found        &&
           horn_3_found        &&
           horn_4_found        &&
           sct_distance_found  &&
           sigma_sct_found     &&
           col_mdistance_found &&
           col_cdistance_found &&
           col_xdistance_found &&
           col_ydistance_found &&
           norm_value_found    &&
           gray_mu_dmax_found  &&
           energy_min_found    &&
           energy_max_found    &&
           l_value_found       &&
           b_value_found       &&
           p_con_found         &&
           distance_con_found  &&
           radius_con_found    &&
#ifdef USE_NU_LATER
           nu_value_found      &&
#endif
           e_mean_con_found       ) read_line = false;

   } // while (read_line)

   // Hyperion sets the rectangular field outline through  real set_score_rectangle(real &, real &, real &, real &);

   // parameter input was successful, now reset collimator openings
   col_x1       = open_x1*col_xdistance/iso_distance;
   col_x2       = open_x2*col_xdistance/iso_distance;
   col_y1       = open_y1*col_ydistance/iso_distance;
   col_y2       = open_y2*col_ydistance/iso_distance;
   col_width_x  = col_x2 - col_x1;
   col_width_y  = col_y2 - col_y1;

   // set auxiliary charged particle spectrum parameters
   e1_con_aux = exp(-e_cut/e_mean_con);
   e2_con_aux = e1_con_aux - exp(-energy_max_aux/e_mean_con);

   // print beam model parameters
   xvmc_message("p_pri:        ",p_pri*100.0,"%",1);
   xvmc_message("pri_distance: ",pri_distance,"cm",0);
   xvmc_message("sigma_pri:    ",sigma_pri,"cm",0);
   xvmc_message("horn_0:       ",horn_0,"",0);
   xvmc_message("horn_1:       ",horn_1,"",0);
   xvmc_message("horn_2:       ",horn_2,"",0);
   xvmc_message("horn_3:       ",horn_3,"",0);
   xvmc_message("horn_4:       ",horn_4,"",0);
   xvmc_message("p_sct:        ",p_sct*100.0,"%",0);
   xvmc_message("sct_distance: ",sct_distance,"cm",0);
   xvmc_message("sigma_sct:    ",sigma_sct,"cm",0);
   xvmc_message("col_mdistance:",col_mdistance,"cm",0);
   xvmc_message("col_cdistance:",col_cdistance,"cm",0);
   xvmc_message("col_xdistance:",col_xdistance,"cm",0);
   xvmc_message("col_ydistance:",col_ydistance,"cm",0);
   xvmc_message("norm_value:   ",norm_value,"",0);
   xvmc_message("gray_mu_dmax: ",gray_mu_dmax,"Gy",0);
   xvmc_message("energy_min:   ",energy_min_aux,"MeV",0);
   xvmc_message("energy_max:   ",energy_max_aux,"MeV",0);
   xvmc_message("l:            ",l_aux,"",0);
   xvmc_message("b:            ",b_aux,"",0);
   xvmc_message("p_con:        ",p_con*100.0,"%",0);
   xvmc_message("distance_con: ",distance_con,"cm",0);
   xvmc_message("radius_con:   ",radius_con,"cm",0);
   xvmc_message("e_mean_con:   ",e_mean_con,"MeV",0);
#ifdef USE_NU
	xvmc_message("nu_value:     ",nu_value,"",0);
#endif
}

// get particle parameters from photon beam with two Gaussian sources
bool p_beam_2gauss_mono::emit(particle_parameters &p,
#ifdef USE_SOBOL
                              sobseq &sobol, int &sobol_count,
#endif
                              ranmar &rndm)
{
   real origin_point_dist; // origin (x=y=z=0) to point distance
                           // (position in X-collimator plane)
   real rx,ry,rp,rr,rphi;  // auxiliary random (Sobol) numbers
   real src_point_dist;    // distance from particle position in the source
                           // plane to position in collimator plane (point)
   real rad_src;           // particle position in source plane (radius)
   real phi_src;           // particle position in source plane (angle)
   real x_src,y_src,z_src; // particle position in source plane (x,y,z)

   real energy_0;          // photon energy before the softening correction

   real rad_pri;           // primary photon position in target plane (radius)
   real phi_pri;           // primary photon position in target plane (angle)
   real x_pri,y_pri,z_pri; // primary photon direction before filter scatter
   real pri_sct_distance;  // distance from primary photon position to the
                           // photon position in the head scatter plane
                           // (simulate Compton)
   real cos_comp;          // Compton scattering angle
   real fac_comp;          // Compton correction factor

   // first of all, we check for charged particle (electron) contamination

   // random (Sobol) number to determine the source particle type
#ifdef USE_SOBOL
   rp = sobol.number(sobol_count); ++sobol_count;
#else
   rp = rndm.number();
#endif
   if (rp <= p_con)
   {
      // this is an electron
      p.type   = ELECTRON;

      // position in electron source plane
#ifdef USE_SOBOL
      rr   = sobol.number(sobol_count); ++sobol_count;
      rphi = sobol.number(sobol_count); ++sobol_count;
#else
      rr = rndm.number(); rphi = rndm.number();
#endif
      // sample initial position from a uniform distribution
      rad_src = radius_con*sqrt(max_of(ZERO,rr));
      phi_src = TWO_PI*rphi;
      x_src   = rad_src*cos(phi_src);
      y_src   = rad_src*sin(phi_src);
      z_src   = distance_con;

      // random (Sobol) numbers to sample initial position
      // in the upper collimator plane
#ifdef USE_SOBOL
      rx = sobol.number(sobol_count); ++sobol_count;
      ry = sobol.number(sobol_count); ++sobol_count;
#else
      rx = rndm.number();
      ry = rndm.number();
#endif

      // sample electron position in the MC starting plane
      // (above the collimators)
      p.pos.z = col_mdistance;

      // X-position in the X-collimator plane
      p.pos.x = col_x1 + col_width_x*rx;
      // X-position in the MC starting plane (above the collimators)
      p.pos.x =
         (p.pos.x-x_src)*(col_mdistance-z_src)/(col_xdistance-z_src)+x_src;

      // Y-position in the Y-collimator plane
      p.pos.y = col_y1 + col_width_y*ry;
      // Y-position in the MC starting plane (above the collimators)
      p.pos.y =
         (p.pos.y-y_src)*(col_mdistance-z_src)/(col_ydistance-z_src)+y_src;

      // origin to point distance
      origin_point_dist =
         sqrt(p.pos.x*p.pos.x + p.pos.y*p.pos.y + p.pos.z*p.pos.z);

      // direction
      p.dir.x = p.pos.x - x_src;
      p.dir.y = p.pos.y - y_src;
      p.dir.z = p.pos.z - z_src;
      src_point_dist =
         sqrt(p.dir.x*p.dir.x + p.dir.y*p.dir.y + p.dir.z*p.dir.z);
      p.dir.x = p.dir.x/src_point_dist;
      p.dir.y = p.dir.y/src_point_dist;
      p.dir.z = p.dir.z/src_point_dist;

      // weight
      p.weight = ONE;

      // sample energy from exponential distribution
      const real EPSI = 1.0e-20;
      p.energy =
         -e_mean_con*log(max_of(EPSI,e1_con_aux-rndm.number()*e2_con_aux));
   }
   else
   {
      // this is a photon, set photon parameters
      p.type = PHOTON;

      // check p_con to avoid division by zero
      if (p_con >= ONE)
         xvmc_error("p_beam_2gauss_mono::emit",
                    "cannot emit photon if contamination probabilty >= 1",8);

      // scale random number to sample primary or head scatter photon
      rp = (rp-p_con)/(ONE-p_con);
      if (rp < p_sct)
      {
         /* simulate photon head scatter */

         // position in head scatter plane
#ifdef USE_SOBOL
         rr   = sobol.number(sobol_count); ++sobol_count;
         rphi = sobol.number(sobol_count); ++sobol_count;
#else
         rr = rndm.number(); rphi = rndm.number();
#endif
         const real EPSI = 1.0e-20;
         rad_src = sigma_sct*sqrt(TWO*max_of(ZERO,real(-log(max_of(EPSI,rr)))));
         phi_src = TWO_PI*rphi;
         x_src   = rad_src*cos(phi_src);
         y_src   = rad_src*sin(phi_src);
         z_src   = sct_distance;

         // random (Sobol) numbers to sample initial position
#ifdef USE_SOBOL
         rx = sobol.number(sobol_count); ++sobol_count;
         ry = sobol.number(sobol_count); ++sobol_count;
#else
         rx = rndm.number();
         ry = rndm.number();
#endif

         // sample photon position in the MC starting plane
         // (above the collimators)
         p.pos.z = col_mdistance;

         // X-position in the X-collimator plane
#ifndef MOVING_JAW
         p.pos.x = col_x1 + col_width_x*rx;
#else
         p.pos.x = col_x1 + col_width_x - col_width_x*rx*rndm.number();
#endif
         // X-position in the MC starting plane (above the collimators)
         p.pos.x =
            (p.pos.x-x_src)*(col_mdistance-z_src)/(col_xdistance-z_src)+x_src;

         // Y-position in the Y-collimator plane
         p.pos.y = col_y1 + col_width_y*ry;
         // Y-position in the MC starting plane (above the collimators)
         p.pos.y =
            (p.pos.y-y_src)*(col_mdistance-z_src)/(col_ydistance-z_src)+y_src;

         // origin to point distance
         origin_point_dist =
            sqrt(p.pos.x*p.pos.x + p.pos.y*p.pos.y + p.pos.z*p.pos.z);

         // direction
         p.dir.x = p.pos.x - x_src;
         p.dir.y = p.pos.y - y_src;
         p.dir.z = p.pos.z - z_src;
         src_point_dist =
            sqrt(p.dir.x*p.dir.x + p.dir.y*p.dir.y + p.dir.z*p.dir.z);
         p.dir.x = p.dir.x/src_point_dist;
         p.dir.y = p.dir.y/src_point_dist;
         p.dir.z = p.dir.z/src_point_dist;

         // energy
         energy_0 = sample_energy(rndm);

         //
         // simulate Compton scatter in the flattening filter
         //

         // primary position in the target plane
         rr = rndm.number(); rphi = rndm.number();
         rad_pri = sigma_pri*sqrt(TWO*max_of(ZERO,real(-log(max_of(EPSI,rr)))));
         phi_pri = TWO_PI*rphi;
         x_pri   = rad_pri*cos(phi_pri);
         y_pri   = rad_pri*sin(phi_pri);
         z_pri   = pri_distance;

         // direction before filter scatter
         x_pri = x_src - x_pri;
         y_pri = y_src - y_pri;
         z_pri = z_src - z_pri;

         // normalize
         pri_sct_distance = sqrt(x_pri*x_pri + y_pri*y_pri + z_pri*z_pri);
         x_pri /= pri_sct_distance;
         y_pri /= pri_sct_distance;
         z_pri /= pri_sct_distance;

         // Compton scattering angle (scalar product)
         cos_comp = p.dir.x*x_pri + p.dir.y*y_pri + p.dir.z*z_pri;

         // Compton correction factor
         fac_comp = ONE+(ONE-cos_comp)*energy_0/EMASS;

         // adjust energy
         p.energy = energy_0/fac_comp;

         // weight
         p.weight = fac_comp;

         // now the brass build-up cap correction:
         p.weight *= (pow(double(energy_0),-0.558)+0.026*energy_0)/
                     (pow(double(p.energy),-0.558)+0.026*p.energy);
      }
      else
      {
         /* simulate primary (target) photons */

         // position in target plane
#ifdef USE_SOBOL
         rr   = sobol.number(sobol_count); ++sobol_count;
         rphi = sobol.number(sobol_count); ++sobol_count;
#else
         rr = rndm.number(); rphi = rndm.number();
#endif
         const real EPSI = 1.0e-20;
         rad_src = sigma_pri*sqrt(TWO*max_of(ZERO,real(-log(max_of(EPSI,rr)))));
         phi_src = TWO_PI*rphi;
         x_src   = rad_src*cos(phi_src);
         y_src   = rad_src*sin(phi_src);
         z_src   = pri_distance;

         // random (Sobol) numbers to sample initial position
#ifdef USE_SOBOL
         rx = sobol.number(sobol_count); ++sobol_count;
         ry = sobol.number(sobol_count); ++sobol_count;
#else
         rx = rndm.number();
         ry = rndm.number();
#endif

         // sample photon position in the MC starting plane
         // (above the collimators)
         p.pos.z = col_mdistance;

         // X-position in the X-collimator plane
#ifndef MOVING_JAW
         p.pos.x = col_x1 + col_width_x*rx;
#else
         p.pos.x = col_x1 + col_width_x - col_width_x*rx*rndm.number();
#endif
         // X-position in the MC starting plane (above the collimators)
         p.pos.x =
            (p.pos.x-x_src)*(col_mdistance-z_src)/(col_xdistance-z_src)+x_src;

         // Y-position in the Y-collimator plane
         p.pos.y = col_y1 + col_width_y*ry;
         // Y-position in the MC starting plane (above the collimators)
         p.pos.y =
            (p.pos.y-y_src)*(col_mdistance-z_src)/(col_ydistance-z_src)+y_src;

         // origin to point distance
         origin_point_dist =
            sqrt(p.pos.x*p.pos.x + p.pos.y*p.pos.y + p.pos.z*p.pos.z);

         // direction
         p.dir.x = p.pos.x - x_src;
         p.dir.y = p.pos.y - y_src;
         p.dir.z = p.pos.z - z_src;
         src_point_dist =
            sqrt(p.dir.x*p.dir.x + p.dir.y*p.dir.y + p.dir.z*p.dir.z);
         p.dir.x = p.dir.x/src_point_dist;
         p.dir.y = p.dir.y/src_point_dist;
         p.dir.z = p.dir.z/src_point_dist;

#ifdef CLIPPING
			// Primary Collimator Clipping
         float zpc = 7.6; // For Varian 21EX
         float lpc = 0.0;
         if (p.dir.z > EPSI) lpc=(zpc-z_src)/p.dir.z;
         float xpc = x_src + p.dir.x * lpc;
         float ypc = y_src + p.dir.y * lpc;
         float rpc = sqrt(xpc*xpc+ypc*ypc);
         if (rpc < zpc*tan(14*ONE_PI/180)){
#endif
         /* estimate the horn effect (primary photons only) */

         // calculate normalized radius
         real rad2_horn = (p.pos.x*p.pos.x+p.pos.y*p.pos.y)/(p.pos.z*p.pos.z);
         real rad1_horn = sqrt(rad2_horn);
         real rad3_horn = rad1_horn*rad2_horn;
         real rad4_horn = rad2_horn*rad2_horn;
         real cor_horn  =
            ONE + rad2_horn*(horn_0 + horn_1*rad1_horn + horn_2*rad2_horn
                                    + horn_3*rad3_horn + horn_4*rad4_horn);
         const real CUT_HORN = 0.5;
         if (cor_horn < CUT_HORN) cor_horn = CUT_HORN;

         /* estimate off-axis softening for primary photons */

         // angle of the position vector to central axis in degrees
         double theta_oas = p.pos.z/origin_point_dist;
         theta_oas = min_of(theta_oas,1.0);
         theta_oas = max_of(theta_oas,-1.0);
         theta_oas = acos(theta_oas)*180.0/ONE_PI;
#ifdef VARIAN
         theta_oas = min_of(14.00,theta_oas); // maximum of mu_cor

			double radius_oas = 100.0*tan(theta_oas*ONE_PI/180);
			// off-axis correction (Yu et al 1997, Med. Phys. 24, pp. 233-239)
         real mu_cor =
            ONE + (0.0028952772+2.4845996E-4*radius_oas)*radius_oas;
#else
         theta_oas = min_of(14.73,theta_oas); // maximum of mu_cor

// Moved by JOKim on Jan 14, 2008
// sample energy and correct due to off-axis softening, i.e.
// in the mono-energetic case the beam is not really mono-energetic
energy_0 = sample_energy(rndm);

  #ifndef NEWOFFAXIS
         // off-axis correction (Tailor et al 1998, Med. Phys. 25, pp. 662-667)
			real mu_cor =
            ONE + (0.00181+(0.00202-0.0000942*theta_oas)*theta_oas)*theta_oas;
  #else

		   real E1 = energy_0;
			if (E1 < 0.75) E1 = 0.75;  //  Supress Low Energy at Off-Axis
			if (E1 > 6.0) E1 = 6.0;
		   real E2 = E1*E1;
		   real E3 = E1*E2;
		   real E4 = E2*E2;

		/*	O3 Fitted
		   real A1 = 7.2552E-4*E4 - 9.7528E-3*E3 + 4.5867E-2*E2 - 9.9351E-2*E1 + 6.8884E-2;
		   real A2 = 9.0826E-5*E4 - 1.6258E-3*E3 + 9.8449E-3*E2 - 2.3305E-2*E1 + 1.4884E-2;
		   real A3 = 2.1611E-5*E3 - 2.2567E-4*E2 + 7.1899E-4*E1 - 5.4195E-4;
		   real S = ((A3*theta_oas + A2)*theta_oas +A1)*theta_oas + ONE;
      		*/

		/* O4 Fitted
		   real A1 = 6.7129E-4*E4 - 8.0021E-3*E3 + 3.2132E-2*E2 - 5.6897E-2*E1 + 3.4334E-2;
		   real A2 =-9.9772E-4*E3 + 1.0326E-2*E2 - 3.3945E-2*E1 + 2.7091E-2;
		   real A3 = 7.8077E-5*E3 - 8.0101E-4*E2 + 2.7207E-3*E1 - 2.2578E-3;
		   real A4 =-2.1047E-6*E3 + 2.1445E-5*E2 - 7.4609E-5*E1 + 6.3957E-5;
		   real S = (((A4*theta_oas + A3)*theta_oas +A2)*theta_oas + A1)*theta_oas+ONE;
      		*/
		// O4 Fitted
		   real A1 = 1.1382E-3*E4 - 1.4172E-2*E3 + 6.0892E-2*E2 - 1.1337E-1*E1 + 6.7807E-2;
		   real A2 =-1.3322E-4*E4 + 8.9809E-4*E3 + 5.0240E-4*E2 - 1.1892E-2*E1 + 1.2362E-2;
		   real A3 = 3.2779E-5*E4 - 3.5839E-4*E3 + 1.2502E-3*E2 - 1.3286E-3*E1 + 1.4793E-4;
		   real A4 =-1.4961E-6*E4 + 1.7744E-5*E3 - 7.1306E-5*E2 + 1.0692E-4*E1 - 4.2867E-5;
		   real S = (((A4*theta_oas + A3)*theta_oas +A2)*theta_oas + A1)*theta_oas+ONE;

  #endif
#endif
         //real mu_cor2 = mu_cor*mu_cor; // mu_cor squared
#ifdef USE_NU
			// printf ("New Nu = %f\n", NEW_NU);
			if (nu_value == 0.0) xvmc_error("p_beam_2gauss_mono::emit","nu_value is zero",8);
			// nu_value = 1.0;
         // cout << " >>>>>> nu_value = " << nu_value * 1.0 << endl;
			real mu_cor2 = pow(double(mu_cor),1.0/nu_value); // mu_cor squared
#else
  #ifndef NEWOFFAXIS
         real mu_cor2 = pow(double(mu_cor),2.22222); // mu_cor squared
  #else
         real mu_cor2 = 1.0/S;
  #endif
#endif
         p.energy = energy_0/mu_cor2;

         // change weight due to off-axis softening and horn correction,
         // the energy fluence should be proportional to cor_horn (dose
         // in air) divided by a correction caused by the energy dependence
         // of the brass attenuation coeffcient
         p.weight  = mu_cor2*cor_horn;

         // now the brass build-up cap correction:
         // p.weight /= pow(mu_cor2,0.628)*exp(0.1*p.energy*(ONE-mu_cor2));
         // p.weight /= pow(mu_cor2,0.584)*exp(0.077*p.energy*(ONE-mu_cor2));
         // p.weight /= pow(mu_cor2,0.584);
         p.weight *= (pow(double(energy_0),-0.558)+0.026*energy_0)/
                     (pow(double(p.energy),-0.558)+0.026*p.energy);

			// gprintf ("%e %e\n", energy_0, p.energy);

#ifdef NEW_NU
			p.weight *= pow(double(energy_0),-NEW_NU/2.2);
#endif

#ifdef CLIPPING
         } else {
         	// Ignore the particle transport because of primary collimator clip
         	p.energy = EPSI;
         	p.weight = EPSI;
         	// printf ("Particle Ignored\n");
			}
#endif
      }
#ifdef SPECTRUM_WEIGHT
		float spect[26] = {0.3057, 1.2322, 1.1644, 1.0258, 0.9752,
								 0.9194, 0.9007, 0.8963, 0.9179, 0.9458,
								 0.9537, 1.0307, 1.0891, 1.1496, 1.2207,
								 1.2639, 1.3803, 1.4861, 1.5523, 1.6304,
								 1.6090, 1.5982, 1.4607, 1.1091, 0.4886,
								 0.0878};
		float eMin = 0.25;
		float eMax = 6.50;
		float delE = 0.25;

		int iE = (int)((energy_0 - eMin)/delE);
		if (energy_0 < eMin) iE = 0;
		if (energy_0 > eMax) iE = 25;

		float sWeight = 1.0;
		if (iE >= 0 && iE < 25) {
			sWeight = spect[iE] + (spect[iE+1] - spect[iE])/delE *(energy_0-(eMin + iE*delE));
		}
		// printf ("%e %d %e %f %e\n", energy_0, iE, p.weight, sWeight, p.weight*sWeight);
		if (sWeight > 0.0) p.weight *= sWeight;
#endif

   }

/*
	printf ("b  %8d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %4d %4d %4d\n",
	     p.type,  p.energy, p.weight,
	     p.pos.x, p.pos.y,  p.pos.z,
	     p.dir.x, p.dir.y,  p.dir.z,
	     p.i.x,   p.i.y,    p.i.z);
*/
   if ( modifier->transport(p,rndm) )
   {
/*
	printf ("m  %8d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %4d %4d %4d\n",
	     p.type,  p.energy, p.weight,
	     p.pos.x, p.pos.y,  p.pos.z,
	     p.dir.x, p.dir.y,  p.dir.z,
	     p.i.x,   p.i.y,    p.i.z);
*/
      // trace particle to the simulation grid (take beam angles into account)
      // return true if the particle hits the calculation cube
      return( trace2cube(p.pos, p.dir, p.i, origin_point_dist, rndm) );
   }

   // there isn't any particle to emit, return false
   return(false);
}

p_beam_2gauss_mono::p_beam_2gauss_mono(const p_beam_2gauss_mono &A)
{
  // beam_model variables
  iso_distance = A.iso_distance;
  type = A.type;
  nominal_energy = A.nominal_energy;
  model_id = A.model_id;
  // p_beam_model variables
  col_mdistance = A.col_mdistance;
  col_cdistance = A.col_cdistance;
  col_xdistance = A.col_xdistance;
  col_ydistance = A.col_ydistance;
  col_x1 = A.col_x1;
  col_x2 = A.col_x2;
  col_y1 = A.col_y1;
  col_y2 = A.col_y2;
  col_width_x = A.col_width_x;
  col_width_y = A.col_width_y;
  //p_beam_2gauss_mono variables
  p_pri = A.p_pri;
  pri_distance = A.pri_distance;
  sigma_pri = A.sigma_pri;
  horn_0 = A.horn_0;
  horn_1 = A.horn_1;
  horn_2 = A.horn_2;
  horn_3 = A.horn_3;
  horn_4 = A.horn_4;
  p_sct = A.p_sct;
  sct_distance = A.sct_distance;
  sigma_sct = A.sigma_sct;
  norm_value = A.norm_value;
  gray_mu_dmax = A.gray_mu_dmax;
  energy_min_aux = A.energy_min_aux;
  energy_max_aux = A.energy_max_aux;
  l_aux = A.l_aux;
  b_aux = A.b_aux;
#ifdef USE_NU
  nu_value = A.nu_value;
#endif
  p_con = A.p_con;
  distance_con = A.distance_con;
  radius_con = A.radius_con;
  e_mean_con = A.e_mean_con;
}
