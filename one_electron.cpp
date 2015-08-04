/*****************************************************************************
 * one_electron.cpp:                                                         *
 *    function:                                                              *
 *       one_electron: simulate single electron history (no repetition)      *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 00/02/24        *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <math.h>
#include "definitions.h"
#include "global.h"
#include "ranmar.h"
#include "xvmc_util.h"
#include "electron_data.h"
#include "moller_XS.h"
#include "brems_XS.h"
#include "portal_dose.h"

// ****************************************
// declare functions and global variables
// ****************************************

// water electron transport data for this beam
extern electron_data *e_h2o;
// Moller cross section for water
extern moller_XS h2o_moller;
// bremsstrahlung cross section for water
extern brems_XS h2o_brems;

// sample multiple scattering angle
bool mscat(const real &, const real &, real &, real &, ranmar &);

// simulate bremsstrahlung photons by KERMA approximation
void kerma_photon(particle_parameters &, ranmar &,
#ifdef CHECK_ENERGY
                  sum_energy_type &,
#endif // CHECK_ENERGY
                  array_3d<double> *, portal_dose *);


// --- Added by JKim 15Nov2010 ------------------------------------

// simulate bremsstrahlung photons by KERMA approximation
void kerma_photon_portal(particle_parameters &, ranmar &,
#ifdef CHECK_ENERGY
                  sum_energy_type &,
#endif // CHECK_ENERGY
                  array_3d<double> *);

// --- End of Adding ----------------------------------------------

// create and trace one electron history
void one_electron(particle_parameters &e, ranmar &rndm,
#ifdef CHECK_ENERGY
                  sum_energy_type &sum_energy,
#endif // CHECK_ENERGY
                  array_3d<double> *batch_dose, portal_dose *portal)

{
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
   if ( (e.i.x < 0) || (e.i.x >= dim.x) ||
        (e.i.y < 0) || (e.i.y >= dim.y) ||
        (e.i.z < 0) || (e.i.z >= dim.z) )
   {
      // count portal image dose
      if (portal != NULL) portal->add(e, rndm); // Added by JOKim 29Oct2011
#ifdef CHECK_ENERGY
      // test energy conservation
      sum_energy.eloss += e.energy*e.weight;
#endif // CHECK_ENERGY
      return;
   }

   next_step =  true;
   num_mfp   = -1.0;
   sub_sh2o  = ZERO;
   sub_loss  = ZERO;

   while (next_step)
   {
   // generate one step

      if (e.energy <= e_cut)
      {
#ifdef ELECTRON_TRACK_FULL_STOP
         // deposit the remaining energy in the present voxel
         batch_dose->matrix[e.i.x][e.i.y][e.i.z] += e.energy*e.weight;
#ifdef CHECK_ENERGY
         // test energy conservation
         sum_energy.edepo += e.energy*e.weight;
#endif // CHECK_ENERGY
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
      else
      {
      // begin: e.energy > e_cut
      //   cout << e.energy << endl;
         // determine energy index
         i_e=int((e.energy-e_cut)*e_h2o->inverse_delta_e);
			// cout << "e.energy = " << e.energy << " " << e_cut << " " <<e_h2o->inverse_delta_e << endl;
         // cout << "i_e = " << i_e << endl;;
         // determine number of mean free paths to the next discrete interaction
         // using the fictitious interaction method (see SLAC-265, p. 16)
         if (num_mfp < ZERO)
         {
            num_mfp = ZERO;
            while (true)
            {
               zeta                   = rndm.number();
               if (zeta == ZERO) zeta = rndm.number();
               num_mfp   -= log(ONE-zeta);
               eloss      = num_mfp/e_h2o->sigma_max;
               new_energy = e.energy - eloss;
               if (new_energy > e_cut)
               {
                  new_i_e = int((new_energy-e_cut)*e_h2o->inverse_delta_e);
                  if (rndm.number() <= e_h2o->sigma_tot[new_i_e]) break;
               }
               else
               {
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
         if (eloss > eloss_max)
         {
            // the maximum step size is reached
            eloss                = eloss_max;
            discrete_interaction = false;
         }
         else
         {
            // there is a discrete interaction at the end of the step
            discrete_interaction = true;
         }

         // determine and check "new_energy"
         new_energy = e.energy - eloss;
			// cout << "eloss = " << eloss << "  " << eloss_max << endl;
			// cout << "new_energy = " << new_energy << endl;
         if (new_energy <= e_cut)
         {
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
// cout << "voxel_size = " << voxel_size.x << " " << voxel_size.y << " " << voxel_size.z << endl;
// cout << "pos = " << e.pos.x << " " << e.pos.y << " " << e.pos.z << endl;
      do
      {
      // we go back to this point if there is a discrete interaction
      // at the and of the step or the present step is the last step,
      // i.e. "next_step=false", the bool variable "go_back" must be
      // "true" in both cases

         // find the X-, Y-, Z-distance to the next voxel boundary
         // check X-direction
         if (e.dir.x > ZERO)
         {
            istep.x = 1;
            step.x  = (voxel_size.x*double(e.i.x+1)-e.pos.x)/e.dir.x;
         }
         else
         {
            if (e.dir.x < ZERO)
            {
               istep.x = -1;
               step.x  = (voxel_size.x*double(e.i.x)-e.pos.x)/e.dir.x;
            }
            else
            {
               istep.x = 0;
               step.x  = HUGE_STEP;
            }
         }

         // check Y-direction
         if (e.dir.y > ZERO)
         {
            istep.y = 1;
            step.y  = (voxel_size.y*double(e.i.y+1)-e.pos.y)/e.dir.y;
         }
         else
         {
            if (e.dir.y < ZERO)
            {
               istep.y = -1;
               step.y  = (voxel_size.y*double(e.i.y)-e.pos.y)/e.dir.y;
            }
            else
            {
               istep.y = 0;
               step.y  = HUGE_STEP;
            }
         }

         // check Z-direction
         if (e.dir.z > ZERO)
         {
            istep.z = 1;
            step.z  = (voxel_size.z*double(e.i.z+1)-e.pos.z)/e.dir.z;
         }
         else
         {
            if (e.dir.z < ZERO)
            {
               istep.z = -1;
               step.z  = (voxel_size.z*double(e.i.z)-e.pos.z)/e.dir.z;
            }
            else
            {
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
         while (repeat_step)
         {
         // the electron will be traced to the next voxel boundary, if the
         // electron reaches the boundary we repeat this step in the next
         // voxel, the bool variable "repeat_step" is "true" in this case
         // cout << "step = " << step.x << " " << step.y << " " << step.z << endl;
            // find the next voxel
            if ( (step.z < step.x) && (step.z < step.y) )
            {
               voxel_step = step.z;
               step.x -= voxel_step;
               step.y -= voxel_step;
               step.z  = voxel_size.z/fabs(e.dir.z);
               inew.z = e.i.z + istep.z;
            }
            else
            {
               if (step.x < step.y)
               {
                  voxel_step = step.x;
                  step.x  = voxel_size.x/fabs(e.dir.x);
                  step.y -= voxel_step;
                  step.z -= voxel_step;
                  inew.x = e.i.x + istep.x;
               }
               else
               {
                  voxel_step = step.y;
                  step.x -= voxel_step;
                  step.y  = voxel_size.y/fabs(e.dir.y);
                  step.z -= voxel_step;
                  inew.y = e.i.y + istep.y;
               }
            }

            // collision stopping power factor of the present voxel
            f_col = dens_scol->matrix[e.i.x][e.i.y][e.i.z]
                  - dens_ccol->matrix[e.i.x][e.i.y][e.i.z]*s_correction;

            // total stopping power factor of the present voxel
            factor =
            (ONE + alpha*dens_srad->matrix[e.i.x][e.i.y][e.i.z])/(ONE + alpha);
            f_tot = f_col*factor;

            // corresponding step size in water
            voxel_sh2o = voxel_step*f_tot;
// cout << "f_col = " << f_col << " " << factor << " " << voxel_sh2o << " " << total_sh2o << " " << voxel_step << endl;
            if (voxel_sh2o >= total_sh2o)
            {
               // the step ends in this voxel
               repeat_step = false;
               path_step  += total_sh2o/f_tot;
               voxel_depo  = total_loss/factor;

               // put the energy in the present voxel
               batch_dose->matrix[e.i.x][e.i.y][e.i.z] +=
                  voxel_depo*e.weight*dose_factor;

               // MS
               sum_chi_cc2 +=
                  total_loss*dens_fchi->matrix[e.i.x][e.i.y][e.i.z];
#ifdef CHECK_ENERGY
               // test energy conservation
               sum_energy.edepo += voxel_depo*e.weight;
               sum_energy.ploss += (total_loss-voxel_depo)*e.weight;
#endif // CHECK_ENERGY
            }
            else
            {
               // the step continues in the next voxel
               voxel_loss  = total_loss*voxel_sh2o/total_sh2o;
               voxel_depo  = voxel_loss/factor;
               total_loss -= voxel_loss;
               total_sh2o -= voxel_sh2o;
               path_step  += voxel_step;

               // put the energy in the present voxel
               batch_dose->matrix[e.i.x][e.i.y][e.i.z] +=
                  voxel_depo*e.weight*dose_factor;

               // MS
               sum_chi_cc2 +=
                  voxel_loss*dens_fchi->matrix[e.i.x][e.i.y][e.i.z];
#ifdef CHECK_ENERGY
               // test energy conservation
               sum_energy.edepo += voxel_depo*e.weight;
               sum_energy.ploss += (voxel_loss-voxel_depo)*e.weight;
#endif // CHECK_ENERGY
					//cout << "inew = " << inew.x << " " << inew.y << " " << inew.z
					//	   << " " << dim.x << " " << dim.y << " " << dim.z << endl;
               if ( (inew.x < 0) || (inew.x >= dim.x) ||
                    (inew.y < 0) || (inew.y >= dim.y) ||
                    (inew.z < 0) || (inew.z >= dim.z) )
               {
                  // the electron leaves the calculation cube, the history
                  // ends here, test energy conservation

      				// count portal image dose
      				if (portal != NULL) portal->add(e, rndm); // Added by JOKim 29Oct2011
#ifdef CHECK_ENERGY
                  sum_energy.eloss +=
                     (total_loss+sub_loss+new_energy)*e.weight;
#endif // CHECK_ENERGY
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
         if (multiple_scattering)
         {  // sample multiple scattering angle between the two sub-steps

            // momentum squared
            p2      = sub_energy*(sub_energy + TWOxEMASS);
            // beta squared, (v/c)^2
            beta2   = p2/(p2 + EMASSxEMASS);
            // number of elastic scatterings
            omega_0 = BC*path/beta2;
            xr      = CHI_CC2*path/p2/beta2/TWO;
            if (save_eloss > DSMALL) xr *= sum_chi_cc2/save_eloss;

            // sample multiple scattering angle (reduced_angle)
#ifdef DEBUG
            if (mscat(omega_0, xr, cos_t, sin_t, rndm))
            {
               xvmc_warning("one_electron","error in mscat",0);
               xvmc_warning("sub_energy",sub_energy,0);
               xvmc_warning("path",path,0);
               xvmc_warning("fchi",dens_fchi->matrix[e.i.x][e.i.y][e.i.z],0);
            }
#else
            mscat(omega_0, xr, cos_t, sin_t, rndm);
#endif // DEBUG
            phi = TWO_PI * rndm.number();
            sin_phi = sin(phi);
            cos_phi = cos(phi);

            // rotate electron direction by the multiple scattering angle
            rotate(cos_t, sin_t, cos_phi, sin_phi, e.dir);

            if (discrete_interaction)
            {
               go_back = true;
               multiple_scattering = false;
               total_sh2o = sub_sh2o;
               total_loss = sub_loss;
               sub_sh2o   = ZERO;
               sub_loss   = ZERO;
            }

            if (!next_step)
            {
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
      if (discrete_interaction)
      {
         // at first, reset number of mean free paths for the next step
         num_mfp = -1.0;

         // determine the type of discrete interaction (Moller or brems. prod.)
         i_e = int((e.energy-e_cut)*e_h2o->inverse_delta_e);
         if (rndm.number() <= e_h2o->p_moller[i_e])
         {  // no Moller interaction if the energy is too small
            // (in this case we jump to the next MS step)
            if ( h2o_moller.interaction(e.energy,   rndm,
                                        new_energy, cos_t_n, sin_t_n,
                                        d.energy,   cos_t_d, sin_t_d) )
            {
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
               one_electron(d,rndm,
#ifdef CHECK_ENERGY
                            sum_energy,
#endif // CHECK_ENERGY
                            batch_dose, portal);

               // direction change of the primary electron
               cos_phi = -cos_phi;
               sin_phi = -sin_phi;

               // rotate primary electron
               rotate(cos_t_n, sin_t_n, cos_phi, sin_phi, e.dir);

               // "e.energy" is now the energy after the Moller interaction
               e.energy = new_energy;
            }
         } // end of Moller interaction
         else
         { // bremsstrahlung production
           // if the energy is larger than TC (corresponds to AP in EGS4)
            if ( h2o_brems.interaction(e.energy,   rndm,
                                       new_energy, cos_t_n, sin_t_n,
                                       p.energy,   cos_t_p, sin_t_p) )
            {
               // "e.energy" is now the energy after bremsstrahlung production
               e.energy = new_energy;
#ifdef CHECK_ENERGY
               sum_energy.ploss += p.energy*e.weight;
#endif // CHECK_ENERGY

               // play Russian Roulette
               if (rndm.number() < P_ROULETTE)
               {
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
#ifdef CHECK_ENERGY
                  sum_energy.ploss -= p.energy*p.weight;
#endif // CHECK_ENERGY
                  // simulate photon
                  kerma_photon(p,rndm,
#ifdef CHECK_ENERGY
                               sum_energy,
#endif // CHECK_ENERGY
                               batch_dose, portal);
               }
            }

         } // end of bremsstrahlung production

      } // end of discrete interaction

   } // end of this step
   return;
}


// --- Added by JKim 15Nov2010 ------------------------------------
// create and trace one electron history
void one_electron_portal(particle_parameters &e, ranmar &rndm,
#ifdef CHECK_ENERGY
                  sum_energy_type &sum_energy,
#endif // CHECK_ENERGY
                  array_3d<double> *batch_dose_portal)
{
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
   if ( (e.i.x < 0) || (e.i.x >= dim_portal.x) ||
        (e.i.y < 0) || (e.i.y >= dim_portal.y) ||
        (e.i.z < 0) || (e.i.z >= dim_portal.z) )
   {
#ifdef CHECK_ENERGY
      // test energy conservation
      sum_energy.eloss += e.energy*e.weight;
#endif // CHECK_ENERGY
      return;
   }

   next_step =  true;
   num_mfp   = -1.0;
   sub_sh2o  = ZERO;
   sub_loss  = ZERO;

   while (next_step)
   {
   // generate one step

      if (e.energy <= e_cut)
      {
#ifdef ELECTRON_TRACK_FULL_STOP
         // deposit the remaining energy in the present voxel
         batch_dose_portal->matrix[e.i.x][e.i.y][e.i.z] += e.energy*e.weight;
#ifdef CHECK_ENERGY
         // test energy conservation
         sum_energy.edepo += e.energy*e.weight;
#endif // CHECK_ENERGY
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
      else
      {
      // begin: e.energy > e_cut
      //   cout << e.energy << endl;
         // determine energy index
         i_e=int((e.energy-e_cut)*e_h2o->inverse_delta_e);
			// cout << "e.energy = " << e.energy << " " << e_cut << " " <<e_h2o->inverse_delta_e << endl;
         // cout << "i_e = " << i_e << endl;;
         // determine number of mean free paths to the next discrete interaction
         // using the fictitious interaction method (see SLAC-265, p. 16)
         if (num_mfp < ZERO)
         {
            num_mfp = ZERO;
            while (true)
            {
               zeta                   = rndm.number();
               if (zeta == ZERO) zeta = rndm.number();
               num_mfp   -= log(ONE-zeta);
               eloss      = num_mfp/e_h2o->sigma_max;
               new_energy = e.energy - eloss;
               if (new_energy > e_cut)
               {
                  new_i_e = int((new_energy-e_cut)*e_h2o->inverse_delta_e);
                  if (rndm.number() <= e_h2o->sigma_tot[new_i_e]) break;
               }
               else
               {
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
         if (eloss > eloss_max)
         {
            // the maximum step size is reached
            eloss                = eloss_max;
            discrete_interaction = false;
         }
         else
         {
            // there is a discrete interaction at the end of the step
            discrete_interaction = true;
         }

         // determine and check "new_energy"
         new_energy = e.energy - eloss;
			// cout << "eloss = " << eloss << "  " << eloss_max << endl;
			// cout << "new_energy = " << new_energy << endl;
         if (new_energy <= e_cut)
         {
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
// cout << "voxel_size = " << voxel_size_portal.x << " " << voxel_size_portal.y << " " << voxel_size_portal.z << endl;
// cout << "pos = " << e.pos.x << " " << e.pos.y << " " << e.pos.z << endl;
      do
      {
      // we go back to this point if there is a discrete interaction
      // at the and of the step or the present step is the last step,
      // i.e. "next_step=false", the bool variable "go_back" must be
      // "true" in both cases

         // find the X-, Y-, Z-distance to the next voxel boundary
         // check X-direction
         if (e.dir.x > ZERO)
         {
            istep.x = 1;
            step.x  = (voxel_size_portal.x*double(e.i.x+1)-e.pos.x)/e.dir.x;
         }
         else
         {
            if (e.dir.x < ZERO)
            {
               istep.x = -1;
               step.x  = (voxel_size_portal.x*double(e.i.x)-e.pos.x)/e.dir.x;
            }
            else
            {
               istep.x = 0;
               step.x  = HUGE_STEP;
            }
         }

         // check Y-direction
         if (e.dir.y > ZERO)
         {
            istep.y = 1;
            step.y  = (voxel_size_portal.y*double(e.i.y+1)-e.pos.y)/e.dir.y;
         }
         else
         {
            if (e.dir.y < ZERO)
            {
               istep.y = -1;
               step.y  = (voxel_size_portal.y*double(e.i.y)-e.pos.y)/e.dir.y;
            }
            else
            {
               istep.y = 0;
               step.y  = HUGE_STEP;
            }
         }

         // check Z-direction
         if (e.dir.z > ZERO)
         {
            istep.z = 1;
            step.z  = (voxel_size_portal.z*double(e.i.z+1)-e.pos.z)/e.dir.z;
         }
         else
         {
            if (e.dir.z < ZERO)
            {
               istep.z = -1;
               step.z  = (voxel_size_portal.z*double(e.i.z)-e.pos.z)/e.dir.z;
            }
            else
            {
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
         while (repeat_step)
         {
         // the electron will be traced to the next voxel boundary, if the
         // electron reaches the boundary we repeat this step in the next
         // voxel, the bool variable "repeat_step" is "true" in this case
         // cout << "step = " << step.x << " " << step.y << " " << step.z << endl;
            // find the next voxel
            if ( (step.z < step.x) && (step.z < step.y) )
            {
               voxel_step = step.z;
               step.x -= voxel_step;
               step.y -= voxel_step;
               step.z  = voxel_size_portal.z/fabs(e.dir.z);
               inew.z = e.i.z + istep.z;
            }
            else
            {
               if (step.x < step.y)
               {
                  voxel_step = step.x;
                  step.x  = voxel_size_portal.x/fabs(e.dir.x);
                  step.y -= voxel_step;
                  step.z -= voxel_step;
                  inew.x = e.i.x + istep.x;
               }
               else
               {
                  voxel_step = step.y;
                  step.x -= voxel_step;
                  step.y  = voxel_size_portal.y/fabs(e.dir.y);
                  step.z -= voxel_step;
                  inew.y = e.i.y + istep.y;
               }
            }

            // collision stopping power factor of the present voxel
            f_col = dens_scol_portal->matrix[e.i.x][e.i.y][e.i.z]
                  - dens_ccol_portal->matrix[e.i.x][e.i.y][e.i.z]*s_correction;

            // total stopping power factor of the present voxel
            factor =
            (ONE + alpha*dens_srad_portal->matrix[e.i.x][e.i.y][e.i.z])/(ONE + alpha);
            f_tot = f_col*factor;

            // corresponding step size in water
            voxel_sh2o = voxel_step*f_tot;
// cout << "f_col = " << f_col << " " << factor << " " << voxel_sh2o << " " << total_sh2o << " " << voxel_step << endl;
            if (voxel_sh2o >= total_sh2o)
            {
               // the step ends in this voxel
               repeat_step = false;
               path_step  += total_sh2o/f_tot;
               voxel_depo  = total_loss/factor;

               // put the energy in the present voxel
               batch_dose_portal->matrix[e.i.x][e.i.y][e.i.z] +=
                  voxel_depo*e.weight*dose_factor;

               // MS
               sum_chi_cc2 +=
                  total_loss*dens_fchi_portal->matrix[e.i.x][e.i.y][e.i.z];
#ifdef CHECK_ENERGY
               // test energy conservation
               sum_energy.edepo += voxel_depo*e.weight;
               sum_energy.ploss += (total_loss-voxel_depo)*e.weight;
#endif // CHECK_ENERGY
            }
            else
            {
               // the step continues in the next voxel
               voxel_loss  = total_loss*voxel_sh2o/total_sh2o;
               voxel_depo  = voxel_loss/factor;
               total_loss -= voxel_loss;
               total_sh2o -= voxel_sh2o;
               path_step  += voxel_step;

               // put the energy in the present voxel
               batch_dose_portal->matrix[e.i.x][e.i.y][e.i.z] +=
                  voxel_depo*e.weight*dose_factor;

               // MS
               sum_chi_cc2 +=
                  voxel_loss*dens_fchi_portal->matrix[e.i.x][e.i.y][e.i.z];
#ifdef CHECK_ENERGY
               // test energy conservation
               sum_energy.edepo += voxel_depo*e.weight;
               sum_energy.ploss += (voxel_loss-voxel_depo)*e.weight;
#endif // CHECK_ENERGY
					//cout << "inew = " << inew.x << " " << inew.y << " " << inew.z
					//	   << " " << dim_portal.x << " " << dim_portal.y << " " << dim_portal.z << endl;
               if ( (inew.x < 0) || (inew.x >= dim_portal.x) ||
                    (inew.y < 0) || (inew.y >= dim_portal.y) ||
                    (inew.z < 0) || (inew.z >= dim_portal.z) )
               {
                  // the electron leaves the calculation cube, the history
                  // ends here, test energy conservation
#ifdef CHECK_ENERGY
                  sum_energy.eloss +=
                     (total_loss+sub_loss+new_energy)*e.weight;
#endif // CHECK_ENERGY
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
         if (multiple_scattering)
         {  // sample multiple scattering angle between the two sub-steps

            // momentum squared
            p2      = sub_energy*(sub_energy + TWOxEMASS);
            // beta squared, (v/c)^2
            beta2   = p2/(p2 + EMASSxEMASS);
            // number of elastic scatterings
            omega_0 = BC*path/beta2;
            xr      = CHI_CC2*path/p2/beta2/TWO;
            if (save_eloss > DSMALL) xr *= sum_chi_cc2/save_eloss;

            // sample multiple scattering angle (reduced_angle)
#ifdef DEBUG
            if (mscat(omega_0, xr, cos_t, sin_t, rndm))
            {
               xvmc_warning("one_electron","error in mscat",0);
               xvmc_warning("sub_energy",sub_energy,0);
               xvmc_warning("path",path,0);
               xvmc_warning("fchi",dens_fchi_portal->matrix[e.i.x][e.i.y][e.i.z],0);
            }
#else
            mscat(omega_0, xr, cos_t, sin_t, rndm);
#endif // DEBUG
            phi = TWO_PI * rndm.number();
            sin_phi = sin(phi);
            cos_phi = cos(phi);

            // rotate electron direction by the multiple scattering angle
            rotate(cos_t, sin_t, cos_phi, sin_phi, e.dir);

            if (discrete_interaction)
            {
               go_back = true;
               multiple_scattering = false;
               total_sh2o = sub_sh2o;
               total_loss = sub_loss;
               sub_sh2o   = ZERO;
               sub_loss   = ZERO;
            }

            if (!next_step)
            {
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
      if (discrete_interaction)
      {
         // at first, reset number of mean free paths for the next step
         num_mfp = -1.0;

         // determine the type of discrete interaction (Moller or brems. prod.)
         i_e = int((e.energy-e_cut)*e_h2o->inverse_delta_e);
         if (rndm.number() <= e_h2o->p_moller[i_e])
         {  // no Moller interaction if the energy is too small
            // (in this case we jump to the next MS step)
            if ( h2o_moller.interaction(e.energy,   rndm,
                                        new_energy, cos_t_n, sin_t_n,
                                        d.energy,   cos_t_d, sin_t_d) )
            {
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
               one_electron_portal(d,rndm,
#ifdef CHECK_ENERGY
                            sum_energy,
#endif // CHECK_ENERGY
                            batch_dose_portal);

               // direction change of the primary electron
               cos_phi = -cos_phi;
               sin_phi = -sin_phi;

               // rotate primary electron
               rotate(cos_t_n, sin_t_n, cos_phi, sin_phi, e.dir);

               // "e.energy" is now the energy after the Moller interaction
               e.energy = new_energy;
            }
         } // end of Moller interaction
         else
         { // bremsstrahlung production
           // if the energy is larger than TC (corresponds to AP in EGS4)
            if ( h2o_brems.interaction(e.energy,   rndm,
                                       new_energy, cos_t_n, sin_t_n,
                                       p.energy,   cos_t_p, sin_t_p) )
            {
               // "e.energy" is now the energy after bremsstrahlung production
               e.energy = new_energy;
#ifdef CHECK_ENERGY
               sum_energy.ploss += p.energy*e.weight;
#endif // CHECK_ENERGY

               // play Russian Roulette
               if (rndm.number() < P_ROULETTE)
               {
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
#ifdef CHECK_ENERGY
                  sum_energy.ploss -= p.energy*p.weight;
#endif // CHECK_ENERGY
                  // simulate photon
                  kerma_photon_portal(p,rndm,
#ifdef CHECK_ENERGY
                               sum_energy,
#endif // CHECK_ENERGY
                               batch_dose_portal);
               }
            }

         } // end of bremsstrahlung production

      } // end of discrete interaction

   } // end of this step
   return;
}
// --- End of Adding ---------------------------------------------
