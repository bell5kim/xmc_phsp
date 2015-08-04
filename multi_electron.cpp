/*****************************************************************************
 * multi_electron.cpp:                                                       *
 *    class member functions for:                                            *
 *       multi_electron: create and simulate electron histories              *
 *                       (history repetition technique)                      *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 00/02/17        *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <math.h>
#include "definitions.h"
#include "global.h"
#include "xvmc_util.h"
#include "electron_data.h"
#include "multi_electron.h"
#include "moller_XS.h"
#include "bhabha_XS.h"
#include "brems_XS.h"

// ****************************************
// declare functions and global variables
// ****************************************

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

// transport delta electron
void one_electron(particle_parameters &, ranmar &,
#ifdef CHECK_ENERGY
                  sum_energy_type &,
#endif // CHECK_ENERGY
                  array_3d<double> *, portal_dose *);

// simulate bremsstrahlung photons by KERMA approximation
void kerma_photon(particle_parameters &, ranmar &,
#ifdef CHECK_ENERGY
                  sum_energy_type &,
#endif // CHECK_ENERGY
                  array_3d<double> *, portal_dose *);

// --- Added by JOKim 15Nov2010 --------------------------------------

// transport delta electron
void one_electron_portal(particle_parameters &, ranmar &,
#ifdef CHECK_ENERGY
                  sum_energy_type &,
#endif // CHECK_ENERGY
                  array_3d<double> *);

// simulate bremsstrahlung photons by KERMA approximation
void kerma_photon_portal(particle_parameters &, ranmar &,
#ifdef CHECK_ENERGY
                  sum_energy_type &,
#endif // CHECK_ENERGY
                  array_3d<double> *);
// --- End of Adding -------------------------------------------------

// create electron history

void multi_electron_create(multi_electron &A, real energy, ranmar &rndm)
{
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

   while (next_step)
   {
   // generate one step

      if (A.n_step >= A.MAX_STEP-1)
      {
         xvmc_warning("multi_electron_create",
                      "too many electron steps (n_step >= MAX_STEP-1)",8);
      }

      if ( (energy <= e_cut) || (A.n_step >= A.MAX_STEP-1) )
      {
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
      else
      {
      // begin: energy > e_cut

         // determine energy index
         i_e=int((energy-e_cut)*e_h2o->inverse_delta_e);

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
               new_energy = energy - eloss;
               if (new_energy > e_cut)
               {
                  new_i_e = int((new_energy-e_cut)*e_h2o->inverse_delta_e);
                  if (rndm.number() <= e_h2o->sigma_tot[new_i_e]) break;
               }
               else
               {
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
         new_energy = energy - eloss;
         if (new_energy <= e_cut)
         {
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
#ifdef DEBUG
      if (mscat(omega_0, xr, cos_z, sin_z, rndm))
      {
         xvmc_warning("multi_electron_create","error in mscat",0);
         xvmc_warning("A.n_step",A.n_step,0);
         xvmc_warning("sub_energy",sub_energy,0);
         xvmc_warning("path",path,0);
         xvmc_warning("A.radiation_loss[A.n_step]",A.radiation_loss[A.n_step],0);
         xvmc_warning("A.radiation_loss[A.n_step-1]",
                       A.radiation_loss[max_of(A.n_step-1,0)],0);
         xvmc_warning("discrete_interaction",discrete_interaction,0);
      }
#else
      mscat(omega_0, xr, cos_z, sin_z, rndm);
#endif // DEBUG

      // MS angle
      A.reduced_angle[A.n_step] = ONE - cos_z;

      // "energy" is now the energy at the end of the MS step
      energy = new_energy;

      // default: no Moller interaction
      A.moller[A.n_step] = false;

      // perform a discrete interaction
      if (discrete_interaction)
      {
         // at first, reset number of mean free paths for the next step
         num_mfp = -1.0;

         // determine the type of discrete interaction (Moller or brems. prod.)
         i_e = int((energy-e_cut)*e_h2o->inverse_delta_e);
         if (rndm.number() <= e_h2o->p_moller[i_e])
         {  // no Moller interaction if the energy is too small
            // (in this case we jump to the next MS step)
            if ( h2o_moller.interaction(energy,     rndm,
                                        new_energy, cos_t_n, sin_t_n,
                                        energy_d,   cos_t_d, sin_t_d) )
            {
               if (A.n_delta >= A.MAX_DELTA)
               {
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
         else
         { // bremsstrahlung production
           // if the energy is larger than TC (corresponds to AP in EGS4)
            if ( h2o_brems.interaction(energy,     rndm,
                                       new_energy, cos_t_n, sin_t_n,
                                       rad_loss,   cos_t_p, sin_t_p) )
            {
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

void multi_electron_trace(multi_electron &A,
                          particle_parameters &e,
                          ranmar &rndm,
#ifdef CHECK_ENERGY
                          sum_energy_type &sum_energy,
#endif // CHECK_ENERGY
                          array_3d<double> *batch_dose,
                          portal_dose      *portal)
{
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

   sub_sh2o = ZERO;
   sub_loss = ZERO;

   for (register int i_step=0; i_step<A.n_step; ++i_step)
   {
   // simulate electron step

#ifdef ELECTRON_TRACK_FULL_STOP
      // deposit the remaining energy in the present voxel
      if (i_step == A.n_step-1)
      {
         total_loss = sub_loss + A.energy_loss[i_step];
         batch_dose->matrix[e.i.x][e.i.y][e.i.z] += total_loss*e.weight;
#ifdef CHECK_ENERGY
         // test energy conservation
         sum_energy.edepo += total_loss*e.weight;
#endif // CHECK_ENERGY
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

      do
      {
      // we go back to this point if there is a discrete interaction
      // at the and of the step or the present step is the last step,
      // i.e. "i_step=A.n_step-1", the bool variable "go_back" must be
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
               if ( (inew.x < 0) || (inew.x >= dim.x) ||
                    (inew.y < 0) || (inew.y >= dim.y) ||
                    (inew.z < 0) || (inew.z >= dim.z) )
               {
                  // the electron leaves the calculation cube, the history
                  // ends here, test energy conservation

                  // count portal image dose
                  if (portal != NULL) portal->add(e, rndm); // Added by JOKim 29Oct2011
#ifdef CHECK_ENERGY
                  sum_energy.eloss += (total_loss+sub_loss)*e.weight;
                  for (int k_step=i_step; k_step<A.n_step; ++k_step)
                  {
                     sum_energy.ploss += A.radiation_loss[k_step]*e.weight;
                     if ( A.moller[k_step] )
                     {
                        sum_energy.eloss +=
                           A.energy_delta[A.i_delta[k_step]]*e.weight;
                     }
                     if (k_step > i_step)
                     {
                        sum_energy.eloss += A.energy_loss[k_step]*e.weight;
                     }
                  }
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

         go_back = false;
         if (multiple_scattering)
         {  // sample multiple scattering angle

            if (save_eloss > DSMALL)
                 temp = A.reduced_angle[i_step]*sum_chi_cc2/save_eloss;
            else temp = ZERO;
            if (temp < DTWO)
            {
               // scattering angles
               cos_t    = ONE - temp;
               sin_t    = sqrt(temp*(ONE+cos_t));
               real phi = TWO_PI * rndm.number(); // azimuthal angle
               sin_phi  = sin(phi);
               cos_phi  = cos(phi);
               sin_z    = (ONE-e.dir.z)*(ONE+e.dir.z);

               // rotate electron
               if (sin_z > DSMALL)
               {
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
               else
               {
                  e.dir.x = sin_t*cos_phi;
                  e.dir.y = sin_t*sin_phi;
                  e.dir.z = cos_t;
               } // if (sin_z > DSMALL)
            }
            else
            {
               temp     = DTWO*rndm.number();
               e.dir.z  = ONE - temp;
               sin_z    = sqrt(temp*(TWO-temp));
               real phi = TWO_PI * rndm.number(); // azimuthal angle
               e.dir.x  = sin_z*cos(phi);
               e.dir.y  = sin_z*sin(phi);
            } // if (temp < DTWO)

            if (A.moller[i_step])
            {
               go_back = true;
               multiple_scattering = false;
               total_sh2o = sub_sh2o;
               total_loss = sub_loss;
               sub_sh2o   = ZERO;
               sub_loss   = ZERO;
            }

            if (i_step == A.n_step-1)
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
      }
      while (go_back);

      // follow brems photons with KERMA approximation
      p.energy = A.radiation_loss[i_step];

      if (p.energy > ZERO)
      {
#ifdef CHECK_ENERGY
         sum_energy.ploss += p.energy*e.weight;
#endif // CHECK_ENERGY
         // play Russian Roulette
         if (rndm.number() < P_ROULETTE)
         {
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

      if ( A.moller[i_step] )
      { // add delta electron contribution

         // get array index
         int i_d = A.i_delta[i_step];

         // calculate moving directions of the primary and delta electrons
         real phi = TWO_PI * rndm.number();
         sin_phi  = sin(phi);
         cos_phi  = cos(phi);
         sin_z    = (ONE-e.dir.z)*(ONE+e.dir.z);
         if ( sin_z > DSMALL )
         {
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
         else
         {
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
         one_electron(d,rndm,
#ifdef CHECK_ENERGY
                      sum_energy,
#endif // CHECK_ENERGY
                      batch_dose, portal);

      } // end of delta electron contribution

   // end of this step
   }

   return;
}

// create positron history
void multi_positron_create(multi_electron &A, real energy, ranmar &rndm)

{
   int    i_e;          // energy index
   int    new_i_e;      // energy index for new_energy

   real   eta;          // random number in [0,1]
   real   zeta;         // random number in [0,1]
   real   num_mfp;      // number of mean free paths
   real   eloss;        // positron energy loss per step
   real   eloss_max;    // maximum positron energy loss per step
   real   new_energy;   // positron energy after one step
   real   sub_energy;   // positron energy after the first sub-step
   real   path;         // positron path for one whole step
   real   p2;           // momentum squared
   real   beta2;        // beta squared, (v/c)^2
   real   omega_0;      // number of elastic scatterings
   real   xr;
   real   cos_z,sin_z;  // multiple scattering angle
   real   sin_t_n;      // sin of discrete scattering angle (prim. positron)
   real   cos_t_n;      // cos of discrete scattering angle (prim. positron)
   real   energy_d;     // delta (Bhabha) electron energy
   real   sin_t_d;      // sin of delta electron scattering angle
   real   cos_t_d;      // cos of delta electron scattering angle
   real   rad_loss;     // bremsstrahlung energy (radiation loss)
   real   sin_t_p;      // sin of bremsstrahlung photon angle
   real   cos_t_p;      // cos of bremsstrahlung photon angle

   // program control
   bool   next_step,discrete_interaction;

// ****************************************
// begin positron transport
// ****************************************

   new_energy =  ZERO;
   next_step  =  true;
   A.n_step     =  0;
   A.n_delta    =  0;
   num_mfp    = -1.0;

   while (next_step)
   {
   // generate one step

      if (A.n_step >= A.MAX_STEP-1)
      {
         xvmc_warning("multi_positron_create",
                      "too many positron steps (n_step >= MAX_STEP-1)",8);
      }

      if ( (energy <= e_cut) || (A.n_step >= A.MAX_STEP-1) )
      {
         next_step                  = false;
         A.energy_loss[A.n_step]    = energy;
         A.radiation_loss[A.n_step] = ZERO;
         A.sin_brems[A.n_step]      = ZERO;
         A.cos_brems[A.n_step]      = ZERO;
         A.dose_cor[A.n_step]       = e_h2o->dose_cor[0];
         eta                        = rndm.number();
         if (eta == ZERO) eta       = rndm.number();
         A.random_number[A.n_step]  = eta;
         path                       = energy/e_h2o->s_res[0];
         A.step_size[A.n_step]      = path;
         A.s_cor[A.n_step]          = e_h2o->s_cor[0];
         A.alpha_r[A.n_step]        = e_h2o->alpha_r[0];
         discrete_interaction       = false;
         sub_energy                 = e_cut;
      }
      else
      {
      // begin: energy > e_cut

         // determine energy index
         i_e=int((energy-e_cut)*e_h2o->inverse_delta_e);

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
               eloss      = num_mfp/e_h2o->sigma_pmx;
               new_energy = energy - eloss;
               if (new_energy > e_cut)
               {
                  new_i_e = int((new_energy-e_cut)*e_h2o->inverse_delta_e);
                  if (rndm.number() <= e_h2o->sigma_pos[new_i_e]) break;
               }
               else
               {
                  // the factor TWO is for numerical reasons only,
                  // it avoids any discrete interaction at the end of the step
                  num_mfp = TWO*energy*e_h2o->sigma_pmx;
                  break;
               }
            }
         }

         // determine energy loss depending on "num_mfp" and "e_h2o->max_loss"
         eloss     = num_mfp/e_h2o->sigma_pmx;
         eloss_max = e_h2o->max_loss[i_e]*energy;
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
         new_energy = energy - eloss;
         if (new_energy <= e_cut)
         {
            eloss                = energy;
            new_energy           = ZERO;
            next_step            = false;
            discrete_interaction = false;
         }

         // the new number of mean free paths
         num_mfp -= eloss*e_h2o->sigma_pmx;

         // use a random number to divide the step into two sub-steps
         eta = rndm.number();
         if (eta == ZERO) eta = rndm.number();
         sub_energy = energy - eloss*eta; // energy after the first sub-step
         if (sub_energy < e_cut) sub_energy = e_cut;
         new_i_e = int((sub_energy-e_cut)*e_h2o->inverse_delta_e);

         // the positron path for the whole step
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

      // sample multiple scattering angle (reduced_angle)
#ifdef DEBUG
      if (mscat(omega_0, xr, cos_z, sin_z, rndm))
      {
         xvmc_warning("multi_positron_create","error in mscat",0);
         xvmc_warning("n_step",A.n_step,0);
         xvmc_warning("sub_energy",sub_energy,0);
         xvmc_warning("path",path,0);
         xvmc_warning("radiation_loss[n_step]",A.radiation_loss[A.n_step],0);
         xvmc_warning("radiation_loss[n_step-1]",
                       A.radiation_loss[max_of(A.n_step-1,0)],0);
         xvmc_warning("discrete_interaction",discrete_interaction,0);
      }
#else
      mscat(omega_0, xr, cos_z, sin_z, rndm);
#endif // DEBUG

      // MS angle
      A.reduced_angle[A.n_step] = ONE - cos_z;

      // "energy" is now the energy at the end of the MS step
      energy = new_energy;

      // default: no Bhabha interaction
      A.bhabha[A.n_step] = false;

      // perform a discrete interaction
      if (discrete_interaction)
      {
         // at first, reset number of mean free paths for the next step
         num_mfp = -1.0;

         // determine the type of discrete interaction (bhabha or brems. prod.)
         i_e = int((energy-e_cut)*e_h2o->inverse_delta_e);
         if (rndm.number() <= e_h2o->p_bhabha[i_e])
         {  // no Bhabha interaction if the energy is too small
            // (in this case we jump to the next MS step)
            if ( h2o_bhabha.interaction(energy,     rndm,
                                        new_energy, cos_t_n, sin_t_n,
                                        energy_d,   cos_t_d, sin_t_d) )
            {
               if (A.n_delta >= A.MAX_DELTA)
               {
                  xvmc_error("multi_positron_create",
                    "too many delta electrons (n_delta >= MAX_DELTA)",8);
               }

               // count Bhabha interaction (delta electron production)
               A.bhabha[A.n_step]        = true;
               A.i_delta[A.n_step]       = A.n_delta;
               A.energy_delta[A.n_delta] = energy_d;

               // store the polar scattering angles
               // (positron and delta electron)
               A.sin_theta[A.n_step]     = sin_t_n;
               A.cos_theta[A.n_step]     = cos_t_n;
               A.sin_theta_d[A.n_delta]  = sin_t_d;
               A.cos_theta_d[A.n_delta]  = cos_t_d;

               // "energy" is now the energy after the Bhabha interaction
               energy = new_energy;

               ++A.n_delta;
            }
         } // end of Bhabha interaction
         else
         { // bremsstrahlung production
           // if the energy is larger than TC (corresponds to AP in EGS4)
            if ( h2o_brems.interaction(energy,     rndm,
                                       new_energy, cos_t_n, sin_t_n,
                                       rad_loss,   cos_t_p, sin_t_p) )
            {
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

// trace pre-calculated positron history through the calculation cube
void multi_positron_trace(multi_electron &A,
                          particle_parameters &e,
                          ranmar &rndm,
#ifdef CHECK_ENERGY
                          sum_energy_type &sum_energy,
#endif // CHECK_ENERGY
                          array_3d<double> *batch_dose,
                          portal_dose      *portal)
{
   real   eta;          // random number in [0,1]
   real   total_sh2o;   // water step size of the first positron sub-step
   real   total_loss;   // energy loss during the first positron sub-step
   real   sub_sh2o;     // water step size of the second positron sub-step
   real   sub_loss;     // energy loss during the second positron sub-step
   real   alpha;        // alpha_r for this step
   real   dose_factor;  // correction for dose, ionization etc.
   real   f_col;        // function f_c(rho) to calculate the collision
                        // stopping power (see eq. (19) in VMC paper I)
   real   f_tot;        // total stopping power factor of the present voxel
   real   factor;       // ratio f_tot/f_col
   real   s_correction; // low energy correction for collision stopping power
   real   path_step;    // positron path length from the beginning of the step
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
// begin positron transport
// ****************************************

   // the positron must be within the calculation cube
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

   sub_sh2o = ZERO;
   sub_loss = ZERO;

   for (register int i_step=0; i_step<A.n_step; ++i_step)
   {
   // simulate positron step

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

      do
      {
      // we go back to this point if there is a discrete interaction
      // at the and of the step or the present step is the last step,
      // i.e. "i_step=n_step-1", the bool variable "go_back" must be
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
         // the positron will be traced to the next voxel boundary, if the
         // positron reaches the boundary we repeat this step in the next
         // voxel, the bool variable "repeat_step" is "true" in this case

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
               if ( (inew.x < 0) || (inew.x >= dim.x) ||
                    (inew.y < 0) || (inew.y >= dim.y) ||
                    (inew.z < 0) || (inew.z >= dim.z) )
               {
                  // the positron leaves the calculation cube, the history
                  // ends here, test energy conservation

      				// count portal image dose
      				if (portal != NULL) portal->add(e, rndm); // Added by JOKim 29Oct2011
#ifdef CHECK_ENERGY
                  sum_energy.eloss += (total_loss+sub_loss)*e.weight;
                  for (int k_step=i_step; k_step<A.n_step; ++k_step)
                  {
                     sum_energy.ploss += A.radiation_loss[k_step]*e.weight;
                     if ( A.bhabha[k_step] )
                     {
                        sum_energy.eloss +=
                           A.energy_delta[A.i_delta[k_step]]*e.weight;
                     }
                     if (k_step > i_step)
                     {
                        sum_energy.eloss += A.energy_loss[k_step]*e.weight;
                     }
                  }
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

         go_back = false;
         if (multiple_scattering)
         {  // sample multiple scattering angle

            if (save_eloss > DSMALL)
                 temp = A.reduced_angle[i_step]*sum_chi_cc2/save_eloss;
            else temp = ZERO;
            if (temp < DTWO)
            {
               // scattering angles
               cos_t    = ONE - temp;
               sin_t    = sqrt(temp*(ONE+cos_t));
               real phi = TWO_PI * rndm.number(); // azimuthal angle
               sin_phi  = sin(phi);
               cos_phi  = cos(phi);
               sin_z    = (ONE-e.dir.z)*(ONE+e.dir.z);

               // rotate positron
               if (sin_z > DSMALL)
               {
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
               else
               {
                  e.dir.x = sin_t*cos_phi;
                  e.dir.y = sin_t*sin_phi;
                  e.dir.z = cos_t;
               } // if (sin_z > DSMALL)
            }
            else
            {
               temp     = DTWO*rndm.number();
               e.dir.z  = ONE - temp;
               sin_z    = sqrt(temp*(TWO-temp));
               real phi = TWO_PI * rndm.number(); // azimuthal angle
               e.dir.x  = sin_z*cos(phi);
               e.dir.y  = sin_z*sin(phi);
            } // if (temp < DTWO)

            if (A.bhabha[i_step])
            {
               go_back = true;
               multiple_scattering = false;
               total_sh2o = sub_sh2o;
               total_loss = sub_loss;
               sub_sh2o   = ZERO;
               sub_loss   = ZERO;
            }

            if (i_step == A.n_step-1)
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
      }
      while (go_back);

      // follow brems photons with KERMA approximation
      p.energy = A.radiation_loss[i_step];

      if (p.energy > ZERO)
      {
#ifdef CHECK_ENERGY
         sum_energy.ploss += p.energy*e.weight;
#endif // CHECK_ENERGY
         // play Russian Roulette
         if (rndm.number() < P_ROULETTE)
         {
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

      if ( A.bhabha[i_step] )
      { // add delta electron contribution

         // get array index
         int i_d = A.i_delta[i_step];

         // calculate moving directions of the positron and delta electron
         real phi = TWO_PI * rndm.number();
         sin_phi  = sin(phi);
         cos_phi  = cos(phi);
         sin_z    = (ONE-e.dir.z)*(ONE+e.dir.z);
         if ( sin_z > DSMALL )
         {
            sin_z = sqrt(sin_z);

            // rotate delta electron
            cos_t    = A.cos_theta_d[i_d];
            sin_t    = A.sin_theta_d[i_d];
            temp     = sin_t/sin_z;
            temp_phi = e.dir.z*cos_phi;
            d.dir.x  = e.dir.x*cos_t + temp*(temp_phi*e.dir.x-e.dir.y*sin_phi);
            d.dir.y  = e.dir.y*cos_t + temp*(temp_phi*e.dir.y+e.dir.x*sin_phi);
            d.dir.z  = e.dir.z*cos_t - sin_z*sin_t*cos_phi;

            // rotate positron
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
         else
         {
            // rotate delta electron
            cos_t    = A.cos_theta_d[i_d];
            sin_t    = A.sin_theta_d[i_d];
            d.dir.x  = sin_t*cos_phi;
            d.dir.y  = sin_t*sin_phi;
            d.dir.z  = cos_t;

            // rotate positron
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
         one_electron(d,rndm,
#ifdef CHECK_ENERGY
                      sum_energy,
#endif // CHECK_ENERGY
                      batch_dose, portal);

      } // end of delta electron contribution

   // end of this step
   }

   return;
}

// --- Added by JOKim 15Nov2010 -------------------------------------

// trace pre-calculated electron history through the calculation cube
void multi_electron_trace_portal(multi_electron &A,
                          particle_parameters &e,
                          ranmar &rndm,
#ifdef CHECK_ENERGY
                          sum_energy_type &sum_energy,
#endif // CHECK_ENERGY
                          array_3d<double> *batch_dose_portal)
{
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

   sub_sh2o = ZERO;
   sub_loss = ZERO;

   for (register int i_step=0; i_step<A.n_step; ++i_step)
   {
   // simulate electron step

#ifdef ELECTRON_TRACK_FULL_STOP
      // deposit the remaining energy in the present voxel
      if (i_step == A.n_step-1)
      {
         total_loss = sub_loss + A.energy_loss[i_step];
         batch_dose_portal->matrix[e.i.x][e.i.y][e.i.z] += total_loss*e.weight;
#ifdef CHECK_ENERGY
         // test energy conservation
         sum_energy.edepo += total_loss*e.weight;
#endif // CHECK_ENERGY
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

      do
      {
      // we go back to this point if there is a discrete interaction
      // at the and of the step or the present step is the last step,
      // i.e. "i_step=A.n_step-1", the bool variable "go_back" must be
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
               if ( (inew.x < 0) || (inew.x >= dim_portal.x) ||
                    (inew.y < 0) || (inew.y >= dim_portal.y) ||
                    (inew.z < 0) || (inew.z >= dim_portal.z) )
               {
                  // the electron leaves the calculation cube, the history
                  // ends here, test energy conservation
#ifdef CHECK_ENERGY
                  sum_energy.eloss += (total_loss+sub_loss)*e.weight;
                  for (int k_step=i_step; k_step<A.n_step; ++k_step)
                  {
                     sum_energy.ploss += A.radiation_loss[k_step]*e.weight;
                     if ( A.moller[k_step] )
                     {
                        sum_energy.eloss +=
                           A.energy_delta[A.i_delta[k_step]]*e.weight;
                     }
                     if (k_step > i_step)
                     {
                        sum_energy.eloss += A.energy_loss[k_step]*e.weight;
                     }
                  }
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

         go_back = false;
         if (multiple_scattering)
         {  // sample multiple scattering angle

            if (save_eloss > DSMALL)
                 temp = A.reduced_angle[i_step]*sum_chi_cc2/save_eloss;
            else temp = ZERO;
            if (temp < DTWO)
            {
               // scattering angles
               cos_t    = ONE - temp;
               sin_t    = sqrt(temp*(ONE+cos_t));
               real phi = TWO_PI * rndm.number(); // azimuthal angle
               sin_phi  = sin(phi);
               cos_phi  = cos(phi);
               sin_z    = (ONE-e.dir.z)*(ONE+e.dir.z);

               // rotate electron
               if (sin_z > DSMALL)
               {
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
               else
               {
                  e.dir.x = sin_t*cos_phi;
                  e.dir.y = sin_t*sin_phi;
                  e.dir.z = cos_t;
               } // if (sin_z > DSMALL)
            }
            else
            {
               temp     = DTWO*rndm.number();
               e.dir.z  = ONE - temp;
               sin_z    = sqrt(temp*(TWO-temp));
               real phi = TWO_PI * rndm.number(); // azimuthal angle
               e.dir.x  = sin_z*cos(phi);
               e.dir.y  = sin_z*sin(phi);
            } // if (temp < DTWO)

            if (A.moller[i_step])
            {
               go_back = true;
               multiple_scattering = false;
               total_sh2o = sub_sh2o;
               total_loss = sub_loss;
               sub_sh2o   = ZERO;
               sub_loss   = ZERO;
            }

            if (i_step == A.n_step-1)
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
      }
      while (go_back);

      // follow brems photons with KERMA approximation
      p.energy = A.radiation_loss[i_step];

      if (p.energy > ZERO)
      {
#ifdef CHECK_ENERGY
         sum_energy.ploss += p.energy*e.weight;
#endif // CHECK_ENERGY
         // play Russian Roulette
         if (rndm.number() < P_ROULETTE)
         {
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

      if ( A.moller[i_step] )
      { // add delta electron contribution

         // get array index
         int i_d = A.i_delta[i_step];

         // calculate moving directions of the primary and delta electrons
         real phi = TWO_PI * rndm.number();
         sin_phi  = sin(phi);
         cos_phi  = cos(phi);
         sin_z    = (ONE-e.dir.z)*(ONE+e.dir.z);
         if ( sin_z > DSMALL )
         {
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
         else
         {
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
         one_electron_portal(d,rndm,
#ifdef CHECK_ENERGY
                      sum_energy,
#endif // CHECK_ENERGY
                      batch_dose_portal);

      } // end of delta electron contribution

   // end of this step
   }

   return;
}

// trace pre-calculated positron history through the calculation cube
void multi_positron_trace_portal(multi_electron &A,
                          particle_parameters &e,
                          ranmar &rndm,
#ifdef CHECK_ENERGY
                          sum_energy_type &sum_energy,
#endif // CHECK_ENERGY
                          array_3d<double> *batch_dose_portal)

{
   real   eta;          // random number in [0,1]
   real   total_sh2o;   // water step size of the first positron sub-step
   real   total_loss;   // energy loss during the first positron sub-step
   real   sub_sh2o;     // water step size of the second positron sub-step
   real   sub_loss;     // energy loss during the second positron sub-step
   real   alpha;        // alpha_r for this step
   real   dose_factor;  // correction for dose, ionization etc.
   real   f_col;        // function f_c(rho) to calculate the collision
                        // stopping power (see eq. (19) in VMC paper I)
   real   f_tot;        // total stopping power factor of the present voxel
   real   factor;       // ratio f_tot/f_col
   real   s_correction; // low energy correction for collision stopping power
   real   path_step;    // positron path length from the beginning of the step
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
// begin positron transport
// ****************************************

   // the positron must be within the calculation cube
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

   sub_sh2o = ZERO;
   sub_loss = ZERO;

   for (register int i_step=0; i_step<A.n_step; ++i_step)
   {
   // simulate positron step

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

      do
      {
      // we go back to this point if there is a discrete interaction
      // at the and of the step or the present step is the last step,
      // i.e. "i_step=n_step-1", the bool variable "go_back" must be
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
         // the positron will be traced to the next voxel boundary, if the
         // positron reaches the boundary we repeat this step in the next
         // voxel, the bool variable "repeat_step" is "true" in this case

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
               if ( (inew.x < 0) || (inew.x >= dim_portal.x) ||
                    (inew.y < 0) || (inew.y >= dim_portal.y) ||
                    (inew.z < 0) || (inew.z >= dim_portal.z) )
               {
                  // the positron leaves the calculation cube, the history
                  // ends here, test energy conservation
#ifdef CHECK_ENERGY
                  sum_energy.eloss += (total_loss+sub_loss)*e.weight;
                  for (int k_step=i_step; k_step<A.n_step; ++k_step)
                  {
                     sum_energy.ploss += A.radiation_loss[k_step]*e.weight;
                     if ( A.bhabha[k_step] )
                     {
                        sum_energy.eloss +=
                           A.energy_delta[A.i_delta[k_step]]*e.weight;
                     }
                     if (k_step > i_step)
                     {
                        sum_energy.eloss += A.energy_loss[k_step]*e.weight;
                     }
                  }
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

         go_back = false;
         if (multiple_scattering)
         {  // sample multiple scattering angle

            if (save_eloss > DSMALL)
                 temp = A.reduced_angle[i_step]*sum_chi_cc2/save_eloss;
            else temp = ZERO;
            if (temp < DTWO)
            {
               // scattering angles
               cos_t    = ONE - temp;
               sin_t    = sqrt(temp*(ONE+cos_t));
               real phi = TWO_PI * rndm.number(); // azimuthal angle
               sin_phi  = sin(phi);
               cos_phi  = cos(phi);
               sin_z    = (ONE-e.dir.z)*(ONE+e.dir.z);

               // rotate positron
               if (sin_z > DSMALL)
               {
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
               else
               {
                  e.dir.x = sin_t*cos_phi;
                  e.dir.y = sin_t*sin_phi;
                  e.dir.z = cos_t;
               } // if (sin_z > DSMALL)
            }
            else
            {
               temp     = DTWO*rndm.number();
               e.dir.z  = ONE - temp;
               sin_z    = sqrt(temp*(TWO-temp));
               real phi = TWO_PI * rndm.number(); // azimuthal angle
               e.dir.x  = sin_z*cos(phi);
               e.dir.y  = sin_z*sin(phi);
            } // if (temp < DTWO)

            if (A.bhabha[i_step])
            {
               go_back = true;
               multiple_scattering = false;
               total_sh2o = sub_sh2o;
               total_loss = sub_loss;
               sub_sh2o   = ZERO;
               sub_loss   = ZERO;
            }

            if (i_step == A.n_step-1)
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
      }
      while (go_back);

      // follow brems photons with KERMA approximation
      p.energy = A.radiation_loss[i_step];

      if (p.energy > ZERO)
      {
#ifdef CHECK_ENERGY
         sum_energy.ploss += p.energy*e.weight;
#endif // CHECK_ENERGY
         // play Russian Roulette
         if (rndm.number() < P_ROULETTE)
         {
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

      if ( A.bhabha[i_step] )
      { // add delta electron contribution

         // get array index
         int i_d = A.i_delta[i_step];

         // calculate moving directions of the positron and delta electron
         real phi = TWO_PI * rndm.number();
         sin_phi  = sin(phi);
         cos_phi  = cos(phi);
         sin_z    = (ONE-e.dir.z)*(ONE+e.dir.z);
         if ( sin_z > DSMALL )
         {
            sin_z = sqrt(sin_z);

            // rotate delta electron
            cos_t    = A.cos_theta_d[i_d];
            sin_t    = A.sin_theta_d[i_d];
            temp     = sin_t/sin_z;
            temp_phi = e.dir.z*cos_phi;
            d.dir.x  = e.dir.x*cos_t + temp*(temp_phi*e.dir.x-e.dir.y*sin_phi);
            d.dir.y  = e.dir.y*cos_t + temp*(temp_phi*e.dir.y+e.dir.x*sin_phi);
            d.dir.z  = e.dir.z*cos_t - sin_z*sin_t*cos_phi;

            // rotate positron
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
         else
         {
            // rotate delta electron
            cos_t    = A.cos_theta_d[i_d];
            sin_t    = A.sin_theta_d[i_d];
            d.dir.x  = sin_t*cos_phi;
            d.dir.y  = sin_t*sin_phi;
            d.dir.z  = cos_t;

            // rotate positron
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
         one_electron_portal(d,rndm,
#ifdef CHECK_ENERGY
                      sum_energy,
#endif // CHECK_ENERGY
                      batch_dose_portal);

      } // end of delta electron contribution

   // end of this step
   }

   return;
}

// --- End of Adding ------------------------------------------------
