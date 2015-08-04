/*****************************************************************************
 * multi_photon.cpp:                                                         *
 *    function:                                                              *
 *       multi_photon: simulate photon history with variance reduction       *
 *                     techniques: interaction forcing, electron history     *
 *                     repetition and Russian Roulette                       *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 00/03/02        *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <math.h>
#include "definitions.h"
#include "global.h"
#include "xvmc_util.h"
#include "photon_data.h"
#include "multi_electron.h"
#include "compton_XS_diff.h"
#include "pair_XS_diff.h"
#include "portal_dose.h"

// ****************************************
// declare functions and global variables
// ****************************************

extern photon_data     *p_h2o;    // water photon transport data for this beam
extern compton_XS_diff  compton;  // differential Compton cross section data
extern pair_XS_diff     dif_pair; // differential pair cross section data

// simulate Compton photons with energy less than "k1_cut"
// by KERMA approximation
void kerma_photon(particle_parameters &, ranmar &,
#ifdef CHECK_ENERGY
                  sum_energy_type &,
#endif // CHECK_ENERGY
                  array_3d<double> *, portal_dose *);

// -- Added by JOKim 15Nov2010  --------------------------------------
void kerma_photon_portal(particle_parameters &, ranmar &,
#ifdef CHECK_ENERGY
                  sum_energy_type &,
#endif // CHECK_ENERGY
                  array_3d<double> *);


// --- End of Adding -------------------------------------------------

// simulate photon history
void multi_photon(particle_parameters &p, int n_repeat, ranmar &rndm,
#ifdef USE_SOBOL
                  real *zeta,
#endif // USE_SOBOL
#ifdef CHECK_ENERGY
                  sum_energy_type &sum_energy,
#endif // CHECK_ENERGY
                  array_3d<double> *batch_dose, portal_dose *portal)
{
#ifdef USE_SOBOL
   // store Sobol random numbers to call "multi_photon" recursively
   real   new_zeta[SOBOL_DIM];
#else
   real   zeta[3];        // array to store random numbers
   for (register int i=0; i<3; ++i) zeta[i] = rndm.number();
#endif // USE_SOBOL
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
   multi_electron electron_comp; // Compton electron
   multi_electron electron_pair; // pair electron
   multi_electron positron_pair; // pair positron
   multi_electron electron_phot; // photo electron

   bool next_interaction,repeat_step;

// ****************************************
// begin photon transport
// ****************************************

   // the photon must be within the calculation cube
   if ( (p.i.x < 0) || (p.i.x >= dim.x) ||
        (p.i.y < 0) || (p.i.y >= dim.y) ||
        (p.i.z < 0) || (p.i.z >= dim.z) )
   {
      // count portal image dose
      if (portal != NULL) portal->add(p, rndm);

#ifdef CHECK_ENERGY
      // test energy conservation
      sum_energy.ploss += p.energy*p.weight;
#endif // CHECK_ENERGY
      return;
   }

   // kill photon if the energy is below p_cut
   if (p.energy <= p_cut)
   {
#ifdef CHECK_ENERGY
      sum_energy.ploss += p.energy*p.weight;
#endif // CHECK_ENERGY
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
   if (p.dir.x > ZERO)
   {
      istep.x = 1;
      step.x  = (voxel_size.x*double(p.i.x+1)-p.pos.x)/p.dir.x;
   }
   else
   {
      if (p.dir.x < ZERO)
      {
         istep.x = -1;
         step.x  = (voxel_size.x*double(p.i.x)-p.pos.x)/p.dir.x;
      }
      else
      {
         istep.x = 0;
         step.x  = HUGE_STEP;
      }
   }

   // check Y-direction
   if (p.dir.y > ZERO)
   {
      istep.y = 1;
      step.y  = (voxel_size.y*double(p.i.y+1)-p.pos.y)/p.dir.y;
   }
   else
   {
      if (p.dir.y < ZERO)
      {
         istep.y = -1;
         step.y  = (voxel_size.y*double(p.i.y)-p.pos.y)/p.dir.y;
      }
      else
      {
         istep.y = 0;
         step.y  = HUGE_STEP;
      }
   }

   // check Z-direction
   if (p.dir.z > ZERO)
   {
      istep.z = 1;
      step.z  = (voxel_size.z*double(p.i.z+1)-p.pos.z)/p.dir.z;
   }
   else
   {
      if (p.dir.z < ZERO)
      {
         istep.z = -1;
         step.z  = (voxel_size.z*double(p.i.z)-p.pos.z)/p.dir.z;
      }
      else
      {
         istep.z = 0;
         step.z  = HUGE_STEP;
      }
   }

   next_interaction = true;
   save_mfp         = ZERO;
   while (next_interaction)
   {
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

      while (repeat_step)
      {
      // the photon will be traced to the next voxel boundary, if the
      // photon reaches the boundary we repeat this step in the next
      // voxel, the bool variable "repeat_step" is "true" in this case

         // determine the interaction probabilities in the present voxel
         prob_phot = prob_phot_h2o*dens_phot->matrix[inew.x][inew.y][inew.z];
         prob_pair = prob_pair_h2o*dens_pair->matrix[inew.x][inew.y][inew.z];
         prob_comp = prob_comp_h2o*dens_comp->matrix[inew.x][inew.y][inew.z];
         prob_tot  = prob_phot + prob_pair + prob_comp;

         // total cross section in this voxel
         mu_tot = mu_tot_h2o*prob_tot;

         // path to the next interaction
         path = num_mfp/mu_tot;

         // find the next voxel
         if ( (step.z < step.x) && (step.z < step.y) )
         {
            voxel_step = step.z;
            if (voxel_step > path)
            {
               voxel_step  = path;
               step.z     -= voxel_step;
               repeat_step = false;
            }
            else
            {
               step.z  = voxel_size.z/fabs(p.dir.z);
               inew.z += istep.z;
            }
            step.x -= voxel_step;
            step.y -= voxel_step;
         }
         else
         {
            if (step.x < step.y)
            {
               voxel_step = step.x;
               if (voxel_step > path)
               {
                  voxel_step  = path;
                  step.x     -= voxel_step;
                  repeat_step = false;
               }
               else
               {
                  step.x  = voxel_size.x/fabs(p.dir.x);
                  inew.x += istep.x;
               }
               step.y -= voxel_step;
               step.z -= voxel_step;
            }
            else
            {
               voxel_step = step.y;
               if (voxel_step > path)
               {
                  voxel_step  = path;
                  step.y     -= voxel_step;
                  repeat_step = false;
               }
               else
               {
                  step.y  = voxel_size.y/fabs(p.dir.y);
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
         if ( (inew.x < 0) || (inew.x >= dim.x) ||
              (inew.y < 0) || (inew.y >= dim.y) ||
              (inew.z < 0) || (inew.z >= dim.z) )
         {
            // the photon is leaving the calculation cube, the history
            // ends here, change weight and count portal image dose
            p.weight *= double(n_repeat-i_repeat)/double(n_repeat);
            if (portal != NULL) portal->add(p, rndm);

            // test energy conservation
#ifdef CHECK_ENERGY
            sum_energy.ploss += p.energy*p.weight;
#endif // CHECK_ENERGY
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
      p.i.x = int(p.pos.x/voxel_size.x);
      p.i.y = int(p.pos.y/voxel_size.y);
      p.i.z = int(p.pos.z/voxel_size.z);

      // check voxel index
      if ( (p.i.x < 0) || (p.i.x >= dim.x) ||
           (p.i.y < 0) || (p.i.y >= dim.y) ||
           (p.i.z < 0) || (p.i.z >= dim.z) )
      {
         // the photon is leaving the calculation cube, the history
         // ends here, change weight and count portal image dose
         p.weight *= double(n_repeat-i_repeat)/double(n_repeat);
         if (portal != NULL) portal->add(p, rndm);

         // test energy conservation
#ifdef CHECK_ENERGY
         sum_energy.ploss += p.energy*p.weight;
#endif // CHECK_ENERGY
         return;
      }

      // find the X-, Y-, Z-distance to the next voxel boundary
      // check X-direction
      if (p.i.x != inew.x)
      {
         if (istep.x == 1)
         {
            step.x  = (voxel_size.x*double(p.i.x+1)-p.pos.x)/p.dir.x;
         }
         else
         {
            if (istep.x == -1)
            {
               step.x  = (voxel_size.x*double(p.i.x)-p.pos.x)/p.dir.x;
            }
            else
            {
               step.x  = HUGE_STEP;
            }
         }
      }

      // check Y-direction
      if (p.i.y != inew.y)
      {
         if (istep.y == 1)
         {
            step.y  = (voxel_size.y*double(p.i.y+1)-p.pos.y)/p.dir.y;
         }
         else
         {
            if (istep.y == -1)
            {
               step.y  = (voxel_size.y*double(p.i.y)-p.pos.y)/p.dir.y;
            }
            else
            {
               step.y  = HUGE_STEP;
            }
         }
      }

      // check Z-direction
      if (p.i.z != inew.z)
      {
         if (istep.z == 1)
         {
            step.z  = (voxel_size.z*double(p.i.z+1)-p.pos.z)/p.dir.z;
         }
         else
         {
            if (istep.z == -1)
            {
               step.z  = (voxel_size.z*double(p.i.z)-p.pos.z)/p.dir.z;
            }
            else
            {
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

      if (zeta[1] <= prob_comp)
      {
         // Compton scattering

         if (first_comp)
         {
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
         multi_electron_trace(electron_comp, e_comp, rndm,
#ifdef CHECK_ENERGY
                              sum_energy,
#endif // CHECK_ENERGY
                              batch_dose, portal);

#ifdef CHECK_ENERGY
         sum_energy.ploss += x_comp_energy*p.weight*delta_eta;
#endif // CHECK_ENERGY
         // with the photon play Russian Roulette
         if (rndm.number() < delta_eta)
         {
            // photon direction
            sin_phi = -sin_phi;
            cos_phi = -cos_phi;
            x_comp.dir = p.dir;
            rotate(cos_t_cx, sin_t_cx, cos_phi, sin_phi, x_comp.dir);

            x_comp.energy = x_comp_energy;
            x_comp.weight = p.weight;
            x_comp.pos    = p.pos;
            x_comp.i      = p.i;
#ifdef CHECK_ENERGY
            sum_energy.ploss -= x_comp.energy*x_comp.weight;
#endif // CHECK_ENERGY
            // simulate Compton photon
            if (x_comp_energy <= k1_cut)
            {
               // KERMA approximation for energies less than "k1_cut"
               kerma_photon(x_comp, rndm,
#ifdef CHECK_ENERGY
                            sum_energy,
#endif // CHECK_ENERGY
                            batch_dose, portal);
            }
            else
            {
#ifdef USE_SOBOL
               // calculate new array of random numbers "new_zeta"
               for (register int i=0; i<SOBOL_DIM; ++i)
               {
                  if (i < (SOBOL_DIM-3)) new_zeta[i] = zeta[i+3];
                  else                   new_zeta[i] = rndm.number();
               }
#endif
               multi_photon(x_comp, n_repeat, rndm,
#ifdef USE_SOBOL
                            new_zeta,
#endif // USE_SOBOL
#ifdef CHECK_ENERGY
                            sum_energy,
#endif // CHECK_ENERGY
                            batch_dose, portal);
            }
         } // end of Russian Roulette

      } // end of Compton interaction
      else
      {
         // pair production or photo effect
         if (zeta[1] < prob_pair)
         {
            // pair production
            if (first_pair)
            {
               dif_pair.interaction(p.energy, zeta[2], rndm,
                                    e_pair_energy, cos_t_pe, sin_t_pe,
                                    p_pair_energy, cos_t_pp, sin_t_pp);

               // create the pair electron history in water
               multi_electron_create(electron_pair, e_pair_energy, rndm);

               // create the pair positron history in water
               multi_positron_create(positron_pair, p_pair_energy, rndm);

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
            multi_electron_trace(electron_pair, e_pair, rndm,
#ifdef CHECK_ENERGY
                                 sum_energy,
#endif // CHECK_ENERGY
                                 batch_dose, portal);

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
            multi_positron_trace(positron_pair, p_pair, rndm,
#ifdef CHECK_ENERGY
                                 sum_energy,
#endif // CHECK_ENERGY
                                 batch_dose, portal);
#ifdef CHECK_ENERGY
            // forget annihilation photons
            sum_energy.ploss += TWOxEMASS*p.weight*delta_eta;
#endif // CHECK_ENERGY

         } // end of pair production
         else
         {
            // photo electric absorption
            if (first_phot)
            {
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
            multi_electron_trace(electron_phot, e_phot, rndm,
#ifdef CHECK_ENERGY
                                 sum_energy,
#endif // CHECK_ENERGY
                                 batch_dose, portal);

         } // end of photo electric absorption

      }

   } // while(next_interaction)

   return;
}

// --- Added by JKim 15Nov2010 ------------------------------
// simulate Compton photons with energy less than "k1_cut"
// by KERMA approximation

// simulate photon history
void multi_photon_portal(particle_parameters &p, int n_repeat, ranmar &rndm,
#ifdef USE_SOBOL
                  real *zeta,
#endif // USE_SOBOL
#ifdef CHECK_ENERGY
                  sum_energy_type &sum_energy,
#endif // CHECK_ENERGY
                  array_3d<double> *batch_dose_portal)
{
#ifdef USE_SOBOL
   // store Sobol random numbers to call "multi_photon" recursively
   real   new_zeta[SOBOL_DIM];
#else
   real   zeta[3];        // array to store random numbers
   for (register int i=0; i<3; ++i) zeta[i] = rndm.number();
#endif // USE_SOBOL
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
   multi_electron electron_comp; // Compton electron
   multi_electron electron_pair; // pair electron
   multi_electron positron_pair; // pair positron
   multi_electron electron_phot; // photo electron

   bool next_interaction,repeat_step;

// ****************************************
// begin photon transport
// ****************************************

   // the photon must be within the calculation cube
   if ( (p.i.x < 0) || (p.i.x >= dim_portal.x) ||
        (p.i.y < 0) || (p.i.y >= dim_portal.y) ||
        (p.i.z < 0) || (p.i.z >= dim_portal.z) )
   {

#ifdef CHECK_ENERGY
      // test energy conservation
      sum_energy.ploss += p.energy*p.weight;
#endif // CHECK_ENERGY
      return;
   }

   // kill photon if the energy is below p_cut
   if (p.energy <= p_cut)
   {
#ifdef CHECK_ENERGY
      sum_energy.ploss += p.energy*p.weight;
#endif // CHECK_ENERGY
      return;
   }
/*
   // Correction of Reaction Rate in Phosphor
   float react_ratio = 1.0;
   float pWeight = p.weight;
   if (p.energy > 0.001) {
      float logE = log10(p.energy);
      //float react_ratio = 1.728389e-2*pow(logE,4.0)+ 3.785169e-1*pow(logE,3.0)
      //                  + 7.304943e-1*logE*logE    - 6.880345e-1*logE
      //                  - 2.929900e-2;
      react_ratio = -5.038339e-2*pow(logE,5.0)
                        -  2.281870e-1*pow(logE,4.0)
                        +  1.383199e-1*pow(logE,3.0)
                        +  9.547950e-1*pow(logE,2.0)
                        -  4.642250e-1*logE
                        -  6.189300e-2;
      react_ratio = pow(10.0, react_ratio);
      // if (react_ratio > 1.85) react_ratio = 1.85;
      p.weight *= react_ratio;
   }
*/
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
   if (p.dir.x > ZERO)
   {
      istep.x = 1;
      step.x  = (voxel_size_portal.x*double(p.i.x+1)-p.pos.x)/p.dir.x;
   }
   else
   {
      if (p.dir.x < ZERO)
      {
         istep.x = -1;
         step.x  = (voxel_size_portal.x*double(p.i.x)-p.pos.x)/p.dir.x;
      }
      else
      {
         istep.x = 0;
         step.x  = HUGE_STEP;
      }
   }

   // check Y-direction
   if (p.dir.y > ZERO)
   {
      istep.y = 1;
      step.y  = (voxel_size_portal.y*double(p.i.y+1)-p.pos.y)/p.dir.y;
   }
   else
   {
      if (p.dir.y < ZERO)
      {
         istep.y = -1;
         step.y  = (voxel_size_portal.y*double(p.i.y)-p.pos.y)/p.dir.y;
      }
      else
      {
         istep.y = 0;
         step.y  = HUGE_STEP;
      }
   }

   // check Z-direction
   if (p.dir.z > ZERO)
   {
      istep.z = 1;
      step.z  = (voxel_size_portal.z*double(p.i.z+1)-p.pos.z)/p.dir.z;
   }
   else
   {
      if (p.dir.z < ZERO)
      {
         istep.z = -1;
         step.z  = (voxel_size_portal.z*double(p.i.z)-p.pos.z)/p.dir.z;
      }
      else
      {
         istep.z = 0;
         step.z  = HUGE_STEP;
      }
   }

   next_interaction = true;
   save_mfp         = ZERO;
   while (next_interaction)
   {
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

      while (repeat_step)
      {
      // the photon will be traced to the next voxel boundary, if the
      // photon reaches the boundary we repeat this step in the next
      // voxel, the bool variable "repeat_step" is "true" in this case

         // determine the interaction probabilities in the present voxel
         prob_phot = prob_phot_h2o*dens_phot_portal->matrix[inew.x][inew.y][inew.z];
         prob_pair = prob_pair_h2o*dens_pair_portal->matrix[inew.x][inew.y][inew.z];
         prob_comp = prob_comp_h2o*dens_comp_portal->matrix[inew.x][inew.y][inew.z];
         prob_tot  = prob_phot + prob_pair + prob_comp;

         // total cross section in this voxel
         mu_tot = mu_tot_h2o*prob_tot;

         // path to the next interaction
         path = num_mfp/mu_tot;

         // find the next voxel
         if ( (step.z < step.x) && (step.z < step.y) )
         {
            voxel_step = step.z;
            if (voxel_step > path)
            {
               voxel_step  = path;
               step.z     -= voxel_step;
               repeat_step = false;
            }
            else
            {
               step.z  = voxel_size_portal.z/fabs(p.dir.z);
               inew.z += istep.z;
            }
            step.x -= voxel_step;
            step.y -= voxel_step;
         }
         else
         {
            if (step.x < step.y)
            {
               voxel_step = step.x;
               if (voxel_step > path)
               {
                  voxel_step  = path;
                  step.x     -= voxel_step;
                  repeat_step = false;
               }
               else
               {
                  step.x  = voxel_size_portal.x/fabs(p.dir.x);
                  inew.x += istep.x;
               }
               step.y -= voxel_step;
               step.z -= voxel_step;
            }
            else
            {
               voxel_step = step.y;
               if (voxel_step > path)
               {
                  voxel_step  = path;
                  step.y     -= voxel_step;
                  repeat_step = false;
               }
               else
               {
                  step.y  = voxel_size_portal.y/fabs(p.dir.y);
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
         if ( (inew.x < 0) || (inew.x >= dim_portal.x) ||
              (inew.y < 0) || (inew.y >= dim_portal.y) ||
              (inew.z < 0) || (inew.z >= dim_portal.z) )
         {
            // the photon is leaving the calculation cube, the history
            // ends here, change weight and count portal image dose
            p.weight *= double(n_repeat-i_repeat)/double(n_repeat);
            // test energy conservation
#ifdef CHECK_ENERGY
            sum_energy.ploss += p.energy*p.weight;
#endif // CHECK_ENERGY
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
      p.i.x = int(p.pos.x/voxel_size_portal.x);
      p.i.y = int(p.pos.y/voxel_size_portal.y);
      p.i.z = int(p.pos.z/voxel_size_portal.z);

      // check voxel index
      if ( (p.i.x < 0) || (p.i.x >= dim_portal.x) ||
           (p.i.y < 0) || (p.i.y >= dim_portal.y) ||
           (p.i.z < 0) || (p.i.z >= dim_portal.z) )
      {
         // the photon is leaving the calculation cube, the history
         // ends here, change weight and count portal image dose
         p.weight *= double(n_repeat-i_repeat)/double(n_repeat);
         // test energy conservation
#ifdef CHECK_ENERGY
         sum_energy.ploss += p.energy*p.weight;
#endif // CHECK_ENERGY
         return;
      }

      // find the X-, Y-, Z-distance to the next voxel boundary
      // check X-direction
      if (p.i.x != inew.x)
      {
         if (istep.x == 1)
         {
            step.x  = (voxel_size_portal.x*double(p.i.x+1)-p.pos.x)/p.dir.x;
         }
         else
         {
            if (istep.x == -1)
            {
               step.x  = (voxel_size_portal.x*double(p.i.x)-p.pos.x)/p.dir.x;
            }
            else
            {
               step.x  = HUGE_STEP;
            }
         }
      }

      // check Y-direction
      if (p.i.y != inew.y)
      {
         if (istep.y == 1)
         {
            step.y  = (voxel_size_portal.y*double(p.i.y+1)-p.pos.y)/p.dir.y;
         }
         else
         {
            if (istep.y == -1)
            {
               step.y  = (voxel_size_portal.y*double(p.i.y)-p.pos.y)/p.dir.y;
            }
            else
            {
               step.y  = HUGE_STEP;
            }
         }
      }

      // check Z-direction
      if (p.i.z != inew.z)
      {
         if (istep.z == 1)
         {
            step.z  = (voxel_size_portal.z*double(p.i.z+1)-p.pos.z)/p.dir.z;
         }
         else
         {
            if (istep.z == -1)
            {
               step.z  = (voxel_size_portal.z*double(p.i.z)-p.pos.z)/p.dir.z;
            }
            else
            {
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

      if (zeta[1] <= prob_comp)
      {
         // Compton scattering

         if (first_comp)
         {
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
         // e_comp.weight *= react_ratio; // Added by JOKim 16Dec2011
         e_comp.pos    = p.pos;
         e_comp.i      = p.i;
         multi_electron_trace_portal(electron_comp, e_comp, rndm,
#ifdef CHECK_ENERGY
                              sum_energy,
#endif // CHECK_ENERGY
                              batch_dose_portal);

#ifdef CHECK_ENERGY
         sum_energy.ploss += x_comp_energy*p.weight*delta_eta;
#endif // CHECK_ENERGY
         // with the photon play Russian Roulette
         if (rndm.number() < delta_eta)
         {
            // photon direction
            sin_phi = -sin_phi;
            cos_phi = -cos_phi;
            x_comp.dir = p.dir;
            rotate(cos_t_cx, sin_t_cx, cos_phi, sin_phi, x_comp.dir);

            x_comp.energy = x_comp_energy;
            x_comp.weight = p.weight;
            x_comp.pos    = p.pos;
            x_comp.i      = p.i;
#ifdef CHECK_ENERGY
            sum_energy.ploss -= x_comp.energy*x_comp.weight;
#endif // CHECK_ENERGY
            // simulate Compton photon
            if (x_comp_energy <= k1_cut)
            {
               // KERMA approximation for energies less than "k1_cut"
               kerma_photon_portal(x_comp, rndm,
#ifdef CHECK_ENERGY
                            sum_energy,
#endif // CHECK_ENERGY
                            batch_dose_portal);
            }
            else
            {
#ifdef USE_SOBOL
               // calculate new array of random numbers "new_zeta"
               for (register int i=0; i<SOBOL_DIM; ++i)
               {
                  if (i < (SOBOL_DIM-3)) new_zeta[i] = zeta[i+3];
                  else                   new_zeta[i] = rndm.number();
               }
#endif
               multi_photon_portal(x_comp, n_repeat, rndm,
#ifdef USE_SOBOL
                            new_zeta,
#endif // USE_SOBOL
#ifdef CHECK_ENERGY
                            sum_energy,
#endif // CHECK_ENERGY
                            batch_dose_portal);
            }
         } // end of Russian Roulette

      } // end of Compton interaction
      else
      {
         // pair production or photo effect
         if (zeta[1] < prob_pair)
         {
            // pair production
            if (first_pair)
            {
               dif_pair.interaction(p.energy, zeta[2], rndm,
                                    e_pair_energy, cos_t_pe, sin_t_pe,
                                    p_pair_energy, cos_t_pp, sin_t_pp);

               // create the pair electron history in water
               multi_electron_create(electron_pair, e_pair_energy, rndm);

               // create the pair positron history in water
               multi_positron_create(positron_pair, p_pair_energy, rndm);

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
            // e_pair.weight *= react_ratio; // Added by JOKim 16Dec2011
            e_pair.pos    = p.pos;
            e_pair.i      = p.i;
            multi_electron_trace_portal(electron_pair, e_pair, rndm,
#ifdef CHECK_ENERGY
                                 sum_energy,
#endif // CHECK_ENERGY
                                 batch_dose_portal);

            // positron direction
            sin_phi = -sin_phi;
            cos_phi = -cos_phi;
            p_pair.dir = p.dir;
            rotate(cos_t_pp, sin_t_pp, cos_phi, sin_phi, p_pair.dir);

            // trace the pair positron
            p_pair.energy = p_pair_energy;
            p_pair.weight = p.weight*delta_eta;
            // p_pair.weight *= react_ratio; // Added by JOKim 16Dec2011
            p_pair.pos    = p.pos;
            p_pair.i      = p.i;
            multi_positron_trace_portal(positron_pair, p_pair, rndm,
#ifdef CHECK_ENERGY
                                 sum_energy,
#endif // CHECK_ENERGY
                                 batch_dose_portal);
#ifdef CHECK_ENERGY
            // forget annihilation photons
            sum_energy.ploss += TWOxEMASS*p.weight*delta_eta;
#endif // CHECK_ENERGY

         } // end of pair production
         else
         {
            // photo electric absorption
            if (first_phot)
            {
               // the photon energy is transfered to the electron
               e_phot_energy = p.energy;

               // create the photo electron history in water
               multi_electron_create(electron_phot, e_phot_energy, rndm);

               first_phot = false;
            }

            // trace the photo electron
            e_phot.energy = e_phot_energy;
            e_phot.weight = p.weight*delta_eta;
            // e_phot.weight *= react_ratio; // Added by JOKim 16Dec2011
            // the photo electron has the direction of th incoming photon
            e_phot.dir    = p.dir;
            e_phot.pos    = p.pos;
            e_phot.i      = p.i;
            multi_electron_trace_portal(electron_phot, e_phot, rndm,
#ifdef CHECK_ENERGY
                                 sum_energy,
#endif // CHECK_ENERGY
                                 batch_dose_portal);

         } // end of photo electric absorption

      }

   } // while(next_interaction)

   return;
}


// --- End of Adding ----------------------------------------










