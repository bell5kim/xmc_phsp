/*****************************************************************************
 * kerma_photon.cpp:                                                         *
 *    function:                                                              *
 *       kerma_photon: simulate single photon by KERMA approximation         *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 00/02/27        *
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
#include "compton_XS_diff.h"
#include "portal_dose.h"

// ****************************************
// declare functions and global variables
// ****************************************

extern photon_data     *p_h2o;   // water photon transport data for this beam
extern compton_XS_diff  compton; // differential Compton cross section data


// simulate photon history by KERMA approximation
void kerma_photon(particle_parameters &p, ranmar &rndm,
#ifdef CHECK_ENERGY
                  sum_energy_type &sum_energy,
#endif // CHECK_ENERGY
                  array_3d<double> *batch_dose, portal_dose *portal)

{
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

   // update voxel indices
   inew = p.i;

   // start transport
   repeat_step = true;
   path_step   = ZERO;

   // store deposited energy
   tot_energy_depo = ZERO;

   while (repeat_step)
   {
   // the photon will be traced to the next voxel boundary, if the
   // photon reaches the boundary we repeat this step in the next
   // voxel, the bool variable "repeat_step" is "true" in this case

      // determine the interaction probabilities in the present voxel
      dens_corr_phot = dens_phot->matrix[p.i.x][p.i.y][p.i.z];
      dens_corr_pair = dens_pair->matrix[p.i.x][p.i.y][p.i.z];
      dens_corr_comp = dens_comp->matrix[p.i.x][p.i.y][p.i.z];
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

      // calculate deposited energy
      mu_en = mu_en_h2o*dens_corr_comp
            + mu_tot_h2o*p_phot_h2o*p.energy*(dens_corr_phot-dens_corr_comp);
      if (p.energy > TWOxEMASS)
      {
         mu_en += mu_tot_h2o*p_pair_h2o*(p.energy-TWOxEMASS)*
                     (dens_corr_pair-dens_corr_comp);
      }
      energy_depo = p.weight*mu_en*voxel_step;
      tot_energy_depo                        += energy_depo;
      batch_dose->matrix[p.i.x][p.i.y][p.i.z] += energy_depo;

#ifdef CHECK_ENERGY
      // check energy conservation
      sum_energy.pdepo += energy_depo;
#endif // CHECK_ENERGY

      if ( (inew.x < 0) || (inew.x >= dim.x) ||
           (inew.y < 0) || (inew.y >= dim.y) ||
           (inew.z < 0) || (inew.z >= dim.z) )
      {
         // the photon is leaving the calculation cube, the history
         // ends here, count portal image dose and
         if (portal != NULL) portal->add(p, rndm);

         // test energy conservation
#ifdef CHECK_ENERGY
         sum_energy.ploss += (p.energy*p.weight-tot_energy_depo);
#endif // CHECK_ENERGY
         return;
      }

      // update voxel indices
      p.i = inew;

   // if "repeat_step=true" we continue the step in the next voxel
   }

   if (rndm.number() < p_comp)
   {
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
      x.i.x = int(x.pos.x/voxel_size.x);
      x.i.y = int(x.pos.y/voxel_size.y);
      x.i.z = int(x.pos.z/voxel_size.z);

      // simulate the new photon
      kerma_photon(x, rndm,
#ifdef CHECK_ENERGY
                   sum_energy,
#endif // CHECK_ENERGY
                   batch_dose, portal);

#ifdef CHECK_ENERGY
      // forget electron energy
      sum_energy.eloss += energy_e*p.weight - tot_energy_depo;
   }

   else
   {
      // photo absorption and pair production are taken into account by
      // the KERMA approximation, only check energy conservation
      sum_energy.ploss += (p.energy*p.weight-tot_energy_depo);
#endif // CHECK_ENERGY
   }

   return;
}

// --- Added by JKIM 15NOV2010  ------------------------------
// simulate photon history by KERMA approximation

void kerma_photon_portal(particle_parameters &p, ranmar &rndm,
                  array_3d<double> *batch_dose_portal)

{
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

// ****************************************
// begin photon transport
// ****************************************

   // the photon must be within the calculation cube
   if ( (p.i.x < 0) || (p.i.x >= dim_portal.x) ||
        (p.i.y < 0) || (p.i.y >= dim_portal.y) ||
        (p.i.z < 0) || (p.i.z >= dim_portal.z) )
   {
      // end of particle history
      return;
   }

   // kill photon if the energy is below p_cut
   if (p.energy <= p_cut)
   {
#ifdef CHECK_ENERGY
      sum_energy.ploss += p.energy*p.weight;
#endif // CHECK_ENERGY


// Added by JOKim 16Dec2010 --------------
#ifdef NO_PORTAL_PCUT
      // calculate deposited energy
      i_e = 0;
      dens_corr_phot = dens_phot_portal->matrix[p.i.x][p.i.y][p.i.z];
      dens_corr_pair = dens_pair_portal->matrix[p.i.x][p.i.y][p.i.z];
      dens_corr_comp = dens_comp_portal->matrix[p.i.x][p.i.y][p.i.z];

      mu_tot_h2o = p_h2o->mu_tot[i_e];
      mu_en_h2o  = p_h2o->mu_en[i_e];
      p_phot_h2o = p_h2o->mu_phot[i_e]/mu_tot_h2o;
      mu_en = mu_en_h2o*dens_corr_comp
            + mu_tot_h2o*p_phot_h2o*p.energy*(dens_corr_phot-dens_corr_comp);
      energy_depo = p.weight*mu_en*react_ratio;
#ifdef PHOSPHOR
      energy_depo *= react_ratio;
#endif
      tot_energy_depo                                += energy_depo;
      batch_dose_portal->matrix[p.i.x][p.i.y][p.i.z] += energy_depo;
#endif
      return;
   }

#ifdef PHOSPHOR
   // Correction of Reaction Rate in Phosphor
   float react_ratio = 1.0;
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
      react_ratio = pow(10.0, react_ratio) * 2.5;
      // if (react_ratio > 1.85) react_ratio = 1.85;
      // p.weight *= react_ratio;

   }
#endif

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

   // update voxel indices
   inew = p.i;

   // start transport
   repeat_step = true;
   path_step   = ZERO;

   // store deposited energy
   tot_energy_depo = ZERO;

   while (repeat_step)
   {
   // the photon will be traced to the next voxel boundary, if the
   // photon reaches the boundary we repeat this step in the next
   // voxel, the bool variable "repeat_step" is "true" in this case

      // determine the interaction probabilities in the present voxel
      dens_corr_phot = dens_phot_portal->matrix[p.i.x][p.i.y][p.i.z];
      dens_corr_pair = dens_pair_portal->matrix[p.i.x][p.i.y][p.i.z];
      dens_corr_comp = dens_comp_portal->matrix[p.i.x][p.i.y][p.i.z];
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

      // calculate deposited energy
      mu_en = mu_en_h2o*dens_corr_comp
            + mu_tot_h2o*p_phot_h2o*p.energy*(dens_corr_phot-dens_corr_comp);
      if (p.energy > TWOxEMASS)
      {
         mu_en += mu_tot_h2o*p_pair_h2o*(p.energy-TWOxEMASS)*
                     (dens_corr_pair-dens_corr_comp);
      }

      energy_depo = p.weight*mu_en*voxel_step;
#ifdef PHOSPHOR
      energy_depo *= react_ratio;  // Added by JOKim 16Dec2010 --------------
#endif
      tot_energy_depo                        += energy_depo;
      batch_dose_portal->matrix[p.i.x][p.i.y][p.i.z] += energy_depo;

      if ( (inew.x < 0) || (inew.x >= dim_portal.x) ||
           (inew.y < 0) || (inew.y >= dim_portal.y) ||
           (inew.z < 0) || (inew.z >= dim_portal.z) )
      {
         // end of particle history
         return;
      }

      // update voxel indices
      p.i = inew;

   // if "repeat_step=true" we continue the step in the next voxel
   }

   if (rndm.number() < p_comp)
   {
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
      x.i.x = int(x.pos.x/voxel_size_portal.x);
      x.i.y = int(x.pos.y/voxel_size_portal.y);
      x.i.z = int(x.pos.z/voxel_size_portal.z);

      // simulate the new photon
      kerma_photon_portal(x, rndm, batch_dose_portal);
   }

   return;
}

// --- End of Adding  ----------------------------------------

