/*****************************************************************************
 * electron_data.cpp:                                                        *
 *    class member functions for:                                            *
 *       electron_data:  electron transport data (function of energy)        *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 00/02/22        *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include "definitions.h"
#include "global.h"
#include "electron_data.h"


/*****************************************************************************
 * member functions for electron_data:                                       *
 *****************************************************************************/

// use (copy) transport data from input array to calculate the array
// of transport parameters for MC simulation (better resolution)

extern moller_XS h2o_moller; // Moller cross section for water
extern bhabha_XS h2o_bhabha; // Bhabha cross section for water
extern brems_XS  h2o_brems;  // bremsstrahlung cross section for water
void electron_data_init(electron_data &A, unsigned int number, electron_data_inp *inp,
                             real energy_max, real energy_mean,
                             result_type result, real film_factor)
{
   register unsigned int i; // loop control variable
   real ee,delta,p,q;
   int  lower,upper,j;
   real p2;        // momentum squared
   real xx;        // inverse screening angle squared
   real beta2;     // beta squared, (v/c)^2
   real ts;        // scattering power
   real s_col_tmp; // temporary collision stopping power value
   real s_rad_tmp; // temporary radiation stopping power value
   real av_e_delta;  // average energy of the delta electron (E > E_cut)
   real av_e_brems;  // average energy of the bremsstrahlung photon (k > TC)
   real sigma_mol;   // total Moller cross section
   real sigma_bha;   // total Bhabha cross section
   real sigma_brems; // total cross section for bremsstrahlung production
   real tot_step_min; // total minimum step size
   real tot_step_max; // total maximum step size
   real step_max;     // maximum step size for the given energy
   real *s_scat;      // temporary scattering power array

   if (energy_max < e_cut)
       xvmc_error("electron_data::electron_data","energy_max < e_cut",8);

   // allocate memory for MC simulation electron data
   A.num       = number;
   A.delta_e   = (1.01*energy_max - e_cut)/double(A.num);
   A.inverse_delta_e = ONE/A.delta_e;
   bool error = false;
   if ( (A.energy    = new real[A.num]) == NULL ) error = true;
   if ( (A.s_res     = new real[A.num]) == NULL ) error = true;
   if ( (A.s_cor     = new real[A.num]) == NULL ) error = true;
   if ( (A.alpha_r   = new real[A.num]) == NULL ) error = true;
   if ( (A.s_tot     = new real[A.num]) == NULL ) error = true;
   if ( (A.dose_cor  = new real[A.num]) == NULL ) error = true;
   if ( (A.max_loss  = new real[A.num]) == NULL ) error = true;
   if ( (A.sigma_tot = new real[A.num]) == NULL ) error = true;
   A.sigma_max = ZERO;
   if ( (A.p_moller  = new real[A.num]) == NULL ) error = true;
   if ( (A.sigma_pos = new real[A.num]) == NULL ) error = true;
   A.sigma_pmx = ZERO;
   if ( (A.p_bhabha  = new real[A.num]) == NULL ) error = true;

   // allocate memory for temporary scattering power array
   if ( (s_scat    = new real[A.num]) == NULL ) error = true;

   if (error)
   {
      xvmc_error("electron_data::electron_data",
                 "cannot allocate memory for electron transport data",8);
   }

   // calculate interpolated parameters
   ee        = e_cut - A.delta_e/TWO;
   for (i=0; i<A.num; ++i)
   {
      ee = ee + A.delta_e;
      get_index(inp->num,inp->energy,ee,lower,upper);
      delta = inp->energy[upper] - inp->energy[lower];
      if (delta > ZERO)
      {
         p = (ee - inp->energy[lower])/delta;
         q = ONE-p;
      }
      else
      {
         p = ONE;
         q = ZERO;
      }

      // energy array
      A.energy[i]  = ee;

      // total stopping power
      A.s_tot[i] = q*inp->s_tot[lower] + p*inp->s_tot[upper];

      // unrestricted collision stopping power
      s_col_tmp = q*inp->s_col[lower] + p*inp->s_col[upper];

      // total Moller cross section
      sigma_mol = h2o_moller.total(ee);

      // calculate restricted collision stopping power: eq. (6) and (7) of
      // Med. Phys. 23 (1996) 445 (VMC paper I)
      // and mean delta electron energy
      av_e_delta = ZERO;
      if (sigma_mol > ZERO)
      {
         // first moment of the Moller cross section
         real first_moment_mol = h2o_moller.first(ee);
         av_e_delta = first_moment_mol/sigma_mol;
         A.s_res[i] = s_col_tmp - first_moment_mol;
      }
      else
      {
         A.s_res[i] = s_col_tmp;
      }

      // total Bhabha cross section
      sigma_bha = h2o_bhabha.total(ee);

      // low energy correction for collision stopping power
      if (ee > ONE) A.s_cor[i] = exp(-(ee-ONE));
      else          A.s_cor[i] = ONE;
      // if (ee > ONE_HALF) s_cor[i] = exp(ONE-0.7*ee);
      // else               s_cor[i] = 1.9155408290139;

      // radiation stopping power
      s_rad_tmp  = q*inp->s_rad[lower] + p*inp->s_rad[upper];

      // ratio s_rad/s_col (alpha_r)
      // A.alpha_r[i] = s_rad_tmp/s_col_tmp;
      A.alpha_r[i] = s_rad_tmp/A.s_res[i];

      // total bremsstrahlung cross section: eq. (6) and (7) of
      // Med. Phys. 23 (1996) 445 (VMC paper I)
      // mean bremstrahlung photon energy
      av_e_brems  = h2o_brems.first(ee)/h2o_brems.total(ee);
      sigma_brems = s_rad_tmp/av_e_brems;

      // total cross section for discrete electron interaction
      A.sigma_tot[i] = sigma_mol + sigma_brems;

      // probability for discrete Moller interaction
      A.p_moller[i] = sigma_mol/A.sigma_tot[i];

      // total cross section for discrete positron interaction
      A.sigma_pos[i] = sigma_bha + sigma_brems;

      // probability for discrete Bhabha interaction
      A.p_bhabha[i] = sigma_bha/A.sigma_pos[i];

      // divide total cross section by the restricted collision stopping power
      // to change from position into the energy space (we sample energy loss
      // instead of the distance to the next discrete interaction)
      A.sigma_tot[i] = A.sigma_tot[i]/A.s_res[i];
      A.sigma_pos[i] = A.sigma_pos[i]/A.s_res[i];

      // calculate the maximum cross section for the fictitious interaction
      // method (SLAC-265, p. 16), later we will calculate the ratios
      if (A.sigma_tot[i] > A.sigma_max) A.sigma_max = A.sigma_tot[i];
      if (A.sigma_pos[i] > A.sigma_pmx) A.sigma_pmx = A.sigma_pos[i];

      // additional dose correction factor depending on the result type:
      // dose to water, film dose, ionization (dose to air) or energy dose,
      // "dose_cor" is used to store the factor
      switch (result)
      {
      case DOSE_TO_H2O:    // dose to water
         A.dose_cor[i] = ONE;
         break;
      case FILM:           // film dose
         A.dose_cor[i] = q*inp->s_photo[lower] + p*inp->s_photo[upper];
         A.dose_cor[i] = A.dose_cor[i]*film_factor;
         break;
      case IONIZATION:     // ionization chamber in water (factor: s_air/s_h2o)
         A.dose_cor[i] = q*inp->s_air[lower] + p*inp->s_air[upper];
         break;
      default:             // energy dose (no correction)
         A.dose_cor[i] = ONE;
         break;
      }

      // re-calculate electron scattering power for linear interpolation
      p2    = ee*(ee+TWOxEMASS);    // momentum squared
      xx    = 4.0*p2/AMU2;          // inverse screening angle squared
      beta2 = p2/(p2+EMASSxEMASS);  // beta squared, (v/c)^2
      // scattering power
      ts    = 4.0*BC*((ONE+ONE/xx)*log(ONE+xx)-ONE)/xx/beta2;
      s_scat[i]=(q*inp->s_scat[lower] + p*inp->s_scat[upper])*ts;
   }

   // Initialize total minimum and maximum step sizes depending on energy,
   // voxel sizes, estepe, etc.,
   //    total minimum step size:
   tot_step_min = voxel_size.x;
   if (voxel_size.y < tot_step_min) tot_step_min = voxel_size.y;
   if (voxel_size.z < tot_step_min) tot_step_min = voxel_size.z;
   tot_step_min = 0.4 * tot_step_min;
   if (tot_step_min > energy_max/ 20.0) tot_step_min = energy_max/20.0;
   if (tot_step_min < energy_max/100.0) tot_step_min = energy_max/100.0;

   //    total maximum step size:
   j = int((energy_max-e_cut) * A.inverse_delta_e);
   tot_step_max = e_step*energy_mean/A.s_res[j];
   if (tot_step_max > TEFF0_600) tot_step_max = TEFF0_600;

   if (tot_step_min > tot_step_max) tot_step_min = tot_step_max;

   xvmc_message("Total minimum step size:",tot_step_min,"cm",1);
   xvmc_message("Total maximum step size:",tot_step_max,"cm",0);

   // calculate maximum step size for the given energy
   ee = e_cut - A.delta_e/TWO;
   for (i=0; i<A.num; ++i)
   {
      ee = ee + A.delta_e;
      step_max=ONE/s_scat[i];
      if (step_max > tot_step_max) step_max=tot_step_max;
      if (step_max < tot_step_min) step_max=tot_step_min;
      // maximum energy loss in units of the electron energy
      A.max_loss[i] = step_max*A.s_res[i]/ee;
      if (A.max_loss[i] > ONE) A.max_loss[i] = ONE;

      // normalize sigma_tot/sigma_pos for the fictitious interaction method
      A.sigma_tot[i] = A.sigma_tot[i]/A.sigma_max;
      A.sigma_pos[i] = A.sigma_pos[i]/A.sigma_pmx;
   }

   // we do not need the scattering power any longer
   delete [] s_scat;
}

//############################################################################
//########################NEWNEWNEWNEWNEWNEW##################################
//############################################################################

electron_data::electron_data(int number,
			     electron_data_inp& inp,
                             real energy_max,
			     real energy_mean,
			     real e_step,
			     real eCutOff,
			     const moller_XS& h2o_moller,  // Moller cross section for water
			     const bhabha_XS& h2o_bhabha,  // Bhabha cross section for water
			     const brems_XS&  h2o_brems,
			     const real_3& voxel_size,
                             result_type result) :
  num(number)
{
   register unsigned int i; // loop control variable
   real ee,delta,p,q;
   int  lower,upper,j;
   real p2;        // momentum squared
   real xx;        // inverse screening angle squared
   real beta2;     // beta squared, (v/c)^2
   real ts;        // scattering power
   real s_col_tmp; // temporary collision stopping power value
   real s_rad_tmp; // temporary radiation stopping power value
   real av_e_delta;  // average energy of the delta electron (E > E_cut)
   real av_e_brems;  // average energy of the bremsstrahlung photon (k > TC)
   real sigma_mol;   // total Moller cross section
   real sigma_bha;   // total Bhabha cross section
   real sigma_brems; // total cross section for bremsstrahlung production
   real tot_step_min; // total minimum step size
   real tot_step_max; // total maximum step size
   real step_max;     // maximum step size for the given energy
   real *s_scat;      // temporary scattering power array

   if (energy_max < eCutOff)
       xvmc_error("electron_data::electron_data","energy_max < eCutOff",8);

   // allocate memory for MC simulation electron data
   delta_e   = (1.01*energy_max - eCutOff)/double(num);
   inverse_delta_e = ONE/delta_e;
   bool error = false;
   if ( (energy    = new real[num]) == NULL ) error = true;
   if ( (s_res     = new real[num]) == NULL ) error = true;
   if ( (s_cor     = new real[num]) == NULL ) error = true;
   if ( (alpha_r   = new real[num]) == NULL ) error = true;
   if ( (s_tot     = new real[num]) == NULL ) error = true;
   if ( (dose_cor  = new real[num]) == NULL ) error = true;
   if ( (max_loss  = new real[num]) == NULL ) error = true;
   if ( (sigma_tot = new real[num]) == NULL ) error = true;
   sigma_max = ZERO;
   if ( (p_moller  = new real[num]) == NULL ) error = true;
   if ( (sigma_pos = new real[num]) == NULL ) error = true;
   sigma_pmx = ZERO;
   if ( (p_bhabha  = new real[num]) == NULL ) error = true;

   // allocate memory for temporary scattering power array
   if ( (s_scat    = new real[num]) == NULL ) error = true;

   if (error)
   {
      xvmc_error("electron_data::electron_data",
                 "cannot allocate memory for electron transport data",8);
   }

   // calculate interpolated parameters
   ee        = eCutOff - delta_e/TWO;
   for (i=0; i<num; ++i)
   {
      ee = ee + delta_e;
      get_index(inp.num,inp.energy,ee,lower,upper);
      delta = inp.energy[upper] - inp.energy[lower];
      if (delta > ZERO)
      {
         p = (ee - inp.energy[lower])/delta;
         q = ONE-p;
      }
      else
      {
         p = ONE;
         q = ZERO;
      }

      // energy array
      energy[i]  = ee;

      // total stopping power
      s_tot[i] = q*inp.s_tot[lower] + p*inp.s_tot[upper];

      // unrestricted collision stopping power
      s_col_tmp = q*inp.s_col[lower] + p*inp.s_col[upper];

      // total Moller cross section
      sigma_mol = h2o_moller.total(ee);

      // calculate restricted collision stopping power: eq. (6) and (7) of
      // Med. Phys. 23 (1996) 445 (VMC paper I)
      // and mean delta electron energy
      av_e_delta = ZERO;
      if (sigma_mol > ZERO)
      {
         // first moment of the Moller cross section
         real first_moment_mol = h2o_moller.first(ee);
         av_e_delta = first_moment_mol/sigma_mol;
         s_res[i] = s_col_tmp - first_moment_mol;
      }
      else
      {
         s_res[i] = s_col_tmp;
      }

      // total Bhabha cross section
      sigma_bha = h2o_bhabha.total(ee);

      // low energy correction for collision stopping power
      if (ee > ONE) s_cor[i] = exp(-(ee-ONE));
      else          s_cor[i] = ONE;
      // if (ee > ONE_HALF) s_cor[i] = exp(ONE-0.7*ee);
      // else               s_cor[i] = 1.9155408290139;

      // radiation stopping power
      s_rad_tmp  = q*inp.s_rad[lower] + p*inp.s_rad[upper];

      // ratio s_rad/s_col (alpha_r)
      // alpha_r[i] = s_rad_tmp/s_col_tmp;
      alpha_r[i] = s_rad_tmp/s_res[i];

      // total bremsstrahlung cross section: eq. (6) and (7) of
      // Med. Phys. 23 (1996) 445 (VMC paper I)
      // mean bremstrahlung photon energy
      av_e_brems  = h2o_brems.first(ee)/h2o_brems.total(ee);
      sigma_brems = s_rad_tmp/av_e_brems;

      // total cross section for discrete electron interaction
      sigma_tot[i] = sigma_mol + sigma_brems;

      // probability for discrete Moller interaction
      p_moller[i] = sigma_mol/sigma_tot[i];

      // total cross section for discrete positron interaction
      sigma_pos[i] = sigma_bha + sigma_brems;

      // probability for discrete Bhabha interaction
      p_bhabha[i] = sigma_bha/sigma_pos[i];

      // divide total cross section by the restricted collision stopping power
      // to change from position into the energy space (we sample energy loss
      // instead of the distance to the next discrete interaction)
      sigma_tot[i] = sigma_tot[i]/s_res[i];
      sigma_pos[i] = sigma_pos[i]/s_res[i];

      // calculate the maximum cross section for the fictitious interaction
      // method (SLAC-265, p. 16), later we will calculate the ratios
      if (sigma_tot[i] > sigma_max) sigma_max = sigma_tot[i];
      if (sigma_pos[i] > sigma_pmx) sigma_pmx = sigma_pos[i];

      // additional dose correction factor depending on the result type:
      // dose to water, film dose, ionization (dose to air) or energy dose,
      // "dose_cor" is used to store the factor
      switch (result)
      {
      case DOSE_TO_H2O:    // dose to water
         dose_cor[i] = ONE;
         break;
//      case FILM:           // film dose
//         dose_cor[i] = q*inp.s_photo[lower] + p*inp.s_photo[upper];
//         dose_cor[i] = dose_cor[i]*film_factor;
//         break;
      case IONIZATION:     // ionization chamber in water (factor: s_air/s_h2o)
         dose_cor[i] = q*inp.s_air[lower] + p*inp.s_air[upper];
         break;
      default:             // energy dose (no correction)
         dose_cor[i] = ONE;
         break;
      }

      // re-calculate electron scattering power for linear interpolation
      p2    = ee*(ee+TWOxEMASS);    // momentum squared
      xx    = 4.0*p2/AMU2;          // inverse screening angle squared
      beta2 = p2/(p2+EMASSxEMASS);  // beta squared, (v/c)^2
      // scattering power
      ts    = 4.0*BC*((ONE+ONE/xx)*log(ONE+xx)-ONE)/xx/beta2;
      s_scat[i]=(q*inp.s_scat[lower] + p*inp.s_scat[upper])*ts;
   }

   // Initialize total minimum and maximum step sizes depending on energy,
   // voxel sizes, estepe, etc.,
   //    total minimum step size:
   tot_step_min = voxel_size.x;
   if (voxel_size.y < tot_step_min) tot_step_min = voxel_size.y;
   if (voxel_size.z < tot_step_min) tot_step_min = voxel_size.z;
   tot_step_min = 0.4 * tot_step_min;
   if (tot_step_min > energy_max/ 20.0) tot_step_min = energy_max/20.0;
   if (tot_step_min < energy_max/100.0) tot_step_min = energy_max/100.0;

   //    total maximum step size:
   j = int((energy_max-eCutOff) * inverse_delta_e);
   tot_step_max = e_step*energy_mean/s_res[j];
   if (tot_step_max > TEFF0_600) tot_step_max = TEFF0_600;

   if (tot_step_min > tot_step_max) tot_step_min = tot_step_max;

   xvmc_message("Total minimum step size:",tot_step_min,"cm",1);
   xvmc_message("Total maximum step size:",tot_step_max,"cm",0);

   // calculate maximum step size for the given energy
   ee = eCutOff - delta_e/TWO;
   for (i=0; i<num; ++i)
   {
      ee = ee + delta_e;
      step_max=ONE/s_scat[i];
      if (step_max > tot_step_max) step_max=tot_step_max;
      if (step_max < tot_step_min) step_max=tot_step_min;
      // maximum energy loss in units of the electron energy
      max_loss[i] = step_max*s_res[i]/ee;
      if (max_loss[i] > ONE) max_loss[i] = ONE;

      // normalize sigma_tot/sigma_pos for the fictitious interaction method
      sigma_tot[i] = sigma_tot[i]/sigma_max;
      sigma_pos[i] = sigma_pos[i]/sigma_pmx;
   }

   // we do not need the scattering power any longer
   delete [] s_scat;
}



// delete data
electron_data::~electron_data()
{
   num = 0;
   delta_e = ZERO;
   inverse_delta_e = ZERO;
   delete [] energy;
   delete [] s_res;
   delete [] s_cor;
   delete [] alpha_r;
   delete [] s_tot;
   delete [] dose_cor;
   delete [] max_loss;
   delete [] sigma_tot;
   sigma_max = ZERO;
   delete [] p_moller;
   delete [] sigma_pos;
   sigma_pmx = ZERO;
   delete [] p_bhabha;
}
