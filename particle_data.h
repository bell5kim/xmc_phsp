#ifndef _PARTICLE_DATA_H_
#define _PARTICLE_DATA_H_

/*****************************************************************************
 * particle_data.h:                                                          *
 *    class declarations:                                                    *
 *       electron_data_inp: input electron transport data                    *
 *       electron_data:     electron transport data (function of energy)     *
 *       photon_data_inp:   input photon transport data                      *
 *       photon_data:       photon transport data (function of energy)       *
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

#include <stdlib.h>
#include "definitions.h"
#include "array_3d.h"

// input electron transport data as a function of energy
class electron_data_inp
{
   friend class electron_data;
   friend void electron_data_init(electron_data &, unsigned int, electron_data_inp *,
                    real, real, result_type, real = 1.0);

   public:
      unsigned int num;   // number of array elements
      real   *energy;     // energy array
      real   *s_col;      // collision stopping power array
      real   *s_rad;      // radiation stopping power array
      real   *s_tot;      // total stopping power array
      real   *s_air;      // air stopping power array
      real   *s_photo;    // photo factor
      real   *s_scat;     // scattering power array

   public:
      electron_data_inp(char *); // input electron transport data from file
      ~electron_data_inp();      // delete data
};

// electron transport data as a function of energy
class electron_data
{
   public:
      unsigned int num;   // number of array elements
      real   delta_e;     // energy resolution
      real   inverse_delta_e; // 1/delta_e
      real   *energy;     // energy array
      real   *s_res;      // restricted collision stopping power array
      real   *s_cor;      // low energy correction for stopping power
      real   *alpha_r;    // array of ratios s_rad/s_col
      real   *s_tot;      // total stopping power array
      real   *dose_cor;   // additional dose correction factor depending on
                          // the plan result type: dose to water, film dose,
                          // ionization (dose to air) or energy dose
      real   *max_loss;   // maximum electron energy loss per step (in units of
                          // the electron energy
      real   *sigma_tot;  // total cross section array for electron interaction
                          // devided by the restricted collision stopping power
                          // and normalized by the maximum (sigma_max)
      real    sigma_max;  // maximum of sigma_tot[i] before normalization
      real   *p_moller;   // probability array for discrete Moller interaction
      real   *sigma_pos;  // total cross section array for positron interaction
                          // devided by the restricted collision stopping power
                          // and normalized by the maximum (sigma_pmx)
      real    sigma_pmx;  // maximum of sigma_pos[i] before normalization
      real   *p_bhabha;   // probability array for discrete Bhabha interaction

   public:
      // calculate transport parameters for MC simulation
      // by copying the input array (parameters from file)
      electron_data() {
	energy = NULL;
        s_res = NULL;
	s_cor = NULL;
        alpha_r = NULL;
        s_tot = NULL;
        dose_cor = NULL;
	max_loss = NULL;
	sigma_tot = NULL;
	p_moller = NULL;
	sigma_pos = NULL;
	p_bhabha = NULL;
      }
      ~electron_data();   // delete data
};

void electron_data_init(electron_data &A, unsigned int, electron_data_inp *,
                    real, real, result_type, real);

// input photon transport data as a function of energy
class photon_data_inp
{
   friend class photon_data;
   friend void photon_data_init(photon_data &, unsigned int, photon_data_inp *, real);

   public:
      unsigned int num;   // number of array elements
      real   *energy;     // energy array
      real   *mu_comp;    // compton cross section array
      real   *mu_pair;    // pair cross section array
      real   *mu_phot;    // photo cross section array
      real   *mu_tot;     // total cross section array
      real   *mu_en;      // array of energy absorption coefficients

   public:
      photon_data_inp(char *); // input photon transport data from file
      ~photon_data_inp();      // delete data
};

// photon transport data as a function of energy
class photon_data
{
   public:
      unsigned int num;   // number of array elements
      real   delta_e;     // energy resolution
      real   inverse_delta_e; // 1/delta_e
      real   *energy;     // energy array
      real   *mu_comp;    // compton cross section array
      real   *mu_pair;    // pair cross section array
      real   *mu_phot;    // photo cross section array
      real   *mu_tot;     // total cross section array
      real   *mu_en;      // array of energy absorption coefficients

   public:
      // calculate transport parameters for MC simulation
      // by copying the input array (parameters from file)
      photon_data() {
	energy = NULL;
	mu_comp = NULL;
	mu_pair = NULL;
	mu_phot = NULL;
	mu_tot = NULL;
	mu_en = NULL;
      }
      ~photon_data();     // delete data
};

void photon_data_init(photon_data &A, unsigned int, photon_data_inp *, real);

#endif /* _PARTICLE_DATA_H_ */
