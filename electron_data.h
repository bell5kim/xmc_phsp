#ifndef _ELECTRON_DATA_H_
#define _ELECTRON_DATA_H_

/*****************************************************************************
 * electron_data.h:                                                          *
 *    class declarations and inline member functions for:                    *
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
#include "electron_data_inp.h"
#include "moller_XS.h"
#include "bhabha_XS.h"
#include "brems_XS.h"

// electron transport data as a function of energy
class electron_data
{
 public:

  electron_data(int tableSize,
		electron_data_inp& DataContainer,
		real energy_max,
		real energy_mean,
		real e_step,
		real eCutOff,
		const moller_XS& h2o_moller,  // Moller cross section for water
		const bhabha_XS& h2o_bhabha,  // Bhabha cross section for water
		const brems_XS&  h2o_brems,
		const real_3& voxel_size,
		result_type result);

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
  electron_data() {   // TODO : private and not used by anybody.
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

void electron_data_init(electron_data &A, unsigned int,
                        electron_data_inp *,
                        real, real, result_type, real);

#endif /* _ELECTRON_DATA_H_ */
