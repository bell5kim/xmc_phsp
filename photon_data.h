#ifndef _PHOTON_DATA_H_
#define _PHOTON_DATA_H_

/*****************************************************************************
 * photon_data.h:                                                            *
 *    class declarations and inline member functions for:                    *
 *       photon_data:    photon transport data (function of energy)          *
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
#include "photon_data_inp.h"
#include "compton_XS_total.h"
#include "pair_XS_total.h"

// photon transport data as a function of energy
class photon_data
{
 public:
  photon_data(int nTableSize,
	      photon_data_inp& inp,
	      real energy_max,
	      real p_cut,
	      const compton_XS_total& tot_comp, // total compton cross section for water
	      const pair_XS_total& tot_pair);

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

void photon_data_init(photon_data &A, unsigned int,
                      photon_data_inp *, real);

#endif /* _PHOTON_DATA_H_ */
