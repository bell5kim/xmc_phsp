#ifndef _PHOTON_DATA_INP_H_
#define _PHOTON_DATA_INP_H_

/*****************************************************************************
 * photon_data_inp.h:                                                        *
 *    class declarations and inline member functions for:                    *
 *       photon_data_inp: input photon transport data                        *
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

class photon_data; // To fix higher version of g++

// input photon transport data as a function of energy
class photon_data_inp
{
  friend class photon_data;
  friend void  photon_data_init(photon_data &, unsigned int,
				photon_data_inp *, real);

  friend photon_data_inp* PhotonDataFactory(const char* filename);

 public:

  photon_data_inp(int nTableSize);

 protected:

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

#endif /* _PHOTON_DATA_INP_H_ */
