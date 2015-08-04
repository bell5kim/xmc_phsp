#ifndef _ELECTRON_DATA_INP_H_
#define _ELECTRON_DATA_INP_H_

/*****************************************************************************
 * electron_data_inp.h:                                                      *
 *    class declarations and inline member functions for:                    *
 *       electron_data_inp: input electron transport data                    *
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

class electron_data; // To fix higher g++ version problem

// input electron transport data as a function of energy
class electron_data_inp
{
  friend class electron_data;
  friend void  electron_data_init(electron_data &, unsigned int,
				  electron_data_inp *,
				  real, real, result_type, real = 1.0);

  friend electron_data_inp* ElectronDataFactory(const char* filename);

 public:

  electron_data_inp(int nTableSize);

 protected:

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

#endif /* _ELECTRON_DATA_INP_H_ */
