#ifndef _PAIR_XS_TOTAL_H_
#define _PAIR_XS_TOTAL_H_

/*****************************************************************************
 * pair_XS_total.h:                                                          *
 *    class declarations and inline member functions for:                    *
 *     class pair_XS_total: total pair cross section                         *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 14.12.1999      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************
#include <string>
#include "definitions.h"

/*****************************************************************************
 * class pair_XS_total:                                                      *
 *    total pair cross section data                                          *
 *****************************************************************************/

class PairTotData
{
  friend PairTotData* PairTotFactory(const std::string& file_name);
  friend class pair_XS_total;     // TODO goes away with inheritance

 public:

  PairTotData(int nEnergy, real energyMin, real energyMax);
  virtual ~PairTotData();

 protected:

  int      n_energy;       // number of bins for photon energy
  real     energy_min;     // minimum energy
  real     energy_max;     // maximum energy
  real*    sigma;          // pointer to total cross section data

 private:

  bool invariant() {return true;}       // TODO

  PairTotData();                               // not implemented
  PairTotData(const PairTotData&);             // not implemented
  PairTotData& operator=(const PairTotData&);  // not implemented
};

//###################################################################################

class pair_XS_total
{
 private:

  int      n_energy;       // number of bins for photon energy
  real     energy_min;     // minimum energy
  real     energy_max;     // maximum energy
  real     log_energy_min; // ln(minimum energy)
  real     log_energy_max; // ln(maximum energy)
  real     inv_delta;      // inverse bin size for photon energy (log scale)
  real    *sigma;          // pointer to total cross section data

 public:
  // default constructor (no total pair cross section data)
  pair_XS_total(void) :
    n_energy(0),
    sigma(NULL)
    { }

  // construct pair cross section by reading the NIST file
  pair_XS_total(char *file_name) { sigma = NULL; nist(file_name); }

  // input total pair data from NIST (XCOM) file
  void nist(char *);

  // input total pair data from file
  void read(char *);

  pair_XS_total(const PairTotData& );
  ~pair_XS_total(void);

  // get total pair cross section for the specified energy
  real get(real) const;
};

#endif /* _PAIR_XS_TOTAL_H_ */
