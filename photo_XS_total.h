#ifndef _PHOTO_XS_TOTAL_H_
#define _PHOTO_XS_TOTAL_H_

/*****************************************************************************
 * photo_XS_total.h:                                                         *
 *    class declarations and inline member functions for:                    *
 *     class photo_XS_total: total photo-absorption cross section            *
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

#include "definitions.h"

/*****************************************************************************
 * class photo_XS_total:                                                     *
 *    total photoelectric absorption cross section                           *
 *****************************************************************************/

class photo_XS_total
{
   private:
      int      n_energy;       // number of bins for photon energy
      real     energy_min;     // minimum energy
      real     energy_max;     // maximum energy
      real     log_energy_min; // ln(minimum energy)
      real     log_energy_max; // ln(maximum energy)
      real     inv_delta;      // inverse bin size for photon energy (log scale)
      real    *log_sigma;      // pointer to total cross section (log scale)

   public:
      // default constructor (no total photo cross section data)
      photo_XS_total(void) { n_energy = 0; log_sigma = NULL; }

      // construct photo cross section by reading the NIST file
      photo_XS_total(char *file_name) { log_sigma = NULL; nist(file_name); }

      // delete total photo data
      ~photo_XS_total(void);

      // input total photo data from NIST (XCOM) file
      void nist(char *);

      // get total photo cross section for the specified energy
      real get(real);
};

#endif /* _PHOTO_XS_TOTAL_H_ */
