#ifndef _ELECTRON_TRANSPORT_DATA_H_
#define _ELECTRON_TRANSPORT_DATA_H_

/*****************************************************************************
 * electron_transport_data.h:                                                *
 *    class declarations and inline member functions for:                    *
 *       electron_transport_data: electron stopping powers and ranges        *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 06.12.2001      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include "definitions.h"

/*****************************************************************************
 * class electron_transport_data:                                            *
 *    electron stopping powers and ranges                                    *
 *****************************************************************************/

class electron_transport_data
{
   private:
      int      n_energy;       // number of bins for electron energy
      real     energy_min;     // minimum energy
      real     energy_max;     // maximum energy
      real     log_energy_min; // ln(minimum energy)
      real     log_energy_max; // ln(maximum energy)
      real     inv_delta;      // inverse bin size for electron energy (log)
      real    *log_dedx;       // pointer to total stopping power (log scale)
      real    *log_csda;       // pointer to CSDA range (log scale)

   public:
      // default constructor (no electron stopping powers and ranges)
      electron_transport_data(void) {
         n_energy = 0; log_dedx = NULL; log_csda = NULL; }

      // construct electron transport data by reading the ESTAR file
      electron_transport_data(char *file_name) {
         log_dedx = NULL; log_csda = NULL; estar(file_name); }

      // delete electron transport data data
      ~electron_transport_data(void);

      // input electron transport data from ESTAR (NIST, ICRU) file
      void estar(char *);

      // get total stopping power for the specified energy
      real get_dedx(real);

      // get continuous slowing down (CSDA) range for the specified energy
      real get_csda(real);
};

#endif /* _ELECTRON_TRANSPORT_DATA_H_ */
