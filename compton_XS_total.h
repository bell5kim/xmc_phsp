#ifndef _COMPTON_XS_TOTAL_H_
#define _COMPTON_XS_TOTAL_H_

/*****************************************************************************
 * compton_XS_total.h:                                                       *
 *     class compton_XS_total: total Compton cross section                   *
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
 * class compton_XS_total:                                                   *
 *    total compton cross section                                            *
 *****************************************************************************/

class compton_XS_total
{
   public:
      // construct the Compton cross section using a given electron density
      //    default: water
      compton_XS_total(real edens=3.343) : factorc(edens * ONE_PIxR0xR0) { }

      // construct the Compton cross section by reading the electron density
      // from the cross section file
      // attention: we do not use the Compton data from this file!
      compton_XS_total(char *file_name);

      // initialize the material dependent cross section factor
      //   edens: electron density in units of 10^23 cm^-3
      void init(real edens) { factorc = edens * ONE_PIxR0xR0; }

      // get total Compton cross section for the specified energy
      real get(real energy) const;

   private:
      // material dependent factor (n_e x pi x r_0 x r_0)
      real     factorc;
};

#endif /* _COMPTON_XS_TOTAL_H_ */
