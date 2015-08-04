#ifndef _MC_MATERIAL_H_
#define _MC_MATERIAL_H_

/*****************************************************************************
 * MC_material.h:                                                            *
 *    class declarations and inline member functions for:                    *
 *       MC_material:    total and differential cross sections as well as    *
 *                       electron transport data for one specific material   *
 *                                                                           *
 * Copyright (C) 2003    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 14.07.2003      *
 *                                                                           *
 *****************************************************************************/

#include "definitions.h"
#include "electron_transport_data.h"
#include "compton_XS_diff.h"
#include "compton_XS_total.h"
#include "pair_XS_diff.h"
#include "pair_XS_total.h"
#include "photo_XS_total.h"

// *********************************************************************
// class MC_material: total and differential cross sections as well as
//                    electron transport data for one specific material
// *********************************************************************

class MC_material
{
   public:
      // allocate material data by the ESTAR (NIST, ICRU) electron
      // stopping power file, the NIST photon cross section file,
      // the differential Compton and pair cross section files
      MC_material(char *material_estar, char *material_nist,
                  char *compton_file,   char *pair_file);

      // delete material data
      ~MC_material(void);

      // Compton cross section
      compton_XS_total *tot_comp;
      compton_XS_diff  *compton;

      // pair cross section
      pair_XS_total *tot_pair;
      pair_XS_diff  *pair;

      // total photo cross section
      photo_XS_total *tot_phot;

      // electron transport data (stopping powers and ranges)
      electron_transport_data *e_data;
};

#endif /* _MC_MATERIAL_H_ */
