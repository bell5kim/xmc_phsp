/*****************************************************************************
 * MC_material.cpp:                                                          *
 *    class member functions for:                                            *
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

// ****************************************
// includes
// ****************************************

#include <new>
using namespace std;

#include "MC_material.h"

// ****************************************
// member functions of class MC_material
// ****************************************

// allocate material data by the ESTAR (NIST, ICRU) electron
// stopping power file, the NIST photon cross section file,
// the differential Compton and pair cross section files
MC_material::MC_material(char *material_estar, char *material_nist,
                         char *compton_file,   char *pair_file)
{
   // create total Compton cross section for the specified material
   tot_comp = NULL;
   if ( (tot_comp = new (nothrow) compton_XS_total(material_nist)) == NULL )
   {
      xvmc_error("MC_material::MC_material",
                 "cannot create total Compton cross section",8);
   }

   // create differential Compton cross section using the Compton data file
   compton = NULL;
   if ( (compton = new (nothrow) compton_XS_diff(compton_file)) == NULL )
   {
      xvmc_error("MC_material::MC_material",
                 "cannot create differential Compton cross section",8);
   }

   // create total pair cross section for the specified material
   tot_pair = NULL;
   if ( (tot_pair = new (nothrow) pair_XS_total(material_nist)) == NULL )
   {
      xvmc_error("MC_material::MC_material",
                 "cannot create total pair cross section",8);
   }

   // create differential pair cross section using the pair data file
   pair = NULL;
   if ( (pair = new (nothrow) pair_XS_diff(pair_file)) == NULL )
   {
      xvmc_error("MC_material::MC_material",
                 "cannot create differential pair cross section",8);
   }

   // create total photo cross section for the specified material
   tot_phot = NULL;
   if ( (tot_phot = new (nothrow) photo_XS_total(material_nist)) == NULL )
   {
      xvmc_error("MC_material::MC_material",
                 "cannot create total photo cross section",8);
   }

   // create electron transport data for the specified material
   e_data = NULL;
   if ( (e_data = new (nothrow) electron_transport_data(material_estar))==NULL )
   {
      xvmc_error("MC_material::MC_material",
                 "cannot create electron transport data",8);
   }
}

// delete material data
MC_material::~MC_material(void)
{
   // delete Compton cross section
   delete tot_comp; tot_comp = NULL;
   delete compton;  compton  = NULL;

   // delete pair cross section
   delete tot_pair; tot_pair = NULL;
   delete pair;     pair     = NULL;

   // delete photo cross section
   delete tot_phot; tot_phot = NULL;

   // delete electron transport data
   delete e_data; e_data = NULL;
}
