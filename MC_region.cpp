/*****************************************************************************
 * MC_region.cpp:                                                            *
 *    class member functions for:                                            *
 *       MC_region:      geometrical region defined by an array of planes    *
 *                                                                           *
 * Copyright (C) 2001    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 31.08.2001      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <new>
using namespace std;

#include "MC_region.h"

// ****************************************
// member functions of class MC_region
// ****************************************

// allocate new region with a given number of surface planes,
// a material specified by the ESTAR (NIST, ICRU) electron
// stopping power file, the NIST photon cross section file,
// the differential Compton and pair cross section files
MC_region::MC_region(unsigned n_planes,
                     char *material_estar, char *material_nist,
                     char *compton_file,   char *pair_file)
{
   num_planes = n_planes;

   // reference point set to zero
   p_ref.x = ZERO;
   p_ref.y = ZERO;
   p_ref.z = ZERO;

   // a cross section data base is not used,
   // we create the cross section data below
   material = NULL;

   // create array of pointers to the surface planes
   surface = NULL;
   if (num_planes > 0)
   {
      if ((surface = new (nothrow) MC_plane*[num_planes]) == NULL )
      {
         xvmc_error("MC_region::MC_region",
                    "cannot allocate memory for surface pointer array",8);
      }
      for (register unsigned i=0; i<num_planes; ++i) surface[i] = NULL;
   }

   // create total Compton cross section for the specified material
   tot_comp = NULL;
   if ( (tot_comp = new (nothrow) compton_XS_total(material_nist)) == NULL )
   {
      xvmc_error("MC_region::MC_region",
                 "cannot create total Compton cross section",8);
   }

   // create differential Compton cross section using the Compton data file
   compton = NULL;
   if ( (compton = new (nothrow) compton_XS_diff(compton_file)) == NULL )
   {
      xvmc_error("MC_region::MC_region",
                 "cannot create differential Compton cross section",8);
   }

   // create total pair cross section for the specified material
   tot_pair = NULL;
   if ( (tot_pair = new (nothrow) pair_XS_total(material_nist)) == NULL )
   {
      xvmc_error("MC_region::MC_region",
                 "cannot create total pair cross section",8);
   }

   // create differential pair cross section using the pair data file
   pair = NULL;
   if ( (pair = new (nothrow) pair_XS_diff(pair_file)) == NULL )
   {
      xvmc_error("MC_region::MC_region",
                 "cannot create differential pair cross section",8);
   }

   // create total photo cross section for the specified material
   tot_phot = NULL;
   if ( (tot_phot = new (nothrow) photo_XS_total(material_nist)) == NULL )
   {
      xvmc_error("MC_region::MC_region",
                 "cannot create total photo cross section",8);
   }

   // create electron transport data for the specified material
   e_data = NULL;
   if ( (e_data = new (nothrow) electron_transport_data(material_estar))==NULL )
   {
      xvmc_error("MC_region::MC_region",
                 "cannot create electron transport data",8);
   }
}

// allocate new region with a given number of surface planes,
// and a cross section data base for one specific material
MC_region::MC_region(unsigned n_planes, MC_material *ini_material)
{
   num_planes = n_planes;

   // reference point set to zero
   p_ref.x = ZERO;
   p_ref.y = ZERO;
   p_ref.z = ZERO;

   // a cross section data base is used
   material = ini_material;

   // create array of pointers to the surface planes
   surface = NULL;
   if (num_planes > 0)
   {
      if ((surface = new (nothrow) MC_plane*[num_planes]) == NULL )
      {
         xvmc_error("MC_region::MC_region",
                    "cannot allocate memory for surface pointer array",8);
      }
      for (register unsigned i=0; i<num_planes; ++i) surface[i] = NULL;
   }

   // total Compton cross section for the specified material
   tot_comp = material->tot_comp;

   // differential Compton cross section using the Compton data file
   compton = material->compton;

   // total pair cross section for the specified material
   tot_pair = material->tot_pair;

   // differential pair cross section using the pair data file
   pair = material->pair;

   // total photo cross section for the specified material
   tot_phot = material->tot_phot;

   // electron transport data for the specified material
   e_data = material->e_data;
}

// delete region
MC_region::~MC_region(void)
{
   // delete surface planes
   if (num_planes > 0)
   {
      for (register unsigned i=0; i<num_planes; ++i)
      {
         if (surface[i] != NULL) { delete surface[i]; surface[i] = NULL; }
      }
   }

   // delete pointer array
   delete [] surface; surface = NULL;

   // free memory for cross section data if no data base
   // for the material (MC_material) was used
   if (material == NULL)
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
   else
   {
      material = NULL;
   }

   // set everything to zero
   num_planes = 0;
   p_ref.x = ZERO;
   p_ref.y = ZERO;
   p_ref.z = ZERO;
}
