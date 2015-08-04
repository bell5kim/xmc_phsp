#ifndef _MC_REGION_H_
#define _MC_REGION_H_

/*****************************************************************************
 * MC_region.h:                                                              *
 *    class declarations and inline member functions for:                    *
 *       MC_region:      geometrical region defined by an array of planes    *
 *                                                                           *
 * Copyright (C) 2001    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 29.08.2001      *
 *                                                                           *
 *****************************************************************************/

#include "definitions.h"
#include "MC_material.h"
#include "MC_plane.h"

// *******************************************************************
// class MC_region: geometrical region defined by an array of planes
// *******************************************************************

class MC_region
{
   public:
      // allocate new region with a given number of surface planes,
      // a material specified by the ESTAR (NIST, ICRU) electron
      // stopping power file, the NIST photon cross section file,
      // the differential Compton and pair cross section files
      MC_region(unsigned n_planes,
                char *material_estar, char *material_nist,
                char *compton_file,   char *pair_file);

      // allocate new region with a given number of surface planes,
      // and a cross section data base for one specific material
      MC_region(unsigned n_planes, MC_material *ini_material);

      // delete region
      ~MC_region(void);

      // get number of planes defining the region
      unsigned  get_num_planes(void) { return(num_planes); }

      // get the reference point
      real_3 get_p_ref(void) { return(p_ref); }

      // calculate ralationship (<0, =0 or >0) of a point to the region
      // +1: the point is outside the region
      //  0: the point is located at one of the surface planes
      // -1: the point is inside
      // inline int relationship(const real_3 &);

      // calculate distance of a particle to the region surface by taking
      // into account every surface plane and the present particle plane given
      // by the plane pointer, this pointer also returns the target plane
      inline real distance(const particle_parameters &, MC_plane *&);

      // array of pointers to plane defining the surface of the region
      MC_plane    **surface;

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

   protected:
      // number of planes
      unsigned     num_planes;

      // to identify the inside of the region we need a reference point within
      // the region to calculate the reference point to surface plane
      // relationships
      real_3       p_ref;

      // pointer to the cross section data base for this material,
      // if NULL, the region was constructed using the cross section
      // file names
      MC_material *material;
};

// *******************************************
// inline member functions of class MC_region
// *******************************************

// calculate distance of a particle to the region surface by taking
// into account every surface plane and the present particle plane given
// by the plane pointer, this pointer also returns the target plane
inline real MC_region::distance(const particle_parameters &p,
                                MC_plane *&present_plane)
{
   // constants
   const real MC_INFINITY = 1.0e+10;

   // if the number of planes is 0, i.e. if the region is an exterior region,
   // return infinity distance and plane pointer to NULL
   if (num_planes == 0)
   {
      present_plane = NULL;
      return(MC_INFINITY);
   }

   // the number of planes is larger than 0
   real       dist_min  = MC_INFINITY;
   MC_plane  *new_plane = NULL;
   for (register unsigned i=0; i<num_planes; ++i)
   {
      // old version: if (surface[i] != present_plane)
      // old version: {
      // calculate linear distance of the particle to plane i
      real distance = surface[i]->distance(p,present_plane);
      if (distance < dist_min)
      {
         dist_min  = distance;
         new_plane = surface[i];
      }
      // old version: }
   }
   present_plane = new_plane;
   return(dist_min);
}

#endif /* _MC_REGION_H_ */
