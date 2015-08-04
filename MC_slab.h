#ifndef _MC_SLAB_H_
#define _MC_SLAB_H_

/*****************************************************************************
 * MC_slab.h:                                                                *
 *    class declarations and inline member functions for:                    *
 *       MC_slab:        slab defined by two parallel planes                 *
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
#include "MC_region.h"

// ****************************************************************
// class MC_slab: slab defined by two parallel planes
// ****************************************************************

class MC_slab : public MC_region
{
   private:
      // initialize slab by two existing parallel planes
      void init(MC_plane *, MC_plane *);

   public:
      // define slab by two existing parallel planes and
      // the cross section file names
      MC_slab(MC_plane *plane0, MC_plane *plane1,
              char *material_estar, char *material_nist,
              char *compton_file,   char *pair_file)
         : MC_region(2,material_estar,material_nist,compton_file,pair_file)
              { init(plane0,plane1); }

      // define slab by two existing parallel planes and
      // a cross section data base for the slab material
      MC_slab(MC_plane *plane0, MC_plane *plane1, MC_material *material)
         : MC_region(2,material) { init(plane0,plane1); }

      // define slab perpendicular to the x, y or z axis and
      // the cross section file names
      MC_slab(const axis &xyz, const real &d0, const real &d1,
              char *material_estar, char *material_nist,
              char *compton_file,   char *pair_file);

      // define slab perpendicular to the x, y or z axis and
      // a cross section data base for the slab material
      MC_slab(const axis &xyz, const real &d0, const real &d1,
              MC_material *material);

      // define slab by four points in 3D space
      // the first three points define one plane, the fourth point is to
      // define the second plane parallel to plane one,
      // use cross section file names
      MC_slab(const real_3 &p0, const real_3 &p1,
              const real_3 &p2, const real_3 &p3,
              char *material_estar, char *material_nist,
              char *compton_file,   char *pair_file);

      // define slab by four points in 3D space
      // the first three points define one plane, the fourth point is to
      // define the second plane parallel to plane one,
      // use cross section data base for the slab material
      MC_slab(const real_3 &p0, const real_3 &p1,
              const real_3 &p2, const real_3 &p3,
              MC_material *material);
};

#endif /* _MC_SLAB_H_ */
