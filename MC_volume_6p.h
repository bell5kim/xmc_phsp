#ifndef _MC_VOLUME_6P_H_
#define _MC_VOLUME_6P_H_

/*****************************************************************************
 * MC_volume_6p.h:                                                           *
 *    class declarations and inline member functions for:                    *
 *       MC_volume_6p:   region (volume) defined by 6 planes                 *
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
// class MC_volume_6p: region (volume) defined by 6 planes
// ****************************************************************

class MC_volume_6p : public MC_region
{
   private:
      // initialize volume by 6 planes
      void init(MC_plane *, MC_plane *, MC_plane *,
                MC_plane *, MC_plane *, MC_plane *);

   public:
      // define volume by 6 planes and the cross section file names
      MC_volume_6p(MC_plane *plane0, MC_plane *plane1,  MC_plane *plane2,
                   MC_plane *plane3, MC_plane *plane4,  MC_plane *plane5,
                   const real_3 &p_inp,
                   char *material_estar, char *material_nist,
                   char *compton_file,   char *pair_file)
         : MC_region(6,material_estar,material_nist,compton_file,pair_file) {
              p_ref.x = p_inp.x; p_ref.y = p_inp.y; p_ref.z = p_inp.z;
              init(plane0,plane1,plane2,plane3,plane4,plane5); }

      // define volume by 6 planes and a cross section
      // data base for the volume material
      MC_volume_6p(MC_plane *plane0, MC_plane *plane1,  MC_plane *plane2,
                   MC_plane *plane3, MC_plane *plane4,  MC_plane *plane5,
                   const real_3 &p_inp, MC_material *material)
         : MC_region(6,material) {
              p_ref.x = p_inp.x; p_ref.y = p_inp.y; p_ref.z = p_inp.z;
              init(plane0,plane1,plane2,plane3,plane4,plane5); }
};

#endif /* _MC_VOLUME_6P_H_ */
