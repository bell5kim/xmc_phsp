#ifndef _MC_DOUBLESLAB_H_
#define _MC_DOUBLESLAB_H_

/*****************************************************************************
 * MC_doubleslab.h:                                                          *
 *    class declarations and inline member functions for:                    *
 *       MC_doubleslab:  two slabs defined by 3 parallel planes              *
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
#include "MC_slab.h"
#include "MC_object.h"

// ****************************************************************
// class MC_doubleslab: two slabs defined by 3 parallel planes
// ****************************************************************

class MC_doubleslab : public MC_object
{
   private:
      // initialize doubleslab by three existing parallel planes
      // and two slab materials
      void init(MC_plane *, MC_plane *, MC_plane *, char *, char *);

   public:
      // define doubleslab by three existing parallel planes and two materials
      MC_doubleslab(MC_plane *plane0, MC_plane *plane1, MC_plane *plane2,
                    char *material0,  char *material1)
         : MC_object(3,2) { init(plane0,plane1,plane2,material0,material1); }

      // define double slab perpendicular to the x, y or z axis
      // and two materials
      MC_doubleslab(const axis &, const real &, const real &, const real &,
                    char *,  char *);

      // define doubleslab by five points in 3D space and two materials,
      // the first three points define one plane, the two other points are to
      // define two further planes parallel to plane one
      MC_doubleslab(const real_3 &, const real_3 &, const real_3 &,
                    const real_3 &, const real_3 &,
                    char *,  char *);
};

#endif /* _MC_DOUBLESLAB_H_ */
