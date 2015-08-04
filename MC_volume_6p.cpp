/*****************************************************************************
 * MC_volume_6p.cpp:                                                         *
 *    class member functions for:                                            *
 *       MC_volume_6p:   region (volume) defined by 6 planes                 *
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

#include "MC_volume_6p.h"

// ****************************************
// member functions of class MC_volume_6p
// ****************************************

// initialize volume by 6 planes
void MC_volume_6p::init(MC_plane *plane0, MC_plane *plane1, MC_plane *plane2,
                        MC_plane *plane3, MC_plane *plane4, MC_plane *plane5)
{
   // assign plane pointers
   surface[0] = plane0;
   surface[1] = plane1;
   surface[2] = plane2;
   surface[3] = plane3;
   surface[4] = plane4;
   surface[5] = plane5;
}
