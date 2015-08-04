#ifndef _MULTI_ELECTRON_H_
#define _MULTI_ELECTRON_H_

/*****************************************************************************
 * multi_electron.h                                                          *
 *    class declaration for                                                  *
 *       multi_electron: create and simulate electron histories              *
 *                       (history repetition technique)                      *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 00/02/16        *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include "definitions.h"
#include "ranmar.h"
#include "global.h"
#include "portal_dose.h"

// ****************************************
// class multi_electron
// ****************************************

class multi_electron
{
 public:

  friend void multi_electron_create(multi_electron &, real, ranmar &);
  friend void multi_electron_trace(multi_electron &,
                                   particle_parameters &, ranmar &,
#ifdef CHECK_ENERGY
                 sum_energy_type &,
#endif // CHECK_ENERGY
                 array_3d<double> *, portal_dose *);

  friend void multi_positron_create(multi_electron &, real, ranmar &);
  friend void multi_positron_trace(multi_electron &,
                                   particle_parameters &, ranmar &,
#ifdef CHECK_ENERGY
                 sum_energy_type &,
#endif // CHECK_ENERGY
                 array_3d<double> *, portal_dose *);
// --- Added by JOKim 15Nov2010 --------------------------------------
  friend void multi_electron_trace_portal(multi_electron &,
                                   particle_parameters &, ranmar &,
#ifdef CHECK_ENERGY
                 sum_energy_type &,
#endif // CHECK_ENERGY
                 array_3d<double> *);

  friend void multi_positron_trace_portal(multi_electron &,
                                   particle_parameters &, ranmar &,
#ifdef CHECK_ENERGY
                 sum_energy_type &,
#endif // CHECK_ENERGY
                 array_3d<double> *);
// --- End of Adding -------------------------------------------------
  friend class xvMC;

 private:
      // electron step data
      enum    {MAX_STEP=128};           // array size (maximum step number)
      int     n_step;                   // number of electron steps
      real    energy_loss[MAX_STEP];    // energy loss per step
      real    radiation_loss[MAX_STEP]; // bremsstrahlung loss per step
      real    sin_brems[MAX_STEP];      // sin of bremsstrahlung photon angle
      real    cos_brems[MAX_STEP];      // cos of bremsstrahlung photon angle
      real    dose_cor[MAX_STEP];       // correction for dose, ionization etc.
      real    random_number[MAX_STEP];  // divide the step into 2 sub-steps
      real    step_size[MAX_STEP];      // the length of this step
      real    s_cor[MAX_STEP];          // low energy correction for s_col
      real    alpha_r[MAX_STEP];        // alpha_r for this step
      real    reduced_angle[MAX_STEP];  // "reduced" MS (mult. scat.) angle
      real    sin_theta[MAX_STEP];      // direction change of the primary el.
      real    cos_theta[MAX_STEP];      // according to a discrete interaction
      int     i_delta[MAX_STEP];        // index of the delta electron
      bool    moller[MAX_STEP];         // true if there is a Moller interaction
                                        // (delta electron) during the step
      bool    bhabha[MAX_STEP];         // true if there is a Bhabha interaction
                                        // (delta electron) during the step

      // stack for delta (secondary) electron data
      enum    {MAX_DELTA=100};         // array size (maximum number)
      int     n_delta;                 // number of delta electrons in stack
      real    energy_delta[MAX_DELTA]; // delta electron energy
      real    sin_theta_d[MAX_DELTA];  // sin of the longitudinal scat. angle
      real    cos_theta_d[MAX_DELTA];  // cos of the longitudinal scat. angle
};

// create one electron history in water
void multi_electron_create(multi_electron &, real, ranmar &);

// trace pre-calculated electron history through the calculation cube
void multi_electron_trace(multi_electron &, particle_parameters &, ranmar &,
#ifdef CHECK_ENERGY
			  sum_energy_type &,
#endif // CHECK_ENERGY
			  array_3d<double> *, portal_dose *);

// create one electron history in water
void multi_positron_create(multi_electron &, real, ranmar &);

// trace pre-calculated electron history through the calculation cube
void multi_positron_trace(multi_electron &, particle_parameters &, ranmar &,
#ifdef CHECK_ENERGY
			  sum_energy_type &,
#endif // CHECK_ENERGY
			  array_3d<double> *, portal_dose *);

// --- Added by JOKim 15Nov2010 --------------------------------------
// trace pre-calculated electron history through the calculation cube
void multi_electron_trace_portal(multi_electron &, particle_parameters &, ranmar &,
#ifdef CHECK_ENERGY
			  sum_energy_type &,
#endif // CHECK_ENERGY
			  array_3d<double> *);

// trace pre-calculated electron history through the calculation cube
void multi_positron_trace_portal(multi_electron &, particle_parameters &, ranmar &,
#ifdef CHECK_ENERGY
			  sum_energy_type &,
#endif // CHECK_ENERGY
			  array_3d<double> *);
// --- End of Adding -------------------------------------------------

#endif  /* _MULTI_ELECTRON_H_ */
