#ifndef _XVMC_H_
#define _XVMC_H_

/*****************************************************************************
 * xvmc.h: define global variables                                           *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 99/12/13        *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include "definitions.h"
#include "global.h"
#include "treatment_plan.h"
#include "beam_model.h"
#include "e_beam_model.h"
#include "p_beam_model.h"

// ****************************************
// define global variables
// ****************************************

// kinetic electron cut-off energy and low energy threshold for delta electron
// production (corresponds to ECUT-electron mass and AE-electron mass in EGS4)
real e_cut        = 0.5;        // (MeV)
real TWOe_cut     = e_cut*TWO;  // 2*e_cut (MeV)
real e_step       = 0.12;       // electron energy step size (default: 12%)
real p_cut        = 0.05;       // photon cut-off energy (MeV)
// KERMA approximation for primary photons with energy below k0_cut
real k0_cut       = 1.0;
// KERMA approximation for secondary photons with energy below k1_cut
real k1_cut       = 2.0;

extern treatment_plan plan;     // treatment plan parameters

// pointer to the beam model(s)
extern beam_model   *linac;
extern e_beam_model *elinac;
extern p_beam_model *plinac;

// pointer to the photon and electron beam modifiers
extern beam_modifier *pbeam_modifier;
extern beam_modifier *ebeam_modifier;

// ****************************************
// declare functions
// ****************************************

// general initializations
void init_xvmc(const char *, const char *);

// beam initializations
void init_beam(beam_core *);

// dose calculation algorithm
void calc_dose(beam_core *);

// read dose matrix file
void read_dose_file(const array_3d<float>  *, const array_3d<float>  *,
                    float &, int_3 &, float &,
                    const char *, const char *, const char *);

// evaluate and save dose distributions
void evaluate(treatment_plan &, beam_core *);
void evaluate(treatment_plan &);

#endif /* _XVMC_H_ */
