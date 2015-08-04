#ifndef _GLOBAL_H_
#define _GLOBAL_H_

/*****************************************************************************
 * global.h: global constants, variables and functions                       *
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

#include <math.h>
#include "definitions.h"
#include "array_3d.h"
#include <cstring>

// ****************************************
// define global constants
// ****************************************

const real MINUS_ONE = -1.0;
const real ZERO      =  0.0;
const real ONE_HALF  =  0.5;
const real ONE       =  1.0;
const real TWO       =  2.0;
const real THREE     =  3.0;
const real SQRT2     =  sqrt(TWO);
const real HUGE_STEP =  100000.0;      // a huge particle step size

const real ONE_PI  = 3.14159265358979324;   // PI
const real TWO_PI  = 6.28318530717958648;   // 2*PI
const real HALF_PI = 1.57079632679489662;   // PI/2

const real SSD_0   = 100.0; // default source surface distance = 100 cm
                            // (the beam sizes are defined there)

const real TC     = 0.01;      // low energy threshold (kinetic energy in MeV)
      // for bremsstrahlung production (corresponds to AP in EGS4)

const real AMU2      = 7.5673E-5; // constant to calculate the screening angle
                                  // (screening mass squared)
const real BC        = 8778.2;    //
const real CHI_CC2   = 0.6626;    //
const real TEFF0_600 = 1.5954;    // electron step size restriction because of
                                  // the path length correction model

const real EMASS        = 0.5110034;   // electron mass
const real TWOxEMASS    = 1.0220068;   // 2*EMASS
const real EMASSxEMASS  = 0.2611245;   // EMASS*EMASS

// pi x electron radius squared in units of 10^(-23) cm^2
const real ONE_PIxR0xR0 = 0.024943591807;

#ifdef USE_SOBOL
// dimension of the Sobol' sequence generator
const int  SOBOL_DIM = 12;
#endif // USE_SOBOL

// ****************************************
// declare global variables
// ****************************************

// number of parallel processes (threads) to calculate dose,
// equal to the number of active processors in the system
extern int n_process;

// path to the input and output files (XVMC_WORK)
extern char *inout_path;

// kinetic electron cut-off energy and low energy threshold for delta electron
// production (corresponds to ECUT-electron mass and AE-electron mass in EGS4)
extern real e_cut;           // (MeV)
extern real TWOe_cut;        // 2*e_cut (MeV)

extern real e_step;          // electron energy step size (default: 12%)
extern real p_cut;           // photon cut-off energy (MeV)
// KERMA approximation for primary photons with energy below k0_cut
extern real k0_cut;
// KERMA approximation for secondary photons with energy below k1_cut
extern real k1_cut;

extern real_3 ref_point;     // reference point

extern real_3            voxel_size; // voxel sizes (cm)
extern int_3             dim;        // cube dimensions
extern real_3            cube_size;  // cube sizes (cm)
// density matrix
extern array_3d<float>  *density;    // density matrix
// density correction for the collison stopping power (eq. (19) in VMC paper I)
extern array_3d<float>  *dens_scol;
// low energy correction for "dens_scol"
extern array_3d<float>  *dens_ccol;
// density correction for radiation stopping power
// (eq. (20) and (24) in VMC paper I)
extern array_3d<float>  *dens_srad;
// density correction for multiple scattering (see VMC paper II)
extern array_3d<float>  *dens_fchi;
// density correction for Compton cross section = electron density
// (eq. (10) in XVMC paper I)
extern array_3d<float>  *dens_comp;
// density correction for pair cross section,
// similar to "dens_srad" (Feynman crossing)
extern array_3d<float>  *dens_pair;
// density correction for photo cross section
extern array_3d<float>  *dens_phot;

// dose matrix for present beam
extern array_3d<float>  *beam_dose;
// dose error matrix for present beam
extern array_3d<float>  *beam_error;
// total dose matrix
extern array_3d<float>  *sum_dose;
// total dose error matrix
extern array_3d<float>  *sum_error;

// Added by JOKim 15NOV2010 -------------------------
extern real_3            voxel_size_portal; // voxel sizes (cm)
extern int_3             dim_portal;        // cube dimensions
extern real_3            cube_size_portal;  // cube sizes (cm)

// density matrix
extern array_3d<float>  *density_portal;    // density matrix
// density correction for the collison stopping power (eq. (19) in VMC paper I)
extern array_3d<float>  *dens_scol_portal;
// low energy correction for "dens_scol"
extern array_3d<float>  *dens_ccol_portal;
// density correction for radiation stopping power
// (eq. (20) and (24) in VMC paper I)
extern array_3d<float>  *dens_srad_portal;
// density correction for multiple scattering (see VMC paper II)
extern array_3d<float>  *dens_fchi_portal;
// density correction for Compton cross section = electron density
// (eq. (10) in XVMC paper I)
extern array_3d<float>  *dens_comp_portal;
// density correction for pair cross section,
// similar to "dens_srad" (Feynman crossing)
extern array_3d<float>  *dens_pair_portal;
// density correction for photo cross section
extern array_3d<float>  *dens_phot_portal;

// dose matrix for present beam
extern array_3d<float>  *beam_dose_portal;
// dose error matrix for present beam
extern array_3d<float>  *beam_error_portal;
// total dose matrix
extern array_3d<float>  *sum_dose_portal;
// total dose error matrix
extern array_3d<float>  *sum_error_portal;
// End of Added -------------------------------------

#ifdef CHECK_ENERGY
extern sum_energy_type tot_energy; // test energy conservation
#endif // CHECK_ENERGY

// ****************************************
// declare global functions
// ****************************************

// error and message handlers
void xvmc_error(const char *, const char *, const int);
void xvmc_warning(const char *, const char *, const int);
void xvmc_warning(const char *, const bool, const int);
void xvmc_warning(const char *, const int, const int);
void xvmc_warning(const char *, const real, const int);
void xvmc_message(const char *, const int);
void xvmc_message(const char *, const char *, const int);
void xvmc_message(const char *, const int , const char *, const int);
void xvmc_message(const char *, const long, const char *, const int);
void xvmc_message(const char *, const float, const char *, const int);
void xvmc_message(const char *, const double, const char *, const int);
void xvmc_message(const char *, const float, const char *, const float,
                  const char *, const int);
void xvmc_message(const char *, const double, const char *, const double,
                  const char *, const int);
void xvmc_message(const char *, const float, const char *, const float,
                  const char *, const float, const char *, const int);
void xvmc_message(const char *, const double, const char *, const double,
                  const char *, const double, const char *, const int);
#ifdef CHECK_ENERGY
void energy_message(const char *, const double, const char *, const int);
#endif // CHECK_ENERGY
void get_index(int, real *, real, int &, int &);
real interpolate(int, real *, real *, real);
void rotate(real, real, real_3 &);
void rotate(real, real, real, real, real_3 &);
void rotate(float, float, float, float, float &, float &, float &);
void rotate(double, double, double, double, double &, double &, double &);
char *get_file_path(const char *, const char *);

// ****************************************
// inline functions
// ****************************************

inline float  min_of(float a, float b) { return( a <= b ? a : b); }
inline double min_of(double a, double b) { return( a <= b ? a : b); }
inline int    min_of(int  a, int  b) { return( a <= b ? a : b); }

inline float  max_of(float a, float b) { return( a >= b ? a : b); }
inline double max_of(double a, double b) { return( a >= b ? a : b); }
inline int    max_of(int  a, int  b) { return( a >= b ? a : b); }

#endif /* _GLOBAL_H_ */
