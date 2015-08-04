#ifndef _INTERACTIONS_H_
#define _INTERACTIONS_H_

/*****************************************************************************
 * interactions.h:                                                           *
 *    class declaration:                                                     *
 *     class electron_transport_data: electron stopping powers and ranges    *
 *     class moller_XS:               Moller cross section                   *
 *     class bhabha_XS:               Bhabha cross section                   *
 *     class brems_XS:                bremsstrahlung cross section           *
 *     class compton_XS_diff:         differential Compton cross section     *
 *     class compton_XS_total:        total Compton cross section            *
 *     class pair_XS_diff:            differential pair cross section        *
 *     class pair_XS_total:           total pair cross section               *
 *     class photo_XS_total:          total photo-absorption cross section   *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 14.12.1999      *
 *    class electron_transport_data                       MF 06.12.2001      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include "definitions.h"
#include "rndm.h"

#define MOLLER_EXACT

/*****************************************************************************
 * class electron_transport_data:                                            *
 *    electron stopping powers and ranges                                    *
 *****************************************************************************/

class electron_transport_data
{
   private:
      int      n_energy;       // number of bins for electron energy
      real     energy_min;     // minimum energy
      real     energy_max;     // maximum energy
      real     log_energy_min; // ln(minimum energy)
      real     log_energy_max; // ln(maximum energy)
      real     inv_delta;      // inverse bin size for electron energy (log)
      real    *log_dedx;       // pointer to total stopping power (log scale)
      real    *log_csda;       // pointer to CSDA range (log scale)

   public:
      // default constructor (no electron stopping powers and ranges)
      electron_transport_data(void) {
         n_energy = 0; log_dedx = NULL; log_csda = NULL; }

      // construct electron transport data by reading the ESTAR file
      electron_transport_data(char *file_name) {
         log_dedx = NULL; log_csda = NULL; estar(file_name); }

      // delete electron transport data data
      ~electron_transport_data(void);

      // input electron transport data from ESTAR (NIST, ICRU) file
      void estar(char *);

      // get total stopping power for the specified energy
      real get_dedx(real);

      // get continuous slowing down (CSDA) range for the specified energy
      real get_csda(real);
};

/*****************************************************************************
 * class moller_XS:                                                          *
 *    Moller cross section                                                   *
 *****************************************************************************/

class moller_XS
{
   private:
      real     factor1;    // material dependent factor (n_e x pi x r_0 x r_0)
      real     factor2;    // additional factor: 2m/(beta**2), beta = v/c
      real     t;          // total electron energy, kinetic (e) + rest (EMASS)
      real     xi1,xi1p;   // see SLAC-265 page 58
      real     c1,c2;      // see SLAC-265 page 58
   public:
      // construct the material dependent cross section factor
      //   default: water, edens: electron density in units of 10^23 cm^-3
      moller_XS(real edens=3.343) { factor1 = edens * ONE_PIxR0xR0; }

      // calculate total Moller cross section for the specified energy
      real total(real &);

      // calculate first moment for the given energy
      real first(real &);

      // sample Moller interaction parameters
      inline bool interaction(real, ranmar &,
                              real &, real &, real &,
                              real &, real &, real &);
};

// sample Moller interaction parameters
inline bool moller_XS::interaction(real energy, ranmar  &rndm,
                             real &energy_n, real &cos_t_n, real &sin_t_n,
                             real &energy_d, real &cos_t_d, real &sin_t_d)
{
   // if the energy is too small we don't have a Moller interaction
   if (energy <= TWOe_cut)
   {
      energy_n = energy;
      cos_t_n  = ONE;
      sin_t_n  = ZERO;
      energy_d = ZERO;
      cos_t_d  = ZERO;
      sin_t_d  = ONE;
      return(false);
   }

#ifdef MOLLER_EXACT
   // sample secondary electron energy from the exact Moller cross section
   // using a method described in the EGSnrc manual (PIRS-701)
   t        = energy + EMASS;
   real c2e = EMASS*(energy+t)/(energy*t*t); // c2/energy
   real c3  = ONE/(TWO*t*t);                 // 1/(2*t*t)
   // maximum to normalize the rejection function
   real gmax = 3.0*energy*energy/(2.0*t*t);

   real rej_weight; // rejection weight
   do
   {
      // pick a random number
      real eta = rndm.number();

      // sample delta energy from (1/energy_d)^2 distribution
      energy_d = (energy-e_cut)*e_cut/((energy-e_cut)*(ONE-eta)+e_cut*eta);

      // calculate rejection weight (function)
      rej_weight = (ONE+(c3*energy_d-c2e)*energy_d)/gmax;
   }
   while (rndm.number() > rej_weight);

   // the secondary electron is the electron with smaller energy
   energy_d = min_of(energy_d,energy-energy_d);
#else
   // sample energy from approximated Moller cross section (1/E^2)
   real eta = rndm.number();
   energy_d = energy*e_cut/(energy*(ONE-eta)+TWOe_cut*eta);
#endif /* MOLLER_EXACT */

   // new energy of the primary electron
   energy_n = energy - energy_d;

   // calculate polar scattering angle of the delta electron
   cos_t_d = energy_d*(energy+TWOxEMASS)/energy/(energy_d+TWOxEMASS);
   if (cos_t_d < ONE)
   {
      sin_t_d = sqrt(ONE-cos_t_d);
      cos_t_d = sqrt(cos_t_d);
   }
   else
   {
      sin_t_d = ZERO;
      cos_t_d = ONE;
   }

   // calculate polar scattering angle of the primary electron
   cos_t_n = energy_n*(energy+TWOxEMASS)/energy/(energy_n+TWOxEMASS);
   if (cos_t_n < ONE)
   {
      sin_t_n = sqrt(ONE-cos_t_n);
      cos_t_n = sqrt(cos_t_n);
   }
   else
   {
      sin_t_n = ZERO;
      cos_t_n = ONE;
   }

   return(true);
}

/*****************************************************************************
 * class bhabha_XS:                                                          *
 *    Bhabha cross section                                                   *
 *****************************************************************************/

class bhabha_XS
{
   private:
      // material dependent factor (n_e x 2 x pi x r_0 x r_0 x m)
      real factor0;

      real tau,yy,inv_beta2,b0,b1,b2,b3,b4,xi,inv_xi;

   public:
      // construct the material dependent cross section factor
      //   default: water, edens: electron density in units of 10^23 cm^-3
      bhabha_XS(real edens=3.343) {
         factor0 = edens * TWOxEMASS * ONE_PIxR0xR0; }

      // calculate total Bhabha cross section for the specified energy
      real total(real &);

      // sample Bhabha interaction parameters
      inline bool interaction(real, ranmar &,
                              real &, real &, real &,
                              real &, real &, real &);
};

// sample Bhabha interaction parameters
inline bool bhabha_XS::interaction(real energy, ranmar  &rndm,
                             real &energy_n, real &cos_t_n, real &sin_t_n,
                             real &energy_d, real &cos_t_d, real &sin_t_d)
{
   // if the energy is too small we don't have a Bhabha interaction
   if (energy <= TWOe_cut)
   {
      energy_n = energy;
      cos_t_n  = ONE;
      sin_t_n  = ZERO;
      energy_d = ZERO;
      cos_t_d  = ZERO;
      sin_t_d  = ONE;
      return(false);
   }

   // sample energy from approximated Bhabha cross section (1/E^2)
   real eta = rndm.number();
   energy_d = energy*e_cut/(energy*(ONE-eta)+TWOe_cut*eta);

   // new energy of the positron
   energy_n = energy - energy_d;

   // calculate polar scattering angle of the delta electron
   cos_t_d = energy_d*(energy+TWOxEMASS)/energy/(energy_d+TWOxEMASS);
   if (cos_t_d < ONE)
   {
      sin_t_d = sqrt(ONE-cos_t_d);
      cos_t_d = sqrt(cos_t_d);
   }
   else
   {
      sin_t_d = ZERO;
      cos_t_d = ONE;
   }

   // calculate polar scattering angle of the positron
   cos_t_n = energy_n*(energy+TWOxEMASS)/energy/(energy_n+TWOxEMASS);
   if (cos_t_n < ONE)
   {
      sin_t_n = sqrt(ONE-cos_t_n);
      cos_t_n = sqrt(cos_t_n);
   }
   else
   {
      sin_t_n = ZERO;
      cos_t_n = ONE;
   }

   return(true);
}

/*****************************************************************************
 * class brems_XS:                                                           *
 *    bremsstrahlung cross section                                           *
 *****************************************************************************/

class brems_XS
{
   private:
      real     factorb;    // material dependent factor

   public:
      // construct the material dependent cross section factor
      //   default: water
      brems_XS(real rad_length=36.0863) { factorb = 4.0/(3.0*rad_length); }

      // calculate total bremsstrahlung cross section for the specified energy
      real total(real &);

      // calculate first moment for the given energy
      real first(real &);

      // sample bremsstrahlung photon energy
      inline real sample(real &, ranmar &);

      // sample bremsstrahlung interaction parameters
      inline bool interaction(real, ranmar &,
                              real &, real &, real &,
                              real &, real &, real &);
};

// sample bremsstrahlung photon energy
inline real brems_XS::sample(real &e, ranmar &rndm)
{
   if (e <= TC) return(ZERO);

#ifdef BREMS_EXACT
   real e2         = e*e;     // energy squared
#endif
   real e_brems    = ZERO;    // photon energy
   real rej_weight = ONE;     // rejection weight
   do
   {
      e_brems    = TC*exp(rndm.number()*log(e/TC));
#ifdef BREMS_EXACT
      rej_weight = ONE - e_brems/e + 0.75*e_brems*e_brems/e2;
#else
      rej_weight = ONE - e_brems/e;
#endif
   }
   while (rndm.number() > rej_weight);

   return(e_brems);
}

// sample bremsstrahlung interaction parameters
inline bool brems_XS::interaction(real energy, ranmar &rndm,
                                  real &energy_n, real &cos_t_n, real &sin_t_n,
                                  real &energy_x, real &cos_t_x, real &sin_t_x)
{
   if (energy <= TC)
   {
      energy_n = energy;
      cos_t_n  = ONE;
      sin_t_n  = ZERO;
      energy_x = ZERO;
      cos_t_x  = ONE;
      sin_t_x  = ZERO;
      return(false);
   }

#ifdef BREMS_EXACT
   real e2         = energy*energy;     // energy squared
#endif
   real rej_weight = ONE;               // rejection weight
   do
   {
      energy_x   = TC*exp(rndm.number()*log(energy/TC));
#ifdef BREMS_EXACT
      rej_weight = ONE - energy_x/energy + 0.75*energy_x*energy_x/e2;
#else
      rej_weight = ONE - energy_x/energy;
#endif
   }
   while (rndm.number() > rej_weight);

   // new energy of the primary electron
   energy_n = energy - energy_x;

   // the direction of the electron is unchanged (as in EGS4)
   cos_t_n  = ONE;
   sin_t_n  = ZERO;

   // fixed photon angle (as in EGS4)
   sin_t_x = sin(EMASS/energy);
   cos_t_x = sqrt(max_of(ZERO,(ONE-sin_t_x)*(ONE+sin_t_x)));

   return(true);
}

/*****************************************************************************
 * class compton_XS_diff:                                                    *
 *    differential compton cross section                                     *
 *****************************************************************************/

class compton_XS_diff
{
   private:
      int      ne_in;        // number of bins for initial photon energy
      int      n_bin;        // number of bins for cross section data
      real     energy_min;   // minimum energy
      real     energy_max;   // maximum energy
      real     inv_delta;    // inverse bin size for photon energy (log scale)
      real     distance;     // distance between first and last data bin
      real    *x_data;       // pointer to cross section data
      real   **x_matrix;     // pointer to cross section matrix (direct access)
   public:
      // default constructor (no Compton data)
      compton_XS_diff(void) {
         ne_in = 0; n_bin = 0; x_data = NULL; x_matrix = NULL; }

      // construct differential Compton cross section by reading the data file
      compton_XS_diff(char *file_name)
      {
         ne_in = 0; n_bin = 0; x_data = NULL; x_matrix = NULL;
         read(file_name);
      }

      // input Compton data from file
      void read(char *);

      // free memory
      ~compton_XS_diff(void) {
         delete [] x_matrix; delete [] x_data; }

      // sample Compton interaction parameters
      inline void interaction(real, real, ranmar &,
                              real &, real &, real &,
                              real &, real &, real &);
};

// sample Compton interaction parameters
inline void compton_XS_diff::interaction(real energy, real rnno, ranmar &rndm,
                                 real &energy_x, real &cos_t_x, real &sin_t_x,
                                 real &energy_e, real &cos_t_e, real &sin_t_e)
{
   int   ie_in,i_bin;    // array indices
   real  fe_in,f_bin;    // auxiliary variables to determine array indices

   // determine initial energy array index
   if (energy <= energy_min)
   {
      ie_in = 0;
   }
   else
   {
      if (energy < energy_max)
      {
         fe_in  = log(energy/energy_min)*inv_delta;
         ie_in  = int(fe_in);
         fe_in -= double(ie_in);
         if (rndm.number() < fe_in) ++ie_in;
      }
      else
      {
         ie_in = ne_in - 1;
      }
   }

   // determine array index for the scattered photon energy
   f_bin  = distance*rnno;
   i_bin  = int(f_bin);
   f_bin -= double(i_bin);

   // the minimum energy of the scattered photon
   real k_min = energy/(ONE + TWO*energy/EMASS);

   // determine Compton photon contribution (interpolate)
   real x = (ONE-f_bin)*x_matrix[ie_in][i_bin]
                + f_bin*x_matrix[ie_in][i_bin+1];

   // energy of the scattered photon
   energy_x = k_min*exp(x*log(energy/k_min));

   // scattering angle
   cos_t_x = ONE + EMASS/energy - EMASS/energy_x;
   sin_t_x = sqrt(max_of(ZERO,(ONE-cos_t_x)*(ONE+cos_t_x)));

   // electron energy
   energy_e = energy - energy_x;

   // electron scattering angle
   real p1 = energy_e/energy;
   real p2 = energy_e*(energy_e+TWOxEMASS);
   if (p2 <= ZERO)
   {
      cos_t_e = ZERO;
      sin_t_e = ONE;
   }
   else
   {
      cos_t_e = (energy+EMASS)*p1/sqrt(p2);
      sin_t_e = sqrt(max_of(ZERO,(ONE-cos_t_e)*(ONE+cos_t_e)));
   }

   return;
}

/*****************************************************************************
 * class compton_XS_total:                                                   *
 *    total compton cross section                                            *
 *****************************************************************************/

class compton_XS_total
{
   private:
      real     factorc;     // material dependent factor (n_e x pi x r_0 x r_0)

   public:
      // initialize the material dependent cross section factor
      //   edens: electron density in units of 10^23 cm^-3
      void init(real edens) { factorc = edens * ONE_PIxR0xR0; }

      // construct the Compton cross section using a given electron density
      //    default: water
      compton_XS_total(real edens=3.343) { init(edens); }

      // construct the Compton cross section by reading the electron density
      // from the cross section file
      // attention: we do not use the Compton data from this file!
      compton_XS_total(char *);

      // get total Compton cross section for the specified energy
      real get(real);
};

/*****************************************************************************
 * class pair_XS_diff:                                                       *
 *    differential pair cross section                                        *
 *****************************************************************************/

class pair_XS_diff
{
   private:
      int      ne_in;        // number of bins for initial photon energy
      int      n_bin;        // number of bins for cross section data
      real     energy_min;   // minimum energy
      real     energy_max;   // maximum energy
      real     delta;        // bin size for initial photon energy (log scale)
      real     distance;     // distance between first and last data bin
      real    *x_data;       // pointer to cross section data
      real    *f_data;       // pointer to cross section data
      real   **x_matrix;     // pointer to cross section matrix (direct access)
      real   **f_matrix;     // pointer to cross section matrix (direct access)
   public:
      // default constructor (no pair data)
      pair_XS_diff(void)
      {
         ne_in = 0; n_bin = 0;
         x_data = NULL; f_data = NULL; x_matrix = NULL; f_matrix = NULL;
      }

      // construct differential pair cross section by reading the data file
      pair_XS_diff(char *file_name)
      {
         x_data = NULL; f_data = NULL; x_matrix = NULL; f_matrix = NULL;
         read(file_name);
      }

      // input pair data from file
      void read(char *);

      // free memory
      ~pair_XS_diff(void)
      {
         delete [] x_matrix; delete [] x_data;
         delete [] f_matrix; delete [] f_data;
      }

      // sample pair production parameters
      inline void interaction(real, real, ranmar &,
                              real &, real &, real &,
                              real &, real &, real &);
};

// sample pair production parameters
inline void pair_XS_diff::interaction(real energy, real rnno, ranmar &rndm,
                              real &energy_e, real &cos_t_e, real &sin_t_e,
                              real &energy_p, real &cos_t_p, real &sin_t_p)
{
   int   ie_in,i_bin;    // array indices
   real  fe_in,f_bin;    // auxiliary variables to determine array indices
   real  x;              // positron contribution
   real  tmp,dx;         // auxiliary variables
   const real EPSI=0.01;

   // sum of electron and positron energy
   real energy_sum = energy - TWOxEMASS;
   if (energy_sum < ZERO)
   {
      xvmc_error("pair_XS_diff::interaction",
                 "photon energy < 2 electron mass",8);
   }

#ifdef PAIR_ENERGY_EGS4
   const real BPAR   = 0.853769;
   const real ELIMIT = 2.1;
   real  e_small,xi,eta0,eta1,eta2,eta3;
   if (energy < ELIMIT)
   {
      e_small = EMASS;
   }
   else
   {
      do
      {
         eta0 = rndm.number();
         if (eta0 < BPAR)
         {
            xi = rndm.number();
         }
         else
         {
            eta1 = rndm.number();
            eta2 = rndm.number();
            eta3 = rndm.number();
            xi   = ONE-max_of(eta1,max_of(eta2,eta3));
         }
         e_small = ONE_HALF*xi*energy;
      }
      while (e_small < EMASS);
   }
   if (rnno < ONE_HALF)
   {
      energy_e = e_small-EMASS;
      energy_p = energy-e_small-EMASS;
   }
   else
   {
      energy_p = e_small-EMASS;
      energy_e = energy-e_small-EMASS;
   }
#else
   // determine array index
   f_bin  = distance*rnno;
   i_bin  = int(f_bin);
   f_bin -= double(i_bin);

   // determine positron contribution x
   if ((energy <= energy_min) || (energy >= energy_max))
   {
      ie_in = 0;
      if (energy >= energy_max) ie_in = ne_in - 1;
      tmp = f_matrix[ie_in][i_bin+1]/f_matrix[ie_in][i_bin] - ONE;
      if (fabs(tmp) < EPSI)
      {
         x = (ONE-f_bin)*x_matrix[ie_in][i_bin]
                 + f_bin*x_matrix[ie_in][i_bin+1];
      }
      else
      {
         dx = x_matrix[ie_in][i_bin+1] - x_matrix[ie_in][i_bin];
         x  = x_matrix[ie_in][i_bin]
            + (sqrt(ONE+f_bin*tmp*(tmp+TWO))-ONE)*dx/tmp;
      }
   }
   else
   {
      fe_in = log(energy/energy_min)/delta;
      ie_in  = int(fe_in);
      fe_in -= double(ie_in);
      if (rndm.number() < fe_in) ++ie_in;
      tmp = f_matrix[ie_in][i_bin+1]/f_matrix[ie_in][i_bin] - ONE;
      if (fabs(tmp) < EPSI)
      {
         x = (ONE-f_bin)*x_matrix[ie_in][i_bin]
                 + f_bin*x_matrix[ie_in][i_bin+1];
      }
      else
      {
         dx = x_matrix[ie_in][i_bin+1] - x_matrix[ie_in][i_bin];
         x  = x_matrix[ie_in][i_bin]
            + (sqrt(ONE+f_bin*tmp*(tmp+TWO))-ONE)*dx/tmp;
      }
   }

   // output energies
   energy_p = x*energy_sum;
   energy_e = (ONE-x)*energy_sum;
#endif

   // sample pair angles
#ifdef PAIR_ANGLE_EGSNRC
   // leading order angular distribution (see EGSnrc manual)
   real p2,beta,eta; // momentum squared, velocity, random number

   // polar electron angle
   p2      = energy_e*(energy_e+TWOxEMASS);
   beta    = sqrt(p2/(p2+EMASSxEMASS));
   eta     = TWO*rndm.number()-ONE;
   cos_t_e = (beta+eta)/(beta*eta+ONE);
   sin_t_e = sqrt(max_of(ZERO,(ONE-cos_t_e)*(ONE+cos_t_e)));

   // polar positron angle
   p2      = energy_p*(energy_p+TWOxEMASS);
   beta    = sqrt(p2/(p2+EMASSxEMASS));
   eta     = TWO*rndm.number()-ONE;
   cos_t_p = (beta+eta)/(beta*eta+ONE);
   sin_t_p = sqrt(max_of(ZERO,(ONE-cos_t_p)*(ONE+cos_t_p)));
#else
   // simple estimate, fixed angle as in EGS4
   sin_t_e = sin(EMASS/energy);
   cos_t_e = sqrt(max_of(ZERO,(ONE-sin_t_e)*(ONE+sin_t_e)));
   sin_t_p = sin_t_e;
   cos_t_p = cos_t_e;
#endif

   return;
}

/*****************************************************************************
 * class pair_XS_total:                                                      *
 *    total pair cross section data                                          *
 *****************************************************************************/

class pair_XS_total
{
   private:
      int      n_energy;       // number of bins for photon energy
      real     energy_min;     // minimum energy
      real     energy_max;     // maximum energy
      real     log_energy_min; // ln(minimum energy)
      real     log_energy_max; // ln(maximum energy)
      real     inv_delta;      // inverse bin size for photon energy (log scale)
      real    *sigma;          // pointer to total cross section data

   public:
      // default constructor (no total pair cross section data)
      pair_XS_total(void) { n_energy = 0; sigma = NULL; }

      // construct pair cross section by reading the NIST file
      pair_XS_total(char *file_name) { sigma = NULL; nist(file_name); }

      // delete total pair data
      ~pair_XS_total(void);

      // input total pair data from NIST (XCOM) file
      void nist(char *);

      // input total pair data from file
      void read(char *);

      // get total pair cross section for the specified energy
      real get(real);
};

/*****************************************************************************
 * class photo_XS_total:                                                     *
 *    total photoelectric absorption cross section                           *
 *****************************************************************************/

class photo_XS_total
{
   private:
      int      n_energy;       // number of bins for photon energy
      real     energy_min;     // minimum energy
      real     energy_max;     // maximum energy
      real     log_energy_min; // ln(minimum energy)
      real     log_energy_max; // ln(maximum energy)
      real     inv_delta;      // inverse bin size for photon energy (log scale)
      real    *log_sigma;      // pointer to total cross section (log scale)

   public:
      // default constructor (no total photo cross section data)
      photo_XS_total(void) { n_energy = 0; log_sigma = NULL; }

      // construct photo cross section by reading the NIST file
      photo_XS_total(char *file_name) { log_sigma = NULL; nist(file_name); }

      // delete total photo data
      ~photo_XS_total(void);

      // input total photo data from NIST (XCOM) file
      void nist(char *);

      // get total photo cross section for the specified energy
      real get(real);
};

#endif /* _INTERACTIONS_H_ */
