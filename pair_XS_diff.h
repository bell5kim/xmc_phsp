#ifndef _PAIR_XS_DIFF_H_
#define _PAIR_XS_DIFF_H_

/*****************************************************************************
 * pair_XS_diff.h:                                                           *
 *    class declarations and inline member functions for:                    *
 *     class pair_XS_diff: differential pair cross section                   *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 14.12.1999      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************
#include <string>
#include "definitions.h"
#include "ranmar.h"

/*****************************************************************************
 * class pair_XS_diff:                                                       *
 *    differential pair cross section                                        *
 *****************************************************************************/

class PairDiffData
{
  friend PairDiffData* PairDiffFactory(const std::string& file_name);
  friend class pair_XS_diff;     // TODO goes away with inheritance

 public:

  PairDiffData(int nEin, int nBin, real energyMin, real energyMax);
  virtual ~PairDiffData();

 protected:

  int      ne_in;        // number of bins for initial photon energy
  int      n_bin;        // number of bins for cross section data
  real     energy_min;   // minimum energy
  real     energy_max;   // maximum energy
  real*    x_data;       // pointer to cross section data
  real*    f_data;       // pointer to cross section data
  real**   x_matrix;     // pointer to cross section matrix (direct access)
  real**   f_matrix;     // pointer to cross section matrix (direct access)

 private:

  bool invariant() {return true;}                // TODO

  PairDiffData();                                // not implemented
  PairDiffData(const PairDiffData&);             // not implemented
  PairDiffData& operator=(const PairDiffData&);  // not implemented
};


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

  pair_XS_diff(const PairDiffData& );
  ~pair_XS_diff();

  // sample pair production parameters
  inline void interaction(real energy,
                          real rnno,
                          ranmar &rndm,
                          real &energy_e,
                          real &cos_t_e,
                          real &sin_t_e,
                          real &energy_p,
                          real &cos_t_p,
                          real &sin_t_p) const;
};

// sample pair production parameters
inline void pair_XS_diff::interaction(real energy, real rnno, ranmar &rndm,
				      real &energy_e, real &cos_t_e, real &sin_t_e,
				      real &energy_p, real &cos_t_p, real &sin_t_p) const
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

#endif /* _PAIR_XS_DIFF_H_ */
