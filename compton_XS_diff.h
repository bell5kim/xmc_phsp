#ifndef _COMPTON_XS_DIFF_H_
#define _COMPTON_XS_DIFF_H_

/*****************************************************************************
 * compton_XS_diff.h:                                                        *
 *    class declarations and inline member functions for:                    *
 *     class compton_XS_diff: differential Compton cross section             *
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
 * class compton_XS_diff:                                                    *
 *    differential compton cross section                                     *
 *****************************************************************************/

class ComptonDiffData
{
  friend ComptonDiffData* ComptonFactory(const std::string& file_name);
  friend class compton_XS_diff;     // TODO goes away with inheritance

 public:

  ComptonDiffData(int nEin, int nBin, real energyMin, real energyMax);
  virtual ~ComptonDiffData();

 protected:

  int      ne_in;        // number of bins for initial photon energy
  int      n_bin;        // number of bins for cross section data
  real     energy_min;   // minimum energy
  real     energy_max;   // maximum energy
  real*    x_data;       // pointer to cross section data
  real**   x_matrix;     // pointer to cross section matrix (direct access)

 private:

  bool invariant() {return true;}              // TODO

  ComptonDiffData();                               // not implemented
  ComptonDiffData(const ComptonDiffData&);             // not implemented
  ComptonDiffData& operator=(const ComptonDiffData&);  // not implemented
};


class compton_XS_diff // TODO better : private ComptonDiffData, conflicts with other constructors below
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


 public:

  compton_XS_diff(const ComptonDiffData& );
  ~compton_XS_diff(void);

  // sample Compton interaction parameters
  inline void interaction(real energy,
                          real rnno,
                          ranmar &rndm,
                          real &energy_x,
                          real &cos_t_x,
                          real &sin_t_x,
                          real &energy_e,
                          real &cos_t_e,
                          real &sin_t_e) const;
};

// sample Compton interaction parameters
inline void compton_XS_diff::interaction(real energy, real rnno, ranmar &rndm,
                                 real &energy_x, real &cos_t_x, real &sin_t_x,
                                 real &energy_e, real &cos_t_e, real &sin_t_e) const
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

#endif /* _COMPTON_XS_DIFF_H_ */
