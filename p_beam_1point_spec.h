#ifndef _P_BEAM_1POINT_SPEC_H_
#define _P_BEAM_1POINT_SPEC_H_

/*****************************************************************************
 * p_beam_1point_spec.h:                                                     *
 *    class declarations and inline member functions for:                    *
 *       p_beam_1point_spec: point source photon beam, energy spectrum       *
 *                                                                           *
 * Copyright (C) 2002    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 31.01.2002      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************
 
#include "p_beam_1point_mono.h"

// point source photon beam, energy spectrum
class p_beam_1point_spec : public p_beam_1point_mono
{
   public:
      p_beam_1point_spec(const beam_core *this_beam)
         : p_beam_1point_mono(this_beam)
      {
         find_base_data(this_beam);
         get_base_data();
         close_base_file();
      }

      real average_energy(void) const { return(energy_mean); }
      real maximum_energy(void) const { return(energy_max); }

      ~p_beam_1point_spec(void) { delete [] energy_bin; energy_bin = NULL;
                                  delete [] dprob;  dprob  = NULL;
                                  delete [] cprob;  cprob  = NULL; }

   protected:
      // parameters for energy spectrum
      int       n_bin;         // number of bins for the energy spectrum
      real      energy_min;    // minimum energy for spectrum (MeV)
      real     *energy_bin;    // upper energy bin limit
      real     *dprob;         // differential probability distribution
      real     *cprob;         // cumulative probability distribution
      double    norm_factor;   // normalization factor
      double    energy_mean;   // mean energy (MeV)
      real      energy_max;    // maximum energy for spectrum (MeV)

      // get base data from file, initialize energy spectrum
      void      get_base_data(void);

      // sample energy from spectrum
      real      sample_energy(ranmar &);
};

// sample energy from spectrum
inline real p_beam_1point_spec::sample_energy(ranmar &rndm)
{
   // we need two random numbers
   real rndm1 = rndm.number();
   real rndm2 = rndm.number();

   // select the energy bin using the first random number
   register int i = 0;
   while ((cprob[i]<rndm1) && (i<(n_bin-1))) ++i;

   // energy_bin[i-1] <= energy <= energy_bin[i]
   real energy_low;
   if (i==0) energy_low = energy_min;
   else      energy_low = energy_bin[i-1];

   // sample energy within the bin using the second random number
   real energy = energy_low + (energy_bin[i]-energy_low)*rndm2;
   return(energy);
}

#endif  /* _P_BEAM_1POINT_SPEC_H_ */
