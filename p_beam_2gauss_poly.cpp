/*****************************************************************************
 * p_beam_2gauss_poly.cpp:                                                   *
 *    class member functions for:                                            *
 *       p_beam_2gauss_poly: two Gaussian sources, poly-energetic photons    *
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

#include "p_beam_2gauss_poly.h"

// ***********************************************
// member functions of class p_beam_2gauss_poly
// ***********************************************

// initialize energy spectrum, calculate Ep, <E> and s:
// the energy E will be sampled from the distribution p(E) = E^l*exp(-b*E),
// therefore, the most probable energy (energy_prob) is given by Ep = l/b and
// the average energy (energy_mean) is estimated by <E> = (l+1)/b (this formula
// is exactly true only for Emin=0 (energy_min) and Emax=infinity (energy_max))
void p_beam_2gauss_poly::init_energy(void)
{
   if (energy_min_aux <= energy_max_aux)
   {
      energy_min = energy_min_aux;
      energy_max = energy_max_aux;
      if (b_aux > ZERO)
      {
         l = l_aux;
         b = b_aux;
         s = sqrt(TWO*l+ONE);
         energy_prob = l/b;
         energy_mean = (l+ONE)/b; // this is just an estimate
      }
      else
      {
         xvmc_error("p_beam_2gauss_poly::init_energy",
                    "spectrum b value <= 0",8);
      }
   }
   else
   {
      xvmc_error("p_beam_2gauss_poly::init_energy",
                 "minimum spectrum energy > maximum spectrum energy",8);
   }
}

p_beam_2gauss_poly::p_beam_2gauss_poly(const p_beam_2gauss_poly &A)
  :p_beam_2gauss_mono(A)
{
  //p_beam_2gauss_poly variables
  energy_min = A.energy_min;
  energy_max = A.energy_max;
  energy_prob = A.energy_prob;
  energy_mean = A.energy_mean;
  l = A.l;
  b = A.b;
  s = A.s;
};
