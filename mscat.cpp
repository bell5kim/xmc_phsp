/*****************************************************************************
 * mscat.cpp:                                                                *
 *    function:                                                              *
 *       mscat:          sample multiple scattering angle                    *
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

#include "definitions.h"
#include "global.h"
#include "ranmar.h"

/*****************************************************************************
 * function mscat: sample multiple scattering angle                          *
 *    input:  omega_0: number of elastic scatterings (0 < omega_0 < 100000)  *
 *            xr:      cos_t = 1 - xr*y^2                                    *
 *    output: cos_t, sin_t - cos(theta) and sin(theta)                       *
 *    return: true if there is an error, false without any error             *
 *****************************************************************************/

bool mscat(const real &omega_0, const real &xr, real &cos_t, real &sin_t,
           ranmar &rndm)
{
   const int max_reject = 100;

   const real FIVE    =  5.0;
   const real FIFTEEN = 15.0;

   const real TINY = 1.0e-35;

   const real PAR1 = -0.072027;
   const real PAR2 =  1.674324;
   const real PAR3 = -0.0897635;
   const real PAR4 =  0.0065706;
   const real PAR5 = -0.000190272;

   const real W2_D0[] = {  1.00000,  0.91579,  0.84469,  0.77406,  0.70108,
                           0.62883,  0.54964,  0.47394,  0.40317,  0.34305,
                           0.29720,  0.24660,  0.20917,  0.17809,  0.15481 };
   const real W2_D1[] = { -0.08421, -0.07110, -0.07063, -0.07298, -0.07225,
                          -0.07919, -0.07570, -0.07077, -0.06012, -0.04585,
                          -0.05060, -0.03743, -0.03108, -0.02328, -0.02328 };
   const real W2_G0[] = {  0.00000,  0.05851,  0.12931,  0.19984,  0.27079,
                           0.34009,  0.41786,  0.49502,  0.57008,  0.63527,
                           0.68549,  0.74016,  0.78039,  0.81347,  0.83812 };
   const real W2_G1[] = {  0.05851,  0.07080,  0.07053,  0.07095,  0.06930,
                           0.07777,  0.07716,  0.07506,  0.06519,  0.05022,
                           0.05467,  0.04023,  0.03308,  0.02465,  0.02465 };

   int     n_reject,i_omega;
   real    w_2,w_3,eta,blc,b,xx,amu,amu2,d_omega,s_omega,alpha;

   n_reject = 0;
   if (omega_0 > FIFTEEN)
   {
   // start: omega_0 > 15
      blc = log(omega_0);
      b=PAR1+blc*(PAR2+blc*(PAR3+blc*(PAR4+blc*PAR5)));
      w_2 = ONE/(0.9499+10.75199*blc);
      do
      {
         ++n_reject;
         if (n_reject > max_reject)
         {
            xvmc_warning("mscat","too many rejections",1);
            xvmc_warning("omega_0",omega_0,0);
            xvmc_warning("xr",xr,0);
            cos_t = TWO*rndm.number() - ONE;
            sin_t = max_of((ONE-cos_t)*(ONE+cos_t),ZERO);
            sin_t = sqrt(sin_t);
            return(true);
         }
         xx = rndm.number();
         if (xx < w_2)
         {
            amu2  = 8.0*b;
            eta   = max_of(rndm.number(),TINY);
            cos_t = ONE - xr*amu2*(ONE-eta)/eta;
         }
         else
         {
            w_3 = w_2 + ONE/(0.3341+0.26164*blc);
            if (xx < w_3)
            {
               amu2  = 1.5944*b;
               eta   = max_of(rndm.number(),rndm.number());
               cos_t = ONE - xr*amu2*(ONE/eta - ONE);
            }
            else
            {
               eta   = max_of(rndm.number(),TINY);
               cos_t = ONE + xr*b*log(eta);
            }
         }
         sin_t=(ONE-cos_t)*(ONE+cos_t);
      }
      while (sin_t < ZERO);
   // end: omega_0 > 15
   }
   else
   {
   // start: omega_0 <= 15
      if (omega_0 < FIVE)
      {
         w_2 = exp(-omega_0);
         if (rndm.number() < w_2)
         {
            cos_t = ONE;
            sin_t = ZERO;
            return(false);
         }
      }
      do
      {
         ++n_reject;
         if (n_reject > max_reject)
         {
            xvmc_warning("mscat","too many rejections",1);
            xvmc_warning("omega_0",omega_0,0);
            xvmc_warning("xr",xr,0);
            cos_t = TWO*rndm.number() - ONE;
            sin_t = max_of((ONE-cos_t)*(ONE+cos_t),ZERO);
            sin_t = sqrt(sin_t);
            return(true);
         }
         i_omega = int(omega_0);
         d_omega = omega_0 - float(i_omega);
         w_2 = W2_D0[i_omega] + d_omega*W2_D1[i_omega];
         s_omega = sqrt(omega_0);
         xx = rndm.number();
         if (xx < w_2)
         {
            amu   = ONE/s_omega + 0.301*s_omega + 0.0055*omega_0*omega_0;
            amu2  = amu*amu;
            eta   = max_of(rndm.number(),TINY);
            cos_t = ONE - xr*amu2*(ONE-eta)/eta;
         }
         else
         {
            w_3 = w_2 + W2_G0[i_omega] + d_omega*W2_G1[i_omega];
            if (xx < w_3)
            {
               amu   = 3.15/s_omega + 0.8 + 0.05*omega_0;
               amu2  = amu*amu;
               eta   = max_of(rndm.number(),rndm.number());
               cos_t = ONE - xr*amu2*(ONE/eta - ONE);
            }
            else
            {
               alpha = (ONE+0.301*omega_0)/s_omega + 0.0055*omega_0*omega_0;
               alpha = alpha*alpha/8.0;
               eta   = max_of(rndm.number(),TINY);
               cos_t = ONE + xr*alpha*log(eta);
            }
         }
         sin_t=(ONE-cos_t)*(ONE+cos_t);
      }
      while (sin_t < ZERO);
   // end: omega_0 <= 15
   }
   sin_t = sqrt(sin_t);
   return(false);
}
