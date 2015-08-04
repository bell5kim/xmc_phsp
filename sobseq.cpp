/*****************************************************************************
 * sobseq.cpp:                                                               *
 *    class member functions for:                                            *
 *       sobseq: quasi-random sequences                                      *
 *****************************************************************************/

#include <iostream>
using namespace std;

#include "definitions.h"
#include "global.h"
#include "sobseq.h"

/*****************************************************************************
 * class sobseq: generates quasi-random sequences                            *
 *               (adapted from "Numerical Recipes")                          *
 *                                                                           *
 * member functions:                                                         *
 *    init:    initialize sobol sequence                                     *
 *    sobseq:  constructors                                                  *
 *    re_init: re-initialize sobol sequence                                  *
 *    next:    generates next set of sobol numbers                           *
 *    number:  returns a sobol number from present set                       *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 99/12/08        *
 *                                                                           *
 *****************************************************************************/

void sobseq::init(const unsigned int n)
{
   const unsigned int MAXDEG = 5;
   const unsigned int ip[MAXDIM]   = {0,1,1,2,1,4,2,4,7,11,13,14};
   const unsigned int mdeg[MAXDIM] = {1,2,3,3,4,4,5,5,5,5,5,5};
   const unsigned int ivs[MAXDIM*MAXDEG] = {
              1,  1,  1,  1,  1,  1,       1,  1,  1,  1,  1,  1,
              3,  1,  3,  3,  1,  1,       3,  3,  1,  1,  3,  3,
              5,  7,  7,  3,  3,  5,       7,  5,  3,  3,  5,  7,
             15, 11,  5, 15, 13,  9,       9,  7,  5,  7,  9,  9,
             17, 39,  7, 49,  0,  0,       5, 11, 13, 17, 19, 31};
   int           kk;
   unsigned long i, ipp, *iu[MAXBIT];

   for (unsigned int k=0; k < MAXDIM*MAXBIT; ++k)
   {
      if ( k < MAXDIM*MAXDEG ) iv[k] = ivs[k];
      else                     iv[k] = 0;
   }

//   for (int j=0, int k=0; j < MAXBIT; ++j, k += MAXDIM)
   kk=0;
   for (unsigned int j=0; j < MAXBIT; ++j)
   {
      iu[j] = &iv[kk];
      kk += MAXDIM;
   }

   for (unsigned int k=0; k < MAXDIM; ++k)
   {
      ix[k] = 0;
      for (unsigned int j=0; j < mdeg[k]; ++j) iu[j][k] <<= (MAXBIT-1-j);
      for (int j=mdeg[k]; j < MAXBIT; ++j)
      {
         ipp = ip[k];
         i   = iu[j-mdeg[k]][k];
         i ^= (i >> mdeg[k]);
         for (int l = mdeg[k]-1; l >= 1; --l)
         {
            if (ipp & 1) i ^= iu[j-l][k];
            ipp >>= 1;
         }
         iu[j][k]=i;
      }
   }
   fac=1.0/(1L << MAXBIT);
   in =0;
   if (n > MAXDIM)
   {
      dim = MAXDIM;
      cerr << " WARNING: maximum dimension for sobol sequence is "
           << int(MAXDIM) << "!" << endl;
   }
   else if (n < 1)
   {
      dim = 1;
      cerr << " WARNING: minimum dimension for sobol sequence is 1!" << endl;
   }
   else dim = n;
}

sobseq::sobseq(void)
{
   init(MAXDIM);
}

sobseq::sobseq(const unsigned int n)
{
   init(n);
}

void sobseq::re_init(void)
{
   init(dim);
}

void sobseq::re_init(const unsigned int n)
{
   init(n);
}

#ifndef RANDOM_INLINE
void sobseq::next(void)
{
   unsigned long im = in;
   int           j;

   for (j=0; j < MAXBIT; ++j)
   {
      if (!(im & 1)) break;
      im >>= 1;
   }
   if (j >= MAXBIT) xvmc_error("function sobseq::next","MAXBIT too small",8);
   im =j*MAXDIM;
   for (unsigned int k=0; k < dim; ++k)
   {
      ix[k] ^= iv[im+k];
      x[k]   = ix[k]*fac;
   }
   ++in;
}

real sobseq::number(const unsigned int i)
{
   if (i < dim) return(x[i]);
   else xvmc_error("function sobseq::number", "dimension exceeded",8);
   return(0.0);
}
#endif
