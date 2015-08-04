#ifndef _SOBSEQ_H_
#define _SOBSEQ_H_

/*****************************************************************************
 * sobseq.h:                                                                 *
 *     class declarations and inline member functions for:                   *
 *       class sobseq: quasi-random sequences                                *
 *****************************************************************************/

#include "definitions.h"
#include "global.h"

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
 *****************************************************************************/

class sobseq
{
   private:
      enum           {MAXBIT = 30, MAXDIM = 12};
      unsigned long  in, ix[MAXDIM], iv[MAXDIM*MAXBIT];
      real           fac, x[MAXDIM];
      unsigned int   dim;
      void           init(const unsigned int);
   public:
      sobseq(void);
      sobseq(const unsigned int);
      void re_init(void);
      void re_init(const unsigned int);
      void next(void);
      real number(const unsigned int);
};

#ifdef RANDOM_INLINE
inline void sobseq::next(void)
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

inline real sobseq::number(const unsigned int i)
{
   if (i < dim) return(x[i]);
   else xvmc_error("function sobseq::number", "dimension exceeded",8);
   return(0.0);
}
#endif /* RANDOM_INLINE */

#endif /* _SOBSEQ_H_ */
