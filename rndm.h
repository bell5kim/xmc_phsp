#ifndef _RNDM_H_
#define _RNDM_H_

/*****************************************************************************
 * rndm.h:                                                                   *
 *    class declaration:                                                     *
 *       class ranmar: uncorrelated random sequences                         *
 *       class sobseq: quasi-random sequences                                *
 *    function declaration:                                                  *
 *       ranmar_test: test random number generator, class ranmar             *
 *****************************************************************************/

#include "definitions.h"
#include "global.h"

/*****************************************************************************
 * class ranmar: generates uncorrelated random number sequences              *
 *               (adapted from CERNLIB FORTRAN function RANMAR)              *
 *                                                                           *
 * member functions:                                                         *
 *    init:    initialize random sequence                                    *
 *    ranmar:  constructors                                                  *
 *    re_init: re-initialize random sequence                                 *
 *    number:  returns a random number (inline function if                   *
 *                RANDOM_INLINE is defined                                   *
 *    status:  returns the generator state, 2 integers (type ranmar_state)   *
 *****************************************************************************/

struct ranmar_state
{
   int i1;
   int i2;
};

class ranmar
{
   private:
      enum  {SIZE=97};
      real  u[SIZE],c,cd,cm;
      int   na1,na2,na3,nb1,i1,i2;
      void  init(int, int, int, int);
   public:
      ranmar(void);
      ranmar(const int, const int, const int, const int);
      void re_init(void);
      void re_init(const int, const int, const int, const int);
      real number(void);
      ranmar_state status(void);
};

#ifdef RANDOM_INLINE
inline real ranmar::number(void)
{
   real r = u[i1] - u[i2];
   if (r < 0.0) r += 1.0;
   u[i1] = r;
   --i1;
   if (i1 < 0) i1=SIZE-1;
   --i2;
   if (i2 < 0) i2=SIZE-1;
   c -= cd;
   if (c < 0.0) c += cm;
   r -= c;
   if (r < 0.0) r += 1.0;
   return(r);
}
#endif /* RANDOM_INLINE */

/*****************************************************************************
 * function ranmar_test: test random number generator, class ranmar          *
 *****************************************************************************/

void ranmar_test(const int);

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

#endif /* _RNDM_H_ */
