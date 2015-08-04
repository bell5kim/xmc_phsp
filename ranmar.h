#ifndef _RANMAR_H_
#define _RANMAR_H_

/*****************************************************************************
 * ranmar.h:                                                                 *
 *    class declarations and inline member functions for:                    *
 *       class ranmar: uncorrelated random sequences                         *
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

#endif /* _RANMAR_H_ */
