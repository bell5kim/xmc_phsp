/*****************************************************************************
 * ranmar.cpp:                                                               *
 *    class member functions for:                                            *
 *       ranmar: uncorrelated random sequences                               *
 *    function:                                                              *
 *       ranmar_test: test random number generator, class ranmar             *
 *****************************************************************************/

#include <iostream>
using namespace std;

#include "definitions.h"
#include "global.h"
#include "ranmar.h"

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
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 99/12/07        *
 *                                                                           *
 *****************************************************************************/

void ranmar::init(int ma1, int ma2, int ma3, int mb1)
{
   na1 = ma1;
   na2 = ma2;
   na3 = ma3;
   nb1 = mb1;
   i1  = SIZE-1;
   i2  = 32;

   for (int i=0; i<SIZE; ++i)
   {
      real s = 0.0;
      real t = 0.5;
      for (int j=0; j<24; ++j)
      {
         int mat = (((ma1*ma2) % 179)*ma3) % 179;
         ma1 = ma2;
         ma2 = ma3;
         ma3 = mat;
         mb1 = (53*mb1+1) % 169;
         if (((mb1*mat) % 64) >= 32) s += t;
         t = 0.5*t;
      }
      u[i] = s;
   }
   c  =   362436.0/16777216.0;
   cd =  7654321.0/16777216.0;
   cm = 16777213.0/16777216.0;

   return;
}

ranmar::ranmar(void)
{
   init(12,34,56,78);
}

ranmar::ranmar(const int ma1, const int ma2, const int ma3, const int mb1)
{
   init(ma1,ma2,ma3,mb1);
}

void ranmar::re_init(void)
{
   init(na1,na2,na3,nb1);
}

void ranmar::re_init(const int ma1, const int ma2, const int ma3, const int mb1)
{
   init(ma1,ma2,ma3,mb1);
}

#ifndef RANDOM_INLINE
real ranmar::number(void)
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
#endif

ranmar_state ranmar::status(void)
{
   ranmar_state state;
   state.i1 = i1;
   state.i2 = i2;
   return(state);
}

/*****************************************************************************
 * function ranmar_test: test random number generator, class ranmar          *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 99/12/08        *
 *                                                                           *
 *****************************************************************************/

void ranmar_test(const int output)
{
   ranmar  test;
   const real expected[6]={ 6533892.0 , 14220222.0 ,  7275067.0,
                            6172232.0 ,  8354498.0 , 10633180.0};
   real       calculated[6],difference[6];
   real       total_difference=0.0;

   for (int i=0; i<20000; ++i) test.number();

   for (int i=0; i<6; ++i)
   {
      calculated[i] = 4096.0*(4096.0*test.number());
      difference[i] = calculated[i]-expected[i];
      total_difference += difference[i];
   }

   if ((output != 0) || (total_difference != 0.0))
   {
      cout << "\t  === TEST OF THE RANDOM-GENERATOR ===\n"
           << "\t" << "EXPECTED VALUE"
           << "\t" << "CALCULATED VALUE"
           << "\t" << "DIFFERENCE\n";
      for (int i=0; i<6; ++i)
      {
         cout << "\t  "   << expected[i]
              << "\t  "   << calculated[i]
              << "\t\t  " << difference[i] << endl;
      }
      cout << "\t  === END OF TEST ===\n";
   }
}
