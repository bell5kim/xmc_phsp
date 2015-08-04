/*****************************************************************************
 * xvmc_util.cpp: utility functions                                          *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 99/12/13        *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <stdlib.h>

#include "definitions.h"
#include "global.h"

// find upper and lower array indices
void get_index(int n, real *x_array, real x, int &lower, int &upper)
{
   int m;

   if ( x_array[0] < x_array[n-1] )
   {
      // x_array[0] <= x_array[1] <= ... <= x_array[n-1]
      if ( x <= x_array[0] )
      {
         lower = 0;
         upper = 1;
         return;
      }
      if ( x >= x_array[n-1] )
      {
         lower = n-2;
         upper = n-1;
         return;
      }
      // x_array[0] < x < x_array[n-1]
      lower = 0;
      upper = n-1;
      while ( upper > lower+1 )
      {
         m = (lower+upper)/2;
         if (x < x_array[m]) upper = m;
         else                lower = m;
      }
      lower = upper - 1;
      return;
   }
   else
   {
      // x_array[n-1] <= ... <= x_array[1] <= x_array[0]
      if ( x <= x_array[n-1] )
      {
         lower = n-2;
         upper = n-1;
         return;
      }
      if ( x >= x_array[0] )
      {
         lower = 0;
         upper = 1;
         return;
      }
      // x_array[n-1] < x < x_array[0]
      lower = n-1;
      upper = 0;
      while ( upper < lower-1 )
      {
         m = (lower+upper)/2;
         if (x < x_array[m]) upper = m;
         else                lower = m;
      }
      upper = lower;
      lower = upper - 1;
      return;
   }
}

#ifndef INTERPOLATE_INLINE
// find function value f(x) by interpolation
// (extrapolation if x < x_array[0] or x > x_array[n-1])
real interpolate(int n, real *x_array, real *f_array, real x)
{
   int  lower,upper;    // array indices
   real delta,p,q;

   get_index(n,x_array,x,lower,upper);

   delta = x_array[upper]-x_array[lower];
   if (delta != ZERO) p = (x_array[upper]-x)/delta;
   else               p = ONE;
   q = ONE-p;
   return(p*f_array[lower]+q*f_array[upper]);
}
#endif   // #ifndef INTERPOLATE_INLINE

#ifndef ROTATE_INLINE
// 3D rotation of vector "dir" by two angles "theta" and "phi"
void rotate(real theta, real phi, real_3 &dir)
{
   const real ONE     =  1.0;
   const real EPSILON =  1.0e-10;
   real temp,temp_p,temp_x,temp_y,temp_x1,temp_y1;

   real cos_t = cos(theta);
   real sin_t = sin(theta);
   real cos_p = cos(phi);
   real sin_p = sin(phi);

   real sin_z = (ONE-dir.z)*(ONE+dir.z);
   if (sin_z > EPSILON)
   {
      sin_z   = sqrt(sin_z);
      temp    = sin_t/sin_z;
      temp_p  = dir.z * cos_p;
      temp_x  = dir.x * cos_t;
      temp_y  = dir.y * cos_t;
      temp_x1 = temp_p*dir.x - dir.y*sin_p;
      temp_y1 = temp_p*dir.y + dir.x*sin_p;
      dir.x   = temp*temp_x1 + temp_x;
      dir.y   = temp*temp_y1 + temp_y;
      dir.z   = dir.z*cos_t - sin_z*sin_t*cos_p;
   }
   else
   {
      dir.x = sin_t*cos_p;
      dir.y = sin_t*sin_p;
      dir.z = cos_t;
   }
}

// 3D rotation of vector "dir" by two angles "theta" and "phi"
void rotate(real cos_t, real sin_t, real cos_fi, real sin_fi,
            real_3 &dir)
{
   const real ONE     =  1.0;
   const real EPSILON =  1.0e-10;
   real sin_z,temp,temp_fi,temp_x,temp_y,temp_x1,temp_y1;

   sin_z = (ONE-dir.z)*(ONE+dir.z);
   if (sin_z > EPSILON)
   {
      sin_z   = sqrt(sin_z);
      temp    = sin_t/sin_z;
      temp_fi = dir.z * cos_fi;
      temp_x  = dir.x * cos_t;
      temp_y  = dir.y * cos_t;
      temp_x1 = temp_fi*dir.x - dir.y*sin_fi;
      temp_y1 = temp_fi*dir.y + dir.x*sin_fi;
      dir.x   = temp*temp_x1 + temp_x;
      dir.y   = temp*temp_y1 + temp_y;
      dir.z   = dir.z*cos_t - sin_z*sin_t*cos_fi;
   }
   else
   {
      dir.x=sin_t*cos_fi;
      dir.y=sin_t*sin_fi;
      dir.z=cos_t;
   }
}

// 3D rotation of vector "(cos_x,cos_y,cos_z)" by two angles "theta" and "phi"
// single precision (float)
void rotate(float cos_t, float sin_t, float cos_fi, float sin_fi,
            float &cos_x, float &cos_y, float &cos_z)
{
   const float ONE     =  1.0;
   const float EPSILON =  1.0e-10;
   float sin_z,temp,temp_fi,temp_x,temp_y,temp_x1,temp_y1;

   sin_z = (ONE-cos_z)*(ONE+cos_z);
   if (sin_z > EPSILON)
   {
      sin_z   = sqrt(sin_z);
      temp    = sin_t/sin_z;
      temp_fi = cos_z * cos_fi;
      temp_x  = cos_x * cos_t;
      temp_y  = cos_y * cos_t;
      temp_x1 = temp_fi*cos_x - cos_y*sin_fi;
      temp_y1 = temp_fi*cos_y + cos_x*sin_fi;
      cos_x   = temp*temp_x1 + temp_x;
      cos_y   = temp*temp_y1 + temp_y;
      cos_z   = cos_z*cos_t - sin_z*sin_t*cos_fi;
   }
   else
   {
      cos_x=sin_t*cos_fi;
      cos_y=sin_t*sin_fi;
      cos_z=cos_t;
   }
}

// 3D rotation of vector "(cos_x,cos_y,cos_z)" by two angles "theta" and "phi"
// double precision (double)
void rotate(double cos_t, double sin_t, double cos_fi, double sin_fi,
            double &cos_x, double &cos_y, double &cos_z)
{
   const double ONE     =  1.0;
   const double EPSILON =  1.0e-10;
   double sin_z,temp,temp_fi,temp_x,temp_y,temp_x1,temp_y1;

   sin_z = (ONE-cos_z)*(ONE+cos_z);
   if (sin_z > EPSILON)
   {
      sin_z   = sqrt(sin_z);
      temp    = sin_t/sin_z;
      temp_fi = cos_z * cos_fi;
      temp_x  = cos_x * cos_t;
      temp_y  = cos_y * cos_t;
      temp_x1 = temp_fi*cos_x - cos_y*sin_fi;
      temp_y1 = temp_fi*cos_y + cos_x*sin_fi;
      cos_x   = temp*temp_x1 + temp_x;
      cos_y   = temp*temp_y1 + temp_y;
      cos_z   = cos_z*cos_t - sin_z*sin_t*cos_fi;
   }
   else
   {
      cos_x=sin_t*cos_fi;
      cos_y=sin_t*sin_fi;
      cos_z=cos_t;
   }
}
#endif   // #ifndef ROTATE_INLINE

// swap byte order
void swap_bytes( float &f )
{
   unsigned int  *i = (unsigned int *) &f;
   unsigned char c1,c2,c3,c4;

   c4 = (unsigned char)   *i;
   c3 = (unsigned char) ( *i >>  8 );
   c2 = (unsigned char) ( *i >> 16 );
   c1 = (unsigned char) ( *i >> 24 );
   *i =   ((unsigned int) c1)         +
        ( ((unsigned int) c2) <<  8 ) +
        ( ((unsigned int) c3) << 16 ) +
        ( ((unsigned int) c4) << 24 );
   return;
}
