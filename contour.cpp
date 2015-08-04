/*****************************************************************************
 * contour.cpp:                                                              *
 *    class member functions for:                                            *
 *       contour:    contour defined by a set of (x,y) points                *
 *                                                                           *
 * Copyright (C) 2001    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 14.06.2001      *
 *                                                                           *
 *****************************************************************************/

#include <string.h>
#include "definitions.h"
#include "global.h"
#include "contour.h"

// ****************************************
// member functions of class contour
// ****************************************

// define a contour by num points
contour::contour(const int num)
{
   // the contour is assumed to be open when the contour is created
   closed = false;

   // number of contour points
   num_points = num;

   // By definition the contour must be closed, i.e. the first and the
   // last point must coincide. However, if the input contour isn't closed,
   // we just add one further point. Therefore, at this point the array size
   // is a little bit larger than the number of contour points.
   array_size = num_points + 1;

   // allocate memory
   if ( (x_con = new real[array_size]) == NULL)
   {
      xvmc_error("contour::contour",
                 "cannot allocate memory for x-coordinates",8);
   }
   if ( (y_con = new real[array_size]) == NULL)
   {
      xvmc_error("contour::contour",
                 "cannot allocate memory for y-coordinates",8);
   }

   // initialize contour points
   for (register int i=0; i<array_size; ++i)
   {
      x_con[i] = ZERO;
      y_con[i] = ZERO;
   }

   // initialize minimum and maximum of x and y coordinates
   min_x =  1.0e+10;
   max_x = -1.0e+10;
   min_y =  1.0e+10;
   max_y = -1.0e+10;

   // allocate memory for the intersection point arrays,
   // the array size cannot be larger than the contour point arrays
   if ( (x_inter = new real[array_size]) == NULL )
   {
      xvmc_error("contour::contour",
                 "cannot allocate memory for unsorted intersection points",8);
   }
   if ( (x_min = new real[array_size]) == NULL )
   {
      xvmc_error("contour::contour",
                 "cannot allocate memory for intersection pair point x_min",8);
   }
   if ( (x_max = new real[array_size]) == NULL )
   {
      xvmc_error("contour::contour",
                 "cannot allocate memory for intersection pair point x_max",8);
   }
}

// delete contour
contour::~contour(void)
{
   delete [] x_con;
   delete [] y_con;
   delete [] x_inter;
   delete [] x_min;
   delete [] x_max;
}

// change one contour point, if this is impossible -> return false
bool contour::change_point(int i_new, real x_new, real y_new)
{
   if (i_new < 0)           return(false);
   if (i_new >= num_points) return(false);

   // the contour is assumed to be open when a point is added
   closed = false;

   // new point
   x_con[i_new] = x_new;
   y_con[i_new] = y_new;

   // adjust minimum and maximum of x and y coordinates
   min_x = min_of(min_x,x_new);
   max_x = max_of(max_x,x_new);
   min_y = min_of(min_y,y_new);
   max_y = max_of(max_y,y_new);

   // one contour point is changed
   return(true);
}

// close contour, if this is impossible -> return false
bool contour::close(void)
{
   if ( (x_con[0] == x_con[num_points-1]) &&
        (y_con[0] == y_con[num_points-1])    )
   {
      // the contour is already closed, change closed to true and return true
      closed = true;
      return(true);
   }

   if (num_points >= array_size)
   {
      // cannot close contour, return false
      xvmc_warning("contour::close","cannot close contour",1);
      closed = false;
      return(false);
   }

   // add one contour point to close contour
   ++num_points;
   x_con[num_points-1] = x_con[0];
   y_con[num_points-1] = y_con[0];

   // change closed to true and return true
   closed = true;
   return(true);
}

// get number of contour points
int contour::get_num(void)
{
   if (!closed)
   {
      xvmc_warning("contour::get_num","contour not closed",1);
   }
   return(num_points);
}

// get x-coordinate of point i
real contour::get_x(int i)
{
   if (!closed)
   {
      xvmc_warning("contour::get_x","contour not closed",1);
   }

   if (i < 0)
   {
      xvmc_warning("contour::get_x","i < 0, return zero",1);
      return(ZERO);
   }

   if (i >= num_points)
   {
      xvmc_warning("contour::get_x","i >= num_points, return zero",1);
      return(ZERO);
   }

   return(x_con[i]);
}

// get y-coordinate of point i
real contour::get_y(int i)
{
   if (!closed)
   {
      xvmc_warning("contour::get_y","contour not closed",1);
   }

   if (i < 0)
   {
      xvmc_warning("contour::get_y","i < 0, return zero",1);
      return(ZERO);
   }

   if (i >= num_points)
   {
      xvmc_warning("contour::get_y","i >= num_points, return zero",1);
      return(ZERO);
   }

   return(y_con[i]);
}
