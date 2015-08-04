#ifndef _CONTOUR_H_
#define _CONTOUR_H_

/*****************************************************************************
 * contour.h:                                                                *
 *    class declarations and inline member functions for:                    *
 *       class contour:  contour defined by a set of (x,y) points            *
 *                                                                           *
 * Copyright (C) 2001    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 12.06.2001      *
 *                                                                           *
 *****************************************************************************/

#include "definitions.h"
#include "global.h"

// contour defined by a set of (x,y) points
class contour
{
   public:
      contour(const int);      // define a contour by a given number of points
      ~contour(void);          // delete contour
      // change one contour point, if this is impossible -> return false
      bool change_point(int, real, real);
      // close contour, if this is impossible -> return false
      bool close(void);
      int  get_num(void);       // get number of contour points
      real get_x(int);         // get x-coordinate of point i
      real get_y(int);         // get y-coordinate of point i
      // return true if a point (x0,y0) is located inside the contour,
      // note that the contour must be closed to run this function
      bool inside(real, real);
   private:
      bool    closed;          // indicates an open or closed contour
      int     num_points;      // number of contour points
      int     array_size;      // size of x and y arrays
      real   *x_con;           // array of contour points (x-coordinate)
      real   *y_con;           // array of contour points (y-coordinate)
      real    min_x;           // minimum of all x-coordinates
      real    max_x;           // maximum of all x-coordinates
      real    min_y;           // minimum of all y-coordinates
      real    max_y;           // maximum of all y-coordinates
      real   *x_inter;         // array to store intersection points
      real   *x_min,*x_max;    // arrays to store pairs of intersection points
};

// return true if a point (x0,y0) is located inside the contour,
// note that the contour must be closed to run this function
inline bool contour::inside(real x0, real y0)
{
   if (!closed)
   {
      xvmc_warning("contour::inside","contour not closed",1);
      return(false);
   }

   // test, whether (x0,y0) is outside the rectangle,
   // defined by min_x,max_x,min_y,max_y (this saves compuation time)
   if (x0 < min_x) return(false);
   if (x0 > max_x) return(false);
   if (y0 < min_y) return(false);
   if (y0 > max_y) return(false);

   // find the intersections between the line y=y0 and the contour
   int n_inter = 0; // number of intersections
   real f_old      = y_con[0]            - y0;
   real f_very_old = y_con[num_points-2] - y0;
   for (register int i=1; i<num_points; ++i)
   {
      real f_new =  y_con[i] - y0;

      // two intersection points are calculated, if the contour touches
      // the line y=y0
      if (f_old == ZERO)
      {
         if (f_very_old*f_new >= ZERO)
         {
           x_inter[n_inter] = x_con[i-1]; ++n_inter;
         }
      }

      // there is an intersection, if f_new is zero
      if (f_new == ZERO)
      {
         x_inter[n_inter] = x_con[i]; ++n_inter;
      }
      else
      {
         // a change of the sign of f indicates an intersection
         if (f_old*f_new < ZERO)
         {
            real f_ratio = f_old/(f_old-f_new);
            x_inter[n_inter] = x_con[i-1] + f_ratio*(x_con[i]-x_con[i-1]);
            ++n_inter;
         }
      }
      f_very_old = f_old;
      f_old      = f_new;
   }

   // determine the number of intersection pairs x_min,x_max
   if ( (n_inter%2) != 0)
      xvmc_error("contour::inside","odd number of intersection points",8);
   int n_pair = n_inter/2;

   // sort the intersection points and create pairs x_min,x_max
   real large = max_x + 10.0; // just a value larger than max_x
   for (register int i=0; i<n_pair; ++i)
   {
      // find x_min
      int j_min = 0;
      x_min[i]  = large;
      for (register int j=0; j<n_inter; ++j)
      {
         if (x_inter[j] < x_min[i])
         {
            x_min[i] = x_inter[j];
            j_min    = j;
         }
      }
      // fill array element with a large value
      x_inter[j_min] = large + 10;

      // find x_max
      int j_max  = 0;
      x_max[i]   = large;
      for (register int j=0; j<n_inter; ++j)
      {
         if (x_inter[j] < x_max[i])
         {
            x_max[i] = x_inter[j];
            j_max    = j;
         }
      }
      // fill array element with a large value
      x_inter[j_max] = large + 10;
   }

   // check for all intersection pairs whether x_min <= x0 <= x_max
   for (register int i=0; i<n_pair; ++i)
   {
      if ( (x_min[i] <= x0) && (x_max[i] >= x0) ) return(true);
   }

   // the point (x0,y0) is outside
   return(false);
}

#endif /* _CONTOUR_H_ */
