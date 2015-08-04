/*****************************************************************************
 * xvmc_aux.cpp: auxiliary functions                                         *
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

#include <string.h>

#include "definitions.h"
#include "global.h"

// calculate intersection plane of a line given by "start" and "dir"
// with calculation cube defined by the global variable "cube_size"
int intersection_plane(const real_3 &start, const real_3 &dir, real_3 &point)
{
   real   temp;

   // check z direction
   if (dir.z > ZERO)
   {
      point.z = ZERO;
      temp    = -start.z/dir.z;
      point.x =  start.x + temp*dir.x;
      point.y =  start.y + temp*dir.y;
      if (point.x >= ZERO && point.x <= cube_size.x &&
          point.y >= ZERO && point.y <= cube_size.y) return(1);
   }
   if (dir.z < ZERO)
   {
      point.z = cube_size.z;
      temp    = (cube_size.z-start.z)/dir.z;
      point.x =  start.x + temp*dir.x;
      point.y =  start.y + temp*dir.y;
      if (point.x >= ZERO && point.x <= cube_size.x &&
          point.y >= ZERO && point.y <= cube_size.y) return(6);
   }

   // the line does not go trough the upper or lower side
   // therefore, check x direction
   if (dir.x > ZERO)
   {
      point.x = ZERO;
      temp    = -start.x/dir.x;
      point.z =  start.z + temp*dir.z;
      point.y =  start.y + temp*dir.y;
      if (point.z >= ZERO && point.z <= cube_size.z &&
          point.y >= ZERO && point.y <= cube_size.y) return(2);
   }
   if (dir.x < ZERO)
   {
      point.x = cube_size.x;
      temp    = (cube_size.x-start.x)/dir.x;
      point.z =  start.z + temp*dir.z;
      point.y =  start.y + temp*dir.y;
      if (point.z >= ZERO && point.z <= cube_size.z &&
          point.y >= ZERO && point.y <= cube_size.y) return(3);
   }

   // the line does not go trough the upper, lower, left or right side
   // therefore, check y direction
   if (dir.y > ZERO)
   {
      point.y = ZERO;
      temp    = -start.y/dir.y;
      point.x =  start.x + temp*dir.x;
      point.z =  start.z + temp*dir.z;
      if (point.x >= ZERO && point.x <= cube_size.x &&
          point.z >= ZERO && point.z <= cube_size.z) return(4);
   }
   if (dir.y < ZERO)
   {
      point.y = cube_size.y;
      temp    = (cube_size.y-start.y)/dir.y;
      point.x =  start.x + temp*dir.x;
      point.z =  start.z + temp*dir.z;
      if (point.x >= ZERO && point.x <= cube_size.x &&
          point.z >= ZERO && point.z <= cube_size.z) return(5);
   }

   // the line doesn't touch the calculation cube
   return(0);
}

// trace a line from "start" to "stop" through the calculation cube
// and return the effective length (equivalent or radiological depth
// in water, approximated by the Compton length)
// the point "start" may be outside the cube, "stop" must be within the cube,
// "length" returns the geometrical length inside the cube
real trace_line(const real_3 &start, const real_3 &stop, real &length)
{
   // "stop" must be within the cube limits
   if ( (stop.x < ZERO) || (stop.x > cube_size.x) ||
        (stop.y < ZERO) || (stop.y > cube_size.y) ||
        (stop.z < ZERO) || (stop.z > cube_size.z) )
   {
      xvmc_error("trace_line",
                 "the end point is outside the calculation cube",8);
   }

   // new starting position if "start" is outside the cube limits
   real_3 pos;

   // check whether "start" is inside or outside
   if ( (start.x < ZERO) || (start.x > cube_size.x) ||
        (start.y < ZERO) || (start.y > cube_size.y) ||
        (start.z < ZERO) || (start.z > cube_size.z) )
   {
      // "start" is outside, find the intersection point with the
      // cube surface

      // difference vector and distance
      real_3 dss;
      dss.x = stop.x - start.x;
      dss.y = stop.y - start.y;
      dss.z = stop.z - start.z;
      real distance = sqrt(dss.x*dss.x + dss.y*dss.y + dss.z*dss.z);

      // normalize difference vector to get the direction cosines
      dss.x /= distance;
      dss.y /= distance;
      dss.z /= distance;

      // find intersection point of line with cube surface,
      if (!intersection_plane(start,dss,pos))
      {
         // in principle, this is impossible because "start" is outside
         // and "stop" is within the cube limits
         xvmc_error("trace_line",
                    "the line doesn't touch the calculation cube",8);
      }
   }
   else
   {
      // "start" is inside, "start" and "pos" are identical
      pos.x = start.x;
      pos.y = start.y;
      pos.z = start.z;
   }

   // new difference vector
   real_3 dps;
   dps.x = stop.x - pos.x;
   dps.y = stop.y - pos.y;
   dps.z = stop.z - pos.z;

   // geometrical length inside the cube
   length = sqrt(dps.x*dps.x + dps.y*dps.y + dps.z*dps.z);

   // voxel index for the starting position "pos"
   int_3 ipos;
   ipos.x = int(pos.x/voxel_size.x);
   ipos.y = int(pos.y/voxel_size.y);
   ipos.z = int(pos.z/voxel_size.z);
   if (ipos.x < 0) ipos.x = 0;
   if (ipos.y < 0) ipos.y = 0;
   if (ipos.z < 0) ipos.z = 0;
   if (ipos.x >= dim.x) ipos.x = dim.x-1;
   if (ipos.y >= dim.y) ipos.y = dim.y-1;
   if (ipos.z >= dim.z) ipos.z = dim.z-1;

   // voxel index for "stop"
   int_3 istop;
   istop.x = int(stop.x/voxel_size.x);
   istop.y = int(stop.y/voxel_size.y);
   istop.z = int(stop.z/voxel_size.z);
   if (istop.x < 0) istop.x = 0;
   if (istop.y < 0) istop.y = 0;
   if (istop.z < 0) istop.z = 0;
   if (istop.x >= dim.x) istop.x = dim.x-1;
   if (istop.y >= dim.y) istop.y = dim.y-1;
   if (istop.z >= dim.z) istop.z = dim.z-1;

   // the effective length is trivial if the two points are within
   // the same voxel
   if ( (ipos.x == istop.x) && (ipos.y == istop.y) && (ipos.z == istop.z) )
      return(length*dens_comp->matrix[ipos.x][ipos.y][ipos.z]);

   // direction cosines
   real_3 dir;
   dir.x = dps.x/length;
   dir.y = dps.y/length;
   dir.z = dps.z/length;

   // find the X-, Y-, Z-distance to the next voxel boundary
   real_3 step;   // step sizes to the next x-, y-, z-voxel boundaries
   int_3  istep;  // step direction to the next voxel boundary (-1,0,1

   // check X-direction
   if (dir.x > ZERO)
   {
      istep.x = 1;
      step.x  = (voxel_size.x*double(ipos.x+1)-pos.x)/dir.x;
   }
   else
   {
      if (dir.x < ZERO)
      {
         istep.x = -1;
         step.x  = (voxel_size.x*double(ipos.x)-pos.x)/dir.x;
      }
      else
      {
         istep.x = 0;
         step.x  = HUGE_STEP;
      }
   }

   // check Y-direction
   if (dir.y > ZERO)
   {
      istep.y = 1;
      step.y  = (voxel_size.y*double(ipos.y+1)-pos.y)/dir.y;
   }
   else
   {
      if (dir.y < ZERO)
      {
         istep.y = -1;
         step.y  = (voxel_size.y*double(ipos.y)-pos.y)/dir.y;
      }
      else
      {
         istep.y = 0;
         step.y  = HUGE_STEP;
      }
   }

   // check Z-direction
   if (dir.z > ZERO)
   {
      istep.z = 1;
      step.z  = (voxel_size.z*double(ipos.z+1)-pos.z)/dir.z;
   }
   else
   {
      if (dir.z < ZERO)
      {
         istep.z = -1;
         step.z  = (voxel_size.z*double(ipos.z)-pos.z)/dir.z;
      }
      else
      {
         istep.z = 0;
         step.z  = HUGE_STEP;
      }
   }

   // new voxel index
   int_3 inew = ipos;

   // step length in the present voxel
   real voxel_step = ZERO;

   // accumulated total geometrical length
   real tot_length = ZERO;

   // accumulated effective length
   real eff_length = ZERO;

   // start tracing
   bool repeat = true;
   while (repeat)
   {
      if ( (ipos.x < 0) || (ipos.x >= dim.x) ||
           (ipos.y < 0) || (ipos.y >= dim.y) ||
           (ipos.z < 0) || (ipos.z >= dim.z)    )
      {
         // the cube boundary has been reached, i.e. there must be an
         // error because the "stop" point is within the cube
         xvmc_error("trace_line","the cube boundary is reached",8);
      }
      else
      {
         if ( (ipos.x == istop.x) &&
              (ipos.y == istop.y) &&
              (ipos.z == istop.z)    )
         {
            // the "stop" voxel has been reached
            repeat = false;
            voxel_step  = length - tot_length;
            eff_length += voxel_step*dens_comp->matrix[ipos.x][ipos.y][ipos.z];
         }
         else
         {
            // find the next voxel
            if ( (step.z < step.x) && (step.z < step.y) )
            {
               voxel_step = step.z;
               step.x -= voxel_step;
               step.y -= voxel_step;
               step.z  = voxel_size.z/fabs(dir.z);
               inew.z  = ipos.z + istep.z;
            }
            else
            {
               if (step.x < step.y)
               {
                  voxel_step = step.x;
                  step.x  = voxel_size.x/fabs(dir.x);
                  step.y -= voxel_step;
                  step.z -= voxel_step;
                  inew.x  = ipos.x + istep.x;
               }
               else
               {
                  voxel_step = step.y;
                  step.x -= voxel_step;
                  step.y  = voxel_size.y/fabs(dir.y);
                  step.z -= voxel_step;
                  inew.y  = ipos.y + istep.y;
               }
            }

            // update total end effective length
            tot_length += voxel_step;
            eff_length += voxel_step*dens_comp->matrix[ipos.x][ipos.y][ipos.z];

            // update voxel index
            ipos = inew;
         }
      }
   } // while (repeat)

   return(eff_length);
}

// get interaction data file path for the specified file name and extension
char *get_file_path(const char *file_name, const char *file_ext)
{
   // get xvmc home directory name
   char *env_value = getenv("XVMC_HOME");
   if (env_value == NULL)
   {
      xvmc_error("get_file_path",
                 "environment variable XVMC_HOME not set",8);
   }

   // calculate string length for file path
   int string_slength = strlen(env_value) + strlen(file_name)
                      + strlen(file_ext)  + 10;

   // allocate memory for file path
   char *file_path = NULL;
   if ( (file_path = new char[string_slength]) == NULL )
   {
      xvmc_error("get_file_path",
                 "cannot allocate memory for file path",8);
   }

   // copy xvmc home directory name into file path string
   strcpy(file_path, env_value);

   // add file path, name and extension
   strcat(file_path, "/dat/");
   strcat(file_path, file_name);
   strcat(file_path, ".");
   strcat(file_path, file_ext);

   // return file path
   return(file_path);
}
