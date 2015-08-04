/*****************************************************************************
 * beam_model.cpp:                                                           *
 *    class member functions for:                                            *
 *       beam_model:         abstract base class for different beam models   *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 11.02.2000      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <math.h>
#include <ctype.h>
#include <string.h>

#include <sstream>
using namespace std;

#include <assert.h>

#include "definitions.h"
#include "global.h"
#include "xvmc_util.h"
#include "beam_core.h" // defines ONLY EMPTY class beam_core for Hyperion
#include "beam_model.h"

// ****************************************
// declare functions and global variables
// ****************************************

// ****************************************
// member functions of class beam_model
// ****************************************

// construct beam model from beam parameters
beam_model::beam_model(const beam_core *this_beam)
{
   iso_distance       = 100.0; // default 100cm
   type               = this_beam->type;
   nominal_energy     = this_beam->energy; 

   // the next two parameters can be defined only for beam models based
   // on measurements under reference conditions in water,
   // therefore these default settings have no meaning
   norm_value         = 0.01;
   gray_mu_dmax       = 0.01;

   iso_center.x       = this_beam->iso_center.x;
   iso_center.y       = this_beam->iso_center.y;
   iso_center.z       = this_beam->iso_center.z;
   gantry_rotation    = this_beam->gantry_rotation;
   start_gantry_angle = this_beam->start_gantry_angle*ONE_PI/180.0;
   stop_gantry_angle  = this_beam->stop_gantry_angle*ONE_PI/180.0;
   // estimate an average gantry angle in the case of gantry rotation
   if (gantry_rotation)
   {
      delta_gantry_angle =  stop_gantry_angle - start_gantry_angle;
      gantry_angle       = (stop_gantry_angle + start_gantry_angle)/TWO;
   }
   else
   {
      stop_gantry_angle  = start_gantry_angle;
      delta_gantry_angle = ZERO;
      gantry_angle       = start_gantry_angle;
   }
   table_angle        = this_beam->table_angle*ONE_PI/180.0;
   collimator_angle   = this_beam->collimator_angle*ONE_PI/180.0;
   open_x1            = this_beam->open_x1;
   open_x2            = this_beam->open_x2;
   open_y1            = this_beam->open_y1;
   open_y2            = this_beam->open_y2;
   cos_alpha          = cos(gantry_angle);
   sin_alpha          = sin(gantry_angle);
   cos_beta           = cos(table_angle);
   sin_beta           = sin(table_angle);
   cos_gamma          = cos(collimator_angle);
   sin_gamma          = sin(collimator_angle);
   beam_dir.x         = -sin_alpha*cos_beta;
   beam_dir.y         =  sin_alpha*sin_beta;
   beam_dir.z         =  cos_alpha;
   origin.x           = iso_center.x - iso_distance*beam_dir.x;
   origin.y           = iso_center.y - iso_distance*beam_dir.y;
   origin.z           = iso_center.z - iso_distance*beam_dir.z;
   theta_x1           = atan(open_x1/iso_distance);
   theta_x2           = atan(open_x2/iso_distance);
   theta_y1           = atan(open_y1/iso_distance);
   theta_y2           = atan(open_y2/iso_distance);
   theta_width_x      = theta_x2 - theta_x1;
   theta_width_y      = theta_y2 - theta_y1;

   num_warnings       = 0;
}

// open file and search for base data entry
void beam_model::find_base_data(const beam_core *this_beam)
{
   extern char *base_path;      // path to the base data files
   char        *base_file_name; // base data file name
   char         line[81] = "";  // lines to read from file
  
   if (this_beam->base_key == NULL)
   {
      xvmc_error("beam_model::find_base_data",
                 "base data file key not specified",8);
   }
  
   // allocate memory for base data file name
   if ( (base_file_name = new char[
                                 strlen(base_path)           // file path
                               + strlen(this_beam->base_key) // file name
                               + 5                           // file extension
      ]) == NULL)
   {
      xvmc_error("beam_model::find_base_data",
                 "cannot allocate memory for base data file name",8);
   }
  
   // create file name
   strcpy(base_file_name, base_path);
   strcat(base_file_name, this_beam->base_key);
   strcat(base_file_name, ".bdt");
  
   xvmc_message("Opening base data file:",base_file_name,1);
  
   // opening file
   base_file.open(base_file_name,ios::in);
   if (!base_file)
   {
      xvmc_error("beam_model::find_base_data",
                 "cannot open base data file",8);
   }
  
   // delete file name
   delete [] base_file_name;
   base_file_name = NULL;
  
   // find "BASE-DATA-FILE-VERSION"
   bool read_line = true;
   while (read_line)
   {
      if (base_file.eof())
      {
         xvmc_error("beam_model::find_base_data",
                    "BASE-DATA-FILE-VERSION entry not found",8);
      }
      base_file.getline(line,sizeof(line));
      istringstream line_stream(line);
      char keyword[81] = "";
      line_stream >> keyword;
      if (!strcmp(keyword,"BASE-DATA-FILE-VERSION:"))
      {
         char version[4] = "";
         line_stream >> version;
         if (!strcmp(version,"1.4")) read_line = false;
         else xvmc_error("beam_model::find_base_data",
                         "wrong base data file version",8);
      }
   }
  
   // here it is clear that the base data file version is correct
   // now we try to find a basa data entry for the present beam
  
   // define bool variables
   bool base_data_found      = false;
   bool particle_type_found  = false;
   bool nominal_energy_found = false;
   bool beam_model_id_found  = false;
   // for electron beams we also need an applicator
   bool applicator_found     = false;

   // read further lines
   read_line = true;
   while (read_line)
   {
      if (base_file.eof())
      {
         xvmc_error("beam_model::find_base_data",
                    "base data for this beam not found",8);
      }
      base_file.getline(line,sizeof(line));
      istringstream line_stream(line);
      char keyword[81] = "";
      line_stream >> keyword;
      if (!strcmp(keyword,"BASE-DATA-ENTRY"))
      {
         base_data_found      = true;
         particle_type_found  = false;
         nominal_energy_found = false;
         // for mono-energetic beams we take the nominal energy
         if (this_beam->model_id <= 0) nominal_energy_found = true;
         beam_model_id_found  = false;
         // for photon beams we don't need to search for an applicator
         if (this_beam->type == ELECTRON) applicator_found = false;
         else                             applicator_found = true;
      }
      else
      {
         if (base_data_found)
         {
            if (!strcmp(keyword,"PARTICLE-TYPE:"))
            {
               char particle_type[10];
               line_stream >> particle_type;
               // convert string to lower case
               for (register unsigned int i=0; i<strlen(particle_type); ++i)
                  particle_type[i] = tolower(particle_type[i]);
               if (this_beam->type == PHOTON)
               {
                  if (!strncmp(particle_type,"photon",6))
                     particle_type_found  = true;
               }
               if (this_beam->type == ELECTRON)
               {
                  if (!strncmp(particle_type,"electron",8))
                     particle_type_found  = true;
               }
            }
            if (!strcmp(keyword,"NOMINAL-ENERGY:"))
            {
               real nominal_energy;
               line_stream >> nominal_energy;
               if (fabs(this_beam->energy-nominal_energy) < 0.01)
                  nominal_energy_found = true;
            }
            if (!strcmp(keyword,"BEAM-MODEL-ID:"))
            {
               int beam_model_id;
               line_stream >> beam_model_id;
               if (this_beam->model_id == beam_model_id)
                  beam_model_id_found = true;
            }
            if (!strcmp(keyword,"APPLICATOR:"))
            {
               char applicator[10];
               line_stream >> applicator;
               if (this_beam->type == ELECTRON)
               {
                  if (this_beam->applicator == NULL)
                  {
                     xvmc_error("beam_model::find_base_data",
                            "an electron applicator has not been defined",8);
                  }
                  else if (!strncmp(this_beam->applicator,applicator,8))
                  {
                     applicator_found  = true;
                  }
               }
            }
            if (!strcmp(keyword,"BEGIN-PARAMETERS"))
            {
               if ( particle_type_found  &&
                    nominal_energy_found &&
                    beam_model_id_found  &&
                    applicator_found        )
               {
                  // stop file search because an entry
                  // for this beam was detected
                  read_line = false;
               }
               else
               {
                  // this entry is not for the present beam
                  // therefore the bool variables are reset
                  base_data_found      = false;
                  particle_type_found  = false;
                  nominal_energy_found = false;
                  // for mono-energetic beams we take the nominal energy
                  if (this_beam->model_id <= 0) nominal_energy_found = true;
                  beam_model_id_found  = false;
                  // for photon beams we don't need to search for an applicator
                  if (this_beam->type == ELECTRON) applicator_found = false;
                  else                             applicator_found = true;
               }
            }
         } // if (base_data_found)
      } // if (!strcmp(keyword,"BASE-DATA-ENTRY"))
   } // while (read_line)
}

// trace particle to the simulation grid (take beam angles into account)
bool beam_model::trace2cube(real_3 &pos, real_3 &dir, int_3 &i,
                            real origin_point_dist, ranmar &rndm)
{
   real       temp;
   real_3     escape;             // linac escaping position
   // check for collimator-patient collisions
   const real DENSITY_THRESHOLD = 0.1;
   // maximum number of collision warnings
   const int  MAX_WARNINGS = 100;

   // in the case of gantry rotation sample new gantry angle
   if (gantry_rotation)
   {
      gantry_angle = start_gantry_angle+rndm.number()*delta_gantry_angle;
      cos_alpha    = cos(gantry_angle);
      sin_alpha    = sin(gantry_angle);
      beam_dir.x   = -sin_alpha*cos_beta;
      beam_dir.y   =  sin_alpha*sin_beta;
      beam_dir.z   =  cos_alpha;
      origin.x     = iso_center.x - iso_distance*beam_dir.x;
      origin.y     = iso_center.y - iso_distance*beam_dir.y;
      origin.z     = iso_center.z - iso_distance*beam_dir.z;
   }
// printf ("a %f %f %f p %f %f %f opd %f\n", origin.x, origin.y, origin.z, pos.x, pos.y, pos.z, origin_point_dist); 

   // normalize position vector
   pos.x  = pos.x/origin_point_dist;
   pos.y  = pos.y/origin_point_dist;
   pos.z  = pos.z/origin_point_dist;

   // rotate position vector by the collimator angle
   temp  = pos.x*sin_gamma;
   pos.x = pos.x*cos_gamma + pos.y*sin_gamma;
   pos.y = -temp           + pos.y*cos_gamma;

   // rotate direction vector by the collimator angle
   temp  = dir.x*sin_gamma;
   dir.x = dir.x*cos_gamma + dir.y*sin_gamma;
   dir.y = -temp           + dir.y*cos_gamma;

   // rotate position vector by the gantry angle
   temp  = pos.x*sin_alpha;
   pos.x = pos.x*cos_alpha - pos.z*sin_alpha;
   pos.z = temp            + pos.z*cos_alpha;

   // rotate direction vector by the gantry angle
   temp  = dir.x*sin_alpha;
   dir.x = dir.x*cos_alpha - dir.z*sin_alpha;
   dir.z = temp            + dir.z*cos_alpha;

   // rotate position vector by the table angle
   temp  = pos.x*sin_beta;
   pos.x = pos.x*cos_beta + pos.y*sin_beta;
   pos.y = -temp          + pos.y*cos_beta;

   // rotate direction vector by the table angle
   temp  = dir.x*sin_beta;
   dir.x = dir.x*cos_beta + dir.y*sin_beta;
   dir.y = -temp          + dir.y*cos_beta;

   // calculate the accelerator escaping position of this particle
   pos.x = origin.x + origin_point_dist*pos.x;
   pos.y = origin.y + origin_point_dist*pos.y;
   pos.z = origin.z + origin_point_dist*pos.z;
// printf ("o %f %f %f p %f %f %f opd %f\n", origin.x, origin.y, origin.z, pos.x, pos.y, pos.z, origin_point_dist); 
   // test whether this particle is inside the calculation cube
   if ( (pos.x > ZERO) && (pos.x < cube_size.x) &&
        (pos.y > ZERO) && (pos.y < cube_size.y) &&
        (pos.z > ZERO) && (pos.z < cube_size.z) )
   {
      // the particle is within the cube, calculate voxel indices and return
      i.x = int(pos.x/voxel_size.x);
      i.y = int(pos.y/voxel_size.y);
      i.z = int(pos.z/voxel_size.z);

      // check for collimator-patient collisions
      if ( density->matrix[i.x][i.y][i.z] > DENSITY_THRESHOLD )
      {
         // the particle is probably in the patient, print a warning
         if (num_warnings < MAX_WARNINGS)
         {
            ++num_warnings;
            xvmc_warning("beam_model::trace2cube",
                         "collimator within patient contour",0);
         }
      }
      // now return
      return(true);

   }

   // the particle is outside the calculation cube, therefore
   // find the surface point where this particle hits the cube

   escape.x = pos.x;
   escape.y = pos.y;
   escape.z = pos.z;

   // check in +z direction
   if (dir.z > ZERO)
   {
      temp  = escape.z/dir.z;
      pos.x = escape.x - temp*dir.x;
      pos.y = escape.y - temp*dir.y;
      if ( (pos.x > ZERO) && (pos.x < cube_size.x) &&
           (pos.y > ZERO) && (pos.y < cube_size.y) )
      {
         pos.z = ZERO;
         i.x = int(pos.x/voxel_size.x);
         i.y = int(pos.y/voxel_size.y);
         i.z = 0;
         return(true);
      }
   }

   // check in -z direction
   if (dir.z < ZERO)
   {
      temp  = (cube_size.z-escape.z)/dir.z;
      pos.x = escape.x + temp*dir.x;
      pos.y = escape.y + temp*dir.y;
      if ( (pos.x > ZERO) && (pos.x < cube_size.x) &&
           (pos.y > ZERO) && (pos.y < cube_size.y) )
      {
         pos.z = cube_size.z;
         i.x = int(pos.x/voxel_size.x);
         i.y = int(pos.y/voxel_size.y);
         i.z = dim.z-1;
         return(true);
      }
   }

   // the particle does not go trough the upper/lower side
   // check in +x direction
   if (dir.x > ZERO)
   {
      temp  = escape.x/dir.x;
      pos.y = escape.y - temp*dir.y;
      pos.z = escape.z - temp*dir.z;
      if ( (pos.y > ZERO) && (pos.y < cube_size.y) &&
           (pos.z > ZERO) && (pos.z < cube_size.z) )
      {
         pos.x = ZERO;
         i.x = 0;
         i.y = int(pos.y/voxel_size.y);
         i.z = int(pos.z/voxel_size.z);
         return(true);
      }
   }

   // check in -x direction
   if (dir.x < ZERO)
   {
      temp  = (cube_size.x-escape.x)/dir.x;
      pos.y = escape.y + temp*dir.y;
      pos.z = escape.z + temp*dir.z;
      if ( (pos.y > ZERO) && (pos.y < cube_size.y) &&
           (pos.z > ZERO) && (pos.z < cube_size.z) )
      {
         pos.x = cube_size.x;
         i.x = dim.x-1;
         i.y = int(pos.y/voxel_size.y);
         i.z = int(pos.z/voxel_size.z);
         return(true);
      }
   }

   // the particle does not go trough the left/right side
   // check in +y direction
   if (dir.y > ZERO)
   {
      temp  = escape.y/dir.y;
      pos.x = escape.x - temp*dir.x;
      pos.z = escape.z - temp*dir.z;
      if ( (pos.x > ZERO) && (pos.x < cube_size.x) &&
           (pos.z > ZERO) && (pos.z < cube_size.z) )
      {
         pos.y = ZERO;
         i.x = int(pos.x/voxel_size.x);
         i.y = 0;
         i.z = int(pos.z/voxel_size.z);
         return(true);
      }
   }

   // check in -y direction
   if (dir.y < ZERO)
   {
      temp  = (cube_size.y-escape.y)/dir.y;
      pos.x = escape.x + temp*dir.x;
      pos.z = escape.z + temp*dir.z;
      if ( (pos.x > ZERO) && (pos.x < cube_size.x) &&
           (pos.z > ZERO) && (pos.z < cube_size.z) )
      {
         pos.y = cube_size.y;
         i.x = int(pos.x/voxel_size.x);
         i.y = dim.y-1;
         i.z = int(pos.z/voxel_size.z);
         return(true);
      }
   }

   // now, it is clear that this history never hits the calculation cube
   return(false);
}
