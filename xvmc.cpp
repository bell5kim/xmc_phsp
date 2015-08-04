/*****************************************************************************
 * xvmc.cpp: main program for XVMC 1.0                                       *
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
#include <stdio.h>

#include <iostream>
using namespace std;

#include "definitions.h"
#include "global.h"
#include "xvmc.h"

int main(int argc, char *argv[])
{
   char    *patient_name;       // patient or phantom name
   char    *plan_name;          // plane name

   if (argc != 3)
   {
      cerr << "XVMC> " << endl;
      cerr << "XVMC> Usage: xvmc <patient name> <plan name>" << endl;
      cerr << "XVMC> " << endl;
      exit(8);
   }
   patient_name = &argv[1][0];
   plan_name    = &argv[2][0];

   // do general initializations (read data, generate phantom. etc. )
   init_xvmc(patient_name, plan_name);

   // calculate plan or read dose distribution from file
   if (plan.calculate)
   {
      while (plan.total_number > 0)
      {
         // initialize present beam
         init_beam(plan.first_beam);
         if (plan.first_beam->n_history == 0)
         {
            // read dose distribution for this beam from file
            xvmc_message("==============================================",1);
            xvmc_message("= Reading 3D dose and error files of beam:",
                          plan.first_beam->id,"",0);
            xvmc_message("==============================================",0);

            // input parameters
            float dose_max = ZERO;
            float cpu_time = ZERO;
            int_3 d_max; d_max.x = 0; d_max.y = 0; d_max.z = 0;

            // generate file names
            int slength = strlen(inout_path)+10; // length of file name strings
            char *hed_name = new char[slength];  // dose matrix header file
            char *d3d_name = new char[slength];  // 3D dose data file
            char *e3d_name = new char[slength];  // 3D dose error file

            strcpy(hed_name,inout_path);
            strcpy(d3d_name,inout_path);
            strcpy(e3d_name,inout_path);

            int  file_id = plan.first_beam->id + 1;
            if ((0<file_id) && (file_id < 100))
            {
               char file_id_string[3];
               sprintf(file_id_string,"%02u",file_id);
               strcat(hed_name,file_id_string); strcat(hed_name,".hed");
               strcat(d3d_name,file_id_string); strcat(d3d_name,".d3d");
               strcat(e3d_name,file_id_string); strcat(e3d_name,".e3d");
            }
            else
            {
               if ((99<file_id) && (file_id < 1000))
               {
                  char file_id_string[4];
                  sprintf(file_id_string,"%03u",file_id);
                  strcat(hed_name,file_id_string); strcat(hed_name,".hed");
                  strcat(d3d_name,file_id_string); strcat(d3d_name,".d3d");
                  strcat(e3d_name,file_id_string); strcat(e3d_name,".e3d");
               }
               else
               {
                  xvmc_error("main","file Id out of range",8);
               }
            }

            // read beam dose distribution from file
            read_dose_file(beam_dose,beam_error,dose_max,d_max,
                           cpu_time,hed_name,d3d_name,e3d_name);

            plan.first_beam->cpu_time = cpu_time;

            // free memory
            delete [] hed_name;
            delete [] d3d_name;
            delete [] e3d_name;

         }
         else
         {
            calc_dose(plan.first_beam);
         }
         evaluate(plan,plan.first_beam);
         plan.erase_first();
         // delete accelerator head model(s) for this beam
         if (linac != NULL)
         {
            delete linac;
            linac  = NULL;
            elinac = NULL;
            plinac = NULL;
         }
      }
   }
   else
   {
      xvmc_message("=========================================",1);
      xvmc_message("= Reading total 3D dose and error files =",0);
      xvmc_message("=========================================",0);

      // input parameters
      float dose_max = ZERO;
      float cpu_time = ZERO;
      int_3 d_max; d_max.x = 0; d_max.y = 0; d_max.z = 0;

      // generate file names
      int slength = strlen(inout_path)+10;   // length of file name strings
      char *hed_name = new char[slength];    // dose matrix header file
      char *d3d_name = new char[slength];    // 3D dose data file
      char *e3d_name = new char[slength];    // 3D dose error file

      strcpy(hed_name,inout_path);
      strcpy(d3d_name,inout_path);
      strcpy(e3d_name,inout_path);

      strcat(hed_name,"00.hed");
      strcat(d3d_name,"00.d3d");
      strcat(e3d_name,"00.e3d");

      // read total dose distribution from file
      read_dose_file(sum_dose,sum_error,dose_max,d_max,
                     cpu_time,hed_name,d3d_name,e3d_name);

      plan.cpu_time = cpu_time;
      plan.dose_max = dose_max;
      plan.d_max.x  = d_max.x;
      plan.d_max.y  = d_max.y;
      plan.d_max.z  = d_max.z;

      // free memory
      delete [] hed_name;
      delete [] d3d_name;
      delete [] e3d_name;
   }

   // evaluate total dose distribution
   evaluate(plan);

   // free memeory of the photon beam modifier
   if (pbeam_modifier != NULL)
   {
      delete pbeam_modifier;
      pbeam_modifier = NULL;
   }

   // free memeory of the electron beam modifier
   if (ebeam_modifier != NULL)
   {
      delete ebeam_modifier;
      ebeam_modifier = NULL;
   }

   // free memory of 3D density and dose cubes
   delete density;    density    = NULL;
   delete dens_scol;  dens_scol  = NULL;
   delete dens_ccol;  dens_ccol  = NULL;
   delete dens_srad;  dens_srad  = NULL;
   delete dens_fchi;  dens_fchi  = NULL;
   delete dens_comp;  dens_comp  = NULL;
   delete dens_pair;  dens_pair  = NULL;
   delete dens_phot;  dens_phot  = NULL;
   delete beam_dose;  beam_dose  = NULL;
   delete beam_error; beam_error = NULL;
   delete sum_dose;   sum_dose   = NULL;
   delete sum_error;  sum_error  = NULL;

   xvmc_message("===========================",1);
   xvmc_message("= END OF DOSE CALCULATION =",0);
   xvmc_message("===========================",0);
   xvmc_message("",0);
}
