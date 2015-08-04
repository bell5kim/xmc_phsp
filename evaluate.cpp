/*****************************************************************************
 * evaluate.cpp:                                                             *
 *    functions evaluate(plan,beam): evaluate and save beam dose             *
 *              evaluate(plan):      evaluate and save plan dose             *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 00/12/06        *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "definitions.h"
#include "global.h"
#include "treatment_plan.h"
#include "beam_model.h"

// ****************************************
// declare functions and global variables
// ****************************************

// write dose matrix file
void write_dose_file(const array_3d<float> *, const array_3d<float> *,
                     const float, const int_3, const float,
                     const char *, const char *, const char *);

// --- Added by JOKim 16Nov2010 ----------------------------------
// write portal dose matrix file
void write_portal_file(const array_3d<float> *, const array_3d<float> *,
                     const float, const int_3, const float,
                     const char *, const char *, const char *, const char *);
// --- End of Adding 16Nov2010 -----------------------------------

// write dose plane
void write_plane(const array_3d<float> *,  const array_3d<float> *,
                 const plane_parameters *, const char *);

// write dose profile
void write_profile(const array_3d<float> *,    const array_3d<float> *,
                   const profile_parameters *, const char *);

// ****************************************
// function evaluate(plan,beam)
// ****************************************

void evaluate(treatment_plan &plan, beam_core *this_beam)
{
   // calculate 3D dose error array, determine dose maximum and
   // add dose and dose error to sum dose and sum error arrays
   float dose,error;
   float dose_max = ZERO;
   int_3 d_max        = {0, 0, 0}; // indices of the dose maximum voxel
   float dose_max_portal = ZERO;   // Added by JOKim 16Nov2010
   int_3 d_max_portal = {0, 0, 0}; // Added by JOKim 16Nov2010

   // dose unit conversion factor(s)
   // for absolute dose calc. we convert dose units from Gy/fluence -> Gy/MUs
   float  conv = ONE;
   extern beam_model *linac;
#ifdef MONITORBSF   
   if (plan.result == ABS_DOSE || plan.result == ABS_DOSE_MBSF)
#else
   if (plan.result == ABS_DOSE)
#endif   
   {
      conv = linac->get_gray()/linac->get_norm();
   }
   float  con2 = conv*conv;

   for (register int k=0; k<dim.z; ++k)
   {
      for (register int j=0; j<dim.y; ++j)
      {
         for (register int i=0; i<dim.x; ++i)
         {
            dose     = beam_dose->matrix[i][j][k];
            if (dose > dose_max)
            {
               dose_max = dose;
               d_max.x  = i;
               d_max.y  = j;
               d_max.z  = k;
            }
            if (this_beam->n_history == 0)
            {
               // we have read the dose error matrix from file
               error = beam_error->matrix[i][j][k];
            }
            else
            {
               // convert average dose squared to standard deviation (error)
               error = float(this_beam->n_batch)*beam_error->matrix[i][j][k];
               error = sqrt(fabs(error-dose*dose)/float(this_beam->n_batch-1));
               beam_error->matrix[i][j][k] = error;
            }
            sum_dose->matrix[i][j][k] += this_beam->weight*dose*conv;
            sum_error->matrix[i][j][k] +=
               this_beam->weight*this_beam->weight*error*error*con2;
         }
      }
   }
// --- Added by JOKim 16Nov2010 ------------------------------------
 
   for (register int k=0; k<dim_portal.z; ++k)
   {
      for (register int j=0; j<dim_portal.y; ++j)
      {
         for (register int i=0; i<dim_portal.x; ++i)
         {
            dose     = beam_dose_portal->matrix[i][j][k];
            if (dose > dose_max_portal)
            {
               dose_max_portal = dose;
               d_max_portal.x  = i;
               d_max_portal.y  = j;
               d_max_portal.z  = k;
            }
            // convert average dose squared to standard deviation (error)
            error = float(this_beam->n_batch)*beam_error_portal->matrix[i][j][k];
            error = sqrt(fabs(error-dose*dose)/float(this_beam->n_batch-1));
            beam_error_portal->matrix[i][j][k] = error;

            sum_dose_portal->matrix[i][j][k] += this_beam->weight*dose*conv;
            sum_error_portal->matrix[i][j][k] +=
               this_beam->weight*this_beam->weight*error*error*con2;
         }
      }
   }
// --- End of Adding 16Nov2010 -------------------------------------
   // get cpu_time
   float cpu_time = this_beam->cpu_time;
   plan.cpu_time += cpu_time;

   //
   // dose at the reference point
   // 

   // calculate voxel indices for x interpolation
   int lower_x = int((ref_point.x-voxel_size.x/TWO)/voxel_size.x);
   int upper_x = lower_x+1;
   if (lower_x <      0) { lower_x =       0; upper_x =       0; }
   if (upper_x >= dim.x) { lower_x = dim.x-1; upper_x = dim.x-1; }

   // interpolation factors
   real px = real(upper_x)+ONE_HALF-ref_point.x/voxel_size.x;
   real qx = ONE-px;

   // calculate voxel indices for y interpolation
   int lower_y = int((ref_point.y-voxel_size.y/TWO)/voxel_size.y);
   int upper_y = lower_y+1;
   if (lower_y <      0) { lower_y =       0; upper_y =       0; }
   if (upper_y >= dim.y) { lower_y = dim.y-1; upper_y = dim.y-1; }

   // interpolation factors
   real py = real(upper_y)+ONE_HALF-ref_point.y/voxel_size.y;
   real qy = ONE-py;

   // calculate voxel indices for z interpolation
   int lower_z = int((ref_point.z-voxel_size.z/TWO)/voxel_size.z);
   int upper_z = lower_z+1;
   if (lower_z <      0) { lower_z =       0; upper_z =       0; }
   if (upper_z >= dim.z) { lower_z = dim.z-1; upper_z = dim.z-1; }

   // interpolation factors
   real pz = real(upper_z)+ONE_HALF-ref_point.z/voxel_size.z;
   real qz = ONE-pz;

   // dose at the reference point
   float ref_dose  = beam_dose->matrix[lower_x][lower_y][lower_z]*px*py*pz
                   + beam_dose->matrix[lower_x][lower_y][upper_z]*px*py*qz
                   + beam_dose->matrix[lower_x][upper_y][lower_z]*px*qy*pz
                   + beam_dose->matrix[lower_x][upper_y][upper_z]*px*qy*qz
                   + beam_dose->matrix[upper_x][lower_y][lower_z]*qx*py*pz
                   + beam_dose->matrix[upper_x][lower_y][upper_z]*qx*py*qz
                   + beam_dose->matrix[upper_x][upper_y][lower_z]*qx*qy*pz
                   + beam_dose->matrix[upper_x][upper_y][upper_z]*qx*qy*qz;

   // error at the reference point
   float ref_error = beam_error->matrix[lower_x][lower_y][lower_z]*px*py*pz
                   + beam_error->matrix[lower_x][lower_y][upper_z]*px*py*qz
                   + beam_error->matrix[lower_x][upper_y][lower_z]*px*qy*pz
                   + beam_error->matrix[lower_x][upper_y][upper_z]*px*qy*qz
                   + beam_error->matrix[upper_x][lower_y][lower_z]*qx*py*pz
                   + beam_error->matrix[upper_x][lower_y][upper_z]*qx*py*qz
                   + beam_error->matrix[upper_x][upper_y][lower_z]*qx*qy*pz
                   + beam_error->matrix[upper_x][upper_y][upper_z]*qx*qy*qz;

   // output
   xvmc_message("Total CPU time:",cpu_time,"s",1);
#ifdef MONITORBSF   
   if (plan.result == ABS_DOSE || plan.result == ABS_DOSE_MBSF)
#else
   if (plan.result == ABS_DOSE)
#endif      
   {
      float mon_units = this_beam->weight*plan.sum_monunits*plan.num_fractions;
      xvmc_message("Maximum dose:  ",dose_max*conv*mon_units,"Gy",0);
      xvmc_message("Reference dose:",ref_dose*conv*mon_units,"+/-",
                                     ref_error*conv*mon_units,"Gy",1);
      xvmc_message("               ",ref_error/ref_dose*100.0,"(%)",0);
   }
   else
   {
      xvmc_message("Maximum dose:  ",dose_max,"10^-10 Gy cm^2",0);
      xvmc_message("Reference dose:",ref_dose,"+/-",
                                     ref_error,"10^-10 Gy cm^2",1);
   }

   // determine the average statistical uncertainty
   // for all voxels with D > D_max/2
   double average_error = ZERO;
   long   num_voxel     = 0;
   float  dose_limit    = dose_max/TWO;
   for (register int k=0; k<dim.z; ++k)
   {
      for (register int j=0; j<dim.y; ++j)
      {
         for (register int i=0; i<dim.x; ++i)
         {
            dose  = beam_dose->matrix[i][j][k];
            error = beam_error->matrix[i][j][k];
            if (dose > dose_limit)
            {
               average_error += error*error/(dose*dose);
               ++num_voxel;
            }
         }
      }
   }
   // output analysis, 0.888 is taken from the EGS4 timing benchmark test
   average_error = sqrt(average_error/double(num_voxel));
   xvmc_message("Average ICCR error:      ",average_error*100.0,"%",1);
#ifdef ICCR_TEST
   //real cpu_500  = cpu_time*0.888;
   real cpu_500  = cpu_time*1.7;
   real eff_500  = ONE/(cpu_500*average_error*average_error);
   real cpu_two  = 2500.0/eff_500;
   xvmc_message("CPU time for PIII 500MHz:",cpu_500,"s",0);
   xvmc_message("Efficiency:              ",eff_500,"s^(-1)",0);
   xvmc_message("CPU time for 2% error:   ",cpu_two,"s",0);
#else
   real efficiency  = ONE/(cpu_time*average_error*average_error);
   real cpu_two     = 2500.0/efficiency;
   xvmc_message("Efficiency:              ",efficiency,"s^(-1)",0);
   xvmc_message("CPU time for 2% error:   ",cpu_two,"s",0);
#endif // ICCR_TEST

   // write 3D dose matrix for each beam if required
   if (plan.out3d_beam)
   {
      xvmc_message("Creating 3D dose and error files for this beam",1);
      xvmc_message("==============================================",0);

      // generate file names and save data
      int slength = strlen(inout_path)+10;   // length of file name strings
      char *hed_name = new char[slength];    // dose matrix header file
      char *d3d_name = new char[slength];    // 3D dose data file
      char *e3d_name = new char[slength];    // 3D dose error file

      strcpy(hed_name,inout_path);
      strcpy(d3d_name,inout_path);
      strcpy(e3d_name,inout_path);

      int  file_id = this_beam->id + 1;
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
            xvmc_error("evaluate(plan,beam)","file Id out of range",8);
         }
      }

      // write beam dose distribution to file
      write_dose_file(beam_dose,beam_error,dose_max,d_max,
                      cpu_time,hed_name,d3d_name,e3d_name);

      // free memory
      delete [] hed_name;
      delete [] d3d_name;
      delete [] e3d_name;

      // write portal dose matrix for each beam if required
      if (this_beam->portal != NULL)
      {
	 xvmc_message("Creating portal dose and error files for this beam",1);
	 xvmc_message("==================================================",0);

	 // generate file names and save data
	 int slength = strlen(inout_path)+10;   // length of file name strings
	 char *pdh_name = new char[slength];    // portal dose matrix header file
	 char *pdi_name = new char[slength];    // portal dose data file
	 char *pde_name = new char[slength];    // portal dose error file
	 char *bmp_name = new char[slength];    // portal dose BMP file

	 strcpy(pdh_name,inout_path);
	 strcpy(pdi_name,inout_path);
	 strcpy(pde_name,inout_path);
	 strcpy(bmp_name,inout_path);

	 int  file_id = this_beam->id + 1;
	 if ((0<file_id) && (file_id < 100))
	 {
            char file_id_string[3];
            sprintf(file_id_string,"%02u",file_id);
            strcat(pdh_name,file_id_string); strcat(pdh_name,".pdh");
            strcat(pdi_name,file_id_string); strcat(pdi_name,".pdi");
            strcat(pde_name,file_id_string); strcat(pde_name,".pde");
            strcat(bmp_name,file_id_string); strcat(bmp_name,".bmp");
	 }
	 else
	 {
            if ((99<file_id) && (file_id < 1000))
            {
               char file_id_string[4];
               sprintf(file_id_string,"%03u",file_id);
               strcat(pdh_name,file_id_string); strcat(pdh_name,".pdh");
               strcat(pdi_name,file_id_string); strcat(pdi_name,".pdi");
               strcat(pde_name,file_id_string); strcat(pde_name,".pde");
               strcat(bmp_name,file_id_string); strcat(bmp_name,".bmp");
            }
            else
            {
               xvmc_error("evaluate(plan,beam)","file Id out of range",8);
            }
	 }

	 // write beam dose distribution to file
	 write_portal_file(beam_dose_portal, beam_error_portal,
                         dose_max_portal, d_max_portal, plan.portal->distance,
                	 pdh_name, pdi_name, pde_name, bmp_name);

	 // free memory
	 delete [] pdh_name;
	 delete [] pdi_name;
	 delete [] pde_name;
	 delete [] bmp_name;
      }
   }
}

// ****************************************
// function evaluate(plan)
// ****************************************

void evaluate(treatment_plan &plan)
{
   // calculate sum dose error array and determine dose maximum
   float dose;
   float dose_max = ZERO;
   int_3 d_max    = {0, 0, 0}; // indices of the dose maximum voxel
   float dose_max_portal = ZERO;   // Added by JOKim 16Nov2010
   int_3 d_max_portal = {0, 0, 0}; // Added by JOKim 16Nov2010

   if (plan.calculate)
   {
      // for abs dose calc. we multiply with the number of MUs
      float mon_units = ONE;
#ifdef MONITORBSF   
   if (plan.result == ABS_DOSE || plan.result == ABS_DOSE_MBSF)
#else
   if (plan.result == ABS_DOSE)
#endif            
      mon_units = plan.sum_monunits * plan.num_fractions;
      for (register int k=0; k<dim.z; ++k)
      {
         for (register int j=0; j<dim.y; ++j)
         {
            for (register int i=0; i<dim.x; ++i)
            {
               dose     = sum_dose->matrix[i][j][k]*mon_units;
               if (dose > dose_max)
               {
                  dose_max = dose;
                  d_max.x  = i;
                  d_max.y  = j;
                  d_max.z  = k;
               }
               sum_dose->matrix[i][j][k]  = dose;
               sum_error->matrix[i][j][k] =
                  sqrt(sum_error->matrix[i][j][k])*mon_units;
            }
         }
      }
// --- Added by JOKim 16Nov2010 --------------------------------
      for (register int k=0; k<dim_portal.z; ++k)
      {
         for (register int j=0; j<dim_portal.y; ++j)
         {
            for (register int i=0; i<dim_portal.x; ++i)
            {
               dose     = sum_dose_portal->matrix[i][j][k]*mon_units;
               if (dose > dose_max_portal)
               {
                  dose_max_portal = dose;
                  d_max_portal.x  = i;
                  d_max_portal.y  = j;
                  d_max_portal.z  = k;
               }
               sum_dose_portal->matrix[i][j][k]  = dose;
               sum_error_portal->matrix[i][j][k] =
                  sqrt(sum_error_portal->matrix[i][j][k])*mon_units;
            }
         }
      }
// --- End of Adding 16Nov2010 ---------------------------------

      // set total dose maximum
      plan.dose_max = dose_max;
      plan.d_max.x  = d_max.x;
      plan.d_max.y  = d_max.y;
      plan.d_max.z  = d_max.z;
   }
   else
   {
      // get total dose maximum
      dose_max = plan.dose_max;
      d_max.x  = plan.d_max.x;
      d_max.y  = plan.d_max.y;
      d_max.z  = plan.d_max.z;
   }

   // get cpu_time
   float cpu_time = plan.cpu_time;

   //
   // dose at the reference point
   // 

   // calculate voxel indices for x interpolation
   int lower_x = int((ref_point.x-voxel_size.x/TWO)/voxel_size.x);
   int upper_x = lower_x+1;
   if (lower_x <      0) { lower_x =       0; upper_x =       0; }
   if (upper_x >= dim.x) { lower_x = dim.x-1; upper_x = dim.x-1; }

   // interpolation factors
   real px = real(upper_x)+ONE_HALF-ref_point.x/voxel_size.x;
   real qx = ONE-px;

   // calculate voxel indices for y interpolation
   int lower_y = int((ref_point.y-voxel_size.y/TWO)/voxel_size.y);
   int upper_y = lower_y+1;
   if (lower_y <      0) { lower_y =       0; upper_y =       0; }
   if (upper_y >= dim.y) { lower_y = dim.y-1; upper_y = dim.y-1; }

   // interpolation factors
   real py = real(upper_y)+ONE_HALF-ref_point.y/voxel_size.y;
   real qy = ONE-py;

   // calculate voxel indices for z interpolation
   int lower_z = int((ref_point.z-voxel_size.z/TWO)/voxel_size.z);
   int upper_z = lower_z+1;
   if (lower_z <      0) { lower_z =       0; upper_z =       0; }
   if (upper_z >= dim.z) { lower_z = dim.z-1; upper_z = dim.z-1; }

   // interpolation factors
   real pz = real(upper_z)+ONE_HALF-ref_point.z/voxel_size.z;
   real qz = ONE-pz;

   // dose at the reference point
   float ref_dose  = sum_dose->matrix[lower_x][lower_y][lower_z]*px*py*pz
                   + sum_dose->matrix[lower_x][lower_y][upper_z]*px*py*qz
                   + sum_dose->matrix[lower_x][upper_y][lower_z]*px*qy*pz
                   + sum_dose->matrix[lower_x][upper_y][upper_z]*px*qy*qz
                   + sum_dose->matrix[upper_x][lower_y][lower_z]*qx*py*pz
                   + sum_dose->matrix[upper_x][lower_y][upper_z]*qx*py*qz
                   + sum_dose->matrix[upper_x][upper_y][lower_z]*qx*qy*pz
                   + sum_dose->matrix[upper_x][upper_y][upper_z]*qx*qy*qz;

   // error at the reference point
   float ref_error = sum_error->matrix[lower_x][lower_y][lower_z]*px*py*pz
                   + sum_error->matrix[lower_x][lower_y][upper_z]*px*py*qz
                   + sum_error->matrix[lower_x][upper_y][lower_z]*px*qy*pz
                   + sum_error->matrix[lower_x][upper_y][upper_z]*px*qy*qz
                   + sum_error->matrix[upper_x][lower_y][lower_z]*qx*py*pz
                   + sum_error->matrix[upper_x][lower_y][upper_z]*qx*py*qz
                   + sum_error->matrix[upper_x][upper_y][lower_z]*qx*qy*pz
                   + sum_error->matrix[upper_x][upper_y][upper_z]*qx*qy*qz;

   // output
   xvmc_message("Total CPU time:",cpu_time,"s",1);
#ifdef MONITORBSF   
   if (plan.result == ABS_DOSE || plan.result == ABS_DOSE_MBSF)
#else
   if (plan.result == ABS_DOSE)
#endif      
   {
      xvmc_message("Maximum Dose:  ",dose_max,"Gy",0);
      xvmc_message("Reference Dose:",ref_dose,"+/-",ref_error,"Gy",1);
      xvmc_message("               ",ref_error/ref_dose*100.0,"(%)",0);
   }
   else
   {
      xvmc_message("Maximum Dose:  ",dose_max,"10^-10 Gy cm^2",0);
      xvmc_message("Reference Dose:",ref_dose,"+/-",
                                     ref_error,"10^-10 Gy cm^2",1);
   }

   // length of file name strings for dose output
   int slength = strlen(inout_path)+10;

   // *******************
   // output 3D dose cube
   // *******************

   // write 3D total dose matrix if required
   if (plan.out3d_plan)
   {
      xvmc_message("Creating total 3D dose and error files",1);
      xvmc_message("======================================",0);

      // generate file names
      char *hed_name = new char[slength];    // dose matrix header file
      char *d3d_name = new char[slength];    // 3D dose data file
      char *e3d_name = new char[slength];    // 3D dose error file

      strcpy(hed_name,inout_path);
      strcpy(d3d_name,inout_path);
      strcpy(e3d_name,inout_path);

      strcat(hed_name,"00.hed");
      strcat(d3d_name,"00.d3d");
      strcat(e3d_name,"00.e3d");

      // save 3D data
      write_dose_file(sum_dose,sum_error,dose_max,d_max,
                      cpu_time,hed_name,d3d_name,e3d_name);

      // free memory
      delete [] hed_name;
      delete [] d3d_name;
      delete [] e3d_name;

      // write portal dose matrix for plan if required
      if (plan.portal != NULL)
      {
	 xvmc_message("Creating portal dose and error files for this beam",1);
	 xvmc_message("==================================================",0);

	 // generate file names and save data
	 int slength = strlen(inout_path)+10;   // length of file name strings
	 char *pdh_name = new char[slength];    // portal dose matrix header file
	 char *pdi_name = new char[slength];    // portal dose data file
	 char *pde_name = new char[slength];    // portal dose error file
	 char *bmp_name = new char[slength];    // portal dose BMP file

	 strcpy(pdh_name,inout_path);
	 strcpy(pdi_name,inout_path);
	 strcpy(pde_name,inout_path);
	 strcpy(bmp_name,inout_path);

	 strcat(pdh_name,"00.pdh");
         strcat(pdi_name,"00.pdi");
         strcat(pde_name,"00.pde");
         strcat(bmp_name,"00.bmp");

	 // write portal dose distribution for this beam to file
         float mon_units = plan.sum_monunits * plan.num_fractions;
	 //write_portal_file(pdh_name,pdi_name,
         //                  pde_name,bmp_name,
         //                  mon_units, plan.portal->distance);

	 write_portal_file(sum_dose_portal,sum_error_portal,
                         dose_max_portal,d_max_portal,plan.portal->distance,
                	 pdh_name,pdi_name,pde_name,bmp_name);

	 // free memory
	 delete [] pdh_name;
	 delete [] pdi_name;
	 delete [] pde_name;
	 delete [] bmp_name;
      }
   }

   // ********************************
   // output 2D dose profiles (planes)
   // ********************************

   // plane file name
   char *pln_name = NULL;

   // number of xy, xz and yz planes
   int num_pln_xy = 0;
   int num_pln_xz = 0;
   int num_pln_yz = 0;

   // output xy, xz and yz planes
   while (plan.first_plane != NULL)
   {
      plane_parameters *plane = plan.first_plane;

      // generate file name
      if (pln_name == NULL)
      {
         xvmc_message("Creating dose planes",1);
         xvmc_message("====================",0);
         pln_name = new char[slength];
      }
      strcpy(pln_name,inout_path);
      switch (plane->type)
      {
      case XY_PLANE:
         strcat(pln_name,".xy");
         char file_xyn[5];
         sprintf(file_xyn,"%u",num_pln_xy);
         strcat(pln_name,file_xyn);
         ++num_pln_xy;
         break;
      case XZ_PLANE:
         strcat(pln_name,".xz");
         char file_xzn[5];
         sprintf(file_xzn,"%u",num_pln_xz);
         strcat(pln_name,file_xzn);
         ++num_pln_xz;
         break;
      case YZ_PLANE:
         strcat(pln_name,".yz");
         char file_yzn[5];
         sprintf(file_yzn,"%u",num_pln_yz);
         strcat(pln_name,file_yzn);
         ++num_pln_yz;
         break;
      default:
         xvmc_error("evaluate(plan)","unknown plane type",8);
         break;
      }
      
      write_plane(sum_dose,sum_error,plane,pln_name);

      plan.first_plane = plane->next;
      if (plane == plan.last_plane) plan.last_plane = NULL;
      delete plane;
   }

   // free memory for dose plane file name
   if (pln_name != NULL) delete [] pln_name;

   // ***********************
   // output 1D dose profiles
   // ***********************

   // profile file name
   char *prf_name = NULL;

   // count x, y and z profiles
   int_3 num_prf;
   num_prf.x = 0;
   num_prf.y = 0;
   num_prf.z = 0;

   // output x, y and z profiles
   while (plan.first_profile != NULL)
   {
      profile_parameters *profile = plan.first_profile;

      // generate file name
      if (prf_name == NULL)
      {
         xvmc_message("Creating profiles",1);
         xvmc_message("=================",0);
         prf_name = new char[slength];
      }
      strcpy(prf_name,inout_path);
      switch (profile->type)
      {
      case X_PROFILE:
         strcat(prf_name,".px");
         char file_pxn[5];
         sprintf(file_pxn,"%u",num_prf.x);
         strcat(prf_name,file_pxn);
         ++num_prf.x;
         break;
      case Y_PROFILE:
         strcat(prf_name,".py");
         char file_pyn[5];
         sprintf(file_pyn,"%u",num_prf.y);
         strcat(prf_name,file_pyn);
         ++num_prf.y;
         break;
      case Z_PROFILE:
         strcat(prf_name,".pz");
         char file_pzn[5];
         sprintf(file_pzn,"%u",num_prf.z);
         strcat(prf_name,file_pzn);
         ++num_prf.z;
         break;
      default:
         xvmc_error("evaluate(plan)","unknown profile type",8);
         break;
      }
      
      write_profile(sum_dose,sum_error,profile,prf_name);

      plan.first_profile = profile->next;
      if (profile == plan.last_profile) plan.last_profile = NULL;
      delete profile;
   }

   // free memory for dose profile file name
   if (prf_name != NULL) delete [] prf_name;

   return;
}
