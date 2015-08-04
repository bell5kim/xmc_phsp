/*****************************************************************************
 * io_routines.cpp: functions for data input and output                      *
 *    function read_inp_file:     read XVMC input file                       *
 *    function read_density_file: read density matrix file                   *
 *    function read_dose_file:    read dose matrix file                      *
 *    function write_dose_file:   write dose matrix file                     *
 *    function write_plane:       write dose plane                           *
 *    function write_profile:     write dose profile                         *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 99/12/27        *
 *    changed class multi_leaf                            MF 06.11.2001      *
 *    added portal dose image                             MF 03/03/13        *
 *    check of valid block types and keywords             MF 03/07/17        *
 *    dynamic MLC mode implemented                        MF 22.10.2003      *
 *    physical jaw implemented                            JK 29.01.2008      *
 *    Change Line Length from 81 to 256                                      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
using namespace std;

#include <math.h>

#include "definitions.h"
#include "global.h"
#include "treatment_plan.h"

// ****************************************
// declare global variables
// ****************************************

void   swap_bytes( float &); // swap bytes
extern treatment_plan plan;  // treatment plan parameters

// read XVMC input file
void read_inp_file(const char *inp_name,
                   const char *hed_name,
                   const char *dmx_name)
{
   // open input file
   xvmc_message("Reading file:",inp_name,1);
   ifstream inp_file(inp_name,ios::in);
   if (!inp_file)
   {
      xvmc_error("read_inp_file","cannot open parameter file",8);
   }

   // set default values
   bool create_phantom = false;

   // read data cards
   char rline[256]    = "";
   char rtype         = '\0';
   char blocktype[17] = "";
   char rkeyword[17]  = "";
   char rvalue[64]    = "";
   bool read_cards = true;
   bool recognize_block = false;
   bool recognize_key   = false;
   while (read_cards)
   {
      inp_file.getline(rline,sizeof(rline));
      rtype = rline[0];
      if (rtype == '*')
      {
         for (register int i=0; i<16; ++i) blocktype[i] = rline[i+1];
         recognize_block = false;
         blocktype[16] = '\0';
         if (!strcmp(blocktype,"END-INPUT      |")) {
            recognize_block = true; read_cards     = false; }
         if (!strcmp(blocktype,"PHANTOM        |")) {
            recognize_block = true; create_phantom = true; }
         if (!strcmp(blocktype,"GLOBAL-DATA    |")) {
            recognize_block = true; plan.calculate = true; }
         if (!strcmp(blocktype,"BEAM-PARAMETERS|")) {
            recognize_block = true; plan.add_beam(0.0); }

         if (!recognize_block)
         {
            xvmc_warning("read_inp_file","block type not recognized",1);
            xvmc_message("unknown block type:",blocktype,1);
         }
      } // if (rtype == '*')

      if (rtype == '-')
      {
         for (register int i=0; i<16; ++i) rkeyword[i] = rline[i+1];
         rkeyword[16] = '\0';
         recognize_key = false;
         for (register int i=0; i<63; ++i) rvalue[i] = rline[i+17];
         rvalue[63] = '\0';
         istringstream values(rvalue);

         if (!strcmp(blocktype,"PHANTOM        |"))
         {

            if (!strcmp(rkeyword,"VOXELSIZE      |"))
            {
               recognize_key = true;
               values >> voxel_size.x >> voxel_size.y >> voxel_size.z;
               cube_size.x  = voxel_size.x * float(dim.x);
               cube_size.y  = voxel_size.y * float(dim.y);
               cube_size.z  = voxel_size.z * float(dim.z);
            }
            if (!strcmp(rkeyword,"DIMENSION      |"))
            {
               recognize_key = true;
               values >> dim.x >> dim.y >> dim.z;
               cube_size.x  = voxel_size.x * float(dim.x);
               cube_size.y  = voxel_size.y * float(dim.y);
               cube_size.z  = voxel_size.z * float(dim.z);
               if (density == NULL)
               {
                  if ( (density = new array_3d<float>(dim, ONE)) == NULL )
                  {
                     xvmc_error("read_inp_file",
                                "cannot create density->matrix",8);
                  }
               }
               else
               {
                  xvmc_error("read_inp_file",
                             "density->matrix already exists",8);
               }
            }
            if (!strcmp(rkeyword,"CHANGE-DENSITY |"))
            {
               recognize_key = true;
               int   il,iu,jl,ju,kl,ku;
               float rho;
               values >> il >> iu >> jl >> ju >> kl >> ku >> rho;
               if (density != NULL)
               {
                  for (register int i=il-1; i<iu; ++i)
                  {
                     for (register int j=jl-1; j<ju; ++j)
                     {
                        for (register int k=kl-1; k<ku; ++k)
                        {
                           density->matrix[i][j][k] = rho;
                        }
                     }
                  }
               }
               else
               {
                  xvmc_error("read_inp_file",
                             "density->matrix does not exists",8);
               }
            }

         } // if (!strcmp(blocktype,"PHANTOM        |"))

         if (!strcmp(blocktype,"GLOBAL-DATA    |"))
         {

            if (!strcmp(rkeyword,"E-CUTOFF       |"))
            {
               recognize_key = true;
               values >> e_cut;
               TWOe_cut = e_cut*TWO;
            }
            if (!strcmp(rkeyword,"ESTEPE         |"))
            {
               recognize_key = true;
               values >> e_step;
            }
            if (!strcmp(rkeyword,"P-CUTOFF       |"))
            {
               recognize_key = true;
               values >> p_cut;
            }
            if (!strcmp(rkeyword,"P-CUTOFF-KERMA |"))
            {
               recognize_key = true;
               values >> k1_cut;
               k0_cut = ZERO;
            }
            if (!strcmp(rkeyword,"P0-CUTOFF-KERMA|"))
            {
               recognize_key = true;
               values >> k0_cut;
            }
            if (!strcmp(rkeyword,"P1-CUTOFF-KERMA|"))
            {
               recognize_key = true;
               values >> k1_cut;
            }
            if (!strcmp(rkeyword,"WRITE-3D-DOSE  |"))
            {
               recognize_key = true;
               int write_3d_dose;
               values >> write_3d_dose;
               switch (write_3d_dose)
               {
               case -1:
                  plan.calculate  = false;
                  plan.out3d_beam = false;
                  plan.out3d_plan = false;
                  break;
               case 0:
                  plan.calculate  = true;
                  plan.out3d_beam = false;
                  plan.out3d_plan = false;
                  break;
               case 1:
                  plan.calculate  = true;
                  plan.out3d_beam = false;
                  plan.out3d_plan = true;
                  break;
               case 2:
                  plan.calculate  = true;
                  plan.out3d_beam = true;
                  plan.out3d_plan = true;
                  break;
               case 3:
                  plan.calculate  = true;
                  plan.out3d_beam = true;
                  plan.out3d_plan = false;
                  break;
               default:
                  xvmc_error("read_inp_file","unknown 3D dose output type",8);
                  break;
               }
            }
            if (!strcmp(rkeyword,"XY-PLANE       |"))
            {
               recognize_key = true;
               // create xy plane
               plane_parameters *plane = NULL;
               if ( (plane = new plane_parameters) == NULL )
               {
                  xvmc_error("read_inp_file",
                             "cannot create xy plane for output",8);
               }
               plane->type = XY_PLANE;   // xy plane
               plane->next = NULL;       // pointer to the next plane
               plane->pos  = ZERO;       // z coordinate of the plane
               values >> plane->pos;

               // assign pointers
               if (plan.first_plane == NULL)
               {
                  plan.first_plane = plane;
               }

               if (plan.last_plane == NULL)
               {
                  plan.last_plane = plane;
               }
               else
               {
                  plan.last_plane->next = plane;
                  plan.last_plane       = plane;
               }
            }
            if (!strcmp(rkeyword,"XZ-PLANE       |"))
            {
               recognize_key = true;
               // create xz plane
               plane_parameters *plane = NULL;
               if ( (plane = new plane_parameters) == NULL )
               {
                  xvmc_error("read_inp_file",
                             "cannot create xz plane for output",8);
               }
               plane->type = XZ_PLANE;   // xz plane
               plane->next = NULL;       // pointer to the next plane
               plane->pos  = ZERO;       // y coordinate of the plane
               values >> plane->pos;

               // assign pointers
               if (plan.first_plane == NULL)
               {
                  plan.first_plane = plane;
               }

               if (plan.last_plane == NULL)
               {
                  plan.last_plane = plane;
               }
               else
               {
                  plan.last_plane->next = plane;
                  plan.last_plane       = plane;
               }
            }
            if (!strcmp(rkeyword,"YZ-PLANE       |"))
            {
               recognize_key = true;
               // create yz plane
               plane_parameters *plane = NULL;
               if ( (plane = new plane_parameters) == NULL )
               {
                  xvmc_error("read_inp_file",
                             "cannot create yz plane for output",8);
               }
               plane->type = YZ_PLANE;   // yz plane
               plane->next = NULL;       // pointer to the next plane
               plane->pos  = ZERO;       // x coordinate of the plane
               values >> plane->pos;

               // assign pointers
               if (plan.first_plane == NULL)
               {
                  plan.first_plane = plane;
               }

               if (plan.last_plane == NULL)
               {
                  plan.last_plane = plane;
               }
               else
               {
                  plan.last_plane->next = plane;
                  plan.last_plane       = plane;
               }
            }
            if (!strcmp(rkeyword,"X-PROFILE      |"))
            {
               recognize_key = true;
               // create x profile
               profile_parameters *profile = NULL;
               if ( (profile = new profile_parameters) == NULL )
               {
                  xvmc_error("read_inp_file",
                             "cannot create x-profile",8);
               }
               profile->type  = X_PROFILE;  // x profile
               profile->next  = NULL;       // pointer to the next profile
               profile->pos.x = ZERO;       // x coordinate unnecessary
               values >> profile->pos.y >> profile->pos.z;

               // assign pointers
               if (plan.first_profile == NULL)
               {
                  plan.first_profile = profile;
               }

               if (plan.last_profile == NULL)
               {
                  plan.last_profile = profile;
               }
               else
               {
                  plan.last_profile->next = profile;
                  plan.last_profile       = profile;
               }
            }
            if (!strcmp(rkeyword,"Y-PROFILE      |"))
            {
               recognize_key = true;
               // create y profile
               profile_parameters *profile = NULL;
               if ( (profile = new profile_parameters) == NULL )
               {
                  xvmc_error("read_inp_file",
                             "cannot create y-profile",8);
               }
               profile->type  = Y_PROFILE;  // y profile
               profile->next  = NULL;       // pointer to the next profile
               profile->pos.y = ZERO;       // y coordinate unnecessary
               values >> profile->pos.x >> profile->pos.z;

               // assign pointers
               if (plan.first_profile == NULL)
               {
                  plan.first_profile = profile;
               }

               if (plan.last_profile == NULL)
               {
                  plan.last_profile = profile;
               }
               else
               {
                  plan.last_profile->next = profile;
                  plan.last_profile       = profile;
               }
            }
            if ( (!strcmp(rkeyword,"Z-PROFILE      |")) ||
                 (!strcmp(rkeyword,"DEPTH-DOSE     |"))    )
            {
               recognize_key = true;
               // create z profile (or depth dose)
               profile_parameters *profile = NULL;
               if ( (profile = new profile_parameters) == NULL )
               {
                  xvmc_error("read_inp_file",
                             "cannot create z-profile",8);
               }
               profile->type  = Z_PROFILE;  // z profile
               profile->next  = NULL;       // pointer to the next profile
               profile->pos.z = ZERO;       // z coordinate unnecessary
               values >> profile->pos.x >> profile->pos.y;

               // assign pointers
               if (plan.first_profile == NULL)
               {
                  plan.first_profile = profile;
               }

               if (plan.last_profile == NULL)
               {
                  plan.last_profile = profile;
               }
               else
               {
                  plan.last_profile->next = profile;
                  plan.last_profile       = profile;
               }
            }
            if (!strcmp(rkeyword,"PORTAL-DOSE    |"))
            {
               recognize_key = true;
               // define and read portal image parameters
               double       distance = ZERO;
               unsigned int dim_x    = 0;
               unsigned int dim_y    = 0;
               unsigned int dim_z    = 0;
               double       res_x    = ZERO;
               double       res_y    = ZERO;
               double       res_z    = ZERO;
               int          bmp_size = 0;
               values >> distance >> dim_x >> dim_y >> dim_z
                                  >> res_x >> res_y >> res_z
                                  >> bmp_size;
               
               // create portal dose image
               plan.portal = NULL;
               if ( (plan.portal = new portal_dose(distance,dim_x,dim_y,dim_z,
                                                   res_x,res_y,res_z,
                                                   bmp_size)) == NULL )
               {
                  xvmc_error("read_inp_file",
                             "cannot create portal dose image",8);
               }
               cube_size_portal.x  = voxel_size_portal.x * float(dim_portal.x);
               cube_size_portal.y  = voxel_size_portal.y * float(dim_portal.y);
               cube_size_portal.z  = voxel_size_portal.z * float(dim_portal.z);
               plan.portal->distance = distance;
               if (density_portal == NULL)
               {
                  if ( (density_portal = new array_3d<float>(dim_portal, ONE)) == NULL )
                  {
                     xvmc_error("read_inp_file",
                                "cannot create density_portal->matrix",8);
                  }
               }
               else
               {
                  xvmc_error("read_inp_file",
                             "density_portal->matrix already exists",8);
               }
            }
            if (!strcmp(rkeyword,"DOSE-TYPE      |"))
            {
               recognize_key = true;
               int dose_type = 0;
               values >> dose_type;
               switch (dose_type)
               {
               case 0:
                  plan.result = DOSE;
                  break;
               case 1:
                  plan.result = IONIZATION;
                  break;
               case 2:
                  plan.result = FILM;
                  values >> plan.film_factor;
                  break;
               case 3:
                  plan.result = DOSE_TO_H2O;
                  break;
               case 4:
                  plan.result = ABS_DOSE;
                  break;
#ifdef MONITORBSF		  
               case 5:  
	          // Same as ABS_DOSE plus Monitor Backscatter Factor Correction
                  plan.result = ABS_DOSE_MBSF;
                  break;
#endif		  
               default:
                  xvmc_error("read_inp_file","unknown dose type",8);
                  break;
               }
            }
            if (!strcmp(rkeyword,"NUM-FRACTIONS  |"))
            {
               recognize_key = true;
               values >> plan.num_fractions;
            }
            if (!strcmp(rkeyword,"REFERENCE-POINT|"))
            {
               recognize_key = true;
               values >> ref_point.x >> ref_point.y >> ref_point.z;
            }
            if (!strcmp(rkeyword,"RANDOM-SET     |"))
            {
               recognize_key = true;
               values >> plan.ini_rndm1 >> plan.ini_rndm2
                      >> plan.ini_rndm3 >> plan.ini_rndm4;
            }

         } // if (!strcmp(blocktype,"GLOBAL-DATA    |"))

         if (!strcmp(blocktype,"BEAM-PARAMETERS|"))
         {
            beam_core *this_beam = plan.last_beam;

            if (!strcmp(rkeyword,"BEAM-WEIGHT    |"))
            {
               recognize_key = true;
               real beam_weight; // new beam weight in % or # of MUs
               values >> beam_weight;
               // subtract old beam weight
               plan.sum_weight -= this_beam->weight;
               // set new beam weight
               this_beam->weight = beam_weight;
               // add new beam weight
               plan.sum_weight += this_beam->weight;
            }
            if (!strcmp(rkeyword,"DEVICE-TYPE    |"))
            {
               recognize_key = true;
               int device_type; // (old) VMC beam model numbers
               values >> device_type;
               if (device_type >= 100)
               {
                  // this is a photon beam
                  this_beam->type = PHOTON;

                  // determine beam model id
                  switch (device_type)
                  {
                  case 100:
                     // mono-energetic point source
                     this_beam->model_id = 0;
                     break;
                  case 101:
                     // poly-energetic point source (not implemented)
                     xvmc_error("read_inp_file",
                       "poly-energetic photon point source not implemented",8);
                     break;
                  case 102:
                     // mono-energetic beam with two Gaussian sources
                     this_beam->model_id = -1;
                     break;
                  case 103:
                     // poly-energetic beam with two Gaussian sources
                     this_beam->model_id = 1;
                     break;
                  case 121:
                     // point source photon beam, energy spectrum
                     this_beam->model_id = 21;
                     break;
                  case 131:
                     // TO DO:: phase space from phase space file
                     this_beam->model_id = 31;
                     break;
                  default:
                     xvmc_error("read_inp_file",
                                "this photon beam model is not implemented",8);
                     break;
                  }
               }
               else
               {
                  // this is an electron beam
                  this_beam->type = ELECTRON;

                  // determine beam model id
                  switch (device_type)
                  {
                  case 0:
                     // mono-energetic point source
                     this_beam->model_id = 0;
                     break;
                  case 1:
                     // poly-energetic point source
                     this_beam->model_id = 1;
                     break;
                  case 2:
                     // triple electron source model, mono-energetic
                     this_beam->model_id = -3;
                     break;
                  case 3:
                     // triple electron source model, poly-energetic
                     this_beam->model_id = 3;
                     break;
                  case -4:
                     // triple electron source model, mono-energetic
                     // with bump correction
                     this_beam->model_id = -4;
                     break;
                  case 4:
                     // triple electron source model, poly-energetic
                     // with bump correction
                     this_beam->model_id = 4;
                     break;
                  default:
                     xvmc_error("read_inp_file",
                            "this electron beam model is not implemented",8);
                     break;
                  }
               }
            }
            if (!strcmp(rkeyword,"DEVICE-KEY     |"))
            {
               recognize_key = true;
               char string[64];
               values >> string;
               int length = strlen(string)+1;
               if (this_beam->base_key == NULL)
               {
                  if ( (this_beam->base_key = new char[length]) == NULL )
                  {
                     xvmc_error("read_inp_file",
                                "cannot create base key",8);
                  }
                  strcpy(this_beam->base_key,string);
               }
               else
               {
                  xvmc_error("read_inp_file",
                             "base key already defined",8);
               }
            }
            if (!strcmp(rkeyword,"APPLICATOR     |"))
            {
               recognize_key = true;
               char string[64];
               values >> string;
               int length = strlen(string)+1;
               if (this_beam->applicator == NULL)
               {
                  if ( (this_beam->applicator = new char[length]) == NULL )
                  {
                     xvmc_error("read_inp_file",
                                "cannot create applicator",8);
                  }
                  strcpy(this_beam->applicator,string);
               }
               else
               {
                  xvmc_error("read_inp_file",
                             "applicator already defined",8);
               }
            }
            if (!strcmp(rkeyword,"EVENT-NUMBER   |"))
            {
               recognize_key = true;
               values >> this_beam->n_history
                      >> this_beam->n_repeat
                      >> this_beam->n_rotate
                      >> this_beam->n_batch;
               // "n_batch" must be an integer multiple of "n_process"
               int n_test = (this_beam->n_batch-1)/n_process + 1;
               int new_batch = n_test*n_process;
               if (this_beam->n_batch != new_batch)
               {
                  this_beam->n_batch = new_batch;
                  xvmc_warning("read_inp_file",
                               "number of batches changed",1);
               }
            }
            if (!strcmp(rkeyword,"NOMINAL-ENERGY |"))
            {
               recognize_key = true;
               values >> this_beam->energy;
            }
            if (!strcmp(rkeyword,"ISOCENTER      |"))
            {
               recognize_key = true;
               values >> this_beam->iso_center.x
                      >> this_beam->iso_center.y
                      >> this_beam->iso_center.z;
            }
            if (!strcmp(rkeyword,"GANTRY-ANGLE   |"))
            {
               recognize_key = true;
               char torf;
               values >> torf;
               if (torf == 'T')
               {
                  this_beam->gantry_rotation = true;
                  values >> this_beam->start_gantry_angle
                         >> this_beam->stop_gantry_angle;
               }
               else
               {
                  this_beam->gantry_rotation = false;
                  values >> this_beam->start_gantry_angle;
               }
            }
            if (!strcmp(rkeyword,"TABLE-ANGLE    |"))
            {
               recognize_key = true;
               values >> this_beam->table_angle;
            }
            if (!strcmp(rkeyword,"COLL-ANGLE     |"))
            {
               recognize_key = true;
               values >> this_beam->collimator_angle;
            }
            if (!strcmp(rkeyword,"COLL-WIDTH-X   |"))
            {
               recognize_key = true;
               real collimator_width_x;
               values >> collimator_width_x;
               this_beam->open_x1 = -collimator_width_x/2.0;
               this_beam->open_x2 =  collimator_width_x/2.0;
            }
            if (!strcmp(rkeyword,"COLL-WIDTH-Y   |"))
            {
               recognize_key = true;
               real collimator_width_y;
               values >> collimator_width_y;
               this_beam->open_y1 = -collimator_width_y/2.0;
               this_beam->open_y2 =  collimator_width_y/2.0;
            }
            if (!strcmp(rkeyword,"COLL-LEFT-X    |"))
            {
               recognize_key = true;
               values >> this_beam->open_x1;
            }
            if (!strcmp(rkeyword,"COLL-RIGHT-X   |"))
            {
               recognize_key = true;
               values >> this_beam->open_x2;
            }
            if (!strcmp(rkeyword,"COLL-LEFT-Y    |"))
            {
               recognize_key = true;
               values >> this_beam->open_y1;
            }
            if (!strcmp(rkeyword,"COLL-RIGHT-Y   |"))
            {
               recognize_key = true;
               values >> this_beam->open_y2;
            }
            if (!strcmp(rkeyword,"IRREGULAR-FIELD|"))
            {
               recognize_key = true;
               int num_points;
               values >> num_points;
               if ( this_beam->irregular == NULL)
               {
                  if ((this_beam->irregular=new contour(num_points))==NULL)
                  {
                     xvmc_error("read_inp_file",
                                "cannot create irregular beam opening",8);
                  }
               }
               else
               {
                  xvmc_error("read_inp_file",
                             "irregular beam opening already defined",8);
               }
               // input contour points
               for (register int i=0; i<num_points; ++i)
               {
                  char pline[256] = "";
                  inp_file.getline(pline,sizeof(pline));
                  istringstream points(pline);
                  // read point (x,y) and put into contour
                  real x,y;
                  points >> x >> y;
                  if (!this_beam->irregular->change_point(i,x,y))
                  {
                     xvmc_error("read_inp_file",
                                "cannot change irregular beam contour point",8);
                  }
               }
               // close contour
               if (!this_beam->irregular->close())
               {
                  xvmc_error("read_inp_file",
                             "cannot close irregular beam contour",8);
               }
            }
            if ( (!strcmp(rkeyword,"SIMPLE-MLC     |")) ||
                 (!strcmp(rkeyword,"SIMPLE-DMLC    |")) ||
                 (!strcmp(rkeyword,"DBLFOCUS-MLC   |")) ||
                 (!strcmp(rkeyword,"DBLFOCUS-DMLC  |")) ||
                 (!strcmp(rkeyword,"RNDFOCUS-MLC   |")) ||
                 (!strcmp(rkeyword,"RNDFOCUS-DMLC  |")) ||
                 (!strcmp(rkeyword,"ELEKTA-MLC     |")) ||
                 (!strcmp(rkeyword,"ELEKTA-DMLC    |")) ||
                 (!strcmp(rkeyword,"VARIAN-MLC     |")) ||
                 (!strcmp(rkeyword,"VARIAN-DMLC    |"))    )
            {
               recognize_key = true;

               // determine MLC type and mode
               mlc_type type = SIMPLE_MLC;
               mlc_mode mode = STATIC_MLC;
               if   (!strcmp(rkeyword,"SIMPLE-MLC     |"))
               {
                  type = SIMPLE_MLC;
                  mode = STATIC_MLC;
               }
               if   (!strcmp(rkeyword,"SIMPLE-DMLC    |"))
               {
                  type = SIMPLE_MLC;
                  mode = DYNAMIC_MLC;
               }
               if   (!strcmp(rkeyword,"DBLFOCUS-MLC   |"))
               {
                  type = DBLFOCUS_MLC;
                  mode = STATIC_MLC;
               }
               if   (!strcmp(rkeyword,"DBLFOCUS-DMLC  |"))
               {
                  type = DBLFOCUS_MLC;
                  mode = DYNAMIC_MLC;
               }
               if   (!strcmp(rkeyword,"RNDFOCUS-MLC   |"))
               {
                  type = RNDFOCUS_MLC;
                  mode = STATIC_MLC;
               }
               if   (!strcmp(rkeyword,"RNDFOCUS-DMLC  |"))
               {
                  type = RNDFOCUS_MLC;
                  mode = DYNAMIC_MLC;
               }
               if   (!strcmp(rkeyword,"ELEKTA-MLC     |"))
               {
                  type = ELEKTA_MLC;
                  mode = STATIC_MLC;
               }
               if   (!strcmp(rkeyword,"ELEKTA-DMLC    |"))
               {
                  type = ELEKTA_MLC;
                  mode = DYNAMIC_MLC;
               }
               if   (!strcmp(rkeyword,"VARIAN-MLC     |"))
               {
                  type = VARIAN_MLC;
                  mode = STATIC_MLC;
               }
               if   (!strcmp(rkeyword,"VARIAN-DMLC    |"))
               {
                  type = VARIAN_MLC;
                  mode = DYNAMIC_MLC;
               }

               int  num_pairs;       // number of leaf pairs
               char mlc_xy;          // leaf moving direction X or Y
               char material[40];    // MLC material name
               real z_upper,z_lower; // upper and lower leaf limits
               real z_curve,r_curve; // parameters for leaf curvature
               values >> num_pairs >> mlc_xy   >> material
                      >> z_upper   >> z_lower >> z_curve >> r_curve;
               if ( this_beam->mlc == NULL)
               {
                  if ( (this_beam->mlc=new multi_leaf(type,mode,
                              mlc_xy,num_pairs,material,
                              &z_upper,&z_lower,&z_curve,&r_curve)) == NULL )
                  {
                     xvmc_error("read_inp_file",
                                "cannot create multi-leaf collimator",8);
                  }
               }
               else
               {
                  xvmc_error("read_inp_file",
                             "multi-leaf collimator already defined",8);
               }
               // input leaf positions
               for (register int i=0; i<num_pairs; ++i)
               {
                  char pline[256] = "";
                  inp_file.getline(pline,sizeof(pline));
                  istringstream pairs(pline);

                  if (mode == STATIC_MLC)
                  {
                     // static MLC mode
                     // read leaf width, left and rigth positions, put into MLC
                     real width,left,right;  // defined at the iso-center plane
                     pairs >> width >> left >> right;
                     if (!this_beam->mlc->change_width(i,width))
                     {
                        xvmc_error("read_inp_file",
                                   "cannot change leaf width",8);
                     }
                     if (!this_beam->mlc->change_pair(i,left,right))
                     {
                        xvmc_error("read_inp_file",
                                   "cannot change leaf positions",8);
                     }
                     if (!this_beam->mlc->change_start_pair(i,left,right))
                     {
                        xvmc_error("read_inp_file",
                                   "cannot change leaf starting positions",8);
                     }
                     if (!this_beam->mlc->change_stop_pair(i,left,right))
                     {
                        xvmc_error("read_inp_file",
                                   "cannot change leaf stopping positions",8);
                     }
                  }
                  else
                  {
                     // dynamic MLC mode
                     // read leaf width, left and rigth starting and stopping
                     // positions defined at the iso-center plane, put into MLC
                     real width,start_left,start_right,stop_left,stop_right;
                     pairs >> width >> start_left >> start_right
                                    >> stop_left  >> stop_right;
                     if (!this_beam->mlc->change_width(i,width))
                     {
                        xvmc_error("read_inp_file",
                                   "cannot change leaf width",8);
                     }
                     if (!this_beam->mlc->change_pair(i,start_left,
                                                        start_right))
                     {
                        xvmc_error("read_inp_file",
                                   "cannot change leaf positions",8);
                     }
                     if (!this_beam->mlc->change_start_pair(i,start_left,
                                                              start_right))
                     {
                        xvmc_error("read_inp_file",
                                   "cannot change leaf starting positions",8);
                     }
                     if (!this_beam->mlc->change_stop_pair(i,stop_left,
                                                             stop_right))
                     {
                        xvmc_error("read_inp_file",
                                   "cannot change leaf stopping positions",8);
                     }
                  }
               }
            }
/*
            if ( (!strcmp(rkeyword,"SIMPLE-JAW     |")) ||
                 (!strcmp(rkeyword,"SIMPLE-DJAW    |")) ||
                 (!strcmp(rkeyword,"DBLFOCUS-JAW   |")) ||
                 (!strcmp(rkeyword,"DBLFOCUS-DJAW  |")) ||
                 (!strcmp(rkeyword,"REAL-JAW       |")) ||
                 (!strcmp(rkeyword,"REAL-DJAW      |"))    )
            {
               recognize_key = true;

               // determine MLC type and mode
               jaw_type type = SIMPLE_JAW;
               jaw_mode mode = STATIC_JAW;
               if   (!strcmp(rkeyword,"SIMPLE-JAW     |"))
               {
                  type = SIMPLE_JAW;
                  mode = STATIC_JAW;
               }
               if   (!strcmp(rkeyword,"SIMPLE-DJAW    |"))
               {
                  type = SIMPLE_JAW;
                  mode = DYNAMIC_JAW;
               }
               if   (!strcmp(rkeyword,"DBLFOCUS-JAW   |"))
               {
                  type = DBLFOCUS_JAW;
                  mode = STATIC_JAW;
               }
               if   (!strcmp(rkeyword,"DBLFOCUS-DJAW  |"))
               {
                  type = DBLFOCUS_JAW;
                  mode = DYNAMIC_JAW;
               }
               if   (!strcmp(rkeyword,"REAL-MLC       |"))
               {
                  type = REAL_JAW;
                  mode = STATIC_JAW;
               }
               if   (!strcmp(rkeyword,"REAL-DJAW      |"))
               {
                  type = REAL_JAW;
                  mode = DYNAMIC_JAW;
               }

               char jaw_xy;          // jaw moving direction X or Y
               char material[40];    // MLC material name
               real open_left,open_right; // open_left and open_right
               real z_upper,z_lower; // upper and lower jaw limits
               values >>  jaw_xy   >> material >> z_upper   >> z_lower;

               axis x_jaw = X;
               // check jaw moving direction (X or Y)
               switch (jaw_xy)
               {
               case 'x':    // X jaw
               case 'X':    // X jaw
                  x_jaw = X;
                  break;
               case 'y':    // Y jaw
               case 'Y':    // Y jaw
                  x_jaw = Y;
                  break;
               default:
                  xvmc_error("io_routine::read_inp_file","there are only X or Y JAWs",8);
                  break;
               }

               if ( this_beam->jaw == NULL)
               {
                  if ( (this_beam->jaw=new MC_jaw(x_jaw,type,mode,
                              open_left,open_right,
                              z_upper,z_lower,material,"air")) == NULL )
                  {
                     xvmc_error("read_inp_file",
                                "cannot create physical collimator",8);
                  }
               }
               else
               {
                  xvmc_error("read_inp_file",
                             "physical collimator already defined",8);
               }
            } // if ( (!strcmp(rkeyword,"SIMPLE-JAW     |")) ||
*/
         } // if (!strcmp(blocktype,"BEAM-PARAMETERS|"))

         if ( (!strcmp(blocktype,"PHANTOM        |")) ||
              (!strcmp(blocktype,"GLOBAL-DATA    |")) ||
              (!strcmp(blocktype,"BEAM-PARAMETERS|"))    )
         {
            if (!recognize_key)
            {
               xvmc_warning("read_inp_file","keyword not recognized",1);
               xvmc_message("unknown keyword:",rkeyword,1);
            }
         }

      } // if (rtype == '-')

   } // while (read_cards)

   // close input file
   inp_file.close();

   // normalize beam weights
   plan.normalize_weights();

   // print plan parameters
   xvmc_message("patient_id:        ",plan.patient_id,1);
#ifndef NO_PHANTOM_INFO
	xvmc_message("phantom x dim:     ",dim.x,"",0);
	xvmc_message("phantom y dim:     ",dim.y,"",0);
	xvmc_message("phantom z dim:     ",dim.z,"",0);
	xvmc_message("voxel x size:      ",voxel_size.x,"cm",0);
	xvmc_message("voxel y size:      ",voxel_size.y,"cm",0);
	xvmc_message("voxel z size:      ",voxel_size.z,"cm",0);
#endif
   xvmc_message("plan_id:           ",plan.plan_id,0);
   xvmc_message("number of beams:   ",plan.total_number,"",0);
   xvmc_message("sum of weights:    ",plan.sum_weight,"",0);
   xvmc_message("number of MUs:     ",plan.sum_monunits,"",0);
   xvmc_message("result type:       ",int(plan.result),"",0);
   xvmc_message("film_factor:       ",plan.film_factor,"",0);

   // print beam parameters
   beam_core *this_beam = plan.first_beam;
   while (this_beam != NULL)
   {
      xvmc_message("id:                ",this_beam->id,"",1);
      xvmc_message("type:              ",int(this_beam->type),"",0);
      xvmc_message("model_id:          ",this_beam->model_id,"",0);
      xvmc_message("weight:            ",this_beam->weight,"",0);
      if (this_beam->base_key != NULL)
         xvmc_message("base_key:          ",this_beam->base_key,0);
      if (this_beam->applicator != NULL)
         xvmc_message("applicator:        ",this_beam->applicator,0);
      xvmc_message("n_history:         ",this_beam->n_history,"",0);
      xvmc_message("n_repeat:          ",this_beam->n_repeat,"",0);
      xvmc_message("n_rotate:          ",this_beam->n_rotate,"",0);
      xvmc_message("n_batch:           ",this_beam->n_batch,"",0);
      xvmc_message("energy:            ",this_beam->energy,"MeV",0);
      xvmc_message("iso_center.x:      ",this_beam->iso_center.x,"cm",0);
      xvmc_message("iso_center.y:      ",this_beam->iso_center.y,"cm",0);
      xvmc_message("iso_center.z:      ",this_beam->iso_center.z,"cm",0);
      xvmc_message("gantry_rotation:   ",
                   int(this_beam->gantry_rotation),"",0);
      xvmc_message("start_gantry_angle:",
                   this_beam->start_gantry_angle,"deg",0);
      xvmc_message("stop_gantry_angle: ",
                   this_beam->stop_gantry_angle,"deg",0);
      xvmc_message("table_angle:       ",this_beam->table_angle,"deg",0);
      xvmc_message("collimator_angle:  ",
                   this_beam->collimator_angle,"deg",0);
      xvmc_message("open_x1:           ",this_beam->open_x1,"cm",0);
      xvmc_message("open_x2:           ",this_beam->open_x2,"cm",0);
      xvmc_message("open_y1:           ",this_beam->open_y1,"cm",0);
      xvmc_message("open_y2:           ",this_beam->open_y2,"cm",0);
      if (this_beam->irregular != NULL)
      {
         int num_points = this_beam->irregular->get_num();
         xvmc_message("== start of irregular beam opening ==",0);
         xvmc_message("  number of points:     ",num_points,"",0);
         for (register int i=0; i<num_points; ++i)
         {
            xvmc_message(
               "  x_con[i]:",this_beam->irregular->get_x(i),
               "  y_con[i]:",this_beam->irregular->get_y(i),"",0);
         }
         xvmc_message("== end of irregular beam opening ==",0);
      }

      // MLC type
      if (this_beam->mlc != NULL)
      {
         switch (this_beam->mlc->get_type())
         {
         case DBLFOCUS_MLC:
            xvmc_message("== double focussing MLC ==",0);
            break;
         case RNDFOCUS_MLC:
            xvmc_message("== focussing MLC with curved leaf ends ==",0);
            break;
         case ELEKTA_MLC:
            xvmc_message("== ELEKTA MLC ==",0);
            break;
         default:
            xvmc_message("== simple MLC ==",0);
            break;
         }

         // MLC mode
         switch (this_beam->mlc->get_mode())
         {
         case DYNAMIC_MLC:
            xvmc_message("  MLC mode:              ","dynamic",0);
            break;
         default:
            xvmc_message("  MLC mode:              ","static",0);
            break;
         }

         int num_pairs  = this_beam->mlc->get_num();
         char mlc_xy[2] = " ";
         mlc_xy[0] = this_beam->mlc->get_xy();
         xvmc_message("  leaf moving direction: ",mlc_xy,0);
         char *material = this_beam->mlc->get_material();
         if (material != NULL)
         {
            xvmc_message("  MLC material:          ",material,0);
         }
         real *z_upper = this_beam->mlc->get_upper();
         if (z_upper != NULL)
         {
            xvmc_message("  upper MLC limit:       ",*z_upper,"",0);
         }
         real *z_lower = this_beam->mlc->get_lower();
         if (z_lower != NULL)
         {
            xvmc_message("  lower MLC limit:       ",*z_lower,"",0);
         }
         real *z_curve = this_beam->mlc->get_z_center_curve();
         if (z_curve != NULL)
         {
            xvmc_message("  leaf curve z-position: ",*z_curve,"",0);
         }
         real *r_curve = this_beam->mlc->get_radius_curve();
         if (r_curve != NULL)
         {
            xvmc_message("  leaf curve radius:     ",*r_curve,"",0);
         }
         xvmc_message("  number of leaf pairs:  ",num_pairs,"",0);
         for (register int i=0; i<num_pairs; ++i)
         {
            real lperp =
             (this_beam->mlc->get_lperp(i+1)+this_beam->mlc->get_lperp(i))/TWO;
            xvmc_message(
               "  leaf:    ",real(i),
               "  width:   ",this_beam->mlc->get_width(i),
               "  position:",lperp,"",0);
         }
         if (this_beam->mlc->get_mode() == DYNAMIC_MLC)
         {
            for (register int k=0; k<num_pairs; ++k)
            {
               xvmc_message(
                  "  leaf:   ",real(k),
                  "  start left:   ",this_beam->mlc->get_start_left(k),
                  "  stop  left:   ",this_beam->mlc->get_stop_left(k),"",0);

            }
            for (register int k=0; k<num_pairs; ++k)
            {
               xvmc_message(
                  "  leaf:   ",real(k),
                  "  start right: ",this_beam->mlc->get_start_right(k),
                  "  stop rigth:  ",this_beam->mlc->get_stop_right(k),"",0);
            }
         }
         else
         {
            for (register int k=0; k<num_pairs; ++k)
            {
               xvmc_message(
                  "  leaf:   ",real(k),
                  "  left:   ",this_beam->mlc->get_left(k),
                  "  rigth:  ",this_beam->mlc->get_right(k),"",0);
            }
         }
         xvmc_message("== end of multi-leaf collimator ==",0);
      }
      this_beam = this_beam->next();
   }

   // write phantom density matrix to file
   if (create_phantom)
   {
      float *out_data = NULL;     // pointer to density output data

      // allocate memory for density output data
      if ( (out_data = new float[dim.x*dim.y*dim.z]) == NULL )
      {
         xvmc_error("read_inp_file",
                    "cannot allocate memory for density output data",8);
      }

      // write density data to output memory
      float rho;
      float *out_count = out_data;
      for (register int k=0; k<dim.z; ++k)
      {
         for (register int j=0; j<dim.y; ++j)
         {
            for (register int i=0; i<dim.x; ++i)
            {
               rho = density->matrix[i][j][k];
               swap_bytes(rho);
               *out_count = rho;
               ++out_count;
            }
         }
      }

      // write header file
      xvmc_message("Writing file:",hed_name,1);
      ofstream hed_file(hed_name,ios::out);
      if (!hed_file)
      {
         xvmc_error("read_inp_file","cannot open density header file",8);
      }
      hed_file.setf( ios::fixed, ios::floatfield );

      hed_file << setprecision(6);
      hed_file << "VOXELSIZE      |" << setw(11) << voxel_size.x
                                     << setw(11) << voxel_size.y
                                     << setw(11) << voxel_size.z << endl;
      hed_file << "DIMENSION      |" << setw(11) << dim.x
                                     << setw(11) << dim.y
                                     << setw(11) << dim.z << endl;
      hed_file << "END-INPUT      |" << endl;

      hed_file.close();

      // write density data file
      xvmc_message("Writing file:",dmx_name,0);
#ifdef OSF1
      FILE *dmx_file=fopen( dmx_name,"wb" );
      fwrite( out_data, sizeof(float), dim.x*dim.y*dim.z, dmx_file );
      fclose( dmx_file );
#else
      ofstream dmx_file(dmx_name,ios::out|ios::binary);
      char *out_dat_memory = (char *) out_data;
      dmx_file.write( out_dat_memory, dim.x*dim.y*dim.z*sizeof(float));
      dmx_file.close();
#endif

      // free memory
      delete [] out_data; out_data = NULL;
   }
}

// read density matrix file
void read_density_file(const char *hed_name, const char *dmx_name)
{
   // set initial values
   voxel_size.x = ZERO;
   voxel_size.y = ZERO;
   voxel_size.z = ZERO;
   dim.x        = 0;
   dim.y        = 0;
   dim.z        = 0;
   cube_size.x  = voxel_size.x * float(dim.x);
   cube_size.y  = voxel_size.y * float(dim.y);
   cube_size.z  = voxel_size.z * float(dim.z);

   // open header file
   xvmc_message("Reading file:",hed_name,1);
   ifstream hed_file(hed_name,ios::in);
   if (!hed_file)
   {
      xvmc_error("read_density_file","cannot open density header file",8);
   }

   // read data cards
   char keyword[17],bar;
   bool read_cards = true;
   while (read_cards)
   {
      hed_file >> keyword;

      if (!strcmp(keyword,"END-INPUT"))
      {
         read_cards = false;
      }
      else
      {
         if (!strcmp(keyword,"VOXELSIZE"))
         {
            hed_file >> bar >> voxel_size.x >> voxel_size.y >> voxel_size.z;
         }
         if (!strcmp(keyword,"DIMENSION"))
         {
            hed_file >> bar >> dim.x >> dim.y >> dim.z;
         }
      }
   }

   // close header file
   hed_file.close();

   // check voxel size
   if ( (voxel_size.x == ZERO) ||
        (voxel_size.y == ZERO) ||
        (voxel_size.z == ZERO) )
   {
      xvmc_error("read_density_file","zero voxel size",8);
   }

   // check cube dimensions
   if ( (dim.x == 0) || (dim.y == 0) || (dim.z == 0) )
   {
      xvmc_error("read_density_file","zero cube dimension",8);
   }

   // set cube size
   cube_size.x  = voxel_size.x * float(dim.x);
   cube_size.y  = voxel_size.y * float(dim.y);
   cube_size.z  = voxel_size.z * float(dim.z);

   // allocate memory for density input data
   float *inp_data = NULL;     // pointer to density input data
   if ( (inp_data = new float[dim.x*dim.y*dim.z]) == NULL )
   {
      xvmc_error("read_density_file",
                 "cannot allocate memory for density input data",8);
   }

   // read density data file
   xvmc_message("Reading file:",dmx_name,0);
#ifdef OSF1
   FILE *dmx_file=fopen( dmx_name,"rb" );
   fread( inp_data, sizeof(float), dim.x*dim.y*dim.z, dmx_file );
   fclose( dmx_file );
#else
   ifstream dmx_file(dmx_name,ios::in|ios::binary);
   char *inp_dat_memory = (char *) inp_data;
   dmx_file.read( inp_dat_memory, dim.x*dim.y*dim.z*sizeof(float));
   dmx_file.close();
#endif

   // create density matrix with rho = 1.0 g/cm^3
   if ( (density = new array_3d<float>(dim, ONE)) == NULL ) {
      xvmc_error("read_density_file","cannot construct density->matrix",8); }

   // write density data to output memory
   float rho;
   float *inp_count = inp_data;
   for (register int k=0; k<dim.z; ++k)
   {
      for (register int j=0; j<dim.y; ++j)
      {
         for (register int i=0; i<dim.x; ++i)
         {
            rho = *inp_count;
            swap_bytes(rho);
            density->matrix[i][j][k] = rho;
            ++inp_count;
         }
      }
   }

   // free memory
   delete [] inp_data; inp_data = NULL;
}

// read dose matrix file
void read_dose_file(const array_3d<float>  *inp_dose,
                    const array_3d<float>  *inp_error,
                    float &dose_max, int_3 &d_max, float &cpu_time,
                    const char *hed_name,
                    const char *dat_name, const char *err_name)
{
   // define input parameters
   real_3 inp_voxel_size;
   int_3  inp_dim;

   // initialize input parameters
   inp_voxel_size.x = ZERO;
   inp_voxel_size.y = ZERO;
   inp_voxel_size.z = ZERO;
   inp_dim.x        = 0;
   inp_dim.y        = 0;
   inp_dim.z        = 0;
   dose_max         = ZERO;
   d_max.x          = 0;
   d_max.y          = 0;
   d_max.z          = 0;
   cpu_time         = ZERO;

   // open header file
   xvmc_message("Reading file:",hed_name,1);
   ifstream hed_file(hed_name,ios::in);
   if (!hed_file)
   {
      xvmc_error("read_dose_file","cannot open dose header file",8);
   }

   // read data cards
   char keyword[17],bar;
   bool read_cards = true;
   while (read_cards)
   {
      hed_file >> keyword;

      if (!strcmp(keyword,"END-INPUT"))
      {
         read_cards = false;
      }
      else
      {
         if (!strcmp(keyword,"VOXELSIZE"))
         {
            hed_file >> bar >> inp_voxel_size.x
                            >> inp_voxel_size.y
                            >> inp_voxel_size.z;
         }
         if (!strcmp(keyword,"DIMENSION"))
         {
            hed_file >> bar >> inp_dim.x >> inp_dim.y >> inp_dim.z;
         }
         if (!strcmp(keyword,"DOSE-MAXIMUM"))
         {
            hed_file >> bar >> dose_max;
         }
         if (!strcmp(keyword,"POS-MAXIMUM"))
         {
            hed_file >> bar >> d_max.x >> d_max.y >> d_max.z;
         }
         if (!strcmp(keyword,"CPU-TIME"))
         {
            hed_file >> bar >> cpu_time;
         }
      }
   }

   // close header file
   hed_file.close();

   // check voxel size
   if ( (inp_voxel_size.x == ZERO) ||
        (inp_voxel_size.y == ZERO) ||
        (inp_voxel_size.z == ZERO) )
   {
      xvmc_error("read_dose_file","zero voxel size",8);
   }
   if ( (inp_voxel_size.x != voxel_size.x) ||
        (inp_voxel_size.y != voxel_size.y) ||
        (inp_voxel_size.z != voxel_size.z)    )
   {
      xvmc_error("read_dose_file","inconsistent voxel sizes",8);
   }

   // check cube dimensions
   if ( (inp_dim.x == 0) || (inp_dim.y == 0) || (inp_dim.z == 0) )
   {
      xvmc_error("read_dose_file","zero cube dimension",8);
   }
   if ( (inp_dim.x != dim.x) || (inp_dim.y != dim.y) || (inp_dim.z != dim.z) )
   {
      xvmc_error("read_dose_file","inconsistent voxel dimensions",8);
   }

   // allocate memory for dose input data
   float *inp_data = NULL;     // pointer to dose input data
   if ( (inp_data = new float[dim.x*dim.y*dim.z]) == NULL )
   {
      xvmc_error("read_dose_file",
                 "cannot allocate memory for input data",8);
   }

   // read dose data file
   xvmc_message("Reading file:",dat_name,0);
#ifdef OSF1
   FILE *dat_file=fopen( dat_name,"rb" );
   fread( inp_data, sizeof(float), dim.x*dim.y*dim.z, dat_file );
   fclose( dat_file );
#else
   ifstream dat_file(dat_name,ios::in|ios::binary);
   if (!dat_file)
   {
      xvmc_error("read_dose_file","cannot open dose data file",8);
   }
   char *inp_dat_memory = (char *) inp_data;
   dat_file.read( inp_dat_memory, dim.x*dim.y*dim.z*sizeof(float));
   dat_file.close();
#endif

   // write dose data to input dose cube
   float  dose;
   float *inp_count = inp_data;
   for (register int k=0; k<dim.z; ++k)
   {
      for (register int j=0; j<dim.y; ++j)
      {
         for (register int i=0; i<dim.x; ++i)
         {
            dose = *inp_count;
            swap_bytes(dose);
            inp_dose->matrix[i][j][k] = dose;
            ++inp_count;
         }
      }
   }

   // read dose error file
   xvmc_message("Reading file:",err_name,0);
   bool err_file_exist = false;
#ifdef OSF1
   FILE *err_file=fopen( err_name,"rb" );
   fread( inp_data, sizeof(float), dim.x*dim.y*dim.z, err_file );
   fclose( err_file );
#else
   ifstream err_file(err_name,ios::in|ios::binary);
   if (err_file)
   {
      err_file_exist = true;
      char *inp_err_memory = (char *) inp_data;
      err_file.read( inp_err_memory, dim.x*dim.y*dim.z*sizeof(float));
      err_file.close();
   }
   else
   {
      xvmc_warning("read_dose_file","dose error file does not exist",1);
   }
#endif

   // if error file exist: write dose error data to input dose error cube
   float error;
   inp_count = inp_data;
   for (register int k=0; k<dim.z; ++k)
   {
      for (register int j=0; j<dim.y; ++j)
      {
         for (register int i=0; i<dim.x; ++i)
         {
            if (err_file_exist)
            {
               error = *inp_count;
               swap_bytes(error);
               inp_error->matrix[i][j][k] = error;
            }
            else
            {
               inp_error->matrix[i][j][k] = ZERO;
            }
            ++inp_count;
         }
      }
   }

   // free memory and return
   delete [] inp_data; inp_data = NULL;
   return;
}

// write dose matrix file
void write_dose_file(const array_3d<float>  *out_dose,
                     const array_3d<float>  *out_error,
                     const float dose_max, const int_3 d_max,
                     const float cpu_time, const char *hed_name,
                     const char *dat_name, const char *err_name)
{
   float *out_data = NULL;     // pointer to dose output data
   float *err_data = NULL;     // pointer to dose error output data

   // allocate memory for dose output data
   if ( (out_data = new float[dim.x*dim.y*dim.z]) == NULL )
   {
      xvmc_error("write_dose_file",
                 "cannot allocate memory for dose output data",8);
   }

   // allocate memory for dose error output data
   if ( (err_data = new float[dim.x*dim.y*dim.z]) == NULL )
   {
      xvmc_error("write_dose_file",
                 "cannot allocate memory for dose error output data",8);
   }

   // write dose and error data to output memory
   float dose,error;
   float *out_count = out_data;
   float *err_count = err_data;
   for (register int k=0; k<dim.z; ++k)
   {
      for (register int j=0; j<dim.y; ++j)
      {
         for (register int i=0; i<dim.x; ++i)
         {
            dose  = out_dose->matrix[i][j][k];
            error = out_error->matrix[i][j][k];
            swap_bytes(dose);      swap_bytes(error);
            *out_count = dose;     *err_count = error;
            ++out_count;           ++err_count;
         }
      }
   }

   // write header file
   xvmc_message("Writing file:",hed_name,1);
   ofstream hed_file(hed_name,ios::out);
   if (!hed_file)
   {
      xvmc_error("write_dose_file","cannot open file",8);
   }
   hed_file.setf( ios::fixed, ios::floatfield );
   hed_file << setprecision(6);
   hed_file << "VOXELSIZE      |" << setw(11) << voxel_size.x
                                  << setw(11) << voxel_size.y
                                  << setw(11) << voxel_size.z << endl;
   hed_file << "DIMENSION      |" << setw(11) << dim.x
                                  << setw(11) << dim.y
                                  << setw(11) << dim.z << endl;
   hed_file << setprecision(6);
   hed_file << "DOSE-MAXIMUM   |" << setw(15) << dose_max << endl;
   hed_file << "CPU-TIME       |" << setw(15) << cpu_time << endl;
   hed_file << "POS-MAXIMUM    |" << setw(11) << d_max.x+1
                                  << setw(11) << d_max.y+1
                                  << setw(11) << d_max.z+1 << endl;
   hed_file << "END-INPUT      |" << endl;

   hed_file.close();

   // write dose data file
   xvmc_message("Writing file:",dat_name,0);
#ifdef OSF1
   FILE *dat_file=fopen( dat_name,"wb" );
   fwrite( out_data, sizeof(float), dim.x*dim.y*dim.z, dat_file );
   fclose( dat_file );
#else
   ofstream dat_file(dat_name,ios::out|ios::binary);
   char *out_dat_memory = (char *) out_data;
   dat_file.write( out_dat_memory, dim.x*dim.y*dim.z*sizeof(float));
   dat_file.close();
#endif

   // write dose error file
   xvmc_message("Writing file:",err_name,0);
#ifdef OSF1
   FILE *err_file=fopen( err_name,"wb" );
   fwrite( err_data, sizeof(float), dim.x*dim.y*dim.z, err_file );
   fclose( err_file );
#else
   ofstream err_file(err_name,ios::out|ios::binary);
   char *out_err_memory = (char *) err_data;
   err_file.write( out_err_memory, dim.x*dim.y*dim.z*sizeof(float));
   err_file.close();
#endif

   // free memory
   delete [] out_data; out_data = NULL;
   delete [] err_data; err_data = NULL;
}
// --- Added by JOKim 16Nov2010 -------------------------------------
// write dose matrix file
void write_portal_file(const array_3d<float>  *out_dose,
                       const array_3d<float>  *out_error,
                       const float dose_max, const int_3 d_max, float distance,
                       const char *hed_name, const char *dat_name, 
                       const char *err_name, const char *bmp_name)
{
   float *out_data = NULL;     // pointer to dose output data
   float *err_data = NULL;     // pointer to dose error output data

   // allocate memory for dose output data
   if ( (out_data = new float[dim_portal.x*dim_portal.y*dim_portal.z]) == NULL )
   {
      xvmc_error("write_portal_file",
                 "cannot allocate memory for dose output data",8);
   }

   // allocate memory for dose error output data
   if ( (err_data = new float[dim_portal.x*dim_portal.y*dim_portal.z]) == NULL )
   {
      xvmc_error("write_portal_file",
                 "cannot allocate memory for dose error output data",8);
   }

   // write dose and error data to output memory
   float dose,error;
   float *out_count = out_data;
   float *err_count = err_data;
   for (register int k=0; k<dim_portal.z; ++k)
   {
      for (register int j=0; j<dim_portal.y; ++j)
      {
         for (register int i=0; i<dim_portal.x; ++i)
         {
            dose  = out_dose->matrix[i][j][k];
            error = out_error->matrix[i][j][k];
            swap_bytes(dose);      swap_bytes(error);
            *out_count = dose;     *err_count = error;
            ++out_count;           ++err_count;
         }
      }
   }

   // write header file
   xvmc_message("Writing Portal File:",hed_name,1);
   ofstream hed_file(hed_name,ios::out);
   if (!hed_file)
   {
      xvmc_error("write_portal_file","cannot open file",8);
   }
   hed_file.setf( ios::fixed, ios::floatfield );
   hed_file << setprecision(6);
   hed_file << "VOXELSIZE      |" << setw(11) << voxel_size_portal.x
                                  << setw(11) << voxel_size_portal.y
                                  << setw(11) << voxel_size_portal.z << endl;
   hed_file << "DIMENSION      |" << setw(11) << dim_portal.x
                                  << setw(11) << dim_portal.y
                                  << setw(11) << dim_portal.z << endl;
   hed_file << "DISTANCE       |" << setw(11) << distance << endl;
   hed_file << setprecision(6);
   hed_file << "DOSE-MAXIMUM   |" << setw(15) << dose_max << endl;
   // hed_file << "CPU-TIME       |" << setw(15) << cpu_time << endl;
   hed_file << "POS-MAXIMUM    |" << setw(11) << d_max.x+1
                                  << setw(11) << d_max.y+1
                                  << setw(11) << d_max.z+1 << endl;
   hed_file << "END-INPUT      |" << endl;

   hed_file.close();

   // write dose data file
   xvmc_message("Writing file:",dat_name,0);
#ifdef OSF1
   FILE *dat_file=fopen( dat_name,"wb" );
   fwrite( out_data, sizeof(float), dim.x*dim.y*dim.z, dat_file );
   fclose( dat_file );
#else
   ofstream dat_file(dat_name,ios::out|ios::binary);
   char *out_dat_memory = (char *) out_data;
   dat_file.write( out_dat_memory, dim.x*dim.y*dim.z*sizeof(float));
   dat_file.close();
#endif

   // write dose error file
   xvmc_message("Writing Portal File:",err_name,0);
#ifdef OSF1
   FILE *err_file=fopen( err_name,"wb" );
   fwrite( err_data, sizeof(float), dim.x*dim.y*dim.z, err_file );
   fclose( err_file );
#else
   ofstream err_file(err_name,ios::out|ios::binary);
   char *out_err_memory = (char *) err_data;
   err_file.write( out_err_memory, dim.x*dim.y*dim.z*sizeof(float));
   err_file.close();
#endif

   // free memory
   delete [] out_data; out_data = NULL;
   delete [] err_data; err_data = NULL;
}
// write dose plane
void write_plane(const array_3d<float>  *out_dose,
                 const array_3d<float>  *out_error,
                 const plane_parameters *plane,
                 const char *pln_name)
{
   // print message and open file
   xvmc_message("Writing file:",pln_name,1);
   ofstream pln_file(pln_name,ios::out);
   if (!pln_file)
   {
      xvmc_error("write_plane","cannot open file",8);
   }
   pln_file.setf( ios::scientific );
   pln_file << setprecision(4);

   // misc variables
   int  lower,upper;
   real p,q;
   real dose,error;

   switch (plane->type)
   {
   case XY_PLANE:
      // print header
      pln_file << "#  Z:     " << setw(15) << plane->pos << endl;
      pln_file << "#   " << "   X           " << "   Y           "
                         << "   DOSE        " << "   ERROR" << endl;

      // calculate voxel indices for z interpolation
      lower = int((plane->pos-voxel_size.z/TWO)/voxel_size.z);
      upper = lower+1;
      if (lower <      0) { lower =       0; upper =       0; }
      if (upper >= dim.z) { lower = dim.z-1; upper = dim.z-1; }

      // interpolation factors
      p = real(upper)+ONE_HALF-plane->pos/voxel_size.z;
      q = ONE-p;

      // y loop
      for (register int j=0; j<dim.y; ++j)
      {
         // calculate y coordinate
         real pos_y = voxel_size.y*(real(j)+0.5);

         // x loop
         for (register int i=0; i<dim.x; ++i)
         {
            // calculate x coordinate
            real pos_x = voxel_size.x*(real(i)+0.5);

            // calculate dose
            dose = out_dose->matrix[i][j][lower]*p
                 + out_dose->matrix[i][j][upper]*q;

            // calculate error
            error = out_error->matrix[i][j][lower]*p
                  + out_error->matrix[i][j][upper]*q;

            // print results
            pln_file << setw(15) << pos_x
                     << setw(15) << pos_y
                     << setw(15) << dose
                     << setw(15) << error << endl;
         }
         pln_file <<  endl;
      }
      break;
   case XZ_PLANE:
      // print header
      pln_file << "#  Y:     " << setw(15) << plane->pos << endl;
      pln_file << "#   " << "   X           " << "   Z           "
                         << "   DOSE        " << "   ERROR" << endl;

      // calculate voxel indices for y interpolation
      lower = int((plane->pos-voxel_size.y/TWO)/voxel_size.y);
      upper = lower+1;
      if (lower <      0) { lower =       0; upper =       0; }
      if (upper >= dim.y) { lower = dim.y-1; upper = dim.y-1; }

      // interpolation factors
      p = real(upper)+ONE_HALF-plane->pos/voxel_size.y;
      q = ONE-p;

      // z loop
      for (register int k=0; k<dim.z; ++k)
      {
         // calculate z coordinate
         real pos_z = voxel_size.z*(real(k)+0.5);

         // x loop
         for (register int i=0; i<dim.x; ++i)
         {
            // calculate x coordinate
            real pos_x = voxel_size.x*(real(i)+0.5);

            // calculate dose
            dose = out_dose->matrix[i][lower][k]*p
                 + out_dose->matrix[i][upper][k]*q;

            // calculate error
            error = out_error->matrix[i][lower][k]*p
                  + out_error->matrix[i][upper][k]*q;

            // print results
            pln_file << setw(15) << pos_x
                     << setw(15) << pos_z
                     << setw(15) << dose
                     << setw(15) << error << endl;
         }
         pln_file <<  endl;
      }
      break;
   case YZ_PLANE:
      // print header
      pln_file << "#  X:     " << setw(15) << plane->pos << endl;
      pln_file << "#   " << "   Y           " << "   Z           "
                         << "   DOSE        " << "   ERROR" << endl;

      // calculate voxel indices for x interpolation
      lower = int((plane->pos-voxel_size.x/TWO)/voxel_size.x);
      upper = lower+1;
      if (lower <      0) { lower =       0; upper =       0; }
      if (upper >= dim.x) { lower = dim.x-1; upper = dim.x-1; }

      // interpolation factors
      p = real(upper)+ONE_HALF-plane->pos/voxel_size.x;
      q = ONE-p;

      // z loop
      for (register int k=0; k<dim.z; ++k)
      {
         // calculate z coordinate
         real pos_z = voxel_size.z*(real(k)+0.5);

         // y loop
         for (register int j=0; j<dim.y; ++j)
         {
            // calculate y coordinate
            real pos_y = voxel_size.y*(real(j)+0.5);

            // calculate dose
            dose = out_dose->matrix[lower][j][k]*p
                 + out_dose->matrix[upper][j][k]*q;

            // calculate error
            error = out_error->matrix[lower][j][k]*p
                  + out_error->matrix[upper][j][k]*q;

            // print results
            pln_file << setw(15) << pos_y
                     << setw(15) << pos_z
                     << setw(15) << dose
                     << setw(15) << error << endl;
         }
         pln_file <<  endl;
      }
      break;
   default:
      xvmc_error("write_plane","unknown plane type",8);
      break;
   }

   // close file and return
   pln_file.close();
   return;
}

// write dose profile
void write_profile(const array_3d<float>  *out_dose,
                   const array_3d<float>  *out_error,
                   const profile_parameters *profile,
                   const char *prf_name)
{
   // print message and open file
   xvmc_message("Writing file:",prf_name,1);
   ofstream prf_file(prf_name,ios::out);
   if (!prf_file)
   {
      xvmc_error("write_profile","cannot open file",8);
   }
   prf_file.setf( ios::scientific );
   prf_file << setprecision(4);

   // misc variables
   int  lower1,upper1,lower2,upper2;
   real p1,q1,p2,q2;
   real position,dose,error;

   switch (profile->type)
   {
   case X_PROFILE:
      // print header
      prf_file << "#  Y:     " << setw(15) << profile->pos.y << endl;
      prf_file << "#  Z:     " << setw(15) << profile->pos.z << endl;
      prf_file << "#  X" << endl;

      // calculate voxel indices for y interpolation
      lower1 = int((profile->pos.y-voxel_size.y/TWO)/voxel_size.y);
      upper1 = lower1+1;
      if (lower1 <      0) { lower1 =       0; upper1 =       0; }
      if (upper1 >= dim.y) { lower1 = dim.y-1; upper1 = dim.y-1; }

      // interpolation factors
      p1 = real(upper1)+ONE_HALF-profile->pos.y/voxel_size.y;
      q1 = ONE-p1;

      // calculate voxel indices for z interpolation
      lower2 = int((profile->pos.z-voxel_size.z/TWO)/voxel_size.z);
      upper2 = lower2+1;
      if (lower2 <      0) { lower2 =       0; upper2 =       0; }
      if (upper2 >= dim.z) { lower2 = dim.z-1; upper2 = dim.z-1; }

      // interpolation factors
      p2 = real(upper2)+ONE_HALF-profile->pos.z/voxel_size.z;
      q2 = ONE-p2;

      // x loop
      for (register int i=0; i<dim.x; ++i)
      {
         // calculate x-position
         position = voxel_size.x*(real(i)+0.5);
#ifdef REF_POINT
       	// Reference Point Shift
      	position -= ref_point.x;
#endif
         // calculate dose
         dose = out_dose->matrix[i][lower1][lower2]*p1*p2
              + out_dose->matrix[i][lower1][upper2]*p1*q2
              + out_dose->matrix[i][upper1][lower2]*q1*p2
              + out_dose->matrix[i][upper1][upper2]*q1*q2;

         // calculate error
         error = out_error->matrix[i][lower1][lower2]*p1*p2
               + out_error->matrix[i][lower1][upper2]*p1*q2
               + out_error->matrix[i][upper1][lower2]*q1*p2
               + out_error->matrix[i][upper1][upper2]*q1*q2;

         // print results
         prf_file << setw(15) << position
                  << setw(15) << dose
                  << setw(15) << error << endl;
      }
      break;
   case Y_PROFILE:
      // print header
      prf_file << "#  X:     " << setw(15) << profile->pos.x << endl;
      prf_file << "#  Z:     " << setw(15) << profile->pos.z << endl;
      prf_file << "#  Y" << endl;

      // calculate voxel indices for x interpolation
      lower1 = int((profile->pos.x-voxel_size.x/TWO)/voxel_size.x);
      upper1 = lower1+1;
      if (lower1 <      0) { lower1 =       0; upper1 =       0; }
      if (upper1 >= dim.x) { lower1 = dim.x-1; upper1 = dim.x-1; }

      // interpolation factors
      p1 = real(upper1)+ONE_HALF-profile->pos.x/voxel_size.x;
      q1 = ONE-p1;

      // calculate voxel indices for z interpolation
      lower2 = int((profile->pos.z-voxel_size.z/TWO)/voxel_size.z);
      upper2 = lower2+1;
      if (lower2 <      0) { lower2 =       0; upper2 =       0; }
      if (upper2 >= dim.z) { lower2 = dim.z-1; upper2 = dim.z-1; }

      // interpolation factors
      p2 = real(upper2)+ONE_HALF-profile->pos.z/voxel_size.z;
      q2 = ONE-p2;

      // y loop
      for (register int j=0; j<dim.y; ++j)
      {
         // calculate y position
         position = voxel_size.y*(real(j)+0.5);
#ifdef REF_POINT
       	// Reference Point Shift
      	position -= ref_point.y;
#endif
         // calculate dose
         dose = out_dose->matrix[lower1][j][lower2]*p1*p2
              + out_dose->matrix[lower1][j][upper2]*p1*q2
              + out_dose->matrix[upper1][j][lower2]*q1*p2
              + out_dose->matrix[upper1][j][upper2]*q1*q2;

         // calculate error
         error = out_error->matrix[lower1][j][lower2]*p1*p2
               + out_error->matrix[lower1][j][upper2]*p1*q2
               + out_error->matrix[upper1][j][lower2]*q1*p2
               + out_error->matrix[upper1][j][upper2]*q1*q2;

         // print results
         prf_file << setw(15) << position
                  << setw(15) << dose
                  << setw(15) << error << endl;
      }
      break;
   case Z_PROFILE:
      // print header
      prf_file << "#  X:     " << setw(15) << profile->pos.x << endl;
      prf_file << "#  Y:     " << setw(15) << profile->pos.y << endl;
      prf_file << "#  Z" << endl;

      // calculate voxel indices for x interpolation
      lower1 = int((profile->pos.x-voxel_size.x/TWO)/voxel_size.x);
      upper1 = lower1+1;
      if (lower1 <      0) { lower1 =       0; upper1 =       0; }
      if (upper1 >= dim.x) { lower1 = dim.x-1; upper1 = dim.x-1; }

      // interpolation factors
      p1 = real(upper1)+ONE_HALF-profile->pos.x/voxel_size.x;
      q1 = ONE-p1;

      // calculate voxel indices for y interpolation
      lower2 = int((profile->pos.y-voxel_size.y/TWO)/voxel_size.y);
      upper2 = lower2+1;
      if (lower2 <      0) { lower2 =       0; upper2 =       0; }
      if (upper2 >= dim.y) { lower2 = dim.y-1; upper2 = dim.y-1; }

      // interpolation factors
      p2 = real(upper2)+ONE_HALF-profile->pos.y/voxel_size.y;
      q2 = ONE-p2;

      // z loop
      for (register int k=0; k<dim.z; ++k)
      {
         // calculate z position
         position = voxel_size.z*(real(k)+0.5);

         // calculate dose
         dose = out_dose->matrix[lower1][lower2][k]*p1*p2
              + out_dose->matrix[lower1][upper2][k]*p1*q2
              + out_dose->matrix[upper1][lower2][k]*q1*p2
              + out_dose->matrix[upper1][upper2][k]*q1*q2;

         // calculate error
         error = out_error->matrix[lower1][lower2][k]*p1*p2
               + out_error->matrix[lower1][upper2][k]*p1*q2
               + out_error->matrix[upper1][lower2][k]*q1*p2
               + out_error->matrix[upper1][upper2][k]*q1*q2;

         // print results
         prf_file << setw(15) << position
                  << setw(15) << dose
                  << setw(15) << error << endl;
      }
      break;
   default:
      xvmc_error("write_profile","unknown profile type",8);
      break;
   }

   // close file and return
   prf_file.close();
   return;
}
