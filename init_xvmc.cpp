/*****************************************************************************
 * init_xvmc.cpp:                                                            *
 *    function init_xvmc: read and initialize data                           *
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
#include <math.h>

#include <iostream>
#include <stdio.h>
using namespace std;

#include <unistd.h>
#include <sys/utsname.h>
#include <time.h>
#include <assert.h>

#include "definitions.h"
#include "global.h"
#include "electron_data_inp.h"
#include "photon_data_inp.h"
#include "moller_XS.h"
#include "bhabha_XS.h"
#include "brems_XS.h"
#include "compton_XS_diff.h"
#include "compton_XS_total.h"
#include "pair_XS_diff.h"
#include "pair_XS_total.h"

#include "treatment_plan.h"

// ****************************************
// define global variables
// ****************************************

// number of parallel processes (threads) to calculate dose,
// equal to the number of active processors in the system
int n_process = 1;

electron_data_inp *e_inp = NULL; // input array for electron transport data
photon_data_inp   *p_inp = NULL; // input array for photon transport data

// Moller cross section for water
moller_XS          h2o_moller(e_cut);

// Bhabha cross section for water
bhabha_XS          h2o_bhabha(e_cut);

// bremsstrahlung cross section forwater
brems_XS           h2o_brems(TC);

compton_XS_diff    compton;    // differential Compton cross section data
compton_XS_total   tot_comp;   // total Compton cross section for water
pair_XS_diff       dif_pair;   // differential pair cross section data
pair_XS_total      tot_pair;   // total pair cross section data for water

char  *inout_path;   // path to the input and output files
char  *base_path;    // path to the base data files

treatment_plan plan;          // treatment plan parameters

real_3            ref_point  = {0.0, 0.0, 0.0}; // reference point

real_3            voxel_size = {0.0, 0.0, 0.0}; // cube voxel sizes (cm)
int_3             dim        = {0, 0, 0};       // cube dimensions
real_3            cube_size  = {0.0, 0.0, 0.0}; // cube sizes (cm)

// density matrix
array_3d<float>  *density   = NULL;
// density correction for the collison stopping power (eq. (19) in VMC paper I)
array_3d<float>  *dens_scol = NULL;
// low energy correction for "dens_scol"
array_3d<float>  *dens_ccol = NULL;
// density correction for radiation stopping power
// (eq. (20) and (24) in VMC paper I)
array_3d<float>  *dens_srad = NULL;
// density correction for multiple scattering (see VMC paper II)
array_3d<float>  *dens_fchi = NULL;
// density correction for Compton cross section = electron density
// (eq. (10) in XVMC paper I)
array_3d<float>  *dens_comp = NULL;
// density correction for pair cross section,
// similar to "dens_srad" (Feynman crossing)
array_3d<float>  *dens_pair = NULL;
// density correction for photo cross section
array_3d<float>  *dens_phot = NULL;

// dose matrix for present beam
array_3d<float>  *beam_dose  = NULL;
// dose error matrix for present beam
array_3d<float>  *beam_error = NULL;
// total dose matrix
array_3d<float>  *sum_dose   = NULL;
// total dose error matrix
array_3d<float>  *sum_error  = NULL;

// --- Added by JKim 15Nov2010 ----------------------------
real_3    voxel_size_portal = {0.0, 0.0, 0.0}; // cube voxel sizes (cm)
int_3     dim_portal        = {0, 0, 0};       // cube dimensions
real_3    cube_size_portal  = {0.0, 0.0, 0.0}; // cube sizes (cm)

// density matrix
array_3d<float>  *density_portal   = NULL;
// density correction for the collison stopping power (eq. (19) in VMC paper I)
array_3d<float>  *dens_scol_portal = NULL;
// low energy correction for "dens_scol"
array_3d<float>  *dens_ccol_portal = NULL;
// density correction for radiation stopping power
// (eq. (20) and (24) in VMC paper I)
array_3d<float>  *dens_srad_portal = NULL;
// density correction for multiple scattering (see VMC paper II)
array_3d<float>  *dens_fchi_portal = NULL;
// density correction for Compton cross section = electron density
// (eq. (10) in XVMC paper I)
array_3d<float>  *dens_comp_portal = NULL;
// density correction for pair cross section,
// similar to "dens_srad" (Feynman crossing)
array_3d<float>  *dens_pair_portal = NULL;
// density correction for photo cross section
array_3d<float>  *dens_phot_portal = NULL;

// dose matrix for present beam
array_3d<float>  *beam_dose_portal  = NULL;
// dose error matrix for present beam
array_3d<float>  *beam_error_portal = NULL;
// total dose matrix
array_3d<float>  *sum_dose_portal   = NULL;
// total dose error matrix
array_3d<float>  *sum_error_portal  = NULL;
// --- End of Adding 

#ifdef CHECK_ENERGY
sum_energy_type tot_energy; // test energy conservation
#endif // CHECK_ENERGY

// ****************************************
// declare functions
// ****************************************

// file input and output
void read_inp_file(const char *, const char *, const char *);
void read_density_file(const char *, const char *);

// Table inits: shared with hyperion version
bool InitCrossSectionTables();

void init_xvmc(const char *patient_name, const char *plan_name)
{
   char *env_value;              // environment variable
   int   slength;                // length of a string
   char *water_file;             // water cross section file name
   char *compton_file;           // diff. Compton cross section data file
   char *pair_file;              // diff. pair cross section data file
   char *tot_pair_file;          // total pair cross section data file
   char *xvmc_inp_file;          // XVMC input file name
   char *dens_hed_file;          // density matrix header file (ASCII)
   char *dens_dat_file;          // density matrix data file (binary)

   /**************************************************
    * get system information
    **************************************************/

   time_t      *nulltime = NULL;
   time_t       timer;
   struct  tm  *br_time;
   timer = time (nulltime);
   br_time = localtime (&timer);
   cout << "XVMC>" << endl;
   cout << "XVMC> Date:                 " << asctime(br_time);

   struct utsname uts;
   uname(&uts);
   xvmc_message("Operating system:    ", uts.sysname,1);
   xvmc_message("Network node:        ", uts.nodename,0);
   xvmc_message("System release:      ", uts.release,0);
   xvmc_message("System version:      ", uts.version,0);
   xvmc_message("Machine type:        ", uts.machine,0);
#ifdef NEW_NU
   printf ("**********  Modification:   Nu Value = %f ***********\n", NEW_NU);
#endif
#ifdef SPECTRUM_WEIGHT
   printf ("**********  Modification:   Spectrum Weighting is applied  ***********\n");
#endif
#ifdef PHSP_WRITE
   printf ("**********  Modification:   Phase Space Writing is applied  ***********\n");
#endif
#ifdef PHSP_READ
   printf ("PHSP_READ> **********  Phase Space Reading Mode  ***********\n");
#endif
#ifdef VARIAN
   printf ("**********  Modification:   Varian Attenuation acording to Yu  ***********\n");
#endif

   /**************************************************
    * get the number of active processors
    **************************************************/

#ifdef ONE_PROCESS
   n_process = 1;
#else
   n_process = sysconf(_SC_NPROCESSORS_ONLN);
#endif
   if (n_process < 1) n_process = 1;
   xvmc_message("Number of processors:", n_process,"",1);

   /**************************************************
    * read water cross section data from files
    **************************************************/

   // get xvmc home directory name
   env_value = getenv("XVMC_HOME");
   if (env_value == NULL)
   {
      xvmc_error("init_xvmc","environment variable XVMC_HOME not set",8);
   }

   // calculate string length for file names
   slength       = strlen(env_value) + 64;

   // allocate memory for file names
   water_file    = new char[slength];
   compton_file  = new char[slength];
   pair_file     = new char[slength];
   tot_pair_file = new char[slength];

   // copy xvmc home directory name into file name strings
   strcpy(water_file,    env_value);
   strcpy(compton_file,  env_value);
   strcpy(pair_file,     env_value);
   strcpy(tot_pair_file, env_value);

   // add file names
   strcat(water_file,    "/dat/water.data");
   strcat(compton_file,  "/dat/compton.data");
   strcat(pair_file,     "/dat/pair.data");
   strcat(tot_pair_file, "/dat/tot_pair.data");

   // input electron transport parameters from file
   if ((e_inp = new electron_data_inp(water_file)) == NULL)
   {
      xvmc_error("init_xvmc","cannot read electron transport data from file",8);
   }

   // input photon transport parameters from file
   if ((p_inp = new photon_data_inp(water_file)) == NULL)
   {
      xvmc_error("init_xvmc","cannot read photon transport data from file",8);
   }

   // input differential Compton cross section data
   compton.read(compton_file);

   // input differential pair cross section data
   dif_pair.read(pair_file);

   // input total pair cross section data
   tot_pair.read(tot_pair_file);

   // free memory for file names
   delete [] water_file;
   delete [] compton_file;
   delete [] pair_file;
   delete [] tot_pair_file;

   /**************************************************
    * generate names for input-output files
    **************************************************/

   // calculate string length for file names
   if (strlen(patient_name) > strlen(plan_name)) {
      slength       = 2 * strlen(patient_name) + 20; }
   else {
      slength       = 2 * strlen(plan_name) + 20; }

   // get xvmc working directory name
   env_value = getenv("XVMC_WORK");
   if (env_value != NULL)
   {
      // calculate new string length for file names
      slength       = slength + strlen(env_value);

      // allocate memory for file names
      inout_path    = new char[slength];
      dens_hed_file = new char[slength];
      dens_dat_file = new char[slength];
      base_path     = new char[slength];

      // copy xvmc working directory name into file name strings
      strcpy(inout_path,    env_value);
      strcpy(dens_hed_file, env_value);
      strcpy(dens_dat_file, env_value);
      strcpy(base_path,     env_value);
   }
   else
   {  // if XVMC_WORK has not been set use present directory

      // allocate memory for file names
      inout_path    = new char[slength];
      dens_hed_file = new char[slength];
      dens_dat_file = new char[slength];
      base_path     = new char[slength];

      // set xvmc working directory
      strcpy(inout_path,    ".");
      strcpy(dens_hed_file, ".");
      strcpy(dens_dat_file, ".");
      strcpy(base_path,     ".");
   }

   // add file names (without extension)
   strcat(inout_path,    "/");
   strcat(inout_path,    patient_name);
   strcat(inout_path,    "/");
   strcat(inout_path,    plan_name);

   strcat(dens_hed_file, "/");
   strcat(dens_hed_file, patient_name);
   strcat(dens_hed_file, "/");
   strcat(dens_hed_file, patient_name);
   strcat(dens_hed_file, ".hed");  // density matrix header file (ASCII)

   strcat(dens_dat_file, "/");
   strcat(dens_dat_file, patient_name);
   strcat(dens_dat_file, "/");
   strcat(dens_dat_file, patient_name);
   strcat(dens_dat_file, ".dmx");  // density matrix data file (binary)

   strcat(base_path,     "/dat/");

   /**************************************************
    * read input file: *.vmc and density data
    **************************************************/

   // for patient initialize treatment plan without beams
   plan.init(patient_name, plan_name);

   // generate name for input file
   slength = strlen(inout_path)+10;
   xvmc_inp_file  = new char[slength];
   strcpy(xvmc_inp_file, inout_path);
   strcat(xvmc_inp_file, ".vmc");

   read_inp_file(xvmc_inp_file, dens_hed_file, dens_dat_file);
      // the density matrix file names are needed to save
      // the phantom density distribution (if it is created)

   // if the density matrix does not exist by now, we try to read it from file
   if (density == NULL)
   {
      read_density_file(dens_hed_file, dens_dat_file);
   }

   // delete file names
   delete [] xvmc_inp_file;
   delete [] dens_hed_file;
   delete [] dens_dat_file;

// --- Added by JOKim 20Nov2010  ------------------------
   // Portal Contruction -----------------------
//   for (register int k=10; k<12; ++k) {
//      for (register int j=0; j<dim_portal.y; ++j) {
//      	 for (register int i=0; i<dim_portal.x; ++i) {
//             density_portal->matrix[i][j][k]  = 2.0;
//      	 }
//      }
//   }
// --- End of Adding --- 20Nov2010 ----------------------
   // print transport parameters
   xvmc_message("electron step size:",e_step,"",1);
   xvmc_message("electron cut-off:  ",e_cut,"",0);
   xvmc_message("photon cut-off:    ",p_cut,"",0);
   xvmc_message("pri. KERMA cut-off:",k0_cut,"",0);
   xvmc_message("sec. KERMA cut-off:",k1_cut,"",0);

   // print reference point
   xvmc_message("ref_point.x:",ref_point.x,"cm",1);
   xvmc_message("ref_point.y:",ref_point.y,"cm",0);
   xvmc_message("ref_point.z:",ref_point.z,"cm",0);

   /**************************************************
    * calculate cross section density corrections
    **************************************************/

   // initialize electron cross section density corrections
   if ( (dens_scol = new array_3d<float>(dim, ONE)) == NULL ) {
      xvmc_error("init_xvmc","cannot construct dens_scol",8); }

   if ( (dens_ccol = new array_3d<float>(dim, ONE)) == NULL ) {
      xvmc_error("init_xvmc","cannot construct dens_ccol",8); }

   if ( (dens_srad = new array_3d<float>(dim, ONE)) == NULL ) {
      xvmc_error("init_xvmc","cannot construct dens_srad",8); }

   if ( (dens_fchi = new array_3d<float>(dim, ONE)) == NULL ) {
      xvmc_error("init_xvmc","cannot construct dens_fchi",8); }

   // initialize photon cross section arrays with negative value, therefore
   // the fit functions to calculate the corrections from density are used
   // Compton effect
   if (dens_comp == NULL)
   {
      if ( (dens_comp = new array_3d<float>(dim, MINUS_ONE)) == NULL ) {
         xvmc_error("init_xvmc","cannot construct dens_comp",8); }
   }

   // pair production
   if (dens_pair == NULL)
   {
      if ( (dens_pair = new array_3d<float>(dim, MINUS_ONE)) == NULL ) {
         xvmc_error("init_xvmc","cannot construct dens_pair",8); }
   }

   // photo effect
   if (dens_phot == NULL)
   {
      if ( (dens_phot = new array_3d<float>(dim, MINUS_ONE)) == NULL ) {
         xvmc_error("init_xvmc","cannot construct dens_phot",8); }
   }

// --- Added by JOKim 15Nov2010 -----------------------------------
   // initialize electron cross section density corrections
   if ( (dens_scol_portal = new array_3d<float>(dim, ONE)) == NULL ) {
      xvmc_error("init_xvmc","cannot construct dens_scol_portal",8); }

   if ( (dens_ccol_portal = new array_3d<float>(dim, ONE)) == NULL ) {
      xvmc_error("init_xvmc","cannot construct dens_ccol_portal",8); }

   if ( (dens_srad_portal = new array_3d<float>(dim, ONE)) == NULL ) {
      xvmc_error("init_xvmc","cannot construct dens_srad_portal",8); }

   if ( (dens_fchi_portal = new array_3d<float>(dim, ONE)) == NULL ) {
      xvmc_error("init_xvmc","cannot construct dens_fchi_portal",8); }

   // initialize photon cross section arrays with negative value, therefore
   // the fit functions to calculate the corrections from density are used
   // Compton effect
   if (dens_comp_portal == NULL)
   {
      if ( (dens_comp_portal = new array_3d<float>(dim, MINUS_ONE)) == NULL ) {
         xvmc_error("init_xvmc","cannot construct dens_comp_portal",8); }
   }

   // pair production
   if (dens_pair_portal == NULL)
   {
      if ( (dens_pair_portal = new array_3d<float>(dim, MINUS_ONE)) == NULL ) {
         xvmc_error("init_xvmc","cannot construct dens_pair_portal",8); }
   }

   // photo effect
   if (dens_phot_portal == NULL)
   {
      if ( (dens_phot_portal = new array_3d<float>(dim, MINUS_ONE)) == NULL ) {
         xvmc_error("init_xvmc","cannot construct dens_phot_portal",8); }
   }
// --- End of Adding  ---------------------------------------------

   InitCrossSectionTables();

   /**************************************************
    * initialize arrays for dose and dose error
    **************************************************/

   xvmc_message("Initialize dose arrays",1);

   if ( (beam_dose = new array_3d<float>(dim, ZERO)) == NULL ) {
     xvmc_error("init_xvmc","cannot construct beam_dose",8); }

   if ( (beam_error = new array_3d<float>(dim, ZERO)) == NULL ) {
     xvmc_error("init_xvmc","cannot construct beam_error",8); }

   if ( (sum_dose = new array_3d<float>(dim, ZERO)) == NULL ) {
      xvmc_error("init_xvmc","cannot construct sum_dose",8); }

   if ( (sum_error = new array_3d<float>(dim, ZERO)) == NULL ) {
      xvmc_error("init_xvmc","cannot construct sum_error",8); }    

// --- Added by JOKim 15Nov2010 -----------------------------------
   if ( (beam_dose_portal = new array_3d<float>(dim, ZERO)) == NULL ) {
     xvmc_error("init_xvmc","cannot construct beam_dose_portal",8); }

   if ( (beam_error_portal = new array_3d<float>(dim, ZERO)) == NULL ) {
     xvmc_error("init_xvmc","cannot construct beam_error_portal",8); }

   if ( (sum_dose_portal = new array_3d<float>(dim, ZERO)) == NULL ) {
      xvmc_error("init_xvmc","cannot construct sum_dose_portal",8); }

   if ( (sum_error_portal = new array_3d<float>(dim, ZERO)) == NULL ) {
      xvmc_error("init_xvmc","cannot construct sum_error_portal",8); }    

// --- End of Adding  ---------------------------------------------

}
