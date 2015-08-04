#ifndef _PORTAL_DOSE_H_
#define _PORTAL_DOSE_H_

/*****************************************************************************
 * portal_dose.h:                                                            *
 *    class declaration:                                                     *
 *       portal_dose:    portal dose image (dummy)                           *
 *                                                                           *
 * Copyright (C) 2003    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 03/03/13        *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include "definitions.h"
#include "global.h"
#include "array_2d.h"
#include "ranmar.h"  // Added by JOK 12/25/2009

// ****************************************
// class multi_electron
// ****************************************

class MULTI_ELECTRON {
 public:

  friend void multi_electron_create(MULTI_ELECTRON &, real, ranmar &);
  friend void multi_electron_trace(MULTI_ELECTRON &,
                                   particle_parameters &, ranmar &); 

//  friend void multi_positron_create(MULTI_ELECTRON &, real, ranmar &);
//  friend void multi_positron_trace(MULTI_ELECTRON &,
//                                   particle_parameters &, ranmar &);

      // electron step data
      enum    {MAX_STEP=128};           // array size (maximum step number)
      int     n_step;                   // number of electron steps
      real    energy_loss[MAX_STEP];    // energy loss per step
      real    radiation_loss[MAX_STEP]; // bremsstrahlung loss per step
      real    sin_brems[MAX_STEP];      // sin of bremsstrahlung photon angle
      real    cos_brems[MAX_STEP];      // cos of bremsstrahlung photon angle
      real    dose_cor[MAX_STEP];       // correction for dose, ionization etc.
      real    random_number[MAX_STEP];  // divide the step into 2 sub-steps
      real    step_size[MAX_STEP];      // the length of this step
      real    s_cor[MAX_STEP];          // low energy correction for s_col
      real    alpha_r[MAX_STEP];        // alpha_r for this step
      real    reduced_angle[MAX_STEP];  // "reduced" MS (mult. scat.) angle
      real    sin_theta[MAX_STEP];      // direction change of the primary el.
      real    cos_theta[MAX_STEP];      // according to a discrete interaction
      int     i_delta[MAX_STEP];        // index of the delta electron
      bool    moller[MAX_STEP];         // true if there is a Moller interaction
                                        // (delta electron) during the step
      bool    bhabha[MAX_STEP];         // true if there is a Bhabha interaction
                                        // (delta electron) during the step

      // stack for delta (secondary) electron data
      enum    {MAX_DELTA=100};         // array size (maximum number)
      int     n_delta;                 // number of delta electrons in stack
      real    energy_delta[MAX_DELTA]; // delta electron energy
      real    sin_theta_d[MAX_DELTA];  // sin of the longitudinal scat. angle
      real    cos_theta_d[MAX_DELTA];  // cos of the longitudinal scat. angle


};


// ****************************************************************
// class portal_dose: portal dose image
// ****************************************************************

class portal_dose
{
   public:
      // define portal dose image plane by the distance to the target
      // the image matrix dimension, the image matrix resolution and
      // the BMP image size
      portal_dose(double, unsigned int, unsigned int, double, double, int);

      // destructor
      ~portal_dose(void);

      // move the portal dose plane to the correct place using the
      // coordinates of the target and the beam direction vector
      void move(real_3, real, real);

      // add portal dose contribution of one photon
      void add(particle_parameters p, ranmar rndm); // Added ranmar rndm by JOK 12-25-2009
      //{ r
         // portal_exit(); 
	 // return(false); 
      //}

      // update dose and dose error for one batch
      void update_batch(double);

      // clear batch dose matrix
      void reset_batch_dose_portal(void);  
         
/*// Added by JOKim 25NOV2009 --------------------------------------------------
      { 
          portal_exit(); 	  
      }
// End of Added -------------------------------------------------------------- */

      // clear beam dose and error matrices
      void reset_beam_dose_portal(void) 
      { 
         // portal_exit(); 
      }

      // write portal dose files
      void write_files(const char *, const char *,
                       const char *, const char *, int);
		       
      // void write_bmp(const char *, int, int, int,int, BYTE *);

      // Added by JOKim 25NOV2009 --------------------------------------------------
      double        distance;  // Portal Dose Plane Distance from Isocenter
      int           bmp_size;  // BMP Size
		
      int_3         p_dim;          // Portal Imager Dimension
      real_3        p_size;         // Portal Imager Resolution, Pixel Size
      real_3        p_size_rot;     // Pixel Size after Rotation
      real_3        beam_dir;       // Beam Direction Vector
      real_3        p_pos;          // Portal Imager Center Position
      real_3        p_dir;          // Portal Imager Center Direction Vector
      real_3        x1, x2, y1, y2; // Portal Imager Corners
      real_3        min_limit;      // Portal Imager Minimum Limit in Phantom Coordinates
      real_3        max_limit;      // Portal Imager Maximum Limit in Phantom Coordinates
      real          iso_distance;   // Isocenter Distance 
      real          cos_alpha, sin_alpha, cos_beta, sin_beta;
      real_3        normal;         // Normal Vector of Portal Imager
      real	    dist2zero;  // Distance to Origin 
		
      real          portal_gantry_angle; // Gantry Angle Rotation
      real          portal_table_angle;  // Table Angle Rotation
     
      array_2d<float> *batch_dose_portal; 
      array_2d<float> *beam_dose;
      array_2d<float> *beam_portal_error;
          
      // End of Added --------------------------------------------------------------

   private:
      void portal_exit(void) {
         xvmc_error("portal_dose",
            "In this version of XVMC portal images are not implemented",8); }
      void portal_kerma(particle_parameters &p, ranmar &);
      void portal_photon(particle_parameters &p, ranmar &);
      void multi_electron_create(MULTI_ELECTRON &, real, ranmar &);
      void multi_electron_trace(MULTI_ELECTRON &, particle_parameters &, ranmar &);
      void one_electron(particle_parameters &, ranmar &);
};


#endif /* _PORTAL_DOSE_H_ */
