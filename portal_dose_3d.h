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
#include "array_3d.h"

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
      void add(particle_parameters p); 
      // void add(particle_parameters p, ranmar rndm); // Added ranmar rndm by JOK 12-25-2009
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
      real	        dist2zero;  // Distance to Origin 
		
      real          portal_gantry_angle; // Gantry Angle Rotation
      real          portal_table_angle;  // Table Angle Rotation
     
      array_3d<float> *batch_portal_dose; 
      array_3d<float> *beam_portal_dose;
      array_3d<float> *beam_portal_error;
          
      // End of Added --------------------------------------------------------------

   private:
      void portal_exit(void) {
         xvmc_error("portal_dose",
            "In this version of XVMC portal images are not implemented",8); }
};


#endif /* _PORTAL_DOSE_H_ */
