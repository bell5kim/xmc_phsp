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

#include "ranmar.h"  // Added by JOK 12/25/2009

// ****************************************************************
// class portal_dose: portal dose image
// ****************************************************************

class portal_dose
{
   public:
      // define portal dose image plane by the distance to the target
      // the image matrix dimension, the image matrix resolution and
      // the BMP image size
      portal_dose(double, unsigned int, unsigned int, unsigned int, 
                  double, double, double, int);

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
      void update_batch(double, double);

      // clear batch dose matrix
      void reset_batch_dose_portal(int);  
         
/*// Added by JOKim 25NOV2009 --------------------------------------------------
      { 
          portal_exit(); 	  
      }
// End of Added -------------------------------------------------------------- */

      // clear beam dose and error matrices
      void reset_beam_dose_portal(void) 
      { 
	 memset((void *) beam_dose_portal->data, 0, beam_dose_portal->size*sizeof(float));
	 memset((void *) beam_error_portal->data, 0, beam_error_portal->size*sizeof(float));
         // portal_exit(); 
      }

      // write portal dose files
      void write_beam_portal(const char *, const char *,
                       const char *, const char *, int, float, real);
		       
      // void write_bmp(const char *, int, int, int,int, BYTE *);

      // Added by JOKim 25NOV2009 --------------------------------------------------
      double        distance;  // Portal Dose Plane Distance from Isocenter
      int           bmp_size;  // BMP Size

      int           nZ_portal;  // number of layers in portal imager

      int           n_repeat_portal;  // number of repeatation in portal imager

      real_3        beam_dir;        // Beam Direction Vector
      real_3        beam_origin;     // Beam Origin in Portal Imager Coordinates
      real_3        portal_center;   // Portal Imager Center Position
      real_3        portal_distance; // Portal Imager Distance from Beam Origin
      real          iso_distance;    // Isocenter Distance 
      real          cos_alpha, sin_alpha, cos_beta, sin_beta;
      real_3        r_origin;        // Origin after rotated 
		
      real          portal_gantry_angle; // Gantry Angle Rotation
      real          portal_table_angle;  // Table Angle Rotation

      array_3d<double> *batch_dose_portal;


      // End of Added --------------------------------------------------------------

   private:
      void portal_exit(void) {
         xvmc_error("portal_dose",
            "In this version of XVMC portal images are not implemented",8); }
};


#endif /* _PORTAL_DOSE_H_ */
