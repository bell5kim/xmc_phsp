#include <math.h>
#include <string>
#include <assert.h>
#include "definitions.h"
#include "global.h"

//#####################################################################################################

bool InitCrossSectionTables()
{
  assert(density!=NULL);
  assert(dens_scol!=NULL);
  assert(dens_ccol!=NULL);
  assert(dens_srad!=NULL);
  assert(dens_fchi!=NULL);
  assert(dens_comp!=NULL);
  assert(dens_pair!=NULL);
  assert(dens_phot!=NULL);

  register int k, j, i;
  for (i=0; i<dim.x; ++i) {
    for (j=0; j<dim.y; ++j) {
      for (k=0; k<dim.z; ++k) {
	// get mass density
	real rho = density->matrix[i][j][k];

	// density correction for collison stopping power
	if (rho > 0.821) {
	  dens_scol->matrix[i][j][k] =
	    pow(double(rho),-0.227)+0.038*(rho-1.0);
	}
	else
	  {
	    dens_scol->matrix[i][j][k] = 1.039;
	  }
	if (fabs(rho-1.0) < 0.02) dens_scol->matrix[i][j][k] = 1.0;
	// this is ICRP bone from ICRU 37
	// if (fabs(rho-1.85) < 0.02) dens_scol->matrix[i][j][k] = 0.912;
	// this is LUNG521ICRU
	// if (fabs(rho-0.26) < 0.02) dens_scol->matrix[i][j][k] = 1.055;
	// this is AL521ICRU
	if (fabs(rho-2.70) < 0.02) dens_scol->matrix[i][j][k] = 0.832;

	// density correction for radiation stopping power
	// (eq. (20) in VMC paper I)
	if (rho > 0.9)
	  {
	    dens_srad->matrix[i][j][k] = 1.13 + 0.56*log(rho-0.3);
            }
	else
	  {
	    dens_srad->matrix[i][j][k] = 1.049-0.228*rho;
	  }
	if (fabs(rho-1.0) < 0.02) dens_srad->matrix[i][j][k] = 1.0;
	// this is ICRP bone from ICRU 37
	// if (fabs(rho-1.85) < 0.02) dens_srad->matrix[i][j][k] = 1.35;

	// we change to f-tilde, ratio of f_r to f_c
	// (eq. (24) in VMC paper I)
	dens_srad->matrix[i][j][k] /= dens_scol->matrix[i][j][k];

	// mass  --> linear collision stopping power
	dens_scol->matrix[i][j][k] *= rho;

	// low energy correction for "dens_scol"
	if (rho > 1.1)
	  {
	    dens_ccol->matrix[i][j][k] = 0.012*rho-0.013;
	  }
	else
	  {
	    dens_ccol->matrix[i][j][k] = 0.05553-0.0617*rho;
	  }
	// water
	if (fabs(rho-1.0) < 0.1) dens_ccol->matrix[i][j][k] = 0.0;
	// this is LUNG521ICRU
	// if (fabs(rho-0.26) < 0.02) dens_ccol->matrix[i][j][k] = 0.0;
	// this is AL521ICRU
	if (fabs(rho-2.70) < 0.02) dens_ccol->matrix[i][j][k] = 0.03;
	// mass  --> linear collision stopping power
	dens_ccol->matrix[i][j][k] *= rho;

	// density correction for multiple scattering
	if (rho > 1.2)
	  {
	    dens_fchi->matrix[i][j][k] = 0.916+1.047*log(rho);
	  }
	else
	  {
	    if (rho > 1.138)
	      {
		dens_fchi->matrix[i][j][k] = 1.107;
	      }
	    else
	      {
		if (rho > 0.92)
		  {
		    dens_fchi->matrix[i][j][k] = -0.352+1.282*rho;
                  }
		else
                  {
		    dens_fchi->matrix[i][j][k] = 0.993-0.18*rho;
                  }
	      }
	  }
	// water
	if (fabs(rho-1.0) < 0.02)  dens_fchi->matrix[i][j][k] = 1.0;
	// this is AL521ICRU
	if (fabs(rho-2.70) < 0.02) dens_fchi->matrix[i][j][k] = 2.055;

	// density correction for Compton cross section (electron density)
	if (dens_comp->matrix[i][j][k] < 0.0)
	  {
	    if (rho > 1.0)
	      {
		dens_comp->matrix[i][j][k] = 0.85*rho+0.15;
	      }
	    else
	      {
		//dens_comp->matrix[i][j][k] = rho;
		dens_comp->matrix[i][j][k] = 0.99*rho+0.01*rho*rho;
	      }
	  }
	// this is AL521ICRU
	if (fabs(rho-2.70) < 0.02) dens_comp->matrix[i][j][k] = 2.3197;

	// density correction for pair cross section
	if (dens_pair->matrix[i][j][k] < 0.0)
	  {
	    if (rho > 0.92)
	      {
		dens_pair->matrix[i][j][k] = 1.15 + 0.66*log(rho-0.3);
	      }
	    else
	      {
		dens_pair->matrix[i][j][k] = 1.049-0.228*rho;
	      }
	    if (fabs(rho-1.0) < 0.02) dens_pair->matrix[i][j][k] = 1.0;
	  }
	// this is AL521ICRU
	if (fabs(rho-2.70) < 0.02) dens_pair->matrix[i][j][k] = 1.638;
	// mass attenuation --> linear attenuation
	dens_pair->matrix[i][j][k] *= rho;

	// density correction for photo cross section
	if (dens_phot->matrix[i][j][k] < 0.0)
	  {
	    if (rho > 1.1)
	      {
		dens_phot->matrix[i][j][k] = 1.0 + 8.0*sqrt(rho-1.1);
	      }
	    else
	      {
		if (rho >= 1.0)
		  {
		    dens_phot->matrix[i][j][k] = 1.0;
		  }
		else
		  {
		    if (rho > 0.9)
		      {
			    dens_phot->matrix[i][j][k] = 1.0-6.0*(1.0-rho);
		      }
		    else
		      {
			dens_phot->matrix[i][j][k] = 1.07; // this is lung
		      }
		  }
	      }
	  }
	// mass attenuation --> linear attenuation
	dens_phot->matrix[i][j][k] *= rho;

      }  // for (register int k=0; k<dim.z; ++k)
    }  // for (register int j=0; j<dim.y; ++j)
  }  // for (register int i=0; i<dim.x; ++i)


// --- Added by JKim 15NOV2010 -----------------------------------
  if (density_portal != NULL) {
     assert(density_portal!=NULL);
     assert(dens_scol_portal!=NULL);
     assert(dens_ccol_portal!=NULL);
     assert(dens_srad_portal!=NULL);
     assert(dens_fchi_portal!=NULL);
     assert(dens_comp_portal!=NULL);
     assert(dens_pair_portal!=NULL);
     assert(dens_phot_portal!=NULL);

     for (i=0; i<dim_portal.x; ++i) {
       for (j=0; j<dim_portal.y; ++j) {
	 for (k=0; k<dim_portal.z; ++k) {
	   // get mass density
	   real rho = density_portal->matrix[i][j][k];

	   // density correction for collison stopping power
	   if (rho > 0.821) {
	     dens_scol_portal->matrix[i][j][k] =
	       pow(double(rho),-0.227)+0.038*(rho-1.0);
	   }
	   else
	     {
	       dens_scol_portal->matrix[i][j][k] = 1.039;
	     }
	   if (fabs(rho-1.0) < 0.02) dens_scol_portal->matrix[i][j][k] = 1.0;
	   // this is ICRP bone from ICRU 37
	   // if (fabs(rho-1.85) < 0.02) dens_scol_portal->matrix[i][j][k] = 0.912;
	   // this is LUNG521ICRU
	   // if (fabs(rho-0.26) < 0.02) dens_scol_portal->matrix[i][j][k] = 1.055;
	   // this is AL521ICRU
	   if (fabs(rho-2.70) < 0.02) dens_scol_portal->matrix[i][j][k] = 0.832;

	   // density correction for radiation stopping power
	   // (eq. (20) in VMC paper I)
	   if (rho > 0.9)
	     {
	       dens_srad_portal->matrix[i][j][k] = 1.13 + 0.56*log(rho-0.3);
               }
	   else
	     {
	       dens_srad_portal->matrix[i][j][k] = 1.049-0.228*rho;
	     }
	   if (fabs(rho-1.0) < 0.02) dens_srad_portal->matrix[i][j][k] = 1.0;
	   // this is ICRP bone from ICRU 37
	   // if (fabs(rho-1.85) < 0.02) dens_srad_portal->matrix[i][j][k] = 1.35;

	   // we change to f-tilde, ratio of f_r to f_c
	   // (eq. (24) in VMC paper I)
	   dens_srad_portal->matrix[i][j][k] /= dens_scol_portal->matrix[i][j][k];

	   // mass  --> linear collision stopping power
	   dens_scol_portal->matrix[i][j][k] *= rho;

	   // low energy correction for "dens_scol"
	   if (rho > 1.1)
	     {
	       dens_ccol_portal->matrix[i][j][k] = 0.012*rho-0.013;
	     }
	   else
	     {
	       dens_ccol_portal->matrix[i][j][k] = 0.05553-0.0617*rho;
	     }
	   // water
	   if (fabs(rho-1.0) < 0.1) dens_ccol_portal->matrix[i][j][k] = 0.0;
	   // this is LUNG521ICRU
	   // if (fabs(rho-0.26) < 0.02) dens_ccol_portal->matrix[i][j][k] = 0.0;
	   // this is AL521ICRU
	   if (fabs(rho-2.70) < 0.02) dens_ccol_portal->matrix[i][j][k] = 0.03;
	   // mass  --> linear collision stopping power
	   dens_ccol_portal->matrix[i][j][k] *= rho;

	   // density correction for multiple scattering
	   if (rho > 1.2)
	     {
	       dens_fchi_portal->matrix[i][j][k] = 0.916+1.047*log(rho);
	     }
	   else
	     {
	       if (rho > 1.138)
		 {
		   dens_fchi_portal->matrix[i][j][k] = 1.107;
		 }
	       else
		 {
		   if (rho > 0.92)
		     {
		       dens_fchi_portal->matrix[i][j][k] = -0.352+1.282*rho;
                     }
		   else
                     {
		       dens_fchi_portal->matrix[i][j][k] = 0.993-0.18*rho;
                     }
		 }
	     }
	   // water
	   if (fabs(rho-1.0) < 0.02)  dens_fchi_portal->matrix[i][j][k] = 1.0;
	   // this is AL521ICRU
	   if (fabs(rho-2.70) < 0.02) dens_fchi_portal->matrix[i][j][k] = 2.055;

	   // density correction for Compton cross section (electron density)
	   if (dens_comp_portal->matrix[i][j][k] < 0.0)
	     {
	       if (rho > 1.0)
		 {
		   dens_comp_portal->matrix[i][j][k] = 0.85*rho+0.15;
		 }
	       else
		 {
		   //dens_comp_portal->matrix[i][j][k] = rho;
		   dens_comp_portal->matrix[i][j][k] = 0.99*rho+0.01*rho*rho;
		 }
	     }
	   // this is AL521ICRU
	   if (fabs(rho-2.70) < 0.02) dens_comp_portal->matrix[i][j][k] = 2.3197;

	   // density correction for pair cross section
	   if (dens_pair_portal->matrix[i][j][k] < 0.0)
	     {
	       if (rho > 0.92)
		 {
		   dens_pair_portal->matrix[i][j][k] = 1.15 + 0.66*log(rho-0.3);
		 }
	       else
		 {
		   dens_pair_portal->matrix[i][j][k] = 1.049-0.228*rho;
		 }
	       if (fabs(rho-1.0) < 0.02) dens_pair_portal->matrix[i][j][k] = 1.0;
	     }
	   // this is AL521ICRU
	   if (fabs(rho-2.70) < 0.02) dens_pair_portal->matrix[i][j][k] = 1.638;
	   // mass attenuation --> linear attenuation
	   dens_pair_portal->matrix[i][j][k] *= rho;

	   // density correction for photo cross section
	   if (dens_phot_portal->matrix[i][j][k] < 0.0)
	     {
	       if (rho > 1.1)
		 {
		   dens_phot_portal->matrix[i][j][k] = 1.0 + 8.0*sqrt(rho-1.1);
		 }
	       else
		 {
		   if (rho >= 1.0)
		     {
		       dens_phot_portal->matrix[i][j][k] = 1.0;
		     }
		   else
		     {
		       if (rho > 0.9)
			 {
			       dens_phot_portal->matrix[i][j][k] = 1.0-6.0*(1.0-rho);
			 }
		       else
			 {
			   dens_phot_portal->matrix[i][j][k] = 1.07; // this is lung
			 }
		     }
		 }
	     }
	   // mass attenuation --> linear attenuation
	   dens_phot_portal->matrix[i][j][k] *= rho;

	 }  // for (register int k=0; k<dim.z; ++k)
       }  // for (register int j=0; j<dim.y; ++j)
     }  // for (register int i=0; i<dim.x; ++i)
   }
// --- End of Adding ---------------------------------------------

  return true;
}

