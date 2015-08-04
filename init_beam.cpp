/*****************************************************************************
 * init_beam.cpp:                                                            *
 *    function init_beam: initialize beam data                               *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 2000/01/12      *
 *    new beam model classes and class names              MF 2000/12/04      *
 *    added portal dose image                             MF 03/03/13        *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <stdlib.h>
#include <math.h>
#include "definitions.h"
#include "global.h"
#include "treatment_plan.h"
#include "beam_model.h"
#include "e_beam_model.h"
#include "e_beam_1point_mono.h"
#include "e_beam_1point_poly.h"
#include "e_beam_triple_mono.h"
#include "e_beam_triple_poly.h"
#include "e_beam_bumpy_mono.h"
#include "e_beam_bumpy_poly.h"
#include "p_beam_model.h"
#include "p_beam_1point_mono.h"
#include "p_beam_1point_spec.h"
#include "p_beam_2gauss_mono.h"
#include "p_beam_2gauss_poly.h"
#include "irregular_modifier.h"
#include "simple_mlc.h"
#include "real_mlc.h"
#include "electron_data.h"
#include "photon_data.h"
#include "electron_data_inp.h"
#include "photon_data_inp.h"

// ****************************************
// declare functions and global variables
// ****************************************

// treatment plan parameters
extern treatment_plan plan;

// input arrays for electron and photon transport parameters
extern electron_data_inp *e_inp;
extern photon_data_inp   *p_inp;

// arrays of electron and photon transport parameters with high resolution
// for linear interpolation during the simulation
electron_data *e_h2o = NULL;
photon_data   *p_h2o = NULL;

// pointer to the beam model(s)
beam_model    *linac  = NULL;
e_beam_model  *elinac = NULL;
p_beam_model  *plinac = NULL;

// pointer to the photon and electron beam modifiers
beam_modifier *pbeam_modifier = NULL;
beam_modifier *ebeam_modifier = NULL;

// calculate intersection plane for central beam axis
int intersection_plane(const real_3 &, const real_3 &, real_3 &);

// calculate surface area and edge points of the beam
real beam_area(const beam_model *, const int, real_3 *);

// ****************************************
// function init_beam
// ****************************************

void init_beam(beam_core *this_beam)
{
   real_3   beam_dir;        // central beam axis direction cosine
   real_3   origin;          // beam origin position 
   int      i_plane;         // intersection plane for central beam axis
   real_3   i_point;         // intersection point
   real     surface_area;    // beam surface area
   real_3   edges[4];        // edge points

   xvmc_message("======================",1);
   xvmc_message("= Initialize beam:",this_beam->id,"",0);
   xvmc_message("======================",0);

   // define and initialize the accelerator head model
   if (linac != NULL)
   {
      delete linac;
      linac  = NULL;
      elinac = NULL;
      plinac = NULL;
   }

   switch (this_beam->type)
   {
   case PHOTON:
      switch (this_beam->model_id)
      {
      case  0: // mono-energetic point source beam
         if ( (plinac = new p_beam_1point_mono(this_beam)) == NULL )
         {
            xvmc_error("init_beam",
                       "cannot create mono-energetic photon point source",8);
         }
         break;
      case -1: // mono-energetic beam with two Gaussian sources
         if ( (plinac = new p_beam_2gauss_mono(this_beam)) == NULL )
         {
            xvmc_error("init_beam",
                       "cannot create mono-energetic photon 2-Gauss source",8);
         }
         break;
      case  1: // poly-energetic beam with two Gaussian sources
         if ( (plinac = new p_beam_2gauss_poly(this_beam)) == NULL )
         {
            xvmc_error("init_beam",
                       "cannot create poly-energetic photon 2-Gauss source",8);
         }
         break;
      case 21: // point source photon beam, energy spectrum
         if ( (plinac = new p_beam_1point_spec(this_beam)) == NULL )
         {
            xvmc_error("init_beam",
                       "cannot create spectral photon point source",8);
         }
         break;
/* Added by JK on 03/08/2015 -- To do for phase space 
      case 31: // point source photon beam, energy spectrum
         if ( (plinac = new phase_space_read(this_beam)) == NULL )
         {
            xvmc_error("init_beam",
                       "cannot read phase space source",8);
         }
         break;
*/
      default:
         xvmc_error("init_beam","beam model not implemented",8);
         break;
      }
      linac = plinac;
      break;
   case ELECTRON:
      switch (this_beam->model_id)
      {
      case  0: // mono-energetic point source beam
         if ( (elinac = new e_beam_1point_mono(this_beam)) == NULL )
         {
            xvmc_error("init_beam",
                       "cannot create mono-energetic electron point source",8);
         }
         break;
      case  1: // poly-energetic point source beam
         if ( (elinac = new e_beam_1point_poly(this_beam)) == NULL )
         {
            xvmc_error("init_beam",
                       "cannot create poly-energetic electron point source",8);
         }
         break;
      case -3: // triple electron source model, mono-energetic
         if ( (elinac = new e_beam_triple_mono(this_beam)) == NULL )
         {
            xvmc_error("init_beam",
                 "cannot create mono-energetic triple electron source",8);
         }
         break;
      case  3: // triple electron source model, poly-energetic
         if ( (elinac = new e_beam_triple_poly(this_beam)) == NULL )
         {
            xvmc_error("init_beam",
                 "cannot create poly-energetic triple electron source",8);
         }
         break;
      case -4: // triple electron source model, mono-energetic
               // with bump correction
         if ( (elinac = new e_beam_bumpy_mono(this_beam)) == NULL )
         {
            xvmc_error("init_beam",
                 "cannot create mono-energetic bumpy electron source",8);
         }
         break;
      case  4: // triple electron source model, poly-energetic
               // with bump correction
         if ( (elinac = new e_beam_bumpy_poly(this_beam)) == NULL )
         {
            xvmc_error("init_beam",
                 "cannot create poly-energetic bumpy electron source",8);
         }
         break;
      default:
         xvmc_error("init_beam","beam model not implemented",8);
         break;
      }
      linac = elinac;
      break;
   default:
      xvmc_error("init_beam","particle beam type not implemented",8);
      break;
   }

   // beam modifiers for photon beams
   if (plinac != NULL)
   {
      // define and initialize irregular photon beam modifier
      if ( this_beam->irregular != NULL)
      {
         if (pbeam_modifier == NULL)
         {
            // delete empty beam modifier
            if (plinac->modifier != NULL) delete plinac->modifier;

            // there is no beam modifier, create irregular beam modifier
            // based on the irregular beam contour at the iso-center plane
            if ( (pbeam_modifier =
               new irregular_modifier(this_beam->irregular,
                                      plinac->get_iso_distance())) == NULL )
            {
               xvmc_error("init_beam",
                          "cannot create irregular beam modifier",8);
            }
         }
         else
         {
            // there is already a modifier, delete old beam modifier
            delete pbeam_modifier; pbeam_modifier = NULL;

            // create irregular beam modifier based on the irregular beam
            // contour at the iso-center plane
            if ( (pbeam_modifier =
               new irregular_modifier(this_beam->irregular,
                                      plinac->get_iso_distance())) == NULL )
            {
               xvmc_error("init_beam",
                          "cannot create new irregular beam modifier",8);
            }
         }

         // assign modifier pointers
         plinac->modifier = pbeam_modifier;
      }

      // define and initialize multi-leaf collimator
      if ( this_beam->mlc != NULL)
      {
         if (pbeam_modifier == NULL)
         {
            // delete empty beam modifier
            if (plinac->modifier != NULL) delete plinac->modifier;

            // there is no beam modifier, create MLC based on the leaf
            // positions at the iso-center plane
            switch (this_beam->mlc->get_type())
            {
            case DBLFOCUS_MLC:
            case RNDFOCUS_MLC:
            case ELEKTA_MLC:
            case VARIAN_MLC:
               xvmc_message("Creating real MLC",1);
               if ( (pbeam_modifier=new real_mlc(this_beam->mlc,
                                        plinac->get_iso_distance())) == NULL )
               {
                  xvmc_error("init_beam",
                             "cannot create real 3D Monte Carlo MLC model",8);
               }
               break;
            default:
               xvmc_message("Creating simple MLC",1);
               if ( (pbeam_modifier=new simple_mlc(this_beam->mlc,
                                        plinac->get_iso_distance())) == NULL )
               {
                  xvmc_error("init_beam",
                             "cannot create simple MLC",8);
               }
               break;
            }
         }
         else
         {
            // there is already a modifier, compare beam modifier types
            switch (this_beam->mlc->get_type())
            {
            case DBLFOCUS_MLC:
            case RNDFOCUS_MLC:
            case ELEKTA_MLC:
            case VARIAN_MLC:
               if (pbeam_modifier->get_modifier_type() == MODIFIER_REAL_MLC)
               {
                  // the beam modifier types coincide, change leaf positions
                  real_mlc *r_mlc = (real_mlc *) pbeam_modifier;
                  xvmc_message("Setting mode of real MLC",1);
                  r_mlc->set_mode(this_beam->mlc->get_mode());
                  xvmc_message("Changing leaf positions of real MLC",1);
                  r_mlc->change_leaf_positions(this_beam->mlc);
               }
               else
               {
                  // inconsistent beam modifier types, delete old beam modifier
                  delete pbeam_modifier; pbeam_modifier = NULL;

                  // create new MLC based on the leaf positions at the
                  // iso-center plane
                  xvmc_message("Creating new real MLC",1);
                  if ( (pbeam_modifier=new real_mlc(this_beam->mlc,
                                plinac->get_iso_distance())) == NULL )
                  {
                     xvmc_error("init_beam",
                         "cannot create new real 3D Monte Carlo MLC model",8);
                  }
               }
               break;
            default:
               if (pbeam_modifier->get_modifier_type() == MODIFIER_SIMPLE_MLC)
               {
                  // the beam modifier types coincide, change leaf positions
                  simple_mlc *s_mlc = (simple_mlc *) pbeam_modifier;
                  xvmc_message("Setting mode of simple MLC",1);
                  s_mlc->set_mode(this_beam->mlc->get_mode());
                  xvmc_message("Changing leaf positions of simple MLC",1);
                  s_mlc->change_leaf_positions(this_beam->mlc);
               }
               else
               {
                  // inconsistent beam modifier types, delete old beam modifier
                  delete pbeam_modifier; pbeam_modifier = NULL;

                  // create new MLC based on the leaf positions at the
                  // iso-center plane
                  xvmc_message("Creating new simple MLC",1);
                  if ( (pbeam_modifier=new simple_mlc(this_beam->mlc,
                                plinac->get_iso_distance())) == NULL )
                  {
                     xvmc_error("init_beam",
                                "cannot create new simple MLC",8);
                  }
               }
               break;
            }
         }

         // assign modifier pointers
         plinac->modifier = pbeam_modifier;
      }
   }

   // beam modifiers for electron beams
   if (elinac != NULL)
   {
      // define and initialize irregular electron beam modifier
      if ( this_beam->irregular != NULL)
      {
         // source to applicator distance
         real SAD = elinac->get_app_distance();

         if (ebeam_modifier == NULL)
         {
            // delete empty beam modifier
            if (elinac->modifier != NULL) delete elinac->modifier;

            // there is no beam modifier, create irregular beam modifier
            // based on the irregular beam contour at the iso-center plane
            if ( (ebeam_modifier =
               new irregular_modifier(this_beam->irregular,
                                      elinac->get_iso_distance(),SAD))==NULL )
            {
               xvmc_error("init_beam",
                          "cannot create irregular beam modifier",8);
            }
         }
         else
         {
            // there is already a modifier, delete old beam modifier
            delete ebeam_modifier; ebeam_modifier = NULL;

            // create irregular beam modifier based on the irregular beam
            // contour at the iso-center plane
            if ( (ebeam_modifier =
               new irregular_modifier(this_beam->irregular,
                                      elinac->get_iso_distance(),SAD))==NULL )
            {
               xvmc_error("init_beam",
                          "cannot create new irregular beam modifier",8);
            }
         }

         // assign modifier pointers
         elinac->modifier = ebeam_modifier;
      }
   }

   // direction cosine of the central beam axis
   beam_dir = linac->direction();

   // beam origin position
   origin = linac->origin_position();

   // move the portal dose plane to the correct place
   if (plan.portal != NULL) plan.portal->move(origin,linac->get_gantry_angle(),
                                                     linac->get_table_angle());

   // assign pointer to the portal dose image
   this_beam->portal = plan.portal;

   // fix intersection plane of central beam axis with calculation cube
   i_plane = intersection_plane(origin,beam_dir,i_point);
   if ((i_plane < 1) || (i_plane > 6))
   {
      xvmc_error("init_beam","beam axis doesn't touch calculation cube",8);
   }

   // calculate surface area and edge points of the beam
   surface_area = beam_area(linac,i_plane,edges);

   xvmc_message("Average beam area on the calculation cube surface:",
                surface_area,"cm^2",1);

   real eff_beam_width = sqrt(surface_area);  // effective beam width
   xvmc_message("  Effective beam width:",eff_beam_width,"cm",0);

   // if n_history is less than zero estimate history number depending on
   // the desired statistical accuracy given by (-n_history/10)%
   if (this_beam->n_history < 0)
   {
      // temporary variables for history and repetition numbers
      real  tmp_history = float(this_beam->n_history);
      real  tmp_repeat  = float(this_beam->n_repeat);
      // final sigma (standard deviation) squared
      const real sigma2 = tmp_history*tmp_history/100.0;
      real voxel_area = ZERO;       // voxel area depending on beam direction
      real voxel_size_depth = ZERO; // voxel size in beam direction
      real cos_depth = ZERO;        // beam depth direction cosine
      real factor = ZERO;           // additional factor

      // calculate voxel_area depending on beam direction
      switch (i_plane)
      {
      case 1:
      case 6: // Z
         voxel_area = voxel_size.x * voxel_size.y;
         voxel_size_depth = voxel_size.z;
         cos_depth=fabs(beam_dir.z);
         break;
      case 2:
      case 3: // X
         voxel_area = voxel_size.y * voxel_size.z;
         voxel_size_depth = voxel_size.x;
         cos_depth=fabs(beam_dir.x);
         break;
      case 4:
      case 5: // Y
         voxel_area = voxel_size.x * voxel_size.z;
         voxel_size_depth = voxel_size.y;
         cos_depth=fabs(beam_dir.y);
         break;
      default:
         xvmc_error("init_beam","beam axis doesn't touch calculation cube",16);
      }

      // estimate number of histories and repetitions
      if (this_beam->type == ELECTRON)
      // electron beam
      {
         tmp_repeat  = surface_area/5.0 + 0.5;
         tmp_history = 50000.0*cos_depth/voxel_area/sigma2;
      }
      else
      // photon beam
      {
         if (voxel_size_depth <= 0.333) {
            tmp_repeat = 75.0; }
         else
         {
            if (voxel_size_depth <= 2.5) {
               tmp_repeat = 25.0/voxel_size_depth; }
            else
            {
               tmp_repeat = 10.0;
            }
         }
         factor = surface_area/5.0 + 0.5;
         factor = factor*50000.0;
         tmp_history = factor*cos_depth/voxel_area/sigma2;
      }

      // print out history number
      this_beam->n_history = int(tmp_history);
      xvmc_message("Estimated number of particle histories:",
                    this_beam->n_history,"",1);

      // print out repetition number
      if (tmp_repeat < 1.0) this_beam->n_repeat = 1;
      else                  this_beam->n_repeat = int(tmp_repeat);
      xvmc_message("  repetitions:",
                    this_beam->n_repeat,"",0);

      // estimate additional repetition number for gantry rotation
      if (this_beam->gantry_rotation)
      {
         this_beam->n_rotate = 1
              + int(fabs(this_beam->stop_gantry_angle
                       - this_beam->start_gantry_angle)/90.0);
         xvmc_message("  additional repetitions for gantry rotation:",
                       this_beam->n_rotate,"",0);
      }

   } // end of history number estimate

   // initialize dose and dose error array for this beam
   // *beam_dose  = ZERO;
   // *beam_error = ZERO;
    memset((void *) beam_dose->data, 0, beam_dose->size*sizeof(float));
    memset((void *) beam_error->data, 0, beam_error->size*sizeof(float));
    
   // initialize portal dose and portal error matrices for this beam
   if (this_beam->portal != NULL) this_beam->portal->reset_beam_dose_portal();

   // initialize the arrays for linear interpolation
   // photon transport parameters
   if (p_h2o == NULL)
   {
      if ( (p_h2o = new photon_data) == NULL )
      {
         xvmc_error("init_beam","cannot create photon transport data",8);
      }
      photon_data_init(*p_h2o, 10000, p_inp, linac->maximum_energy());
   }
   else
   {
      xvmc_error("init_beam","photon transport data already created",8);
   }
   
   // electron transport parameters
   if (e_h2o == NULL)
   {
      if ( (e_h2o = new electron_data) == NULL )
      {
         xvmc_error("init_beam","cannot create electron transport data",8);
      }
      real maximum_energy = linac->maximum_energy();
      real average_energy = linac->average_energy();
				
		maximum_energy = 25;
		average_energy = 1;
				
      // for electron beams we initialize the electron transport data
      // using the maximum spectrum energy only
      if (linac->get_type() == ELECTRON) average_energy = maximum_energy;
      electron_data_init(*e_h2o, 32768, e_inp,
			  maximum_energy, average_energy,
			  plan.result, plan.film_factor);
/*		
		cout << "# " << maximum_energy << " " << average_energy << " " << e_step << " " << e_cut << endl;
		cout << "# " << voxel_size.x << " " << voxel_size.y << " " << voxel_size.z << endl;
		cout << "# " << e_h2o->inverse_delta_e << " " << e_h2o->sigma_max << " " << endl;
		for (int i=0; i<10000; i++) {
			cout << i/e_h2o->inverse_delta_e << " " << e_h2o->s_res[i] << " " << e_h2o->s_cor[i] << " " << e_h2o->alpha_r[i]
					<< " " << e_h2o->max_loss[i] << " " << e_h2o->dose_cor[i]
					<< " " << e_h2o->p_moller[i] << " " << e_h2o->sigma_tot[i] << endl;			
		}
*/		
   }
   else
   {
      xvmc_error("init_beam","electron transport data already created",8);
   }
   
   // test energy conservation
#ifdef CHECK_ENERGY
   tot_energy.init();
#endif // CHECK_ENERGY

}

// **************************************************
// function beam_area
// calculate surface area and edge points of the beam
// **************************************************

real beam_area(const beam_model *linac, const int i_plane, real_3 *edges)
{
   real    tan_x[4],tan_y[4];     // corner points
   real    cos_alpha,sin_alpha;   // gantry angle
   real    cos_beta,sin_beta;     // table angle
   real    cos_gamma,sin_gamma;   // collimator angle
   real_3  dir;                   // direction cosine for corner points
   real_3  origin;                // beam origin position
   real    s1,s2;

   // point 1
   tan_x[0] = linac->open_x1/SSD_0;
   tan_y[0] = linac->open_y1/SSD_0;
   // point 2
   tan_x[1] = linac->open_x2/SSD_0;
   tan_y[1] = linac->open_y1/SSD_0;
   // point 3
   tan_x[2] = linac->open_x2/SSD_0;
   tan_y[2] = linac->open_y2/SSD_0;
   // point 4
   tan_x[3] = linac->open_x1/SSD_0;
   tan_y[3] = linac->open_y2/SSD_0;

   // beam angles
   cos_alpha = linac->cos_alpha;
   sin_alpha = linac->sin_alpha;
   cos_beta  = linac->cos_beta;
   sin_beta  = linac->sin_beta;
   cos_gamma = linac->cos_gamma;
   sin_gamma = linac->sin_gamma;

   // origin position
   origin.x = linac->origin.x;
   origin.y = linac->origin.y;
   origin.z = linac->origin.z;

   // loop over points
   for (register int i=0; i<4; ++i)
   {
      // calculate direction cosine for corner point
      real temp = sqrt(ONE + tan_x[i]*tan_x[i] + tan_y[i]*tan_y[i]);
      dir.x = tan_x[i]/temp;
      dir.y = tan_y[i]/temp;
      dir.z = ONE/temp;

      // rotate by the collimator angle
      temp  = dir.x*sin_gamma;
      dir.x = dir.x*cos_gamma + dir.y*sin_gamma;
      dir.y = -temp           + dir.y*cos_gamma;

      // rotate by the gantry angle
      temp  = dir.x*sin_alpha;
      dir.x = dir.x*cos_alpha - dir.z*sin_alpha;
      dir.z = temp            + dir.z*cos_alpha;

      // rotate by the table angle
      temp  = dir.x*sin_beta;
      dir.x = dir.x*cos_beta + dir.y*sin_beta;
      dir.y = -temp          + dir.y*cos_beta;

      // calculate corner points
      switch (i_plane)
      {
      case 1:
         edges[i].z =  ZERO;
         temp       = -origin.z/dir.z;
         edges[i].x =  origin.x + temp*dir.x;
         edges[i].y =  origin.y + temp*dir.y;
         break;
      case 6:
         edges[i].z =  cube_size.z;
         temp       =  (cube_size.z-origin.z)/dir.z;
         edges[i].x =  origin.x + temp*dir.x;
         edges[i].y =  origin.y + temp*dir.y;
         break;
      case 2:
         edges[i].x =  ZERO;
         temp       = -origin.x/dir.x;
         edges[i].y =  origin.y + temp*dir.y;
         edges[i].z =  origin.z + temp*dir.z;
         break;
      case 3:
         edges[i].x =  cube_size.x;
         temp       =  (cube_size.x-origin.x)/dir.x;
         edges[i].y =  origin.y + temp*dir.y;
         edges[i].z =  origin.z + temp*dir.z;
         break;
      case 4:
         edges[i].y =  ZERO;
         temp       = -origin.y/dir.y;
         edges[i].x =  origin.x + temp*dir.x;
         edges[i].z =  origin.z + temp*dir.z;
         break;
      case 5:
         edges[i].y =  cube_size.y;
         temp       =  (cube_size.y-origin.y)/dir.y;
         edges[i].x =  origin.x + temp*dir.x;
         edges[i].z =  origin.z + temp*dir.z;
         break;
      default:
         xvmc_error("beam_area","wrong plane number",8);
      } // switch (i_plane)
   } // for (register int i=0; i<4; ++i)

   // calculate the beam surface area
   switch (i_plane)
   {
   case 1:
   case 6: // planes perpendicular to z
      s1=fabs((edges[0].x-edges[1].x)*(edges[0].y-edges[2].y)-
              (edges[0].x-edges[2].x)*(edges[0].y-edges[1].y));
      s2=fabs((edges[0].x-edges[2].x)*(edges[0].y-edges[3].y)-
              (edges[0].x-edges[3].x)*(edges[0].y-edges[2].y));
      return((s1+s2)/TWO);
   case 2:
   case 3: // planes perpendicular to x
      s1=fabs((edges[0].y-edges[1].y)*(edges[0].z-edges[2].z)-
              (edges[0].y-edges[2].y)*(edges[0].z-edges[1].z));
      s2=fabs((edges[0].y-edges[2].y)*(edges[0].z-edges[3].z)-
              (edges[0].y-edges[3].y)*(edges[0].z-edges[2].z));
      return((s1+s2)/TWO);
   case 4:
   case 5: // planes perpendicular to y
      s1=fabs((edges[0].x-edges[1].x)*(edges[0].z-edges[2].z)-
              (edges[0].x-edges[2].x)*(edges[0].z-edges[1].z));
      s2=fabs((edges[0].x-edges[2].x)*(edges[0].z-edges[3].z)-
              (edges[0].x-edges[3].x)*(edges[0].z-edges[2].z));
      return((s1+s2)/TWO);
   default:
      xvmc_error("beam_area","wrong plane number",16);
   } // switch (i_plane)

   return(ZERO);
}

