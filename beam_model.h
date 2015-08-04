#ifndef _BEAM_MODEL_H_
#define _BEAM_MODEL_H_

/*****************************************************************************
 * beam_model.h:                                                             *
 *    class declarations for different accelerator head models:              *
 *       beam_model:          abstract base class for beam models            *
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

#include <fstream>
using namespace std;

#include "definitions.h"
#include "ranmar.h"
#include "sobseq.h"
#include "beam_core.h"

// abstract base class for different photon and electron beam models
class beam_model
{
  friend real  beam_area(const beam_model *, const int, real_3 *);

 public:

  // construct beam model from beam parameters
  beam_model(const beam_core *);

  beam_model(particle_type,
	     int model_id,
	     real nominalEnergy,
	     real iso_distance,
	     real norm_value,
	     real grayPerMuAtDmax) :
    type(type),
    model_id(model_id),
    nominal_energy(nominalEnergy),
    iso_distance(iso_distance),
    norm_value(norm_value),
    gray_mu_dmax(grayPerMuAtDmax)
    {}

  virtual ~beam_model() {}

  particle_type  get_type() const         { return type; }
  int            getModelID() const       { return model_id;}
  virtual real   getNominalEnergy() const { return nominal_energy; }
  virtual real   average_energy() const = 0;   // get average energy of the spectrum (MeV)
  virtual real   maximum_energy() const = 0;   // get maximum energy of the spectrum (MeV)
  real           get_iso_distance() const { return iso_distance; } // get source to iso-center distance in CM!!
  // set source to iso-center distance in CM!!
  //real set_iso_distance(real id) { return(iso_distance = id);}   // Gate to Hell
  // get dose normalization parameters                             // ... or Texas, whichever is worse.
  real           get_norm() const         { return norm_value; }
  real           get_gray() const         { return gray_mu_dmax; }

  // get particle parameters from the beam
  virtual bool   emit(particle_parameters &,
#ifdef USE_SOBOL
		      sobseq &, int &,
#endif
		      ranmar &) = 0;


 protected:

  particle_type type;               // PHOTON or ELECTRON beam
  int           model_id;           // beam model id
  real          nominal_energy;     // nominal beam energy (MeV)
  real          iso_distance;       // origin to iso-center distance (cm)
  // normalization value = maximum depth dose per initial fluence
  // in water for the reference beam
  real          norm_value;
  // Gray per monitor unit at Dmax of the reference beam in water
  real          gray_mu_dmax;

 public:

  // get beam angles
  real get_gantry_angle(void) const { return(gantry_angle); }
  real get_table_angle(void) const { return(table_angle); }
  // get beam direction cosine
  real_3 direction(void) { return(beam_dir); }
  // get beam origin position
  real_3 origin_position(void) { return(origin); }
#ifdef PHSP_READ
  // trace particle to the simulation grid (take beam angles into account)
  bool trace2cube(real_3 &, real_3 &, int_3 &, real, ranmar &);
#endif
 protected:

  beam_model() {}
      // base data file
      ifstream      base_file;

      real_3        iso_center;         // iso-center position (cm)
      bool          gantry_rotation;    // true for gantry rotation
      real          start_gantry_angle; // start gantry angle (rad)
      real          stop_gantry_angle;  // stop gantry angle for rotation (rad)
      real          delta_gantry_angle; // gantry angle difference for rot.
      real          gantry_angle;       // average gantry angle
      real          table_angle;        // table angle (rad)
      real          collimator_angle;   // collimator angle (rad)
      real          open_x1;            // beam limit X1 (cm, iso-center plane)
      real          open_x2;            // beam limit X2 (cm, iso-center plane)
      real          open_y1;            // beam limit Y1 (cm, iso-center plane)
      real          open_y2;            // beam limit Y2 (cm, iso-center plane)
      real          cos_alpha;          // cos(gantry_angle)
      real          sin_alpha;          // sin(gantry_angle)
      real          cos_beta;           // cos(table_angle)
      real          sin_beta;           // sin(table_angle)
      real          cos_gamma;          // cos(collimator_angle)
      real          sin_gamma;          // sin(collimator_angle)
      real_3        beam_dir;           // beam direction cosine
      real_3        origin;             // beam origin (target) position
      real          theta_x1;           // minimum collimator angle (X1)
      real          theta_x2;           // maximum collimator angle (X2)
      real          theta_y1;           // minimum collimator angle (Y1)
      real          theta_y2;           // maximum collimator angle (Y2)
      real          theta_width_x;      // theta_x2 - theta_x1
      real          theta_width_y;      // theta_y2 - theta_y1
#ifndef PHSP_READ
      // trace particle to the simulation grid (take beam angles into account)
      bool trace2cube(real_3 &, real_3 &, int_3 &, real, ranmar &);
#endif
      // number of collimator-patient collision warnings
      int  num_warnings;
      // open file and search for base data entry
      void find_base_data(const beam_core *);

      // close base data file
      void close_base_file(void) {
         base_file.close();
      }
};

#endif  /* _BEAM_MODEL_H_ */
