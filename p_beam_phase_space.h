#ifndef _P_BEAM_PHASE_SPACE_H_
#define _P_BEAM_PHASE_SPACE_H_

/*****************************************************************************
 * p_beam_phase_space.h:                                                     *
 *    class declarations and inline member functions for:                    *
 *       p_beam_phase_space: particle sources from BEAM Mode2 Phase Space    *
 *                                                                           *
 * Copyright (C) 2002    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 31.01.2002      *
 *    Modified from p_beam_2gauss_mono.h                  JK 07-31-2015      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include "p_beam_model.h"

// phase space beam model
class p_beam_phase_space : public p_beam_model
{
public:
    p_beam_phase_space(const beam_core *this_beam)
        : p_beam_model(this_beam)
    {
        find_base_data(this_beam);
        get_base_data();
        close_base_file();
    }

    // that was reduplicated
    p_beam_phase_space(const p_beam_phase_space &);
    virtual real average_energy(void) const
    {
        return(nominal_energy);
    }
    virtual real maximum_energy(void) const
    {
        return(nominal_energy);
    }

    // get particle parameters from the beam
    bool emit(particle_parameters &,
#ifdef USE_SOBOL
              sobseq &, int &,
#endif
              ranmar &);

    FILE *open(int i_process);

protected:
    real_3    iso_center;
    long      positronCount=0, negatronCount=0, photonCount=0;
    long      primaryCount=0, scatterCount=0;
    real      pETotal=0.0, eETotal=0.0;
    real      beamArea;
    // real      col_cdistance = 44.999;  // Default Position of Collimator

    int       nTotal, nPhoton;
    int       nPhsp = 0, nCycl = 0;
    float     mxPhotEng, mnElectEng, nSrcPart;


    // photon head scatter model parameters
    real      p_pri;         // probability for primary photons
    real      pri_distance;  // distance from z=0 to primary source (cm)
    real      sigma_pri;     // sigma of Gaussian for primary photons (cm)
    real      horn_0;        // horn correction parameter 0
    real      horn_1;        // horn correction parameter 1
    real      horn_2;        // horn correction parameter 2
    real      horn_3;        // horn correction parameter 3
    real      horn_4;        // horn correction parameter 4
    real      p_sct;         // head scatter probability (p_sct=1-p_pri)
    real      sct_distance;  // distance from z=0 to scatter source (cm)
    real      sigma_sct;     // sigma of Gaussian for scattered photons (cm)

    // parameters for charged particle contamination,
    real      p_con;         // probability for electron contamination
    real      distance_con;  // distance from z=0 to electron source (cm)
    real      radius_con;    // radius of the electron source (cm)
    real      e_mean_con;    // average energy of charged particles

    // auxiliary parameters for the charged particle energy spectrum
    real      e1_con_aux,e2_con_aux;

    // further parameters from base data file needed by the derived class
    real      energy_min_aux,energy_max_aux,l_aux,b_aux;

    // NU Parameter for Off-axis Beam Hardening and Softening
    real      nu_value;      // nu value for off-axis spectrum

    // get base data from file
    void      get_base_data();

    // here "sample_energy" just returns the nominal energy,
    // however, derived beam models will overload this function
    virtual real sample_energy(ranmar &rndm)
    {
        return(nominal_energy);
    }
};

#endif  /* _P_BEAM_PHASE_SPACE_H_ */
