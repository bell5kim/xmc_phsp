/*****************************************************************************
 * p_beam_phase_space.cpp:                                                     *
 *    class declarations and inline member functions for:                    *
 *       p_beam_phase_space: particle sources from BEAM Mode2 Phase Space    *
 *                                                                           *
 * Copyright (C) 2002    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 31.01.2002      *
 *    Modified from p_beam_2gauss_mono.cpp                JK Feb, 1, 2008    *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************
#include <iostream>
#include <sstream>

#include "p_beam_phase_space.h"

// ***********************************************
// member functions of class p_beam_phase_space
// ***********************************************

// base data file input
void p_beam_phase_space::get_base_data()
{
    // define bool variables
    bool p_pri_found         = false;
    bool pri_distance_found  = false;
    bool sigma_pri_found     = false;
    bool horn_0_found        = false;
    bool horn_1_found        = false;
    bool horn_2_found        = false;
    bool horn_3_found        = false;
    bool horn_4_found        = false;
    bool sct_distance_found  = false;
    bool sigma_sct_found     = false;
    bool col_mdistance_found = false;
    bool col_cdistance_found = false;
    bool col_xdistance_found = false;
    bool col_ydistance_found = false;
    bool norm_value_found    = false;
    bool gray_mu_dmax_found  = false;
    bool energy_min_found    = false;
    bool energy_max_found    = false;
    bool l_value_found       = false;
    bool b_value_found       = false;
    bool p_con_found         = false;
    bool distance_con_found  = false;
    bool radius_con_found    = false;
    bool e_mean_con_found    = false;
#ifdef USE_NU
    bool nu_value_found      = false;
    nu_value = 0.45;         // Defualt Value
#endif

    // read lines
    char line[81]  = "";    // lines to read from base data file
    bool read_line = true;
    while (read_line)
    {
        if (base_file.eof())
        {
            xvmc_error("p_beam_phase_space::get_base_data",
                       "end of file reached, base data incomplete",8);
        }
        base_file.getline(line,sizeof(line));
        istringstream line_stream(line);
        char keyword[81] = "";
        line_stream >> keyword;
        if (!strcmp(keyword,"END-PARAMETERS"))
        {
            xvmc_error("p_beam_phase_space::get_base_data",
                       "end of parameter entry, base data incomplete",8);
        }
        else
        {
            if (!strcmp(keyword,"PRIMARY-PHOTONS:"))
            {
                line_stream >> p_pri;
                p_sct = ONE-p_pri;
                p_pri_found = true;
            }
            if (!strcmp(keyword,"PRIMARY-DIST:"))
            {
                line_stream >> pri_distance;
                pri_distance_found = true;
            }
            if (!strcmp(keyword,"PRIMARY-SIGMA:"))
            {
                line_stream >> sigma_pri;
                sigma_pri_found = true;
            }
            if (!strcmp(keyword,"PRIMARY-HORN0:"))
            {
                line_stream >> horn_0;
                horn_0_found = true;
            }
            if (!strcmp(keyword,"PRIMARY-HORN1:"))
            {
                line_stream >> horn_1;
                horn_1_found = true;
            }
            if (!strcmp(keyword,"PRIMARY-HORN2:"))
            {
                line_stream >> horn_2;
                horn_2_found = true;
            }
            if (!strcmp(keyword,"PRIMARY-HORN3:"))
            {
                line_stream >> horn_3;
                horn_3_found = true;
            }
            if (!strcmp(keyword,"PRIMARY-HORN4:"))
            {
                line_stream >> horn_4;
                horn_4_found = true;
            }
            if (!strcmp(keyword,"SCATTER-DIST:"))
            {
                line_stream >> sct_distance;
                sct_distance_found = true;
            }
            if (!strcmp(keyword,"SCATTER-SIGMA:"))
            {
                line_stream >> sigma_sct;
                sigma_sct_found = true;
            }
            if (!strcmp(keyword,"COLM-DIST:"))
            {
                line_stream >> col_mdistance;
                col_mdistance_found = true;
            }
            if (!strcmp(keyword,"COLC-DIST:"))
            {
                line_stream >> col_cdistance;
                col_cdistance_found = true;
            }
            if (!strcmp(keyword,"COLX-DIST:"))
            {
                line_stream >> col_xdistance;
                col_xdistance_found = true;
            }
            if (!strcmp(keyword,"COLY-DIST:"))
            {
                line_stream >> col_ydistance;
                col_ydistance_found = true;
            }
            if (!strcmp(keyword,"NORM-VALUE:"))
            {
                line_stream >> norm_value;
                norm_value_found = true;
            }
            if (!strcmp(keyword,"GY/MU-DMAX:"))
            {
                line_stream >> gray_mu_dmax;
                gray_mu_dmax_found = true;
            }
            if (!strcmp(keyword,"ENERGY-MIN:"))
            {
                line_stream >> energy_min_aux;
                energy_min_found = true;
            }
            if (!strcmp(keyword,"ENERGY-MAX:"))
            {
                line_stream >> energy_max_aux;
                energy_max_found = true;
            }
            if (!strcmp(keyword,"L-VALUE:"))
            {
                line_stream >> l_aux;
                l_value_found = true;
            }
            if (!strcmp(keyword,"B-VALUE:"))
            {
                line_stream >> b_aux;
                b_value_found = true;
            }
            if (!strcmp(keyword,"CHARGED-PARTICLES:"))
            {
                line_stream >> p_con;
                p_con_found = true;
            }
            if (!strcmp(keyword,"CHARGED-DIST:"))
            {
                line_stream >> distance_con;
                distance_con_found = true;
            }
            if (!strcmp(keyword,"CHARGED-RADIUS:"))
            {
                line_stream >> radius_con;
                radius_con_found = true;
            }
            if (!strcmp(keyword,"CHARGED-E-MEAN:"))
            {
                line_stream >> e_mean_con;
                e_mean_con_found = true;
            }
#ifdef USE_NU
            if (!strcmp(keyword,"NU-VALUE:"))
            {
                line_stream >> nu_value;
                nu_value_found = true;
            }
#endif
        } // if (!strcmp(keyword,"END-PARAMETERS"))

        if ( p_pri_found         &&
                pri_distance_found  &&
                sigma_pri_found     &&
                horn_0_found        &&
                horn_1_found        &&
                horn_2_found        &&
                horn_3_found        &&
                horn_4_found        &&
                sct_distance_found  &&
                sigma_sct_found     &&
                col_mdistance_found &&
                col_cdistance_found &&
                col_xdistance_found &&
                col_ydistance_found &&
                norm_value_found    &&
                gray_mu_dmax_found  &&
                energy_min_found    &&
                energy_max_found    &&
                l_value_found       &&
                b_value_found       &&
                p_con_found         &&
                distance_con_found  &&
                radius_con_found    &&
#ifdef USE_NU_LATER
                nu_value_found      &&
#endif
                e_mean_con_found       ) read_line = false;

    } // while (read_line)

    // parameter input was successful, now reset collimator openings
    col_x1       = open_x1*col_xdistance/iso_distance;
    col_x2       = open_x2*col_xdistance/iso_distance;
    col_y1       = open_y1*col_ydistance/iso_distance;
    col_y2       = open_y2*col_ydistance/iso_distance;
    col_width_x  = col_x2 - col_x1;
    col_width_y  = col_y2 - col_y1;

    // print beam model parameters
    xvmc_message("p_pri:        ",p_pri*100.0,"%",1);
    xvmc_message("pri_distance: ",pri_distance,"cm",0);
    xvmc_message("sigma_pri:    ",sigma_pri,"cm",0);
    xvmc_message("horn_0:       ",horn_0,"",0);
    xvmc_message("horn_1:       ",horn_1,"",0);
    xvmc_message("horn_2:       ",horn_2,"",0);
    xvmc_message("horn_3:       ",horn_3,"",0);
    xvmc_message("horn_4:       ",horn_4,"",0);
    xvmc_message("p_sct:        ",p_sct*100.0,"%",0);
    xvmc_message("sct_distance: ",sct_distance,"cm",0);
    xvmc_message("sigma_sct:    ",sigma_sct,"cm",0);
    xvmc_message("col_mdistance:",col_mdistance,"cm",0);
    xvmc_message("col_cdistance:",col_cdistance,"cm",0);
    xvmc_message("col_xdistance:",col_xdistance,"cm",0);
    xvmc_message("col_ydistance:",col_ydistance,"cm",0);
    xvmc_message("norm_value:   ",norm_value,"",0);
    xvmc_message("gray_mu_dmax: ",gray_mu_dmax,"Gy",0);
    xvmc_message("energy_min:   ",energy_min_aux,"MeV",0);
    xvmc_message("energy_max:   ",energy_max_aux,"MeV",0);
    xvmc_message("l:            ",l_aux,"",0);
    xvmc_message("b:            ",b_aux,"",0);
    xvmc_message("p_con:        ",p_con*100.0,"%",0);
    xvmc_message("distance_con: ",distance_con,"cm",0);
    xvmc_message("radius_con:   ",radius_con,"cm",0);
    xvmc_message("e_mean_con:   ",e_mean_con,"MeV",0);
#ifdef USE_NU
    xvmc_message("nu_value:     ",nu_value,"",0);
#endif
}

// get particle parameters from phase space
bool p_beam_phase_space::emit(particle_parameters &p,
#ifdef USE_SOBOL
                              sobseq &sobol, int &sobol_count,
#endif
                              ranmar &rndm)
{
    real origin_point_dist; // origin (x=y=z=0) to point distance
    // (position in X-collimator plane)
    real rx,ry,rp,rr,rphi;  // auxiliary random (Sobol) numbers
    real src_point_dist;    // distance from particle position in the source
    // plane to position in collimator plane (point)
    real rad_src;           // particle position in source plane (radius)
    real phi_src;           // particle position in source plane (angle)
    real x_src,y_src,z_src; // particle position in source plane (x,y,z)

    real energy_0;          // photon energy before the softening correction

    real rad_pri;           // primary photon position in target plane (radius)
    real phi_pri;           // primary photon position in target plane (angle)
    real x_pri,y_pri,z_pri; // primary photon direction before filter scatter
    real pri_sct_distance;  // distance from primary photon position to the
    // photon position in the head scatter plane
    // (simulate Compton)
    real cos_comp;          // Compton scattering angle
    real fac_comp;          // Compton correction factor

    // first of all, we check for charged particle (electron) contamination

    // random (Sobol) number to determine the source particle type
#ifdef USE_SOBOL
    rp = sobol.number(sobol_count);
    ++sobol_count;
#else
    rp = rndm.number();
#endif
    // this is a photon, set photon parameters
    p.type = PHOTON;
    /* simulate primary (target) photons */

    // sample photon position in the MC starting plane
    // (above the collimators)
    p.pos.z = col_mdistance;

    // X-position in the X-collimator plane
#ifndef MOVING_JAW
    p.pos.x = col_x1 + col_width_x*rx;
#else
    p.pos.x = col_x1 + col_width_x - col_width_x*rx*rndm.number();
#endif
    // X-position in the MC starting plane (above the collimators)
    p.pos.x =
        (p.pos.x-x_src)*(col_mdistance-z_src)/(col_xdistance-z_src)+x_src;

    // Y-position in the Y-collimator plane
    p.pos.y = col_y1 + col_width_y*ry;
    // Y-position in the MC starting plane (above the collimators)
    p.pos.y =
        (p.pos.y-y_src)*(col_mdistance-z_src)/(col_ydistance-z_src)+y_src;

    // origin to point distance
    origin_point_dist =
        sqrt(p.pos.x*p.pos.x + p.pos.y*p.pos.y + p.pos.z*p.pos.z);

    // direction
    p.dir.x = p.pos.x - x_src;
    p.dir.y = p.pos.y - y_src;
    p.dir.z = p.pos.z - z_src;
    src_point_dist =
        sqrt(p.dir.x*p.dir.x + p.dir.y*p.dir.y + p.dir.z*p.dir.z);
    p.dir.x = p.dir.x/src_point_dist;
    p.dir.y = p.dir.y/src_point_dist;
    p.dir.z = p.dir.z/src_point_dist;

    // Moved by JOKim on Jan 14, 2008
    // sample energy and correct due to off-axis softening, i.e.
    // in the mono-energetic case the beam is not really mono-energetic
    energy_0 = sample_energy(rndm);

    // p.weight  = mu_cor2*cor_horn;
    p.weight  = 1.0;

    if ( modifier->transport(p,rndm) )
    {
        // trace particle to the simulation grid (take beam angles into account)
        // return true if the particle hits the calculation cube
        return( trace2cube(p.pos, p.dir, p.i, origin_point_dist, rndm) );
    }

// there isn't any particle to emit, return false
    return(false);
}

p_beam_phase_space::p_beam_phase_space(const p_beam_phase_space &A)
{
    // beam_model variables
    iso_distance = A.iso_distance;
    type = A.type;
    nominal_energy = A.nominal_energy;
    model_id = A.model_id;
    // p_beam_model variables
    col_mdistance = A.col_mdistance;
    col_cdistance = A.col_cdistance;
    col_xdistance = A.col_xdistance;
    col_ydistance = A.col_ydistance;
    col_x1 = A.col_x1;
    col_x2 = A.col_x2;
    col_y1 = A.col_y1;
    col_y2 = A.col_y2;
    col_width_x = A.col_width_x;
    col_width_y = A.col_width_y;
    //p_beam_phase_space variables
    p_pri = A.p_pri;
    pri_distance = A.pri_distance;
    sigma_pri = A.sigma_pri;
    horn_0 = A.horn_0;
    horn_1 = A.horn_1;
    horn_2 = A.horn_2;
    horn_3 = A.horn_3;
    horn_4 = A.horn_4;
    p_sct = A.p_sct;
    sct_distance = A.sct_distance;
    sigma_sct = A.sigma_sct;
    norm_value = A.norm_value;
    gray_mu_dmax = A.gray_mu_dmax;
    energy_min_aux = A.energy_min_aux;
    energy_max_aux = A.energy_max_aux;
    l_aux = A.l_aux;
    b_aux = A.b_aux;
#ifdef USE_NU
    nu_value = A.nu_value;
#endif
    p_con = A.p_con;
    distance_con = A.distance_con;
    radius_con = A.radius_con;
    e_mean_con = A.e_mean_con;
}

FILE* p_beam_phase_space::open(int i_process)
{
    char sPhspFile[256];
    sprintf(sPhspFile, "phspRead%d.bin", i_process);

    printf ("PHSP_READ> Open Phase Space File at %f cm: %s for CPU %d\n", col_cdistance, sPhspFile, i_process);

    FILE *istrm = fopen(sPhspFile,"rb"); /* input stream */
    if(istrm==NULL)
    {
        printf("\n ERROR: opening input %s \n", sPhspFile);
        exit(-1);
    }
    rewind (istrm);

    char ModeRW[20];
    if (fread(ModeRW, sizeof(char), 5, istrm) == 5)
    {
        ModeRW[5] = '\0';
    }
    else
    {
        printf ("ERROR: Can not read MODE_RW\n");
        exit(-1);
    }

    int intArray[10];
    if(fread(intArray, sizeof(int), 2, istrm)!=2)
    {
        printf("ERROR: reading intArray from header");
        exit(-1);
    }

    nTotal  = intArray[0];
    nPhoton = intArray[1];

    float floatArray[20];
    if (fread(floatArray, sizeof(float), 3, istrm)!=3)
    {
        printf("ERROR: reading floatArray from header");
        exit(-1);
    }

    mxPhotEng  = floatArray[0];
    mnElectEng = floatArray[1];
    nSrcPart   = floatArray[2];

    // for Temporary Debug
    printf ("PHSP_READ> nTotal = %d nPhoton= %d mxPhotEng = %E  mnElectEng = %E  nSrcPart = %E\n",
            nTotal, nPhoton, mxPhotEng,mnElectEng,nSrcPart);

    /* padding */
    char charArray[200];
    if (fread(charArray, sizeof(char), 7, istrm)!=7)
    {
        printf("ERROR: reading mode from header\n");
        exit(-1);
    }
    charArray[7] = '\0';

    return (istrm);
}
