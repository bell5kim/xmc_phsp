/*****************************************************************************
 * p_beam_1point_spec.cpp:                                                   *
 *    class member functions for:                                            *
 *       p_beam_1point_spec: point source photon beam, energy spectrum       *
 *                                                                           *
 * Copyright (C) 2002    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 31.01.2002      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <sstream>
#include "p_beam_1point_spec.h"

// ***********************************************
// member functions of class p_beam_1point_spec
// ***********************************************

// get base data from file, initialize energy spectrum
void p_beam_1point_spec::get_base_data()
{
   // read lines
   char line[81]  = "";    // lines to read from base data file
   bool read_line = true;
   while (read_line)
   {
      if (base_file.eof())
      {
         xvmc_error("p_beam_1point_spec::get_base_data",
                    "end of file reached, base data incomplete",8);
      }
      base_file.getline(line,sizeof(line));
      istringstream line_stream(line);
      char keyword[81] = "";
      line_stream >> keyword;
      if (!strcmp(keyword,"END-PARAMETERS"))
      {
         xvmc_error("p_beam_1point_spec::get_base_data",
                    "end of parameter entry, base data incomplete",8);
      }
      else
      {
         if (!strcmp(keyword,"ENERGY-PD:"))
         {
            line_stream >> n_bin >> energy_min;

            // allocate memory for the spectrum arrays
            if ( (energy_bin  = new real[n_bin]) == NULL )
            {
               xvmc_error("p_beam_1point_spec::get_base_data",
                  "cannot allocate memory for energy bin array",8);
            }
            if ( (dprob  = new real[n_bin]) == NULL )
            {
               xvmc_error("p_beam_1point_spec::get_base_data",
                  "cannot allocate memory for differential distribution",8);
            }
            if ( (cprob  = new real[n_bin]) == NULL )
            {
               xvmc_error("p_beam_1point_spec::get_base_data",
                  "cannot allocate memory for cumulative distribution",8);
            }

            // input spectral data
            char comma;
            norm_factor = ZERO;
            for (register int i=0; i<n_bin; ++i)
            {
               base_file.getline(line,sizeof(line));
               istringstream line_stream(line);
               line_stream >> energy_bin[i] >> comma >> dprob[i];
               norm_factor += dprob[i];
            }

            // end of spectrum data
            read_line = false;
         }
      } // if (!strcmp(keyword,"END-PARAMETERS"))

   } // while (read_line)

   // normalize and calculate cumulative energy distribution
   dprob[0]   /= norm_factor;
   cprob[0]    = dprob[0];
   energy_mean = (energy_min+energy_bin[0])*dprob[0]/TWO;
   for (register int i=1; i<n_bin; ++i)
   {
      dprob[i]    /= norm_factor;
      cprob[i]     = cprob[i-1]+dprob[i];
      energy_mean += (energy_bin[i-1]+energy_bin[i])*dprob[i]/TWO;
   }
   energy_max = energy_bin[n_bin-1];

   // print beam model parameters
   xvmc_message("n_bin:        ",n_bin,"",1);
   xvmc_message("energy_min:   ",energy_min,"MeV",0);
   xvmc_message("energy_mean:  ",energy_mean,"MeV",0);
   xvmc_message("energy_max:   ",energy_max,"MeV",0);
}
