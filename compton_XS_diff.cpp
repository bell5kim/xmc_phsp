/*****************************************************************************
 * compton_XS_diff.cpp:                                                      *
 *    class member functions for:                                            *
 *       compton_XS_diff: differential compton cross section                 *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 14.12.1999      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <fstream>
#include <string>
using namespace std;

#include "compton_XS_diff.h"

/*****************************************************************************
 * member functions for class compton_XS_diff:                               *
 *    differential compton cross section                                     *
 *****************************************************************************/

//##################################################################################
// NEW data container for compton_XS_diff
//##################################################################################

ComptonDiffData::ComptonDiffData(int nEin, int nBin, real energyMin, real energyMax) :
  ne_in(nEin),
  n_bin(nBin),
  energy_min(energyMin),
  energy_max(energyMax),
  x_data(NULL),
  x_matrix(NULL)
{
  x_data = new real[ne_in*n_bin];
  x_matrix = new real*[ne_in];
  // init pointer array for direct access
  int i;
  for (i=0; i<ne_in; ++i) {
    x_matrix[i] = x_data + i*n_bin;
  }
}

//##################################################################################

ComptonDiffData::~ComptonDiffData()
{
  delete[] x_data;
  delete[] x_matrix;
}

//##################################################################################

ComptonDiffData* ComptonFactory(const std::string& file_name)
{
  real energy_min, energy_max;
  int ne_in, n_bin;

  // open data file
  ifstream file(file_name.c_str(), ios::in);
  if (!file) {
    xvmc_error("ComptonDiffData::read","cannot open file",8); }

   // read minimum and maximum energy
   file >> energy_min >> energy_max;

   // read the number of bins
   file >> ne_in >> n_bin;

   ComptonDiffData* pcxs = new ComptonDiffData(ne_in, n_bin, energy_min, energy_max);

   // read data from file
   int i, j;
   for (i=0; i<ne_in; ++i) {
     for (j=0; j<n_bin; ++j) {
       file >> pcxs->x_matrix[i][j];
     }
   }

   // close file
   file.close();

   if(pcxs->invariant() == false) {
     delete pcxs;
     return NULL;
   }

   return pcxs;
}

//##################################################################################

compton_XS_diff::compton_XS_diff(const ComptonDiffData& A) :
  ne_in(A.ne_in),
  n_bin(A.n_bin),
  energy_min(A.energy_min),
  energy_max(A.energy_max),
  x_data(NULL),
  x_matrix(NULL)
{
  // calculate inverse bin size for initial photon energy (log scale)
  inv_delta = float(ne_in-1)/log(energy_max/energy_min);
  // calculate distance between first and last cross section data bin
  distance = float(n_bin - 1);

  x_data = new real[ne_in*n_bin];
  x_matrix = new real*[ne_in];
  // init pointer array for direct access
  int i, j;
  for (i=0; i<ne_in; ++i) {
    x_matrix[i] = x_data + i*n_bin;
  }
  for (i=0; i<ne_in; ++i) {
    for (j=0; j<n_bin; ++j) {
      x_matrix[i][j] = A.x_matrix[i][j];
    }
  }
}

//##################################################################################

compton_XS_diff::~compton_XS_diff(void)
{
  delete [] x_matrix;
  delete [] x_data;
}

//##################################################################################

// read data from file
void compton_XS_diff::read(char *file_name)
{
   // open data file
   ifstream file(file_name, ios::in);
   if (!file) {
      xvmc_error("compton_XS_diff::read","cannot open file",8); }

   // read minimum and maximum energy
   file >> energy_min >> energy_max;

   // read the number of bins
   file >> ne_in >> n_bin;

   // calculate inverse bin size for initial photon energy (log scale)
   inv_delta = float(ne_in-1)/log(energy_max/energy_min);

   // calculate distance between first and last cross section data bin
   distance = float(n_bin - 1);

   // test memory allocation
   if ( (x_data != NULL) || (x_matrix != NULL) )
   {
      xvmc_error("compton_XS_diff::read",
                 "memory for compton data already allocated",8);
   }

   // allocate memory for compton data
   if ( (x_data = new real[ne_in*n_bin]) == NULL )
   {
      xvmc_error("compton_XS_diff::read",
                 "cannot allocate memory for compton data",8);
   }

   // allocate memory for pointer matrix
   if ( (x_matrix = new real*[ne_in]) == NULL )
   {
      xvmc_error("compton_XS_diff::read",
                 "cannot allocate memory for pointer matrix",8);
   }

   // calculate pointer matrix
   for (register int i=0; i<ne_in; ++i) {
      x_matrix[i] = x_data + i*n_bin; }

   // read data from file
   for (register int i=0; i<ne_in; ++i) {
      for (register int j=0; j<n_bin; ++j) {
         file >> x_matrix[i][j]; } }

   // close file
   file.close();
}
