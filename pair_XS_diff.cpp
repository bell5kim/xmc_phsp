/*****************************************************************************
 * pair_XS_diff.cpp:                                                         *
 *    class member functions for:                                            *
 *       pair_XS_diff:   differential pair cross section                     *
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
using namespace std;

#include "pair_XS_diff.h"

/*****************************************************************************
 * member functions for class pair_XS_diff:                                  *
 *    differential pair cross section                                        *
 *****************************************************************************/

//##################################################################################
// NEW data container for compton_XS_diff
//##################################################################################

PairDiffData::PairDiffData(int nEin, int nBin, real energyMin, real energyMax) :
  ne_in(nEin),
  n_bin(nBin),
  energy_min(energyMin),
  energy_max(energyMax),
  x_data(NULL),
  x_matrix(NULL)
{
  x_data = new real[ne_in*n_bin];
  x_matrix = new real*[ne_in];
  f_data = new real[ne_in*n_bin];
  f_matrix = new real*[ne_in];

  // init pointer array for direct access
  int i;
  for (i=0; i<ne_in; ++i) {
    x_matrix[i] = x_data + i*n_bin;
    f_matrix[i] = f_data + i*n_bin;
  }
}

//##################################################################################

PairDiffData::~PairDiffData()
{
  delete[] x_data;
  delete[] x_matrix;
  delete[] f_data;
  delete[] f_matrix;
}

//##################################################################################

PairDiffData* PairDiffFactory(const std::string& file_name)
{
  real energy_min, energy_max;
  int ne_in, n_bin;

  // open data file
  ifstream file(file_name.c_str(), ios::in);
  if (!file) {
    xvmc_error("pair_XS_diff::read","cannot open file",8); }

  // read minimum and maximum energy
  file >> energy_min >> energy_max;

  // read the number of bins
  file >> ne_in >> n_bin;

  PairDiffData* pcxs = new PairDiffData(ne_in, n_bin, energy_min, energy_max);

  int i, j;

  // read data from file
  for (i=0; i<ne_in; ++i) {
    for (j=0; j<n_bin; ++j) {
      file >> pcxs->x_matrix[i][j];
    }
    for (j=0; j<n_bin; ++j) {
      file >> pcxs->f_matrix[i][j];
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

pair_XS_diff::pair_XS_diff(const PairDiffData& A) :
  ne_in(A.ne_in),
  n_bin(A.n_bin),
  energy_min(A.energy_min),
  energy_max(A.energy_max),
  x_data(NULL),
  f_data(NULL),
  x_matrix(NULL),
  f_matrix(NULL)
{
  // calculate bin size for initial photon energy (log scale)
  delta = log(energy_max/energy_min)/float(ne_in-1);
  // calculate distance between first and last cross section data bin
  distance = float(n_bin - 1);

  x_data = new real[ne_in*n_bin];
  x_matrix = new real*[ne_in];
  f_data = new real[ne_in*n_bin];
  f_matrix = new real*[ne_in];

  // init pointer array for direct access
  int i, j;
  for (int i=0; i<ne_in; ++i) {
    x_matrix[i] = x_data + i*n_bin;
    f_matrix[i] = f_data + i*n_bin;
  }

  for (i=0; i<ne_in; ++i) {
    for (j=0; j<n_bin; ++j) {
      x_matrix[i][j] = A.x_matrix[i][j];
      f_matrix[i][j] = A.f_matrix[i][j];
    }
  }
}

//##################################################################################

pair_XS_diff:: ~pair_XS_diff(void)
{
  delete [] x_matrix;
  delete [] x_data;
  delete [] f_matrix;
  delete [] f_data;
}

//##################################################################################
// read data from file
void pair_XS_diff::read(char *file_name)
{
   // open data file
   ifstream file(file_name, ios::in);
   if (!file) {
      xvmc_error("pair_XS_diff::read","cannot open file",8); }

   // read minimum and maximum energy
   file >> energy_min >> energy_max;

   // read the number of bins
   file >> ne_in >> n_bin;

   // calculate bin size for initial photon energy (log scale)
   delta = log(energy_max/energy_min)/float(ne_in-1);

   // calculate distance between first and last cross section data bin
   distance = float(n_bin - 1);

   // test memory allocation
   if ( (x_data != NULL) || (x_matrix != NULL) ||
        (f_data != NULL) || (f_matrix != NULL) )
   {
      xvmc_error("pair_XS_diff::read",
                 "memory for pair data already allocated",8);
   }

   // allocate memory for pair data
   if ( (x_data = new real[ne_in*n_bin]) == NULL )
   {
      xvmc_error("pair_XS_diff::read",
                 "cannot allocate memory for pair x_data",8);
   }

   if ( (f_data = new real[ne_in*n_bin]) == NULL )
   {
      xvmc_error("pair_XS_diff::read",
                 "cannot allocate memory for pair f_data",8);
   }

   // allocate memory for pointer matrix
   if ( (x_matrix = new real*[ne_in]) == NULL )
   {
      xvmc_error("pair_XS_diff::read",
                 "cannot allocate memory for pointer x_matrix",8);
   }

   if ( (f_matrix = new real*[ne_in]) == NULL )
   {
      xvmc_error("pair_XS_diff::read",
                 "cannot allocate memory for pointer f_matrix",8);
   }

   // calculate pointer matrix
   for (register int i=0; i<ne_in; ++i)
   {
      x_matrix[i] = x_data + i*n_bin;
      f_matrix[i] = f_data + i*n_bin;
   }

   // read data from file
   for (register int i=0; i<ne_in; ++i)
   {
      for (register int j=0; j<n_bin; ++j) {
         file >> x_matrix[i][j]; }
      for (register int j=0; j<n_bin; ++j) {
         file >> f_matrix[i][j]; }
   }

   // close file
   file.close();
}

