/*****************************************************************************
 * array_3d.cpp:                                                             *
 *    class member functions for:                                            *
 *       array_3d: 3D matrix for various types (template class)              *
 *                 (for density distribution, dose distribution etc.         *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 00/02/22        *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <string.h>
#include "definitions.h"
#include "global.h"
#include "array_3d.h"

#ifndef OSF1
// template instantiation of different "array_3d" types
template class array_3d<float>;
template class array_3d<double>;
#endif

// initialize empty matrix
template <class type>
void array_3d<type>::init(void)
{
   dimension.x = 0;
   dimension.y = 0;
   dimension.z = 0;
   size = 0;
   data   = NULL;
   matrix = NULL;
}

// allocate memory for matrix
template <class type>
void array_3d<type>::get_memory(const int_3 &new_dim)
{
   // array must be empty
   if (size != 0)
      xvmc_error("array_3d::get_memory", "array not empty", 8);

   // test new dimensions
   if ((new_dim.x <= 0) || (new_dim.y <= 0) || (new_dim.z <= 0))
      xvmc_error("array_3d::get_memory", "wrong array dimension", 8);

   // update matrix dimension
   dimension.x = new_dim.x;
   dimension.y = new_dim.y;
   dimension.z = new_dim.z;

   // calculate new size
   size = dimension.x*dimension.y*dimension.z;

   // auxiliary variables
   int dim_yz = dimension.y*dimension.z;
   int dim_z  = dimension.z;

   // allocate memory for new array data
   if ( (data = new type[size]) == NULL )
   {
      xvmc_error("array_3d::get_memory",
                 "cannot allocate memory for array data",8);
   }

   // allocate memory for pointer matrix (x dimension)
   if ( (matrix = new type**[dimension.x]) == NULL )
   {
      xvmc_error("array_3d::get_memory",
                 "cannot allocate memory for matrix (x dimension)",8);
   }

   // calculate pointer matrix
   for (register int i=0; i<dimension.x; ++i)
   {
      // allocate memory for pointer matrix (y dimension)
      if ( (matrix[i] = new type*[dimension.y]) == NULL )
      {
         xvmc_error("array_3d::get_memory",
                    "cannot allocate memory for matrix (y dimension)",8);
      }

      for (register int j=0; j<dimension.y; ++j)
      {
         // assign pointers
         matrix[i][j] = data + i*dim_yz + j*dim_z; 
      }
   }
}

// initialize and fill array with a single value
template <class type>
void array_3d<type>::fill(const int_3 &new_dim, type value)
{
   // free memory
   if (size != 0) clear();

   // get new memory
   get_memory(new_dim);

   // fill array with value
   for (register unsigned int i=0; i<size; ++i)
      { *(data + i) = value; }
}

// clear array and free memory
template <class type>
void array_3d<type>::clear(void)
{
   // array is empty
   if (size == 0) return;

   // free old data memory
   if (data != NULL) { delete [] data; data = NULL; }

   // free old matrix memory
   if (dimension.x > 0)
   {
      for (register int i=0; i<dimension.x; ++i)
      {
         if (matrix[i] != NULL) { delete [] matrix[i]; matrix[i] = NULL; }
      }
   }
   if (matrix != NULL) { delete [] matrix; matrix = NULL; }

   dimension.x = 0;
   dimension.y = 0;
   dimension.z = 0;
   size = 0;
}

// copy data of another array
template <class type>
void array_3d<type>::copy(const array_3d<type> &old)
{
   // free memory
   if (size != 0) clear();

   if (old.size != 0)
   {
      // get new memory
      get_memory(old.dimension);

      // fill array
      for (register unsigned int i=0; i<size; ++i)
         { *(data + i) = *(old.data + i); }
   }
}

// assign data of another array
template <class type>
array_3d<type>& array_3d<type>::operator = (const array_3d<type> &old)
{
   copy(old);
   return(*this);
}

// fill matrix with value
template <class type>
array_3d<type>& array_3d<type>::operator = (type value)
{
  if(value == ZERO) {
    memset((void *) data, 0, size*sizeof(type));
  }
  else {
    for (register unsigned int i=0; i<size; ++i) { 
      *(data + i) = value; 
    }
  }
  return(*this);
}
