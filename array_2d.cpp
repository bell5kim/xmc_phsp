/*****************************************************************************
 * array_2d.cpp:                                                             *
 *    class member functions for:                                            *
 *       array_2d: 2D matrix for various types (template class)              *
 *                 (for portal images etc.)                                  *
 *                                                                           *
 * Copyright (C) 2003    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding (adapted from array_3d)              MF 03/03/13        *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include "global.h"
#include "array_2d.h"

#ifndef OSF1
// template instantiation of different "array_2d" types
template class array_2d<float>;
template class array_2d<double>;
#endif

// initialize empty matrix
template <class type>
void array_2d<type>::init(void)
{
   dim_x  = 0;
   dim_y  = 0;
   size   = 0;
   data   = NULL;
   matrix = NULL;
}

// allocate memory for matrix
template <class type>
void array_2d<type>::get_memory(unsigned int new_dim_x, unsigned int new_dim_y)
{
   // array must be empty
   if (size != 0)
      xvmc_error("array_2d::get_memory","array not empty", 8);

   // test new dimensions
   if ((new_dim_x == 0) || (new_dim_y == 0))
      xvmc_error("array_2d::get_memory","cannot create array of zero size", 8);

   // update matrix dimension
   dim_x = new_dim_x;
   dim_y = new_dim_y;

   // calculate new size
   size = dim_x*dim_y;

   // allocate memory for new array data
   data = NULL;
   if ( (data = new type[size]) == NULL )
   {
      xvmc_error("array_2d::get_memory",
                 "cannot allocate memory for array data",8);
   }

   // allocate memory for pointer matrix
   matrix = NULL;
   if ( (matrix = new type*[dim_x]) == NULL )
   {
      xvmc_error("array_2d::get_memory",
                 "cannot allocate memory for pointer matrix",8);
   }

   // calculate pointer matrix
   for (register unsigned int i=0; i<dim_x; ++i)
   {
      // assign pointers
      matrix[i] = data + i*dim_y; 
   }
}

// initialize and fill array with a single value
template <class type>
void array_2d<type>::fill(unsigned int new_dim_x, unsigned int new_dim_y,
                          type value)
{
   // free memory
   if (size != 0) clear();

   // get new memory
   get_memory(new_dim_x,new_dim_y);

   // fill array with value
   for (register unsigned int i=0; i<size; ++i)
      { *(data + i) = value; }
}

// clear array and free memory
template <class type>
void array_2d<type>::clear(void)
{
   // array is empty
   if (size == 0) return;

   // free old data memory
   if (data != NULL) { delete [] data; data = NULL; }

   // free old matrix memory
   if (matrix != NULL) { delete [] matrix; matrix = NULL; }

   dim_x = 0;
   dim_y = 0;
   size  = 0;
}

// copy data of another array
template <class type>
void array_2d<type>::copy(const array_2d<type> &old)
{
   // free memory
   if (size != 0) clear();

   if (old.size != 0)
   {
      // get new memory
      get_memory(old.dim_x,old.dim_y);

      // fill array
      for (register unsigned int i=0; i<size; ++i)
         { *(data + i) = *(old.data + i); }
   }
}

// assign data of another array
template <class type>
array_2d<type>& array_2d<type>::operator = (const array_2d<type> &old)
{
   copy(old);
   return(*this);
}

// fill matrix with value
template <class type>
array_2d<type>& array_2d<type>::operator = (type value)
{
  for (register unsigned int i=0; i<size; ++i)
  { 
     *(data + i) = value; 
  }
  return(*this);
}
