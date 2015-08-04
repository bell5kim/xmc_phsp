#ifndef _ARRAY_3D_H_
#define _ARRAY_3D_H_

/*****************************************************************************
 * array_3d.h:                                                               *
 *    class declaration:                                                     *
 *       array_3d:            3D matrix for various types (template class)   *
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

#include <stdlib.h>
#include "definitions.h"

// 3D matrix for density distribution, dose distribution etc.
template <class type>
class array_3d
{
   public:
      int_3              dimension;      // matrix dimensions
      unsigned int       size;           // dim.x * dim.y * dim.z
      type              *data;           // pointer to 3d data
      type            ***matrix;         // pointer to 3d matrix
      // default constructor returns empty array
      array_3d(void) { init(); }
      // 2nd constructor returns non-empty array filled with value
      array_3d(int_3 &new_dim, type value) {
         init(); fill(new_dim, value); }
      // copy constructor
      array_3d(const array_3d<type> &old) {
         init(); copy(old); }
      // destructor
      ~array_3d(void) { clear(); }
      array_3d<type>& operator = (const array_3d<type> &);
      array_3d<type>& operator = (type value);
      
   private:
      // initialize empty matrix
      void init(void);
      // allocate memory for matrix
      void get_memory(const int_3 &);
      // initialize and fill array with value
      void fill(const int_3 &, type);
      // clear array and free memory
      void clear(void);
      // copy data of another array
      void copy(const array_3d<type> &);
};

#endif /* _ARRAY_3D_H_ */
