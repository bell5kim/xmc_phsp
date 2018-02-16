#ifndef _ARRAY_2D_H_
#define _ARRAY_2D_H_

/*****************************************************************************
 * array_2d.h:                                                               *
 *    class declaration:                                                     *
 *       array_2d:       2D matrix for various types (template class)        *
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

#include "definitions.h"

// 2D matrix for portal images etc.
template <class type>
class array_2d
{
   public:
// Modified by JOKim 25NOV2009 Moved from private --------------
      // matrix dimensions
      unsigned int       dim_x,dim_y;

      // size of the 2d data array (dim_x * dim_y)
      unsigned int       size;

      // pointer to 2d data
      type              *data;
// End of Modified      ----------------------------------------
      // pointer to 2d matrix
      type            **matrix;

      // default constructor returns empty array
      array_2d(void) { init(); }

      // 2nd constructor returns non-empty array filled with value
      array_2d(unsigned int new_dim_x, unsigned int new_dim_y, type value) {
         init(); fill(new_dim_x, new_dim_y, value); }

      // copy constructor
      array_2d(const array_2d<type> &old) {
         init(); copy(old); }

      // destructor
      ~array_2d(void) { clear(); }

      // assignment operators
      array_2d& operator = (const array_2d<type> &);
      array_2d& operator = (type value);

   private:
   
      // initialize empty matrix
      void init(void);

      // allocate memory for matrix
      void get_memory(unsigned int, unsigned int);

      // initialize and fill array with value
      void fill(unsigned int, unsigned int, type);

      // clear array and free memory
      void clear(void);

      // copy data of another array
      void copy(const array_2d<type> &);
};

#endif /* _ARRAY_2D_H_ */
