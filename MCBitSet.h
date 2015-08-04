#ifndef _MCBITSET_H_
#define _MCBITSET_H_

/*****************************************************************************
 * MCBitSet.h:                                                               *
 *    class declarations and inline member functions for:                    *
 *       MCBitSet:       manipulate bit vectors                              *
 *                                                                           *
 * Copyright (C) 2001    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 08.11.2001      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include <new>
using namespace std;

#include "definitions.h"
#include "global.h"
//void xvmc_error(const char *, const char *, const int);

// ****************************************************************
// class MCBitSet: manipulate bit vectors
// ****************************************************************

class MCBitSet
{
   public:
      // create bit set by a number of bits
      MCBitSet(unsigned num) { init(num); }

      // copy constructor
      MCBitSet(const MCBitSet &orig) {
         init(orig.number); copy(orig); }

      // delete bit set
      inline ~MCBitSet(void);

      // clear the whole bit set
      inline void clear(void)
         { for (register unsigned i=0; i<num_words; ++i) data[i] = 0x00; }

      // get number of bits
      unsigned get_number(void) const { return(number); }

      // clear one bit (set to 0)
      inline void clear(unsigned);

      // set one bit (set to 1)
      inline void set(unsigned);

      // invert one bit
      inline void invert(unsigned);

      // test one bit
      inline char test(unsigned);

      // bitwise AND
      inline MCBitSet& bit_and(const MCBitSet &, const MCBitSet &);

      // copy bit set
      MCBitSet& operator = (const MCBitSet &orig)
         {  copy(orig); return(*this); }

      // bitwise AND
      MCBitSet& operator &= (const MCBitSet &orig) {
         for (register unsigned i=0; i<num_words; ++i) {
            data[i] &= orig.data[i];
         }
         return(*this);
      }

   private:
      // number of bits
      unsigned       number;

      // data array to store the bits
      char *data;

      // number of data elements
      unsigned       num_words;

      // initialize bit set by a number of bits
      inline void init(unsigned);

      // copy bit set
      inline void copy(const MCBitSet &);

      // operators with access to private members
      friend MCBitSet operator &  (const MCBitSet &, const MCBitSet &);
      friend int      operator == (const MCBitSet &, const MCBitSet &);
};

// delete bit set
inline MCBitSet::~MCBitSet(void)
{
   if (data != NULL)
   {
      delete [] data;
      data = NULL;
   }
}

// initialize bit set by a number of bits
inline void MCBitSet::init(unsigned num)
{
#ifdef DEBUG
   // check number of bits
   if (num < 1) xvmc_error("MCBitSet::init",
                           "number of bits must be larger than 0",8);
#endif

   // set number of bits
   number = num;

#ifdef DEBUG
   // calculate word size in bits (must be 8)
   unsigned word_size = 8*sizeof(char);
   if (word_size != 8) xvmc_error("MCBitSet::init",
                                  "type char is not 8 bits long",8);
#endif

   // calculate number of data elements
   num_words = (number-1)/8 + 1;

   // allocate memory for the new bit data array
#ifdef DEBUG
   data = NULL;
   if ( (data = new (nothrow) char[num_words]) == NULL )
   {
      xvmc_error("MCBitSet::init",
                 "cannot allocate memory for bit data",8);
   }
#else
   data = new char[num_words];
#endif

   // set all bits to zero
   clear();
}

// copy bit set
inline void MCBitSet::copy(const MCBitSet &orig)
{
#ifdef DEBUG
   // check number of bits
   if (number != orig.get_number())
      xvmc_error("MCBitSet::copy","the number of bits are different",8);
#endif

   // copy data
   for (register unsigned i=0; i<num_words; ++i) data[i] = orig.data[i];
}

// clear one bit (set to 0)
inline void MCBitSet::clear(unsigned n)
{
#ifdef DEBUG
   if (n >= number) xvmc_error("MCBitSet::clear","cannot clear this bit",8);
#endif
   data[n/8] &= ~(0x80 >> (n%8));
}

// set one bit (set to 1)
inline void MCBitSet::set(unsigned n)
{
#ifdef DEBUG
   if (n >= number) xvmc_error("MCBitSet::set","cannot set this bit",8);
#endif
   data[n/8] |= (0x80 >> (n%8));
}

// invert one bit
inline void MCBitSet::invert(unsigned n)
{
#ifdef DEBUG
   if (n >= number) xvmc_error("MCBitSet::invert","cannot invert this bit",8);
#endif
   data[n/8] ^= (0x80 >> (n%8));
}

// test one bit
inline char MCBitSet::test(unsigned n)
{
#ifdef DEBUG
   if (n >= number) xvmc_error("MCBitSet::test","cannot test this bit",8);
#endif
   return(data[n/8] & (0x80 >> (n%8)));
}

// bitwise AND
inline MCBitSet& MCBitSet::bit_and(const MCBitSet &s1, const MCBitSet &s2)
{
#ifdef DEBUG
   // check number of bits
   if ( (number != s1.get_number()) || (number != s2.get_number()) )
      xvmc_error("MCBitSet::bit_and","the number of bits are different",8);
#endif

   // create and copy data
   for (register unsigned i=0; i<num_words; ++i)
      data[i] = s1.data[i] & s2.data[i];

   return(*this);
}

// bitwise AND
inline MCBitSet operator & (const MCBitSet &s1, const MCBitSet &s2)
{
   unsigned num = s1.get_number();

#ifdef DEBUG
   // check number of bits
   if (num != s2.get_number())
      xvmc_error("MCBitSet::operator &","the number of bits are different",8);
#endif

   // create new bit set for the result
   MCBitSet result(num);

   // create and copy data
   for (register unsigned i=0; i<result.num_words; ++i)
      result.data[i] = s1.data[i] & s2.data[i];

   return(result);
}

// equality operator
inline int operator == (const MCBitSet &s1, const MCBitSet &s2)
{
#ifdef DEBUG
   // check number of bits
   if (s1.number != s2.number)
      xvmc_error("MCBitSet::operator ==","the number of bits are different",8);
#endif

   // compare
   int result = 1;
   for (register unsigned i=0; i<s1.num_words; ++i)
   {
      result = ( result && (s1.data[i] == s2.data[i]) );
   }

   return(result);
}

#endif /* _MCBITSET_H_ */
