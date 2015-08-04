/*****************************************************************************
 * MC_object.cpp:                                                            *
 *    class member functions for:                                            *
 *       MC_object:      geometrical object defined by regions and planes    *
 *                                                                           *
 * Copyright (C) 2001    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 31.08.2001      *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************

#include "MC_object.h"

// ****************************************
// member functions of class MC_object
// ****************************************

// allocate new object
MC_object::MC_object(unsigned n_planes, unsigned n_regions)
{
   // initialize number of planes and regions
   num_planes  = n_planes;
   num_regions = n_regions;

   // initialize pointers
   starting_plane = NULL;
   final_plane    = NULL;
   separator      = NULL;
   piece          = NULL;
   bit_mask       = NULL;
   bit_pattern    = NULL;
 
   if (num_planes > 0)
   {
      // create array of pointers to the separator planes
      if ( (separator = new (nothrow) MC_plane*[num_planes]) == NULL )
      {
         xvmc_error("MC_object::MC_object",
                    "cannot allocate memory for separator pointer array",8);
      }
      for (register unsigned ip=0; ip<num_planes; ++ip) separator[ip] = NULL;
   }

   if (num_regions > 0)
   {
      // create array of pointers to the pieces of the object
      if ( (piece = new (nothrow) MC_region*[num_regions]) == NULL )
      {
         xvmc_error("MC_object::MC_object",
                    "cannot allocate memory for piece pointer array",8);
      }

      // create array of bit mask pointers
      if ( (bit_mask = new (nothrow) MCBitSet*[num_regions]) == NULL )
      {
         xvmc_error("MC_object::MC_object",
                    "cannot allocate memory for bit mask pointer array",8);
      }

      // create array of bit pattern pointers
      if ( (bit_pattern = new (nothrow) MCBitSet*[num_regions]) == NULL )
      {
         xvmc_error("MC_object::MC_object",
                    "cannot allocate memory for bit pattern pointer array",8);
      }

      // create bit masks and bit patterns, initialize
      for (register unsigned i_region = 0; i_region < num_regions; ++i_region)
      {
         bit_mask[i_region]    = NULL;
         if ( (bit_mask[i_region] = new (nothrow)
               MCBitSet(num_planes)) == NULL )
         {
            xvmc_error("MC_object::MC_object",
                       "cannot allocate memory for bit masks",8);
         }
         bit_mask[i_region]->clear();

         bit_pattern[i_region] = NULL;
         if ( (bit_pattern[i_region] = new (nothrow)
               MCBitSet(num_planes)) == NULL )
         {
            xvmc_error("MC_object::MC_object",
                       "cannot allocate memory for bit patterns",8);
         }
         bit_pattern[i_region]->clear();

         piece[i_region] = NULL;
      }
   }
}

// delete object
MC_object::~MC_object(void)
{
   // delete separator planes
   if (num_planes > 0)
   {
      for (register unsigned ip=0; ip<num_planes; ++ip)
      {
         if (separator[ip] != NULL)
            { delete separator[ip]; separator[ip] = NULL; }
      }
      delete [] separator;
   }

   if (num_regions > 0)
   {
      for (register unsigned ir=0; ir<num_regions; ++ir)
      {
         // delete piece regions
         if (piece[ir] != NULL)
         {
            // all planes are already deleted, therefore we assign
            // "NULL" to the region surface plane pointers before
            // deleting the piece regions
            for (register unsigned iq=0; iq<piece[ir]->get_num_planes(); ++iq)
               {  piece[ir]->surface[iq] = NULL; }
            delete piece[ir]; piece[ir] = NULL;
         }

         // delete bit masks
         if (bit_mask[ir] != NULL)
         {
            delete bit_mask[ir]; bit_mask[ir] = NULL;
         }

         // delete bit patterns
         if (bit_pattern[ir] != NULL)
         {
            delete bit_pattern[ir]; bit_pattern[ir] = NULL;
         }
      }
      delete [] piece;

      // delete array of bit mask pointers
      delete [] bit_mask;

      // delete array of bit pattern pointers
      delete [] bit_pattern;
   }
 
   // reset pointers
   starting_plane = NULL;
   final_plane    = NULL;
   separator      = NULL;
   piece          = NULL;
   bit_mask       = NULL;
   bit_pattern    = NULL;

   // set everything to zero
   num_planes  = 0;
   num_regions = 0;
}

// set plane indices for all planes of the object,
// set bit masks of all regions by checking the surface plane indices and
// set bit patterns of all regions by checking the relationships of the
// region reference points to the corresponding surface planes
void MC_object::set_bits(void)
{
   // first of all, set plane indices
   for (register unsigned i_plane = 0; i_plane < num_planes; ++i_plane)
   {
      separator[i_plane]->index = i_plane;
   }

   for (register unsigned i_region = 0; i_region < num_regions; ++i_region)
   {
      // determine and set the bit mask for each region
      for (register unsigned i = 0; i < piece[i_region]->get_num_planes(); ++i)
      {
         int plane_index = piece[i_region]->surface[i]->index;
         if (plane_index < 0)
         {
            xvmc_error("MC_object::set_bits","plane index not set",8);
         }
         else
         {
            bit_mask[i_region]->set(plane_index);
         }
      }

      // check the reference points relationships to the surface planes
      // of this corresponding region
      for (register unsigned k = 0; k < num_planes; ++k)
      {
         if (bit_mask[i_region]->test(k))
         {
            // get the reference point of this region
            real_3 p_ref = piece[i_region]->get_p_ref();

            // this plane is a surface plane of the region,
            // check the reference point to plane relationship
            if      (separator[k]->relationship(p_ref) < 0)
               bit_pattern[i_region]->set(k);
            else if (separator[k]->relationship(p_ref) > 0)
               bit_pattern[i_region]->clear(k);
            else
            {
               xvmc_error("MC_object::set_bits",
                          "the reference point is located at the plane",8);
            }
         }
         else
         {
            // this plane has nothing to do with the region, clear bit
            bit_pattern[i_region]->clear(k);
         }
      }
   }
}
