#ifndef _MC_OBJECT_H_
#define _MC_OBJECT_H_

/*****************************************************************************
 * MC_object.h:                                                              *
 *    class declarations and inline member functions for:                    *
 *       MC_object:      geometrical object defined by regions and planes    *
 *                                                                           *
 * Copyright (C) 2001    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 29.08.2001      *
 *                                                                           *
 *****************************************************************************/

#include <new>
using namespace std;

#define SECONDARY_MODIFIER_ELECTRONS

//#define CPU_WATCH
#ifdef CPU_WATCH
#include "cpu_watch.h"
#endif

#include "definitions.h"
#include "MCBitSet.h"
#include "xvmc_util.h"
#include "ranmar.h"
#include "MC_plane.h"
#include "MC_region.h"

// *******************************************************************
// class MC_object: geometrical object defined by regions and planes
// *******************************************************************

class MC_object
{
   public:
      MC_object(unsigned, unsigned);  // allocate new object
      virtual ~MC_object(void);       // delete object

      // get number of planes defining the object
      unsigned  get_num_planes(void) { return(num_planes); }

      // get bit mask of one region
      MCBitSet  get_bit_mask(const unsigned i_region)
                   { return(*bit_mask[i_region]); }

      // get bit pattern of one region
      MCBitSet  get_bit_pattern(const unsigned i_region)
                   { return(*bit_pattern[i_region]); }

      // get number of regions defining the object
      unsigned  get_num_regions(void) { return(num_regions); }

      // get particle starting plane of the object
      MC_plane *get_starting_plane(void) { return(starting_plane); }

      // get final particle transport plane of the object
      MC_plane *get_final_plane(void) { return(final_plane); }

      // estimate the region index of a point, default: -1 == no estimation
      virtual int estimate_region_index(const real_3 &p0) { return(-1); }

      // calculate particle bit pattern by determining the plane relationships,
      // only the bit associated with the specified plane will be set explicitly
      inline MCBitSet calc_bit_pattern(const particle_parameters &,
                                       const MC_plane *, const bool &);

      // determine the region index of the particle by comparing the
      // the masked particle bit pattern with the region bit patterns,
      // return -1 if the particle is outside the object,
      // the integer input variable can be used for an index estimation
      inline int get_region_index(const MCBitSet &, const int &);

      // transport photon through the object, the plane pointer points to
      // the starting plane of the photon (it must not be the NULL pointer),
      // if the photon survives (or one daughter particle) -> return true,
      // if no particle survives  -> return false,
      // the particle_parameters struct returns at most one particle
      inline bool primary_photon(particle_parameters &,
                                 MC_plane *, ranmar &);

      // transport electron through the object, the plane pointer points to
      // the starting plane of the electron (it must not be the NULL pointer),
      // if the electron survives (and/or daughter particles) -> return true,
      // if no particle survives  -> return false,
      // this is just a Continuous Slowing Down Approximation (CSDA),
      // the random number generator is reserved for later use
      inline bool primary_electron(particle_parameters &,
                                   MC_plane *, ranmar &);

      // transport secondary electron through the object,
      // the bit pattern and the region index of the electron must be known,
      // if the electron survives (and/or daughter particles) -> return true,
      // if no particle survives  -> return false,
      // this is just a Continuous Slowing Down Approximation (CSDA),
      // the random number generator is reserved for later use
      inline bool secondary_electron(particle_parameters &,
                                     MCBitSet &, int &, ranmar &);

      // array of pointers to plane defining the object
      MC_plane    **separator;

      // array of pointers to region
      MC_region   **piece;

   protected:
      // number of planes
      unsigned      num_planes;

      // number of regions, including the outer region with index 0
      unsigned      num_regions;

      // pointer to the starting plane
      // (generally the particle transport begins at this plane)
      MC_plane     *starting_plane;

      // pointer to the final plane
      // (the transport will be stopped if the particle crosses this plane)
      MC_plane     *final_plane;

      // set plane indices for all planes of the object,
      // set bit masks of all regions by checking the surface plane indices and
      // set bit patterns of all regions by checking the relationships of the
      // region reference points to the corresponding surface planes
      void set_bits(void);

   private:
      // for each region of this object we define a bit mask,
      // set bit if the corresponding plane is a surface plane of the region,
      // clear bit if the corresponding plane has nothing to do with the plane
      MCBitSet **bit_mask;

      // for each region of this object we define a bit pattern
      // providing the region to plane relationships,
      // set bit to 1 (or true) if the corresponding plane is a surface plane
      // of the region (see bit_mask) and the region is behind this plane,
      // clear bit (set to 0 or to false) if the corresponding plane is a
      // surface plane of the region (see bit_mask) and the region is before
      // this plane, clear also bit if the corresponding plane has nothing
      // to do with the region (see bit_mask)
      MCBitSet **bit_pattern;
};

// *******************************************
// inline member functions of class MC_object
// *******************************************

// calculate particle bit pattern by determining the plane relationships,
// only the bit associated with the specified plane will be set explicitly
inline MCBitSet MC_object::calc_bit_pattern(const particle_parameters &p,
                                            const MC_plane *p_plane,
                                            const bool &p_plane_bit)
{
#ifdef CPU_WATCH
   extern cpu_watch t_cpu_pbp; t_cpu_pbp.start();
#endif

   MCBitSet p_bit_pattern(num_planes);

   for (register unsigned i_plane=0; i_plane<num_planes; ++i_plane)
   {
      if (separator[i_plane] == p_plane)
      {
         // set bit explicitly for plane p_plane
         if (p_plane_bit) p_bit_pattern.set(i_plane);
         else             p_bit_pattern.clear(i_plane);
      }
      else
      {
         // calculate bits for all planes with the exception of plane p_plane
         if (separator[i_plane]->relationship(p.pos) < 0)
         {
            // the particle is behind the plane, set bit
            p_bit_pattern.set(i_plane);
         }
         else
         {
            // the particle is either before or at the plane
            if (separator[i_plane]->relationship(p.pos) > 0)
            {
               // the particle is before the plane, clear bit
               p_bit_pattern.clear(i_plane);
            }
            else
            {
               // the particle is at the plane, therefore we take the particle
               // direction into account, compare with the plane normal vector
               real_3 normal = separator[i_plane]->get_normal();
            
               // scalar product of particle direction and plane normal vectors
               real product =
                  p.dir.x*normal.x + p.dir.y*normal.y + p.dir.z*normal.z;

               if (product > ZERO)
               {
                  // the particle points to the region before the plane,
                  // clear bit
                  p_bit_pattern.clear(i_plane);
               }
               else
               {
                  // the particle points to the region behind the plane,
                  // set bit,
                  // or it will stay at the plane, set also bit because
                  // in this case we define the particle to be behind the plane
                  p_bit_pattern.set(i_plane);
               }
            }
         }
      }
   }


#ifdef CPU_WATCH
   t_cpu_pbp.stop();
#endif

   // return new particle bit pattern
   return(p_bit_pattern);
}

// determine the region index of the particle by comparing the
// the masked particle bit pattern with the region bit patterns,
// return -1 if the particle is outside the object,
// the integer input variable can be used for an index estimation
inline int MC_object::get_region_index(const MCBitSet &p_bit_pattern,
                                       const int      &estimate_i_region)
{
#ifdef CPU_WATCH
   extern cpu_watch t_cpu_ind; t_cpu_ind.start();
#endif

   MCBitSet test_bit_set(p_bit_pattern);

   // test whether we can use the estimated index
   if ( (estimate_i_region > 0) && (estimate_i_region < int(num_regions)) )
   {
      // use estimated index

      // the new region index should be close to the estimated region index,
      // therefore we start the index search from this estimated index
      for (register int i1=estimate_i_region-1, i2=estimate_i_region;
                        ( (0 <= i1) || (i2<int(num_regions)) ); --i1, ++i2)
      {
         if (i1 >= 0)
         {
            if ( ( test_bit_set.bit_and(p_bit_pattern,*bit_mask[i1]) )
                                                        == *bit_pattern[i1] )
            {
#ifdef CPU_WATCH
               t_cpu_ind.stop();
#endif
               return(i1);
            }
         }

         if (i2 < int(num_regions))
         {
            if ( ( test_bit_set.bit_and(p_bit_pattern,*bit_mask[i2]) )
                                                     == *bit_pattern[i2] )
            {
#ifdef CPU_WATCH
               t_cpu_ind.stop();
#endif
               return(i2);
            }
         }
      }
   }
   else
   {
      // don't use estimated index
      for (register unsigned i0=0; i0<num_regions; ++i0)
      {
         if ( ( test_bit_set.bit_and(p_bit_pattern,*bit_mask[i0]) )
                                                     == *bit_pattern[i0] )
         {
#ifdef CPU_WATCH
            t_cpu_ind.stop();
#endif
            return(i0);
         }
      }
   }
 
#ifdef CPU_WATCH
   t_cpu_ind.stop();
#endif

   return(-1);
}

// transport photon through the object, the plane pointer points to
// the starting plane of the photon (it must not be the NULL pointer),
// if the photon survives (or one daughter particle) -> return true,
// if no particle survives  -> return false,
// the particle_parameters struct returns at most one particle
inline bool MC_object::primary_photon(particle_parameters &p,
                                      MC_plane *p_plane, ranmar &rndm)
{
   // at the beginning of this routine the primary photon is still alive
   bool primary_photon_alive = true;

#ifdef SECONDARY_MODIFIER_ELECTRONS
   // we create a linked list of secondary particles,
   // the pointers point to the first, last and new entries in this list
   particle_parameters *p_first = NULL;
   particle_parameters *p_last  = NULL;
   particle_parameters *p_new   = NULL;

   // number of secondary particles in linked list
   int num_secondaries = 0;

   // Compton electron
   particle_parameters e_comp;
   e_comp.type     = ELECTRON;
   e_comp.weight   = p.weight;
   e_comp.i.x      = 0;
   e_comp.i.y      = 0;
   e_comp.i.z      = 0;
   e_comp.next     = NULL;

   // pair electron
   particle_parameters e_pair;
   e_pair.type     = ELECTRON;
   e_pair.weight   = p.weight;
   e_pair.i.x      = 0;
   e_pair.i.y      = 0;
   e_pair.i.z      = 0;
   e_pair.next     = NULL;

   // pair positron
   particle_parameters p_pair;
   p_pair.type     = POSITRON;
   p_pair.weight   = p.weight;
   p_pair.i.x      = 0;
   p_pair.i.y      = 0;
   p_pair.i.z      = 0;
   p_pair.next     = NULL;

   // photo electron
   particle_parameters e_phot;
   e_phot.type     = ELECTRON;
   e_phot.weight   = p.weight;
   e_phot.i.x      = 0;
   e_phot.i.y      = 0;
   e_phot.i.z      = 0;
   e_phot.next     = NULL;
#endif

   // there must be an error if p_plane points to NULL
   if (p_plane == NULL)
   {
      xvmc_error("MC_object::primary_photon",
                 "the starting plane is undefined",8);
   }

   // determine the relationship of the photon to p_plane by taking only the
   // photon direction into account, compare with the plane normal vector
   real_3 normal = p_plane->get_normal();
            
   // scalar product of photon direction and plane normal vector
   real product = p.dir.x*normal.x + p.dir.y*normal.y + p.dir.z*normal.z;

   bool p_plane_bit;
   if (product > ZERO)
   {
      // the photon points to the region before the plane, clear bit
      p_plane_bit = false;
   }
   else
   {
      // the photon points to the region behind the plane, set bit
      // or it will stay at the plane, set also bit because
      // in this case we define the photon to be behind the plane
      p_plane_bit = true;
   }

   // determine the photon bit pattern for all planes
   // including the starting plane (set by bool p_plane_bit)
   MCBitSet p_bit_pattern = calc_bit_pattern(p,p_plane,p_plane_bit);

   // estimate and get region index of the photon
   int old_i_region = estimate_region_index(p.pos);
   int p_i_region   = get_region_index(p_bit_pattern,old_i_region);

   // the photon is outside of the object, return false
   if (p_i_region < 0)
   {
      p.energy = ZERO;
      p.weight = ZERO;
      return(false);
   }

   // number of mean free paths to the next interaction point
   real num_mfp = -log(ONE-rndm.number());

   // now trace photon through the object
   while (p_i_region >= 0)
   {
      // distance to the next plane and next plane
      real distance = piece[p_i_region]->distance(p,p_plane);

      if (p_plane != NULL)
      {
         // total Compton cross section
         real mu_comp = piece[p_i_region]->tot_comp->get(p.energy);

         // total pair cross section
         real mu_pair = piece[p_i_region]->tot_pair->get(p.energy);

         // total photo cross section
         real mu_phot = piece[p_i_region]->tot_phot->get(p.energy);

         // total attenuation coefficient
         real mu_tot = mu_comp + mu_pair + mu_phot;

         // Compton interaction probability
         real prob_comp = mu_comp/mu_tot;

         // pair interaction probability
         real prob_pair = mu_pair/mu_tot;

         // first interval:  Compton interaction
         // second interval: pair production
         // third interval:  photo-electric absorption
         prob_pair += prob_comp;

         // path length to the next interaction point
         real path_length = num_mfp/mu_tot;

         // compare photon path length and distance
         if (path_length < distance)
         {
            // move photon to the interaction site
            p.pos.x += path_length*p.dir.x;
            p.pos.y += path_length*p.dir.y;
            p.pos.z += path_length*p.dir.z;

            // the photon is within the region (not at a surface plane)
            p_plane = NULL;

            // determine interaction type
            real eta_itype = rndm.number();
            if (eta_itype <= prob_comp)
            {
               // Compton scattering
               real rnno = rndm.number();  // just a random number
               real energy_cx;             // Compton photon energy
               real energy_ce;             // Compton electron energy
               real cos_t_cx,sin_t_cx;     // Compton photon scattering angle
               real cos_t_ce,sin_t_ce;     // Compton electron scattering angle
               piece[p_i_region]->compton->interaction(p.energy, rnno, rndm,
                                              energy_cx, cos_t_cx, sin_t_cx,
                                              energy_ce, cos_t_ce, sin_t_ce);

               // azimuthal scattering angle
               real phi = TWO_PI * rndm.number();
               real sin_phi = sin(phi);
               real cos_phi = cos(phi);

#ifdef SECONDARY_MODIFIER_ELECTRONS
               if (energy_ce > e_cut)
               {
                  // simulate Compton electron
                  e_comp.energy = energy_ce;
                  e_comp.pos.x  = p.pos.x;
                  e_comp.pos.y  = p.pos.y;
                  e_comp.pos.z  = p.pos.z;
                  e_comp.dir.x  = p.dir.x;
                  e_comp.dir.y  = p.dir.y;
                  e_comp.dir.z  = p.dir.z;

                  // change electron direction
                  rotate(cos_t_ce, sin_t_ce, cos_phi, sin_phi, e_comp.dir);

                  // set electron region identifiers
                  MCBitSet ec_bit_pattern = p_bit_pattern;
                  int      ec_i_region    = p_i_region;

                  // transport electron
                  if ( secondary_electron(e_comp,
                                          ec_bit_pattern, ec_i_region,
                                          rndm) )
                  {
                     // the electron survived, create new list entry
                     ++num_secondaries;
                     p_new  = NULL;
                     if ( (p_new = new (nothrow) particle_parameters) == NULL )
                     {
                        xvmc_error("MC_object::primary_photon",
                           "cannot allocate memory for Compton electron",8);
                     }

                     // store parameters
                     p_new->type     = e_comp.type;
                     p_new->energy   = e_comp.energy;
                     p_new->weight   = e_comp.weight;
                     p_new->pos.x    = e_comp.pos.x;
                     p_new->pos.y    = e_comp.pos.y;
                     p_new->pos.z    = e_comp.pos.z;
                     p_new->dir.x    = e_comp.dir.x;
                     p_new->dir.y    = e_comp.dir.y;
                     p_new->dir.z    = e_comp.dir.z;
                     p_new->i.x      = e_comp.i.x;
                     p_new->i.y      = e_comp.i.y;
                     p_new->i.z      = e_comp.i.z;
                     p_new->next     = e_comp.next;

                     // assign to p_first if this is the first secondary,
                     // assign to the next list entry otherwise
                     if (p_first == NULL) p_first      = p_new;
                     else                 p_last->next = p_new;

                     // the electron is now the last particle in the list
                     p_last = p_new;
                  }
               }
#endif
               // stop transport if the Compton photon energy is below p_cut
               if (energy_cx <= p_cut)
               {
                  // stop transport
                  p_i_region = -1;

                  // and kill photon
                  primary_photon_alive = false;
               }
               else
               {
                  // rotate direction of the photon
                  rotate(cos_t_cx, sin_t_cx, cos_phi, sin_phi, p.dir);

                  // change energy
                  p.energy = energy_cx;

                  // sample new number of mean free paths
                  num_mfp = -log(ONE-rndm.number());
               }
            }
            else
            {
#ifdef SECONDARY_MODIFIER_ELECTRONS
               // pair production or photoelectric absorption,
               if (eta_itype < prob_pair)
               {
                  // pair production
                  real rnno = rndm.number();  // just a random number
                  real energy_pe;             // pair electron energy
                  real energy_pp;             // pair positron energy
                  real cos_t_pe,sin_t_pe;     // pair electron scattering angle
                  real cos_t_pp,sin_t_pp;     // pair positron scattering angle
                  piece[p_i_region]->pair->interaction(p.energy, rnno, rndm,
                                              energy_pe, cos_t_pe, sin_t_pe,
                                              energy_pp, cos_t_pp, sin_t_pp);

                  // azimuthal scattering angle
                  real phi = TWO_PI * rndm.number();
                  real sin_phi = sin(phi);
                  real cos_phi = cos(phi);

                  if (energy_pe > e_cut)
                  {
                     // simulate pair electron
                     e_pair.energy = energy_pe;
                     e_pair.pos.x  = p.pos.x;
                     e_pair.pos.y  = p.pos.y;
                     e_pair.pos.z  = p.pos.z;
                     e_pair.dir.x  = p.dir.x;
                     e_pair.dir.y  = p.dir.y;
                     e_pair.dir.z  = p.dir.z;

                     // change electron direction
                     rotate(cos_t_pe, sin_t_pe, cos_phi, sin_phi, e_pair.dir);

                     // set pair electron region identifiers
                     MCBitSet ep_bit_pattern = p_bit_pattern;
                     int      ep_i_region    = p_i_region;

                     // transport electron
                     if ( secondary_electron(e_pair,
                                             ep_bit_pattern, ep_i_region,
                                             rndm) )
                     {
                        // the electron survived, create new list entry
                        ++num_secondaries;
                        p_new  = NULL;
                        if ( (p_new = new (nothrow) particle_parameters)==NULL )
                        {
                           xvmc_error("MC_object::primary_photon",
                              "cannot allocate memory for pair electron",8);
                        }

                        // store parameters
                        p_new->type     = e_pair.type;
                        p_new->energy   = e_pair.energy;
                        p_new->weight   = e_pair.weight;
                        p_new->pos.x    = e_pair.pos.x;
                        p_new->pos.y    = e_pair.pos.y;
                        p_new->pos.z    = e_pair.pos.z;
                        p_new->dir.x    = e_pair.dir.x;
                        p_new->dir.y    = e_pair.dir.y;
                        p_new->dir.z    = e_pair.dir.z;
                        p_new->i.x      = e_pair.i.x;
                        p_new->i.y      = e_pair.i.y;
                        p_new->i.z      = e_pair.i.z;
                        p_new->next     = e_pair.next;

                        // assign to p_first if this is the first secondary,
                        // assign to the next list entry otherwise
                        if (p_first == NULL) p_first      = p_new;
                        else                 p_last->next = p_new;

                        // the electron is now the last particle in the list
                        p_last = p_new;
                     }
                  }

                  if (energy_pp > e_cut)
                  {
                     // simulate pair positron
                     p_pair.energy = energy_pp;
                     p_pair.pos.x  = p.pos.x;
                     p_pair.pos.y  = p.pos.y;
                     p_pair.pos.z  = p.pos.z;
                     p_pair.dir.x  = p.dir.x;
                     p_pair.dir.y  = p.dir.y;
                     p_pair.dir.z  = p.dir.z;

                     // change positron direction
                     sin_phi = -sin_phi;
                     cos_phi = -cos_phi;
                     rotate(cos_t_pp, sin_t_pp, cos_phi, sin_phi, p_pair.dir);

                     // set pair positron region identifiers
                     MCBitSet pp_bit_pattern = p_bit_pattern;
                     int      pp_i_region    = p_i_region;

                     // transport positron
                     if ( secondary_electron(p_pair,
                                             pp_bit_pattern, pp_i_region,
                                             rndm) )
                     {
                        // the positron survived, create new list entry
                        ++num_secondaries;
                        p_new  = NULL;
                        if ( (p_new = new (nothrow) particle_parameters)==NULL )
                        {
                           xvmc_error("MC_object::primary_photon",
                              "cannot allocate memory for pair positron",8);
                        }

                        // store parameters
                        p_new->type     = p_pair.type;
                        p_new->energy   = p_pair.energy;
                        p_new->weight   = p_pair.weight;
                        p_new->pos.x    = p_pair.pos.x;
                        p_new->pos.y    = p_pair.pos.y;
                        p_new->pos.z    = p_pair.pos.z;
                        p_new->dir.x    = p_pair.dir.x;
                        p_new->dir.y    = p_pair.dir.y;
                        p_new->dir.z    = p_pair.dir.z;
                        p_new->i.x      = p_pair.i.x;
                        p_new->i.y      = p_pair.i.y;
                        p_new->i.z      = p_pair.i.z;
                        p_new->next     = p_pair.next;

                        // assign to p_first if this is the first secondary,
                        // assign to the next list entry otherwise
                        if (p_first == NULL) p_first      = p_new;
                        else                 p_last->next = p_new;

                        // the positron is now the last particle in the list
                        p_last = p_new;
                     }
                  }

                  // stop transport
                  p_i_region = -1;

                  // because the photon has been killed
                  primary_photon_alive = false;
               }
               else
               {
                  // photoelectric absorption, transform photon into electron
                  if (p.energy > e_cut)
                  {
                     // simulate photo electron
                     e_phot.energy = p.energy;
                     e_phot.pos.x  = p.pos.x;
                     e_phot.pos.y  = p.pos.y;
                     e_phot.pos.z  = p.pos.z;
                     e_phot.dir.x  = p.dir.x;
                     e_phot.dir.y  = p.dir.y;
                     e_phot.dir.z  = p.dir.z;

                     // set photo electron region identifiers
                     MCBitSet ea_bit_pattern = p_bit_pattern;
                     int      ea_i_region    = p_i_region;

                     // transport electron
                     if ( secondary_electron(e_phot,
                                             ea_bit_pattern, ea_i_region,
                                             rndm) )
                     {
                        // the electron survived, create new list entry
                        ++num_secondaries;
                        p_new  = NULL;
                        if ( (p_new = new (nothrow) particle_parameters)==NULL )
                        {
                           xvmc_error("MC_object::primary_photon",
                              "cannot allocate memory for photo electron",8);
                        }

                        // store parameters
                        p_new->type     = e_phot.type;
                        p_new->energy   = e_phot.energy;
                        p_new->weight   = e_phot.weight;
                        p_new->pos.x    = e_phot.pos.x;
                        p_new->pos.y    = e_phot.pos.y;
                        p_new->pos.z    = e_phot.pos.z;
                        p_new->dir.x    = e_phot.dir.x;
                        p_new->dir.y    = e_phot.dir.y;
                        p_new->dir.z    = e_phot.dir.z;
                        p_new->i.x      = e_phot.i.x;
                        p_new->i.y      = e_phot.i.y;
                        p_new->i.z      = e_phot.i.z;
                        p_new->next     = e_phot.next;

                        // assign to p_first if this is the first secondary,
                        // assign to the next list entry otherwise
                        if (p_first == NULL) p_first      = p_new;
                        else                 p_last->next = p_new;

                        // the electron is now the last particle in the list
                        p_last = p_new;
                     }
                  }

                  // stop transport
                  p_i_region = -1;

                  // because the photon has been killed
                  primary_photon_alive = false;
               }
#else
               // pair production or photoelectric absorption,
               // stop transport because secondary electrons are neglected
               p_i_region = -1;

               // and kill photon
               primary_photon_alive = false;
#endif
            }
         }
         else
         {
            // move the photon to this plane
            p.pos.x += distance*p.dir.x;
            p.pos.y += distance*p.dir.y;
            p.pos.z += distance*p.dir.z;

            if (p_plane == final_plane)
            {
               // the photon is at the final plane, stop transport
               p_i_region = -1;
            }
            else
            {
               // reduce the number of mean free paths
               num_mfp -= distance*mu_tot;

               // the photon crosses the plane, check the corresponding bit
               // if the bit is set the new bit will be cleared (false)
               // if the bit is not set the new bit will be set (true)
               // that is, we invert the corresponding bit
               if (p_bit_pattern.test(p_plane->index)) p_plane_bit = false;
               else                                    p_plane_bit = true;

               // set this bit and calculate all other bits
               p_bit_pattern = calc_bit_pattern(p,p_plane,p_plane_bit);

               // determine the new region index of the photon, it should
               // be close to the present region, therefore we start the
               // region search at the present region (index "old_i_region")
               old_i_region = p_i_region;
               p_i_region = get_region_index(p_bit_pattern,old_i_region);
            }
         }
      }
      else
      {
         // the distance to the next plane is infinity
         p_i_region = -1;
      }
   }

#ifdef SECONDARY_MODIFIER_ELECTRONS
   // we sample only one of all primary and secondary particles
   if (primary_photon_alive)
   {
      // the primary photon is still alive
      if (p_first != NULL)
      {
         // there are also secondary particles
         int   num_particles = num_secondaries + 1;
         float sum_particles = float(num_particles);
         int p_index = int(rndm.number()*sum_particles); // int random number
         p_index     = max_of(0,p_index);
         p_index     = min_of(num_particles-1,p_index);
         p_new       = &p;           // initialize loop pointer
         register int i_entry = 0;   // index of list entry
         while (p_new != NULL)
         {
            if (p_index == i_entry)
            {
               // move parameters only for secondary particles,
               // change only the weight if we return the primary photon
               if (p_new == &p)
               {
                  p.weight *= sum_particles;
               }
               else
               {
                  p.type     = p_new->type;
                  p.energy   = p_new->energy;
                  p.weight   = p_new->weight*sum_particles;
                  p.pos.x    = p_new->pos.x;
                  p.pos.y    = p_new->pos.y;
                  p.pos.z    = p_new->pos.z;
                  p.dir.x    = p_new->dir.x;
                  p.dir.y    = p_new->dir.y;
                  p.dir.z    = p_new->dir.z;
               }
            }

            // old particle in linked list
            p_last = p_new;

            // next active particle
            if (p_last == &p)
            {
               p_new = p_first;
            }
            else
            {
               p_new = p_last->next;
               delete  p_last;  // delete old particle
            }

            // increase list entry index
            ++i_entry;
         }
      }

      // there is at least one surviving particle
      return(true);
   }

   if (p_first != NULL)
   {
      // there are only secondary patricles

      // we need an interger random number
      float sum_particles = float(num_secondaries);
      int p_index = int(rndm.number()*sum_particles);
      p_index     = max_of(0,p_index);
      p_index     = min_of(num_secondaries-1,p_index);
      p_new       = p_first;  // initialize loop pointer
      register int i_entry = 0;   // index of list entry
      while (p_new != NULL)
      {
         if (p_index == i_entry)
         {
            // store parameters
            p.type     = p_new->type;
            p.energy   = p_new->energy;
            p.weight   = p_new->weight*sum_particles;
            p.pos.x    = p_new->pos.x;
            p.pos.y    = p_new->pos.y;
            p.pos.z    = p_new->pos.z;
            p.dir.x    = p_new->dir.x;
            p.dir.y    = p_new->dir.y;
            p.dir.z    = p_new->dir.z;
         }

         // old particle in linked list
         p_last = p_new;

         // next active particle
         p_new  = p_last->next;

         // delete old particle
         delete p_last;

         // increase list entry index
         ++i_entry;
      }

      // there is a surviving secondary particle
      return(true);
   }

   // there is neither a primary nor a secondary surviving particle
   p.energy = ZERO;
   p.weight = ZERO;
   return(false);
#else
   if (primary_photon_alive) return(true);

   // there is no surviving particle, set energy and weight to zero
   p.energy = ZERO;
   p.weight = ZERO;
   return(false);
#endif
}

// transport primary electron through the object, the plane pointer points to
// the starting plane of the electron (it must not be the NULL pointer),
// if the electron survives (and/or daughter particles) -> return true,
// if no particle survives  -> return false,
// this is just a Continuous Slowing Down Approximation (CSDA),
// the random number generator is reserved for later use
inline bool MC_object::primary_electron(particle_parameters &e,
                                        MC_plane *e_plane, ranmar &rndm)
{
   // there must be an error if e_plane points to NULL
   if (e_plane == NULL)
   {
      xvmc_error("MC_object::primary_electron",
                 "the starting plane is undefined",8);
   }

   // determine the relationship of the electron to e_plane by taking only the
   // electron direction into account, compare with the plane normal vector
   real_3 normal = e_plane->get_normal();
            
   // scalar product of electron direction and plane normal vector
   real product = e.dir.x*normal.x + e.dir.y*normal.y + e.dir.z*normal.z;

   bool e_plane_bit;
   if (product > ZERO)
   {
      // the electron points to the region before the plane, clear bit
      e_plane_bit = false;
   }
   else
   {
      // the electron points to the region behind the plane, set bit
      // or it will stay at the plane, set also bit because
      // in this case we define the electron to be behind the plane
      e_plane_bit = true;
   }

   // determine the electron bit pattern for all planes
   // including the starting plane (set by bool e_plane_bit)
   MCBitSet e_bit_pattern = calc_bit_pattern(e,e_plane,e_plane_bit);

   // estimate and get region index of the electron
   int old_i_region = estimate_region_index(e.pos);
   int e_i_region   = get_region_index(e_bit_pattern,old_i_region);

   // the electron is outside of the object, return false
   if (e_i_region < 0)
   {
      e.energy = ZERO;
      e.weight = ZERO;
      return(false);
   }


   // now trace electron through the object
   while (e_i_region >= 0)
   {
      // distance to the next plane and next plane
      real distance = piece[e_i_region]->distance(e,e_plane);

      if (e_plane != NULL)
      {
         // CSDA range of the electron for the material of the present region
         real csda_range = piece[e_i_region]->e_data->get_csda(e.energy);

         // compare electron CSDA range and distance
         if (csda_range < distance)
         {
            // kill electron
            return(false);
         }
         else
         {
            // move electron to the next plane
            e.pos.x += distance*e.dir.x;
            e.pos.y += distance*e.dir.y;
            e.pos.z += distance*e.dir.z;

            // estimate the new energy of the electron
            real new_energy = e.energy -
                 distance*piece[e_i_region]->e_data->get_dedx(e.energy);
            if (new_energy <= e_cut) return(false);
            e.energy = new_energy;

            if (e_plane == final_plane)
            {
               // the electron is at the final plane, stop transport
               e_i_region = -1;
            }
            else
            {
               // the electron crosses the plane, check the corresponding bit
               // if the bit is set the new bit will be cleared (false)
               // if the bit is not set the new bit will be set (true)
               // that is, we invert the corresponding bit
               if (e_bit_pattern.test(e_plane->index)) e_plane_bit = false;
               else                                    e_plane_bit = true;

               // set this bit and calculate all other bits
               e_bit_pattern = calc_bit_pattern(e,e_plane,e_plane_bit);

               // determine the new region index of the electron, it should
               // be close to the present region, therefore we start the
               // region search at the present region (index "old_i_region")
               old_i_region = e_i_region;
               e_i_region = get_region_index(e_bit_pattern,old_i_region);
            }
         }
      }
      else
      {
         // the distance to the next plane is infinity
         e_i_region = -1;
      }
   }

   // the electron survived, return true
   return(true);
}

// transport secondary electron through the object,
// the bit pattern and the region index of the electron must be known,
// if the electron survives (and/or daughter particles) -> return true,
// if no particle survives  -> return false,
// this is just a Continuous Slowing Down Approximation (CSDA),
// the random number generator is reserved for later use
inline bool MC_object::secondary_electron(particle_parameters &e,
                                          MCBitSet &e_bit_pattern,
                                          int      &e_i_region,
                                          ranmar   &rndm)
{
   // there must be an error if the secondary electron is outside of the object
   if (e_i_region < 0)
   {
      xvmc_error("MC_object::secondary_electron",
                 "the electron is outside of the object",8);
   }

   MC_plane *e_plane = NULL;
   bool      e_plane_bit;
   int       old_i_region;

   // now trace electron through the object
   while (e_i_region >= 0)
   {
      // distance to the next plane and next plane
      real distance = piece[e_i_region]->distance(e,e_plane);

      if (e_plane != NULL)
      {
         // CSDA range of the electron for the material of the present region
         real csda_range = piece[e_i_region]->e_data->get_csda(e.energy);

         // compare electron CSDA range and distance
         if (csda_range < distance)
         {
            // kill electron
            return(false);
         }
         else
         {
            // move electron to the next plane
            e.pos.x += distance*e.dir.x;
            e.pos.y += distance*e.dir.y;
            e.pos.z += distance*e.dir.z;

            // estimate the new energy of the electron
            real new_energy = e.energy -
                 distance*piece[e_i_region]->e_data->get_dedx(e.energy);
            if (new_energy <= e_cut) return(false);
            e.energy = new_energy;

            if (e_plane == final_plane)
            {
               // the electron is at the final plane, stop transport
               e_i_region = -1;
            }
            else
            {
               // the electron crosses the plane, check the corresponding bit
               // if the bit is set the new bit will be cleared (false)
               // if the bit is not set the new bit will be set (true)
               // that is, we invert the corresponding bit
               if (e_bit_pattern.test(e_plane->index)) e_plane_bit = false;
               else                                    e_plane_bit = true;

               // set this bit and calculate all other bits
               e_bit_pattern = calc_bit_pattern(e,e_plane,e_plane_bit);

               // determine the new region index of the electron, it should
               // be close to the present region, therefore we start the
               // region search at the present region (index "old_i_region")
               old_i_region = e_i_region;
               e_i_region = get_region_index(e_bit_pattern,old_i_region);
            }
         }
      }
      else
      {
         // the distance to the next plane is infinity
         e_i_region = -1;
      }
   }

   // the electron survived, return true
   return(true);
}

#endif /* _MC_OBJECT_H_ */
