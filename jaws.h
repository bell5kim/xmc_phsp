#ifndef _JAWS_H_
#define _JAWS_H_

/*****************************************************************************
 * jaws.h:                                                                   *
 *    class declarations and inline member functions for:                    *
 *       class jaws: nominal jaws (iso-center plane)                         *
 *                                                                           *
 * Copyright (C) 2001    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 15.06.2001      *
 *    dynamic MLC mode implemented                        MF 09.10.2003      *
 *    modified from multi_leaf.h                          JK 29.01.2008      *
 *                                                                           *
 *****************************************************************************/

#include "definitions.h"
#include "global.h"

enum jaw_type {SIMPLE_JAW, DBLFOCUS_JAW, REAL_JAW};
enum jaw_mode {STATIC_JAW, DYNAMIC_JAW};

// nominal jaws (iso-center plane)
class jaws
{
   public:
      // define JAW by type, mode, 
      // leaf moving direction (X or Y),
      // number of jaw pairs (must be ONE), JAW material,
      // upper and lower limits of the JAW
      jaws(const jaw_type  inp_type, const jaw_mode  inp_mode,
           const char      xy,
           const int       num,      const char     *inp_material,
           const real      z_upper,  const real     z_lower ) {
           init(inp_type,inp_mode,xy,num,inp_material,z_upper,z_lower); }

      // delete JAW
      ~jaws(void);

      // reset jaw positions (zero beam opening)
      void reset(void);

      // change width of jaws,
      // if this is impossible -> return false
      bool change_width(real);

      // change position of jaw pair,
      // if this is impossible -> return false
      bool change_jaws(real, real);

      // change starting position of jaw pair,
      // if this is impossible -> return false
      bool change_start_jaws(real, real);

      // change stopping position of jaw pair,
      // if this is impossible -> return false
      bool change_stop_jaws(real, real);

      // get JAW type
      jaw_type get_type(void) { return(type); }

      // get JAW mode
      jaw_mode get_mode(void) { return(mode); }

      // set JAW mode
      void set_mode(jaw_mode new_mode) { mode = new_mode; }

      // get jaw moving direction X or Y
      char get_xy(void) { return(x_jaw ? 'X' : 'Y'); }

      // get JAW material
      char *get_material(void) { return(material); }

      // get upper JAW limit
      real get_upper(void) { return(upper_limit); }

      // get lower JAW limit
      real get_lower(void) { return(lower_limit); }

      // get number of jaw pair (it must be ONE)
      int get_num(void) {  return(num_pairs); }

      // get width of jaw
      real get_width() { return(width); }
      // get perpendicular jaw limit 
      real get_length() { return(length); } 
      // get position of left jaw
      real get_left() { return(left); }
      // get position of right jaw
      real get_right() { return(right); }
      // get starting position of left jaw
      real get_start_left() { return(start_left); }
      // get starting position of right jaw  
      real get_start_right() { return(start_right); }
      // get stopping position of left jaw 
      real get_stop_left() { return(stop_left); }
      // get stopping position of right jaw 
      real get_stop_right(){ return(stop_right); }

      // for a point (x0,y0) get the corresponding jaw pair index,
      // return -1 if the point is outside the JAW
      int get_index(real, real);

      // return true if a point (x0,y0) is located inside the beam limits
      bool inside(real, real);

   private:
      // JAW type (SIMPLE_JAW, DBLFOCUS_JAW, REAL_JAW, etc.)
      jaw_type  type;

      // JAW mode (static or dynamic)
      jaw_mode  mode;

      // jaw moving direction: true  if the jaw are parallel to the X axis
      //                        false if the jaw are parallel to the Y axis
      bool      x_jaw;

      // JAW material
      char *material;

      // upper and lower JAW limits
      real upper_limit, lower_limit;

      int       num_pairs;       // number of jaw pairs
      real      width;           // jaw widths
      real      length;          // perpendicular jaw limits
      real      left;            // left jaw positions
      real      right;           // right jaw positions
      real      min_width;       // jaw position (max. opening)
      real      max_width;       // jaw position (max. opening)
      real      start_left;      // jaw starting positions
      real      start_right;     // jaw starting positions
      real      min_start_left;  // jaw starting position
      real      max_start_right; // jaw starting position
      real      stop_left;       // jaw stopping positions
      real      stop_right;      // jaw stopping positions
      real      min_stop_left;   // minimum left  jaw stopping position
      real      max_stop_right;  // maximum right jaw stopping position

      // minimum and maximum beam openings perpendicular to the jaw direction
      real      min_length,max_length;

      // initialize JAW by type, mode, jaw moving direction (X or Y),
      // number of jaw pairs, JAW material,
      // upper and lower limits of the JAW
      void init(const jaw_type, const jaw_mode,
                const char,
                const int,      const char *,
                const real,     const real);
};

// return true if a point (x0,y0) is located inside the beam limits
inline bool jaws::inside(real x0, real y0)
{
   if (x_jaw) // X jaw
   {
      // test, whether (x0,y0) is outside the rectangle defined by
      // min_width,max_width,min_length,max_length (to save compuation time)
      if (x0 < min_width)   return(false);
      if (x0 > max_width)   return(false);
      if (y0 < min_length)  return(false);
      if (y0 > max_length)  return(false);

      // check whether left < x0 < right
      if ( (left < x0) && (x0 < right) ) return(true);
   }
   else         // Y jaw
   {
      // test, whether (x0,y0) is outside the rectangle defined by
      // min_left,max_right,min_length,max_length (to save compuation time)
      if (x0 < min_length)  return(false);
      if (x0 > max_length)  return(false);
      if (y0 < min_width)   return(false);
      if (y0 > max_width)   return(false);

      // check whether left < y0 < right
      if ( (left < y0) && (y0 < right) ) return(true);
   }

   // the point (x0,y0) is outside
   return(false);
}

#endif /* _JAWS_H_ */
