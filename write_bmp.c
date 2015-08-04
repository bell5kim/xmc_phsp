/* ========================================================================= *
 *                                                                           *
 * write_bmp.c: subroutine to write a BMP file                               *
 *                                                                           *
 * Copyright (C) 1997    Version 2.0: Iwan Kawrakow, Matthias Fippel         *
 * Copyright (C) 1998    Version 2.1: Iwan Kawrakow, Matthias Fippel         *
 *                       Clinic for Radiation Therapy,                       *
 *                       University of Leipzig, Germany                      *
 *                       Germany                                             *
 * Copyright (C) 1998    Version 3.0: Matthias Fippel                        *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * ------------------------------------------------------------------------- *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the Free Software             *
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               *
 *                                                                           *
 * ------------------------------------------------------------------------- *
 *                                                                           *
 *  Contacts:                                                                *
 *                                                                           *
 *  Matthias Fippel                             Iwan Kawrakow                *
 *  Abteilung fuer Medizinische Physik          Ionizing Radiation Standards *
 *  Radiologische Universitaetsklinik           Institute for National       *
 *  Hoppe-Seyler-Str. 3                         Measurement Standards        *
 *  72076 Tuebingen                             National Research Council    *
 *  Germany                                     Ottawa, K1A 0R6, Canada      *
 *  voice: ++49 7071 2985887                    voice: ++1 613 993 2197      *
 *  fax:   ++49 7071 295920                     fax:   ++1 613 952 9865      *
 *  email: msfippel@med.uni-tuebingen.de        email: iwan@irs.phy.nrc.ca   *
 *                                                                           *
 * ------------------------------------------------------------------------- *
 *                                                                           *
 *   Revision:   96/05/09    M. Fippel  Initial coding                       *
 *               97/01/09    M. Fippel  Various changes for Version 2.0      *
 *                                                                           *
 * ========================================================================= */

#include <stdio.h>
#include <stdlib.h>

#ifdef AIX
#  define SWAP_BYTE_ORDER
#endif

#ifdef IRX
#  define SWAP_BYTE_ORDER
#endif

#ifdef HPUX
#  define SWAP_BYTE_ORDER
#endif

#ifdef CONVEX
#  define SWAP_BYTE_ORDER
#endif

typedef char            BYTE;
typedef unsigned short  UINT2;
typedef signed short    SINT2;
typedef unsigned int    UINT4;
typedef signed int      SINT4;

typedef struct
        {  UINT4 bfSize;           /* size of the file */
           UINT2 bfReserved1;      /* reserved = 0  */
           UINT2 bfReserved2;      /* reserved = 0  */
           UINT4 bfOffBits;        /* Offset from BITMAPFILEHEADER */
                                   /* to actual bitmap in file = 1078 */
        } BITMAPFILEHEADER;

typedef struct
        {  UINT4  biSize ;         /* size of BITMAPINFOHEADER = 40 */
           UINT4  biWidth;         /* width of bitmap (pixel) */
           UINT4  biHeight;        /* height of bitmap (pixel) */
           UINT2  biPlanes ;       /* number of color planes = 1 */
           UINT2  biBitCount;      /* number of bits/pixel = 8 */
           UINT4  biCompression;   /* compression type = 0 */
           UINT4  biSizeImage;     /* image size (bytes) */
           UINT4  biXPelsPerMeter; /* x-resolution = 0 */
           UINT4  biYPelsPerMeter; /* y-resolution = 0 */
           UINT4  biClrUsed;       /* = 0 */
           UINT4  biClrImportant;  /* = 0 */
        }
        BITMAPINFOHEADER;

typedef  struct
        {
           BYTE   rgbBlue;
           BYTE   rgbGreen;
           BYTE   rgbRed;
           BYTE   rgbReserved;
        }
            RGBQUAD;


typedef struct
        {
           BITMAPFILEHEADER   bf;
           BITMAPINFOHEADER   bi;
           RGBQUAD            rgbquad[256];
        }BMPHEADER;

/*************************  DECLARE FUNCTIONS  ***************************/

UINT2 mirror_uint2( UINT2 );
UINT4 mirror_uint4( UINT4 );

void write_bmp( const char *bmp_file_name,
                int max_x, int max_y,
                int dwX,   int dwY,
                BYTE *pl_array)
{
UINT2        bmp_type;
UINT4        offset;
UINT4        bmp_x,bmp_y;
UINT4        image_size;
UINT4        file_size;
register int i,j;
UINT4        addr1,addr2;
FILE        *bmp_file;
BMPHEADER   *head;
BYTE        *image;


/*** create bmp header ***/

  offset     = sizeof(BMPHEADER);
  bmp_x      = ((UINT4) dwX + 3)/4 * 4;
  bmp_y      =  (UINT4) dwY;
  image_size = bmp_x * bmp_y;
  file_size  = offset + image_size + 2;

  if ( ( head=malloc(offset) ) == NULL )
  {
    printf("\n ERROR, cannot allocate space for header! \n");
    exit(-1);
  }

  bmp_type                 = 0x4d42;

  head->bf.bfSize          = file_size;
  head->bf.bfReserved1     = (UINT2) 0x0000;
  head->bf.bfReserved2     = (UINT2) 0x0000;
  head->bf.bfOffBits       = offset+2;

  head->bi.biSize          = sizeof(BITMAPINFOHEADER);
  head->bi.biWidth         = (UINT4) dwX;
  head->bi.biHeight        = (UINT4) dwY;
  head->bi.biPlanes        = (UINT2) 1;
  head->bi.biBitCount      = (UINT2) 8;
  head->bi.biCompression   = (UINT4) 0x00;
  head->bi.biSizeImage     = (UINT4) image_size;
  head->bi.biXPelsPerMeter = (UINT4) 0x00;
  head->bi.biYPelsPerMeter = (UINT4) 0x00;
  head->bi.biClrUsed       = (UINT4) 0x00;
  head->bi.biClrImportant  = (UINT4) 0x00;

/*** mirror bytes ***/

#ifdef SWAP_BYTE_ORDER
   bmp_type                 = mirror_uint2( bmp_type );

   head->bf.bfSize          = mirror_uint4( head->bf.bfSize );
   head->bf.bfReserved1     = mirror_uint2( head->bf.bfReserved1 );
   head->bf.bfReserved2     = mirror_uint2( head->bf.bfReserved2 );
   head->bf.bfOffBits       = mirror_uint4( head->bf.bfOffBits );

   head->bi.biSize          = mirror_uint4( head->bi.biSize );
   head->bi.biWidth         = mirror_uint4( head->bi.biWidth );
   head->bi.biHeight        = mirror_uint4( head->bi.biHeight );
   head->bi.biPlanes        = mirror_uint2( head->bi.biPlanes );
   head->bi.biBitCount      = mirror_uint2( head->bi.biBitCount );
   head->bi.biCompression   = mirror_uint4( head->bi.biCompression );
   head->bi.biSizeImage     = mirror_uint4( head->bi.biSizeImage );
   head->bi.biXPelsPerMeter = mirror_uint4( head->bi.biXPelsPerMeter );
   head->bi.biYPelsPerMeter = mirror_uint4( head->bi.biYPelsPerMeter );
   head->bi.biClrUsed       = mirror_uint4( head->bi.biClrUsed );
   head->bi.biClrImportant  = mirror_uint4( head->bi.biClrImportant );
#endif

/*** create color table ***/
  
  for ( i=0; i<256; i++)
  {
    head->rgbquad[i].rgbBlue     = (BYTE) i;
    head->rgbquad[i].rgbGreen    = (BYTE) i;
    head->rgbquad[i].rgbRed      = (BYTE) i;
    head->rgbquad[i].rgbReserved = (BYTE) 0;
  }

/*
  head->rgbquad[246].rgbBlue     = (BYTE) 0;
  head->rgbquad[246].rgbGreen    = (BYTE) 127;
  head->rgbquad[246].rgbRed      = (BYTE) 255;
  head->rgbquad[246].rgbReserved = (BYTE) 0;

  head->rgbquad[247].rgbBlue     = (BYTE) 255;
  head->rgbquad[247].rgbGreen    = (BYTE) 127;
  head->rgbquad[247].rgbRed      = (BYTE) 127;
  head->rgbquad[247].rgbReserved = (BYTE) 0;

  head->rgbquad[248].rgbBlue     = (BYTE) 127;
  head->rgbquad[248].rgbGreen    = (BYTE) 255;
  head->rgbquad[248].rgbRed      = (BYTE) 127;
  head->rgbquad[248].rgbReserved = (BYTE) 0;

  head->rgbquad[249].rgbBlue     = (BYTE) 127;
  head->rgbquad[249].rgbGreen    = (BYTE) 127;
  head->rgbquad[249].rgbRed      = (BYTE) 255;
  head->rgbquad[249].rgbReserved = (BYTE) 0;

  head->rgbquad[250].rgbBlue     = (BYTE) 255;
  head->rgbquad[250].rgbGreen    = (BYTE) 255;
  head->rgbquad[250].rgbRed      = (BYTE) 0;
  head->rgbquad[250].rgbReserved = (BYTE) 0;

  head->rgbquad[251].rgbBlue     = (BYTE) 255;
  head->rgbquad[251].rgbGreen    = (BYTE) 0;
  head->rgbquad[251].rgbRed      = (BYTE) 255;
  head->rgbquad[251].rgbReserved = (BYTE) 0;

  head->rgbquad[252].rgbBlue     = (BYTE) 0;
  head->rgbquad[252].rgbGreen    = (BYTE) 255;
  head->rgbquad[252].rgbRed      = (BYTE) 255;
  head->rgbquad[252].rgbReserved = (BYTE) 0;

  head->rgbquad[253].rgbBlue     = (BYTE) 255;
  head->rgbquad[253].rgbGreen    = (BYTE) 0;
  head->rgbquad[253].rgbRed      = (BYTE) 0;
  head->rgbquad[253].rgbReserved = (BYTE) 0;

  head->rgbquad[254].rgbBlue     = (BYTE) 0;
  head->rgbquad[254].rgbGreen    = (BYTE) 255;
  head->rgbquad[254].rgbRed      = (BYTE) 0;
  head->rgbquad[254].rgbReserved = (BYTE) 0;

  head->rgbquad[255].rgbBlue     = (BYTE) 0;
  head->rgbquad[255].rgbGreen    = (BYTE) 0;
  head->rgbquad[255].rgbRed      = (BYTE) 255;
  head->rgbquad[255].rgbReserved = (BYTE) 0;
*/

/*** open BMP file ***/
  
  if ( ( bmp_file=fopen( bmp_file_name,"wb") ) == NULL )
  {
    printf("\n ERROR, cannot open BMP file! \n");
    exit(-1);
  }

/*** write header to BMP file ***/
  if ( 1 != fwrite( &bmp_type , sizeof(UINT2), 1, bmp_file) )
  {
    printf("\n ERROR, cannot write BMP type! \n");
    exit(-1);
  }

  if ( offset != fwrite( (BYTE *)head, sizeof(BYTE),
                         offset, bmp_file) )
  {
    printf("\n ERROR, cannot write BMP header! \n");
    exit(-1);
  }
  
/*** allocate space for image ***/

  if ( ( image=malloc(image_size) ) == NULL )
  {
    printf("\n ERROR, cannot allocate space for image! \n");
    exit(-1);
  }

/*** write image to memory ***/

  for ( j=0; j< dwY; j++)
  {
    addr1 = bmp_x * j;
    addr2 = max_x * ( dwY - j - 1 );
    for ( i=0; i < (int) bmp_x; i++)
    {
      if ( i < dwX )
      {
        *(image + addr1 + i) = *(pl_array + addr2 + i);
      }
      else
      {
        *(image + addr1 + i) = (BYTE) 0x00;
      }
    }
  }

/*** write image to BMP file ***/

  if ( image_size != fwrite( (BYTE *)image, sizeof(BYTE),
                              image_size, bmp_file) )
  {
    printf("\n ERROR, cannot write image! \n");
    exit(-1);
  }

/*** close BMP file and free memory ***/

  fclose( bmp_file );
  free( head );
  free( image );
  
/*  printf("\n file: %s written! \n", bmp_file_name); */

  return;
  
}

/*-------------------------------------------------------------------------
 *      mirror_uint2
 *------------------------------------------------------------------------*/
UINT2 mirror_uint2( UINT2 vuint2 )
{
  BYTE c1,c2;
  c2 = (BYTE) vuint2;
  c1 = (BYTE) (vuint2 >> 8);
  return( ((UINT2) c1) + (((UINT2) c2) << 8) );
}

/*-------------------------------------------------------------------------
 *      mirror_uint4
 *------------------------------------------------------------------------*/
UINT4 mirror_uint4( UINT4 vuint4 )
{
  BYTE c1,c2,c3,c4;
  c4 = (BYTE) vuint4;
  c3 = (BYTE) (vuint4 >> 8);
  c2 = (BYTE) (vuint4 >> 16);
  c1 = (BYTE) (vuint4 >> 24);
  return(  ((UINT4) c1)        +
          (((UINT4) c2) <<  8) +
          (((UINT4) c3) << 16) +
          (((UINT4) c4) << 24) );
}


