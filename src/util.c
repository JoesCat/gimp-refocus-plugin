/* Refocus plug-in
 * Copyright (C) 1999-2004... Ernst Lippe - (original author)
 * Copyright (C) 2024 Jose Da Silva (updates and improvements)
 *
 * Based on the Convolution Matrix plug-in by Lauri Alanko <la@iki.fi>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <string.h>
#include <libgimp/gimp.h>
#include <libgimp/gimpui.h>
#include "util.h"

gint
floorm (gint a, gint b)
 /* return largest multiple of b that is <= a */
 /*
  & & m = floorm(a,b)
    & a = b*m + r
    &  0 <= r < b
  */
{
#ifdef RLXTEST
  printf("floorm: a/b %d, fl %g\n", a/b, floor ((gdouble) a / b));
#endif
  return (b * floor ((gdouble) a / b));
}

gint
ceilm (gint a, gint b)
 /* return least multiple of b that is >= a */
 /*
    & m = ceilm(a,b)
    & a = b*m - r;
    & m = a/b
    % r = a%b
    & -a = -b*m + r

    & ceilm = (r == 0 ? b*m : (b+1)*m)
  */
{
#ifdef RLXTEST
  printf("ceil: a %d, b %d, -(-a/b) %d,a/b+(a%b != 0 ? 1:0) %d,  fl %g\n",
         a,b,
         -((-a)/b),
          a/b+(a%b != 0 ? 1:0),
         ceil ((gdouble) a / b) );
#endif
  return (b * ceil ((gdouble) a / b));
}


void
copy_rect (guchar * dest_buf, gint dest_x, gint dest_y,
           gint dest_width, gint dest_height,
           guchar * src_buf, gint src_x, gint src_y,
           gint src_width, gint src_height, gint bpp)
{
  gint x_lo, x_hi, y_lo, y_hi, y;

  x_lo = MAX (src_x, dest_x);
  x_hi = MIN (src_x + src_width, dest_x + dest_width);
  y_lo = MAX (src_y, dest_y);
  y_hi = MIN (src_y + src_height, dest_y + dest_height);
  if (x_hi > x_lo)
    {
      for (y = y_lo; y < y_hi; y++)
        {
          memmove (dest_buf +
                   ((y - dest_y) * dest_width + x_lo - dest_x) * bpp,
                   src_buf + ((y - src_y) * src_width + x_lo - src_x) * bpp,
                   (x_hi - x_lo) * bpp);
        };
    };
}

/**
 * tile_width:
 * 
 * This function converts the guint that is returned by
 * gimp_tile_width into a gint. This helps prevent
 * all kind of signed/unsigned problems in the code.
 * 
 * Return value: The width of a tile.
 **/
gint
tile_width (void)
{
  return (gimp_tile_width ());
}

/**
 * tile_height:
 * 
 * This function converts the guint that is returned by
 * gimp_tile_height into a gint. This helps prevent
 * all kind of signed/unsigned problems in the code.
 * 
 * Return value: The height of a tile.
 **/
gint
tile_height (void)
{
  return (gimp_tile_height ());
}

/**
 * get_pixel:
 * @ptr: pointer to pixel_color value in source image.
 * @bpc: color value type as defined from refocus.c:run().
 * This function gets a color pixel from the source image
 * and returns a converted double[0.0..1.0] value.
 **/
gdouble
get_pixel (guchar *ptr, gint bpc)
{
  /* Returned value is in the range of [0.0..1.0]. */
  gdouble ret = 0.0;
  if (bpc == 1) {
    ret += *ptr;
    ret /= 255;
  } else if (bpc == 2) {
    uint16_t *p = (uint16_t *)(ptr);
    ret += *p;
    ret /= 65535;
   } else if (bpc == 4) {
    uint32_t *p = (uint32_t *)(ptr);
    ret += *p;
    ret /= 4294967295;
  } else if (bpc == 8) {
    uint64_t *p = (uint64_t *)(ptr);
    long double lret = 0.0;
    lret += *p;
    lret /= 18446744073709551615UL;
    ret = lret;
  } else if (bpc == -8) {
    double *p = (double *)(ptr);
    ret += *p;
  } else if (bpc == -4) {
    float *p = (float *)(ptr);
    ret += *p;
//} else if (bpc == -2) {
//  half *p = ptr;
//  ret += *p;
  }
  return ret;
}

/**
 * set_pixel:
 * @dest: pointer to destination image to set color pixel.
 * @d: input pixel_color value in range [0.0..1.0].
 * @bpc: color value type as defined from refocus.c:run().
 * This function sets a color pixel in destination image.
 * The input value is a double in the range [0.0..1.0].
 **/
void
set_pixel (guchar *dest, gdouble d, gint bpc)
{
  /* input value is in the range of [0.0..1.0]. */
  if (bpc == 1) {
    *dest = round (d * 255);
  }  else if (bpc == 2) {
    uint16_t *p = (uint16_t *)(dest);
    *p = round (d * 65535);
  } else if (bpc == 4) {
    uint32_t *p = (uint32_t *)(dest);
    *p = round (d * 4294967295);
  } else if (bpc == 8) {
    uint64_t *p = (uint64_t *)(dest);
    *p = roundl (d * 18446744073709551615UL);
  } else if (bpc == -8) {
    double *p = (double *)(dest);
    *p = d;
  } else if (bpc == -4) {
    float *p = (float *)(dest);
    *p = (float)(d);
//} else if (bpc == -2) {
//  half *p = (half *)(dest);
//  *p = d;
  }
  return;
}
