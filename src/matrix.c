/* Refocus plug-in
 * Copyright (C) 1999-2004... Ernst Lippe - (original author)
 * Copyright (C) 2024 Jose Da Silva (updates and improvements)
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

#include <libgimp/gimp.h>
#include <math.h>

#include "matrix.h"
#include "fwlapack.h"
#include "util.h"

static const char *OUT_OF_MEMORY = "Error: Ran out of RAM memory.\n";

static Mat *allocate_matrix (gint nrows, gint ncols) {
  Mat *result;

  if (!(result = g_malloc0 (sizeof (Mat)))) {
    fprintf(stderr, OUT_OF_MEMORY);
  } else {
    result->rows = nrows;
    result->cols = ncols;
    if (!(result->data = g_new0 (double, nrows * ncols))) {
      g_free (result);
      result = NULL;
      fprintf(stderr, OUT_OF_MEMORY);
    }
  }
  return (result);
}

static void finish_and_free_matrix (Mat * mat) {
  g_free (mat->data);
  g_free (mat);
}

static double *mat_eltptr (Mat * mat, const gint r, const gint c) {
#ifdef RF_DEBUG
  g_assert ((r >= 0) && (r < mat->rows));
  g_assert ((c >= 0) && (c < mat->cols));
#endif
  return (&(mat->data[mat->rows * c + r]));
}

static double mat_elt (const Mat * mat, const gint r, const gint c) {
#ifdef RF_DEBUG
  g_assert ((r >= 0) && (r < mat->rows));
  g_assert ((c >= 0) && (c < mat->cols));
#endif
  return (mat->data[mat->rows * c + r]);
}

int init_c_mat (CMat * mat, const gint radius) {
  mat->radius = radius;
  mat->row_stride = 2 * radius + 1;
  if ((mat->data = g_new0 (double, SQR (mat->row_stride))) == NULL) {
    fprintf(stderr, OUT_OF_MEMORY);
    return (0);
  }
  mat->center = mat->data + mat->row_stride * mat->radius + mat->radius;
  return (1);
}

static CMat *allocate_c_mat (const gint radius) {
  CMat *result;

  if ((result = g_malloc0 (sizeof (CMat)))) {
    if (!(init_c_mat (result, radius))) {
      g_free (result);
      fprintf(stderr, OUT_OF_MEMORY);
      result = NULL;
    }
  }
  return (result);
}

void finish_c_mat (CMat * mat) {
  g_free (mat->data);
  mat->data = NULL;
}

static double *c_mat_eltptr (CMat * mat, const gint col, const gint row) {
#ifdef RF_DEBUG
  g_assert ((ABS (row) <= mat->radius) && (ABS (col) <= mat->radius));
#endif
  return (mat->center + mat->row_stride * row + col);
}

static int fill_matrix (CMat * matrix, const gint m,
             double f (const gint, const gint, const double), const double fun_arg) {
  register gint x, y;

  if (!(init_c_mat (matrix, m)))
    return (0);

  // This is a symetric fill when looking at the results.
  // Speed-up this fill by calculating once and copy 4x.
  //for (y = -m; y <= m; y++) {
  //  for (x = -m; x <= m; x++) {
  //    *c_mat_eltptr (matrix, x, y) = f (x, y, fun_arg);
  //  }
  for (y = -m; y <= 0; y++) {
    for (x = -m; x <= 0; x++) {
      *c_mat_eltptr (matrix, -x, -y) =
      *c_mat_eltptr (matrix, -x,  y) =
      *c_mat_eltptr (matrix,  x, -y) =
      *c_mat_eltptr (matrix,  x,  y) = f (x, y, fun_arg);
    }
  }
  return (1);
}

static int fill_matrix2 (CMat * matrix, const gint m,
              double f (const gint, const gint, const double, const double),
              const double fun_arg1, const double fun_arg2) {
  register gint x, y;

  if (!(init_c_mat (matrix, m)))
    return (0);

  // This is a symetric fill when looking at the results.
  // Speed-up this fill by calculating once and copy 4x.
  //for (y = -m; y <= m; y++) {
  //  for (x = -m; x <= m; x++) {
  //    *c_mat_eltptr (matrix, x, y) = f (x, y, fun_arg1, fun_arg2);
  //  }
  for (y = -m; y <= 0; y++) {
    for (x = -m; x <= 0; x++) {
      *c_mat_eltptr (matrix, -x, -y) =
      *c_mat_eltptr (matrix, -x,  y) =
      *c_mat_eltptr (matrix,  x, -y) =
      *c_mat_eltptr (matrix,  x,  y) = f (x, y, fun_arg1, fun_arg2);
    }
  }
  return (1);
}

static double c_mat_elt (const CMat * const mat, const gint col, const gint row) {
#ifdef RF_DEBUG
  g_assert ((ABS (row) <= mat->radius) && (ABS (col) <= mat->radius));
#endif
  return (mat->center[mat->row_stride * row + col]);
}

static void convolve_mat (CMat * result, const CMat * const mata, const CMat * const matb) {
  register gint xr, yr, xa, ya;

  for (yr = -result->radius; yr <= result->radius; yr++) {
    const gint ya_low  = MAX (-mata->radius, yr - matb->radius);
    const gint ya_high = MIN ( mata->radius, yr + matb->radius);
    for (xr = -result->radius; xr <= result->radius; xr++) {
      const gint xa_low  = MAX (-mata->radius, xr - matb->radius);
      const gint xa_high = MIN ( mata->radius, xr + matb->radius);
      register double val = 0.0;

      for (ya = ya_low; ya <= ya_high; ya++) {
        for (xa = xa_low; xa <= xa_high; xa++) {
          val += c_mat_elt (mata, xa, ya) *
                 c_mat_elt (matb, xr - xa, yr - ya);
        }
      }
      *c_mat_eltptr (result, xr, yr) = val;
    }
  }
}

void convolve_star_mat (CMat * result, const CMat * const mata,
                   const CMat * const matb) {
  register gint xr, yr, xa, ya;

  for (yr = -result->radius; yr <= result->radius; yr++) {
    const gint ya_low  = MAX (-mata->radius, -matb->radius - yr);
    const gint ya_high = MIN ( mata->radius,  matb->radius - yr);
    for (xr = -result->radius; xr <= result->radius; xr++) {
      const gint xa_low  = MAX (-mata->radius, -matb->radius - xr);
      const gint xa_high = MIN ( mata->radius,  matb->radius - xr);
      register double val = 0.0;

      for (ya = ya_low; ya <= ya_high; ya++) {
        for (xa = xa_low; xa <= xa_high; xa++) {
          val += c_mat_elt (mata, xa, ya) *
                 c_mat_elt (matb, xr + xa, yr + ya);
        }
      }
      *c_mat_eltptr (result, xr, yr) = val;
    }
  }
}

static gint as_idx (const gint k, const gint l, const gint m) {
  return ((k + m) * (2 * m + 1) + (l + m));
}

static gint as_cidx (const gint k, const gint l) {
  const gint a = MAX (ABS (k), ABS (l));
  const gint b = MIN (ABS (k), ABS (l));
  return ((a * (a + 1)) / 2 + b);
}

#ifdef RF_DEBUG
void print_c_mat (FILE * file, const CMat * const mat) {
  register gint x, y;

  for (y = -mat->radius; y <= mat->radius; y++) {
    for (x = -mat->radius; x <= mat->radius; x++)
      fprintf (file, " %g", c_mat_elt (mat, x, y));
    fprintf (file, "\n");
  }
}

static void print_matrix (FILE * file, Mat * matrix) {
  gint col_idx, row_idx;

  for (row_idx = 0; row_idx < matrix->rows; row_idx++) {
    for (col_idx = 0; col_idx < matrix->cols; col_idx++)
      fprintf (file, " %g", mat_elt (matrix, row_idx, col_idx));
    fprintf (file, "\n");
  }
}
#endif

static Mat *make_s_matrix (CMat * mat, gint m, double noise_factor) {
  const gint mat_size = SQR (2 * m + 1);
  register gint yr, yc, xr, xc;
  Mat *result;

  if (!(result = allocate_matrix (mat_size, mat_size)))
    return (NULL);

  for (yr = -m; yr <= m; yr++) {
    for (xr = -m; xr <= m; xr++) {
      for (yc = -m; yc <= m; yc++) {
        for (xc = -m; xc <= m; xc++) {
          *mat_eltptr (result, as_idx (xr, yr, m), as_idx (xc, yc, m)) =
                       c_mat_elt (mat, xr - xc, yr - yc);
          if ((xr == xc) && (yr == yc)) {
            *mat_eltptr (result, as_idx (xr, yr, m),
                         as_idx (xc, yc, m)) += noise_factor;
          }
        }
      }
    }
  }
  return (result);
}

static Mat *make_s_cmatrix (CMat * mat, gint m, double noise_factor) {
  const gint mat_size = as_cidx (m + 1, 0);
  Mat *result;
  register gint yr, yc, xr, xc;

  if ((result = allocate_matrix (mat_size, mat_size))) {
    for (yr = 0; yr <= m; yr++) {
      for (xr = 0; xr <= yr; xr++) {
        for (yc = -m; yc <= m; yc++) {
          for (xc = -m; xc <= m; xc++) {
            *mat_eltptr (result, as_cidx (xr, yr), as_cidx (xc, yc)) +=
                         c_mat_elt (mat, xr - xc, yr - yc);
            if ((xr == xc) && (yr == yc)) {
                      *mat_eltptr (result, as_cidx (xr, yr),
                                   as_cidx (xc, yc)) += noise_factor;
            }
          }
        }
      }
    }
  }
  return (result);
}


static double correlation (const gint x, const gint y, const double gamma, const double musq) {
  return (musq + pow (gamma, sqrt (SQR (x) + SQR (y))));
}

static Mat *copy_vec (const CMat * const mat, const gint m) {
  Mat *result;
  register gint x, y, index = 0;

  if ((result = allocate_matrix (SQR (2 * m + 1), 1))) {
    for (y = -m; y <= m; y++) {
      for (x = -m; x <= m; x++) {
        *mat_eltptr (result, index, 0) = c_mat_elt (mat, x, y);
        index++;
      }
    }
    g_assert (index == SQR (2 * m + 1));
  }
  return (result);
}

Mat *copy_cvec (const CMat * const mat, const gint m) {
  Mat *result;
  register gint x, y, index = 0;

  if ((result = allocate_matrix (as_cidx (m + 1, 0), 1))) {
    for (y = 0; y <= m; y++) {
      for (x = 0; x <= y; x++) {
        *mat_eltptr (result, index, 0) = c_mat_elt (mat, x, y);
        index++;
      }
    }
    g_assert (index == as_cidx (m + 1, 0));
  }
  return (result);
}

static CMat *copy_cvec2mat (const Mat * const cvec, const gint m) {
  CMat *result;
  register gint x, y;

  if ((result = allocate_c_mat (m))) {
    for (y = -m; y <= m; y++) {
      for (x = -m; x <= m; x++) {
        *c_mat_eltptr (result, x, y) = mat_elt (cvec, as_cidx (x, y), 0);
      }
    }
  }
  return (result);
}

static CMat *copy_vec2mat (const Mat * const cvec, const gint m) {
  CMat *result;
  register gint x, y;

  if ((result = allocate_c_mat (m))) {
    for (y = -m; y <= m; y++) {
      for (x = -m; x <= m; x++) {
        *c_mat_eltptr (result, x, y) = mat_elt (cvec, as_idx (x, y, m), 0);
      }
    }
  }
  return (result);
}

static CMat *compute_g (const CMat * const convolution, const gint m, const double gamma,
                        const double noise_factor, const double musq, const gboolean symmetric) {
  CMat h_conv_ruv, a, corr;
  CMat *result = NULL;
  Mat *b;
  Mat *s;
  int status;

  if (!(init_c_mat (&h_conv_ruv, 3 * m)))
    goto compute_g_err_h_conv;
  if (!(fill_matrix2 (&corr, 4 * m, correlation, gamma, musq)))
    goto compute_g_err_corr;
  convolve_mat (&h_conv_ruv, convolution, &corr);
  if (!(init_c_mat (&a, 2 * m)))
    goto compute_g_err_a;
  convolve_star_mat (&a, convolution, &h_conv_ruv);
  if (symmetric) {
      if ((s = make_s_cmatrix (&a, m, noise_factor)) == NULL)
        goto compute_g_err_s;
      if ((b = copy_cvec (&h_conv_ruv, m)) == NULL)
        goto compute_g_err_b;
  } else {
      if ((s = make_s_matrix (&a, m, noise_factor)) == NULL)
        goto compute_g_err_s;
      if ((b = copy_vec (&h_conv_ruv, m)) == NULL)
        goto compute_g_err_b;
  }
#ifdef RF_DEBUG
  fprintf (stderr, "Convolution:\n");
  print_c_mat (stderr, convolution);
  fprintf (stderr, "h_conv_ruv:\n");
  print_c_mat (stderr, &h_conv_ruv);
  fprintf (stderr, "Value of s:\n");
  print_matrix (stderr, s);
  fprintf (stderr, "\n");
#endif

  g_assert (s->cols == s->rows);
  g_assert (s->rows == b->rows);
  status = dgesv (s->rows, 1, s->data, s->rows, b->data, b->rows);
#ifdef RF_DEBUG
  fprintf (stderr, "dgesv status %d\n\n", status);
#endif

  if (symmetric) {
    result = copy_cvec2mat (b, m);
  } else {
    result = copy_vec2mat (b, m);
  };
#ifdef RF_DEBUG
  fprintf (stderr, "Deconvolution:\n");
  print_c_mat (stderr, result);
  fprintf (stderr, "\n");
#endif

  finish_and_free_matrix (b);
compute_g_err_b:
  finish_and_free_matrix (s);
compute_g_err_s:
  finish_c_mat (&a);
compute_g_err_a:
  finish_c_mat (&corr);
compute_g_err_corr:
  finish_c_mat (&h_conv_ruv);
compute_g_err_h_conv:
  return (result);
}

CMat *compute_g_matrix (const CMat * const convolution, const gint m,
                        const double gamma, const double noise_factor,
                        const double musq, const gboolean symmetric) {
  CMat *g;
  gint r, c;
  double sum = 0.0;

  if ((g = compute_g (convolution, m, gamma, noise_factor, musq, symmetric))) {
    /* Determine sum of array */
    for (r = -g->radius; r <= g->radius; r++) {
      for (c = -g->radius; c <= g->radius; c++) {
        sum += c_mat_elt (g, r, c);
      }
    }
#ifdef RF_DEBUG
  fprintf (stderr, "compute_g_matrix sum %g\n", sum);
#endif
    for (r = -g->radius; r <= g->radius; r++) {
      for (c = -g->radius; c <= g->radius; c++) {
        *c_mat_eltptr (g, r, c) /= sum;
      }
    }
  }
  return (g);
}

int make_gaussian_convolution (const double gradius, CMat * convolution,
                           const gint m) {
  register gint x, y;

  if (!(init_c_mat (convolution, m)))
    return (0);

  if (SQR (gradius) <= 1 / G_MAXFLOAT) {
    for (y = -m; y <= m; y++) {
      for (x = -m; x <= m; x++) {
        *c_mat_eltptr (convolution, x, y) = 0;
      }
    }
    *c_mat_eltptr (convolution, 0, 0) = 1;
  } else {
    const double alpha = log (2.0) / SQR (gradius);
    for (y = -m; y <= m; y++) {
      for (x = -m; x <= m; x++) {
        *c_mat_eltptr (convolution, x, y) =
                       exp (-alpha * (SQR (x) + SQR (y)));
      }
    }
  }
  return (1);
}

static double circle_integral (const double x, const double radius) {
  /* Return the integral of sqrt(radius^2 - z^2) for z = 0 to x */
  double retval;
  if (ABS (radius) <= 0.001) {
    /* Perhaps some epsilon must be added here? */
    retval = 0.0;
  } else {
    const double sin = x / radius;
    const double sq_diff = SQR (radius) - SQR (x);
    /* From a mathematical point of view the following is redundant.
       Numerically they are not equivalent! */
    if ((sq_diff < 0.0) || (sin < -1.0) || (sin > 1.0)) {
      if (sin < 0) {
        retval = (-0.25 * SQR (radius) * G_PI);
      } else {
        retval = (0.25 * SQR (radius) * G_PI);
      }
    } else {
      retval = (0.5 * x * sqrt (sq_diff) + 0.5 * SQR (radius) * asin (sin));
    }
    if (retval < 0.00001)
      retval = 0.0;
  }
  return (retval);
}

double circle_intensity (const gint x, const gint y, const double radius) {

  if (ABS (radius) <= 0.001) {
    return (((x == 0) && (y == 0)) ? 1 : 0);
  } else {
    register double xlo = ABS (x) - 0.5, xhi = ABS (x) + 0.5,
        ylo = ABS (y) - 0.5, yhi = ABS (y) + 0.5;
    register double symmetry_factor = 1, xc1, xc2, retval;

    if (xlo < 0) {
      xlo = 0;
      symmetry_factor *= 2;
    }
    if (ylo < 0) {
      ylo = 0;
      symmetry_factor *= 2;
    }
    if (SQR (xlo) + SQR (yhi) > SQR (radius)) {
      xc1 = xlo;
    } else if (SQR (xhi) + SQR (yhi) > SQR (radius)) {
      xc1 = sqrt (SQR (radius) - SQR (yhi));
    } else {
      xc1 = xhi;
    }
    if (SQR (xlo) + SQR (ylo) > SQR (radius)) {
      xc2 = xlo;
    } else if (SQR (xhi) + SQR (ylo) > SQR (radius)) {
      xc2 = sqrt (SQR (radius) - SQR (ylo));
    } else {
      xc2 = xhi;
    }
    retval = (((yhi - ylo) * (xc1 - xlo) +
              circle_integral (xc2, radius) - circle_integral (xc1, radius) -
              (xc2 - xc1) * ylo) * symmetry_factor / (G_PI * SQR (radius)));
    if (retval < 0.00001)
      retval = 0.0;
    return (retval);
  }
}

int make_circle_convolution (const double radius, CMat * convolution, const gint m) {
  if (!(fill_matrix (convolution, m, circle_intensity, radius)))
    return (0);
  return (1);
}
