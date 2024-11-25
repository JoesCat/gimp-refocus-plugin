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

#if __has_include("refocus-config.h")
#include "refocus-config.h"
#else
#define PLUGIN_VERSION "refocus is local"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <libgimp/gimp.h>
#include <libgimp/gimpui.h>
#include <sys/types.h>
#include <signal.h>
#include <unistd.h>
#include "refocus.h"
#include "prevman.h"
#include "matrix.h"
#include "conv.h"
#include "bdclosure.h"
#include "util.h"

/* No i18n for now */
#define _(x)	x
#define d_(x)	x
#define N_(x)	x

#define SCALE_WIDTH	  150
#define ENTRY_WIDTH	  4

#define MATRIX_MIN	  1.0
#define MATRIX_MAX	  25.0
#define MATRIX_STEP	  1.0
#define RADIUS_MIN	  0.0
#define RADIUS_MAX	  25.0
#define RADIUS_STEP	  0.1
#define ALPHA_MIN	  0.0
#define ALPHA_MAX	  25.0
#define ALPHA_STEP	  0.1
#define GAMMA_MIN	  0.0
#define GAMMA_MAX	  1.0
#define GAMMA_STEP	  0.05
#define NOISE_FACTOR_MIN  0.0
#define NOISE_FACTOR_MAX  1.0
#define NOISE_FACTOR_STEP 0.01

static const char PROCEDURE_NAME[] = "plug-in-refocus";
static const char SHORT_DESCRIPTION[] = "Refocus using FIR Wiener Deconvolution";
static const char LONG_DESCRIPTION[] =
  "The Refocus Gimp plug-in can be used to sharpen images by deblurring focus-blur or gaussian-blur.\n\n"
  "Frequently, when processing images, e.g. scanning photo's or slides, the images can become slightly blurred. "
  "The blurred impression of these images is due to the fact that image pixels are averaged with their neighbors. "
  "Blurred images don't have sharp boundaries and look as though they have been taken with an unfocussed camera.\n\n"
  "This plug-in attempts to 'refocus' the image. "
  "In some cases Refocus can produce better results than plug-ins such as sharpen or unsharp mask.";

typedef struct {
  gint    mat_width;
  gdouble radius;
  gdouble alpha;
  gdouble gamma;
  gdouble noise_factor;
} config;

typedef struct {
  config       *config;
  GimpDrawable *drawable;
  CMat         *matrix;
  gdouble       mat_width_d;
} ref_params;

/* Declare local functions. */
static void query (void);
static void run (const gchar *name,
                 gint nparams,
                 const GimpParam *param,
                 gint *nreturn_vals, GimpParam **return_vals);
static gboolean refocus_dialog (ref_params *params);
static gint doit (ref_params *params);
static int  check_config_values (config *my_config, gboolean reset);
static void refocus_help (const gchar *help_id, gpointer help_data);

GimpPlugInInfo PLUG_IN_INFO = {
  NULL,  /* init_proc  */
  NULL,  /* quit_proc  */
  query, /* query_proc */
  run,   /* run_proc   */
};

gboolean matrix_needs_updating = TRUE;

static const config default_config = {
  5,     /* matrix_width */
  1,     /* radius       */
  0.0,   /* gaussian     */
  0.5,   /* correlation  */
  0.01   /* noise factor */
};

struct {
  GtkWidget *preview;
  GtkWidget *mat_width_entry;
  GtkWidget *radius_entry;
  GtkWidget *alpha_entry;
  GtkWidget *gamma_entry;
  GtkWidget *noise_entry;
  GtkWidget *ok_button;
  GtkWidget *update_preview_button;
} my_widgets;

MAIN ()

static void query (void) {
  static GimpParamDef args[] = {
    {GIMP_PDB_INT32, "run_mode", "Interactive, non-interactive"},
    {GIMP_PDB_IMAGE, "image", "Input image (unused)"},
    {GIMP_PDB_DRAWABLE, "drawable", "Input drawable to modify"},
    {GIMP_PDB_INT32, "mat_size", "Size of matrix"},
    {GIMP_PDB_FLOAT, "radius", "Circle radius"},
    {GIMP_PDB_FLOAT, "gauss", "Parameter for Gaussian convolution"},
    {GIMP_PDB_FLOAT, "correlation", "Correlation"},
    {GIMP_PDB_FLOAT, "noise", "Noise to Signal ratio"},
  };
  static GimpParamDef *return_vals = NULL;
  static gint nreturn_vals = 0;

  gimp_install_procedure (PROCEDURE_NAME,
                          SHORT_DESCRIPTION,
                          LONG_DESCRIPTION,
  /* copyright author  */ "Ernst Lippe, 1999, <ernstl@planet.nl>",
  /* copyright license */ "GPL3+",
  /* build date        */ "2024",
                          "Refocus",
                          "RGB*, GRAY*",
                          GIMP_PLUGIN,
                          G_N_ELEMENTS(args), nreturn_vals, args, return_vals);

  gimp_plugin_menu_register (PROCEDURE_NAME, "<Image>/Filters/Enhance");
}

static void run (const gchar *name, gint n_params, const GimpParam *param,
                 gint *nreturn_vals, GimpParam **return_vals) {
  static GimpParam  values[1];
  config            my_config;
  GimpDrawable     *drawable;
  ref_params        my_params;
  GimpRunMode       run_mode;
  GimpPDBStatusType status;

  /* sleep(30); */ /* Sleep so the debugger can attach to this process */
  *nreturn_vals  = 1;
  *return_vals   = values;
  values[0].type = GIMP_PDB_STATUS;
  status         = GIMP_PDB_SUCCESS;

  if (param[0].type != GIMP_PDB_INT32 || strcmp(name, PROCEDURE_NAME) != 0) {
    values[0].data.d_status = GIMP_PDB_CALLING_ERROR;
    return;
  }

  run_mode = param[0].data.d_int32;

  /* Get the specified drawable */
  my_params.drawable = drawable = gimp_drawable_get (param[2].data.d_drawable);

  /* Make sure that the drawable is gray or RGB color */
  if (!gimp_drawable_is_rgb (drawable->drawable_id) &&
      !gimp_drawable_is_gray (drawable->drawable_id)) {
    g_message("Image is not RGB, RGBA, gray or grayA!");
    status = GIMP_PDB_EXECUTION_ERROR;;
    return;
  }

  /* Although the convolution should work fine with a tiny cache,
     it is made bigger to improve scrolling speed */
  gimp_tile_cache_ntiles (20);

  my_config = default_config;
  my_params.config = &my_config;
  my_params.matrix = NULL;

  switch (run_mode) {
    case GIMP_RUN_WITH_LAST_VALS:
      gimp_get_data (PROCEDURE_NAME, &my_config);
      gimp_ui_init (PROCEDURE_NAME, TRUE);
      break;

    case GIMP_RUN_NONINTERACTIVE:
      if (n_params != 8) {
        g_message("Incorrect number of input parameters!");
        status = GIMP_PDB_EXECUTION_ERROR;;
      } else {
        my_config.mat_width = param[3].data.d_int32;
        my_config.radius = param[4].data.d_float;
        my_config.alpha = param[5].data.d_float;
        my_config.gamma = param[6].data.d_float;
        my_config.noise_factor = param[7].data.d_float;
        if (check_config_values (&my_config, FALSE)) {
          status = GIMP_PDB_EXECUTION_ERROR;;
        }
      }
      break;

    case GIMP_RUN_INTERACTIVE:
      gimp_get_data (PROCEDURE_NAME, &my_config);
      check_config_values (&my_config, TRUE);
      gimp_ui_init (PROCEDURE_NAME, TRUE);
      /* Oh boy. We get to do a dialog box, because we can't really expect
       * the user to set this up with the right values using gdb.
       */
      if (!refocus_dialog (&my_params)) {
        /* The dialog was closed, or something similarly evil happened. */
        status = GIMP_PDB_EXECUTION_ERROR;
      }
      break;

    default:
      break;
  }

  if (status == GIMP_PDB_SUCCESS) {
    if (doit (&my_params))
      status = GIMP_PDB_EXECUTION_ERROR;
    else if (run_mode == GIMP_RUN_INTERACTIVE)
      gimp_set_data (PROCEDURE_NAME, &my_config, sizeof (my_config));

    if (run_mode != GIMP_RUN_NONINTERACTIVE)
      gimp_displays_flush ();

    gimp_drawable_detach (drawable);
  }

  values[0].data.d_status = status;
}

static int update_matrix (ref_params *params) {
  /* Recompute matrix if needed */
  CMat circle, gaussian, convolution;

  if (params->matrix)
    finish_c_mat (params->matrix);

  if (!(make_gaussian_convolution (params->config->alpha, &gaussian,
                             params->config->mat_width)))
    goto update_matrix_mem_err_gaussian;

  if (!(make_circle_convolution (params->config->radius, &circle,
                           params->config->mat_width)))
    goto update_matrix_mem_err_circle;
#ifdef RF_DEBUG
  fprintf (stderr, "mat_width=%d alpha=%f radius=%f gamma=%f noise=%f\n",
           params->config->mat_width, params->config->alpha,
           params->config->radius, params->config->gamma,
           params->config->noise_factor);
  fprintf (stderr, "Gauss\n");
  print_c_mat (stderr, &gaussian);
  fprintf (stderr, "Circle\n");
  print_c_mat (stderr, &circle);
#endif

  /* TODO: must use normal convolution */
  if (!(init_c_mat (&convolution, params->config->mat_width)))
    goto update_matrix_mem_err_convolution;
  convolve_star_mat (&convolution, &gaussian, &circle);
  if (!(params->matrix = compute_g_matrix (&convolution, params->config->mat_width,
                                      params->config->gamma,
                                      params->config->noise_factor, 0.0, TRUE)))
    goto update_matrix_mem_err_matrix;

  finish_c_mat (&convolution);
  finish_c_mat (&gaussian);
  finish_c_mat (&circle);
  return (1);

update_matrix_mem_err_matrix:
  finish_c_mat (&convolution);
update_matrix_mem_err_convolution:
  finish_c_mat (&circle);
update_matrix_mem_err_circle:
  finish_c_mat (&gaussian);
update_matrix_mem_err_gaussian:
  return (0);
}

/* Checks that the configuration is valid for the image type */
static int check_config_values (config *my_config, gboolean reset) {
  if (reset) {
    /* reset value inside limits if out of bounds */
    if (my_config->mat_width < MATRIX_MIN) my_config->mat_width = MATRIX_MIN;
    if (my_config->mat_width > MATRIX_MAX) my_config->mat_width = MATRIX_MAX;
    if (my_config->radius < RADIUS_MIN) my_config->radius = RADIUS_MIN;
    if (my_config->radius > RADIUS_MAX) my_config->radius = RADIUS_MAX;
    if (my_config->alpha < ALPHA_MIN) my_config->alpha = ALPHA_MIN;
    if (my_config->alpha > ALPHA_MAX) my_config->alpha = ALPHA_MAX;
    if (my_config->gamma < GAMMA_MIN) my_config->gamma = GAMMA_MIN;
    if (my_config->gamma > GAMMA_MAX) my_config->gamma = GAMMA_MAX;
    if (my_config->noise_factor < NOISE_FACTOR_MIN)
      my_config->noise_factor = NOISE_FACTOR_MIN;
    if (my_config->noise_factor > NOISE_FACTOR_MAX)
      my_config->noise_factor = NOISE_FACTOR_MAX;
  } else {
    /* report error if parameter is out of bounds */
    if (my_config->mat_width < MATRIX_MIN ||
        my_config->mat_width > MATRIX_MAX ||
        my_config->radius < RADIUS_MIN ||
        my_config->radius > RADIUS_MAX ||
        my_config->alpha < ALPHA_MIN ||
        my_config->alpha > ALPHA_MAX ||
        my_config->gamma < GAMMA_MIN ||
        my_config->gamma > GAMMA_MAX ||
        my_config->noise_factor < NOISE_FACTOR_MIN ||
        my_config->noise_factor > NOISE_FACTOR_MAX) {
      g_message("Bad input parameter!");
      return (1);
    }
  }
#ifdef RF_DEBUG
  fprintf (stderr,
           "my_config{ matrix=%d, radius=%g, alpha=%g, gamma=%g, noise=%g }\n",
           my_config->mat_width, my_config->radius, my_config->alpha,
           my_config->gamma, my_config->noise_factor);
#endif
  return (0);
}

gboolean preview_progress_update_fun (const gpointer data, gdouble arg) {
  return TRUE;
}

static void preview_callback (GtkWidget * widget, ref_params *params) {
  GimpDrawablePreview *preview;
  GimpPreview *ptr;
  gint32       preview_ID;
  gint         image_x, image_y, im_width, im_height, bppImg;
  const Babl  *format;
  TileSource   source;
  TileSink     sink;
  gint         row;
  guchar      *buf;
  gboolean event_is_current = TRUE;
  BDClosure update_progress_closure;

  preview = GIMP_DRAWABLE_PREVIEW (widget);
  ptr = GIMP_PREVIEW (preview);
  gimp_preview_get_position (ptr, &image_x, &image_y);
  gimp_preview_get_size (ptr, &im_width, &im_height);

  preview_ID = gimp_drawable_preview_get_drawable_id (preview);
  format = gimp_drawable_get_format (preview_ID);
  bppImg = babl_format_get_bytes_per_pixel(format);

  params->config->mat_width = round(params->mat_width_d);
  if (!(update_matrix (params))) return;
  tile_source_init_from_drawable (&source, params->drawable, image_x, image_y,
                                  im_width, im_height);
  tile_sink_init_for_preview (&sink, params->drawable, image_x, image_y,
                              im_width, im_height);
  bd_closure_init (&update_progress_closure,
                   preview_progress_update_fun,
                   GINT_TO_POINTER (preview_ID));
  convolve_image (&source, &sink, image_x, image_y,
                  im_width, im_height,
                  TB_BOUNDARY_MIRROR,
                  params->matrix, 2 * params->config->mat_width + 1,
                  &update_progress_closure);
  buf = g_new (guchar, im_height * im_width * bppImg);
  for (row = 0; event_is_current && (row < im_height); row++) {
    tile_sink_get_row (&sink, &buf[row * im_width * bppImg],
                       image_x, image_y + row, im_width);
  };
  gimp_preview_draw_buffer (ptr, buf, im_width * bppImg);
  g_free (buf);
  tile_sink_free_buffers (&sink);
}

static gint refocus_dialog (ref_params *params) {
  GtkWidget *dialog;
  GtkWidget *main_vbox;
  GtkWidget *table;
  GtkWidget *preview;
  GtkObject *entry;
  gboolean  run;

  params->mat_width_d = params->config->mat_width;

  dialog = gimp_dialog_new (_("Refocus"), "refocus", NULL, 0,
                            refocus_help, PROCEDURE_NAME,

                            GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                            GTK_STOCK_OK,     GTK_RESPONSE_OK,

                            NULL);

  main_vbox = gtk_vbox_new (FALSE, 12);
  gtk_container_set_border_width (GTK_CONTAINER (main_vbox), 12);
  gtk_container_add (GTK_CONTAINER (GTK_DIALOG (dialog)->vbox), main_vbox);
  gtk_widget_show (main_vbox);

  preview = gimp_drawable_preview_new_from_drawable_id (params->drawable->drawable_id);
  gtk_box_pack_start (GTK_BOX (main_vbox), preview, TRUE, TRUE, 0);
  gtk_widget_show (preview);
  g_signal_connect (preview, "invalidated",
                    G_CALLBACK (preview_callback),
                    params);

  table = gtk_table_new (2, 5, FALSE);
  gtk_table_set_col_spacings (GTK_TABLE (table), 6);
  gtk_table_set_row_spacings (GTK_TABLE (table), 6);
  gtk_box_pack_start (GTK_BOX (main_vbox), table, FALSE, FALSE, 0);
  gtk_widget_show (table);

  entry = gimp_scale_entry_new (GTK_TABLE (table), 0, 0,
                                _("Matrix Size"), SCALE_WIDTH, ENTRY_WIDTH,
                                params->mat_width_d,
                                MATRIX_MIN, MATRIX_MAX, 1.0, 5.0, 0,
                                TRUE, 0, 0, NULL, NULL);
  g_signal_connect (entry, "value_changed",
                    G_CALLBACK (gimp_double_adjustment_update),
                    &(params->mat_width_d));
  g_signal_connect_swapped (entry, "value_changed",
                            G_CALLBACK (gimp_preview_invalidate), preview);

  entry = gimp_scale_entry_new (GTK_TABLE (table), 0, 1,
                                _("Radius"), SCALE_WIDTH, ENTRY_WIDTH,
                                params->config->radius,
                                RADIUS_MIN, RADIUS_MAX, 0.1, 0.5, 1,
                                TRUE, 0, 0, NULL, NULL);
  g_signal_connect (entry, "value_changed",
                    G_CALLBACK (gimp_double_adjustment_update),
                    &(params->config->radius));
  g_signal_connect_swapped (entry, "value_changed",
                            G_CALLBACK (gimp_preview_invalidate), preview);

  entry = gimp_scale_entry_new (GTK_TABLE (table), 0, 2,
                                _("Gauss"), SCALE_WIDTH, ENTRY_WIDTH,
                                params->config->alpha,
                                ALPHA_MIN, ALPHA_MAX, 0.1, 0.5, 1,
                                TRUE, 0, 0, NULL, NULL);
  g_signal_connect (entry, "value_changed",
                    G_CALLBACK (gimp_double_adjustment_update),
                    &(params->config->alpha));
  g_signal_connect_swapped (entry, "value_changed",
                            G_CALLBACK (gimp_preview_invalidate), preview);

  entry = gimp_scale_entry_new (GTK_TABLE (table), 0, 3,
                                _("Correlation"), SCALE_WIDTH, ENTRY_WIDTH,
                                params->config->gamma,
                                GAMMA_MIN, GAMMA_MAX, 0.05, 0.1, 2,
                                TRUE, 0, 0, NULL, NULL);
  g_signal_connect (entry, "value_changed",
                    G_CALLBACK (gimp_double_adjustment_update),
                    &(params->config->gamma));
  g_signal_connect_swapped (entry, "value_changed",
                            G_CALLBACK (gimp_preview_invalidate), preview);

  entry = gimp_scale_entry_new (GTK_TABLE (table), 0, 4,
                                _("Noise"), SCALE_WIDTH, ENTRY_WIDTH,
                                params->config->noise_factor,
                                NOISE_FACTOR_MIN, NOISE_FACTOR_MAX, 0.01, 0.1, 2,
                                TRUE, 0, 0, NULL, NULL);
  g_signal_connect (entry, "value_changed",
                    G_CALLBACK (gimp_double_adjustment_update),
                    &(params->config->noise_factor));
  g_signal_connect_swapped (entry, "value_changed",
                            G_CALLBACK (gimp_preview_invalidate), preview);
  gtk_widget_show (table);

  gtk_widget_show (dialog);
  run = (gimp_dialog_run (GIMP_DIALOG (dialog)) == GTK_RESPONSE_OK);
  gtk_widget_destroy (dialog);
  g_free (params->matrix);
  return run;
}

static gboolean
gimp_progress_update_fun (gpointer data, gdouble fraction)
{
  return (gimp_progress_update (fraction));
}

static gint doit (ref_params *params) {
  TileSource source;
  TileSink sink;
  gint x, y, width, height;
  BDClosure update_progress_closure;

#ifdef RF_DEBUG
  fprintf (stderr, "doing doit()\n");
#endif

  if (!(gimp_drawable_mask_intersect(params->drawable->drawable_id, &x, &y, &width, &height)))
    return (-1);

  bd_closure_init (&update_progress_closure, gimp_progress_update_fun, NULL);
  gimp_progress_init ("Computing matrix");
  if (!(update_matrix (params)))
    return (-1);
  gimp_progress_init ("Applying convolution");
  tile_source_init_from_drawable (&source, params->drawable, x, y, width, height);
  tile_sink_init_from_drawable (&sink, params->drawable, x, y, width, height);
  convolve_image (&source, &sink, x, y, width, height,
                  TB_BOUNDARY_MIRROR,
                  params->matrix, 2 * params->config->mat_width + 1,
                  &update_progress_closure);
  gimp_drawable_flush (params->drawable);
  gimp_drawable_merge_shadow (params->drawable->drawable_id, TRUE);
  gimp_drawable_update (params->drawable->drawable_id, x, y, width, height);
  g_free (params->matrix);
  return (0);
}

static void refocus_help (const gchar *help_id, gpointer help_data) {
  gimp_message (_(LONG_DESCRIPTION));
}
