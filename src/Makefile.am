## Process this file with automake to produce Makefile.in

AM_CFLAGS = ${CFLAGS} ${CPPFLAGS} ${GTK_CFLAGS} ${GIMP_CFLAGS} ${LAPACK_INCLUDES}

AM_CPPFLAGS = -I${builddir} -I${srcdir} ${GIMP_CFLAGS} -I${includedir}

noinst_HEADERS = bdclosure.h refocus.h matrix.h tilebuf.h conv.h \
		 fwlapack.h prevman.h util.h
EXTRA_DIST = dummy-plugin.c
nodist_EXTRA_DATA = .dep .lib

bin_PROGRAMS = refocus
bindir = $(GIMP_BINDIR)/plug-ins
noinst_PROGRAMS = test-matrix

ldadd_atlas=@LAPACK_LIBS@

refocus_SOURCES = refocus.c conv.c prevman.c util.c \
		  matrix.c tilebuf.c bdclosure.c fwlapack.c
refocus_LDADD = -lm ${ldadd_atlas} ${LIBS} ${GIMP_LIBS} ${GTK_LIBS}
if HAVE_ATLAS
else
refocus_DEPENDENCIES=@BUNDLED_LAPACK_LIB_DIR@/liblapack.a
endif

test_matrix_SOURCES = matrix.c test-matrix.c fwlapack.c
test_matrix_LDADD = -lm ${ldadd_atlas} ${LIBS} ${GLIB_LIBS} ${GTK_LIBS}
if HAVE_ATLAS
else
test_matrix_DEPENDENCIES=@BUNDLED_LAPACK_LIB_DIR@/liblapack.a
endif

@BUNDLED_LAPACK_LIB_DIR@/liblapack.a:
	cd @BUNDLED_LAPACK_LIB_DIR@; gmake
