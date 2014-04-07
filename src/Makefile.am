## Process this file with automake to produce Makefile.in

INCLUDES=@GTK_CFLAGS@ @GIMP_CFLAGS@ @LAPACK_INCLUDES@

bin_PROGRAMS = refocus
noinst_PROGRAMS = test-matrix

ldadd_atlas=@LAPACK_LIBS@

noinst_HEADERS=bdclosure.h refocus.h matrix.h tilebuf.h conv.h fwlapack.h gimppreview.h prevman.h util.h
EXTRA_DIST = dummy-plugin.c

refocus_SOURCES = refocus.c conv.c gimppreview.c prevman.c util.c matrix.c tilebuf.c bdclosure.c fwlapack.c
refocus_LDADD = -lm ${ldadd_atlas} @GIMP_LIBS@
if HAVE_ATLAS
else
refocus_DEPENDENCIES=@BUNDLED_LAPACK_LIB_DIR@/liblapack.a
endif

test_matrix_SOURCES = matrix.c test-matrix.c fwlapack.c
test_matrix_LDADD = -lm ${ldadd_atlas} @GLIB_LIBS@
if HAVE_ATLAS
else
test_matrix_DEPENDENCIES=@BUNDLED_LAPACK_LIB_DIR@/liblapack.a
endif

install-bin: refocus
	@if test "x@GIMPTOOL@" = "x" ; then \
		echo "Warning: configure did not find gimptool" ; \
		echo "install src/refocus manually in your plug-in directory" ; \
		echo "e.g.: cp src/refocus ~/.gimp-1.2/plug-ins" ; \
	else \
		@GIMPTOOL@ --install-bin refocus ; \
	fi

install-admin-bin: refocus
	@if test "x@GIMPTOOL@" = "x" ; then \
		echo "Warning: configure did not find gimptool" ; \
		echo "install src/refocus manually in your plug-in directory" ; \
		echo "e.g.: install src/refocus /usr/local/lib/gimp-1.2/plug-ins" ; \
	else \
		@GIMPTOOL@ --install-admin-bin refocus ; \
	fi

@BUNDLED_LAPACK_LIB_DIR@/liblapack.a:
	cd @BUNDLED_LAPACK_LIB_DIR@; gmake

%.s: %.c
	${CC} ${CFLAGS} -S $< -o $@


