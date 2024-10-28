dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.69)

dnl Name your plug-in here
m4_define([plugin_name], [refocus])

dnl These two define the plug-in version number
m4_define([plugin_major_version], [0])
m4_define([plugin_minor_version], [9])

m4_define([plugin_version], [plugin_major_version.plugin_minor_version])
m4_define([plugin_package_name], [gimp-refocus])
m4_define([plugin_package_home], [https://github.com/JoesCat/gimp-refocus])
m4_define([plugin_package_email], [https://github.com/JoesCat/gimp-refocus/issues])

AC_INIT([plugin_name], [plugin_version], [plugin_package_email],
	[plugin_package_name], [plugin_package_home])

AC_DEFINE(PLUGIN_NAME, ["plugin_name"], [Plug-In name])
AC_DEFINE(PLUGIN_EMAIL, [plugin_package_email], [Plug-In package email])
AC_DEFINE(PACKAGE_NAME, ["plugin_package_name"], [Plug-In Package name])
AC_DEFINE(PLUGIN_VERSION, ["plugin_version"], [Plug-In version])

AC_DEFINE(PLUGIN_MAJOR_VERSION, plugin_major_version, [Plug-In major version])
AC_DEFINE(PLUGIN_MINOR_VERSION, plugin_minor_version, [Plug-In minor version])

AC_CONFIG_SRCDIR([src/refocus.c])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_TARGET
AC_CANONICAL_HOST
AC_CANONICAL_BUILD
AC_USE_SYSTEM_EXTENSIONS
AM_INIT_AUTOMAKE([foreign -Wall])
#--------------------------------------------------------------------------
# automake 1.12 needs AM_PROG_AR but automake < 1.11.2 doesn't recognize it
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

LT_INIT
AC_SUBST([LIBTOOL_DEPS])

#--------------------------------------------------------------------------
# Checks for programs.
AC_PROG_CC
AC_PROG_SED
AC_PROG_LN_S
AC_PROG_MKDIR_P
AC_PATH_PROG([CHMOD],[chmod],[:])
AC_PATH_PROG([STRIP],[strip],[:])
AC_PATH_PROG([GIMP],[gimp-2.10],[:])
AC_PATH_PROG([GIMPTOOL],[gimptool-2.0],[:])
AM_CONFIG_HEADER(refocus-config.h)
AC_PROG_INSTALL
AC_PROG_MAKE_SET

AM_MAINTAINER_MODE

dnl Use -Wall if we have gcc.
changequote(,)dnl
if test "x$GCC" = "xyes"; then
  case " $CFLAGS " in
  *[\ \	]-Wall[\ \	]*) ;;
  *) CFLAGS="$CFLAGS -Wall" ;;
  esac
fi
changequote([,])dnl

#--------------------------------------------------------------------------
# Search for math.h include and math lib (some OSes have -lm built-in).
have_libm=maybe
AC_CHECK_HEADERS([math.h],
  AC_SEARCH_LIBS([cos],[m],[have_libm=yes]))
if test x"${have_libm}" != xyes; then
  AC_MSG_FAILURE([ERROR: Please install the Math library and math.h],[1])
fi
#Check if erf is defined in the math library
have_erf=no
AC_SEARCH_LIBS([erf],[m],[have_erf=yes],
  AC_CHECK_FUNC([erf],[have_erf=yes]))
if test x"${have_erf}" = xyes; then
  AC_DEFINE(HAVE_ERF,["have_erf"],[PMath library has erf])
fi

BUNDLED_LAPACK_LIB_DIR='${top_srcdir}/lib'
AC_SUBST(BUNDLED_LAPACK_LIB_DIR)

have_atlas=no
AC_ARG_WITH([lapack-libs],AS_HELP_STRING([--with-lapack-libs=LIBS]))
if test "x$with_lapack_libs" = "x" ; then
  LAPACK_LIBS='-L${top_srcdir}/lib -latlas -lcblas -llapack'
else
  LAPACK_LIBS="$with_lapack_libs"
  have_atlas=yes
fi
AC_SUBST(LAPACK_LIBS)
AM_CONDITIONAL(HAVE_ATLAS, test x${have_atlas} = xyes)

AC_ARG_WITH(lapack-includes, [  --with-lapack-includes=INCLUDES ])
if test "x$with_lapack_includes" = "x" ; then
  LAPACK_INCLUDES='-I ${top_srcdir}/lib'
else
  LAPACK_INCLUDES="$with_lapack_includes"
fi
AC_SUBST(LAPACK_INCLUDES)

AC_MSG_CHECKING([for some Win32 platform])
case "$host" in
  *-*-mingw*|*-*-cygwin*)
    platform_win32=yes
    ;;
  *)
    platform_win32=no
    ;;
esac
AC_MSG_RESULT([$platform_win32])
AM_CONDITIONAL(PLATFORM_WIN32, test "$platform_win32" = "yes")

AC_MSG_CHECKING([for native Win32])
case "$host" in
  *-*-mingw*)
    os_win32=yes
    PATHSEP=';'
    ;;
  *)
    os_win32=no
    PATHSEP=':'
    ;;
esac
AC_MSG_RESULT([$os_win32])
AC_SUBST(PATHSEP)
AM_CONDITIONAL(OS_WIN32, test "$os_win32" = "yes")

if test "$os_win32" = "yes"; then
  AC_CHECK_PROG(ms_librarian, lib.exe, yes, no)
fi
AM_CONDITIONAL(MS_LIB_AVAILABLE, test x$ms_librarian = xyes)

BUNDLED_LAPACK_LIB_DIR='${top_srcdir}/lib'
AC_SUBST(BUNDLED_LAPACK_LIB_DIR)

have_atlas=no
AC_ARG_WITH(lapack-libs, [  --with-lapack-libs=LIBS ])
if test "x$with_lapack_libs" = "x" ; then
  LAPACK_LIBS='-L${top_srcdir}/lib -latlas -lcblas -llapack'
else
  LAPACK_LIBS="$with_lapack_libs"
  have_atlas=yes
fi
AC_SUBST(LAPACK_LIBS)
AM_CONDITIONAL(HAVE_ATLAS, test x${have_atlas} = xyes)

AC_ARG_WITH(lapack-includes, [  --with-lapack-includes=INCLUDES ])
if test "x$with_lapack_includes" = "x" ; then
  LAPACK_INCLUDES='-I ${top_srcdir}/lib'
else
  LAPACK_INCLUDES="$with_lapack_includes"
fi
AC_SUBST(LAPACK_INCLUDES)

# Check for libgimp/gimp.h include file and libgimp library.
GIMP_CFLAGS=
GIMP_LIBS=
have_libgimp=no
PKG_CHECK_MODULES([GIMP],[gimp-2.0 >= 2.10.0 gimpui-2.0 >= 2.10.0],[have_libgimp=yes])
if test x"${have_libgimp}" != xyes; then
    AC_MSG_FAILURE([ERROR: Please install the developer version of libgimp2.],[1])
fi
AC_SUBST(GIMP_CFLAGS)
AC_SUBST(GIMP_LIBS)
AM_CONDITIONAL([HAVEGIMPTOOL],[test "${GIMPTOOL}"x != x])

GTK_CFLAGS=
GTK_LIBS=
have_libgtk=no
PKG_CHECK_MODULES([GTK],[gtk+-2.0 >= 2.3.0],[have_libgtk=yes])
if test x"${have_libgtk}" != xyes; then
    AC_MSG_FAILURE([ERROR: Please install the developer version of libgtk+2.],[1])
fi
AC_SUBST(GTK_CFLAGS)
AC_SUBST(GTK_LIBS)

CPPFLAGS="$CPPFLAGS -DG_DISABLE_DEPRECATED -DGDK_DISABLE_DEPRECATED -DGTK_DISABLE_DEPRECATED"
CPPFLAGS="${CPPFLAGS} AS_ESCAPE([-I${top_builddir}]) AS_ESCAPE([-I${top_srcdir}])"
AC_SUBST([CPPFLAGS],["${CPPFLAGS}"])

#--------------------------------------------------------------------------
dnl Ensure MSVC-compatible struct packing convention is used when
dnl compiling for Win32 with gcc. GTK+ uses this convention, so we must, too.
dnl What flag to depends on gcc version: gcc3 uses "-mms-bitfields", while
dnl gcc2 uses "-fnative-struct".
if test x"$os_win32" = xyes; then
  if test x"$GCC" = xyes; then
    msnative_struct=''
    AC_MSG_CHECKING([how to get MSVC-compatible struct packing])
    if test -z "$ac_cv_prog_CC"; then
      our_gcc="$CC"
    else
      our_gcc="$ac_cv_prog_CC"
    fi
    case `$our_gcc --version | sed -e 's,\..*,.,' -e q` in
      2.)
	if $our_gcc -v --help 2>/dev/null | grep fnative-struct >/dev/null; then
	  msnative_struct='-fnative-struct'
	fi
	;;
      *)
	if $our_gcc -v --help 2>/dev/null | grep ms-bitfields >/dev/null; then
	  msnative_struct='-mms-bitfields'
	fi
	;;
    esac
    if test x"$msnative_struct" = x ; then
      AC_MSG_RESULT([no way])
      AC_MSG_WARN([produced libraries will be incompatible with prebuilt GTK+ DLLs])
    else
      CFLAGS="$CFLAGS $msnative_struct"
      AC_MSG_RESULT([${msnative_struct}])
    fi
  fi
fi

# This is a check for gtk-doc which you can insert into your configure.in.
# You shouldn't need to change it at all.
##################################################
# Check for gtk-doc.
##################################################

AC_ARG_WITH(html-dir, [  --with-html-dir=PATH path to installed docs ])

if test "x$with_html_dir" = "x" ; then
  HTML_DIR='${datadir}/gtk-doc/html'
else
  HTML_DIR=$with_html_dir
fi

AC_SUBST(HTML_DIR)

AC_CHECK_PROG(GTKDOC, gtkdoc-mkdb, true, false)

gtk_doc_min_version=0.6
if $GTKDOC ; then 
    gtk_doc_version=`gtkdoc-mkdb --version`
    AC_MSG_CHECKING([gtk-doc version ($gtk_doc_version) >= $gtk_doc_min_version])
    if perl <<EOF ; then
      exit (("$gtk_doc_version" =~ /^[[0-9]]+\.[[0-9]]+$/) &&
            ("$gtk_doc_version" >= "$gtk_doc_min_version") ? 0 : 1);
EOF
      AC_MSG_RESULT(yes)
   else
      AC_MSG_RESULT(no)
      GTKDOC=false
   fi
fi

dnl Let people enable the gtk-doc stuff.
AC_ARG_ENABLE(gtk-doc, [  --enable-gtk-doc        Use gtk-doc to build documentation [default=no]], enable_gtk_doc="$enableval", enable_gtk_doc=no)

if test x$enable_gtk_doc = xyes ; then
  if test x$GTKDOC = xtrue ; then
    enable_gtk_doc=yes
  else
    enable_gtk_doc=no 
  fi
fi

AM_CONDITIONAL(ENABLE_GTK_DOC, test x$enable_gtk_doc = xyes)
# End of gtk-doc stuff.

dnl DEBUG
enable_debug=no
AC_ARG_ENABLE([debug],AS_HELP_STRING([--enable-debug],[Enable debugging @<:@default=no@:>]),[enable_debug=yes],[enable_debug=no])
if test x$enable_debug = xyes ; then
  AC_DEFINE(RF_DEBUG,[1],[RF Debug])
  AC_DEFINE(PREVIEW_DEBUG,[1],[Preview Debug])
fi

#--------------------------------------------------------------------------
# Pass variables to MAKEFILE.AM
AC_SUBST([PLUGIN_MAJOR_VERSION],[plugin_major_version])
AC_SUBST([PLUGIN_MINOR_VERSION],[plugin_minor_version])
AC_SUBST([PLUGIN_VERSION],[plugin_version])
AC_SUBST([PLUGIN_PACKAGE_NAME],[plugin_package_email])
AC_SUBST([PLUGIN_EMAIL],[plugin_package_email])
AC_SUBST([HOST],["$host"])

# Put ifndef wrapper on refocus-config.h so we don't call it repeatedly.
AH_TOP([#ifndef REFOCUS_CONFIG_H
#define REFOCUS_CONFIG_H 1])
AH_BOTTOM([

#endif])

#--------------------------------------------------------------------------
AC_CONFIG_FILES([
Makefile
src/Makefile
lib/Makefile
doc/Makefile
gtk-doc/Makefile
])

AC_OUTPUT
