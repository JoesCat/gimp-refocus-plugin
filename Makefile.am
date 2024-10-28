## Process this file with automake to produce Makefile.in

SUBDIRS = src lib doc
if ENABLE_GTK_DOC
SUBDIRS += gtk-doc
endif

EXTRA_DIST = README README.md \
	     img/squirrel.jpg img/squirrel-refocussed.jpg \
	     img/wilber-unsharp.png img/wilber-refocussed.png
nodist_EXTRA_DATA = .git .github .dep .lib
DISTCHECK_CONFIGURE_FLAGS = --disable-silent-rules

strip:
	${STRIP} ${builddir}/src/refocus

if HAVE_ATLAS
dist-hook:
	@echo "**** Don't run make dist when system is configured with atlas."
	@echo "**** Rename lib-atlas, rm config.cache and re-run configure."
	@false
else
dist-hook:
endif

# These use gimptool to install/uninstall refocus in user directory.
# NOTE: gimptool installs refocus directly into the plugin directory
install-bin: ${builddir}/src/refocus
	${GIMPTOOL} --install-bin ${builddir}/src/refocus

uninstall-bin:
	${GIMPTOOL} --uninstall-bin refocus

# These use gimptool to install/uninstall refocus in plugin directory.
# This is for backwards compatibility with v0.9.0, the recommended way
# is to preferrably use 'sudo make install' and 'sudo make uninstall'.
install-admin-bin: ${builddir}/src/refocus
	${GIMPTOOL} --install-admin-bin ${builddir}/src/refocus

uninstall-admin-bin:
	${GIMPTOOL} --uninstall-admin-bin refocus
