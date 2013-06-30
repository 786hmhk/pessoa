#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([oc-sp], [0.61], [BUG-REPORT-ADDRESS])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([cudd-2.5.0/obj/cuddObj.hh])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB

# Documentation

AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN";
   then AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])
fi

#AM_CONDITIONAL([HAVE_DOXYGEN],[test -n "$DOXYGEN"])

if test -n "$DOXYGEN"
then
    DOXYGEN_MAIN_MAKEFILE_TARGETS='
doc:
	@echo "${COL_ON}\nCreating Documentation... \c"; \
	cd $(SRC_DIR)/lib && $(MAKE) $@ > /dev/null 2>&1; \
	echo "done! (Location: /src/lib/doc) ${COL_OFF}\n"

'

    DOXYGEN_MAKEFILE_TARGETS='
doc:
	@doxygen documentation.doxygen.Doxyfile > /dev/null 2>&1

'
else
    DOXYGEN_MAIN_MAKEFILE_TARGETS=''
    DOXYGEN_MAKEFILE_TARGETS=''
fi
AC_SUBST([DOXYGEN_MAKEFILE_TARGETS])
AC_SUBST([DOXYGEN_MAIN_MAKEFILE_TARGETS])


### Platform.

AC_CANONICAL_HOST
AC_DEFINE_UNQUOTED(R_PLATFORM, "${host}",
[Define this to be the canonical name (cpu-vendor-os) of your system.])
AC_DEFINE_UNQUOTED(R_CPU, "${host_cpu}",
[Define this to be the name of the CPU of your system.])
AC_DEFINE_UNQUOTED(R_VENDOR, "${host_vendor}",
[Define this to be the name of the vendor of your system.])
AC_DEFINE_UNQUOTED(R_OS, "${host_os}",
[Define this to be the name of the OS of your system.])



R_PLATFORM="${host}"
AC_SUBST(R_PLATFORM)
R_OS="${host_os}"
AC_SUBST(R_OS)
#linux*
#cygwin*|mingw*|windows*|winnt
#darwin*



case "${host_os}" in
  mingw*|windows*|winnt)
    AC_DEFINE(Win32, 1,
              [Define according to your operating system type.])
    R_OSTYPE="windows"
    ;;
  *)
    AC_DEFINE(Unix, 1,
              [Define according to your operating system type.])
    R_OSTYPE="unix"
   ;;
esac
AC_SUBST(R_OSTYPE)

R_CONFIG_ARGS="${ac_configure_args}"
AC_SUBST(R_CONFIG_ARGS)




## We need to establish suitable defaults for a 64-bit OS
## x86_64|mips64|ppc64|powerpc64|sparc64|s390x
## i686|i386|mips|sparc|ppc|
platform_arch=
case "${host_os}" in
  linux*)
    ## 32bit or 64bit
    case "${host_cpu}" in
      x86_64)
         platform_arch="x86_64"
         CUDD_ICFLAGS="-g -O6 -fPIC"
         CUDD_XCFLAGS="-mtune=native -DHAVE_IEEE_754 -DBSD -DSIZEOF_VOID_P=8 -DSIZEOF_LONG=8"
      ;;
    i386|i686)
         platform_arch="i386"
         CUDD_ICFLAGS="-g -O3 -Wall"
         CUDD_XCFLAGS="-mtune=native -malign-double -DHAVE_IEEE_754 -DBSD"
      ;;
    esac
    ;;
  solaris*)
    ## 
    ;;
  darwin*)
    ## 32bit or 64bit
    case "${host_cpu}" in
      x86_64)
         platform_arch="x86_64"
         CUDD_ICFLAGS="-g -O6 -fPIC"
         CUDD_XCFLAGS="-DPIC -DHAVE_IEEE_754 -DSIZEOF_VOID_P=8 -DSIZEOF_LONG=8 -DHAVE_IEEE_754 -fno-common"
      ;;
    i386|i686)
         platform_arch="i386"
         CUDD_ICFLAGS="-g -O3 -Wall"
         CUDD_XCFLAGS="-m32 -DHAVE_IEEE_754 -fno-common"
      ;;
    esac
    ;;
esac
: ${R_PLATFORM_ARCH="$platform_arch"}
AC_SUBST(R_PLATFORM_ARCH)
AC_SUBST(CUDD_ICFLAGS)
AC_SUBST(CUDD_XCFLAGS)


### MATLAB 

AC_ARG_WITH([matlab-path],
  [AS_HELP_STRING([--with-matlab-path],
    [location of the MATLAB root folder.])],
  MATLAB_DIR=${withval},
  MATLAB_DIR=)

##
if test -n "${MATLAB_DIR}"
then
    AC_MSG_CHECKING(for Matlab software)
else
    AC_MSG_ERROR([Must specify matlab directory using --with-matlab-path=])
fi

##
if test -e ${MATLAB_DIR}/bin/mex
then
    AC_MSG_CHECKING(for MATLAB location ...)
else
    AC_MSG_ERROR([Wrong MATLAB root location! Specify matlab directory using --with-matlab-path=])
fi

##
if test -n "${MATLAB_DIR}"
then
    AC_MSG_CHECKING(for Matlab software)

        ## Linux 32bit MATLAB
	if test "`${MATLAB_DIR}/bin/mex -v 2>&1 | grep 'glnx86'`"; then
            MATLAB_MEXSYSTEM="glnx86"
        fi
        ## Linux 64bit MATLAB
	if test "`${MATLAB_DIR}/bin/mex -v 2>&1 | grep 'glnxa64'`"; then
            MATLAB_MEXSYSTEM="glnxa64"
        fi
        ## MAC OS X 32bit MATLAB
	if test "`${MATLAB_DIR}/bin/mex -v 2>&1 | grep 'maci32'`"; then
            MATLAB_MEXSYSTEM="maci32"
        fi
        ## MAC OS X 64bit MATLAB
	if test "`${MATLAB_DIR}/bin/mex -v 2>&1 | grep 'maci64'`"; then
            MATLAB_MEXSYSTEM="maci64"
        fi

        ## CC
        temp1="`${MATLAB_DIR}/bin/mex -v 2>&1 | grep "CC  "`";
        temp2=`echo ${temp1} | sed -e 's/->//'`;
        temp3=`echo ${temp2} | sed -e 's/^[ \t]*//'`;
        MATLAB_CC="${temp3}";
        ## CFLAGS
        temp1="`${MATLAB_DIR}/bin/mex -v 2>&1 | grep "CFLAGS  "`";
        temp2=`echo ${temp1} | sed -e 's/^[ \t]*//'`;
        MATLAB_CFLAGS="${temp2}";
        ## CDEBUGFLAGS
        temp1="`${MATLAB_DIR}/bin/mex -v 2>&1 | grep "CDEBUGFLAGS  "`";
        temp2=`echo ${temp1} | sed -e 's/^[ \t]*//'`;
        MATLAB_CDEBUGFLAGS="${temp2}";
        ## COPTIMFLAGS
        temp1="`${MATLAB_DIR}/bin/mex -v 2>&1 | grep "COPTIMFLAGS  "`";
        temp2=`echo ${temp1} | sed -e 's/^[ \t]*//'`;
        MATLAB_COPTIMFLAGS="${temp2}";
        ## CLIBS
        temp1="`${MATLAB_DIR}/bin/mex -v 2>&1 | grep "CLIBS  "`";
        temp2=`echo ${temp1} | sed -e 's/^[ \t]*//'`;
        MATLAB_CLIBS="${temp2}";
        ## ARGUMENTS
        temp1="`${MATLAB_DIR}/bin/mex -v 2>&1 | grep "arguments  " | sed -n '1p'`";
        temp2=`echo ${temp1} | sed -e 's/^[ \t]*//'`;
        MATLAB_CCARGUMENTS="${temp2}";
        ## CXX
        temp1="`${MATLAB_DIR}/bin/mex -v 2>&1 | grep "CXX  "`";
        temp2=`echo ${temp1} | sed -e 's/->//'`;
        temp3=`echo ${temp2} | sed -e 's/^[ \t]*//'`;
        MATLAB_CXX="${temp3}";
        ## CXXFLAGS
        temp1="`${MATLAB_DIR}/bin/mex -v 2>&1 | grep "CXXFLAGS  "`";
        temp2=`echo ${temp1} | sed -e 's/^[ \t]*//'`;
        MATLAB_CXXFLAGS="${temp2}";
        ## CXXDEBUGFLAGS
        temp1="`${MATLAB_DIR}/bin/mex -v 2>&1 | grep "CXXDEBUGFLAGS  "`";
        temp2=`echo ${temp1} | sed -e 's/^[ \t]*//'`;
        MATLAB_CXXDEBUGFLAGS="${temp2}";
        ## CXXOPTIMFLAGS
        temp1="`${MATLAB_DIR}/bin/mex -v 2>&1 | grep "CXXOPTIMFLAGS  "`";
        temp2=`echo ${temp1} | sed -e 's/^[ \t]*//'`;
        MATLAB_CXXOPTIMFLAGS="${temp2}";
        ## CXXLIBS
        temp1="`${MATLAB_DIR}/bin/mex -v 2>&1 | grep "CXXLIBS  "`";
        temp2=`echo ${temp1} | sed -e 's/^[ \t]*//'`;
        MATLAB_CXXLIBS="${temp2}";
        ## ARGUMENTS
        temp1="`${MATLAB_DIR}/bin/mex -v 2>&1 | grep "arguments  " | sed -n '2p'`";
        temp2=`echo ${temp1} | sed -e 's/^[ \t]*//'`;
        MATLAB_CXXARGUMENTS="${temp2}";
        ## LD
        temp1="`${MATLAB_DIR}/bin/mex -v 2>&1 | grep "LD  "`";
        temp2=`echo ${temp1} | sed -e 's/->//'`;
        temp3=`echo ${temp2} | sed -e 's/^[ \t]*//'`;
        MATLAB_LD="${temp3}";
        ## LDFLAGS
        temp1="`${MATLAB_DIR}/bin/mex -v 2>&1 | grep "LDFLAGS  "`";
        temp2=`echo ${temp1} | sed -e 's/^[ \t]*//'`;
        MATLAB_LDFLAGS="${temp2}";
        ## LDDEBUGFLAGS
        temp1="`${MATLAB_DIR}/bin/mex -v 2>&1 | grep "LDDEBUGFLAGS  "`";
        temp2=`echo ${temp1} | sed -e 's/^[ \t]*//'`;
        MATLAB_LDDEBUGFLAGS="${temp2}";
        ## LDOPTIMFLAGS
        temp1="`${MATLAB_DIR}/bin/mex -v 2>&1 | grep "LDOPTIMFLAGS  "`";
        temp2=`echo ${temp1} | sed -e 's/^[ \t]*//'`;
        MATLAB_LDOPTIMFLAGS="${temp2}";
        ## LDEXTENSION
        temp1="`${MATLAB_DIR}/bin/mex -v 2>&1 | grep "LDEXTENSION  "`";
        temp2=`echo ${temp1} | sed -e 's/^[ \t]*//'`;
        MATLAB_LDEXTENSION="${temp2}";
        ## ARGUMENTS
        temp1="`${MATLAB_DIR}/bin/mex -v 2>&1 | grep "arguments  " | sed -n '4p'`";
        temp2=`echo ${temp1} | sed -e 's/^[ \t]*//'`;
        MATLAB_LDARGUMENTS="${temp2}";
        ## LDCXX
        temp1="`${MATLAB_DIR}/bin/mex -v 2>&1 | grep "LDCXX  "`";
        temp2=`echo ${temp1} | sed -e 's/->//'`;
        temp3=`echo ${temp2} | sed -e 's/^[ \t]*//'`;
        MATLAB_LDCXX="${temp3}";
        ## LDCXXFLAGS
        temp1="`${MATLAB_DIR}/bin/mex -v 2>&1 | grep "LDCXXFLAGS  "`";
        temp2=`echo ${temp1} | sed -e 's/^[ \t]*//'`;
        MATLAB_LDCXXFLAGS="${temp2}";
        ## LDCXXDEBUGFLAGS
        temp1="`${MATLAB_DIR}/bin/mex -v 2>&1 | grep "LDCXXDEBUGFLAGS  "`";
        temp2=`echo ${temp1} | sed -e 's/^[ \t]*//'`;
        MATLAB_LDCXXDEBUGFLAGS="${temp2}";
        ## LDCXXOPTIMFLAGS
        temp1="`${MATLAB_DIR}/bin/mex -v 2>&1 | grep "LDCXXOPTIMFLAGS  "`";
        temp2=`echo ${temp1} | sed -e 's/^[ \t]*//'`;
        MATLAB_LDCXXOPTIMFLAGS="${temp2}";
        ## LDCXXEXTENSION
        temp1="`${MATLAB_DIR}/bin/mex -v 2>&1 | grep "LDCXXEXTENSION  "`";
        temp2=`echo ${temp1} | sed -e 's/^[ \t]*//'`;
        MATLAB_LDCXXEXTENSION="${temp2}";
        ## ARGUMENTS
        temp1="`${MATLAB_DIR}/bin/mex -v 2>&1 | grep "arguments  " | sed -n '5p'`";
        temp2=`echo ${temp1} | sed -e 's/^[ \t]*//'`;
        MATLAB_LDCXXARGUMENTS="${temp2}";

	## The MATLAB compiler
	MEX="${MATLAB_DIR}/bin/mex"

    AC_MSG_RESULT($MATLAB_DIR)
    AC_SUBST(MATLAB_DIR)
    AC_SUBST(MATLAB_CC)
    AC_SUBST(MATLAB_CFLAGS)
    AC_SUBST(MATLAB_CDEBUGFLAGS)
    AC_SUBST(MATLAB_COPTIMFLAGS)
    AC_SUBST(MATLAB_CLIBS)
    AC_SUBST(MATLAB_CCARGUMENTS)
    AC_SUBST(MATLAB_CXX)
    AC_SUBST(MATLAB_CXXFLAGS)
    AC_SUBST(MATLAB_CXXDEBUGFLAGS)
    AC_SUBST(MATLAB_CXXOPTIMFLAGS)
    AC_SUBST(MATLAB_CXXLIBS)
    AC_SUBST(MATLAB_CXXARGUMENTS)
    AC_SUBST(MATLAB_LD)
    AC_SUBST(MATLAB_LDFLAGS)
    AC_SUBST(MATLAB_LDDEBUGFLAGS)
    AC_SUBST(MATLAB_LDOPTIMFLAGS)
    AC_SUBST(MATLAB_LDEXTENSION)
    AC_SUBST(MATLAB_LDARGUMENTS)
    AC_SUBST(MATLAB_LDCXX)
    AC_SUBST(MATLAB_LDCXXFLAGS)
    AC_SUBST(MATLAB_LDCXXDEBUGFLAGS)
    AC_SUBST(MATLAB_LDCXXOPTIMFLAGS)
    AC_SUBST(MATLAB_LDCXXEXTENSION)
    AC_SUBST(MATLAB_LDCXXARGUMENTS)
    AC_SUBST(MATLAB_MEXSYSTEM)
    AC_SUBST(MEX)
    AM_CONDITIONAL(HAVE_MATLAB, test "xyes" = "xyes" )
else
    AC_MSG_ERROR([Must specify matlab directory using --with-matlab-path=])
fi


###

AC_ARG_WITH([maintainers-makefile-targets],
  [AS_HELP_STRING([--with-maintainers-makefile-targets],
    [enables Makefile targets for maintainers])],
  [MAINTAINERS_MAKEFILE_TARGETS="$withval"],
  [MAINTAINERS_MAKEFILE_TARGETS=""])


if test -n "${MAINTAINERS_MAKEFILE_TARGETS}"
then
    MAINTAINERS_MAKEFILE_TARGETS='

tarname = ${package}
distdir = $(tarname)-${version}

dist: $(distdir).tar.gz

$(distdir).tar.gz: $(distdir)
	tar chof - $(distdir) | gzip -9 -c > $@
	@rm -rf $(distdir)
	@rm -rf temp_dist

$(distdir): FORCE
	@mkdir -p temp_dist
	-@cp * temp_dist
	@cp -r $(CUDD_DIR) temp_dist
	@cp -r $(SRC_DIR) temp_dist
	@cp -r bin temp_dist
	@find temp_dist -type f -name "Makefile" -exec rm -f {} \;
	find temp_dist -type f -name "*~" -exec rm -f {} \;
	@mkdir -p $(distdir)
	@cd temp_dist; \
	cp Makefile.in configure config.h.in config.sub install-sh  config.guess ../$(distdir); \
	cp -r $(CUDD_DIR) ../$(distdir); \
	cp -r $(SRC_DIR) ../$(distdir); \
	cp -r bin ../$(distdir)

FORCE: clean distclean
	-@rm $(distdir).tar.gz >/dev/null 2>&1
	-@rm -rf $(distdir) >/dev/null 2>&1
	-@rm -rf temp_dist

distcheck: $(distdir).tar.gz
	gzip -cd $(distdir).tar.gz | tar xvf -
	cd $(distdir) && ./configure
	cd $(distdir) && $(MAKE) all
	cd $(distdir) && $(MAKE) clean
	cd $(distdir) && $(MAKE) distclean
	rm -rf $(distdir)
	@echo "${COL_ON}\n\n*** Package $(distdir).tar.gz is ready for distribution. ***\n${COL_OFF}"'
fi


AC_SUBST([MAINTAINERS_MAKEFILE_TARGETS])




# Checks for header files.
AC_CHECK_HEADERS([float.h limits.h stdlib.h string.h sys/file.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_OFF_T
AC_C_RESTRICT
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRCOLL
AC_FUNC_WAIT3
AC_CHECK_FUNCS([bzero dup2 gethostname getpagesize gettimeofday memchr memmove memset pow sqrt strchr strcspn strerror strpbrk strrchr strspn strstr])

AC_CONFIG_FILES([Makefile
                 cudd-2.5.0/cudd/Makefile
                 cudd-2.5.0/dddmp/Makefile
                 cudd-2.5.0/epd/Makefile
                 cudd-2.5.0/mnemosyne/Makefile
                 cudd-2.5.0/mtr/Makefile
                 cudd-2.5.0/nanotrav/Makefile
                 cudd-2.5.0/obj/Makefile
                 cudd-2.5.0/st/Makefile
                 cudd-2.5.0/util/Makefile
                 src/Makefile
		 src/lib/Makefile])
AC_OUTPUT
