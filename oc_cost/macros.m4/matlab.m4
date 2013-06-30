dnl  ====================================================
dnl             MATLAB library macros
dnl  ====================================================

dnl
dnl  CHECK_MATLAB - find Matlab libraries
dnl  ----------------------------------------------------
AC_DEFUN([CHECK_MATLAB],
[
AC_ARG_WITH(matlab,
[  --with-matlab=DIR	the directory where Matlab is installed ],
MATLAB_DIR=${withval},
MATLAB_DIR=)

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
    AM_CONDITIONAL(HAVE_MATLAB, test "xyes" = "xyes" )
    have_matlab=yes
else
    AM_CONDITIONAL(HAVE_MATLAB, test "xno" = "xyes" )
    have_matlab=no
fi
])
