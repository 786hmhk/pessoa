dnl  ====================================================
dnl             Check the host platform
dnl  ====================================================

AC_DEFUN([CHECK_PLATFORM],
[


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
      ;;
    i386|i686)
         platform_arch="i386"
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
      ;;
    i386|i686)
         platform_arch="i386"
      ;;
    esac
    ;;
esac
: ${R_PLATFORM_ARCH="$platform_arch"}
AC_SUBST(R_PLATFORM_ARCH)


])
