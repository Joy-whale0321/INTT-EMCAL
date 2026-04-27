AC_INIT(trackcalomatch,[1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE
AC_PROG_CXX
LT_INIT([disable-static])

if test "$GXX" = yes; then
   CXXFLAGS="$CXXFLAGS -Wall -Wextra -Wshadow -Wno-unused-variable -Wno-unused-but-set-variable"
fi

CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT