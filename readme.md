

# PREREQUISITES

## Eigen

## cminpack

See <https://github.com/devernay/cminpack>. The `make` command for cminpack must be passed options as follows
to compile the long double version:

		`make LIBSUFFIX=s CFLAGS="-O3 -g -Wall -Wextra -D__cminpack_long_double__"`

Next, run the command `make install` or `sudo make install`.

In the makefile in the *build* subdirectory of this repository, set the variable `CMINPACKDIR`
to the install path of cminpack.

## SuiteSparse