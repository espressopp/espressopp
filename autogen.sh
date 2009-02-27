#! /bin/sh -x

rm -f config.cache src/acconfig.hpp

aclocal -I build-aux/macros	&& \
autoheader			&& \
automake --copy --add-missing -Wall	&& \
autoconf			&& \
exit 0

exit 1

