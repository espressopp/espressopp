#! /bin/sh -x

rm -f config.cache src/cpp/acconfig.hpp

aclocal -I build-aux/macros	&& \
autoheader			&& \
automake --copy --add-missing	&& \
autoconf			&& \
exit 0

exit 1

