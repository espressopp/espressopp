#! /bin/sh -x

rm -f config.cache acconfig.h

if [ -n "`which libtoolize`" ]; then
    libtoolize --copy
fi

aclocal -I build-aux/macros		&& \
autoheader				&& \
automake --copy --foreign --add-missing	&& \
autoconf				&& \
exit 0

exit 1

