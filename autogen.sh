#! /bin/sh -x

rm -f config.cache acconfig.h

# make directory for autoconf-specific helper scripts
if test ! -d build-aux; then
	mkdir build-aux
fi

if [ -n "`which libtoolize`" ]; then
    libtoolize --copy
fi

aclocal -I config			&& \
autoheader				&& \
automake --copy --foreign --add-missing	&& \
autoconf				&& \
exit 0

exit 1

