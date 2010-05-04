#!/bin/sh
(cd src; bash collect_classes.sh)
exec autoreconf -v -i -Wall
