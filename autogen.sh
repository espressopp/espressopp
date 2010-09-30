#!/bin/sh
(cd src; python collect_classes.py)
exec autoreconf -v -i -Wall
