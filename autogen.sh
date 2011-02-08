#!/bin/sh

echo "############################################"
echo "# We offer cmake as an alternative!        #"
echo "# consider running 'cmake .; make'         #"
echo "# instead of 'autogen.sh; configure; make' #"
echo "############################################"
sleep 1

(cd src; python collect_classes.py)
autoreconf -v -i -Wall

cp configure configure.tmp
cat >configure <<EOF
#!/bin/sh
echo "############################################"
echo "# We offer cmake as an alternative!        #"
echo "# Please consider running 'cmake .; make'  #"
echo "# instead of './configure; make'           #"
echo "############################################"
sleep 1
EOF
cat configure.tmp >> configure
rm configure.tmp
