#! /bin/sh -f

clear

echo
echo "Compiling analysis code:"
make clean
make
ls -tlr lib/
ls -trl bin/


echo 
echo

