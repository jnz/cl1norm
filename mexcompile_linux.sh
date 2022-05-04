#!/bin/bash

echo -n "Please enter path to mex application > "
read -i "/usr/local/MATLAB/R2016b/bin/mex" -e pathmex
echo "You entered: $pathmex"
#check if file exist
if [ ! -f $pathmex ]
then
  echo "File doesnt exists"
  exit
fi

CCPATH="/usr/bin/gcc-4.9"
CPPPATH="/usr/bin/g++-4.9"

if [ ! -f $CCPATH ] || [ ! -f $CPPPATH ]
then
  echo "GCC not found"
  echo "Install GCC 4.9: sudo apt-get install gcc-4.9"
  exit
fi

pushd build

$pathmex GCC='/usr/bin/g++-4.9' -O ../src/cl1norm.cpp

popd
