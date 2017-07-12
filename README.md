# compbio
The Penn State University Comp Bio Computing Framework

## How to build on linux and Mac
```
cd compbio
mkdir build
cd build
ccmake ..
make -j 8

## How to build on Windows
Download and install cygwin. After you do the default installation, you will need to also add vim, gcc (c and g++), git, make, cmake, python (version 3 has been tested)

then it is similar to linux and mac instructions, i.e.:
cd compbio
mkdir build
cd build
ccmake ..
