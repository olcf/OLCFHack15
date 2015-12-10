#!/bin/bash

####  A simple automation script to build GCC offload tools, by Philip LoCascio (c) 2015.
skip_download=0
skip_build=1;
## These are the prerequisites for the building system, incase your 
## distro doesn't include them
#wget ftp://gcc.gnu.org/pub/gcc/infrastructure/gmp-4.3.2.tar.bz2
#wget ftp://gcc.gnu.org/pub/gcc/infrastructure/mpfr-2.4.2.tar.bz2
#wget ftp://gcc.gnu.org/pub/gcc/infrastructure/mpc-0.8.1.tar.gz
## Get this date parameter from http://mirrors.concertpass.com/gcc/snapshots/LATEST-5
date=20151013

## SET THIS FIRST AS ROOT OF BUILD
GCC5ROOT="$HOME/Packages/GCCAUTO/"
export CUDA=/usr/local/cuda
if [ "x$GCC5ROOT" = "x" ]; then
echo "set GCC5ROOT to a writable directory, so the build to start. All executables and libraries will be installed under this."
echo "you also need to set CUDA variable to where your version of CUDA is installed. Currently CUDA=$CUDA"
echo "Skip download of sources == $skip_download , current snapshot is date=$date, Skip build  == $skip_build"
exit 0;
fi

## allow us to skip building for testing
if [ $skip_build == 0 ]; then
## Some subdirs
SRC=$GCC5ROOT/src
BUILD=$GCC5ROOT/build
INSTALL=$GCC5ROOT/install
## some build params
threads=16
LCC="gcc-4.9"
LCXX="g++-4.9"

mkdir -p $GCC5ROOT 
cd $GCC5ROOT

export LD_LIBRARY_PATH=$CUDA/lib64:$LD_LIBRARY_PATH
mkdir -p $SRC
cd $SRC
if [ $skip_download == 0 ];then
git clone https://github.com/MentorEmbedded/nvptx-newlib.git
git clone https://github.com/MentorEmbedded/nvptx-tools.git

wget http://mirrors.concertpass.com/gcc/snapshots/LATEST-5/gcc-5-$date.tar.bz2
tar -xjvf gcc-5-$date.tar.bz2
rm -fr gcc
mv gcc-5-$date gcc
fi

echo "Building nvptx-tools : GPU assembler and linking tools"
mkdir -p $BUILD/nvptx-build
cd $BUILD/nvptx-build
$SRC/nvptx-tools/configure    --prefix=$INSTALL   --target=nvptx-none   --with-cuda-driver-include=$CUDA/include   --with-cuda-driver-lib=$CUDA/lib64   --with-cuda-runtime-include=$CUDA/include   --with-cuda-runtime-lib=$CUDA/lib64   CC="$LCC -m64"   CXX="$LCXX -m64" 
make 
make install

echo "Building GCC with nvptx-none target(offload compiler)"
mkdir -p $BUILD/gcc5-accel
cd $BUILD/gcc5-accel
ln -vs $SRC/nvptx-newlib/newlib $SRC/gcc/newlib
ln -vs . $INSTALL/nvptx-none/usr
target=$($SRC/gcc/config.guess)
$SRC/gcc/configure   --prefix=   --target=nvptx-none   --enable-as-accelerator-for="$target"   --enable-languages=c,c++,fortran,lto   --enable-checking=yes,df,fold,rtl   --disable-multilib   --with-sysroot=$INSTALL/nvptx-none   --with-build-sysroot=$INSTALL/nvptx-none   --with-build-time-tools=$INSTALL/nvptx-none/bin   --disable-sjlj-exceptions   --enable-newlib-io-long-long   CC="$LCC -m64"  CXX="$LCXX -m64"
make 
make DESTDIR=$INSTALL install

mkdir -p $BUILD/gcc5
cd $BUILD/gcc5
$SRC/gcc/configure   --prefix=   --disable-bootstrap   --enable-languages=c,c++,fortran,lto   --disable-multilib   --enable-offload-targets=nvptx-none=$INSTALL   --with-cuda-driver-include=$CUDA/include   CC="$LCC -m64"   CXX="$LCXX -m64" --with-sysroot=
make 
make DESTDIR=$INSTALL install



## end if skip_build
echo "Done building. Now try the tests"
fi
echo "Building the gcc5 wrapper script..."
## build the simple wrapper so testing is inplace.
if [ -f rungcc5.TEMPLATE ]; then
cat rungcc5.TEMPLATE |sed "s#_GCC5ROOT_#$GCC5ROOT#g"| sed "s#_CUDA_HOME_#$CUDA#g" > rungcc5.sh
chmod +x rungcc5.sh
echo "ok if you run $PWD/rungcc5.sh with no arguments, it will give you summary of operation."
else
echo "Cannot find rungcc5.TEMPLATE, you will have to set up the paths yourself!"
fi

