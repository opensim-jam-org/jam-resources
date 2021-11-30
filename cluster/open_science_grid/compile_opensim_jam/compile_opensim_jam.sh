#!/bin/sh
#
#
#

module load GCC/8.3.0
module load CMake/3.16.5

# Lapack
########
mkdir lapack-build lapack-install
git clone https://github.com/Reference-LAPACK/lapack-release.git lapack-source
cd lapack-build
cmake ../lapack-source/ -DCMAKE_INSTALL_PREFIX=../lapack-install -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_SHARED_LIBS=ON
make –j8
make install -j8
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH: ~/lapack-install/lib64
cd ..

# Simbody
#########
mkdir simbody-build simbody-install
git clone https://github.com/simbody/simbody.git simbody-source
#git checkout Simbody-3.6
cd simbody-build
cmake ../simbody-source -D BUILD_EXAMPLES=OFF -D BUILD_STATIC_LIBRARIES=ON -D BUILD_USING_OTHER_LAPACK="~/lapack-install/lib64/liblapack.so;~/lapack-install/lib64/libblas.so" -D CMAKE_INSTALL_PREFIX=../simbody-install
make -j8
ctest -j8
make install -j8

cd ..

# OpenSim
#########
mkdir opensim-build opensim-install opensim-build-dependencies opensim-install-dependencies
git clone https://github.com/opensim-jam-org/opensim-core.git opensim-core


cd opensim-build-dependencies 

cmake -D CMAKE_BUILD_TYPE=Release -D CMAKE_INSTALL_PREFIX=$_CONDOR_SCRATCH_DIR/opensim-install-dependencies -D SUPERBUILD_BTK=OFF -D SUPERBUILD_simbody=OFF ../opensim-core/dependencies

make -j8

cd ../opensim-build

cmake -D CMAKE_BUILD_TYPE=Release -D CMAKE_INSTALL_PREFIX=$_CONDOR_SCRATCH_DIR/opensim-install -D OPENSIM_DEPENDENCIES_DIR=$_CONDOR_SCRATCH_DIR/opensim-install-dependencies -D BUILD_JAVA_WRAPPING=false -D BUILD_PYTHON_WRAPPING=false -D docopt_DIR=$_CONDOR_SCRATCH_DIR/opensim-install-dependencies/docopt/lib64/cmake/docopt -D SIMBODY_HOME=../simbody-install -D OPENSIM_INSTALL_UNIX_FHS=false ../opensim-core

make -j8

ctest -j8

make install -j8

cd ..

#Collect everything in package for HTC distribution
JAM_DIR=opensim-jam
mkdir $JAM_DIR 

SIMBODY=~/simbody-install/lib64
OPENSIM_EXE=~/opensim-install/bin
OPENSIM_LIB=~/opensim-install/sdk/lib
LAPACK=~/lapack-install/lib64

chmod +x $SIMBODY/libSimTKcommon.so
chmod +x $SIMBODY/libSimTKcommon.so.3.8
chmod +x $SIMBODY/libSimTKmath.so
chmod +x $SIMBODY/libSimTKmath.so.3.8
chmod +x $SIMBODY/libSimTKsimbody.so
chmod +x $SIMBODY/libSimTKsimbody.so.3.8

cp $SIMBODY/* $JAM_DIR/
cp $OPENSIM_EXE/* $JAM_DIR/
cp $OPENSIM_LIB/* $JAM_DIR/
cp $LAPACK/libblas.so.3 $JAM_DIR/
cp $LAPACK/liblapack.so.3 $JAM_DIR/

cp ~/opensim-install-dependencies/adol-c/lib64/* $JAM_DIR/
cp ~/opensim-install-dependencies/ipopt/lib/* $JAM_DIR/
cp ~/opensim-install-dependencies/hdf5/lib/* $JAM_DIR/
cp ~/opensim-install-dependencies/colpack/lib64/libColPack.so.0 $JAM_DIR/

cp /lib64/libpthread.so.0 $JAM_DIR/
cp /lib64/librt.so.1 $JAM_DIR/
cp /lib64/libdl.so.2 $JAM_DIR/
cp /software/chtc/eb/GCCcore/8.3.0/lib64/libstdc++.so.6 $JAM_DIR/
cp /lib64/libm.so.6 $JAM_DIR/
cp /software/chtc/eb/GCCcore/8.3.0/lib64/libgcc_s.so.1 $JAM_DIR/
cp /lib64/libc.so.6 $JAM_DIR/
cp /lib64/ld-linux-x86-64.so.2 $JAM_DIR/
cp /software/chtc/eb/GCCcore/8.3.0/lib64/libgfortran.so.5 $JAM_DIR/
cp /software/chtc/eb/GCCcore/8.3.0/lib64/libquadmath.so.0 $JAM_DIR/
cp /software/chtc/eb/zlib/1.2.11-GCCcore-8.3.0/lib/libz.so.1 $JAM_DIR/

tar czf opensim-jam.tar.gz $JAM_DIR

cd $JAM_DIR/

export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH

ldd opensim-cmd



# tar czf opensim-jam-install.tar.gz install 

# tar -czf opensim-core-install.tar.gz opensim-install 

#Copy lapack, simbody, and OpenSim into JAM
# BIN=~/opensim-jam/opensim-jam-release/bin/centos7/
# SIMBODY=~/opensim-jam/opensim-jam-release/opensim/centos7/sdk/Simbody/lib64
# OPENSIM_EXE=~/opensim-jam/opensim-jam-release/opensim/centos7/bin
# OPENSIM_LIB=~/opensim-jam/opensim-jam-release/opensim/centos7/sdk/lib
# LAPACK=~/lapack-install/lib64
# cp -r ~/opensim-install/* ~/opensim-jam/opensim-jam-release/opensim/centos7/

# chmod +x $SIMBODY/libSimTKcommon.so
# chmod +x $SIMBODY/libSimTKcommon.so.3.7
# chmod +x $SIMBODY/libSimTKmath.so
# chmod +x $SIMBODY/libSimTKmath.so.3.7
# chmod +x $SIMBODY/libSimTKsimbody.so
# chmod +x $SIMBODY/libSimTKsimbody.so.3.7

# cp $SIMBODY/* $BIN/
# cp $OPENSIM_EXE/* $BIN/
# cp $OPENSIM_LIB/* $BIN/
# cp $LAPACK/libblas.so.3 $BIN/
# cp $LAPACK/liblapack.so.3 $BIN/

# OpenSim Jam 
#############
# mkdir build build-dependencies build-dependencies/HDF5 install install-dependencies install-dependencies/HDF5
# cd build-dependencies/HDF5/

# cmake -D CMAKE_INSTALL_PREFIX=$_CONDOR_SCRATCH_DIR/install-dependencies/HDF5 -D CMAKE_BUILD_TYPE=Release -D BUILD_SHARED_LIBS=False -D HDF5_BUILD_TOOLS=False ../../opensim-jam/dependencies/HDF5/
# make -j8

# make install -j8

# cd ../../build

# cmake -D CMAKE_INSTALL_PREFIX=$_CONDOR_SCRATCH_DIR/install -D OpenSim_DIR=../opensim-jam/opensim-jam-release/opensim/centos7/cmake -D CMAKE_BUILD_TYPE=Release ../opensim-jam

# make -j8

# make install -j8 

# cd ..

#update OS dependencies (use ldd on jam .exe or .dll to see list and locations)
# cp /lib64/libpthread.so.0 $BIN
# cp /lib64/librt.so.1 $BIN
# cp /lib64/libdl.so.2 $BIN
# cp /software/chtc/eb/GCCcore/8.3.0/lib64/libstdc++.so.6 $BIN
# cp /lib64/libm.so.6 $BIN
# cp /software/chtc/eb/GCCcore/8.3.0/lib64/libgcc_s.so.1 $BIN
# cp /lib64/libc.so.6 $BIN
# cp /lib64/ld-linux-x86-64.so.2 $BIN
# cp /software/chtc/eb/GCCcore/8.3.0/lib64/libgfortran.so.5 $BIN
# cp /software/chtc/eb/GCCcore/8.3.0/lib64/libquadmath.so.0 $BIN
# cp /software/chtc/eb/zlib/1.2.11-GCCcore-8.3.0/lib/libz.so.1 $BIN


# tar czf opensim-jam-install.tar.gz install 