#!/bin/sh

#============
# Run Walking
#============

JAM_DIR=$PWD/shared/opensim-jam

module load GCC/8.3.0
module load CMake/3.16.5

export LD_LIBRARY_PATH=$JAM_DIR:$LD_LIBRARY_PATH
export PATH=$JAM_DIR:$PATH

mkdir results

tar xzf shared.tar.gz
tar xzf shared/opensim-jam.tar.gz -C shared

#COMAK Inverse Kinematics
# $JAM_DIR/opensim-cmd ./inputs/comak_inverse_kinematics_settings.xml
# mv opensim.log results/comak_ik_opensim.log

#COMAK
opensim-cmd run-tool comak_settings.xml
#mv opensim.log results/comak_opensim.log

#Inverse Dynamics
opensim-cmd run-tool inverse_dynamics_settings.xml

#Joint Mechanics Analysis
opensim-cmd run-tool jnt_mech_settings.xml 
#mv opensim.log results/joint_mechanics_opensim.log

#tar czf results.tar.gz results

mv results/joint-mechanics/walking.h5 .
tar czf walking.h5.tar.gz walking.h5
rm opensim.log
rm walking.h5
