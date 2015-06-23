#!/bin/bash -exl

# set case directory
export CASE_DIR=`pwd`

# load modules
module load gnu/4.6.4 openmpi/1.5.5_gnu4.6

# compile MFIX in ./src/
echo "******** Compiling MFIX..."
cd $CASE_DIR
../../../model/make_mfix --dmp --opt=O3 --compiler=gcc --exe=mfix.exe -j

cd $CASE_DIR

# Run case
echo "******** Running simulation..."

mpirun -np 16 mfix.exe \
  mu_g0=0.01 \
  imax=128 jmax=128 \
  nodesi=4 nodesj=4 nodesk=1 #> out.log
mv u_profile.dat u_profile_Re100_S.dat
mv v_profile.dat v_profile_Re100_S.dat
rm {TFM03.*,out.log}

mpirun -np 16 mfix.exe \
  mu_g0=0.0025 \
  ur_fac=0.25,0.5,0.15,0.15,0.5,0.8,1.0,0.5,0.8 \
  imax=128 jmax=128 \
  nodesi=4 nodesj=4 nodesk=1 #> out.log
mv u_profile.dat u_profile_Re400_S.dat
mv v_profile.dat v_profile_Re400_S.dat
rm {TFM03.*,out.log}

# Re=1000 case on a fine mesh.
# Not for regular testing. Takes a long time to converge.
# Norm_G is NOT set to 0 for this case.
#mpirun -np 64 mfix.exe \
#  mu_g0=0.001 \
#  leq_method(1)=3 \
#  ic_u_g(1)=0.001 \
#  imax=512 jmax=512 \
#  nodesi=8 nodesj=8 nodesk=1 #> out.log
#mv u_profile.dat u_profile_Re1000_S.dat
#mv v_profile.dat v_profile_Re1000_S.dat
#rm {TFM03.*,out.log}

#rm $CASE_DIR/mfix.exe

echo "******** Done."

# uncomment the following to generate plots:
#echo "******** Generating plots..."
#python plot_results.py &
