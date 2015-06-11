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
mpirun -np 1 $CASE_DIR/mfix.exe imax=64 jmax=64 nodesi=1 nodesj=1 nodesk=1 > out.log
#rm $CASE_DIR/{TFM03.*,out.log}
#rm $CASE_DIR/mfix.exe

echo "******** Done."

# uncomment the following to generate plots:
#echo "******** Generating plots..."
#python plot_results.py &
