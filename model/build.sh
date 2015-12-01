#!/bin/sh
f2py -c mfix.f90 -m mfix --build-dir blah  *.o dmp_modules/*.o des/*.o des/pic/*.o cartesian_grid/*.o thermochemical/*.o GhdTheory/*.o qmomk/*.o chem/*.o check_data/*.o dqmom/*.o
