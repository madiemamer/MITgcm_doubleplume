#!/bin/bash
printf "Removing old build and runs folders."
rm -r build
rm -r run

printf "Making new directories."
mkdir build
mkdir run

printf "Building the model."
cd build
../../../tools/genmake2 -mods ../code
make depend FFLAGS="-O2 -w -fallow-argument-mismatch"
make FFLAGS="-O2 -w -fallow-argument-mismatch"

printf "Copying and linking to the run folder."
cd ../run
ln -s ../input/* .
cp ../build/mitgcmuv .

printf "Running the model."
./mitgcmuv

