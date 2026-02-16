This repository is dedicated to a university course, TANS (UNITO).
The intent of this project is to simulate a particle interaction in a collider (beam - beam interaction).
The geometry of the collider is cilindrindal. To take the measures we have two concentrical detectors layers with inside a beam pipe (berillium).
The goal is to reconstruct the interaction vertex. 

In the directory ./include/const.h you are able to set the configuration of the code. 


## how to run the code?
you should already have a build directory, for the sake we should firstly delete it. go in the main directory
```
rm -rf build/
mkdir build
cd build
cmake ..
make
./sim
```
in the directory ./data/hist_sim you'll be able to see all the data from the simulation and the plot from the reconstruction.


