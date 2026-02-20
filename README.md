# T.A.N.S. project (UNITO)
Alessandro Pisani, Mariachiara Spesso

### Overview
This repository is dedicated to the university course TANS (UNITO).
The intent of this project is to simulate an high energy particle interaction in a collider (beam - beam interaction) and reconstruct the primary vertex.
The geometry of the collider is cilindrindal. To take the measures we have two concentrical detectors layers with inside a beam pipe (berillium). Multiple scattering in taking into account.

### Requirements
ROOT framework.
The code was written with the 6.36.04 version of the framework.

In the directory ./include/const.h you are able to set the configuration of the code. 


## How to run the code?
To compile the code we uses CMAKE. Therefore we need to create the build directory and compile the code from there: 
let's start from the main director.

```
mkdir build
cd build
cmake ..
make -jN
./sim
```

in the directory ./data/hist_sim you'll be able to see all the data from the simulation and the plot from the reconstruction.

## const.h
From the directory .include/const.h you can change all the ccode const and decide some different option how the code works:
* nEvents : number of event simulated and therefore recostructed (default value = 100000);
* get_data_from_kinem : for the simulation you are able to either take the eta and multeplicity values from uniform distribution of from the ./data/kinem.root file. 
* multiple_scattering_on : ask you if you want to take into account the multiple scattering for the simulation

## Output
If not major issues occured you may see the outcome data in /data/hist_sim.root using the ROOT TBrowser in the following way.
```
root -l
TBrowser b
```


Enjoy

