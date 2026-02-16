#ifndef CONST_H
#define CONST_H

// Vertex constant sigmas
const double X_rms = 0.1; //mm
const double Y_rms = 0.1; //mm
const double Z_rms = 53.; //mm


// BEAM PIPE (berillium)
const double beam_pipe_X0 = 350; //mm
const int beam_pipe_Z=4;
const double beam_pipe_radius = 30; //mm
const double beam_pipe_thickness = 0.8; //mm
const double beam_pipe_lenght = 270; //mm

// LAYER1 
const double layer1_X0 = 93.7; //mm X0 del silicio
const int layer1_Z=14;
const double layer1_radius = 40; //mm
const double layer1_thickness = 0.2; //mm
const double layer1_lenght = 270; //mm

// LAYER2
const double layer2_radius = 70; //mm
const double layer2_lenght = 270; //mm

//MS
const bool multiple_scattering_on = true;

//reconstuction constant
const double zBinWidth = 50; //mm
const double delta_phi = 0.01; //rad

//distribution
const bool distrib_assegnata = true;

#endif



