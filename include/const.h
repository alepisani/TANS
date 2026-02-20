#ifndef CONST_H
#define CONST_H

//---------------------Simulation----------------------

// Simulated events
const int nEvents = 100000;

/** To simulated the interaction you can either generate the eta and multeplicity 
 *  value by usiform distribution or from the ../data/kinem.root file
 *  true --> data from the file
 *  false --> data from uniform distribution
 */
const bool get_data_from_kinem = false;

// Vertex constant sigmas
const double X_rms = 0.1; //mm
const double Y_rms = 0.1; //mm
const double Z_rms = 53.; //mm

// BEAM PIPE (berillium)
const double beam_pipe_X0 = 350; //mm
const int beam_pipe_Z = 4;
const double beam_pipe_radius = 30; //mm
const double beam_pipe_thickness = 0.8; //mm
const double beam_pipe_lenght = 270; //mm

// LAYER1 (silicon)
const double layer1_X0 = 93.7; //mm 
const int layer1_Z = 14;
const double layer1_radius = 40; //mm
const double layer1_thickness = 0.2; //mm
const double layer1_lenght = 270; //mm

// LAYER2
const double layer2_radius = 70; //mm
const double layer2_lenght = 270; //mm

/** Multiple Scattering on beam pipe and layer1
 *  true --> MS is on
 *  false --> MS is off
 */
const bool multiple_scattering_on = true;


//---------------------Reconstruction----------------------


//Binning to create the Z_vertex reconstruction histogram
const int bin_zvtx = 200;

// Delta phi selection for Z_vertex candidate
const double delta_phi = 0.01; //rad

//Half Running window
const double half_window = 0.05; //mm



#endif



