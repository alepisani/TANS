#ifndef CONST_H
#define CONST_H

/**
 * compressionlevel = 1, faster and doesnt care about file sizing
 * compressionlevel = 9, slower but optimizes the file sizing
 * 1,9 are only the lower and upper boundary, you can choose in [1,9]
 */

constexpr int CompressionLevel = 1;


//---------------------Simulation----------------------

// Simulated events
constexpr int nEvents = 1000000;

/** To simulated the interaction you can either generate the eta and multeplicity 
 *  value by usiform distribution or from the ../data/kinem.root file
 *  true --> data from the file
 *  false --> data from uniform distribution
 */
constexpr bool get_data_from_kinem = true;

// Vertex constant sigmas
constexpr double X_rms = 0.1; //mm
constexpr double Y_rms = 0.1; //mm
constexpr double Z_rms = 53.; //mm

// BEAM PIPE (berillium)
constexpr double beam_pipe_X0 = 350; //mm
constexpr int beam_pipe_Z = 4;
constexpr double beam_pipe_radius = 30; //mm
constexpr double beam_pipe_thickness = 0.8; //mm

// LAYER1 (silicon)
constexpr double layer1_X0 = 93.7; //mm 
constexpr int layer1_Z = 14;
constexpr double layer1_radius = 40; //mm
constexpr double layer1_thickness = 0.2; //mm
constexpr double layer1_lenght = 270; //mm

// LAYER2
constexpr double layer2_radius = 70; //mm
constexpr double layer2_lenght = 270; //mm

/** Multiple Scattering on beam pipe and layer1
 *  true --> MS is on
 *  false --> MS is off
 */
constexpr bool multiple_scattering_on = true;

constexpr int noise_mu = 5;

constexpr double smearing_sigma = 0.03; //mm


//---------------------Reconstruction----------------------

// Delta phi selection for Z_vertex candidate
//constexpr double delta_phi = 0.005; //rad
constexpr double delta_phi = 0.01; //rad

//Half Running window
constexpr double half_window = 1; //mm



#endif



