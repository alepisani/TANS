#include "../include/particle.h"
#include "../include/point.h"
#include "../include/const.h"
#include "../include/event.h"
#include "../include/reconstruction.h"
#include "TRandom3.h"
#include <iostream>
#include <cmath>
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"

using namespace std;

ClassImp(reconstruction);

double reconstruction::running_window(){

    double window = 0.1; //mm - reduced to improve resolution
    int maxfrequency = 0;
    int frequency = 0;
    double vertex = 187; //dummy value
    double running_sum = 0;  

    //creation of the window to iterate on each value
    for(auto cand_vtx : vertex_candidate){
        
        double window_in = cand_vtx - window/2.; 
        double window_out = cand_vtx + window/2.;
        
        frequency = 0;        
        running_sum = 0;      
        
        //check the window frequency 
        for(auto cand_vtx2 : vertex_candidate){ 

            if(cand_vtx2 > window_in && cand_vtx2 < window_out){
                frequency++;
                running_sum += cand_vtx2;
            }

        }

        //find the maximum and consider the vertex as the mean of the window
        if(frequency > maxfrequency){
            maxfrequency = frequency;
            vertex = running_sum / frequency;  
        }
    }

    return vertex;
}



double reconstruction::reco_z(event* ev, TH1D* hResidui){
    
    // Clear vectors from previous event
    z_candidates.clear();
    vertex_candidate.clear();
    
    //create an histogram to fill with the z candidates
    TH1D hist_z_vtx("hist_reco", "Ricostruzione Vertice Z; z_{cand} [mm]; Conteggi", 200, -beam_pipe_lenght/2., beam_pipe_lenght/2.);

        double z1, z2, r1, r2, phi1, phi2, z_cand;

        //cout << ev->points_L1.size() << endl;
        for(size_t i = 0; i < ev->points_L1.size(); ++i) {  //loop on L1 points
            
            z1 = ev->points_L1[i].get_z();
            r1 = layer1_radius; 
            phi1 = ev->points_L1[i].get_phi();

            for(size_t j = 0; j < ev->points_L2.size(); ++j) {  //loop on L2 points
                
                z2 = ev->points_L2[j].get_z();
                r2 = layer2_radius;
                phi2 = ev->points_L2[j].get_phi();

                //taking into accout the periodicity
                double dphi = std::abs(phi1 - phi2);  
                if (dphi > M_PI) dphi = 2.0 * M_PI - dphi;

                //selection on phi and geometry
                if(dphi <= delta_phi) {
                    z_cand = (r2 * z1 - r1 * z2) / (r2 - r1); //if the points are compatible with a tracklet, finds its vertex
                    if(abs(z_cand) <= beam_pipe_lenght/2.){
                        hist_z_vtx.Fill(z_cand);
                        z_candidates.push_back(z_cand);
                    } 
                }
            }
    }

    int bin = hist_z_vtx.GetMaximumBin();  //find the bin with the maximum frequency of candidates
    double bin_low = hist_z_vtx.GetBinLowEdge(bin); 
    double bin_high = bin_low + hist_z_vtx.GetBinWidth(bin); 

    //fill a vector with all the z values of the maximum bin
    for(double z : z_candidates){

        if(z >= bin_low && z < bin_high){
            vertex_candidate.push_back(z);
        }

    }

    //find primary vertex
    double z_vtx = running_window();
    
    //fill residuals histogram if provided
    //if(hResidui != nullptr && ev->get_multiplicity() > 10){
    if(hResidui != nullptr){
        double residuo = (z_vtx - ev->get_vertex().get_z())*1000 ;
        hResidui->Fill(residuo);
    }

    return z_vtx;

}





