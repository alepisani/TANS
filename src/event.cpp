#include "../include/event.h"
#include "../include/point.h"
#include "../include/particle.h"
#include "../include/const.h"
#include "../include/reconstruction.h"
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TBranch.h"
#include "TGraphAsymmErrors.h"
#include "TH2D.h"
#include "TClonesArray.h"
#include <algorithm>
#include "TApplication.h" 
#include "TVirtualGeoTrack.h" 
#include "TGeoTrack.h"
#include "TPolyLine3D.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include <iostream>
#include "TF1.h"
#include <sstream>
#include <iomanip>

using namespace std;

event::event(int m, particle prtl):TObject(), multiplicity(m), prtl(prtl){}

void event::setmultiplicity(TH1I* hist_mult) {

    if(!get_data_from_kinem){

        multiplicity = static_cast<int>(gRandom->Uniform(1, 50));

    }

    if(get_data_from_kinem){

        if(hist_mult) {

            multiplicity = hist_mult->GetRandom();

        }
        
    }

}

void event::set_vertex(const point& vtx){

    prtl.set_point(vtx);

}


void event::noise(TClonesArray* HitL1, TClonesArray* HitL2){
    
    //add noise on L1
    int nNoiseL1 = gRandom->Poisson(noise_mu); 
    
    for (int i = 0; i < nNoiseL1; ++i) {
        
        point p;
        
        double phi = gRandom->Uniform(-M_PI, M_PI);
        double z = gRandom->Uniform(-layer1_lenght * 0.5, +layer1_lenght * 0.5);
        double R = layer1_radius; 
        
        p.set_cil_coord(phi, R, z);
        
        new((*HitL1)[counterL1++]) point(p);
        
    }
    
    //add noise on L2
    int nNoiseL2 = gRandom->Poisson(noise_mu); 
    
    for (int i = 0; i < nNoiseL2; ++i) {
        
        point p;
        
        double phi = gRandom->Uniform(-M_PI, M_PI);
        double z = gRandom->Uniform(-layer2_lenght * 0.5, +layer2_lenght * 0.5);
        double R = layer2_radius;
        
        p.set_cil_coord(phi, R, z);
        
        new((*HitL2)[counterL2++]) point(p);
        
    }
    
}

void event::single_event(TClonesArray* VTX, TClonesArray* HitL1, TClonesArray* HitL2, TH1I* hist_mult, TH1D* hist_eta){

    counterVTX = 0; 
    counterL1 = 0; 
    counterL2 = 0;
    
    this->setmultiplicity(hist_mult); 
    vertex.generate_VTX();
    
    if(abs(vertex.get_z()) < beam_pipe_lenght * 0.5) 
        new((*VTX)[counterVTX++]) point(vertex);
    
    int mult = this->get_multiplicity();
    for(int i = 0; i < mult; i++){
        
        //reset particle
        prtl.set_point(vertex);
        prtl.generate_phi();
        prtl.generate_theta(hist_eta);
        
        // BEAM PIPE
        prtl.find_intersection(beam_pipe_radius);
        prtl.multiple_scattering(beam_pipe_Z, beam_pipe_X0, beam_pipe_thickness);
    
        // LAYER 1
        prtl.find_intersection(layer1_radius);
        //smearing only in the ttree, the particle will be transported without smearing
        point& p1 = prtl.get_point(); 
        if(abs(p1.get_z()) < layer1_lenght * 0.5) {
            point* storedL1 = new((*HitL1)[counterL1++]) point(p1);
            // apply the smearing directly in the array
            storedL1->smearing(); 
        }
        prtl.multiple_scattering(layer1_Z, layer1_X0, layer1_thickness);
        
        // LAYER 2
        prtl.find_intersection(layer2_radius);
        point& p2 = prtl.get_point(); 
        if(abs(p2.get_z()) < layer2_lenght * 0.5) {
            point* storedL2 = new((*HitL2)[counterL2++]) point(p2);
            storedL2->smearing();
        }
    }
    
    this->noise(HitL1, HitL2);
}

void event::reset() {
    
    multiplicity = 0;
    counterVTX = 0;
    counterL1 = 0;
    counterL2 = 0;
    vertex.reset();
    prtl.reset();
    
}

std::ostream &operator<<(std::ostream &output, const event & ev) {
    output << "+-----------------------------------evento" << endl;
    output << "| x del vertice: " << ev.vertex.get_x() << endl;
    output << "| y del vertice: " << ev.vertex.get_y() << endl;
    output << "| z del vertice: " << ev.vertex.get_z() << endl;
    output << "| multiplicity: " << ev.multiplicity << endl;
    output << "+-----------------------------------------" << endl;
    return output;
}

