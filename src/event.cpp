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

void event::setup_event(TH1I* hist_mult, TH1D* hist_eta){

    this->setmultiplicity(hist_mult);
    vertex.generate_VTX();
    prtl.generate_phi();
    prtl.generate_theta(hist_eta);
    prtl.set_point(vertex);

}

void event::noise(TClonesArray* HitL1, TClonesArray* HitL2){

    //add noise on L1
    int nNoiseL1 = gRandom->Poisson(noise_mu); 

    for (int i = 0; i < nNoiseL1; ++i) {
        
        point p;

        double phi = gRandom->Uniform(-M_PI, M_PI);
        double z = gRandom->Uniform(-layer1_lenght * 0.5, +layer1_lenght * 0.5);
        double R = 40; //mm

        p.set_cil_coord(phi, R, z);

        new((*HitL1)[counterL1++]) point(p);

    }

    //add noise on L2
    int nNoiseL2 = gRandom->Poisson(noise_mu); 

    for (int i = 0; i < nNoiseL2; ++i) {
        
        point p;

        double phi = gRandom->Uniform(-M_PI, M_PI);
        double z = gRandom->Uniform(-layer2_lenght * 0.5, +layer2_lenght * 0.5);
        double R = 70; //mm

        p.set_cil_coord(phi, R, z);

        new((*HitL2)[counterL2++]) point(p);

    }

}

void event::single_event(TClonesArray* VTX, TClonesArray* HitL1, TClonesArray* HitL2){

    //reset counters for new event
    counterVTX = 0;
    counterL1 = 0;
    counterL2 = 0;
    
    //inizialisation of the event
    vertex.generate_VTX();
    if(abs(vertex.get_z()) < beam_pipe_lenght * 0.5) new((*VTX)[counterVTX++]) point(vertex);
    
    //interaction with beam pipe
    prtl.multiple_scattering(beam_pipe_Z, beam_pipe_X0, beam_pipe_thickness);

    //interaction with layer1
    prtl.find_intersection(layer1_radius);
    prtl.multiple_scattering(layer1_Z, layer1_X0, layer1_thickness);
    prtl.get_point().smearing();
    //save interactions in the ttree only if in acceptance
    if(abs(prtl.get_point().get_z()) < layer1_lenght * 0.5) new((*HitL1)[counterL1++]) point(prtl.get_point());
    
    //interaction with layer2
    prtl.find_intersection(layer2_radius);
    prtl.get_point().smearing();
    if(abs(prtl.get_point().get_z()) < layer2_lenght * 0.5) new((*HitL2)[counterL2++]) point(prtl.get_point());

    this->noise(HitL1, HitL2);
    
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

