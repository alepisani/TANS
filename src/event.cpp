#include "../include/event.h"
#include "../include/point.h"
#include "../include/particle.h"
#include "../include/const.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include <algorithm>
#include "TApplication.h" 
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMaterial.h"
#include "TGeoTube.h"
#include "TColor.h"
#include "TVirtualGeoTrack.h" 
#include "TGeoTrack.h"
#include "TPolyLine3D.h"
#include "TRandom3.h"
#include <iostream>
using namespace std;

event::event(int m, point p, point vtx):TObject(), multiplicity(m), pnt(p), vertex(vtx){}

void event::setmultiplicity() {
    
    multiplicity = static_cast<int>(gRandom->Uniform(1, 50));
    //multiplicity = fMultiplicityHist->GetRandom();
    //multiplicity = 5000; 

    /*TF1 *f = new TF1("f", "gaus", -5, 5);

    TH1D *h = new TH1D("h", "PDF", 1000, -5, 5);
    h->FillRandom("f", 1e6);   // campiona la funzione
    h->Scale(1.0 / h->Integral()); // normalizzazione (opzionale)
    multiplicity = h->GetRandom();*/

}

void event::set_vertex(point vtx){

    vertex = vtx;

}

void event::eventsimulation(){

    point vertex;
    vertex.generate_VTX();
    cout << vertex << endl;
    
    this->set_vertex(vertex);
    this->setmultiplicity();

    for (int i = 0; i < multiplicity; i++){
        
        particle prtl;
        prtl.set_point(vertex);
        prtl.generate_theta();
        prtl.generate_phi();
        //cout << "particle on vertex" << prtl << endl;

        //trasport the particle until Beam Pipe
        prtl.find_intersection(beam_pipe_radius);
        points_BP.push_back(prtl.get_point());
        if (multiple_scattering_on){
        prtl.multiple_scattering(beam_pipe_Z, beam_pipe_X0, beam_pipe_thickness);
        };
        //cout << "particle on BP" << prtl << endl;

        //transport the particle until Layer1
        prtl.find_intersection(layer1_radius);
        //if(abs(prtl.get_point().get_z()) <= layer1_lenght/2){points_L1.push_back(prtl.get_point());}
        points_L1.push_back(prtl.get_point());
        if (multiple_scattering_on){
        prtl.multiple_scattering(layer1_Z, layer1_X0, layer1_thickness);
        };
        //cout << "particle on L1" << prtl << endl;

        //transport the particle until Layer2
        prtl.find_intersection(layer2_radius);
        //if(abs(prtl.get_point().get_z()) <= layer2_lenght/2){points_L2.push_back(prtl.get_point());}
        points_L2.push_back(prtl.get_point());
        //cout << "particle on L2" << prtl << endl;
 
    

    }

}

void event::display_event(){
    
    //inizializzo ambiente
    TGeoManager *geom = new TGeoManager("Hierarchy", "Esempio Tracker");
    
    // definizione materiali
    TGeoMaterial *mat = new TGeoMaterial("Vacuum", 0, 0, 0);
    TGeoMedium *med = new TGeoMedium("Vacuum", 1, mat);
    
    // TOP Volume
    TGeoVolume *top = geom->MakeBox("TOP", med, 100, 100, 100);
    geom->SetTopVolume(top);
    
    // Geometria (Beam pipe e Layer)
    TGeoVolume *beampipe = geom->MakeTube("beam_pipe", med, beam_pipe_radius, beam_pipe_radius + beam_pipe_thickness, beam_pipe_lenght / 2.);
    beampipe->SetLineColor(kGreen);
    top->AddNode(beampipe, 1);
    
    //layer1 interno
    TGeoVolume *layer1 = geom->MakeTube("Interno", med, layer1_radius, layer1_radius + layer1_thickness, layer1_lenght / 2.);
    layer1->SetLineColor(kRed);
    top->AddNode(layer1, 1);
    
    //layer2 esterno
    TGeoVolume *layer2 = geom->MakeTube("Esterno", med, layer2_radius, layer2_radius + 1, layer2_lenght / 2.);
    layer2->SetLineColor(kBlue);
    top->AddNode(layer2, 1);
    
    geom->CloseGeometry();
    
    // Disegno
    top->Draw("ogl");

    //creazione delle tracce
    int n_punti = 2;
    for(int i = 0; i < multiplicity; i++){

        if(abs(points_L1[i].get_z()) <= layer1_lenght/2){

            TPolyLine3D *trkl_vtx_bp = new TPolyLine3D(n_punti);   
            trkl_vtx_bp->SetPoint(0, vertex.get_x(), vertex.get_y(), vertex.get_z());
            trkl_vtx_bp->SetPoint(1, points_L1[i].get_x(), points_L1[i].get_y(), points_L1[i].get_z());
            trkl_vtx_bp->SetLineColor(kGreen);
            trkl_vtx_bp->SetLineWidth(2);
            trkl_vtx_bp->Draw("same");
            
            TPolyLine3D *trkl_bp_l1 = new TPolyLine3D(n_punti);   
            trkl_bp_l1->SetPoint(0, points_BP[i].get_x(), points_BP[i].get_y(), points_BP[i].get_z());
            trkl_bp_l1->SetPoint(1, points_L1[i].get_x(), points_L1[i].get_y(), points_L1[i].get_z());
            trkl_bp_l1->SetLineColor(kGreen);
            trkl_bp_l1->SetLineWidth(2);
            trkl_bp_l1->Draw("same");
            
            
            if(abs(points_L2[i].get_z()) <= layer2_lenght/2){
                
                TPolyLine3D *trkl_l1_l2 = new TPolyLine3D(n_punti);   
                trkl_l1_l2->SetPoint(0, points_L1[i].get_x(), points_L1[i].get_y(), points_L1[i].get_z());
                trkl_l1_l2->SetPoint(1, points_L2[i].get_x(), points_L2[i].get_y(), points_L2[i].get_z());
                trkl_l1_l2->SetLineColor(kGreen);
                trkl_l1_l2->SetLineWidth(2);
                trkl_l1_l2->Draw("same");
                
            }
        
        }
            
    }

}

void event::RunFullSimulation() {
    
    int nEvents = 1000;

    TFile* hfile = new TFile("../data/hist_sim.root", "RECREATE");
    TTree* tree = new TTree("Tree", "Tree simulazione");

    // Buffer per i dati scalari (Vertice e Molteplicità)
    // Usiamo una struct o variabili singole per la Leaf List
    double vtx[3]; // X, Y, Z
    int multiplicity;

    // Buffer per i dati vettoriali (Punti nei vari Layer)
    // TClonesArray riduce l'overhead di allocazione/deallocazione
    TClonesArray* hitsBP = new TClonesArray("point", 200);
    TClonesArray* particle_VTX_BP = new TClonesArray("particle", 200);
    TClonesArray* particle_BP_L1 = new TClonesArray("particle", 200);
    TClonesArray* particle_L1_L2 = new TClonesArray("particle", 200);
    TClonesArray* hitsL1 = new TClonesArray("point", 200);
    TClonesArray* hitsL2 = new TClonesArray("point", 200);

    // Definizione dei Branch
    tree->Branch("Vtx", vtx, "vtxX/D:vtxY/D:vtxZ/D");
    tree->Branch("Mult", &multiplicity, "mult/I");
    tree->Branch("HitsBP", &hitsBP);
    tree->Branch("Particle_VTX_BP", &particle_VTX_BP);
    tree->Branch("Particle_BP_L1", &particle_BP_L1);
    tree->Branch("Particle_L1_L2", &particle_L1_L2);
    tree->Branch("HitsL1", &hitsL1);
    tree->Branch("HitsL2", &hitsL2);

    event evSim;

    for (int iEv = 0; iEv < nEvents; ++iEv) {
        
        point vertex;
        vertex.generate_VTX();
        evSim.set_vertex(vertex);
        evSim.setmultiplicity();
        
        vtx[0] = vertex.get_x();
        vtx[1] = vertex.get_y();
        vtx[2] = vertex.get_z();
        multiplicity = evSim.get_multiplicity();

        for (int iPart = 0; iPart < multiplicity; ++iPart) {
            particle prtl;
            prtl.set_point(vertex);
            prtl.generate_theta();
            prtl.generate_phi();

            //trasporto da vtx a bp
            new((*particle_VTX_BP)[iPart]) particle(prtl);
            prtl.find_intersection(beam_pipe_radius);
            new((*hitsBP)[iPart]) point(prtl.get_point());  
            points_BP.push_back(prtl.get_point());
            prtl.multiple_scattering(beam_pipe_Z, beam_pipe_X0, beam_pipe_thickness);

            //trasporto da bp a l1
            new((*particle_BP_L1)[iPart]) particle(prtl);
            prtl.find_intersection(layer1_radius);
            prtl.get_point().smearing();
            new((*hitsL1)[iPart]) point(prtl.get_point());
            points_L1.push_back(prtl.get_point());
            prtl.multiple_scattering(layer1_Z, layer1_X0, layer1_thickness);

            //trasporto l1 bp a l2
            new((*particle_L1_L2)[iPart]) particle(prtl);
            prtl.find_intersection(layer2_radius);
            prtl.get_point().smearing();
            new((*hitsL2)[iPart]) point(prtl.get_point());
            points_L2.push_back(prtl.get_point());

        }

        double PNoiseL1 = gRandom->Rndm();
        double PNoiseL2 = gRandom->Rndm();

        if(PNoiseL1 < 0.001){
            point noisy_point;
            noisy_point.set_phi(gRandom->Uniform(0, 2 * M_PI));
            noisy_point.set_R(layer1_radius);
            noisy_point.set_z(gRandom->Uniform(-layer1_lenght/2, +layer1_lenght/2));
            new((*hitsL1)[multiplicity]) point(noisy_point);
        }

        if(PNoiseL2 < 0.001){
            point noisy_point;
            noisy_point.set_phi(gRandom->Uniform(0, 2 * M_PI));
            noisy_point.set_R(layer2_radius);
            noisy_point.set_z(gRandom->Uniform(-layer2_lenght/2, +layer2_lenght/2));
            new((*hitsL2)[multiplicity]) point(noisy_point);
        }        

        tree->Fill();

        hitsBP->Clear("C");
        particle_VTX_BP->Clear("C");
        particle_BP_L1->Clear("C");
        particle_L1_L2->Clear("C");
        hitsL1->Clear("C");
        hitsL2->Clear("C");

        if (iEv % 100 == 0) std::cout << "Event " << iEv << " has been simulated." << std::endl;
    }


    hfile->Write();
    hfile->Close();

    delete hitsBP;
    delete hitsL1;
    delete hitsL2;
    delete particle_VTX_BP;
    delete particle_BP_L1;
    delete particle_L1_L2;
    delete hfile;

    std::cout << "Simulation ended: all event has been processed succesfully. Data available on /data/hist_sim.root" << std::endl;

}

std::ostream &operator<<(std::ostream &output, const event & ev) {
    output << "+-----------------------------------evento" << endl;
    output << "| x del vertice: " << ev.pnt.get_x() << endl;
    output << "| y del vertice: " << ev.pnt.get_y() << endl;
    output << "| z del vertice: " << ev.pnt.get_z() << endl;
    output << "| multiplicity: " << ev.multiplicity << endl;
    output << "+-----------------------------------------" << endl;
    return output;
}

