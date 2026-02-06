#include "../include/event.h"
#include "../include/point.h"
#include "../include/particle.h"
#include "../include/const.h"
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
    
    //multiplicity = static_cast<int>(gRandom->Uniform(1, 50));
    //multiplicity = fMultiplicityHist->GetRandom();
    multiplicity = 2; 

}

void event::set_vertex(point vtx){

    vertex = vtx;

}

void event::eventsimulation(){

    point vertex;
    vertex.generate_VTX();
    
    this->set_vertex(vertex);
    this->setmultiplicity();

    for (int i = 0; i < multiplicity; i++){
        
        particle prtl;
        prtl.set_point(vertex);
        prtl.generate_theta();
        prtl.generate_phi();
        cout << "particle on vertex" << prtl << endl;

        //trasport the particle until Beam Pipe
        prtl.find_intersection(beam_pipe_radius);
        points_BP.push_back(prtl.get_point());
        prtl.multiple_scattering(beam_pipe_Z, beam_pipe_X0, beam_pipe_thickness);
        cout << "particle on BP" << prtl << endl;

        //transport the particle until Layer1
        prtl.find_intersection(layer1_radius);
        points_L1.push_back(prtl.get_point());
        prtl.multiple_scattering(layer1_Z, layer1_X0, layer1_thickness);
        cout << "particle on L1" << prtl << endl;

        //transport the particle until Layer2
        prtl.find_intersection(layer2_radius);
        points_L2.push_back(prtl.get_point());
        cout << "particle on L2" << prtl << endl;


    }

}


void event::smearing(){

    //smearing di 30um applicato sull'arco di circonferenza --> sta cambiando il phi del punto
    //30um = 0.03mm

    for(int i = 0; i < points_L1.size(); i++){

        double smearing = gRandom->Gaus(0, 0.03);
        points_L1[i].set_phi(points_L1[i].get_phi() + smearing / points_L1[i].get_R());

    }

    for(int i = 0; i < points_L2.size(); i++){

        double smearing = gRandom->Gaus(0, 0.03);
        points_L2[i].set_phi(points_L2[i].get_phi() + smearing / points_L2[i].get_R());

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

        TPolyLine3D *trkl_vtx_bp = new TPolyLine3D(n_punti);   
        trkl_vtx_bp->SetPoint(0, vertex.get_x(), vertex.get_y(), vertex.get_z());
        trkl_vtx_bp->SetPoint(1, points_L1[i].get_x(), points_L1[i].get_y(), points_L1[i].get_z());
        trkl_vtx_bp->SetLineColor(kBlue);
        trkl_vtx_bp->SetLineWidth(2);
        trkl_vtx_bp->Draw("same");
        
        TPolyLine3D *trkl_bp_l1 = new TPolyLine3D(n_punti);   
        trkl_bp_l1->SetPoint(0, points_BP[i].get_x(), points_BP[i].get_y(), points_BP[i].get_z());
        trkl_bp_l1->SetPoint(1, points_L1[i].get_x(), points_L1[i].get_y(), points_L1[i].get_z());
        trkl_bp_l1->SetLineColor(kBlue);
        trkl_bp_l1->SetLineWidth(2);
        trkl_bp_l1->Draw("same");

        TPolyLine3D *trkl_l1_l2 = new TPolyLine3D(n_punti);   
        trkl_l1_l2->SetPoint(0, points_L1[i].get_x(), points_L1[i].get_y(), points_L1[i].get_z());
        trkl_l1_l2->SetPoint(1, points_L2[i].get_x(), points_L2[i].get_y(), points_L2[i].get_z());
        trkl_l1_l2->SetLineColor(kBlue);
        trkl_l1_l2->SetLineWidth(2);
        trkl_l1_l2->Draw("same");

    }


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






