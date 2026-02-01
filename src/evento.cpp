#include "../include/evento.h"
#include "../include/point.h"
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



//evento::evento() : TObject(), x(0), y(0), z(0), multiplicity(0), rand() {
    //costruttore di default
//}

//evento::evento() {}

//evento::evento(double x, double y, double z, int multiplicity) : TObject(), x(x), y(y), z(z), multiplicity(multiplicity), rand() {
    //costruttore con parametri
//}



void evento::setmultiplicity() {
    multiplicity = static_cast<int>(gRandom->Uniform(1, 5));
}

void evento::display_event(){
    
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
    for(int i = 0; i < trkl_BP_L1.size(); i++){

        TPolyLine3D *trkl = new TPolyLine3D(n_punti);   
        trkl->SetPoint(0, trkl_BP_L1[i].get_point_int().get_x(),trkl_BP_L1[i].get_point_int().get_y(),trkl_BP_L1[i].get_point_int().get_z());
        trkl->SetPoint(1, trkl_BP_L1[i].get_point_ext().get_x(),trkl_BP_L1[i].get_point_ext().get_y(),trkl_BP_L1[i].get_point_ext().get_z());
        trkl->SetLineColor(kRed);
        trkl->SetLineWidth(2);
        trkl->Draw("same");

    }
    
    cout << "_______________" << trkl_BP_L1.size() << "_________" << endl;

}

void evento::generate_vertex(){
    double x = gRandom->Gaus(0,0.1);
    double y = gRandom->Gaus(0,0.1);
    double z = gRandom->Gaus(0,53.);
    vertex.set_point(x,y,z);
}

void evento::event(){
    
    
    //molteplicità --> genera n tracklet
    cout << "AAAAAAAAAAAAAAAAA----------------------------------" << multiplicity << endl;
    for (int i = 0; i < multiplicity; i++){
        
        tracklet trkl_BP_to_L1;
        trkl_BP_to_L1.generate_theta();
        trkl_BP_to_L1.generate_phi();
        trkl_BP_to_L1.generate_eta();
        
        //point point_on_BP;
        //double x_BP = vertex.get_x() + beam_pipe_radius * sin(trkl_BP_to_L1.get_theta()) * sin(trkl_BP_to_L1.get_phi());
        //double y_BP = vertex.get_y() + beam_pipe_radius * cos(trkl_BP_to_L1.get_theta());
        //double z_BP = vertex.get_z() + beam_pipe_radius * sin(trkl_BP_to_L1.get_theta()) * cos(trkl_BP_to_L1.get_phi());
        //point_on_BP.set_point(x_BP, y_BP, z_BP);
        trkl_BP_to_L1.find_beampipe_intersection(vertex.get_x(), vertex.get_y(), vertex.get_z(), trkl_BP_to_L1.get_theta(), trkl_BP_to_L1.get_phi());
        
        trkl_BP_to_L1.set_points(vertex, trkl_BP_to_L1.get_point_ext());
        trkl_BP_L1.push_back(trkl_BP_to_L1);
    }

    for(int i = 0; i < trkl_BP_L1.size(); i++){
        cout << "----" << endl;
        cout << trkl_BP_L1[i] << endl;
        points_BP.push_back(trkl_BP_L1[i].get_point_ext());
    }


    //TODO LIST
    //intersezione tracklet beam pipe
    //crea punti su beam pipe
    //tracklet BP_L1
    //smearing
    //crea punti su L1
    //tracklet L1_L2
    //smearing
    //riempi vettori 

}




std::ostream &operator<<(std::ostream &output, const evento & ev) {
    output << "x del vertice: " << ev.vertex.get_x() << endl;
    output << "y del vertice: " << ev.vertex.get_y() << endl;
    output << "z del vertice: " << ev.vertex.get_z() << endl;
    output << "multiplicity: " << ev.multiplicity << endl;
    return output;
}