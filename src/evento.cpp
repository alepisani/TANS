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
    //multiplicity = static_cast<int>(gRandom->Uniform(1, 50));
    multiplicity = 3; 
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
    for(int i = 0; i < trkl_VTX_BP.size(); i++){

        TPolyLine3D *trkl = new TPolyLine3D(n_punti);   
        trkl->SetPoint(0, trkl_VTX_BP[i].get_point_int().get_x(),trkl_VTX_BP[i].get_point_int().get_y(),trkl_VTX_BP[i].get_point_int().get_z());
        trkl->SetPoint(1, trkl_VTX_BP[i].get_point_ext().get_x(),trkl_VTX_BP[i].get_point_ext().get_y(),trkl_VTX_BP[i].get_point_ext().get_z());
        trkl->SetLineColor(kBlue);
        trkl->SetLineWidth(1);
        trkl->Draw("same");

    }
    for(int i = 0; i < trkl_BP_L1.size(); i++){

        TPolyLine3D *trkl = new TPolyLine3D(n_punti);   
        trkl->SetPoint(0, trkl_BP_L1[i].get_point_int().get_x(),trkl_BP_L1[i].get_point_int().get_y(),trkl_BP_L1[i].get_point_int().get_z());
        trkl->SetPoint(1, trkl_BP_L1[i].get_point_ext().get_x(),trkl_BP_L1[i].get_point_ext().get_y(),trkl_BP_L1[i].get_point_ext().get_z());
        trkl->SetLineColor(kBlue);
        trkl->SetLineWidth(6);
        trkl->Draw("same");

    }
    for(int i = 0; i < trkl_L1_L2.size(); i++){

        TPolyLine3D *trkl = new TPolyLine3D(n_punti);   
        trkl->SetPoint(0, trkl_L1_L2[i].get_point_int().get_x(),trkl_L1_L2[i].get_point_int().get_y(),trkl_L1_L2[i].get_point_int().get_z());
        trkl->SetPoint(1, trkl_L1_L2[i].get_point_ext().get_x(),trkl_L1_L2[i].get_point_ext().get_y(),trkl_L1_L2[i].get_point_ext().get_z());
        trkl->SetLineColor(kRed);
        trkl->SetLineWidth(2);
        trkl->Draw("same");

    }
    
    cout << "_______________" << trkl_VTX_BP.size() << "_________" << endl;
    cout << "_______________" << trkl_BP_L1.size() << "_________" << endl;
    cout << "_______________" << trkl_L1_L2.size() << "_________" << endl;

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
        
        // Tracklet vertex -> Beam Pipe
        tracklet trkl_VTX_to_BP;
        trkl_VTX_to_BP.generate_theta();
        trkl_VTX_to_BP.generate_phi();
        trkl_VTX_to_BP.generate_eta();
        trkl_VTX_to_BP.set_point_int(vertex);  
        trkl_VTX_to_BP.set_point_ext(trkl_VTX_to_BP.find_intersection(beam_pipe_radius));
        points_BP.push_back(trkl_VTX_to_BP.get_point_ext());  //vettore con punti di intersezione con la beam pipe, serve per la ricostruzione
        trkl_VTX_BP.push_back(trkl_VTX_to_BP);
        cout << "vertex   " << vertex << endl;
        cout << "BP   " << trkl_VTX_to_BP.get_point_ext() << endl;
        
        // Tracklet BP -> L1
        tracklet trkl_BP_to_L1;
        //trkl_BP_to_L1.set_theta(trkl_VTX_to_BP.get_theta());
        trkl_BP_to_L1.set_theta(trkl_VTX_to_BP.get_theta()+ trkl_VTX_to_BP.multiple_scattering(beam_pipe_Z, beam_pipe_X0, beam_pipe_thickness));  
        //trkl_BP_to_L1.set_phi(trkl_VTX_to_BP.get_phi());  
        trkl_BP_to_L1.set_phi(trkl_VTX_to_BP.get_phi()+ trkl_VTX_to_BP.multiple_scattering(beam_pipe_Z, beam_pipe_X0, beam_pipe_thickness));    
        trkl_BP_to_L1.set_point_int(trkl_VTX_to_BP.get_point_ext().extend_segment(trkl_VTX_to_BP.get_theta(), trkl_VTX_to_BP.get_phi(), beam_pipe_radius + beam_pipe_thickness));
        trkl_BP_to_L1.set_point_ext(trkl_BP_to_L1.find_intersection(layer1_radius));
        points_L1.push_back(trkl_BP_to_L1.get_point_ext());
        trkl_BP_L1.push_back(trkl_BP_to_L1);  
        
        // Tracklet L1 -> L2
        tracklet trkl_L1_to_L2;
        //trkl_L1_to_L2.set_theta(trkl_VTX_to_BP.get_theta());  
        //trkl_L1_to_L2.set_phi(trkl_VTX_to_BP.get_phi());   
        trkl_L1_to_L2.set_theta(trkl_BP_to_L1.get_theta()+ trkl_BP_to_L1.multiple_scattering(layer1_Z, layer1_X0, layer1_thickness));
        trkl_L1_to_L2.set_phi(trkl_BP_to_L1.get_phi()+ trkl_BP_to_L1.multiple_scattering(layer1_Z, layer1_X0, layer1_thickness));
        trkl_L1_to_L2.set_point_int(trkl_BP_to_L1.get_point_ext().extend_segment(trkl_VTX_to_BP.get_theta(), trkl_VTX_to_BP.get_phi(), layer1_radius + layer1_thickness));
        trkl_L1_to_L2.set_point_ext(trkl_L1_to_L2.find_intersection(layer2_radius));
        points_L2.push_back(trkl_L1_to_L2.get_point_ext());
        trkl_L1_L2.push_back(trkl_L1_to_L2);  
        
        cout << "------------" << endl;
        cout<<"THETA "<< trkl_VTX_to_BP.get_theta() << endl;
        cout<<"THETA dopo BP:  "<< trkl_BP_to_L1.get_theta() << endl;
        cout<<"THETA dopo Layer1:  "<< trkl_L1_to_L2.get_theta() << endl;
        cout << "------------" << endl;
                

    }

}

void evento::smearing(){

    //smearing di 30um applicato sull'arco di circonferenza --> sta cambiando il phi del punto
    //30um = 0.03mm

    for(int i = 0; i < points_L1.size(); i++){

        points_L1[i].set_cilindrical();
        double smearing = gRandom->Gaus(0, 0.03);
        points_L1[i].set_phi(points_L1[i].get_phi() + smearing / points_L1[i].get_R());
        points_L1[i].update_coordinates();

    }

    for(int i = 0; i < points_L2.size(); i++){

        points_L2[i].set_cilindrical();
        double smearing = gRandom->Gaus(0, 0.03);
        points_L2[i].set_phi(points_L2[i].get_phi() + smearing / points_L2[i].get_R());
        points_L2[i].update_coordinates();

    }


}


std::ostream &operator<<(std::ostream &output, const evento & ev) {
    output << "x del vertice: " << ev.vertex.get_x() << endl;
    output << "y del vertice: " << ev.vertex.get_y() << endl;
    output << "z del vertice: " << ev.vertex.get_z() << endl;
    output << "multiplicity: " << ev.multiplicity << endl;
    return output;
}