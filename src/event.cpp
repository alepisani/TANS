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
#include "TF1.h"

using namespace std;

event::event(int m, point p, point vtx):TObject(), multiplicity(m), pnt(p), vertex(vtx){}

void event::setmultiplicity(TH1I* hist_mult) {

    if(!distrib_assegnata){

        multiplicity = static_cast<int>(gRandom->Uniform(1, 50));

    }

    if(distrib_assegnata){
        if(hist_mult) {
            multiplicity = hist_mult->GetRandom();
        } else {
            std::cerr << "Warning: hist_mult is nullptr, using random multiplicity\n";
        }
    }

}

void event::set_vertex(point vtx){

    vertex = vtx;

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


// QUELLA CHE SEGUE E' LA VERSIONE CON L'ISTOGRAMMA DEI RESIDUI!!!!
void event::RunFullSimulation() {
    
    reconstruction reco;

    TFile* hfile = new TFile("../data/hist_sim.root", "RECREATE");
    TTree* tree = new TTree("Tree", "Tree simulazione");
    TH1D* hResidui = new TH1D("hResidui", "Residui; (z_{gen} - z_{reco}) [#mum]; Conteggi", 500, -10*1000, 10*1000);
    TH1D* hGen = new TH1D("hGen", "Eventi Generati; Molteplicita'; Conteggi", 30, 0, 50);
    TH1D* hReco = new TH1D("hReco", "Eventi Ricostruiti; Molteplicita'; Conteggi", 30, 0, 50);
    TH1D* hGenVsZ = new TH1D("hGenvsZ", "Eventi generati vs Z; Conteggi", 30, -beam_pipe_lenght/2., beam_pipe_lenght/2.);
    TH1D* hRecoVsZ = new TH1D("hRecovsZ", "Eventi ricostruiti vs Z; Conteggi", 30, -beam_pipe_lenght/2., beam_pipe_lenght/2.);
    TH2D* h2D = new TH2D("h2D", "Residui vs Mult; Molteplicita'; Residuo [#mum]", 10, 0, 50, 200, -500, 500);
    TH2D* h2D_Z = new TH2D("h2D_Z", "Residui vs Z; z [mm]; Residuo [#mum]", 30, -beam_pipe_lenght/2., beam_pipe_lenght/2., 200, -500, 500);

    // Buffer per i dati scalari (Vertice e Molteplicità)
    // Usiamo una struct o variabili singole per la Leaf List
    double vtx[3]; // X, Y, Z

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

    // Load eta histogram once at the beginning
    TH1D* hist_eta = nullptr;
    TH1I* hist_mult = nullptr;
    TFile* hist_kinem = nullptr;
    if(distrib_assegnata) {
        hist_kinem = TFile::Open("../data/kinem.root");
        if(hist_kinem && !hist_kinem->IsZombie()) {
            hist_eta = (TH1D*)hist_kinem->Get("heta2");
            hist_mult = (TH1I*)hist_kinem->Get("hm");
            if(!hist_eta || !hist_mult) {
                std::cerr << "Errore: istogrammi non trovato nel file\n";
            }
        } else {
            std::cerr << "Errore nell'aprire il file kinem.root\n";
        }
    }

    //event evSim;

    for (int iEv = 0; iEv < nEvents; ++iEv) {
        
        point vertex;
        vertex.generate_VTX();
        this->set_vertex(vertex);
        //cout << vertex << endl;
        this->setmultiplicity(hist_mult);
        
        vtx[0] = vertex.get_x();
        vtx[1] = vertex.get_y();
        vtx[2] = vertex.get_z();
        this->get_multiplicity();

        for (int iPart = 0; iPart < this->get_multiplicity(); ++iPart) {
            particle prtl;
            prtl.set_point(vertex);
            prtl.generate_theta(hist_eta);  // Pass histogram to avoid reopening file
            prtl.generate_phi();

            //trasporto da vtx a bp
            new((*particle_VTX_BP)[iPart]) particle(prtl);
            prtl.find_intersection(beam_pipe_radius);
            new((*hitsBP)[iPart]) point(prtl.get_point());  
            points_BP.push_back(prtl.get_point());
            if (multiple_scattering_on){
                prtl.multiple_scattering(beam_pipe_Z, beam_pipe_X0, beam_pipe_thickness);
            };

            //trasporto da bp a l1
            new((*particle_BP_L1)[iPart]) particle(prtl);
            prtl.find_intersection(layer1_radius);
            prtl.get_point().smearing();
            new((*hitsL1)[iPart]) point(prtl.get_point());
            points_L1.push_back(prtl.get_point());
            if (multiple_scattering_on){
                prtl.multiple_scattering(layer1_Z, layer1_X0, layer1_thickness);
            };

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
        
        //reco
        double z_reco = reco.reco_z(this, hResidui);
        double z_true = vertex.get_z();
        int m = this->get_multiplicity();
        double residuo = (z_true - z_reco)*1000;

        // Riempi sempre il denominatore
        hGen->Fill(m);  //praticamente aggiunge 1 al bin corrispondente alla molteplicità m, quindi conta quanti vertici genera per ogni molteplicità
        hGenVsZ->Fill(z_true);
        h2D->Fill(m, residuo);
        h2D_Z->Fill(z_true, residuo);
        // Stabilisci il criterio di "successo":
        // Se lo scarto è minore di una soglia (es. 1.5 mm), considero il vertice trovato
        if (std::abs(z_reco - z_true) < 1.5) {
            hReco->Fill(m);
            hRecoVsZ->Fill(z_true);
        }

        //devi pulire i vector di punti
        points_BP.clear();
        points_L1.clear();
        points_L2.clear();


        tree->Fill();

        hitsBP->Clear("C");
        particle_VTX_BP->Clear("C");
        particle_BP_L1->Clear("C");
        particle_L1_L2->Clear("C");
        hitsL1->Clear("C");
        hitsL2->Clear("C");

        if (iEv % 100 == 0) std::cout << "Event " << iEv << " has been simulated." << std::endl;
    }
    
    // Clean up eta histogram file
    if(hist_kinem) {
        hist_kinem->Close();
        delete hist_kinem;
    }
    
    hfile->cd();

    //MODIFICHE FATTE ORA    
    TGraphAsymmErrors* gEff = new TGraphAsymmErrors(hReco, hGen, "cp"); //"cp" indica che metodo usare per calcolare gli errori (binomiali in questo caso)
    // 2. Estetica del grafico (punti ed errori)
    gEff->SetTitle("Efficienza di Ricostruzione;Molteplicita';Efficienza");
    gEff->SetMarkerStyle(20);
    gEff->SetMarkerSize(1.0);
    gEff->SetMarkerColor(kAzure+2);
    gEff->SetLineColor(kAzure+2);
    // Se vuoi anche la curva continua (il fit) sopra:
    TF1 *f_logis = new TF1("f_logis", "[0] / (1.0 + exp(-[1]*(x-[2])))", 0, 50);
    f_logis->SetParameters(1, 0.5, 5);
    gEff->Fit(f_logis, "R");
    f_logis->SetLineColor(kRed);
    gEff->Write();

    TGraphAsymmErrors* gEffVsZ = new TGraphAsymmErrors(hRecoVsZ, hGenVsZ, "cp");
    gEffVsZ->SetTitle("Efficienza di Ricostruzione vs Z; z [mm]; Efficienza");
    gEffVsZ->SetMarkerStyle(20);
    gEffVsZ->SetMarkerSize(1.0);
    gEffVsZ->SetMarkerColor(kAzure+2);
    gEffVsZ->SetLineColor(kAzure+2);
    gEffVsZ->Write();

    hResidui->Fit("gaus");
    // Recupera la funzione di fit associata all'istogramma
    TF1 *fitFunc = hResidui->GetFunction("gaus");

    if (fitFunc) {
        double mean  = fitFunc->GetParameter(1); // Media (z_gen - z_reco)
        double sigma = fitFunc->GetParameter(2); // Questa è la tua risoluzione migliorata!
        double sigmaErr = fitFunc->GetParError(2); // Errore sulla sigma

        cout << "Media dei residui CALCOLATA DAL FIT: " << mean << " um" << endl; 
        cout << "Risoluzione del vertice (Sigma) CALCOLATA DAL FIT: " << sigma << " um" << endl;
        cout << "Errore sulla risoluzione: " << sigmaErr << " um" << endl;
    }
    
    cout << "Deviazione standard dei residui: " << hResidui->GetStdDev() << " μm" << std::endl;

    // Questo comando divide l'istogramma 2D in fette (slice) verticali (per bin di molteplicità)
    // e fa un fit gaussiano su ogni fetta automaticamente.
    h2D->FitSlicesY();

    // ROOT crea automaticamente degli istogrammi con i risultati dei fit (altezza, media, sigma, chi2).
    // L'istogramma con suffisso "_2" contiene le SIGMA di ogni fetta.
    TH1D *hRisoluzione = (TH1D*)gDirectory->Get("h2D_2");
    hRisoluzione->SetTitle("Risoluzione vs Molteplicita'; Molteplicita'; #sigma [um]");
    hRisoluzione->Draw();
    hRisoluzione->Write();

    h2D_Z->FitSlicesY();
    TH1D *hRisoluzione_Z = (TH1D*)gDirectory->Get("h2D_Z_2");
    hRisoluzione_Z->SetTitle("Risoluzione vs Z; z [mm]; #sigma [um]");
    hRisoluzione_Z->Write();

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

