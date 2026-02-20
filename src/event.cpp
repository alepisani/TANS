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
#include <iostream>
#include "TF1.h"
#include <sstream>
#include <iomanip>

using namespace std;

event::event(int m, point p, point vtx):TObject(), multiplicity(m), pnt(p), vertex(vtx){}

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

    vertex = vtx;

}

void event::RunFullSimulation() {
    
    gErrorIgnoreLevel = kError;
    reconstruction reco;

    TFile* hfile = new TFile("../data/hist_sim.root", "RECREATE");
    TTree* tree = new TTree("Tree", "Tree simulazione");
    TH1D* hResidui = new TH1D("hResidui", "Residui; (z_{gen} - z_{reco}) [#mum]; Conteggi", 500, -1000, 1000);
    TH1D* hGen = new TH1D("hGen", "Eventi Generati; Molteplicita'; Conteggi", 50, 0, 50);
    TH1D* hReco = new TH1D("hReco", "Eventi Ricostruiti; Molteplicita'; Conteggi", 50, 0, 50);
    TH1D* hGenVsZ = new TH1D("hGenvsZ", "Eventi generati vs Z; Conteggi", 100, -beam_pipe_lenght/2., beam_pipe_lenght/2.);
    TH1D* hRecoVsZ = new TH1D("hRecovsZ", "Eventi ricostruiti vs Z; Conteggi", 100, -beam_pipe_lenght/2., beam_pipe_lenght/2.);
    TH2D* hResVsMult = new TH2D("hResVsMult", "Residui vs Mult; Molteplicita'; Residuo [#mum]", 50, 0, 50, 100, -2000, 2000);
    TH2D* hResVsZtrue = new TH2D("hResVsZtrue", "Residui vs Ztrue; z [mm]; Residuo [#mum]", 10, -beam_pipe_lenght/2., beam_pipe_lenght/2., 100, -2000, 2000);
    TH2D* hResVsZreco = new TH2D("hResVsZreco", "Residui vs Zreco; z [mm]; Residuo [#mum]", 10, -beam_pipe_lenght/2., beam_pipe_lenght/2., 100, -2000, 2000);
    TH1D* hist_z_vtx = new TH1D("hist_reco", "Ricostruzione Vertice Z; z_{cand} [mm]; Conteggi", bin_zvtx, -beam_pipe_lenght/2., beam_pipe_lenght/2.);
    hist_z_vtx->SetDirectory(0);

    //buffer usefull for the branch
    double vtx[3]; // X, Y, Z

    TClonesArray* hitsBP = new TClonesArray("point", 200);
    TClonesArray* particle_VTX_BP = new TClonesArray("particle", 200);
    TClonesArray* particle_BP_L1 = new TClonesArray("particle", 200);
    TClonesArray* particle_L1_L2 = new TClonesArray("particle", 200);
    TClonesArray* hitsL1 = new TClonesArray("point", 200);
    TClonesArray* hitsL2 = new TClonesArray("point", 200);

    tree->Branch("Vtx", vtx, "vtxX/D:vtxY/D:vtxZ/D");
    tree->Branch("Mult", &multiplicity, "mult/I");
    tree->Branch("HitsBP", &hitsBP);
    tree->Branch("Particle_VTX_BP", &particle_VTX_BP);
    tree->Branch("Particle_BP_L1", &particle_BP_L1);
    tree->Branch("Particle_L1_L2", &particle_L1_L2);
    tree->Branch("HitsL1", &hitsL1);
    tree->Branch("HitsL2", &hitsL2);

    TH1D* hist_eta = nullptr;
    TH1I* hist_mult = nullptr;
    TFile* hist_kinem = nullptr;
    if(get_data_from_kinem) {
        hist_kinem = TFile::Open("../data/kinem.root");
        if(hist_kinem && !hist_kinem->IsZombie()) {
            hist_eta = (TH1D*)hist_kinem->Get("heta2");
            hist_mult = (TH1I*)hist_kinem->Get("hm");
            if(!hist_eta || !hist_mult) {
                std::cerr << "Warning: histograms not found in the given file\n";
            }
        } else {
            std::cerr << "Warning: could not open the given file\n";
        }
    }

    //-----------------------------Event_Simulation------------------------------

    for (int iEv = 0; iEv < nEvents; ++iEv) {
        
        point vertex;
        vertex.generate_VTX();
        this->set_vertex(vertex);
        this->setmultiplicity(hist_mult);
        
        vtx[0] = vertex.get_x();
        vtx[1] = vertex.get_y();
        vtx[2] = vertex.get_z();
        int mult = this->get_multiplicity();

        //for each particle create in the collision we evaluate its path
        for (int iPart = 0; iPart < mult; ++iPart) {
            
            particle prtl;
            prtl.set_point(vertex);
            prtl.generate_theta(hist_eta);  
            prtl.generate_phi();

            //trasport of the particle from the vertex to the beam pipe
            new((*particle_VTX_BP)[iPart]) particle(prtl);
            prtl.find_intersection(beam_pipe_radius);
            new((*hitsBP)[iPart]) point(prtl.get_point());  
            points_BP.push_back(prtl.get_point());
            if (multiple_scattering_on){
                prtl.multiple_scattering(beam_pipe_Z, beam_pipe_X0, beam_pipe_thickness);
            };

            //trasport of the particle from the beam pipe to layer1
            new((*particle_BP_L1)[iPart]) particle(prtl);
            prtl.find_intersection(layer1_radius);
            prtl.get_point().smearing();
            new((*hitsL1)[iPart]) point(prtl.get_point());
            points_L1.push_back(prtl.get_point());
            if (multiple_scattering_on){
                prtl.multiple_scattering(layer1_Z, layer1_X0, layer1_thickness);
            };

            //trasport of the particle from the layer1 to layer2
            new((*particle_L1_L2)[iPart]) particle(prtl);
            prtl.find_intersection(layer2_radius);
            prtl.get_point().smearing();
            new((*hitsL2)[iPart]) point(prtl.get_point());
            points_L2.push_back(prtl.get_point());

        }

        //add noisy point on the detectors

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
        
        //close the file where we take the distribution
        if(hist_kinem) {
            hist_kinem->Close();
            delete hist_kinem;
        }
        
        //-----------------------------Reconstruction------------------------------
        
        double z_reco = reco.reco_z(this, hResidui, hist_z_vtx);
        double z_true = vertex.get_z();
        double residuo = (z_true - z_reco)*1000;  //from mm to um
        
        //clear the vector for the next iteration
        points_BP.clear();
        points_L1.clear();
        points_L2.clear();
        hitsBP->Clear("C");
        particle_VTX_BP->Clear("C");
        particle_BP_L1->Clear("C");
        particle_L1_L2->Clear("C");
        hitsL1->Clear("C");
        hitsL2->Clear("C");
        
        //filling the histograms
        hGen->Fill(mult);  
        hGenVsZ->Fill(z_true);
        hResVsMult->Fill(mult, residuo);
        hResVsZtrue->Fill(z_true, residuo);
        hResVsZreco->Fill(z_reco, residuo);


        //SELEZIONE DA RIVEDEREEEEEEEEEEEEEEEEEEEEEEEEEEE!!!!!!!!!!!!!!!!!!!!!!!
        if (std::abs(z_reco - z_true) < 1.5) {
            hReco->Fill(mult);
            hRecoVsZ->Fill(z_true);
        }



        tree->Fill();


        //loading bar
        if (iEv % 100 == 0 || iEv == nEvents - 1) {
           printProgressBar(iEv + 1, nEvents, 30);
        }
        
    }
    
    cout << endl;

    
    hfile->cd();


    //---------------------------------------------------------fill all the distros we're intrensted in

    //plot efficienty vs multeplicity
    if (hReco->GetEntries() > 0 && hGen->GetEntries() > 0) {
        TGraphAsymmErrors* gEff = new TGraphAsymmErrors(hReco, hGen, "cp"); // cp --> errori binomiali
        gEff->SetTitle("Efficiency vs Multiplicity; Multeplicity; Efficiency");
        gEff->SetMarkerStyle(20);
        gEff->SetMarkerSize(1.0);
        gEff->SetMarkerColor(kAzure+2);
        gEff->SetLineColor(kAzure+2);
        gEff->Write("eff_vs_mult");
    }

    //plot efficienty vs reconstructed z_vertex 
    TGraphAsymmErrors* gEffVsZ = new TGraphAsymmErrors(hRecoVsZ, hGenVsZ, "cp");
    gEffVsZ->SetTitle("Efficiency vs Z_{true}; Z [mm]; Efficiecy");
    gEffVsZ->SetMarkerStyle(20);
    gEffVsZ->SetMarkerSize(1.0);
    gEffVsZ->SetMarkerColor(kAzure+2);
    gEffVsZ->SetLineColor(kAzure+2);
    gEffVsZ->Write("eff_vs_zreco");

    if (hResidui->GetEntries() > 0) {
        hResidui->Fit("gaus", "Q");
        TF1 *fitFunc = hResidui->GetFunction("gaus");
        if (fitFunc) {
            double mean     = fitFunc->GetParameter(1);
            double sigma    = fitFunc->GetParameter(2);
            double sigmaErr = fitFunc->GetParError(2);
            cout << endl;
            cout << "+---------------------------------------------" << endl;
            cout << "| Gaussin fit on the residual plot: " << endl;
            cout << "| Mean : " << mean  << " um" << endl;
            cout << "| Dev.Std : (" << sigma  << " ± " << sigmaErr << ") #mum" << endl;
            cout << "+---------------------------------------------" << endl;
        }
    }
    
    /**
     * slice the TH2D to find the plot we want: 
     *  - resolution vs multeplicity
     *  - resolution vs z_reconstructed
     */

    hResVsMult->FitSlicesY(nullptr, 0, -1, 0, "Q");
    TH1D *hRisoluzione = (TH1D*)gDirectory->Get("hResVsMult_2");
    hRisoluzione->SetTitle("Risoluzione vs Molteplicita'; Molteplicita'; #sigma [#mum]");    
    // remove the plot i'm not intrested in
    if(auto h = (TH1D*)gDirectory->Get("hResVsMult_0")) h->SetDirectory(0);
    if(auto h = (TH1D*)gDirectory->Get("hResVsMult_1")) h->SetDirectory(0);
    if(auto h = (TH1D*)gDirectory->Get("hResVsMult_chi2")) h->SetDirectory(0);

    hResVsZtrue->FitSlicesY(nullptr, 0, -1, 0, "Q");
    TH1D *hRisoluzione_Z = (TH1D*)gDirectory->Get("hResVsZtrue_2");
    hRisoluzione_Z->SetTitle("Risoluzione vs Z; z [mm]; #sigma [#mum]");
    // remove the plot i'm not intrested in
    if(auto h = (TH1D*)gDirectory->Get("hResVsZtrue_0")) h->SetDirectory(0);
    if(auto h = (TH1D*)gDirectory->Get("hResVsZtrue_1")) h->SetDirectory(0);
    if(auto h = (TH1D*)gDirectory->Get("hResVsZtrue_chi2")) h->SetDirectory(0);
    
    


    hfile->Write("", TObject::kOverwrite);
    hfile->Close();

    delete hitsBP;
    delete hitsL1;
    delete hitsL2;
    delete particle_VTX_BP;
    delete particle_BP_L1;
    delete particle_L1_L2;
    delete hfile;

    cout << endl; cout << endl;
    std::cout << "Simulation ended: all event has been processed succesfully. \nData available on /data/hist_sim.root" << std::endl;
    cout << endl; cout << endl;

}

void event::printProgressBar(int current, int total, int barWidth = 30) {
    float progress = static_cast<float>(current) / total;
    int pos = static_cast<int>(barWidth * progress);

    std::ostringstream oss;
    oss << "\r[";

    for (int i = 0; i < barWidth; ++i) {
        if (i < pos)
            oss << "█";
        else
            oss << " ";
    }

    oss << "] ";
    oss << std::setw(3) << int(progress * 100.0) << "% ";
    oss << "(events=" << current << "/" << total << ")";

    std::string output = oss.str();
    size_t terminal_width = 80;
    if (output.size() < terminal_width)
        output += std::string(terminal_width - output.size(), ' ');

    std::cout << output << std::flush;
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

