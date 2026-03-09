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
#include "TFile.h"
#include "TTree.h"

using namespace std;

ClassImp(reconstruction);

void reconstruction::reset(){

    HitL1.clear();
    HitL2.clear();

}

void reconstruction::read_TTree(){

    TFile hfile("../data/simulation.root");
    
    TTree *tree = (TTree*)hfile.Get("MC");
    
    TClonesArray *hitsL1 = nullptr;
    TClonesArray *hitsL2 = nullptr;
    tree->SetBranchAddress("HitsL1", &hitsL1);
    tree->SetBranchAddress("HitsL2", &hitsL2);

    for(int nEvent = 0; nEvent < tree->GetEntries() ; nEvent++){
        
        tree->GetEvent(nEvent);
        
        this->reset();
        
        for(int i = 0; i < hitsL1->GetEntries(); i++){

            point* p1 = (point*)hitsL1->At(i);
            HitL1.push_back(*p1);

        }

        for(int i = 0; i < hitsL2->GetEntries(); i++){

            point* p2 = (point*)hitsL2->At(i);
            HitL2.push_back(*p2);

        }

    }

    hfile.Close();

}



void reconstruction::reco(){

    this->read_TTree();
    
    
    
    /**
     * read the ttree
     * fill zcand_hist
     * find the maximum with next bins
     * running window 
     * 
     */

}