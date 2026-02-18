#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H
#include <iostream>
using namespace std;
#include "TObject.h"


class reconstruction : public TObject {

    public:
    reconstruction() : TObject() {}
    
    double reco_z(event*, class TH1D*, class TH1D*);
    double running_window();

    std::vector <double> z_candidates;     //all z found in the tracklet combination
    std::vector <double> vertex_candidate; //z values of the moda

    ClassDef(reconstruction,1);

};



#endif