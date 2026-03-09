#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H
#include "point.h"
#include <iostream>
using namespace std;
#include "TObject.h"


class reconstruction : public TObject {

    public:
    reconstruction() : TObject() {}
    
    void reco();

    //read and fill vector
    void read_TTree();
    
    //fill zcandidates in a hist
    void fill_zcand_hist(TH1D* zcand_hist);
    void running_window(TH1D* zcand_hist);
    //zreco e fill histogramma
    double z_reco();
    
    point get_HitL1(int index) {return HitL1[index];}
    point get_HitL2(int index) {return HitL2[index];}
    void set_HitL1(point&, int);
    void set_HitL2(point&, int);

    //setter e getter dei data member
    void reset();

    private:
    std::vector<point> HitL1;
    std::vector<point> HitL2;



    ClassDef(reconstruction,1);

};



#endif