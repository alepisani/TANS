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
    double running_window();
    
    point get_HitL1(int index) {return HitL1[index];}
    point get_HitL2(int index) {return HitL2[index];}
    std::vector<double> get_zcand() {return z_cand;}

    //setter e getter dei data member
    void reset();

    void analysis();

    private:
    std::vector<point> HitL1;
    std::vector<point> HitL2;
    std::vector<double> z_cand;



    ClassDef(reconstruction,1);

};



#endif