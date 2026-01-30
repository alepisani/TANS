#ifndef TRACKLET_H
#define TRACKLET_H
#include <iostream>
using namespace std;
#include "TObject.h"
#include "point.h"

class tracklet : public TObject {

    public:
    tracklet();          
    virtual ~tracklet(); 
    //void print_tracklet();
    void generate_theta();
    void generate_phi();
    void generate_eta();
    void generate_vertex();
    //double multiple_scattering();

    

    friend std::ostream &operator<<(std::ostream &output, const tracklet &tracklet);

    private:
        point point;
        double theta;
        double eta;
        double phi;  

    ClassDef(tracklet, 1);
};

#endif