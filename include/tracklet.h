#ifndef TRACKLET_H
#define TRACKLET_H
#include <iostream>
#include "TObject.h"
#include "point.h"
using namespace std;

class tracklet : public TObject {

    public:
    tracklet();          
    virtual ~tracklet(); 
    //void print_tracklet();
    void generate_theta();
    void generate_phi();
    void generate_eta();
    double get_theta() const {return theta;}
    double get_phi() const {return phi;}
    double get_eta() const {return eta;}
    auto get_point_int() const {return point_int;}
    auto get_point_ext() const {return point_ext;}
    void set_points(point, point);
    void set_point_ext(point);
    void find_beampipe_intersection();
    //double multiple_scattering();

    

    friend std::ostream &operator<<(std::ostream &output, const tracklet &tracklet);


    private:
        point point_int;
        point point_ext;
        double theta;
        double eta;
        double phi;  

    ClassDef(tracklet, 1);
};

#endif