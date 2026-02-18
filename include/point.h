#ifndef POINT_H
#define POINT_H
#include <iostream>
using namespace std;
#include "TObject.h"

class point : public TObject {

    public:
    point():TObject(), x(0.), y(0.), z(0.), R(0.), phi(0.){}
    point(double, double, double);
    
    double get_x() const {return x;}
    double get_y() const {return y;}
    double get_z() const {return z;}
    double get_R() const {return R;}
    double get_phi() const {return phi;}

    void set_phi(double);
    void set_R(double);
    void set_z(double);
    void set_point(double, double, double);

    void generate_VTX();
    void smearing();
  
    friend std::ostream &operator<<(std::ostream &output, const point &point);

    private:
        //Cartesian coordinates (x,y,z)
        double x;
        double y;
        double z;  
        
        //Polar coordinates (R, phi)
        double R;
        double phi;

    ClassDef(point, 1);
};

#endif