#ifndef POINT_H
#define POINT_H
#include <iostream>
using namespace std;
#include "TObject.h"

class point : public TObject {

    public:
    point();              
    virtual ~point();     
    void set_point(double, double, double);
    //auto get_point() {return point;}
    double get_x() const {return x;}
    double get_y() const {return y;}
    double get_z() const {return z;}
    point extend_segment(double, double, double);
  
    friend std::ostream &operator<<(std::ostream &output, const point &point);

    private:
        double x;
        double y;
        double z;  

    ClassDef(point, 1);
};

#endif