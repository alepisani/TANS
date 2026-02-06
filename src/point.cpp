#include "../include/point.h"
#include "../include/const.h"
#include "TRandom3.h"
#include <cmath>
#include <iostream>
#include <cmath>

using namespace std;


ClassImp(point);

point::point(double a, double b, double c):TObject(),
    x(a), y(b), z(c) 
    {
        R = sqrt(x*x + y*y);
        phi = atan2(y,x);
    }

void point::set_point(double a, double b, double c){
    
    x = a;
    y = b;
    z = c;
    R = sqrt(x*x + y*y);
    phi = atan2(y,x);

}

void point::set_phi(double p){

    phi = p;
    x = R * cos(phi);
    y = R * sin(phi);

}

void point::generate_VTX(){
    
    double x = gRandom->Gaus(0,X_rms);
    double y = gRandom->Gaus(0,Y_rms);
    double z = gRandom->Gaus(0,Z_rms);
    
    this->set_point(x,y,z);


}

std::ostream &operator<<(std::ostream &output, const point &point){
    output << "point coordinate (mm): (x, y, z) = (" << point.x << ", " 
                                    << point.y << ", " 
                                    << point.z << ") (" 
                                    << point.R << ", "  
                                    << point.phi << ")" << endl;
    return output;
}
