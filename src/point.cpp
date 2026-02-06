#include "../include/point.h"
#include <cmath>
#include <iostream>
#include <cmath>

using namespace std;


ClassImp(point);

// COSTRUTTORE
point::point() : x(0), y(0), z(0) {}

// DISTRUTTORE
point::~point() {}

void point::set_point(double a, double b, double c){
    x = a;
    y = b;
    z = c;
}

point point::extend_segment(double theta, double phi, double lenght){

    double c1 = sin(theta) * cos(phi);
    double c2 = sin(theta) * sin(phi);
    double c3 = cos(theta);
    double x0 = this->get_x();
    double y0 = this->get_y();
    double z0 = this->get_z();

    double delta = (x0 * c1 + y0 * c2) * (x0 * c1 + y0 * c2) - (c1 * c1 + c2 * c2) * (x0 * x0 + y0 * y0 - lenght * lenght);
    double t = (-(x0 * c1 + y0 * c2)+sqrt(delta))/(c1 * c1 + c2 * c2);

    double x = x0 + c1 * t;
    double y = y0 + c2 * t;
    double z = z0 + c3 * t;

    this->set_point(x, y, z);
    return *this;

}

std::ostream &operator<<(std::ostream &output, const point &point){
    output << "point coordinate (mm): (" << point.x << ", " 
                                    << point.y << ", " 
                                    << point.z << ")" << endl;
    return output;
}