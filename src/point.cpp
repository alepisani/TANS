#include "../include/point.h"
#include <iostream>
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


std::ostream &operator<<(std::ostream &output, const point &point){
    output << "point coordinate (mm): (" << point.x << ", " 
                                    << point.y << ", " 
                                    << point.z << ")" << endl;
    return output;
}