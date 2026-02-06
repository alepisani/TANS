#include "../include/point.h"
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

std::ostream &operator<<(std::ostream &output, const point &point){
    output << "point coordinate (mm): (x, y, z) = (" << point.x << ", " 
                                    << point.y << ", " 
                                    << point.z << ") (" 
                                    << point.R << ", "  
                                    << point.phi << ")" << endl;
    return output;
}
