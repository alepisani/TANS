#ifndef POINT_H
#define POINT_H
#include <iostream>
using namespace std;
#include "TObject.h"

class point : public TObject {

    public:
    point();              // COSTRUTTORE
    virtual ~point();     // DISTRUTTORE (obbligatorio con ROOT)
    void set_point(double, double, double);
    //auto get_point() {return point;}
  
    friend std::ostream &operator<<(std::ostream &output, const point &point);

    private:
        double x;
        double y;
        double z;  

    ClassDef(point, 1);
};

#endif