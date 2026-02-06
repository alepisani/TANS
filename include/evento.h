#ifndef EVENTO_H
#define EVENTO_H

#include "TObject.h"
#include "./tracklet.h"
#include "./point.h"

class evento : public TObject {

 public:
    // inline default constructor so ROOT/cling finds it when compiling only headers
    evento();
   

    //implementazione in riga di funzioni

    void setmultiplicity();
    int getmultiplicity() const {return multiplicity;};
    void display_event();
    void generate_vertex();
    void event();
    void smearing();

    friend std::ostream &operator<<(std::ostream &output, const evento &ev);

    std::vector<tracklet> trkl_VTX_BP;    //will contain all tracklet from vtx to bp
    std::vector<tracklet> trkl_BP_L1;    //will contain all tracklet from bp to l1
    std::vector<tracklet> trkl_L1_L2;    //will contain all tracklet from l1 to l2
    std::vector<point> points_BP;        //will contain all points on bp
    std::vector<point> points_L1;        //will contain all points on l1
    std::vector<point> points_L2;        //will contain all points on l2
     

  
 private:
    int multiplicity;
    point vertex;

    ClassDef(evento,1); //classe evento per root
};





#endif
