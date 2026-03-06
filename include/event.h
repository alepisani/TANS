#ifndef EVENT_H
#define EVENT_H

#include "TObject.h"
#include "./particle.h"
#include "./point.h"
#include <TClonesArray.h>

class event : public TObject {

public:
   event():TObject(), multiplicity(0), prtl(){}
   event(int, particle);
   
   int get_multiplicity() const {return multiplicity;}
   const point& get_vertex() const {return vertex;}
   
   void setmultiplicity(TH1I* hist_mult);
   void set_vertex(const point& vtx);

   void single_event(TClonesArray* vtx, TClonesArray* hitL1, TClonesArray* hitL2, TH1I* hist_mult, TH1D* hist_eta);
   void noise(TClonesArray* hitL1, TClonesArray* hitL2);

   friend std::ostream &operator<<(std::ostream &output, const event &ev);

   std::vector<point> points_BP;        //will contain all points on bp
   std::vector<point> points_L1;        //will contain all points on l1
   std::vector<point> points_L2;        //will contain all points on l2
     

  
private:
   int multiplicity;
   point vertex;
   particle prtl;
   int counterVTX = 0;
   int counterL1 = 0;
   int counterL2 = 0;

   ClassDef(event,1);
};





#endif
