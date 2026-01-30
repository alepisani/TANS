#ifndef EVENTO_H
#define EVENTO_H

#include "TObject.h"
#include "TRandom3.h"

class evento : public TObject {

 public:
    // inline default constructor so ROOT/cling finds it when compiling only headers
    evento() : TObject(), x(0), y(0), z(0), molteplicita(0), rand() {rand.SetSeed(0);}
    //evento(double z, double x, double y, int molteplicita);

    //evento(const evento& source); //copy constructor
    //virtual ~evento();
    //evento& operator=(const evento& source); //assignment operator

    //implementazione in riga di funzioni

    void setX() {x=rand.Gaus(0,0.1);};
    void setY() {y=rand.Gaus(0,0.1);};
    void setZ() {z=rand.Uniform(0,53);};

    void setMolteplicita() ;

    double getX() const {return x;};
    double getY() const {return y;};            
    double getZ() const {return z;};
    int getMolteplicita() const {return molteplicita;};

    friend std::ostream &operator<<(std::ostream &output, const evento &ev);

  
 private:
    double x;
    double y;
    double z;
    int molteplicita;
    TRandom3 rand; //generatore di numeri casuali

    ClassDef(evento,1); //classe evento per root
};





#endif
