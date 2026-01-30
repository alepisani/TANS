#ifndef TRACCIA_H
#define TRACCIA_H

#include "TObject.h"
#include "TRandom3.h"

class traccia : public TObject { 

    public:
        traccia() : TObject(), eta(0), theta(0), phi(0), rand() {rand.SetSeed(0);}

        void setTheta();
        void setPhi();
        void setEta();
        double getTheta() const {return theta;};
        double getPhi() const {return phi;};
        double getEta() const {return eta;};

        friend std::ostream &operator<<(std::ostream &output, const traccia &track);

    private:
        double eta;
        double theta;
        double phi;
        TRandom3 rand;
        //manca l'impulso?

    ClassDef(traccia, 1);
};



#endif