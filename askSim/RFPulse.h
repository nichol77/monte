///////////////////////////////////////////////////////////////////////////////
/////   RFPulse.h
/////   Describes the direction and characteristics of an RFPulse at a given
/////   point
//////////////////////////////////////////////////////////////////////////////

#ifndef RFPULSE_H
#define RFPULSE_H

#include "TObject.h"
#include "TClonesArray.h"
#include "TVector3.h"


class RFPulse : public TObject {
public:
    RFPulse();
    RFPulse(TVector3 location, TVector3 direction);
    ~RFPulse();

 private:
    TVector3 fLocation;
    TVector3 fDirection;

    //Need something that describes the frequency and phase content here
    //Will probably just be a couple of arrays.
    

    ClassDef(RFPulse,1);
};


#endif //RFPULSE_H
