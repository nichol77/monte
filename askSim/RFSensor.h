///////////////////////////////////////////////////////////////////////////////
/////   RFSensor.h
/////   Describes the location and characteristics of an RFSensor within
/////   an RFStation
//////////////////////////////////////////////////////////////////////////////

#ifndef RFSENSOR_H
#define RFSENSOR_H

#include "TObject.h"
#include "TClonesArray.h"

#include "TVector3.h"

class RFPulse;

class RFSensor : public TObject {
public:
    RFSensor();
    RFSensor(TVector3 relPos);
    ~RFSensor();

    Bool_t testTrigger(RFPulse &pulsey);
    TVector3 getRelPos() {return fRelPos;}


 private:
    TVector3 fRelPos;

/*     Int_t fUseSimple */
/*     Double_t fBandwidth; */
/*     Double_t fCentreFrequency; //Note the spelling */
    

    ClassDef(RFSensor,1);
};


#endif //RFSENSOR_H
