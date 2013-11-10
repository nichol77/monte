///////////////////////////////////////////////////////////////////////////////
/////   NeuGenRunSummary.h
/////   Describes one neutrino generation run
//////////////////////////////////////////////////////////////////////////////

#ifndef NEUGENRUNSUMMARY_H
#define NEUGENRUNSUMMARY_H

#include "TObject.h"
#include "TClonesArray.h"

#include "AskConventions.h"

class NeuGenRunSummary : public TObject {
public:
    NeuGenRunSummary();
    NeuGenRunSummary(Long64_t runId, Long64_t numRolls,
		     Long64_t numInteractions, Double_t volumeRadius,
		     Double_t maxDistToIce, AskCons::WhichFlux_t whichFlux, 
		     Double_t fixedEnergy,
		     Double_t crossSectionScale, UInt_t randomSeed,
		     Int_t doZCut, Double_t zCut);
		     
    ~NeuGenRunSummary();

    
    Long64_t fRunId; //MC run id
    Long64_t fNumRolls; //Number of neutrinos tried
    Long64_t fNumInts; //number that interacted in/near ice
    Double_t fVolRad; //radius of volume over which flux was generated
    Double_t fMaxDistToIce; //maximum distance from ice that an interaction was saved
    AskCons::WhichFlux_t fWhichFlux; //holder for a flux type flag
    Double_t fFixedEnergy;
    Double_t fCrossSectionScale; //holder for a cross section scale value
    UInt_t fRandomSeed;
    Int_t fDoZCut;
    Double_t fZCut;

    ClassDef(NeuGenRunSummary,4);
};


#endif //NEUGENRUNSUMMARY_H
