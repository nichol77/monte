///////////////////////////////////////////////////////////////////////////////
/////   NeuGenRunSummary.h
/////   Describes one neutrino generation run
//////////////////////////////////////////////////////////////////////////////


#include "NeuGenRunSummary.h"


NeuGenRunSummary::NeuGenRunSummary()
{

}

NeuGenRunSummary::NeuGenRunSummary(Long64_t runId, Long64_t numRolls,
				   Long64_t numInteractions, Double_t volumeRadius,
				   Double_t maxDistToIce, AskCons::WhichFlux_t whichFlux,
				   Double_t fixedEnergy,
				   Double_t crossSectionScale,UInt_t randomSeed,Int_t doZCut,
				   Double_t zCut)
   :fRunId(runId),fNumRolls(numRolls),fNumInts(numInteractions),fVolRad(volumeRadius),
    fMaxDistToIce(maxDistToIce),fWhichFlux(whichFlux),fFixedEnergy(fixedEnergy),fCrossSectionScale(crossSectionScale),fRandomSeed(randomSeed),fDoZCut(doZCut),fZCut(zCut)
{


}

NeuGenRunSummary::~NeuGenRunSummary()
{

}
