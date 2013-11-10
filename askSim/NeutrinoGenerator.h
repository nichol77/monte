///////////////////////////////////////////////////////////////////////////////
/////   NeutrinoGenerator.h
/////   Describes one neutrino event in the ice and contains all its
/////   secondary interactions
//////////////////////////////////////////////////////////////////////////////

#ifndef NEUTRINOGENERATOR_H
#define NEUTRINOGENERATOR_H

#include "TObject.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TH2.h"
#include "TMath.h"

#include "AskConventions.h"
#include "NeutrinoEvent.h"
#include "NeuGenRunSummary.h"
#include "UsefulPhysics.h"
#include "WorldModel.h"


class NeutrinoGenerator : public TObject {
public:
  NeutrinoGenerator(Long64_t runId,  AskCons::WhichFlux_t whichFlux, Double_t maxDistToIce=0, Int_t doZCut=0,Double_t zCut=0);
    NeutrinoGenerator();
    ~NeutrinoGenerator();

    int getSurfacePointAndDirection(TVector3 &surfPos, TVector3 &intDir);
    NeutrinoEvent *getNextNeutrino();
    NeuGenRunSummary *getRunSummary();


    //Below should really be private... but for now
    inline Double_t pickRandomCosThetaOnSphere() 
       {return -1+2*gRandom->Rndm();}
    inline Double_t pickRandomCosThetaHalfPi()
       { return gRandom->Rndm();}
    inline Double_t pickRandomPhiTwoPi()
       { return TMath::TwoPi()*gRandom->Rndm();}
    
    void pickRandomThetaPhiOnSphere(Double_t &theta, Double_t &phi);  
    void pickRandomDowngoingDirection(Double_t &theta, Double_t &phi);
    void pickRandomThetaPhiOnSphere(Double_t &theta, Double_t &phi,Double_t &costheta);
    void pickRandomDowngoingDirection(Double_t &theta, Double_t &phi,Double_t &costheta);
    void pickPointAndDirectionOnVolumeInGlobalCoords(TVector3 &intPos, TVector3 &intDir);
    int getSurfacePointAndDirectionFromVolumePoint(TVector3 intPos, TVector3 &intDir, TVector3 &surfPos);

    void pickPointAndDirectionOnEarth(TVector3 &intPos, TVector3 &intDir);
    int setPointAndDirectionOnEarth();
    int getIceThicknessesAndLengths(Double_t &iceStart, Double_t &iceEnd);

    void setRandomSeed(UInt_t seed);
    void setFixedEnergy(Double_t fixedEnergy) {
       //Takes in eV sets it as GeV
       fFixedEnergy=fixedEnergy/1e9;
       fWhichFlux=AskCons::kMonoEnergetic;
    }


private:
    
    //Essentially temporary variables
    TVector3 fSurfPos;
    TVector3 fIntDir;
    TVector3 fIntPoint;
    TVector3 fOppPoint;
    TVector3 fRadialDir;
    
    Double_t fStartRadius;
    Double_t fEndRadius;
    Double_t fIntRadius;
    Double_t fSurfCosTheta;
    Double_t fChordLength;
    Double_t fStartChordLength;
    Double_t fEndChordLength;
    Double_t fStartThickness;
    Double_t fEndThickness;



    Int_t fDoZCut;    
    Double_t fZCut;
    Long64_t fRunId;
    Double_t fVolRad;
    Long64_t fRollCounter;
    Long64_t fIntCounter;
    Double_t fMaxDistToIce;
    UInt_t fRandomSeed;
    AskCons::WhichFlux_t fWhichFlux;
    Double_t fFixedEnergy;
    


    UsefulPhysics *fUseful;
    WorldModel *fWorldModel;
    TH2F *fHistLatLonSurf;
    
    ClassDef(NeutrinoGenerator,1);
};


#endif //NEUTRINOGENERATOR_H
