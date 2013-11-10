///////////////////////////////////////////////////////////////////////////////
/////   NeutrinoEvent.h
/////   Describes one neutrino event in the ice and contains all its
/////   secondary interactions
//////////////////////////////////////////////////////////////////////////////

#ifndef NEUTRINOEVENT_H
#define NEUTRINOEVENT_H

#include "TObject.h"
#include "TClonesArray.h"

#include "AskConventions.h"
#include "InteractionInfo.h"
#include "TVector3.h"


class NeutrinoEvent : public TObject {
public:
    NeutrinoEvent();
    NeutrinoEvent(Long64_t runId, Long64_t nuId, 
		  TVector3 surfPos,
		  Double_t energy, AskCons::ParticleType_t flavour,
		  AskCons::MaterialType_t matType,
		  TVector3 intPos, TVector3 nuDir,
		  AskCons::InteractionType_t intType,
		  Double_t elasticity);
    ~NeutrinoEvent();


    void AddInteraction(Double_t intEnergy, Double_t afterLength,
			AskCons::MaterialType_t matType,
			AskCons::InteractionType_t intType,
			Double_t showerEnergy);

    AskCons::ParticleType_t getCurrentParticleType(){ 
	InteractionInfo *previous = (InteractionInfo*)fInts->Last();
	return previous->getFinalParticleType();}
    
    Double_t getCurrentParticleEnergy(){ 
	InteractionInfo *previous = (InteractionInfo*)fInts->Last();	
	return previous->getFinalParticleEnergy();}

    TClonesArray *getInteractions()
       {return fInts;}

    TVector3 getSurfPos()
       { return TVector3(fSurfPos);}


    Long64_t fRunId; //MC Run Id
    Long64_t fNuId; //MC Neutrino Roll Number
    Double_t initialEnergy; //Neutrino energy 
    AskCons::ParticleType_t nuFlavour; //Neutrino flavour
    Int_t fNumInts; //Number of interactions
    Double_t fSurfPos[3]; // Position on surface (in m)
    TClonesArray *fInts; // Interaction info
    static TClonesArray *fgInts;

    ClassDef(NeutrinoEvent,4);
};


#endif //NEUTRINOEVENT_H
