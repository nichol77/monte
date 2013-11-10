///////////////////////////////////////////////////////////////////////////////
/////   InteractionInfo.h
/////   Simple description of a stochastic energy loss interaction
//////////////////////////////////////////////////////////////////////////////

#ifndef INTERACTIONINFO_H
#define INTERACTIONINFO_H

#include "TObject.h"
#include "TVector3.h"

#include "AskConventions.h"


class InteractionInfo : public TObject {
public:
    InteractionInfo();
    InteractionInfo(AskCons::ParticleType_t ipType,
		    AskCons::MaterialType_t matType,
		    TVector3 &iLoc,
		    Double_t theTime,
		    TVector3 &iDir,
		    Double_t energy,
		    AskCons::InteractionType_t intKind,
		    Double_t fracEM, Double_t fracHad, Double_t fracFinal,
		    AskCons::ParticleType_t fpType);

    TVector3 getDirection()
       {return TVector3(intDir);}
    TVector3 getLocation() 
	{return TVector3(intLocation);} 
    AskCons::ParticleType_t getFinalParticleType()
	{ return finalParticleType;}
    Double_t getIntTime()
	{return intTime;}
    Double_t getIntEnergy() 
	{return intEnergy;}
    Double_t getFinalParticleEnergy() 
	{return intEnergy*fracFinalParticle;}


    AskCons::ParticleType_t initialParticleType; //Particle that interacted

    AskCons::MaterialType_t intMaterial; //Material the interaction occured in
    //    TVector3 intLocation;
    Double_t intLocation[3]; // Interaction point
    Double_t intTime; // t=0 is the time of neutrino interaction
    //    TVector3 intDir;
    Double_t intDir[3]; // direction of travel at interaction point
    Double_t intEnergy; // energy of interacting particle
    AskCons::InteractionType_t intType; // type of interaction

    Double_t fracEMShower; // fraction of energy as em shower
    Double_t fracHadShower; // fraction of energy as hadronic shower
    Double_t fracFinalParticle;  // The three fracs should add to unity   
    
    AskCons::ParticleType_t finalParticleType; //Final particle type if one was created

    ClassDef(InteractionInfo,1);
};


#endif //INTERACTIONINFO_H
