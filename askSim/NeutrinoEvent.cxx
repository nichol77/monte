///////////////////////////////////////////////////////////////////////////////
/////   NeutrinoEvent.cxx
/////   Contains information about the various interactions (primary and
/////   secondary)
///////////////////////////////////////////////////////////////////////////////
#include "NeutrinoEvent.h"

#include <iostream>

#include "PhysicalConstants.h"

ClassImp(NeutrinoEvent)

TClonesArray *NeutrinoEvent::fgInts = 0;



NeutrinoEvent::NeutrinoEvent()
    :initialEnergy(0),
     nuFlavour(AskCons::kNothing),fNumInts(0)
{
//Default (zero) constructor 
   if(!fgInts)
    fgInts = new TClonesArray("InteractionInfo", 10);
   fInts=fgInts;
}

NeutrinoEvent::NeutrinoEvent(Long64_t runId, Long64_t nuId, 
			     TVector3 surfPos, Double_t energy, 
			     AskCons::ParticleType_t flavour,		
			     AskCons::MaterialType_t matType,
			     TVector3 intPos, TVector3 nuDir, 
			     AskCons::InteractionType_t intType,
			     Double_t elasticity)
    :fRunId(runId), fNuId(nuId), initialEnergy(energy),
     nuFlavour(flavour),fNumInts(0)
{
//Assignment constructor
   if(!fgInts)
      fgInts = new TClonesArray("InteractionInfo",10);
   fInts=fgInts;
   fSurfPos[0]=surfPos.X();
   fSurfPos[1]=surfPos.Y();
   fSurfPos[2]=surfPos.Z();
    Double_t fEM=0,fHad=0,fPart=0;
    AskCons::ParticleType_t particle=flavour;
    switch(intType) {
	case AskCons::kNC:
	    fHad=elasticity;
	    fPart=1-elasticity;
	    break;
	case AskCons::kCC:
	    fHad=elasticity;
	    fPart=1-elasticity;
	    particle=AskCons::getCCPartner(particle);
	    if(abs(flavour)==AskCons::kElectron) {
		fEM=1-elasticity;
		fPart=0;
		particle=AskCons::kNothing;
	    }
	    break;
	default:
	    std::cerr << "Don't understand interaction type: "
		      << intType << std::endl;
	    break;
    }


    TClonesArray &interactions = *fInts;
    //    InteractionInfo *intInfo = 
    new(interactions[fNumInts++]) InteractionInfo(flavour,matType,intPos,0,nuDir,
						  energy,intType,
						  fEM,fHad,fPart,particle);
    //    std::cout << fInts->Last() << "\n";
    
}


NeutrinoEvent::~NeutrinoEvent()
{
//Default destructor 
    //May need to loop over array and delete things
   //std::cout << "NeutrinoEvent::~NeutrinoEvent" << std::endl;

   //   for(Int_t i=0;i<1+fInts->GetLast();i++){
   //      fInts->At(i)->Clear();
   //   }
   //   fInts->Clear();
   fInts->Clear("C");
   //   delete fInts;
}


void NeutrinoEvent::AddInteraction(Double_t intEnergy, Double_t afterLength,		
				   AskCons::MaterialType_t matType,
				   AskCons::InteractionType_t intType,
				   Double_t showerEnergy)
{
    
    Double_t fEM=0,fHad=0,fPart=0;
    InteractionInfo *previous=(InteractionInfo*)fInts->Last();
    AskCons::ParticleType_t intParticle=previous->getFinalParticleType();
    AskCons::ParticleType_t finalParticle=intParticle;
    TVector3 theDir=previous->getDirection();
    TVector3 intPos=previous->getLocation();
    theDir.SetMag(afterLength);
    intPos+=theDir;
    theDir.SetMag(1);
    Double_t theTime=previous->getIntTime()+(afterLength/c_SI);
    switch(intType) {
	case AskCons::kBrem:
	    fHad=0;
	    fEM=showerEnergy;
	    fPart=intEnergy-showerEnergy;
	    break;
	    
	case AskCons::kPair:
	    fHad=0;
	    fEM=showerEnergy;
	    fPart=intEnergy-showerEnergy;
	    break;
	    
	case AskCons::kPhoto:
	    fEM=0;
	    fHad=showerEnergy;
	    fPart=intEnergy-showerEnergy;
	    break;
	    
	case AskCons::kKnockOn:
	    fHad=0;
	    fEM=showerEnergy;
	    fPart=intEnergy-showerEnergy;
	    break;

	case AskCons::kCC:
	    fHad=showerEnergy;
	    fPart=intEnergy-showerEnergy;
	    finalParticle=AskCons::getCCPartner(intParticle);
	    break;

	case AskCons::kDecay:
	    fHad=showerEnergy;
	    fPart=intEnergy-showerEnergy;
	    finalParticle=AskCons::getCCPartner(intParticle);
	    break;
	default:
	    std::cerr << "Don't understand interaction type: "
		      << intType << std::endl;
	    break;
    }
    fEM/=intEnergy;
    fHad/=intEnergy;
    fPart/=intEnergy;
    
    //     std::cerr << "AddInteraction:\t" << intParticle << "\t" << intEnergy
    // 	      << "\t" << fEM << "\t" << fHad << "\t" << fPart << std::endl;
    TClonesArray &interactions = *fInts;
    new(interactions[fNumInts++]) InteractionInfo(intParticle,matType,
    						  intPos,theTime,
    						  theDir,
    						  intEnergy,intType,
    						  fEM,fHad,fPart,
    						  finalParticle);  
}

