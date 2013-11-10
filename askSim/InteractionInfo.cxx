///////////////////////////////////////////////////////////////////////////////
/////   InteractionInfo.cxx
/////   Contains information about the various interactions (primary and
/////   secondary)
///////////////////////////////////////////////////////////////////////////////
#include "InteractionInfo.h"

#include <iostream>

ClassImp(InteractionInfo)

InteractionInfo::InteractionInfo()
{
//Default (zero) constructor 
  
}

InteractionInfo::InteractionInfo(AskCons::ParticleType_t ipType,
				 AskCons::MaterialType_t matType,
				 TVector3 &iLoc,
				 Double_t theTime,
				 TVector3 &iDir,
				 Double_t energy,
				 AskCons::InteractionType_t intKind,
				 Double_t fracEM, Double_t fracHad, 
				 Double_t fracFinal,
				 AskCons::ParticleType_t fpType)
    :initialParticleType(ipType),
     intMaterial(matType),
     intTime(theTime),
     intEnergy(energy),
     intType(intKind),
     fracEMShower(fracEM),
     fracHadShower(fracHad),
     fracFinalParticle(fracFinal),
     finalParticleType(fpType)
{
   intDir[0]=iDir.X();
   intDir[1]=iDir.Y();
   intDir[2]=iDir.Z();
   intLocation[0]=iLoc.X();
   intLocation[1]=iLoc.Y();
   intLocation[2]=iLoc.Z();
   
//     std::cout << "II:\t" << finalParticleType << "\t" << fracEMShower
// 	      << "\t" << fracHadShower << "\t" << fracFinalParticle << std::endl
// 	;
// Assignment Constructor
}
