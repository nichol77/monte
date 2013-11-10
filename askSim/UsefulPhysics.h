///////////////////////////////////////////////////////////////////////////////
/////   UsefulPhysics.h
/////   Describes one neutrino event in the ice and contains all its
/////   secondary interactions
//////////////////////////////////////////////////////////////////////////////

#ifndef USEFULPHYSICS_H
#define USEFULPHYSICS_H

#include "TObject.h"
#include "TRandom3.h"
#include "TF1.h"

#include "AskConventions.h"


class UsefulPhysics : public TObject {
public:
  UsefulPhysics();
  UsefulPhysics(double sigmaFactor, double massNucleus);
  ~UsefulPhysics();
  
  //Instance generator
  static UsefulPhysics*  Instance(double sigmaFactor=1, double massNucleus=1.66E-27);

  static double Gety();
  static AskCons::InteractionType_t GetCurrent();
  static AskCons::ParticleType_t GetNuFlavor();
  void FillGzkFlux();
  double PickGzkEnergy();
  double PickGzkIronEnergy();


 protected:
   static UsefulPhysics *fgInstance; 
   

private: 



    ClassDef(UsefulPhysics,1);
};


#endif //USEFULPHYSICS_H
