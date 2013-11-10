///////////////////////////////////////////////////////////////////////////////
/////   UsefulPhysics.cxx
/////   Contains information about the various interactions (primary and
/////   secondary)
///////////////////////////////////////////////////////////////////////////////
#include "UsefulPhysics.h"

#include <iostream>
#include "interactions.h"

//#include "TF1.h"
//#include "TRandom.h"
#include "TMath.h"

using namespace std;

ClassImp(UsefulPhysics)
   //static thingies and functions that only exist here
UsefulPhysics*  UsefulPhysics::fgInstance = 0;
  



  namespace UPThings {
    double fSigmaFactor=1;
    double fMassNucleus=1.66E-27;
    
    double energy[12];
    double EdNdEdAdt[12]; //flux of incident neutrinos vs. energy
    double maxflux;
    double dNdEdAdt[12];
  }

UsefulPhysics::UsefulPhysics(double sigmaFactor, double massNucleus)
{
//Default (zero) constructor 
  UPThings::fSigmaFactor=sigmaFactor;
  UPThings::fMassNucleus=massNucleus;
}


UsefulPhysics::UsefulPhysics()
{
//Default (zero) constructor 

}

UsefulPhysics::~UsefulPhysics()
{
//Default Destructor
}


//______________________________________________________________________________
UsefulPhysics*  UsefulPhysics::Instance(double sigmaFactor, double massNucleus)
{
   //static function
  if(fgInstance)
    return fgInstance;

  fgInstance = new UsefulPhysics(sigmaFactor,massNucleus);
  return fgInstance;
  //   return (fgInstance) ? (AnitaEventCalibrator*) fgInstance : new AnitaEventCalibrator();
}


// The interaction
double UsefulPhysics::Gety() {
  // THIS IS A ROUGH PARAMETRIZATION OF PLOT 6 FROM 
  //  Ghandhi,Reno,Quigg,Sarcevic  hep-ph/9512364
  //  (the curves are not in their later article.)
  //  There is also a slow energy dependence.
  
  float rnd;
  float x = 0;
  const double R1=0.36787944;  // 1/e
  const double R2=0.63212056;  // 1-r1
  
  // generate according to Ghandi fig. 6 
  // adjust exponent until looks like the curve
  //  and has right mean.
  //  (Note this is not the fcn, but the inverse of the integral...)
  
  rnd = gRandom->Rndm(1); // (0,1)
  x=pow(-log(R1+rnd*R2),2.5); 
  
  return x;   
}

//Pick Interaction
AskCons::InteractionType_t UsefulPhysics::GetCurrent() {
  // choose CC or NC
  //  get from ratios in Ghandi etal paper
  // updated for the CTEQ6-DIS parton distribution functions
  double rnd=gRandom->Rndm();
  if (rnd<=0.6865254) // 10^18 eV - 10^21 eV (use this one for ANITA)
    //if (rnd<=0.6893498) // 10^17 eV - 10^20 eV (use this one for SalSA)
    return AskCons::kCC;
  else
    return AskCons::kNC; 

  return AskCons::kCC;
} //GetCurrent

AskCons::ParticleType_t UsefulPhysics::GetNuFlavor() {
    // pick a neutrino type, flavor ratio 1:1:1
  double rnd=gRandom->Rndm();
  AskCons::ParticleType_t particle=AskCons::kNue;
  if (rnd<=(1./3.)) {  
    particle=AskCons::kNue;
  } //if
  else if(rnd<=(2./3.)) { 
    particle=AskCons::kNumu;
  } //else if
  else if(rnd<=(1.)) { 
    particle=AskCons::kNutau;
  } //else if
  else
    cout << "unable to pick nu flavor\n";

  //Adding this 'cause I don't know any better
  rnd=gRandom->Rndm();
  if(rnd>0.666666)
     particle=AskCons::ParticleType_t(-1*int(particle));
  return particle;

} //GetNuFlavor

void UsefulPhysics::FillGzkFlux()
{

  //double E2dNdEdAdt[12]; //log(brightness)

  double Emuons[12]; // E^2 dN/dE/dA/dt for neutrinos that are produced as muon neutrinos or muon antineutrinos.
  double Eelectrons[12];// E^2 dN/dE/dA/dt for neutrinos that are produced as electron neutrinos or muon antineutrinos.
  
  for (int i=0;i<12;i++) {
    UPThings::energy[i]=16.+((double)i)/2.;
    Emuons[i]=-30.;
    Eelectrons[i]=-30.;
  } //for
  // what I was using previously, from upper curve of ANITA proposal
//    E2dNdEdAdt[0]=-9.6; // 16.
//    E2dNdEdAdt[1]=-8.9; // 16.5
//    E2dNdEdAdt[2]=-8.1; // 17.
//    E2dNdEdAdt[3]=-7.5; // 17.5
//    E2dNdEdAdt[4]=-7.2; // 18.
//    E2dNdEdAdt[5]=-6.8; // 18.5
//    E2dNdEdAdt[6]=-6.7; // 19
//    E2dNdEdAdt[7]=-6.8; // 19.5
//    E2dNdEdAdt[8]=-7.2; // 20.
//    E2dNdEdAdt[9]=-7.5; // 20.5
//    E2dNdEdAdt[10]=-8.2; // 21.0
//    E2dNdEdAdt[11]=-9.1; // 21.5

  // electron component of Figure 4 of ES&S
  // astro-ph/0101216
  Eelectrons[0]=-17.2; // 16.
  Eelectrons[1]=-17.35; // 16.5
  Eelectrons[2]=-17.2; // 17.
  Eelectrons[3]=-17.1; // 17.5
  Eelectrons[4]=-17.2; // 18.
  Eelectrons[5]=-17.5; // 18.5
  Eelectrons[6]=-18.0; // 19
  Eelectrons[7]=-18.5; // 19.5
  Eelectrons[8]=-19.4; // 20.
  Eelectrons[9]=-30.; // 20.5 punt
  Eelectrons[10]=-30.; // 21.0 punt
  Eelectrons[11]=-30.; // 21.5 punt

  // muon component of Figure 4 of ES&S
  // astro-ph/0101216
//    Emuons[0]=-17.8; // 16.
//    Emuons[1]=-17.4; // 16.5
//    Emuons[2]=-17.; // 17.
//    Emuons[3]=-16.75; // 17.5
//    Emuons[4]=-16.9; // 18.
//    Emuons[5]=-17.2; // 18.5
//    Emuons[6]=-17.7; // 19
//    Emuons[7]=-18.3; // 19.5
//    Emuons[8]=-19.1; // 20.
//    Emuons[9]=-30.; // 20.5 punt
//    Emuons[10]=-30.; // 21.0 punt
//    Emuons[11]=-30.; // 21.5 punt

  // lower curve of Figure 9 of ES&S
  // astro-ph/0101216
//    Emuons[0]=-17.1;  //16.
//    Emuons[1]=-16.6;  //16.5
//    Emuons[2]=-16.3;  //17.
//    Emuons[3]=-16.2; // 17.5
//    Emuons[4]=-16.4; // 18.
//    Emuons[5]=-16.7; // 18.5
//    Emuons[6]=-17.3; // 19
//    Emuons[7]=-17.95; // 19.5
//    Emuons[8]=-18.85; // 20.
//    Emuons[9]=-19.9; // 20.5 punt
//    Emuons[10]=-30.; // 21.0 punt
//    Emuons[11]=-30.; // 21.5 punt


  // upper curve in Figure 9 of ES&S
  // astro-ph/0101216
  Emuons[0]=-16.85;  //16.
  Emuons[1]=-16.4;  //16.5
  Emuons[2]=-16.05;  //17.
  Emuons[3]=-16.; // 17.5
  Emuons[4]=-16.15; // 18.
  Emuons[5]=-16.5; // 18.5
  Emuons[6]=-17.1; // 19
  Emuons[7]=-17.7; // 19.5
  Emuons[8]=-18.65; // 20.
  Emuons[9]=-19.75; // 20.5 punt
  Emuons[10]=-30.; // 21.0 punt
  Emuons[11]=-30.; // 21.5 punt

  for (int i=0;i<12;i++) {
    UPThings::EdNdEdAdt[i]=pow(10.,Eelectrons[i])+pow(10.,Emuons[i]);
//      cout << "UPThings::EdNdEdAdt is " << UPThings::EdNdEdAdt[i] << "\n";
//      cout << "energy is " << UPThings::energy[i] << "\n";
//      cout << "log(UPThings::EdNdEdAdt) is " << log10(UPThings::EdNdEdAdt[i]) << "\n";
  }

  //  for (int i=0;i<12;i++) {
  //UPThings::EdNdEdAdt[i]=pow(10,E2dNdEdAdt[i]-(UPThings::energy[i]-9.));
//      cout << "UPThings::EdNdEdAdt is " << UPThings::EdNdEdAdt[i] << "\n";
//      cout << "energy is " << UPThings::energy[i] << "\n";
//      cout << "log(UPThings::EdNdEdAdt) is " << log10(UPThings::EdNdEdAdt[i]) << "\n";
  //} //for

  
   UPThings::maxflux=TMath::MaxElement(12,UPThings::EdNdEdAdt);

    for (int i=0;i<12;i++) {
      UPThings::EdNdEdAdt[i]=UPThings::EdNdEdAdt[i]/UPThings::maxflux;
    } //for

}




double UsefulPhysics::PickGzkEnergy()
   //Returns energy in GeV
{
  static int doneInit=0;
  if(!doneInit) {
    FillGzkFlux();
    doneInit=1;
  }
  
  double thisenergy=16.;
  double thisflux=2.;
  double max=1.;
  int energybin=0;
  while(thisflux>max) {
    // pick an energy  
    thisenergy=gRandom->Rndm()*(TMath::MaxElement(12,UPThings::energy)-TMath::MinElement(12,UPThings::energy));
    energybin=(int)(thisenergy/0.5);
    max=UPThings::EdNdEdAdt[energybin];
    thisflux=gRandom->Rndm();
  } //while
  return pow(10.,thisenergy+TMath::MinElement(12,UPThings::energy))/1e9;
}

double UsefulPhysics::PickGzkIronEnergy()
   //Returns energy in GeV
{


  double energy[8];
 
  for (int i=0;i<8;i++) {
    energy[i]=12.+((double)i);
  } //for
  
  UPThings::EdNdEdAdt[0]=-17.1; // 12.
  UPThings::EdNdEdAdt[1]=-16.5; // 13.
  UPThings::EdNdEdAdt[2]=-16.1; // 14.
  UPThings::EdNdEdAdt[3]=-16.5; // 15.
  UPThings::EdNdEdAdt[4]=-17.3; // 16.
  UPThings::EdNdEdAdt[5]=-16.8; // 17.
  UPThings::EdNdEdAdt[6]=-17.2; // 18.
  UPThings::EdNdEdAdt[7]=-19.; // 19.


  for (int i=0;i<12;i++) {
    UPThings::EdNdEdAdt[i]=pow(10,UPThings::EdNdEdAdt[i]);
  } //for
  
  double maxflux=TMath::MaxElement(8,UPThings::EdNdEdAdt);

  for (int i=0;i<8;i++) {
    UPThings::EdNdEdAdt[i]=UPThings::EdNdEdAdt[i]/maxflux;
  } //for

  // now throw at a dartboard.
  
  double thisenergy=16.;
  double thisflux=2.;
  double max=1.;
  int energybin=0;
  while(thisflux>max) {
    // pick an energy  
    thisenergy=gRandom->Rndm()*(TMath::MaxElement(8,energy)-TMath::MinElement(8,energy));
    energybin=(int)(thisenergy);
    max=UPThings::EdNdEdAdt[energybin];
    thisflux=gRandom->Rndm(); 
  } //while
  
  for (int i=0;i<8;i++) {
    UPThings::EdNdEdAdt[i]=UPThings::EdNdEdAdt[i]*maxflux;
  } //for 

  return pow(10.,thisenergy+TMath::MinElement(8,energy))/1e9;;

}


