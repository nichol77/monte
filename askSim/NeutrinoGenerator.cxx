///////////////////////////////////////////////////////////////////////////////
/////   NeutrinoGenerator.cxx
/////   Contains information about the various interactions (primary and
/////   secondary)
///////////////////////////////////////////////////////////////////////////////
#include "NeutrinoGenerator.h"
#include "UsefulPhysics.h"
#include "NeutrinoCrossSectionModel.h"
#include "AskGeom.h"



#include <iostream> 
#include "interactions.h"



#include "TF1.h"
#include "TRandom.h"
#include "TMath.h"

using namespace std;

ClassImp(NeutrinoGenerator)

NeutrinoGenerator::NeutrinoGenerator()
   :fRunId(0),fVolRad(AskGeom::R_EARTH),fRollCounter(0),fIntCounter(0),fMaxDistToIce(0),fUseful(0),fHistLatLonSurf(0)
{
//Default (zero) constructor 
}


NeutrinoGenerator::~NeutrinoGenerator()
{
//Default Destructor
}

NeutrinoGenerator::NeutrinoGenerator(Long64_t runId,  AskCons::WhichFlux_t whichFlux, Double_t maxDistToIce, Int_t doZCut,Double_t zCut)
:fRunId(runId), fVolRad(AskGeom::R_EARTH),fRollCounter(0),fIntCounter(0),fMaxDistToIce(maxDistToIce),fWhichFlux(whichFlux),fFixedEnergy(0),fUseful(0)
{
  if(doZCut) {
    fDoZCut=doZCut;
    fZCut=zCut;
  }
  fRandomSeed=gRandom->GetSeed();
}

void NeutrinoGenerator::setRandomSeed(UInt_t seed) {
      gRandom->SetSeed(seed);
      fRandomSeed=gRandom->GetSeed();
}



void NeutrinoGenerator::pickRandomThetaPhiOnSphere(Double_t &theta, Double_t &phi)
{

  phi=pickRandomPhiTwoPi();
  theta=TMath::ACos(pickRandomCosThetaOnSphere());

}

void NeutrinoGenerator::pickRandomThetaPhiOnSphere(Double_t &theta, Double_t &phi, Double_t &costheta)
{
  //Just a hack for now doesn't pick theta

  phi=pickRandomPhiTwoPi();
  costheta=pickRandomCosThetaOnSphere();
  //  theta=TMath::ACos(costheta);
}





void NeutrinoGenerator::pickRandomDowngoingDirection(Double_t &theta, Double_t &phi)
{

  phi=pickRandomPhiTwoPi();
  theta=TMath::ACos(gRandom->Rndm());

}

void NeutrinoGenerator::pickRandomDowngoingDirection(Double_t &theta, Double_t &phi, Double_t &costheta)
  //Just a hack for now doesn't pick theta
{

  phi=2*TMath::Pi()*gRandom->Rndm();
  costheta=pickRandomCosThetaHalfPi();
  //  theta=TMath::ACos(costheta);

}

void NeutrinoGenerator::pickPointAndDirectionOnVolumeInGlobalCoords(TVector3 &intPos, TVector3 &intDir)
{

  Double_t thetaV,phiV;
  Double_t thetaMom,phiMom;
  
  pickRandomThetaPhiOnSphere(thetaV,phiV);
  pickRandomDowngoingDirection(thetaMom,phiMom);
  intDir.SetX(-1*TMath::Cos(phiMom)*TMath::Sin(thetaMom));
  intDir.SetY(-1*TMath::Sin(phiMom)*TMath::Sin(thetaMom));
  intDir.SetZ(-1*TMath::Cos(thetaMom));
  
  intDir.RotateZ(phiV);
  intDir.RotateX(thetaV);


  double rNew=TMath::Sqrt(AskGeom::R_EARTH*AskGeom::R_EARTH + fVolRad*fVolRad - 2*AskGeom::R_EARTH*fVolRad*TMath::Cos(TMath::Pi()-thetaV));
  double thetaNew=TMath::ASin(TMath::Sin(TMath::Pi()-thetaV)*fVolRad/rNew);
  double phiNew=phiV;
  
  intPos.SetZ(rNew*TMath::Cos(thetaNew));
  intPos.SetX(rNew*TMath::Sin(thetaNew)*TMath::Cos(phiNew));
  intPos.SetY(rNew*TMath::Sin(thetaNew)*TMath::Sin(phiNew));
  
}

void NeutrinoGenerator::pickPointAndDirectionOnEarth(TVector3 &intPos, TVector3 &intDir)
{
  setPointAndDirectionOnEarth();
  intPos=fSurfPos;
  intDir=fIntDir;



}

int NeutrinoGenerator::setPointAndDirectionOnEarth() {

  Double_t phiV; //thetaV;
  Double_t costv,sintv,cospv,sinpv;
  Double_t thetaMom,phiMom;
  Double_t costm,sintm,cospm,sinpm;
  
  costv=pickRandomCosThetaOnSphere();  
  fSurfCosTheta=costv;
  if(fDoZCut) {
    if(costv*AskGeom::R_EARTH < fZCut)
      return 0;
  }
  phiV=pickRandomPhiTwoPi();   
  pickRandomDowngoingDirection(thetaMom,phiMom,costm);
  //thetaMom not set
  //thetaV not set

  //Calcuate sins and coss
  cospm=TMath::Cos(phiMom);
  //  sinpm=TMath::Sin(phiMom);
  sinpm=AskGeom::getSinFromCos(cospm);
  //  costm=TMath::Cos(thetaMom);
  //  sintm=TMath::Sin(thetaMom);
  sintm=AskGeom::getSinFromCos(costm);

  //Set initial direction
  fIntDir.SetX(cospm*sintm);
  fIntDir.SetY(sinpm*sintm);
  fIntDir.SetZ(costm);
  

  //Calculate more sins and coss
  cospv=TMath::Cos(phiV);
  //  sinpv=TMath::Sin(phiV);
  sinpv=AskGeom::getSinFromCos(cospv);
  if(phiV>TMath::Pi()) 
    sinpv*=-1;
  //  costv=TMath::Cos(thetaV);
  //  sintv=TMath::Sin(thetaV);
  sintv=AskGeom::getSinFromCos(costv);


  //Rotate initial direction to be in usual frame
  //  fIntDir.RotateZ(phiV);
  //  fIntDir.RotateX(thetaV);
  AskGeom::fastRotateZ(fIntDir,sinpv,cospv);
  AskGeom::fastRotateX(fIntDir,sintv,costv);


  //Set initial position
  fSurfPos.SetX(AskGeom::R_EARTH*sintv*cospv);
  fSurfPos.SetY(AskGeom::R_EARTH*sintv*sinpv);
  fSurfPos.SetZ(AskGeom::R_EARTH*costv);


  //Work out which dir is in
  //Hack because I can't do geometry
  Double_t testVal=TMath::Sqrt((fSurfPos.X()+fIntDir.X())*(fSurfPos.X()+fIntDir.X())+
			       (fSurfPos.Y()+fIntDir.Y())*(fSurfPos.Y()+fIntDir.Y())+
			       (fSurfPos.Z()+fIntDir.Z())*(fSurfPos.Z()+fIntDir.Z()));
  if(testVal>AskGeom::R_EARTH)
    fIntDir*=-1;
  

  return 1;
}

int NeutrinoGenerator::getSurfacePointAndDirectionFromVolumePoint(TVector3 intPos, TVector3 &intDir, TVector3 &surfPos) 
{
  Double_t b=intDir.X()*intPos.X() + intDir.Y()*intPos.Y() + intDir.Z()*intPos.Z();
  Double_t c = intPos.X()*intPos.X() + intPos.Y()*intPos.Y() + intPos.Z()*intPos.Z() - AskGeom::R_EARTH*AskGeom::R_EARTH;
  if(b*b < c) 
    return 0;


  Double_t l1=-1*b + TMath::Sqrt(b*b - c);
  Double_t l2=-1*b - TMath::Sqrt(b*b - c);
  
  
  if(intPos.Mag2()> AskGeom::R_EARTH*AskGeom::R_EARTH) {
    //Start outside take l1
    //return 0;
    //    cout << "Outside:\t" << l1 << "\t" << l2 << endl;
    intDir*=-1;
    l1*=-1;
    surfPos.SetX(intPos.X() + l1*intDir.X());
    surfPos.SetY(intPos.Y() + l1*intDir.Y());
    surfPos.SetZ(intPos.Z() + l1*intDir.Z());
    //    cout << surfPos.Mag() << endl;
  }
  else {
    //Start inside take l2
    //    cout << "Inside:\t" << l1 << "\t" << l2 <<  endl;
    //    cout << l2 << endl;
    //    return 0;
    surfPos.SetX(intPos.X() + l2*intDir.X());
    surfPos.SetY(intPos.Y() + l2*intDir.Y());
    surfPos.SetZ(intPos.Z() + l2*intDir.Z());
    //    cout << surfPos.Mag() << endl;
  }
  
  return 1;
}


int NeutrinoGenerator::getSurfacePointAndDirection(TVector3 &surfPos, TVector3 &intDir)
{
  TVector3 intPos;
  pickPointAndDirectionOnVolumeInGlobalCoords(intPos, intDir);
  return getSurfacePointAndDirectionFromVolumePoint(intPos,intDir,surfPos);
}



int NeutrinoGenerator::getIceThicknessesAndLengths(Double_t &iceStart, Double_t &iceEnd)
{
  //  if(fRollCounter>1350)
  //    std::cout << "getIceThicknessesAndLengths\t" << fRollCounter << "\n";
  // TVector3 radialDir=fSurfPos.Unit();
  //Assume magnitude is radius of earth
  //  fRadialDir.SetX(fSurfPos.X()/AskGeom::R_EARTH);
  //  fRadialDir.SetY(fSurfPos.Y()/AskGeom::R_EARTH);
  //  fRadialDir.SetZ(fSurfPos.Z()/AskGeom::R_EARTH);

  

   //First find chord across world
//   Double_t chordTheta=TMath::Pi()-fRadialDir.Angle(fIntDir);
   Double_t cosChordTheta=AskGeom::getCosThetaBewteenVectors(fSurfPos,fIntDir);
   fChordLength=TMath::Abs(2*AskGeom::R_EARTH*cosChordTheta);

   fOppPoint.SetX(fSurfPos.X()+fChordLength*fIntDir.X());
   fOppPoint.SetY(fSurfPos.Y()+fChordLength*fIntDir.Y());
   fOppPoint.SetZ(fSurfPos.Z()+fChordLength*fIntDir.Z());

   //Check if we are anywhere near ice
   if((fSurfPos.Z()/AskGeom::R_EARTH<0.85) && (fOppPoint.Z()/AskGeom::R_EARTH<0.85)) {
     //     if(fRollCounter>1350)       
     //       std::cout << "\t end getIceThicknessesAndLengths\t" << fRollCounter << "\n";
     return 0;
   }
   
   //   if(fRollCounter>1350)
   //     std::cout << "aa\t" << fChordLength << std::endl;
   //Now if we are try and find actual surface points
   fStartRadius=fWorldModel->getSurface(fSurfPos);
   //   if(fRollCounter>1350)
   //     std::cout << "sr\t" << fStartRadius << std::endl;
   fSurfPos.SetMag(fStartRadius);
   fChordLength=TMath::Abs(2*fStartRadius*cosChordTheta);


   //Refind opposite point
   fOppPoint.SetX(fSurfPos.X()+fChordLength*fIntDir.X());
   fOppPoint.SetY(fSurfPos.Y()+fChordLength*fIntDir.Y());
   fOppPoint.SetZ(fSurfPos.Z()+fChordLength*fIntDir.Z());


   Double_t fEndRadius=fWorldModel->getSurface(fOppPoint);
   fEndChordLength=TMath::Abs(2*fEndRadius*cosChordTheta);
   Double_t endDistToEarth=(fChordLength-fEndChordLength)/2;
   
   //   if(fRollCounter>1350)
   //     std::cout << "b\t" << fChordLength << "\t" << endDistToEarth << std::endl;
   //Now move end points
   //Seems a bit dodgy but there we are
   fChordLength-=endDistToEarth;
   fOppPoint.SetX(fOppPoint.X()-endDistToEarth*fIntDir.X());
   fOppPoint.SetY(fOppPoint.Y()-endDistToEarth*fIntDir.Y());
   fOppPoint.SetZ(fOppPoint.Z()-endDistToEarth*fIntDir.Z());


   //Time to check for ice depth
   iceStart=0,iceEnd=0;

   fStartThickness=0;
   fEndThickness=0;
   //In Northern hemisphere so no ice
   //Placing a cut at about 30 degrees (ie. only latitudes
   if(fSurfPos.Z()/AskGeom::R_EARTH>0.9)
      fStartThickness=fWorldModel->getIceThickness(fSurfPos);
   if(fOppPoint.Z()/AskGeom::R_EARTH>0.9)
      fEndThickness=fWorldModel->getIceThickness(fOppPoint);
   
   
   //   if(fRollCounter>1350)
   //     std::cout << "c\t" << fStartThickness << "\t" << fEndThickness << std::endl;

   if(fStartThickness<0.1 && fEndThickness<0.1) {
     //Neither end near ice give up
     //     if(fRollCounter>1350)
     //       std::cout << "\t end getIceThicknessesAndLengths\t" << fRollCounter << "\n";
     return 0;
   }
   
   //This is only supposed to be an estimate of where the ice is
   if(fStartThickness>0 && fEndThickness>0) {
      //Start in ice end in ice
      iceStart=0;
      iceEnd=fChordLength*AskCons::getDensityForMaterial(AskCons::kIce);
   } 
   else if(fEndThickness>0) {
     cosChordTheta=AskGeom::getCosThetaBewteenVectors(fOppPoint,fIntDir);
     Double_t stepLength=TMath::Abs(fEndThickness/cosChordTheta);

     //Start in rock end in ice
     iceStart=(fChordLength-stepLength)*AskCons::getDensityForMaterial(AskCons::kRock);
     iceEnd=iceStart+(stepLength*AskCons::getDensityForMaterial(AskCons::kIce));
   }
   else if(fStartThickness>0) {
     cosChordTheta=AskGeom::getCosThetaBewteenVectors(fSurfPos,fIntDir);     
     Double_t stepLength=TMath::Abs(fStartThickness/cosChordTheta);
//       //Start in ice end in rock
      iceStart=0;
      iceEnd=stepLength*AskCons::getDensityForMaterial(AskCons::kIce);
   }
   //   if(fRollCounter>1350)
   //     std::cout << "\t end getIceThicknessesAndLengths\t" << fRollCounter << "\n";

   return 1;

//    if(fStartThickness<1 && fEndThickness>0) {
//       Double_t fEndThickness=fWorldModel->getIceThickness(fOppPoint.X(),fOppPoint.Y());
//       Double_t stepSize=TMath::Abs(fEndThickness/cosChordTheta);
      
//       TVector3 newPoint=fOppPoint-stepSize*fIntDir;
//       Double_t newThickness=fWorldModel->getIceThickness(newPoint.X(),newPoint.Y());
//       Int_t stepCount=0;
//       Double_t stepLength=stepSize;
//       //      std::cout << "Step Ice:\t" << fStartThickness << "\t" << newThickness << "\t" << AskGeom::R_EARTH-newPoint.Mag() << std::endl;
//       Double_t newStep=newThickness-(AskGeom::R_EARTH-newPoint.Mag());
//       Int_t converged=1;
//       while(TMath::Abs(newThickness-(AskGeom::R_EARTH-newPoint.Mag()))>50) {
// 	 stepLength+=newStep;
// 	 newPoint=newPoint-newStep*fIntDir;
// 	 newThickness=fWorldModel->getIceThickness(newPoint.X(),newPoint.Y());
// 	 //	 if(stepCount>10) {
// 	 //	    std::cout << "Step Ice:\t" << stepCount << "\t" << fStartThickness << "\t" << newThickness << "\t" << AskGeom::R_EARTH-newPoint.Mag() << std::endl;
// 	 //	 }
// 	 stepCount++;
// 	 newStep=newThickness-(AskGeom::R_EARTH-newPoint.Mag());
// 	 if(stepCount>100) {
// 	    converged=0;
// 	    break;
// 	 }
//       }
//       if(!converged) {
// 	 //Just step linearly
// 	 stepSize=stepSize/2.;
// 	 newPoint=fOppPoint-stepSize*fIntDir;
// 	 newThickness=fWorldModel->getIceThickness(newPoint.X(),newPoint.Y());
// 	 stepLength=stepSize;
// 	 stepCount=0;
// 	 //For downgoing guys just step linearly
// 	 stepSize=100;

// 	 while(newThickness-(AskGeom::R_EARTH-newPoint.Mag())>50 && newThickness>0.1) {
// 	    stepLength+=stepSize;
// 	    newPoint=newPoint-stepSize*fIntDir;
// 	    newThickness=fWorldModel->getIceThickness(newPoint.X(),newPoint.Y());
// // 	    //	 if(stepCount>10) {
// // 	    if(fRollCounter==554291) {
// // 	       std::cout << "Step Ice:\t" << stepCount << "\t" << stepLength << "\t" << fChordLength << "\t"  << fStartThickness << "\t" << newThickness << "\t" << AskGeom::R_EARTH-newPoint.Mag() << std::endl;
// // 	    }
// 	    stepCount++;
// 	    if(stepCount==100)
// 	       stepSize*=5;
// 	    if(stepCount>5000) {
// 	       std::cout << "Didn't converge\n";
// 	       std::cout << fRollCounter << "\n";
// 	       std::cout << stepCount << "\t" << stepSize <<"\t" << stepLength 
// 			 << "\t"<< newThickness << "\t" << AskGeom::R_EARTH-newPoint.Mag() << std::endl;
// 	       std::cout << fStartThickness << "\t" << fEndThickness << "\t" << fChordLength << "\n";
// 	       std::cout << fSurfPos.X() << "\t" << fSurfPos.Y() << "\t" << fSurfPos.Z()
// 			 << "\t" << fSurfPos.Mag() << "\n";
// 	       std::cout << fIntDir.X() << "\t" << fIntDir.Y() << "\t" << fIntDir.Z()
// 			 << "\t" << fIntDir.Mag() << "\n";
// 	       std::cout << fOppPoint.X() << "\t" << fOppPoint.Y() << "\t" << fOppPoint.Z()
// 		   << "\t" << fOppPoint.Mag() << "\n";
// 	       break;
// 	    }
// 	 }
//       }
//       //      std::cout << stepCount << "\t" << stepSize <<"\t" << stepLength << std::endl;
      
//       //Start in rock end in ice
//       iceStart=(fChordLength-stepLength)*AskCons::getDensityForMaterial(AskCons::kRock);
//       iceEnd=iceStart+(stepLength*AskCons::getDensityForMaterial(AskCons::kIce));
      
//    }
//    else if(fStartThickness>0 && fEndThickness>0) {
//       //Start in ice end in ice
//       iceStart=0;
//       iceEnd=fChordLength*AskCons::getDensityForMaterial(AskCons::kIce);
//    }
//    else if(fStartThickness>0 && fEndThickness<1) {
      
//       Double_t fStartThickness=fWorldModel->getIceThickness(fSurfPos.X(),fSurfPos.Y());
//       Double_t stepSize=TMath::Abs(fStartThickness/cosChordTheta);
      
//       TVector3 newPoint=fOppPoint+stepSize*fIntDir;
//       Double_t newThickness=fWorldModel->getIceThickness(newPoint.X(),newPoint.Y());
//       Int_t stepCount=0;
//       Double_t stepLength=stepSize;
//       //      std::cout << "Step Ice:\t" << fStartThickness << "\t" << newThickness << "\t" << AskGeom::R_EARTH-newPoint.Mag() << std::endl;
//       Double_t newStep=newThickness-(AskGeom::R_EARTH-newPoint.Mag());

//       Int_t converged=1;
//       while(TMath::Abs(newThickness-(AskGeom::R_EARTH-newPoint.Mag()))>50 && newThickness) {
// 	 stepLength+=newStep;
// 	 newPoint=newPoint+newStep*fIntDir;
// 	 newThickness=fWorldModel->getIceThickness(newPoint.X(),newPoint.Y());
// 	 //	 if(stepCount>10) {
// 	 //	    std::cout << "Step Ice:\t" << stepCount << "\t" << fStartThickness << "\t" << newThickness << "\t" << AskGeom::R_EARTH-newPoint.Mag() << std::endl;
// 	 //	 }
// 	 stepCount++;
// 	 newStep=newThickness-(AskGeom::R_EARTH-newPoint.Mag());
// 	 if(stepCount>100) {
// 	    converged=0;
// 	    break;
// 	 }
//       }
//       if(!converged) {
// 	 //For downgoing guys just step linearly
// 	 stepSize=stepSize/2;
// 	 newPoint=fSurfPos+stepSize*fIntDir;
// 	 newThickness=fWorldModel->getIceThickness(newPoint.X(),newPoint.Y());
// 	 stepLength=stepSize;
// 	 stepCount=0;
// 	 stepSize=100;

// // 	 if(fRollCounter==277) {
// // 	   std::cout << fRollCounter << "\t" << stepSize << "\t" << fSurfPos.Mag()
// // 		     << "\t"<< newPoint.Mag() << "\n";
// // 	   std::cout << fSurfPos.X() << "\t" << fSurfPos.Y() << "\t" << fSurfPos.Z() << "\t" << fSurfPos.Mag() << "\n";
// // 	   std::cout << fOppPoint.X() << "\t" << fOppPoint.Y() << "\t" << fOppPoint.Z() << "\t" << fOppPoint.Mag() << "\n";
// // 	 }
// 	 //For downgoing guys just step linearly

// 	 while(newThickness-(AskGeom::R_EARTH-newPoint.Mag())>50 && newThickness>0.1) {
// 	    stepLength+=stepSize;
// 	    newPoint=newPoint+stepSize*fIntDir;
// 	    newThickness=fWorldModel->getIceThickness(newPoint.X(),newPoint.Y());
// 	//     if(fRollCounter==277) {
// // 	      //	      std::cout << "Step Ice:\t" << stepCount << "\t" << stepLength << "\t" << fStartThickness << "\t" << newThickness << "\t" << AskGeom::R_EARTH-newPoint.Mag() << std::endl;
// // 	    }
// 	    stepCount++;
// 	    if(stepCount==100)
// 	       stepSize*=5;
// 	    if(stepCount>1000) {
// 	       std::cout << "Didn't converge\n";
// 	       std::cout << fRollCounter << "\n";
// 	       std::cout << stepCount << "\t" << stepSize <<"\t" << stepLength 
// 			 << "\t"<< newThickness << "\t" << AskGeom::R_EARTH-newPoint.Mag() << std::endl;
// 	       std::cout << fStartThickness << "\t" << fEndThickness << "\t" << fChordLength << "\n";
// 	       std::cout << fSurfPos.X() << "\t" << fSurfPos.Y() << "\t" << fSurfPos.Z()
// 			 << "\t" << fSurfPos.Mag() << "\n";
// 	       std::cout << fIntDir.X() << "\t" << fIntDir.Y() << "\t" << fIntDir.Z()
// 			 << "\t" << fIntDir.Mag() << "\n";
// 	       std::cout << fOppPoint.X() << "\t" << fOppPoint.Y() << "\t" << fOppPoint.Z()
// 		   << "\t" << fOppPoint.Mag() << "\n";

// 	       break;
// 	    }
// 	 }
// 	 //	 std::cout << stepCount << "\t" << stepSize <<"\t" << stepLength 
// 	 //		   << "\t"<< newThickness << "\t" << AskGeom::R_EARTH-newPoint.Mag() << std::endl;
//       }
//       //Start in ice end in rock
//       iceStart=0;
//       iceEnd=stepLength*AskCons::getDensityForMaterial(AskCons::kIce);
//    }
}

NeutrinoEvent *NeutrinoGenerator::getNextNeutrino()
{
  fRollCounter++;
   NeutrinoEvent *theEvent=0;
   if(!fUseful)
      fUseful=UsefulPhysics::Instance();
   if(!fWorldModel)
      fWorldModel=WorldModel::Instance();
   //   if(!fHistLatLonSurf)
   //      fHistLatLonSurf= new TH2F("histLatLonSurf","Latitude v Longitude of Neutrino Entry",360,-180,180,180,-90,90);

   //   TVector3 surfPos;
   //   TVector3 intDir;
   
   //   pickPointAndDirectionOnEarth(surfPos,intDir);
   int isValid=setPointAndDirectionOnEarth();
   if(!isValid) {
     return 0;
   }

   //   int isValid=getSurfacePointAndDirection(fSurfPos, fIntDir);
   //   if(!isValid) {
   //      return 0;
   //   }

   //Find if it is near ice
   Double_t iceStart=0,iceEnd=0;
   isValid=getIceThicknessesAndLengths(iceStart,iceEnd);
   if(!isValid) {

        return 0;
   }
   

   //Now get energy and flavour
   AskCons::ParticleType_t particle=fUseful->GetNuFlavor();
   Double_t energy=fFixedEnergy;
   if(fWhichFlux==AskCons::kStandardGzk) 
     energy=fUseful->PickGzkEnergy();

   //Now find interaction length and roll for actual length
   double intLengthWater=NeutrinoCrossSectionModel::getInteractionLength(energy,particle,AskCons::kWater);
   double actualLengthWater=NeutrinoCrossSectionModel::getDistanceToInteraction(intLengthWater);

   //So if interaction is between iceStart & iceEnd interaction is in ice
   //if it is within Xm of the ice we will keep it around for secondary interactions.

   //   std::cout << actualLengthWater << "\t" << iceStart << "\t" << iceEnd << std::endl;

   if((actualLengthWater>iceStart && actualLengthWater<iceEnd) || 
      (actualLengthWater<iceEnd && actualLengthWater+fMaxDistToIce>iceStart)) {
      //Got an in-ice or near-ice interaction
      AskCons::InteractionType_t intType=fUseful->GetCurrent();
      Double_t elasticity=fUseful->Gety();
      
      Double_t actualLength=0;
      if(actualLengthWater>iceStart) {
	 actualLength=iceStart/AskCons::getDensityForMaterial(AskCons::kRock);
	 actualLength+=(actualLengthWater-iceStart)/AskCons::getDensityForMaterial(AskCons::kIce);
      }
      else {
	 actualLength=actualLengthWater/AskCons::getDensityForMaterial(AskCons::kRock);
      }
	    
      fIntPoint.SetX(fSurfPos.X() + actualLength*fIntDir.X());
      fIntPoint.SetY(fSurfPos.Y() + actualLength*fIntDir.Y());
      fIntPoint.SetZ(fSurfPos.Z() + actualLength*fIntDir.Z());
      fIntRadius=fWorldModel->getSurface(fIntPoint);
      if(fIntPoint.Mag()>fIntRadius) {
	//Ooops
	return 0;
      }
      Double_t iceThickness=fWorldModel->getIceThickness(fIntPoint);
      AskCons::MaterialType_t matType = AskCons::kRock;
      if(iceThickness>0 && (fIntRadius-fIntPoint.Mag())<iceThickness)
	matType=AskCons::kIce;
      
      if(fIntPoint.Z()<3e6) {
	 std::cout << "Unphysical Point\n";
	 std::cout << fSurfPos.X() << "\t" << fSurfPos.Y() << "\t" << fSurfPos.Z()
		   << "\t" << fSurfPos.Mag() << "\n";
	 std::cout << fIntDir.X() << "\t" << fIntDir.Y() << "\t" << fIntDir.Z()
		   << "\t" << fIntDir.Mag() << "\n";
	 std::cout << fIntPoint.X() << "\t" << fIntPoint.Y() << "\t" << fIntPoint.Z()
		   << "\t" << fIntPoint.Mag() << "\n";
	 std::cout << fOppPoint.X() << "\t" << fOppPoint.Y() << "\t" << fOppPoint.Z()
		   << "\t" << fOppPoint.Mag() << "\n";
	 std::cout << actualLength<< "\t" << actualLengthWater << "\t" << iceStart << "\t" << iceEnd << std::endl;
	 std::cout << fChordLength << "\t" << fStartThickness << "\t" << fEndThickness << "\n";
      }
	 
	    

   
      theEvent = new NeutrinoEvent(fRunId,fRollCounter-1,fSurfPos,energy,particle,matType,
      				   fIntPoint,fIntDir,intType,elasticity);
      fIntCounter++;

   }
   return theEvent;
}

NeuGenRunSummary *NeutrinoGenerator::getRunSummary()
{
  NeuGenRunSummary *runSummary 
    = new NeuGenRunSummary(fRunId,fRollCounter,fIntCounter,fVolRad,fMaxDistToIce,
			   fWhichFlux,fFixedEnergy,1,fRandomSeed,fDoZCut,fZCut);
  return runSummary;
}
