///////////////////////////////////////////////////////////////////////////////
/////   NeutrinoCrossSectionModel.h
/////  Namespace for holding the cross section related code and calls
//////////////////////////////////////////////////////////////////////////////


#include "NeutrinoCrossSectionModel.h"
#include "UtilityFuncs.h"

#include "TMath.h"
#include "TRandom3.h"
#include "TF1.h"

#include <iostream>


//For interaction length thingy
double NeutrinoCrossSectionModel::intLengthFromProb(double *p, double *par) {
   double L=par[0];
   return -1*L*TMath::Log(1-p[0]);
}

namespace NeutrinoCrossSectionModel {
   TF1 *fIntLength = new TF1("fIntLengthTemp",intLengthFromProb,0,1,1);      

   Double_t fScaleFactor=1;
   WhichCrossection_t fWhichCross=kDefault;   


   //Add lookup tables here
   Int_t fNumNUSigmaPoints=34;
   Double_t fNuEnergy[34]={10,25,60,100,250,600,1000,2500,6000,10000,25000,60000,100000,250000,600000,1e+06,2.5e+06,6e+06,1e+07,2.5e+07,6e+07,1e+08,2.5e+08,6e+08,1e+09,2.5e+09,6e+09,1e+10,2.5e+10,6e+10,1e+11,2.5e+11,6e+11,1e+12};
   Double_t fNuSigmaCC[34]={7.988e-38,1.932e-37,4.45e-37,7.221e-37,1.728e-36,3.964e-36,6.399e-36,1.472e-35,3.096e-35,4.617e-35,8.824e-35,1.514e-34,2.022e-34,3.255e-34,4.985e-34,6.342e-34,9.601e-34,1.412e-33,1.749e-33,2.554e-33,3.63e-33,4.436e-33,6.283e-33,8.699e-33,1.049e-32,1.466e-32,2.01e-32,2.379e-32,3.289e-32,4.427e-32,5.357e-32,7.32e-32,9.927e-32,1.179e-31};
   Double_t fNuSigmaNC[34]={2.492e-38,6.033e-38,1.391e-37,2.261e-37,5.43e-37,1.255e-36,2.039e-36,4.781e-36,1.035e-35,1.575e-35,3.139e-35,5.615e-35,7.667e-35,1.28e-34,2.017e-34,2.6e-34,4.018e-34,6.001e-34,7.482e-34,1.104e-33,1.581e-33,1.939e-33,2.763e-33,3.837e-33,4.641e-33,6.49e-33,8.931e-33,1.066e-32,1.465e-32,1.995e-32,2.377e-32,3.247e-32,4.377e-32,5.196e-32};
   Double_t fNuSigmaTot[34]={1.048e-37,2.535e-37,5.841e-37,9.482e-37,2.271e-36,5.219e-36,8.438e-36,1.95e-35,4.131e-35,6.192e-35,1.196e-34,2.076e-34,2.789e-34,4.535e-34,7.002e-34,8.942e-34,1.362e-33,2.012e-33,2.497e-33,3.658e-33,5.211e-33,6.375e-33,9.046e-33,1.254e-32,1.513e-32,2.115e-32,2.903e-32,3.445e-32,4.754e-32,6.422e-32,7.734e-32,1.057e-31,1.43e-31,1.699e-31};
   Double_t fNubarEnergy[34]={10,25,60,100,250,600,1000,2500,6000,10000,25000,60000,100000,250000,600000,1e+06,2.5e+06,6e+06,1e+07,2.5e+07,6e+07,1e+08,2.5e+08,6e+08,1e+09,2.5e+09,6e+09,1e+10,2.5e+10,6e+10,1e+11,2.5e+11,6e+11,1e+12};
   Double_t fNubarSigmaCC[34]={3.936e-38,9.726e-38,2.287e-37,3.747e-37,9.154e-37,2.153e-36,3.542e-36,8.548e-36,1.922e-35,3.008e-35,6.355e-35,1.199e-34,1.683e-34,2.909e-34,4.667e-34,6.051e-34,9.365e-34,1.393e-33,1.734e-33,2.542e-33,3.622e-33,4.43e-33,6.278e-33,8.696e-33,1.05e-32,1.464e-32,2.011e-32,2.406e-32,3.286e-32,4.481e-32,5.335e-32,7.306e-32,9.854e-32,1.165e-31};
   Double_t fNubarSigmaNC[34]={1.381e-38,3.403e-38,7.982e-38,1.307e-37,3.193e-37,7.531e-37,1.243e-36,3.026e-36,6.896e-36,1.091e-35,2.358e-35,4.57e-35,6.515e-35,1.158e-34,1.901e-34,2.493e-34,3.929e-34,5.93e-34,7.423e-34,1.1e-33,1.578e-33,1.937e-33,2.762e-33,3.836e-33,4.641e-33,6.489e-33,8.931e-33,1.066e-32,1.465e-32,1.995e-32,2.377e-32,3.247e-32,4.377e-32,5.195e-32};
   Double_t fNubarSigmaTot[34]={5.317e-38,1.313e-37,3.085e-37,5.054e-37,1.235e-36,2.906e-36,4.785e-36,1.157e-35,2.612e-35,4.099e-35,8.713e-35,1.656e-34,2.334e-34,4.067e-34,6.568e-34,8.544e-34,1.329e-33,1.986e-33,2.476e-33,3.642e-33,5.2e-33,6.367e-33,9.04e-33,1.253e-32,1.514e-32,2.113e-32,2.904e-32,3.472e-32,4.751e-32,6.476e-32,7.712e-32,1.055e-31,1.423e-31,1.685e-31};
}

   
void NeutrinoCrossSectionModel::setWhichCrossection(NeutrinoCrossSectionModel::WhichCrossection_t whichCross)
{
   fWhichCross=whichCross;
}

void NeutrinoCrossSectionModel::setScaleFactor(Double_t scaleFactor) {
   fScaleFactor=scaleFactor;
}

Double_t NeutrinoCrossSectionModel::getCrossSection(Double_t energy, AskCons::ParticleType_t particle) {
   Double_t isAnti=0;
   if(particle<0)
      isAnti=1;
   

   switch(fWhichCross) {
   case kRenoCteq6:
      return getRenoCteq6CrossSection(energy);
   case kDefault:
   default:
      return funcTotCrossSec(&energy,&isAnti);
   }
   return funcTotCrossSec(&energy,&isAnti);
}



Double_t NeutrinoCrossSectionModel::getInteractionLength(Double_t energy, AskCons::ParticleType_t particle, AskCons::MaterialType_t mat) {
   
   Double_t isAnti=0;
   if(particle<0)
      isAnti=1;
   double density=AskCons::getDensityForMaterial(mat);
   double params[2]={isAnti,density};   
   return 1e3*funcNuTotInteractionLength(&energy,params); //convert from km to m
}
    

double NeutrinoCrossSectionModel::getDistanceToInteraction(double intLength) 
{
   fIntLength->SetParameter(0,intLength);
   return fIntLength->Eval(gRandom->Rndm());      
}


Double_t NeutrinoCrossSectionModel::getRenoCteq6CrossSection(double energy) {
  // calculate cross section
   Double_t sigma=0;
  if (energy<1.2E15) {
    std::cout <<  "Need a parameterization for this energy region.\n";
    return 0;
  } //if
  else {
    // fit to cross sections calculated by M.H. Reno using the same method as Gandhi et al, but with the CTEQ6-DIS parton distribution functions instead of the CTEQ4-DIS distribution functions
    sigma=(2.501E-39)*pow(energy/1.E9,0.3076)*fScaleFactor; // 10^18 eV - 10^21 eV (use this one for ANITA)
    //sigma=(1.2873E-39)*pow(energy/1.E9,0.33646)*UPThings::fScaleFactor; // 10^17 eV - 10^20 eV (use this one for SalSA)
  } //if
  // interaction length in kg/m^2  
  return sigma;
} //GetRenoCteq6CrossSection


Double_t NeutrinoCrossSectionModel::funcCCCrossSec(Double_t *energy, Double_t *par)
{
    //Energy in GeV
    //Cross-section in cm^2
    //This is the cross-section for neutrino-nucleon scattering, so must multiply by A to get the neutrino nucleus cross-section
    Int_t isAntiNeutrino=int(par[0]+0.1);
 //    if(isAntiNeutrino)
// 	return (5.52e-36)*TMath::Power(energy[0],0.363);
//     return (5.53e-36)*TMath::Power(energy[0],0.363);
    if(isAntiNeutrino)
	return UtilityFuncs::linearInterpolation(energy[0],fNubarEnergy,fNubarSigmaCC,fNumNUSigmaPoints);
    return UtilityFuncs::linearInterpolation(energy[0],fNuEnergy,fNuSigmaCC,fNumNUSigmaPoints);
    
}

Double_t NeutrinoCrossSectionModel::funcNCCrossSec(Double_t *energy, Double_t *par)
{ 
    //Energy in GeV
    //Cross-section in cm^2
    Int_t isAntiNeutrino=int(par[0]+0.1);
//     if(isAntiNeutrino)
// 	return (2.29e-36)*TMath::Power(energy[0],0.363);
//     return (2.31e-36)*TMath::Power(energy[0],0.363);
    if(isAntiNeutrino)
	return UtilityFuncs::linearInterpolation(energy[0],fNubarEnergy,fNubarSigmaNC,fNumNUSigmaPoints);
    return UtilityFuncs::linearInterpolation(energy[0],fNuEnergy,fNuSigmaNC,fNumNUSigmaPoints);
}

Double_t NeutrinoCrossSectionModel::funcTotCrossSec(Double_t *energy, Double_t *par)
{
    //Energy in GeV
    //Cross-section in cm^2
    Int_t isAntiNeutrino=int(par[0]+0.1);
//     if(isAntiNeutrino)
// 	return (7.80e-36)*TMath::Power(energy[0],0.363);
//     return (7.84e-36)*TMath::Power(energy[0],0.363);
    if(isAntiNeutrino)
	return UtilityFuncs::linearInterpolation(energy[0],fNubarEnergy,fNubarSigmaTot,fNumNUSigmaPoints);
    return UtilityFuncs::linearInterpolation(energy[0],fNuEnergy,fNuSigmaTot,fNumNUSigmaPoints);
}

Double_t NeutrinoCrossSectionModel::funcNuCCInteractionLength(Double_t *energy, Double_t *par)
{
    Double_t density=par[1];
    return (1e-2)/(funcCCCrossSec(energy,par)*Avogadro*density);
}

Double_t NeutrinoCrossSectionModel::funcNuNCInteractionLength(Double_t *energy, Double_t *par)
{
    Double_t density=par[1];
    return (1e-2)/(funcNCCrossSec(energy,par)*Avogadro*density);
}


Double_t NeutrinoCrossSectionModel::funcNuTotInteractionLength(Double_t *energy, Double_t *par)
{

    //There are factors of A which cancel from the top and bottom, as
    // the interaction length is usually A/(N_A * rho * sigma), but
    // sigma is A * funcCCCrossSec
    Double_t density=par[1];
    return (1*1e-2)/(funcTotCrossSec(energy,par)*Avogadro*density);
}


Double_t NeutrinoCrossSectionModel::funcNuCCInteractionLengthInCM(Double_t *energy, Double_t *par)
{
    //There are factors of A which cancel from the top and bottom, as
    // the interaction length is usually A/(N_A * rho * sigma), but
    // sigma is A * funcCCCrossSec
    Double_t density=par[1];
    return 1./(funcCCCrossSec(energy,par)*Avogadro*density);
}

Double_t NeutrinoCrossSectionModel::funcNuNCInteractionLengthInCM(Double_t *energy, Double_t *par)
{
    //There are factors of A which cancel from the top and bottom, as
    // the interaction length is usually A/(N_A * rho * sigma), but
    // sigma is A * funcNCCrossSec
    Double_t density=par[1];
    return 1./(funcNCCrossSec(energy,par)*Avogadro*density);
}


Double_t NeutrinoCrossSectionModel::funcNuTotInteractionLengthInCM(Double_t *energy, Double_t *par)
{
    //There are factors of A which cancel from the top and bottom, as
    // the interaction length is usually A/(N_A * rho * sigma), but
    // sigma is A * funcTotCrossSec
    Double_t density=par[1];
    return 1./(funcTotCrossSec(energy,par)*Avogadro*density);
}

