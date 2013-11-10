#include "TGraph.h"
#ifndef ROOT_Rtypes
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Rtypes.h"
#endif
#endif

#include "AskConventions.h"

void pairPlot();
void plotPairBetaAsFunctionOfEnergy();
void bremPlot();
void plotBremBetaAsFunctionOfEnergy();
void plotBremAsFunctionOfNu(); 
//Double_t funcPlot(Double_t *x, Double_t *par);

Double_t funcPhotonudSdnu(Double_t *nu, Double_t *par);
Double_t funcPhotodSdnu(Double_t *nu, Double_t *par);
Double_t funcPairnudSdnu(Double_t *nu, Double_t *par);
Double_t funcSlowPairdSdnu(Double_t *nu, Double_t *par);
Double_t funcPairdSdnu(Double_t *nu, Double_t *par);
Double_t funcBremnudSdnu(Double_t *nu, Double_t *par);
Double_t funcBremdSdnu(Double_t *nu, Double_t *par);
Double_t funcKnockOnnudSdnu(Double_t *nu, Double_t *par);
Double_t funcKnockOndSdnu(Double_t *nu, Double_t *par);


Double_t getBetaForBrem(Double_t *energy, Double_t *par);

Double_t getdEdXCombined(Double_t *energy, Double_t *par);
Double_t getBetaCombined(Double_t *energy, Double_t *par);
Double_t getBetaForPair(Double_t *energy, Double_t *par);
Double_t pairdSdnu(Double_t energy, Double_t nu, Double_t leptonMass,Double_t Z);
Double_t pairRhoPart(Double_t rho, Double_t energy, Double_t nu, Double_t leptonMass, Double_t Z);
Double_t getPairNuMin(Double_t energy);
Double_t getPairNuMax(Double_t energy,Double_t leptonMass, Double_t Z);
Double_t funcPairRhoPart(Double_t *rho, Double_t *par);
Double_t getSigmaGammaN(Double_t nu, Double_t energy);
Double_t softPhotodSdnu(Double_t energy, Double_t nu, Double_t Z, Double_t A, Double_t leptonMass);
Double_t getGz(Double_t z, Double_t Z);
Double_t getHnu(Double_t nu);
Double_t funcsoftnudSdnu(Double_t *nu, Double_t *par);
//Double_t interpolate(Double_t x, Double_t *iX, Double_t *iY, Int_t n);
Double_t funchardnudSdnu(Double_t *nu, Double_t *par);
Double_t hardPhotodSdnu(Double_t energy, Double_t nu, Int_t isATau);
//void getCoeffs();
Double_t linearInterpolation(Double_t x, Double_t *iX, Double_t *iY, Int_t numPoints);
Double_t getBetaFromFulldSdnu(Double_t *energy, Double_t *par);
Double_t bremdSdnu(Double_t energy, Double_t nu, Double_t leptonMass,Double_t Z);
Double_t getBremNuMax(Double_t energy,Double_t leptonMass, Double_t Z);
Double_t knockOndSdNu(Double_t energy, Double_t nu, Double_t leptonMass, Double_t Z);
Double_t getKnockOnNuMax(Double_t energy, Double_t leptonMass);
Double_t funcnuDelta_egamma(Double_t *nu, Double_t *par);
Double_t getDelta_egamma(Double_t energy, Double_t nu, Double_t leptonMass);
Double_t getdeltaX(Double_t beta, Double_t gamma, Double_t X0, Double_t X1, Double_t a, Double_t m, Double_t C);
Double_t ionizationdEdx(Double_t energy, Double_t leptonMass, Int_t material);
Double_t funcIonizationdEdx(Double_t *energy, Double_t *par);
Double_t funcOneOverBeta(Double_t *energy, Double_t *par);
Double_t funcOneOverBetaSlow(Double_t *energy, Double_t *par);
Double_t funcIntegralOneOverBetaSlow(Double_t *energy, Double_t *par);
Double_t funcIntegralOneOverBeta(Double_t *energy, Double_t *par);
Double_t funcDecayRange(Double_t *energy, Double_t *par);
Double_t funcOneOverBetaTau(Double_t *energy, Double_t *par); 
Double_t funcOneOverBetaMuon(Double_t *energy, Double_t *par); 
Double_t funcContinuousEnergyLoss(Double_t *energy, Double_t *par);
Double_t funcMeanFreePath(Double_t *energy, Double_t *par);
Double_t funcOneOverCELMFP(Double_t *energy, Double_t *par);
Double_t funcOneOverCEL(Double_t *energy, Double_t *par);
Double_t funcOneOverCELQuick(Double_t *energy, Double_t *par);

Double_t getBremCrossAboveNuCut(Double_t energy, Double_t isATau, Double_t A, Double_t Z, Double_t nu_cut);
Double_t getPairCrossAboveNuCut(Double_t energy, Double_t isATau, Double_t A, Double_t Z, Double_t nu_cut);
Double_t getPhotoCrossAboveNuCut(Double_t energy, Double_t isATau, Double_t A, Double_t Z, Double_t nu_cut);
Double_t getKnockOnCrossAboveNuCut(Double_t energy, Double_t isATau, Double_t A, Double_t Z, Double_t nu_cut);
Double_t funcCrossAboveNuCut(Double_t *energy, Double_t *par);

//Double_t funcEnergyLossRange(Double_t *energy, Double_t *par);

AskCons::InteractionType_t whichInteraction(Double_t randNum, Double_t energy, Double_t isATau, Double_t material, Double_t nu_cut);

TGraph *computeRightIntegral(TF1 *inputFunc, Double_t maxVal, Double_t minVal, Int_t numPoints);
TGraph *getRightIntegral(TF1 *inputFunc, Double_t energy, Double_t lengthScale, Int_t numPoints);
Double_t getNewEnergy(Double_t startEnergy, Double_t length, Double_t isATau, Double_t material, Double_t nuCut, Int_t quickMethod=1);
Double_t guessNewEnergyStraightLineMethod(TF1 *inputFunc, Double_t energy, Double_t lengthScale);
Double_t funcDecayExponential(Double_t *x, Double_t *par);
Double_t funcDecayRangeCmWithRho(Double_t *energy, Double_t *par);

Double_t newFuncOneOverCEL(Double_t *energy, Double_t *par);
Double_t newFuncCEL(Double_t *energy, Double_t *par);
inline int getNewEnergyIndex(double energy);
