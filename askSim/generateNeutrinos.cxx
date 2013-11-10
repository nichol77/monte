#include <iostream>
#include <fstream>
#include <string>

using namespace std;

#include "NeutrinoEvent.h"
#include "NeutrinoGenerator.h"
#include "NeuGenRunSummary.h"
#include "AskConventions.h"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"



int main(int argc, char **argv) {

   Long64_t numRolls=10000000;
   Int_t runId=100;
   AskCons::WhichFlux_t whichFlux=AskCons::kStandardGzk;
   Int_t setRandomSeed=0;
   Double_t energyExponent=18;
   Double_t maxDistToIce=1000;
   Int_t doZCut=1;
   Double_t zCut=4e6;


   char outputName[FILENAME_MAX];
   
   //Will fix this at some point
   if(argc>=3) {
      runId=atoi(argv[1]);
      numRolls=atoll(argv[2]);
   }
   if(argc>=4) {
      int tempInt=atoi(argv[3]);
      if(tempInt>=1 && tempInt<=2)
	 whichFlux=(AskCons::WhichFlux_t)tempInt;
   }
   if(argc>=5) {
      setRandomSeed=atoi(argv[4]);
   }
   if(argc>=6) {
      energyExponent=atof(argv[5]);
   }
      

   std::cout << "Run " << runId << "\n";
   std::cout << "Rolling " << numRolls << " rolls\n";
   std::cout << "Which Flux " << whichFlux << "\n";
   std::cout << "Setting Random Seed " << setRandomSeed << "\n";
   if(whichFlux==AskCons::kMonoEnergetic) {
         std::cout << "Energy Exponent " << energyExponent << "\n";
   }

   char *outputDir=getenv("NEUGEN_OUTPUT_DIR");
   if(!outputDir)
      sprintf(outputName,"rootFiles/nuInt%d.root",runId);
   else
      sprintf(outputName,"%s/nuInt%d.root",outputDir,runId);
   
   std::cout << "Output file: " << outputName << std::endl;


    //Output file stuff
    NeutrinoEvent *nuEvent=0; 
    NeuGenRunSummary *nuRun=0;
    //    sprintf(outputName,"rootFiles/nuInt%d.root",runId);
    TFile *outFile = new TFile(outputName,"RECREATE");
    TTree *nuIntTree = new TTree("nuIntTree","Tree of Neutrino Interactions");
    nuIntTree->Branch("event","NeutrinoEvent",&nuEvent);
    TTree *runSumTree = new TTree("runSumTree","Run Summary Tree");
    runSumTree->Branch("run","NeuGenRunSummary",&nuRun);

    NeutrinoGenerator *neugen= new NeutrinoGenerator(runId,whichFlux,maxDistToIce,doZCut,zCut);
    if(whichFlux==AskCons::kMonoEnergetic) 
       neugen->setFixedEnergy(TMath::Power(10,energyExponent));
    if(setRandomSeed)
       neugen->setRandomSeed(0);
      

    Long64_t starEvery=numRolls/30;
    if(starEvery==0)starEvery++;
    

    for(Long64_t i=0;i<numRolls;i++) {
       //Need to have some sort of steering file to decide how many neutrinos to 
       // generate, of what energy, using what volume, etc, etc, etc
       if(i%starEvery==0) 
	  std::cerr << "*";



       nuEvent=neugen->getNextNeutrino();
       if(nuEvent) {       
	  nuIntTree->Fill();
	  delete nuEvent;
	  nuEvent=0;
       }
    }
    std::cerr << "\n";
    nuRun=neugen->getRunSummary();
    runSumTree->Fill();
    nuIntTree->AutoSave();
    runSumTree->AutoSave();
    outFile->Close();
}
