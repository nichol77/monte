#define theTree_cxx
#include "theTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <iostream>
using namespace std;

void theTree::Loop(Double_t selCut)
{
//   In a ROOT session, you can do:
//      Root > .L theTree.C
//      Root > theTree t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Int_t nentries = Int_t(fChain->GetEntriesFast());

   Int_t nbytes = 0, nb = 0;
   Int_t lastParticleNumber=0;
   TH1F *histNum = new TH1F("histNum","histNum",200,-0.5,199.5);

   Int_t numEventsAboveCut=0;
   for (Int_t jentry=0; jentry<nentries;jentry++) {
      Int_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(particleNumber!=lastParticleNumber) {
	  histNum->Fill(numEventsAboveCut);
	  numEventsAboveCut=0;
      }
      
      // if (Cut(ientry) < 0) continue;
      if(stepSEL>selCut) numEventsAboveCut++;
      lastParticleNumber=particleNumber;
   }
   histNum->Fill(numEventsAboveCut);
	  
}

void theTree::FillHistos()
{
   if (fChain == 0) return;

   Int_t nentries = Int_t(fChain->GetEntriesFast());

   Int_t nbytes = 0, nb = 0;
   Int_t lastParticleNumber=0;
   Double_t energyCut[8]={1e8,3e8,1e9,3e9,1e10,3e10,1e11,3e11};
   Double_t c_light=299792458;
 
   TH1F *histNums[8];
   TH1F *histDists[8];
   TH1F *histTimes[8];
   TH1F *histIntTypes[8];
   char histName[80];
   char histTitle[80];

   for(int i=0;i<8;i++) {
       sprintf(histTitle,"E > %1.0e GeV",energyCut[i]);
       sprintf(histName,"histNum%d",i);
       histNums[i]= new TH1F(histName,histTitle,200,-0.5,199.5);
       sprintf(histName,"histDist%d",i);
       histDists[i]= new TH1F(histName,histTitle,250,0,50);
       sprintf(histName,"histTime%d",i);
       histTimes[i]= new TH1F(histName,histTitle,250,0,100);
       sprintf(histName,"histIntType%d",i);
       histIntTypes[i]= new TH1F(histName,histTitle,5,-1.5,3.5);
   }

   Int_t numEventsAboveCut[8]={0};
   Double_t lastTotalLengthForSEL[8]={0};
   for (Int_t jentry=0; jentry<nentries;jentry++) {
      Int_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(particleNumber!=lastParticleNumber) {
	  for(int i=0;i<8;i++) {
	      histNums[i]->Fill(numEventsAboveCut[i]);
	      numEventsAboveCut[i]=0;
	      lastTotalLengthForSEL[i]=0;
	  }
      }
      
      // if (Cut(ientry) < 0) continue;
      for(int i=0;i<8;i++) {
	  if(stepSEL>energyCut[i]) {
	      numEventsAboveCut[i]++;
// 	      if(numEventsAboveCut[i]==1) {
// 		  cout << totalLength << "\t" << lastTotalLengthForSEL[i] << endl;
// 	      }
	      histDists[i]->Fill((totalLength-lastTotalLengthForSEL[i])/1e5);
	      histTimes[i]->Fill(1e4*(totalLength-lastTotalLengthForSEL[i])/c_light);
	      histIntTypes[i]->Fill(stepIntType);

	      lastTotalLengthForSEL[i]=totalLength;
	  }
      }
      lastParticleNumber=particleNumber;
   }
   for(int i=0;i<8;i++) 
       histNums[i]->Fill(numEventsAboveCut[i]);
	  
}
