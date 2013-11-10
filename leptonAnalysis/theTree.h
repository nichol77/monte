//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan  4 11:32:08 2005 by ROOT version 4.00/06
// from TTree theTree/Energy Loss Tree
// found on file: muonFile1e9Ice.root
//////////////////////////////////////////////////////////

#ifndef theTree_h
#define theTree_h

#include <iostream>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

using namespace std;

class theTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   Int_t           particleNumber;
   Double_t        isAtau;
   Double_t        A;
   Double_t        Z;
   Double_t        material;
   Double_t        startEnergy;
   Double_t        nu_cut;
   Double_t        totalLength;
   Double_t        totalEnergyLoss;
   Int_t           stepNumber;
   Double_t        stepLength;
   Double_t        stepCEL;
   Double_t        stepSEL;
   Double_t        stepFinalEnergy;
   Double_t        stepIntEnergy;
   Double_t        stepStartEnergy;
   Int_t           stepIntType;
   Double_t        stepLbarMeanFree;
   Double_t        stepLbarDecay;
   Double_t        stepLbarWeak;
   Double_t        stepLApproxMeanFree;
   Double_t        stepLApproxDecay;
   Double_t        stepLApproxWeak;

   // List of branches
   TBranch        *b_particleNumber;   //!
   TBranch        *b_isAtau;   //!
   TBranch        *b_A;   //!
   TBranch        *b_Z;   //!
   TBranch        *b_material;   //!
   TBranch        *b_startEnergy;   //!
   TBranch        *b_nu_cut;   //!
   TBranch        *b_totalLength;   //!
   TBranch        *b_totalEnergyLoss;   //!
   TBranch        *b_stepNumber;   //!
   TBranch        *b_stepLength;   //!
   TBranch        *b_stepCEL;   //!
   TBranch        *b_stepSEL;   //!
   TBranch        *b_stepFinalEnergy;   //!
   TBranch        *b_stepIntEnergy;   //!
   TBranch        *b_stepStartEnergy;   //!
   TBranch        *b_stepIntType;   //!
   TBranch        *b_stepLbarMeanFree;   //!
   TBranch        *b_stepLbarDecay;   //!
   TBranch        *b_stepLbarWeak;   //!
   TBranch        *b_stepLApproxMeanFree;   //!
   TBranch        *b_stepLApproxDecay;   //!
   TBranch        *b_stepLApproxWeak;   //!

   theTree(TTree *tree=0);
   theTree(char *filename);
   ~theTree();
//   Int_t  Cut(Int_t entry);
   Int_t  GetEntry(Int_t entry);
   Int_t  LoadTree(Int_t entry);
   void   Init(TTree *tree);
   void   Loop(Double_t selCut=1e9);
   void   FillHistos();
   Bool_t Notify();
   void   Show(Int_t entry = -1);
};

#endif

#ifdef theTree_cxx
theTree::theTree(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("muonFile1e9Ice.root");
      if (!f) {
         f = new TFile("muonFile1e9Ice.root");
      }
      tree = (TTree*)gDirectory->Get("theTree");

   }
   Init(tree);
}


theTree::theTree(char *filename)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
    TTree *tree=0;
   if (tree == 0) {
       TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
      if (!f) {
         f = new TFile(filename);
      }
      tree = (TTree*)gDirectory->Get("theTree");
      if(!tree) {
	  cout << "Couldn't get theTree from: " << filename << endl;
	  return;
      }
   }
   Init(tree);
}

theTree::~theTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t theTree::GetEntry(Int_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Int_t theTree::LoadTree(Int_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Int_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->IsA() != TChain::Class()) return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void theTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses of the tree
   // will be set. It is normaly not necessary to make changes to the
   // generated code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running with PROOF.

   // Set branch addresses
   if (tree == 0) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("particleNumber",&particleNumber);
   fChain->SetBranchAddress("isAtau",&isAtau);
   fChain->SetBranchAddress("A",&A);
   fChain->SetBranchAddress("Z",&Z);
   fChain->SetBranchAddress("material",&material);
   fChain->SetBranchAddress("startEnergy",&startEnergy);
   fChain->SetBranchAddress("nu_cut",&nu_cut);
   fChain->SetBranchAddress("totalLength",&totalLength);
   fChain->SetBranchAddress("totalEnergyLoss",&totalEnergyLoss);
   fChain->SetBranchAddress("stepNumber",&stepNumber);
   fChain->SetBranchAddress("stepLength",&stepLength);
   fChain->SetBranchAddress("stepCEL",&stepCEL);
   fChain->SetBranchAddress("stepSEL",&stepSEL);
   fChain->SetBranchAddress("stepFinalEnergy",&stepFinalEnergy);
   fChain->SetBranchAddress("stepIntEnergy",&stepIntEnergy);
   fChain->SetBranchAddress("stepStartEnergy",&stepStartEnergy);
   fChain->SetBranchAddress("stepIntType",&stepIntType);
   fChain->SetBranchAddress("stepLbarMeanFree",&stepLbarMeanFree);
   fChain->SetBranchAddress("stepLbarDecay",&stepLbarDecay);
   fChain->SetBranchAddress("stepLbarWeak",&stepLbarWeak);
   fChain->SetBranchAddress("stepLApproxMeanFree",&stepLApproxMeanFree);
   fChain->SetBranchAddress("stepLApproxDecay",&stepLApproxDecay);
   fChain->SetBranchAddress("stepLApproxWeak",&stepLApproxWeak);
   Notify();
}

Bool_t theTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. Typically here the branch pointers
   // will be retrieved. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed.

   // Get branch pointers
   b_particleNumber = fChain->GetBranch("particleNumber");
   b_isAtau = fChain->GetBranch("isAtau");
   b_A = fChain->GetBranch("A");
   b_Z = fChain->GetBranch("Z");
   b_material = fChain->GetBranch("material");
   b_startEnergy = fChain->GetBranch("startEnergy");
   b_nu_cut = fChain->GetBranch("nu_cut");
   b_totalLength = fChain->GetBranch("totalLength");
   b_totalEnergyLoss = fChain->GetBranch("totalEnergyLoss");
   b_stepNumber = fChain->GetBranch("stepNumber");
   b_stepLength = fChain->GetBranch("stepLength");
   b_stepCEL = fChain->GetBranch("stepCEL");
   b_stepSEL = fChain->GetBranch("stepSEL");
   b_stepFinalEnergy = fChain->GetBranch("stepFinalEnergy");
   b_stepIntEnergy = fChain->GetBranch("stepIntEnergy");
   b_stepStartEnergy = fChain->GetBranch("stepStartEnergy");
   b_stepIntType = fChain->GetBranch("stepIntType");
   b_stepLbarMeanFree = fChain->GetBranch("stepLbarMeanFree");
   b_stepLbarDecay = fChain->GetBranch("stepLbarDecay");
   b_stepLbarWeak = fChain->GetBranch("stepLbarWeak");
   b_stepLApproxMeanFree = fChain->GetBranch("stepLApproxMeanFree");
   b_stepLApproxDecay = fChain->GetBranch("stepLApproxDecay");
   b_stepLApproxWeak = fChain->GetBranch("stepLApproxWeak");

   return kTRUE;
}

void theTree::Show(Int_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
/* Int_t theTree::Cut(Int_t entry) */
/* { */
/* // This function may be called from Loop. */
/* // returns  1 if entry is accepted. */
/* // returns -1 otherwise. */
/*    return 1; */
/* } */
#endif // #ifdef theTree_cxx
