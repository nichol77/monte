
void testNuIntFile() {
   gSystem->Load("libAskRay.so");
   testNuIntFile(100);
}



void testNuIntFile(int run) {
   char fileName[180];
   sprintf(fileName,"rootFiles/nuInt%d.root",run);
   TFile *fp = new TFile(fileName);
   TTree *nuIntTree = (TTree*) fp->Get("nuIntTree");
   
   NeutrinoEvent *nuEv=0;
   nuIntTree->Branch("event",&nuEv);
   
   cout << nuIntTree->GetEntries() << endl;

   Int_t numInts=0;
   Double_t xVals[1000000];
   Double_t yVals[1000000];
   Double_t zVals[1000000];
   

   for(int i=0;i<nuIntTree->GetEntries();i++) {
      //      cout << nuEv << endl;
      nuIntTree->GetEntry(i);

   //    xVals[numInts]=nuEv->fSurfPos[0];
//       yVals[numInts]=nuEv->fSurfPos[1];
//       zVals[numInts]=nuEv->fSurfPos[2];
//       numInts++;


      TClonesArray *theInts = nuEv->getInteractions();
      //      cout << theInts.GetLast() << endl;
      for(Int_t j=0;j<theInts->GetLast()+1;j++){	 
	 InteractionInfo *intInfo = (InteractionInfo*) theInts->At(j);
	 //       	 cout << j << "\t" << intInfo->intLocation[0] << endl;
	 xVals[numInts]=intInfo->intLocation[0];
	 yVals[numInts]=intInfo->intLocation[1];
	 zVals[numInts]=intInfo->intLocation[2];
	 numInts++;
      }


   }
   TCanvas *can = new TCanvas("can","can");
   can->Divide(1,3);
   can->cd(1);
   TGraph *gr = new TGraph(numInts,xVals,yVals);
   gr->Draw("ap");
   can->cd(2);
   TGraph *gr2 = new TGraph(numInts,xVals,zVals);
   gr2->Draw("ap");
   can->cd(3);
   TGraph *gr3 = new TGraph(numInts,yVals,zVals);
   gr3->Draw("ap");
   
   
}
