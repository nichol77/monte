
void testDirectionPicker() {
   gSystem->Load("libAskRay.so");

   
   NeutrinoGenerator *genny = new NeutrinoGenerator(100,25000);

   TH1F *histPhi = new TH1F("histPhi","histPhi",360,0,360);
   TH1F *histTheta = new TH1F("histTheta","histTheta",360,0,360);
   TH1F *histCosTheta = new TH1F("histCosTheta","histCosTheta",360,-1,1);

   Double_t phi,theta;
   for(int i=0;i<10000;i++) {
      //      genny->pickRandomThetaPhiOnSphere(theta,phi);
      genny->pickRandomDowngoingDirection(theta,phi);
      histPhi->Fill(phi*TMath::RadToDeg());
      histTheta->Fill(theta*TMath::RadToDeg());
      histCosTheta->Fill(TMath::Cos(theta));
   }

   TCanvas *can = new TCanvas("can","can");
   can->Divide(1,3);
   can->cd(1);
   histPhi->Draw();
   can->cd(2);
   histTheta->Draw();
   can->cd(3);
   histCosTheta->Draw();

}
