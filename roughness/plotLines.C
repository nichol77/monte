

void plotLines() {
  TVector2 intPoint(0,-1000);
  TVector2 balloonPos(600000,37000);
  
  Double_t xPoints[3]={intPoint.X(),0,balloonPos.X()};
  Double_t yPoints[3]={intPoint.Y(),0,balloonPos.Y()};

  TCanvas *can = new TCanvas("can","can",600,600);
  TH1F *framey = can->DrawFrame(-10000,-2000,650000,40000);
 
  Double_t nAir=1;
  Double_t nIce=1.78;
  Double_t speedOfLight=3e8;
  
  TH1F *histTime = new TH1F("histTime","histTime",1000,1.8e6,1.9e6);
  TH2F *histTime2d = new TH2F("histTime2d","histTime2d",100,-10000,10000,1000,1.8e6,1.9e6);

  for(int dist=-10000;dist<10000;dist+=100) {
    TVector2 surfPos(dist,0);
    Double_t time=TMath::Sqrt(TMath::Power(surfPos.X()-intPoint.X(),2) + TMath::Power(surfPos.Y()-intPoint.Y(),2))*nIce;
    time+=TMath::Sqrt(TMath::Power(balloonPos.X()-surfPos.X(),2) + TMath::Power(balloonPos.Y()-surfPos.Y(),2))*nAir;
    xPoints[1]=surfPos.X();
    yPoints[1]=surfPos.Y();
    TGraph *gr = new TGraph(3,xPoints,yPoints);
    gr->Draw("l");

    //    cout << 3*time << "\n";
    histTime->Fill(3*time);
    histTime2d->Fill(dist,3*time);

  }
  TCanvas *canTime = new TCanvas();
  canTime->Divide(1,2);
  canTime->cd(1);
  histTime->Draw("");
  canTime->cd(2);
  histTime2d->Draw("");



}
