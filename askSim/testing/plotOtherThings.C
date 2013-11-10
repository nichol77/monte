

void plotOtherThings() 
{
  //  TF1 *nfunc = new TF1("nfunc",straightNFunc,0,130,2);
  //  nfunc->SetParameters(1.78,(1.35-1.78)/130.);
  //  nfunc->Draw();


  TCanvas *can = new TCanvas("can","can");
  can->Divide(1,2);
  can->cd(1);
  gPad->DrawFrame(0,0,1000,150);
  TF1 *fitty = new TF1("fitty",newCoshBased,0,1000,3);
  fitty->SetLineStyle(1);
  fitty->SetLineWidth(1);

  for(double angle=TMath::Pi()/8;angle<TMath::Pi()/2;angle+=TMath::Pi()/90) {
    //    cout << 180*angle/TMath::Pi() << "\t" << simpleGetXFromY(angle,130) << endl;
    Double_t xSurface=simpleGetXFromY(angle,130);
    fitty->SetRange(0,1000);
    fitty->SetParameters(1.78,(1.35-1.78)/130,angle);
    if(!TMath::IsNaN(xSurface))
      fitty->SetRange(0,xSurface);
    //    fitty->SetParameters(1.35,(1.78-1.35)/130,angle);
    fitty->DrawCopy("same");
  }
  can->cd(2);
  gPad->DrawFrame(0,0,1000,1000);
  TF1 *fitty2 = new TF1("fitty2",timeInFirn,0,1000,3);
  fitty2->SetLineStyle(1);
  fitty2->SetLineWidth(1);

  for(double angle=TMath::Pi()/8;angle<TMath::Pi()/2;angle+=TMath::Pi()/90) {
    //    cout << 180*angle/TMath::Pi() << "\t" << simpleGetXFromY(angle,130) << endl;
    Double_t xSurface=simpleGetXFromY(angle,130);
    fitty2->SetRange(0,1000);
    fitty2->SetParameters(1.78,(1.35-1.78)/130,angle);
    if(!TMath::IsNaN(xSurface))
      fitty2->SetRange(0,xSurface);
    //    fitty2->SetParameters(1.35,(1.78-1.35)/130,angle);
    fitty2->DrawCopy("same");
  }

}

Double_t simpleN(Double_t depth) {
  Double_t par[2]={1.78,(1.35-1.78)/130.};
  return straightNFunc(&depth,par);
}

Double_t straightNFunc(Double_t *z, Double_t *par){
  Double_t a=par[0];
  Double_t b=par[1];

  if(z[0]<130 ) {
    return a + b*z[0];
  }
  //  else return b;
}

Double_t  simpleGetXFromY(Double_t theta0, Double_t y){  
  Double_t par[3]={1.78,(1.35-1.78)/130.,theta0};
  return newCoshBasedInverse(&y,par);
}


Double_t newCoshBased(Double_t *x, Double_t *par) {
  Double_t n0=par[0];
  Double_t m=1*par[1];
  Double_t theta0=par[2];
  Double_t k= n0 *TMath::Sin(theta0);
  Double_t xbar=-1*k*TMath::ACosH(n0/k)/m; 
  Double_t z=-1*(n0-k*TMath::CosH((m/k)*(x[0]-xbar)))/m;
  //  cout << k << "\t" << xbar << "\t" << z << endl;
  return z;
}

Double_t newCoshBasedInverse(Double_t *z, Double_t *par) {
  Double_t n0=par[0];
  Double_t m=1*par[1];
  Double_t theta0=par[2];
  Double_t k= n0 *TMath::Sin(theta0);
  Double_t xbar=-1*k*TMath::ACosH(n0/k)/m; 
  Double_t x=xbar + (k/m)*TMath::ACosH((n0+(m*z[0]))/k);
  //  cout << k << "\t" << xbar << "\t" << z << endl;
  return x;
}

Double_t timeInFirn(Double_t *x, Double_t *par) {
  Double_t n0=par[0];
  Double_t m=1*par[1];
  Double_t theta0=par[2];
  Double_t k= n0 *TMath::Sin(theta0);
  Double_t xbar=-1*k*TMath::ACosH(n0/k)/m; 
  
  Double_t timeVal=(k/(4*m))*(2*x[0]*m +k*TMath::SinH(2*m*(x[0]-xbar)/k) + k*TMath::SinH(2*m*xbar/k));
  return timeVal;
}
