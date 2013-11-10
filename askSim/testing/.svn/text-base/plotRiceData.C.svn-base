
Float_t depths[19]={5,10,15,20,25,30,35,45,55,65,75,85,95,105,115,125,135,145,200};
Float_t ns[19]={1.35,1.38,1.45,1.5,1.45,1.53,1.5,1.55,1.51,1.55,1.55,1.62,1.68,1.73,1.72,1.72,1.75,1.78,1.78};


void plotRiceData() {

  TGraph *gr= new TGraph(19,depths,ns);
  gr->Draw("ap");

  TF1 * curvey = new TF1("curvey",peterFunc,0,200,3);
  curvey->SetParameters(1.78,0.33,30);
  curvey->Draw("same");

  TF1 * liney = new TF1("liney",ryanFunc,0,200,3);
  liney->SetParameters(1.35,1.78,130);
  liney->Draw("same");
}


Double_t peterFunc(Double_t *x, Double_t *par) {
  Double_t a=par[0];
  Double_t b=par[1];
  Double_t c=par[2];

  return a * (1 - b * TMath::Exp(-x[0]/c));
}

Double_t ryanFunc(Double_t *x, Double_t *par){
  Double_t a=par[0];
  Double_t b=par[1];
  Double_t c=par[2];

  if(x[0]<c) {
    return a + (b-a)*x[0]/c;
  }
  else return b;
}
