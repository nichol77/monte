


Double_t getC(Double_t *par, Double_t D) {
  Double_t n0=par[0];  //1.35
  Double_t m=par[1];  //(1.78-1.35)/140
  Double_t z0=par[2]; //antenna position
  
  Double_t val=(D/m);
  val*=TMath::ACosH((n0+m*z0)/D);
  return val;
}

Double_t strangeLeft(Double_t *x, Double_t *par) {
  Double_t n0=par[0];  //1.35
  Double_t m=par[1];  //(1.78-1.35)/140
  Double_t z0=par[2]; //antenna position
  Double_t x1=par[3]; // pulser position
  Double_t z1=par[4]; // pulser position
  Double_t D=x[0];
  //  Double_t C=getC(par,D);
  return (n0+m*z1)/D;
}

Double_t strangeRight(Double_t *x, Double_t *par) {
  Double_t n0=par[0];  //1.35
  Double_t m=par[1];  //(1.78-1.35)/140
  Double_t z0=par[2]; //antenna position
  Double_t x1=par[3]; // pulser position
  Double_t z1=par[4]; // pulser position
  Double_t D=x[0];
  Double_t C=getC(par,D);
  
  Double_t val=TMath::CosH(m*(x1+C)/D);
  //  std::cout << D << "\t" << C << "\t" << val <<"\n";
  return val;
}


Double_t strangeZ(Double_t *x, Double_t *par) {
  Double_t n0=par[0];  //1.35
  Double_t m=par[1];  //(1.78-1.35)/140
  Double_t z0=par[2]; //antenna position
  Double_t x1=par[3]; // pulser position
  Double_t D=x[0];
  Double_t C=getC(par,D);
  Double_t z1=(D/m)*TMath::CosH(m*(x1+C)/D)-n0/m;
  //  std::cout << D << "\t" << C << "\t" << z1 <<"\n";
  return z1;  
}



void solveNumerically()
{
  
  Double_t n0=1.35;
  Double_t m=-1*(1.78-1.35)/140;
  Double_t z0=-100;
  Double_t dMax=n0+m*z0;


  TCanvas *can = new TCanvas("can","can",600,400);
  // TF1 *lefty = new TF1("lefty",strangeLeft,0,dMax,5);
  // lefty->SetParameters(n0,m,z0,40,-30);
  // lefty->SetLineColor(8);
  // lefty->SetLineStyle(1);
  // lefty->SetNpx(1000);
  // lefty->Draw("");

  // TF1 *righty = new TF1("righty",strangeRight,0,dMax,5);
  // righty->SetParameters(n0,m,z0,40,-30);
  // righty->SetLineColor(kViolet);
  // righty->SetLineStyle(1);
  // righty->SetNpx(1000);
  // righty->Draw("same");


  TH1F*framey = can->DrawFrame(0,-50,dMax*1.1,0);
  TF1 *endy = new TF1("endy",strangeZ,0.02,dMax,4);
  endy->SetParameters(n0,m,z0,200);
  endy->SetLineColor(kViolet);
  endy->SetLineStyle(1);
  endy->SetNpx(1000);
  endy->Draw("same");

}
