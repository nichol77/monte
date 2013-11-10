
Float_t depths[19]={5,10,15,20,25,30,35,45,55,65,75,85,95,105,115,125,135,145,200};
Float_t ns[19]={1.35,1.38,1.45,1.5,1.45,1.53,1.5,1.55,1.51,1.55,1.55,1.62,1.68,1.73,1.72,1.72,1.75,1.78,1.78};


void plotRiceData() {
  TCanvas *can = new TCanvas("can","can",600,400);
  TGraph *gr= new TGraph(19,depths,ns);
  gr->SetTitle("RICE Measurements");
  gr->Draw("ap");
  gr->GetXaxis()->SetTitle("Depth (m)");
  gr->GetYaxis()->SetTitle("Refractive Index (n)");

  TF1 * curvey = new TF1("curvey",peterFunc,0,200,3);
  curvey->SetParameters(1.325,0.463,0.0140);
  curvey->SetLineStyle(1);
  curvey->SetLineColor(8);
  curvey->Draw("same");

  TF1 * liney = new TF1("liney",ryanFunc,0,200,3);
  liney->SetParameters(1.35,1.78,140);
  liney->SetLineStyle(1);
  liney->SetLineColor(50);
  liney->Draw("same");

  TF1 * poly = new TF1("poly",ryanFunc2,0,200,3);
  poly->SetParameters(1.35,1.78,180);
  poly->SetLineStyle(1);
  poly->SetLineColor(kViolet);
  //  gr->Fit("poly");
  
  poly->Draw("same");
}


Double_t peterFunc(Double_t *x, Double_t *par) {
  Double_t a=par[0];
  Double_t b=par[1];
  Double_t c=par[2];


  Double_t val=(a  + b *( 1- TMath::Exp(-x[0]*c)));
  //  std::cout << x[0] << "\t" << val << "\n";
  return val;
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


Double_t ryanFunc2(Double_t *x, Double_t *par){
  // d is 150`
  //Require that at x=150 the polynomial is 1.78   ---> a*150^2 + b*150 + 1.35 = 1.78
  //Require that at x=0 the polynomial is 1.35 --> c=1.35
  //Require that at x=150 the differential is zero  --> 2a*150 + b = 0   --> b = -300*a
  //Putting that into the top one gives  ---> a*150^2 -a*300*150 = 1.78-1.35
  //Or a = -1*(n_deep-n_top)/d^2.
  //And b = -2*d*a  = -2*(n_deep-n_top)/d


  Double_t c=par[0];  //n_top
  Double_t d=par[2];  // depth when we get to n_deep
  Double_t e=par[1];  //n_deep

  Double_t a=-1*(e-c)/(d*d);
  Double_t b=-2*d*a;

  //  std::cout << a << "\t" << b << "\t" << c << "\t" << d << "\t" << e << "\n";

  if(x[0]<d) {
    return (a*x[0]*x[0] + -2*d*a*x[0] + c);
  }
  else return e;



}
