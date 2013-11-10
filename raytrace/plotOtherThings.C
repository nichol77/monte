

void plotOtherThings() 
{
  //  TF1 *nfunc = new TF1("nfunc",straightNFunc,0,130,2);
  //  nfunc->SetParameters(1.78,(1.35-1.78)/130.);
  //  nfunc->Draw();


  Double_t startN=1.78;
  Double_t startDepth=130;

  //  Double_t startN=1.77669;
  //  Double_t startDepth=129;

  //  Double_t startN=1.51538;
  //  Double_t startDepth=50;
  Double_t surfaceTheta=1*TMath::DegToRad();
  Double_t missTheta=90*TMath::DegToRad();
  Double_t nextTheta=45*TMath::DegToRad();

  while(1) {
    Double_t xSurface=simpleGetXFromY(nextTheta,startDepth,startN);
    if(TMath::IsNaN(xSurface)) {
      if(nextTheta<missTheta) {
	missTheta=nextTheta;
      }
    }
    else {
      if(nextTheta>surfaceTheta) {
	surfaceTheta=nextTheta;
      }
    }
    if(TMath::Abs(surfaceTheta-missTheta)<1e-6) {
      break;
    }
    nextTheta=0.5*(surfaceTheta+missTheta);
  }
  std::cout << "Miss theta: " << missTheta*TMath::RadToDeg() << "\n";
  std::cout << "Surface theta: " << surfaceTheta*TMath::RadToDeg() << "\n";
    
    



  //  Double_t startN=1.68077;
  //  Double_t startDepth=100;
  
  // TCanvas *can2 = new TCanvas("can2","can2");
  // can2->Divide(1,2);
  // can2->cd(1);
  // gPad->DrawFrame(-400,-150,1000,150);
  //  Double_t theta=63.9665*TMath::DegToRad();
  Double_t theta=135*TMath::DegToRad();
  // //  theta=TMath::Pi()-theta;
  // //Double_t theta=TMath::PiOver2()+0.05;

  Double_t xSurface=simpleGetXFromY(theta,startDepth,startN);
  Double_t xBottom=simpleGetXFromY(theta,startDepth-130,startN);
  Double_t xBottom2=otherSimpleGetXFromY(theta,startDepth-130,startN);
    
  std::cout << "xSurface: " << xSurface << "\n";
  std::cout << "xBottom: " << xBottom << "\n";
  std::cout << "xBottom2: " << xBottom2 << "\n";

  TF1 *fittydZdR = new TF1("fittydZdR",dZdRvsR,-1000,1000,3);
  fittydZdR->SetParameters(startN,(1.35-1.78)/130,theta);
  std::cout << "dZdR: " << fittydZdR->Eval(xBottom) << "\n";
  //  return;
  // TF1 *fittyZ = new TF1("fittyZ",newCoshBased,-1000,1000,3);
  // fittyZ->SetLineStyle(1);
  // fittyZ->SetLineWidth(1);
  // fittyZ->SetNpx(1000);
  // fittyZ->SetRange(-1000,1000);
  // fittyZ->SetParameters(1.68077,(1.35-1.78)/130,theta);
  // fittyZ->DrawCopy("same");
  // can2->cd(2);
  
  // gPad->DrawFrame(-1000,-2,1000,2);
  // fittydZdR->SetLineStyle(1);
  // fittydZdR->SetLineWidth(1);
  // fittydZdR->SetNpx(1000);
  // fittydZdR->SetRange(-1000,1000);
  // 
  // fittydZdR->DrawCopy("same");
  // return;


  TCanvas *can = new TCanvas("can","can");
  can->Divide(1,3);
  can->cd(1);
  gPad->DrawFrame(0,-150,1000,150);
  TF1 *fitty = new TF1("fitty",newCoshBased,0,1000,3);
  fitty->SetLineStyle(1);
  fitty->SetLineWidth(1);
  fitty->SetNpx(1000);




  


  int counter=0;
  //for(double angle=-1*TMath::Pi()/4;angle<0;angle+=TMath::Pi()/90) {
  //   for(double angle=TMath::Pi()/4;angle<TMath::Pi();angle+=TMath::Pi()/16) {
  for(double angle=TMath::Pi()/8;angle<TMath::Pi();angle+=TMath::Pi()/30) {
    //    cout << 180*angle/TMath::Pi() << "\t" << simpleGetXFromY(angle,130) << endl;
    Double_t xSurface=simpleGetXFromY(angle,startDepth,startN);
    Double_t xBottom=simpleGetXFromY(angle,startDepth-130,startN);
    Double_t xBottom2=otherSimpleGetXFromY(angle,startDepth-130,startN);   
    std::cout << "Here: " << TMath::RadToDeg()*angle << "\t" << xBottom << "\t" << xBottom2 << "\t" << xSurface << "\n";
    if(xBottom>0 && (xBottom<xBottom2 || xBottom2<0)) 
      xBottom2=xBottom;
    fitty->SetRange(0,1000);
    fitty->SetParameters(startN,(1.35-1.78)/130,angle);
    if(!TMath::IsNaN(xBottom2) && xBottom2>0) 
      fitty->SetRange(0,xBottom2);
    if(!TMath::IsNaN(xSurface) && xSurface>0)
      fitty->SetRange(0,xSurface);
    //    std::cout << "Here2: " << angle << "\t" << fitty->Eval(300) << "\n";

    //    fitty->SetParameters(1.35,(1.78-1.35)/130,angle);
    fitty->SetLineColor(40+counter);
    counter++;
    fitty->DrawCopy("same");
  }

  can->cd(2);
  gPad->DrawFrame(-0,-10,1000,10);
  TF1 *fitty2 = new TF1("fitty2",dZdRvsR,0,1000,3);
  fitty2->SetLineStyle(1);
  fitty2->SetLineWidth(1);
  counter=0;
  for(double angle=TMath::Pi()/8;angle<TMath::Pi();angle+=TMath::Pi()/30) {
    //  for(double angle=TMath::Pi()/8;angle<TMath::Pi()/2;angle+=TMath::Pi()/30) {
    //    cout << 180*angle/TMath::Pi() << "\t" << simpleGetXFromY(angle,130) << endl;
    Double_t degAngle= TMath::RadToDeg()*angle;
    Double_t xSurface=simpleGetXFromY(angle,startDepth,startN);
    Double_t xBottom=simpleGetXFromY(angle,startDepth-130,startN);
    Double_t xBottom2=otherSimpleGetXFromY(angle,startDepth-130,startN);   
    if(xBottom>0 && (xBottom<xBottom2 || xBottom2<0)) 
      xBottom2=xBottom;
    fitty2->SetRange(0,1000);
    fitty2->SetParameters(startN,(1.35-1.78)/130,angle);
    if(!TMath::IsNaN(xBottom2) && xBottom2>0) 
      fitty2->SetRange(0,xBottom2);
    if(!TMath::IsNaN(xSurface) && xSurface>0)
      fitty2->SetRange(0,xSurface);
    //    fitty2->SetParameters(1.35,(1.78-1.35)/130,angle);
    char title[180];
    sprintf(title,"dZdR #theta=%f",degAngle);
    fitty2->SetTitle(title);
    fitty2->SetLineColor(40+counter);
    counter++;
    fitty2->DrawCopy("same");
  }

  can->cd(3);
  gPad->DrawFrame(0,-1000,1000,1000);
  TF1 *fitty3 = new TF1("fitty3",timeInFirn,0,1000,3);
  fitty3->SetLineStyle(1);
  fitty3->SetLineWidth(1);
  counter=0;
  for(double angle=TMath::Pi()/8;angle<TMath::Pi();angle+=TMath::Pi()/30) {
    //    cout << 180*angle/TMath::Pi() << "\t" << simpleGetXFromY(angle,130) << endl;
    Double_t xSurface=simpleGetXFromY(angle,startDepth,startN);
    fitty3->SetRange(0,1000);
    fitty3->SetParameters(startN,(1.35-1.78)/130,angle);
    if(!TMath::IsNaN(xSurface))
      fitty3->SetRange(0,xSurface);
    //    fitty3->SetParameters(1.35,(1.78-1.35)/130,angle);
    fitty3->SetLineColor(40+counter);
    counter++;
    fitty3->DrawCopy("same");
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

Double_t  simpleGetXFromY(Double_t theta0, Double_t y, Double_t inititalN){  
  Double_t par[3]={inititalN,(1.35-1.78)/130.,theta0};
  return newCoshBasedInverse(&y,par);
}

Double_t  otherSimpleGetXFromY(Double_t theta0, Double_t y, Double_t inititalN){  
  Double_t par[3]={inititalN,(1.35-1.78)/130.,theta0};
  return newCoshBasedInverseOtherRoot(&y,par);
}

Double_t dZdRvsR(Double_t *x, Double_t *par) {
  Double_t n0=par[0];
  Double_t m=1*par[1];
  Double_t theta0=par[2];
  Double_t thisX=x[0];
  Double_t direction=1;
  if(theta0>TMath::PiOver2()) {
    thisX*=-1;
    direction=-1;
  }
  Double_t k= n0 *TMath::Sin(theta0);
  Double_t xbar=-1*k*TMath::ACosH(n0/k)/m; 
  return direction*TMath::SinH((m/k)*(thisX-xbar));
}

Double_t newCoshBased(Double_t *x, Double_t *par) {
  Double_t n0=par[0];
  Double_t m=1*par[1];
  Double_t theta0=par[2];
  Double_t thisX=x[0];
  Double_t thisTheta=TMath::Abs(theta0);
  Double_t direction=1;
  if(thisTheta>TMath::PiOver2()) {
    thisX*=-1;
  }
  Double_t k= n0 *TMath::Sin(thisTheta);
  Double_t xbar=-1*k*TMath::ACosH(n0/k)/m; 
  Double_t z=-1*(n0-k*TMath::CosH((m/k)*(thisX-xbar)))/m;
  //  cout << k << "\t" << xbar << "\t" << z << endl;
  return z;
}

Double_t newCoshBasedInverse(Double_t *z, Double_t *par) {
  Double_t n0=par[0];
  Double_t m=1*par[1];
  Double_t theta0=par[2];
  Double_t direction=1;
  if(TMath::Abs(theta0)>TMath::PiOver2())
    direction=-1;
  Double_t k= n0 *TMath::Sin(theta0);
  Double_t xbar=-1*k*TMath::ACosH(n0/k)/m; 
  Double_t x=xbar + (k/m)*TMath::ACosH((n0+(m*z[0]))/k);
  //  cout << "newCoshBasedInverse: " << n0 << "\t" << theta0 << "\t" << k << "\t" << xbar << "\t" << z << endl;
  return (x*direction);
}


Double_t newCoshBasedInverseOtherRoot(Double_t *z, Double_t *par) {
  Double_t n0=par[0];
  Double_t m=1*par[1];
  Double_t theta0=par[2];
  Double_t direction=1;
  if(TMath::Abs(theta0)>TMath::PiOver2())
    direction=-1;
  Double_t k= n0 *TMath::Sin(theta0);
  Double_t xbar=-1*k*TMath::ACosH(n0/k)/m; 
  Double_t x=xbar - (k/m)*TMath::ACosH((n0+(m*z[0]))/k);
  //  cout << k << "\t" << xbar << "\t" << z << endl;
  return x*direction;
}

Double_t timeInFirn(Double_t *x, Double_t *par) {
  Double_t n0=par[0];
  Double_t m=1*par[1];
  Double_t theta0=TMath::Abs(par[2]);
  Double_t thisX=x[0];
  Double_t posNeg=1;
  if(theta0>TMath::PiOver2()) {
    thisX*=-1;
    posNeg=-1;
  }
      
  Double_t k= n0 *TMath::Sin(theta0);
  Double_t xbar=-1*k*TMath::ACosH(n0/k)/m; 
  
  Double_t timeVal=(k/(4*m))*(2*x[0]*m +k*TMath::SinH(2*m*(thisX-xbar)/k) + k*TMath::SinH(2*m*xbar/k));
  return posNeg*timeVal;
}
