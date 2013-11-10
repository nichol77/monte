TGraph *drawPaths(Double_t startN, Double_t deepN, Double_t startDepth, Double_t firnDepth, Double_t gradN, char *plotTitle);

Int_t MyPalette[256];

void makeAbbyPlot()
{
  //  TF1 *nfunc = new TF1("nfunc",straightNFunc,0,130,2);
  //  nfunc->SetParameters(1.78,(1.35-1.78)/130.);
  //  nfunc->Draw();
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
   
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor color;
  color.InitializeColors();
  Int_t fi=color.CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  for (int i=0;i<255;i++) {
    //     std::cout << i << "\t" << TColor::GetColor(i) << "\n";
    MyPalette[i] = fi+i;
    //     MyPalette[i] = TColor::GetColor(i);
  }

  TGraph *gr[2];

  TCanvas *can = new TCanvas("can","can",800,800);
  can->Divide(1,2);
  can->cd(1);
  {


    Double_t startN=1.35;
    Double_t deepN=1.78;
    Double_t startDepth=5;
    Double_t firnDepth=100;
    Double_t gradN=(startN-deepN)/firnDepth;
    gr[0]=drawPaths(startN,deepN,startDepth,firnDepth,gradN,"");
    
  
  }
  can->cd(2);
  {
    
    Double_t startN=1.35;
    Double_t deepN=1.78;
    Double_t startDepth=120;
    Double_t firnDepth=100;
    Double_t gradN=(startN-deepN)/firnDepth;
    gr[1]=drawPaths(startN,deepN,startDepth,firnDepth,gradN,"");
 
  }

  TCanvas *can2 = new TCanvas("can2","can2");
  //Now some smoothing
  Double_t *xVals = gr[1]->GetX();
  Double_t *yVals = gr[1]->GetY();
  Double_t *yVals2 = gr[0]->GetY();
  for(int i=0;i<gr[1]->GetN();i++) {
    if(xVals[i]<60) yVals2[i]=0;
    if(xVals[i]<225) yVals[i]=0;
    if(xVals[i]>=420) yVals[i]=-100;
  }
  
  gr[0]->SetLineColor(getNiceColour(1));
  gr[0]->Draw("al");
  gr[1]->SetLineColor(getNiceColour(2));
  gr[1]->Draw("l");

  Double_t volFrac[101]={0};

  Double_t vSurf=0;
  Double_t vDrill=0;
  

  for(int i=0;i<gr[1]->GetN();i++) {
    Double_t aSurf=TMath::Abs(TMath::Pi()*xVals[i]*(1000+yVals2[i]));
    Double_t aDrill=TMath::Abs(TMath::Pi()*xVals[i]*(1000+yVals[i]));
    vDrill+=aDrill*10;
    vSurf+=aSurf*10;
    std::cout << xVals[i] << "\t" << vSurf << "\t" << vDrill << "\n";
    if(vDrill>0) {
      volFrac[i]=vSurf/vDrill;
    }
    else {
      volFrac[i]=1;
    }
  }
  TGraph *grFrac = new TGraph(101,xVals,volFrac);
  grFrac->SetLineStyle(1);
  grFrac->SetLineWidth(3);
  grFrac->SetLineColor(kRed);
  grFrac->SetTitle("Deep vs Surface");
  grFrac->Draw("al");
  grFrac->GetXaxis()->SetTitle("Horizontal Distance (m)");
  grFrac->GetXaxis()->SetNoExponent(1);
  grFrac->GetYaxis()->SetTitle("V_{surface}/V_{deep}");
  grFrac->GetYaxis()->SetNoExponent(1); 
  grFrac->GetYaxis()->SetRangeUser(0,1.1);
}





TGraph *drawPaths(Double_t startN, Double_t deepN, Double_t startDepth, Double_t firnDepth, Double_t gradN, char *plotTitle) {
 
  Double_t xValues[101]={0};
  Double_t maxDepth[101]={0};
  for(int i=0;i<101;i++) {
    xValues[i]=i*10;
    maxDepth[i]=-1000;
  }

  TH1F *framey=gPad->DrawFrame(0,-1000,1000,10);
  framey->SetTitle(plotTitle);
  framey->SetXTitle("Distance (m)");
  framey->GetXaxis()->SetNoExponent();
  framey->SetYTitle("Depth (m)");
  framey->GetYaxis()->SetNoExponent();
  TF1 *fitty = new TF1("fitty",newCoshBased,0,1000,5);
  TF1 *fitStraight2 = new TF1("fitStraight","pol1",0,1000);
  TF1 *fitStraight1 = new TF1("fitStraight1","pol1",0,1000);
  fitty->SetLineStyle(1);
  fitty->SetLineWidth(1);
  fitty->SetNpx(200);
  fitStraight2->SetLineStyle(1);
  fitStraight2->SetLineWidth(1);
  fitStraight2->SetNpx(10000);
  fitStraight1->SetLineStyle(1);
  fitStraight1->SetLineWidth(1);
  fitStraight1->SetNpx(10000);

  int counter=0;
  //for(double angle=-1*TMath::Pi()/4;angle<0;angle+=TMath::Pi()/90) {
  //   for(double angle=TMath::Pi()/4;angle<TMath::Pi();angle+=TMath::Pi()/16) {

  Int_t numLines=180;
  Int_t colourJump=255/numLines;

  for(double angle=TMath::Pi()/90;angle<TMath::Pi();angle+=TMath::Pi()/(numLines)) {

    Bool_t drawStraight1=true;
    Bool_t drawStraight2=true;
    Bool_t drawCurved=true;
    Double_t testDepth=startDepth;
    Double_t testX=0;


    if(startDepth>firnDepth) {
      //Below the firn

      Double_t m=1./TMath::Tan(angle);
      Double_t c=-startDepth;
      fitStraight1->SetParameters(c,m); 
      fitStraight1->SetRange(0,1000);
      if(angle<TMath::PiOver2()) {
	//Going up
	testX=TMath::Tan(angle)*(startDepth-firnDepth);
	testDepth=firnDepth;
	fitStraight1->SetRange(0,testX);
      }
      else {
	drawCurved=false;
	drawStraight2=false;
      }  
      std::cout << "First: " << startDepth << "\t" << firnDepth << "\t" << angle*TMath::RadToDeg() << "\t" << testX << "\t" << testDepth << "\t"  << drawCurved << "\n";
    
    }
    else {
      drawStraight1=false;
    }


    if(drawCurved) {
      Double_t xSurface=simpleGetXFromY(angle,testDepth,startN,gradN);
      Double_t xSurface2=otherSimpleGetXFromY(angle,testDepth,startN,gradN);
      Double_t otherBottom=otherSimpleGetXFromY(angle,testDepth-firnDepth,startN,gradN);
      Double_t otherBottom2=otherSimpleGetXFromY(angle,(testDepth-firnDepth)+0.01,startN,gradN);
      Double_t xBottom=simpleGetXFromY(angle,testDepth-firnDepth,startN,gradN);
      if(xBottom<0) xBottom*=-1;
      Double_t xBottom2=simpleGetXFromY(angle,(testDepth-firnDepth)+0.01,startN,gradN);   
      if(xBottom2<0) xBottom2*=-1;
      std::cout << "Here: " << counter << "\t" << TMath::RadToDeg()*angle << "\t" << xBottom << "\t" << otherBottom <<  "\t" << xSurface << "\t" << xSurface2 << "\n";
      
      if(TMath::IsNaN(xSurface) && otherBottom>xBottom) {
	xBottom=otherBottom;
	xBottom2=otherBottom2;
      }
      
      fitty->SetRange(0,1000);
      fitty->SetParameters(testDepth,testX,startN,gradN,angle);
      if(!TMath::IsNaN(xSurface) && xSurface>0) {
	fitty->SetRange(testX,testX+xSurface);
	drawStraight2=false;
      }
      else if(!TMath::IsNaN(xBottom) && xBottom>0) {	
	  fitty->SetRange(testX,xBottom+testX);
	  fitStraight2->SetRange(xBottom+testX,1000);
	  Double_t m=-0.01/(xBottom-xBottom2);
	  Double_t c=-(1*firnDepth) -m*(xBottom+testX);
	  fitStraight2->SetParameters(c,m);      	
      }
    }
    std::cout << "Here2: " << fitty->GetXmin() << "\t" << fitty->GetXmax() << "\n";
    
    //    fitty->SetParameters(1.35,(1.78-1.35)/firnDepth,angle);
    fitty->SetLineColor(MyPalette[counter]);
    fitStraight2->SetLineColor(MyPalette[counter]);
    fitStraight1->SetLineColor(MyPalette[counter]);
    
    counter+=colourJump;
    char grTitle[180];
    sprintf(grTitle,"Angle %d",(int)angle*TMath::RadToDeg());
    fitty->SetTitle(grTitle);
    fitStraight->SetTitle(grTitle);
    if(drawStraight1) fitStraight1->DrawCopy("same");
    if(drawCurved) fitty->DrawCopy("same");
    if(drawStraight2) fitStraight2->DrawCopy("same");

    for(int i=0;i<101;i++) {
      if(drawStraight1 && fitStraight1->GetXmin()<=xValues[i] && fitStraight1->GetXmax()>=xValues[i]) {
	Double_t testVal=fitStraight1->Eval(xValues[i]);
	if(testVal>maxDepth[i]) 
	  maxDepth[i]=testVal;
      }
      if(drawStraight2 && fitStraight2->GetXmin()<=xValues[i] && fitStraight2->GetXmax()>=xValues[i]) {
	Double_t testVal=fitStraight2->Eval(xValues[i]);
	if(testVal>maxDepth[i]) 
	  maxDepth[i]=testVal;
      }
      if(drawCurved && fitty->GetXmin()<=xValues[i] && fitty->GetXmax()>=xValues[i]) {
	Double_t testVal=fitty->Eval(xValues[i]);
	if(testVal>maxDepth[i]) 
	  maxDepth[i]=testVal;
      }      
    }
  }
  TGraph *gr = new TGraph(101,xValues,maxDepth);
  return gr;
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

Double_t  simpleGetXFromY(Double_t theta0, Double_t y, Double_t inititalN, Double_t gradN){  
  Double_t par[3]={inititalN,gradN,theta0};
  return newCoshBasedInverse(&y,par);
}

Double_t  otherSimpleGetXFromY(Double_t theta0, Double_t y, Double_t inititalN,Double_t gradN){  
  Double_t par[3]={inititalN,gradN,theta0};
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
  //      fitty->SetParameters(testDepth,testX,startN,gradN,angle);
  Double_t z0=par[0];
  Double_t x0=par[1];
  Double_t n0=par[2];
  Double_t m=1*par[3];
  Double_t theta0=par[4];
  Double_t thisX=x[0]-x0;
  Double_t thisTheta=TMath::Abs(theta0);
  Double_t direction=1;
  if(thisTheta>TMath::PiOver2()) {
    thisX*=-1;
  }
  Double_t k= n0 *TMath::Sin(thisTheta);
  Double_t xbar=-1*k*TMath::ACosH(n0/k)/m; 
  Double_t z=(-z0)-1*(n0-k*TMath::CosH((m/k)*(thisX-xbar)))/m;
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
