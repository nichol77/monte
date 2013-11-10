

void plotRayPaths() 
{
  //  TF1 *nfunc = new TF1("nfunc",straightNFunc,0,130,2);
  //  nfunc->SetParameters(1.78,(1.35-1.78)/130.);
  //  nfunc->Draw();
  plotRayPaths(1000);
  
}


void plotRayPaths(Int_t depth) {

  TCanvas *can = new TCanvas("can","can",800,800);
  can->Divide(1,2);
  can->cd(1);
  char histTitle[180];
  sprintf(histTitle,"Ray Tracing -- Interaction Depth %d m",depth);
  TH1F *framey=gPad->DrawFrame(0,-1.1*depth,5000,0,histTitle);
  framey->GetXaxis()->SetTitle("Horizontal Distance (m)");
  framey->GetYaxis()->SetTitle("Depth (m)");


  TF1 *fitty = new TF1("fitty",rayPath,0,5000,4);
  fitty->SetLineStyle(1);
  fitty->SetLineWidth(1);
  fitty->SetNpx(100);

  Double_t intDepth=depth;

  fitty->SetParameters(1.78,(1.35-1.78)/130,TMath::Pi()/4,depth);
  //  fitty->Draw("same");
  for(double angle=TMath::Pi()/90;angle<TMath::Pi()/2;angle+=TMath::Pi()/45) {
    //    cout << 180*angle/TMath::Pi() << "\n";// << simpleGetXFromY(angle,130) << endl;
    Int_t angDeg= 180*angle/TMath::Pi();
    //    if(angDeg<30)
    //      fitty->SetNpx(4000);
    fitty->SetParameters(1.78,(1.35-1.78)/130,angle,depth);
    fitty->DrawCopy("same");
  }
  
  TLine *liney = new TLine();
  liney->SetLineStyle(2);
  liney->SetLineWidth(2);
  liney->SetLineColor(38);
  liney->DrawLine(0,-130,5000,-130);

  can->cd(2);
  TF1 *timey = new TF1("timey",timePath,0,5000,4);
  TH1F *framey=gPad->DrawFrame(0,0,5000,10000,histTitle);
  framey->GetXaxis()->SetTitle("Horizontal Distance (m)");
  framey->GetYaxis()->SetTitle("Optical Distance (m)");
  
  timey->SetParameters(1.78,(1.35-1.78)/130,TMath::Pi()/4,depth);
  timey->SetLineWidth(1);
  timey->SetLineStyle(1);
  //  timey->Draw("same");
  for(double angle=TMath::Pi()/90;angle<TMath::Pi()/2;angle+=TMath::Pi()/45) {
    timey->SetParameters(1.78,(1.35-1.78)/130,angle,depth);
    timey->DrawCopy("same");
  }
  

}

Double_t simpleYAsFunctionOfX(Double_t intDepth, Double_t intTheta, Double_t x)
{
  Double_t par[4]={1.78,(1.35-1.78)/130,intTheta,intDepth};
  return rayPath(&x,par);
}



Double_t rayPath(Double_t *x, Double_t *par)
{
  //Parameters are
  //0 == n_ice
  //1 == m firn (1.35-1.78)/130
  //2 == initial theta
  //3 == initial depth
  Double_t depth=par[3];
  Double_t theta=par[2];

  Double_t y=0;
  Double_t xf=0;
  if(depth>130) 
    xf=TMath::Tan(theta)*(depth-130);  
  if(x[0]<xf) {
    //Still in straight line section
    y=x[0]/TMath::Tan(theta);
    y-=depth;
  }
  else {
    Double_t xSurface=simpleGetXFromY(theta,130);        
    Double_t xprime=x[0]-xf;
    //    cout << xSurface << "\t" << xprime << "\n";
    if(TMath::IsNaN(xSurface)) {
      //Doesn't reach surface
      //still need to check for when it leaves the firn

      Double_t testy=0;
      Double_t outOfFirnX=newCoshBasedInverseOtherRoot(&testy,par);
      if(xprime<outOfFirnX) {
	//We're in the Firn
	y=newCoshBased(&xprime,par);
	y-=130;
      }
      else {
	//We've been bent out the firn
	y=-130-(x[0]-(xf+outOfFirnX))/TMath::Tan(theta);       
      }
    }
    else if(xprime<xSurface) {
      //Hasn't Reached Surface
      y=newCoshBased(&xprime,par);
      y-=130;
    }
    else {
      //Reached Surface
      //So we're going to try and add a reflection in here
      if(xprime<2*xSurface) {
	//Bounced off the surface and are on the way back down through the Firn
	Double_t equivX=(2*xSurface)-xprime;	
	y=newCoshBased(&equivX,par);
	y-=130;
      }
      else {
	//Bounced off the surface and have left the Firn	
	y=-130-(x[0]-(xf+2*xSurface))/TMath::Tan(theta);      
      }
      //      y=0;
    }
  }
  //  cout << depth << "\t" << theta << "\t" << xf << "\t"<<x[0] << "\n";
  return y;
}



Double_t timePath(Double_t *x, Double_t *par)
{
  //Parameters are
  //0 == n_ice
  //1 == m firn (1.35-1.78)/130
  //2 == initial theta
  //3 == initial depth
  Double_t n_ice=par[0];
  Double_t depth=par[3];
  Double_t theta=par[2];

  Double_t t=0;
  Double_t tf=0;
  Double_t xf=0;
  if(depth>130) {
    xf=TMath::Tan(theta)*(depth-130);  
    tf=n_ice*xf/TMath::Sin(theta);
  }
  if(x[0]<xf) {
    //Still in straight line section
    t=n_ice*x[0]/TMath::Sin(theta);
    //    cout << tf << "\t" << t << "\t" << x[0] << "\t" << xf<< "\n";
  }
  else {
    Double_t xSurface=simpleGetXFromY(theta,130);        
    Double_t xprime=x[0]-xf;
    //    cout << xSurface << "\t" << xprime << "\n";
    if(TMath::IsNaN(xSurface)) {
      //Doesn't reach surface
      //still need to check for when it leaves the firn

      Double_t testy=0;
      Double_t outOfFirnX=newCoshBasedInverseOtherRoot(&testy,par);
      if(xprime<outOfFirnX) {
	//We're in the Firn
	t=tf+timeInFirn(&xprime,par);
      }
      else {
	//We've been bent out the firn
	t=tf+timeInFirn(&outOfFirnX,par)+n_ice*(x[0]-(xf+outOfFirnX))/TMath::Sin(theta);       
      }
    }
    else if(xprime<xSurface) {
      //Hasn't Reached Surface
      t=tf+timeInFirn(&xprime,par);
    }
    else {
      //Reached Surface
      //So we're going to try and add a reflection in here
      Double_t tSurface=timeInFirn(&xSurface,par);
      if(xprime<2*xSurface) {
	//Bounced off the surface and are on the way back down through the Firn
	Double_t equivX=(2*xSurface)-xprime;	
	t=tf+tSurface+(tSurface-timeInFirn(&equivX,par));
      }
      else {
	//Bounced off the surface and have left the Firn	
	t=tf+2*tSurface+n_ice*(x[0]-(xf+2*xSurface))/TMath::Sin(theta);      
      }
      //      y=0;
    }
  }
  //  cout << depth << "\t" << theta << "\t" << xf << "\t"<<x[0] << "\n";
  return t;
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

//Double_t  simpleGetXFromY(Double_t theta0, Double_t y){  
//  Double_t par[3]={1.78,(1.35-1.78)/130.,theta0};
//  return exactMethod(&y,par);
//}

  
Double_t exactMethod(Double_t *y, Double_t *par) {
  Double_t n0=par[0];
  Double_t m=par[1];
  Double_t theta0=par[2];
  Double_t z=n0 +m *y[0] ;
  Double_t k= n0 *TMath::Sin(theta0);
  
  
  Double_t firstLog=TMath::Log(z + TMath::Sqrt(z*z-k*k));
  Double_t secondLog=TMath::Log(n0 + TMath::Sqrt(n0*n0-k*k));
  Double_t scaleFactor=(k/m);
  //  cout << scaleFactor << "\t" << firstLog << "\t" << secondLog << "\n";
  return scaleFactor*(firstLog-secondLog);
}

Double_t exactMethodOtherRoot(Double_t *y, Double_t *par) {
  Double_t n0=par[0];
  Double_t m=par[1];
  Double_t theta0=par[2];
  Double_t z=n0 +m *y[0] ;
  Double_t k= n0 *TMath::Sin(theta0);
  
  
  Double_t firstLog=TMath::Log(z - TMath::Sqrt(z*z-k*k));
  Double_t secondLog=TMath::Log(n0 + TMath::Sqrt(n0*n0-k*k));
  Double_t scaleFactor=(k/m);
  //  cout << scaleFactor << "\t" << firstLog << "\t" << secondLog << "\n";
  return scaleFactor*(firstLog-secondLog);
}

Double_t exactInverse(Double_t *x, Double_t *par) {
  Double_t n0=par[0];
  Double_t m=par[1];
  Double_t theta0=par[2];
  Double_t k= n0 *TMath::Sin(theta0);
  //  Double_t expA=TMath::Exp((m*x[0]/k) + TMath::Log(n0+ TMath::Sqrt(n0*n0-k*k)));
  //  Double_t z=(0.5/expA)*(expA*expA+k*k);
  Double_t z=n0*TMath::CosH(m*x[0]/k)+TMath::Sqrt(n0*n0 -k*k)*TMath::SinH(m*x[0]/k);
  Double_t y=(z-n0)/m ;
  return y;
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

Double_t newCoshBasedInverseOtherRoot(Double_t *z, Double_t *par) {
  Double_t n0=par[0];
  Double_t m=1*par[1];
  Double_t theta0=par[2];
  Double_t k= n0 *TMath::Sin(theta0);
  Double_t xbar=-1*k*TMath::ACosH(n0/k)/m; 
  Double_t x=xbar - (k/m)*TMath::ACosH((n0+(m*z[0]))/k);
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
