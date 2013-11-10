
TRandom3 fRandom;

const double rEarth=6378.1e3;
void eventDirGenerator() {


  double xVals[10000],yVals[10000],zVals[10000];

  double rV=25000;
  double thetaV,phiV;

  double thetaMom,phiMom;

  int numPoints=0;

  TCanvas * can = new TCanvas("can","can");
  TH1F *framey = can->DrawFrame(-7e6,-7e6,+7e6,+7e6);

  TEllipse *elipsey = new TEllipse(0,0,6378.1e3,6378.1e3);
  elipsey->SetLineColor(8);
  elipsey->SetLineWidth(3);			
  elipsey->Draw();//Ellipse(0,0,6378.1e3,6378.1e3);

  TLine *liney = new TLine();
  liney->SetLineColor(9);
  liney->SetLineWidth(1);
  liney->SetLineStyle(2);

  for(int i=0;i<10000;i++) {
    pickRandomThetaPhiOnSphere(phiV,thetaV);
    pickRandomDowngoingDirection(phiMom,thetaMom);
    double dxp=-1*TMath::Cos(phiMom)*TMath::Sin(thetaMom);
    double dyp=-1*TMath::Sin(phiMom)*TMath::Sin(thetaMom);
    double dzp=-1*TMath::Cos(thetaMom);

    TVector3 veccy(dxp,dyp,dzp);
    veccy.RotateZ(phiV);
    veccy.RotateX(thetaV);

//     double cosphi=TMath::Cos(phiV);
//     double sinphi=TMath::Sin(phiV);
//     double costheta=TMath::Cos(thetaV);
//     double sintheta=TMath::Sin(thetaV);

//     double dx=cosphi*dxp-sinphi*dyp;
//     double dy=costheta*sinphi*dxp + costheta*cosphi*dyp - sintheta*dzp;
//     double dz=sintheta*sinphi*dxp + sintheta*cosphi*dxp + costheta*dzp;

      
//    cout << veccy.X() << "\t" << veccy.Y() << "\t" << veccy.Z() << endl;
    //    cout << dx << "\t" << dy << "\t" << dz << endl;
    //    cout << veccy.Mag() << endl;


    //cout << phiV << "\t" << thetaV << endl;
    //    cout <<rEarth*rEarth << "\t" << rV*rV << "\t" <<  2*rEarth*rV*TMath::Cos(TMath::Pi()-thetaV) << endl;
    double rNew=TMath::Sqrt(rEarth*rEarth + rV*rV - 2*rEarth*rV*TMath::Cos(TMath::Pi()-thetaV));
    double thetaNew=TMath::ASin(TMath::Sin(TMath::Pi()-thetaV)*rV/rNew);
    double phiNew=phiV;
    
//     if(rNew>rEarth) {
//       cout << "Above:\t" <<  endl;
//     }
//     else {
//       cout << "Below:\t" <<  endl;
//     }
  

    double z=rNew*TMath::Cos(thetaNew);
    double x=rNew*TMath::Sin(thetaNew)*TMath::Cos(phiNew);
    double y=rNew*TMath::Sin(thetaNew)*TMath::Sin(phiNew);

    TVector3 intPos(x,y,z);
    TVector3 surfPos;
    
    int isValid=getSurfacePoint(intPos,veccy,surfPos);
    //    cout << isValid << endl;


    double x2=surfPos.X();
    double y2=surfPos.Y();
    double z2=surfPos.Z();



    if(isValid) {
      double x3=surfPos.X()+veccy.X()*1e6;
      double z3=surfPos.Z()+veccy.Z()*1e6;
      liney->DrawLine(x2,z2,x3,z3);
      xVals[numPoints]=surfPos.X();
      yVals[numPoints]=surfPos.Y();
      zVals[numPoints]=surfPos.Z();
      
      numPoints++;
    }
  }
  cout << numPoints << endl;
   TGraph *grxz= new TGraph(numPoints,xVals,yVals);
   grxz->Draw("p");
  



}


void pickRandomThetaPhiOnSphere(Double_t &phi, Double_t &theta) {
  phi=2*TMath::Pi()*fRandom.Rndm();
  theta=TMath::Pi()*fRandom.Rndm();
}

void pickRandomDowngoingDirection(Double_t &phi, Double_t &theta) {
  phi=2*TMath::Pi()*fRandom.Rndm();
  theta=0.5*TMath::Pi()*fRandom.Rndm();
}

int getSurfacePoint(TVector3 intPos, TVector3 &intDir, TVector3 &surfPos) {
  Double_t b=intDir.X()*intPos.X() + intDir.Y()*intPos.Y() + intDir.Z()*intPos.Z();
  Double_t c = intPos.X()*intPos.X() + intPos.Y()*intPos.Y() + intPos.Z()*intPos.Z() - rEarth*rEarth;
  if(b*b < c) 
    return 0;


  Double_t l1=-1*b + TMath::Sqrt(b*b - c);
  Double_t l2=-1*b - TMath::Sqrt(b*b - c);
  
  
  if(intPos.Mag2()> rEarth*rEarth) {
    //Start outside take l1
    //return 0;
    //    cout << "Outside:\t" << l1 << "\t" << l2 << endl;
    intDir*=-1;
    l1*=-1;
    surfPos.SetX(intPos.X() + l1*intDir.X());
    surfPos.SetY(intPos.Y() + l1*intDir.Y());
    surfPos.SetZ(intPos.Z() + l1*intDir.Z());
    //    cout << surfPos.Mag() << endl;
  }
  else {
    //Start inside take l2
    //    cout << "Inside:\t" << l1 << "\t" << l2 <<  endl;
    //    cout << l2 << endl;
    //    return 0;
    surfPos.SetX(intPos.X() + l2*intDir.X());
    surfPos.SetY(intPos.Y() + l2*intDir.Y());
    surfPos.SetZ(intPos.Z() + l2*intDir.Z());
    //    cout << surfPos.Mag() << endl;
  }
  
  return 1;
}
