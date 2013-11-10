

const double rEarth=6378.1e3;

void plotThisOne() {
   gSystem->Load("libAskRay.so");


  TCanvas * can = new TCanvas("can","can");
  TH1F *framey = can->DrawFrame(-7e6,-7e6,+7e6,+7e6);

  TEllipse *elipsey = new TEllipse(0,0,6378.1e3,6378.1e3);
  elipsey->SetLineColor(8);
  elipsey->SetLineWidth(3);			
  elipsey->Draw();//Ellipse(0,0,6378.1e3,6378.1e3);

  Double_t gz[361],gx[361];
  Int_t count=0;
  for(Double_t theta=-180;theta<=180;theta+=1) {
     
     Double_t radius=AskGeom::getGeoidFromTheta(theta*TMath::DegToRad());
     gz[count]=radius*TMath::Cos(theta*TMath::DegToRad());
      gx[count]=radius*TMath::Sin(theta*TMath::DegToRad());      
      
      //      cout << theta << "\t" << radius << "\t" << gz[count] << "\t" << gx[count] << endl;
      count++;
      
  }
  TGraph *geoid = new TGraph(count,gx,gz);
  geoid->Draw("l");
  
  TLine *liney = new TLine();
  liney->SetLineColor(9);
  liney->SetLineWidth(1);
  liney->SetLineStyle(2);

  Double_t point1[3]={-294476, 489656,  6.33461e+06};
  Double_t point2[3]={-297716, 487929,  6.34309e+06};

  liney->DrawLine(point1[1],point1[2],point2[1],point2[2]);

}

