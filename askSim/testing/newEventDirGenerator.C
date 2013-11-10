
void newEventDirGenerator() {
   gSystem->Load("libAskRay.so");

   double xVals[10000],yVals[10000],zVals[10000];
   
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



   NeutrinoGenerator *genny = new NeutrinoGenerator(202,25000);

   TVector3 surfPos,intDir;
   int numPoints=0;
   for(int i=0;i<5;i++) {
      genny->pickPointAndDirectionOnEarth(surfPos,intDir);

      xVals[i]=surfPos.X();
      yVals[i]=surfPos.Y();
      zVals[i]=surfPos.Z();

      //      if(!isValid) continue;
      //      cout << surfPos.Mag() << endl;
      //      cout  << TMath::RadToDeg()*AskGeom::getTheta(surfPos) << "\t" <<  TMath::RadToDeg()*AskGeom::getPhi(surfPos) << endl;

      TVector3 radDir=surfPos.Unit();
      //      intDir*=-1;


      Double_t chordLength=TMath::Abs(2*AskGeom::R_EARTH*TMath::Cos(radDir.Angle(intDir)));

      TVector3 oppPoint=surfPos+chordLength*intDir;
      TVector3 otherPoint=surfPos-chordLength*intDir;

      cout << surfPos.Mag() << "\t" << oppPoint.Mag() << "\t" << otherPoint.Mag() << endl;

      liney->SetLineColor(9);     
      liney->DrawLine(surfPos.Y(),surfPos.Z(),oppPoint.Y(),oppPoint.Z());
      numPoints++;
   }
   TGraph *gr = new TGraph(numPoints,yVals,zVals);
   gr->Draw("p");

}
