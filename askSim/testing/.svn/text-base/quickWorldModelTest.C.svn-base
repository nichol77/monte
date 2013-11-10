#include "RyansConventions.h"


void quickWorldModelTest() {
    gSystem->Load("lib/libPropagation.so");
//    BedmapTable *iceTable = new BedmapTable (RyansCons::kBathymetry);
    WorldModel *myWorld = new WorldModel();
    TH2F *iceHist = new TH2F("iceHist","Ice Thickness",1000,-3e6,3e6,1000,-3e6,3e6);
    TH2F *iceHist2 = new TH2F("iceHist2","Ice Thickness (Bedmap)",1000,-3e6,3e6,1000,-3e6,3e6);
    TH2F *surfaceHist = new TH2F("surfaceHist","Surface",1000,-3e6,3e6,1000,-3e6,3e6);
//    int goodFlag;
    for(int binx=1;binx<=iceHist->GetNbinsX();binx++) {
	Double_t x=iceHist->GetXaxis()->GetBinCenter(binx);
	for(int biny=1;biny<=iceHist->GetNbinsY();biny++) {
	    Double_t y=iceHist->GetYaxis()->GetBinCenter(biny);
	    Double_t value=myWorld->getCrust2Data(x,y,RyansCons::kThicknessOfIce);
	    iceHist->SetBinContent(binx,biny,value);
	    value=myWorld->getIceThickness(x,y);
	    iceHist2->SetBinContent(binx,biny,value);
	    value=myWorld->getSurface(x,y);
	    surfaceHist->SetBinContent(binx,biny,value/1e6);
	    if(value>1e7 || value<1e6) cout << x << "\t" << y << "\t" << value << endl;
	}
    }
    TCanvas *can = new TCanvas("can","can");
    can->Divide(1,2);
    can->cd(1);
    iceHist->Draw("colz");
    can->cd(2);
    iceHist2->Draw("colz");
	
    TCanvas *can2 = new TCanvas("can2","can2");
    surfaceHist->Draw("colz");
    
}
