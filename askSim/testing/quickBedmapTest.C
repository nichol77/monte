//#include "../AskConventions.h"

#define NBINS 500

void quickBedmapTest() {
    gSystem->Load("libAskRay.so");
    BedmapTable *iceTable = new BedmapTable (AskCons::kIceThickness);
//     BedmapTable *groundBedTable = new BedmapTable (AskCons::kGroundBed);
//     cout << "Here" << endl;
//     BedmapTable *waterDepthTable = new BedmapTable (AskCons::kWaterDepth);
//     BedmapTable *bathTable = new BedmapTable (AskCons::kBathymetry);
//     BedmapTable *surfaceTable = new BedmapTable (AskCons::kSurfaceHeight);
//     BedmapTable *bedElevationTable = new BedmapTable (AskCons::kBedElevation);


    TH2F *iceHist = new TH2F("iceHist","Ice Thickness",NBINS,-3e6,3e6,NBINS,-3e6,3e6);
//     TH2F *groundBedHist = new TH2F("groundBedHist","Ground Bed",NBINS,-3e6,3e6,NBINS,-3e6,3e6);
//     TH2F *waterDepthHist = new TH2F("waterDepthHist","Water Depth",NBINS,-3e6,3e6,NBINS,-3e6,3e6);
//     TH2F *bathHist = new TH2F("bathHist","Bathymetry",NBINS,-3e6,3e6,NBINS,-3e6,3e6);
//     TH2F *surfaceHist = new TH2F("surfaceHist","Surface Height",NBINS,-3e6,3e6,NBINS,-3e6,3e6);
//     TH2F *bedElevationHist = new TH2F("bedElevHist","Bed Elevation",NBINS,-3e6,3e6,NBINS,-3e6,3e6);
//     TH2F *compositeHist = new TH2F("compositeHist","Composite",NBINS,-3e6,3e6,NBINS,-3e6,3e6);
    int goodFlag;
    for(int binx=1;binx<=iceHist->GetNbinsX();binx++) {
	Double_t x=iceHist->GetXaxis()->GetBinCenter(binx);
//	cout << binx << endl;
	for(int biny=1;biny<=iceHist->GetNbinsY();biny++) {
	    Double_t y=iceHist->GetYaxis()->GetBinCenter(biny);

	    Double_t ice=0;
	    Double_t water=0;
	    Double_t elev=0;
	    Double_t bath=0;

	    Double_t value=iceTable->getValue(x,y,goodFlag);
	    if(goodFlag) {
		iceHist->SetBinContent(binx,biny,value);
		ice=value;
	    }
//	    value=groundBedTable->getValue(x,y,goodFlag);
// 	    if(goodFlag) groundBedHist->SetBinContent(binx,biny,value);
// 	    value=waterDepthTable->getValue(x,y,goodFlag);
// 	    if(goodFlag) {
// 		waterDepthHist->SetBinContent(binx,biny,value);
// 		water=value;
// 	    }
// 	    value=bathTable->getValue(x,y,goodFlag);
// 	    if(goodFlag) {
// 		bathHist->SetBinContent(binx,biny,value);
// 		bath=value;
// 	    }
// 	    value=surfaceTable->getValue(x,y,goodFlag);
// 	    if(goodFlag) surfaceHist->SetBinContent(binx,biny,value);
// 	    value=bedElevationTable->getValue(x,y,goodFlag);
// 	    if(goodFlag) {
// 		bedElevationHist->SetBinContent(binx,biny,value);
// 		elev=value;
// 	    }
// 	    compositeHist->SetBinContent(binx,biny,ice+water+elev);
	}
    }
//     gStyle->SetPadRightMargin(0.2);
//     TCanvas *can = new TCanvas("can","can");
//     can->Divide(2,3);
//     can->cd(1);
//     iceHist->Draw("colz");
//     can->cd(2);
//     groundBedHist->Draw("colz");
//     can->cd(3);
//     waterDepthHist->Draw("colz");
//     can->cd(4);
//     bathHist->Draw("colz");
//     can->cd(5);
//     surfaceHist->Draw("colz");
//     can->cd(6);
//     bedElevationHist->Draw("colz");
    
    TCanvas *can2 = new TCanvas("can2","can2");
    iceHist->Draw("colz");

}
