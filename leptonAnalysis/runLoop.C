char *getParticleName(int isATau);
char *getParticleNameCaps(int isATau);
TCanvas *drawBangGapCanvas(int isATau, char *energy);
TCanvas *drawBangGapTimeCanvas(int isATau, char *energy);
TCanvas *drawBangIntTypeCanvas(int isATau, char *energy);
TCanvas *drawFirstBangGapTimeCanvas(int isATau, char *energy);
TCanvas *drawFirstBangIntTypeCanvas(int isATau, char *energy);
TCanvas *drawLastIntTypeCanvas(int isATau, char *energy);
TCanvas *drawNumBangsCanvas(int isATau, char *energy);

char muonName[]="muon";
char MuonName[]="Muon";
char tauName[]="tau";
char TauName[]="Tau";


void runLoop() {
    int isATau=0;
    char energy[]="1e12";
    char histFilename[80];
    sprintf(histFilename,"newestHist%sFile%sIce.root",getParticleName(isATau),energy);    
    
    gStyle->SetStatH(0.3);
    gStyle->SetStatW(0.3);
    Long_t fileId,fileFlags,fileModTime;
    Long64_t fileSize;
    Int_t retVal=gSystem->GetPathInfo(histFilename,&fileId,&fileSize,&fileFlags,&fileModTime);
    if(retVal) {
	cout << "Need to make histos" << endl;
	char filename[80];	
	gSystem->CompileMacro("analysisTree.C","k");
	sprintf(filename,"newest%sFile%sIce.root",getParticleNameCaps(isATau),energy);
	analysisTree *t =new analysisTree(filename);
	t->FillHistos(histFilename);
    }
    TFile *histFile = new TFile(histFilename);
    TCanvas *canNum=drawNumBangsCanvas(isATau,energy);
    TCanvas *canGap=drawBangGapCanvas(isATau,energy);
    TCanvas *canGapTime=drawBangGapTimeCanvas(isATau,energy);
    TCanvas *canFirstGap=drawFirstBangGapCanvas(isATau,energy);
    TCanvas *canFirstGapTime=drawFirstBangGapTimeCanvas(isATau,energy);
    TCanvas *canIntType=drawBangIntTypeCanvas(isATau,energy);
    TCanvas *canLastIntType=drawLastIntTypeCanvas(isATau,energy);
    
    char outputFilename[180];
    sprintf(outputFilename,"piccies/%sNumBangs%s",
	    getParticleName(isATau),energy);
    canNum->Update();
    canNum->Modified();
    PrintEPSandGIF(canNum,outputFilename);
    sprintf(outputFilename,"piccies/%sGapDistance%s",
	    getParticleName(isATau),energy);
    canGap->Update();
    canGap->Modified();
    PrintEPSandGIF(canGap,outputFilename);
    sprintf(outputFilename,"piccies/%sGapTime%s",
	    getParticleName(isATau),energy);
    canGapTime->Update();
    canGapTime->Modified();
    PrintEPSandGIF(canGapTime,outputFilename);
    sprintf(outputFilename,"piccies/%sFirstGapDistance%s",
	    getParticleName(isATau),energy);
    canFirstGap->Update();
    canFirstGap->Modified();
    PrintEPSandGIF(canFirstGap,outputFilename);
    sprintf(outputFilename,"piccies/%sFirstGapTime%s",
	    getParticleName(isATau),energy);
    canFirstGapTime->Update();
    canFirstGapTime->Modified();
    PrintEPSandGIF(canFirstGapTime,outputFilename);
    sprintf(outputFilename,"piccies/%sInteractionType%s",
	    getParticleName(isATau),energy);
    canIntType->Update();
    canIntType->Modified();
    PrintEPSandGIF(canIntType,outputFilename);
    sprintf(outputFilename,"piccies/%sParticleDestiny%s",
	    getParticleName(isATau),energy);
    canLastIntType->Update();
    canLastIntType->Modified();
    PrintEPSandGIF(canLastIntType,outputFilename);
    

}

TCanvas *drawNumBangsCanvas(int isATau, char *energy) {
    char canTitle[180];
    sprintf(canTitle,"Number of 'Bangs' Per Track -- %ss with Starting Energy %s GeV",getParticleNameCaps(isATau),energy); 
    TCanvas *canNum = new TCanvas("canNum","Bangs per Track",800,800);
    TPaveLabel *pl = new TPaveLabel(0.1,0.96,0.9,0.99,canTitle,"br NDC");
    pl->SetBorderSize(0);
    pl->SetFillColor(0);
    pl->SetFillStyle(0);
    pl->Draw();
    TPad *subCanNum = new TPad("subCanNum","",0,0,1,0.95);
    subCanNum->Draw();
    subCanNum->cd();
    subCanNum->Divide(2,4);
    subCanNum->Update();
    

    TH1F *histNums[8];
    char histName[80];           
    for(int i=0;i<8;i++) {	
	sprintf(histName,"histNum%d",i);
	histNums[i]=(TH1F*) gROOT->FindObject(histName);
	subCanNum->cd(i+1);
	subCanNum->SetTopMargin(0.12);
	histNums[i]->SetLineWidth(2);
	histNums[i]->SetLineColor(getNiceColour(i));	
	histNums[i]->Draw();
    }
         
    for(int i=0;i<8;i++) {
//	cout << "Doing title "<< i+1 << endl;
	subCanNum->cd(i+1);
	sortOutTitle(0.07);
    }
    return canNum;
}


TCanvas *drawBangGapCanvas(int isATau, char *energy) {
    char canTitle[180];
    sprintf(canTitle,"Distance Between 'Bangs' -- %ss with Starting Energy %s GeV",getParticleNameCaps(isATau),energy); 
    TCanvas *canGap = new TCanvas("canGap","Bang Gap",800,800);
    TPaveLabel *pl = new TPaveLabel(0.1,0.96,0.9,0.99,canTitle,"br NDC");
    pl->SetBorderSize(0);
    pl->SetFillColor(0);
    pl->SetFillStyle(0);
    pl->Draw();
    TPad *subCanGap = new TPad("subCanGap","",0,0,1,0.95);
    subCanGap->Draw();
    subCanGap->cd();
    subCanGap->Divide(2,4);
    subCanGap->Update();
    

    TLine *liney = new TLine();
    TH1F *histDists[8];
    char histName[80];           
    for(int i=0;i<8;i++) {	
	sprintf(histName,"histDist%d",i);
	histDists[i]=(TH1F*) gROOT->FindObject(histName);
	subCanGap->cd(i+1);
	histDists[i]->SetLineWidth(2);
	histDists[i]->SetLineColor(getNiceColour(i));
	histDists[i]->GetXaxis()->SetTitle("Distance (km)");
	histDists[i]->Draw();
	gPad->Update();
// 	Double_t yMin=histDists[i]->GetYaxis()->GetBinLowEdge(1);
// 	Double_t yMax= gPad->GetFrame()->GetY2();
	
// 	cout << yMin << "\t" << yMax << endl;
// 	liney->DrawLine(histDists[i]->GetMean(),0,histDists[i]->GetMean(),yMax);
	if(histDists[i]->GetEntries())
	    gPad->SetLogy();
// 	gPad->Modified();
// 	gPad->Update();
// 	Double_t yMin=0;//Or use GetY1 like below
// 	Double_t yMax=TMath::Power(10,gPad->GetFrame()->GetY2());
	
// //	cout << yMin << "\t" << yMax << endl;
// 	liney->DrawLine(histDists[i]->GetMean(),0,histDists[i]->GetMean(),yMax);
	
    }
         
    for(int i=0;i<8;i++) {
//	cout << "Doing title "<< i+1 << endl;
	subCanGap->cd(i+1);
	sortOutTitle(0.07);
    }
    return canGap;
}

TCanvas *drawFirstBangGapCanvas(int isATau, char *energy) {
    char canTitle[180];
    sprintf(canTitle,"Distance Until First 'Bang' -- %ss with Starting Energy %s GeV",getParticleNameCaps(isATau),energy); 
    TCanvas *canFirstGap = new TCanvas("canFirstGap","First Bang Gap",800,800);
    TPaveLabel *pl = new TPaveLabel(0.1,0.96,0.9,0.99,canTitle,"br NDC");
    pl->SetBorderSize(0);
    pl->SetFillColor(0);
    pl->SetFillStyle(0);
    pl->Draw();
    TPad *subCanFirstGap = new TPad("subCanFirstGap","",0,0,1,0.95);
    subCanFirstGap->Draw();
    subCanFirstGap->cd();
    subCanFirstGap->Divide(2,4);
    subCanFirstGap->Update();
    

    TLine *liney = new TLine();
    TH1F *histDists[8];
    char histName[80];           
    for(int i=0;i<8;i++) {	
	sprintf(histName,"histFirstDist%d",i);
	histDists[i]=(TH1F*) gROOT->FindObject(histName);
	subCanFirstGap->cd(i+1);
	histDists[i]->SetLineWidth(2);
	histDists[i]->SetLineColor(getNiceColour(i));
	histDists[i]->GetXaxis()->SetTitle("Distance (km)");
	histDists[i]->Draw();
	gPad->Update();
// 	Double_t yMin=histDists[i]->GetYaxis()->GetBinLowEdge(1);
// 	Double_t yMax= gPad->GetFrame()->GetY2();
	
// 	cout << yMin << "\t" << yMax << endl;
// 	liney->DrawLine(histDists[i]->GetMean(),0,histDists[i]->GetMean(),yMax);
	if(histDists[i]->GetEntries())
	    gPad->SetLogy();
// 	gPad->Modified();
// 	gPad->Update();
// 	Double_t yMin=0;//Or use GetY1 like below
// 	Double_t yMax=TMath::Power(10,gPad->GetFrame()->GetY2());
	
// //	cout << yMin << "\t" << yMax << endl;
// 	liney->DrawLine(histDists[i]->GetMean(),0,histDists[i]->GetMean(),yMax);
	
    }
         
    for(int i=0;i<8;i++) {
//	cout << "Doing title "<< i+1 << endl;
	subCanFirstGap->cd(i+1);
	sortOutTitle(0.07);
    }
    return canFirstGap;
}


TCanvas *drawBangGapTimeCanvas(int isATau, char *energy) {
    char canTitle[180];
    sprintf(canTitle,"Time Between 'Bangs' -- %ss with Starting Energy %s GeV",getParticleNameCaps(isATau),energy); 
    TCanvas *canGapTime = new TCanvas("canGapTime","Bang GapTime",800,800);
    TPaveLabel *pl = new TPaveLabel(0.1,0.96,0.9,0.99,canTitle,"br NDC");
    pl->SetBorderSize(0);
    pl->SetFillColor(0);
    pl->SetFillStyle(0);
    pl->Draw();
    TPad *subCanGapTime = new TPad("subCanGapTime","",0,0,1,0.95);
    subCanGapTime->Draw();
    subCanGapTime->cd();
    subCanGapTime->Divide(2,4);
    subCanGapTime->Update();
    

    TH1F *histTimes[8];
    char histName[80];           
    for(int i=0;i<8;i++) {	
	sprintf(histName,"histTime%d",i);
	histTimes[i]=(TH1F*) gROOT->FindObject(histName);
	subCanGapTime->cd(i+1);
	subCanGapTime->SetTopMargin(0.12);
	histTimes[i]->SetLineWidth(2);
	histTimes[i]->SetLineColor(getNiceColour(i));
	histTimes[i]->GetXaxis()->SetTitle("Time (#mus)");
	histTimes[i]->Draw();
	if(histTimes[i]->GetEntries())
	    gPad->SetLogy();
    }
         
    for(int i=0;i<8;i++) {
//	cout << "Doing title "<< i+1 << endl;
	subCanGapTime->cd(i+1);
	sortOutTitle(0.07);
    }
    return canGapTime;
}


TCanvas *drawFirstBangGapTimeCanvas(int isATau, char *energy) {
    char canTitle[180];
    sprintf(canTitle,"Time Until First 'Bang' -- %ss with Starting Energy %s GeV",getParticleNameCaps(isATau),energy); 
    TCanvas *canFirstGapTime = new TCanvas("canFirstGapTime","Bang FirstGapTime",800,800);
    TPaveLabel *pl = new TPaveLabel(0.1,0.96,0.9,0.99,canTitle,"br NDC");
    pl->SetBorderSize(0);
    pl->SetFillColor(0);
    pl->SetFillStyle(0);
    pl->Draw();
    TPad *subCanFirstGapTime = new TPad("subCanFirstGapTime","",0,0,1,0.95);
    subCanFirstGapTime->Draw();
    subCanFirstGapTime->cd();
    subCanFirstGapTime->Divide(2,4);
    subCanFirstGapTime->Update();
    

    TH1F *histTimes[8];
    char histName[80];           
    for(int i=0;i<8;i++) {	
	sprintf(histName,"histFirstTime%d",i);
	histTimes[i]=(TH1F*) gROOT->FindObject(histName);
	subCanFirstGapTime->cd(i+1);
	subCanFirstGapTime->SetTopMargin(0.12);
	histTimes[i]->SetLineWidth(2);
	histTimes[i]->SetLineColor(getNiceColour(i));
	histTimes[i]->GetXaxis()->SetTitle("Time (#mus)");
	histTimes[i]->Draw();
	if(histTimes[i]->GetEntries())
	    gPad->SetLogy();
    }
         
    for(int i=0;i<8;i++) {
//	cout << "Doing title "<< i+1 << endl;
	subCanFirstGapTime->cd(i+1);
	sortOutTitle(0.07);
    }
    return canFirstGapTime;
}


TCanvas *drawBangIntTypeCanvas(int isATau, char *energy) {
    char canTitle[180];
    sprintf(canTitle,"Interaction Type of 'Bangs' -- %ss with Starting Energy %s GeV",getParticleNameCaps(isATau),energy); 
    TCanvas *canIntType = new TCanvas("canIntType","Bang Interaction Type",800,800);
    TPaveLabel *pl = new TPaveLabel(0.1,0.96,0.9,0.99,canTitle,"br NDC");
    pl->SetBorderSize(0);
    pl->SetFillColor(0);
    pl->SetFillStyle(0);
    pl->Draw();
    TPad *subCanIntType = new TPad("subCanIntType","",0,0,1,0.95);
    subCanIntType->Draw();
    subCanIntType->cd();
    subCanIntType->Divide(2,4);
    subCanIntType->Update();
    
    gStyle->SetOptStat(0);
    TH1F *histIntTypes[8];
    char histName[80];           
    for(int i=0;i<8;i++) {	
	sprintf(histName,"histIntType%d",i);
	histIntTypes[i]=(TH1F*) gROOT->FindObject(histName);
	subCanIntType->cd(i+1);
	gPad->SetTopMargin(0.12);
//	gPad->SetBottomMargin(0.2);
	histIntTypes[i]->SetLineWidth(3);
	histIntTypes[i]->SetLineColor(getNiceColour(i));
//	histIntTypes[i]->GetXaxis()->SetTitle("IntType (#mus)");
	histIntTypes[i]->GetXaxis()->SetBinLabel(2,"Bremsstrahlung");
	histIntTypes[i]->GetXaxis()->SetBinLabel(3,"Pair Production");
	histIntTypes[i]->GetXaxis()->SetBinLabel(4,"Photonuclear");
	histIntTypes[i]->GetXaxis()->SetLabelSize(0.09);
	histIntTypes[i]->GetXaxis()->SetLabelOffset(0.02);
	histIntTypes[i]->Draw();
	if(histIntTypes[i]->GetEntries())
	    gPad->SetLogy();
    }
         
    for(int i=0;i<8;i++) {
//	cout << "Doing title "<< i+1 << endl;
	subCanIntType->cd(i+1);
	sortOutTitle(0.07);
    }
//    gStyle->SetOptStat(1110);
    return canIntType;
}


TCanvas *drawLastIntTypeCanvas(int isATau, char *energy) {
    char histTitle[180];
    sprintf(histTitle,"Particle Destiny -- %ss with Starting Energy %s GeV",getParticleNameCaps(isATau),energy); 
    TCanvas *canLastIntType = new TCanvas("canLastIntType","Final Interaction Type",600,400);

    
    gStyle->SetOptStat(0);
    TH1F *histLastInteraction;
    histLastInteraction=(TH1F*) gROOT->FindObject("histLastInteraction");
    gPad->SetTopMargin(0.12);
//	gPad->SetBottomMargin(0.2);
    histLastInteraction->SetLineWidth(3);
    histLastInteraction->SetLineColor(getNiceColour(1));
//	histLastInteraction->GetXaxis()->SetTitle("IntType (#mus)");
    histLastInteraction->GetXaxis()->SetBinLabel(2,"Energy Threshold");
    histLastInteraction->GetXaxis()->SetBinLabel(4,"Decay");
    histLastInteraction->GetXaxis()->SetBinLabel(6,"Weak (CC)");
    histLastInteraction->GetXaxis()->SetLabelSize(0.06);
    histLastInteraction->GetXaxis()->SetLabelOffset(0.02);
    histLastInteraction->SetTitle(histTitle);
    histLastInteraction->Draw();
    gPad->SetLogy();
    gPad->SetGridy();
  
    sortOutTitle(0.06);

    return canLastIntType;
}

 

char *getParticleName(int isATau)
{    
    if(isATau) return tauName;
    return muonName;
}

char *getParticleNameCaps(int isATau)
{
    if(isATau) return TauName;
    return MuonName;
}
