

char muonName[]="muon";
char MuonName[]="Muon";
char tauName[]="tau";
char TauName[]="Tau";


void plotFinalEnergy()  {

    
    gStyle->SetStatH(0.3);
    gStyle->SetStatY(0.88);
    gStyle->SetStatW(0.3);
    char canTitle[180];
    sprintf(canTitle,"Energy at Terminal Interaction (10000 tries)"); 
    TCanvas *canFinalEnergy = new TCanvas("canFinalEnergy","Final Energy",800,800);
    TPaveLabel *pl = new TPaveLabel(0.1,0.96,0.9,0.99,canTitle,"br NDC");
    pl->SetBorderSize(0);
    pl->SetFillColor(0);
    pl->SetFillStyle(0);
    pl->Draw();
    TPad *subCanFinalEnergy = new TPad("subCanFinalEnergy","",0,0,1,0.95);
    subCanFinalEnergy->Draw();
    subCanFinalEnergy->cd();
    subCanFinalEnergy->Divide(2,4);
    subCanFinalEnergy->Update();
    


//    gStyle->SetOptStat(0);
    char fileName[80];   
    char histTitle[80];

    char theEnergies[4][5]={"1e9","1e10","1e11","1e12"};

    int theColours[2][3]={{50,42,46},{40,30,38}};
    
    for(int isATau=0;isATau<=1;isATau++) {
	for(int i=0;i<4;i++) {		 
	    sprintf(fileName,
		    "newest%sFile%sIce.root",getParticleNameCaps(isATau),theEnergies[i]);
	    sprintf(histTitle,"%s -- %s GeV",getParticleNameCaps(isATau),theEnergies[i]);
	    TFile *fp = new TFile(fileName);
	    TH1F *histEnergy = 
		new TH1F("histEnergy",histTitle,100,7.5,12.5);
	    TH1F *histEnergy2 = 
		new TH1F("histEnergy2","Last Energy",100,7.5,12.5);
	    TH1F *histEnergy3 = 
		new TH1F("histEnergy3","Last Energy",100,7.5,12.5);

	    TTree *theTree = (TTree*) fp->Get("theTree");
	    cout << theTree->GetEntries() << endl;
	    subCanFinalEnergy->cd((2*i)+isATau+1);
	    gPad->SetTopMargin(0.12);
//	gPad->SetBottomMargin(0.2);
	    histEnergy->SetLineWidth(3);
	    histEnergy->SetLineColor(theColours[isATau][0]);
	    histEnergy2->SetLineWidth(3);
	    histEnergy2->SetLineColor(theColours[isATau][1]);
	    histEnergy3->SetLineWidth(3);
	    histEnergy3->SetLineColor(theColours[isATau][2]);
	    theTree->Draw("log10(stepIntEnergy)>>histEnergy","stepIntType>=4");
//	histEnergy->GetXaxis()->SetTitle("IntType (#mus)");
// 	histEnergy->GetXaxis()->SetBinLabel(2,"Bremsstrahlung");
// 	histEnergy->GetXaxis()->SetBinLabel(3,"Pair Production");
// 	histEnergy->GetXaxis()->SetBinLabel(4,"Photonuclear");
// 	histEnergy->GetXaxis()->SetLabelSize(0.09);
// 	histEnergy->GetXaxis()->SetLabelOffset(0.02);
 	    theTree->Draw("log10(stepIntEnergy)>>histEnergy2","stepIntType==4");
 	    theTree->Draw("log10(stepIntEnergy)>>histEnergy3","stepIntType==5");

	    histEnergy->DrawCopy();
	    if(histEnergy2->GetEntries())
		histEnergy2->DrawCopy("same");
	    if(histEnergy3->GetEntries())
		histEnergy3->DrawCopy("same");
	    
	    if(histEnergy->GetEntries())
		gPad->SetLogy();

	    if(i==0) {
		TLegend *leg = new TLegend(0.7,0.2,0.9,0.6);
		leg->SetFillColor(0);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->AddEntry(histEnergy2,"Decay","l");
		leg->AddEntry(histEnergy3,"Weak (CC)","l");
		leg->AddEntry(histEnergy,"Either","l");
		leg->Draw("same");
	    }
	}
    } 
    for(int i=0;i<8;i++) {
//	cout << "Doing title "<< i+1 << endl;
	subCanFinalEnergy->cd(i+1);
	sortOutTitle(0.07);
    }
//    gStyle->SetOptStat(1110);

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
