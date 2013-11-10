

char muonName[]="muon";
char MuonName[]="Muon";
char tauName[]="tau";
char TauName[]="Tau";


void plotSummaryGraphs()  {

    
    gStyle->SetStatH(0.3);
    gStyle->SetStatY(0.88);
    gStyle->SetStatW(0.3);
    gStyle->SetTitleOffset(0.7,"y");
    TCanvas *canBangs = new TCanvas("canBangs","Final Energy",800,800);
    TPaveLabel *pl = new TPaveLabel(0.1,0.96,0.9,0.99,
				    "Number of 'Bangs'","br NDC");
    pl->SetBorderSize(0);
    pl->SetFillColor(0);
    pl->SetFillStyle(0);
    pl->Draw();
    TPad *subCanBangs = new TPad("subCanBangs","",0,0,1,0.95);
    subCanBangs->Draw();
    subCanBangs->cd();
    subCanBangs->Divide(2,4);
    subCanBangs->Update();
    
//    gStyle->SetOptStat(0);
    char fileName[80];   
    char histTitle[80];

    char theEnergies[4][5]={"1e9","1e10","1e11","1e12"};
    Double_t theFourEnergies[4]={1e9,1e10,1e11,1e12};
    char energyCut[8][5]={"1e8","3e8","1e9","3e9","1e10","3e10","1e11","3e11"};


//    int theColours[2][3]={{50,42,46},{40,30,38}};   
    
    Double_t numBangsMean[2][8][4]={0};
    Double_t numBangsRMS[2][8][4]={0};
    Double_t energyErrs[4]={0};

    char histName[180];
     
    for(int isATau=0;isATau<=1;isATau++) {
	for(int i=0;i<4;i++) {		 
	    sprintf(fileName,
		    "newestHist%sFile%sIce.root",getParticleName(isATau),theEnergies[i]);
	    
	    TFile fp(fileName,"OLD");
	    for(int j=0;j<8;j++) {
		sprintf(histName,"histNum%d",j);
		TH1F *histTemp = (TH1F*) fp.Get(histName);
		if(histTemp) {
		    numBangsMean[isATau][j][i]=histTemp->GetMean();
		    numBangsRMS[isATau][j][i]=histTemp->GetRMS();
		}
	    }
	}
    }

    Double_t myColours[2]={50,30};
    char graphName[180];
    TGraphErrors *grs[2];
    for(int i=0;i<8;i++) {
	subCanBangs->cd(i+1);
	gPad->SetLogx();
// 	gPad->SetLogy();
	for(int isATau=1;isATau>=0;isATau--) {	    
	    grs[isATau] 
		= new TGraphErrors(4,theFourEnergies,numBangsMean[isATau][i],
				   energyErrs,numBangsRMS[isATau][i]);
	    
	    grs[isATau]->SetLineColor(myColours[isATau]);
	    grs[isATau]->SetLineWidth(3);
	    if(isATau) {
		sprintf(graphName,"E_{cut} > %s GeV",energyCut[i]);
		grs[isATau]->SetTitle(graphName);
		grs[isATau]->Draw("al");
		grs[isATau]->GetXaxis()->SetTitle("Energy (GeV)");
		grs[isATau]->GetYaxis()->SetTitle("Num Bangs");
	    }
	    else
		grs[isATau]->Draw("l");
	}
	if(i==0) {
	    TLegend *leg = new TLegend(0.15,0.65,0.3,0.85);
	    leg->SetFillColor(0);
	    leg->SetFillStyle(0);
	    leg->SetBorderSize(0);
	    leg->AddEntry(grs[0],"Muons","l");
	    leg->AddEntry(grs[1],"Taus","l");
	    leg->Draw();
	}
 
    }
    for(int i=0;i<8;i++) {
	subCanBangs->cd(i+1);
	sortOutTitle(0.07);
    }

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
