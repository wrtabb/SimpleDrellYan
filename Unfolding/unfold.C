#include "include/ToyModel.hh"
#include "include/Unfolding.hh"

TString inputDistribution = "../data/unfoldingHistograms";

vector<TString> histograms = {
	"hInvMassData",	
	"hInvMassReco",	
	"hInvMassBack",	
	"hInvMassDressed",	
	"hInvMassMatrix"
};
void unfold()
{
	inputDistribution += "EE";
	inputDistribution += ".root";
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	TFile*loadFile = new TFile(inputDistribution);
	int nHistograms = histograms.size();
	
	TH1F*hData   = (TH1F*)loadFile->Get(histograms.at(0));
	TH1F*hReco   = (TH1F*)loadFile->Get(histograms.at(1));
	TH1F*hBack   = (TH1F*)loadFile->Get(histograms.at(2));
	TH1F*hTrue   = (TH1F*)loadFile->Get(histograms.at(3));
	TH2F*hMatrix = (TH2F*)loadFile->Get(histograms.at(4));

	TH1F*hBack0 = (TH1F*)hBack->Clone("hBack0");
	int nBins = hBack0->GetNbinsX();
	for(int i=0;i<nBins;i++){
		// hBack0 will be blank for a closure test
		hBack0->SetBinContent(i,0);
	}

	Unfold*unfold = new Unfold();
	Unfold::RegType regType = unfold->NO_REG;

	TH1F*hUnfoldedClosure;
	TH1F*hUnfolded;
	hUnfoldedClosure = unfold->unfoldTUnfold(regType,hReco,hBack0,hTrue,hMatrix);
	//hUnfolded = unfold->unfoldTUnfold(regType,hData,hBack,hTrue,hMatrix);
	hUnfolded->SetMarkerStyle(25);
	hUnfolded->SetMarkerColor(kBlue+2);
	hUnfolded->SetLineColor(kBlue+2);

	bool logPlot = true;
	TCanvas*c1 = unfold->plotUnfolded("c1","Unfold Data",hData,hTrue,hUnfoldedClosure,logPlot);

	TH2F*hResponse = unfold->makeResponseMatrix(hMatrix);
	TCanvas*c2 = new TCanvas("c2","",0,0,1000,1000);
	c2->SetGrid();
	c2->SetLogy();
	c2->SetLogx();
	hResponse->Draw("colz");

//	TH1F*hUnfoldInv;
//	hUnfoldInv = unfold->unfoldInversion(hReco,hTrue,hMatrix);
//	TCanvas*c3 = unfold->plotUnfolded("c3","Unfold Reco Inv",hReco,hTrue,hUnfolded,logPlot);
}
