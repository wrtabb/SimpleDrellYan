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

	Unfold*unfold = new Unfold();
	Unfold::RegType regType = unfold->NO_REG;

	hReco->Add(hBack);
	TH1F*hUnfoldClosure = unfold->unfoldTUnfold(regType,hReco,hBack,hTrue,hMatrix);
	TH1F*hUnfold = unfold->unfoldTUnfold(regType,hData,hBack,hTrue,hMatrix);

	bool logPlot = true;
	TCanvas*c1 = 
		unfold->plotUnfolded("c1","Closure Test",hReco,hTrue,hUnfoldClosure,logPlot);
	TCanvas*c2 = 
		unfold->plotUnfolded("c2","Data Unfold",hData,hTrue,hUnfold,logPlot);
	c1->SaveAs("../plots/unfoldedEE_ClosureTest.png");
	c2->SaveAs("../plots/unfoldedEE_Data.png");
/*
	TH2F*hResponse = unfold->makeResponseMatrix(hMatrix);
	TCanvas*c2 = new TCanvas("c2","",0,0,1000,1000);
	c2->SetGrid();
	c2->SetLogy();
	c2->SetLogx();
	hResponse->Draw("colz");
*/
//	TH1F*hUnfoldInv;
//	hUnfoldInv = unfold->unfoldInversion(hReco,hTrue,hMatrix);
//	TCanvas*c3 = unfold->plotUnfolded("c3","Unfold Reco Inv",hReco,hTrue,hUnfolded,logPlot);
}
