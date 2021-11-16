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
	hTrue->SetMarkerColor(kRed+2);
	TCanvas*c1 = 
		unfold->plotUnfolded("c1","Closure Test",hReco,hTrue,hUnfoldClosure,logPlot);
	TCanvas*c2 = 
		unfold->plotUnfolded("c2","Data Unfold",hData,hTrue,hUnfold,logPlot);

	TH2F*hResponse = unfold->makeResponseMatrix(hMatrix);
	double condition = unfold->GetConditionNumber(hResponse);
cout << "Condition number of response matrix: " << condition << endl;
	TCanvas*c3 = new TCanvas("c3","",0,0,1000,1000);
	c3->SetGrid();
	c3->SetLogy();
	c3->SetLogx();
	hResponse->Draw("colz");

	c1->SaveAs("../plots/unfoldedEE_ClosureTest.png");
	c2->SaveAs("../plots/unfoldedEE_Data_NoFakes.png");
	c3->SaveAs("../plots/unfoldedEE_ResponseMatrix.png");
}
