#include "include/ToyModel.hh"
#include "include/Unfolding.hh"

TString inputDistribution = "data/unfoldingHistograms.root";
vector<TString> histograms = {
	"hInvMassData",	
	"hInvMassReco",	
	"hInvMassBack",	
	"hInvMassDressed",	
	"hInvMassMatrix"
};
void unfold()
{
	TFile*loadFile = new TFile(inputDistribution);
	int nHistograms = histograms.size();
	
	TH1F*hData   = (TH1F*)loadFile->Get(histograms.at(0));
	TH1F*hReco   = (TH1F*)loadFile->Get(histograms.at(1));
	TH1F*hBack   = (TH1F*)loadFile->Get(histograms.at(2));
	TH1F*hTrue   = (TH1F*)loadFile->Get(histograms.at(3));
	TH2F*hMatrix = (TH2F*)loadFile->Get(histograms.at(4));

	Unfold*unfold = new Unfold();
	Unfold::RegType regType = unfold->NO_REG;

	TH1F*hUnfoldedClosure;
	TH1F*hUnfolded;
	hUnfoldedClosure = unfold->unfoldTUnfold(regType,hReco,hTrue,hMatrix);
	hUnfolded = unfold->unfoldTUnfold(regType,hReco,hTrue,hMatrix);
	hUnfolded->SetMarkerStyle(25);
	hUnfolded->SetMarkerColor(kBlue+2);
	hUnfolded->SetLineColor(kBlue+2);

	bool logPlot = true;
	TCanvas*c1 = unfold->plotUnfolded("c1","Closure Test",hReco,hTrue,hUnfoldedClosure,logPlot);
}
