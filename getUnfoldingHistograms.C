
vector<TString> fileList_DYLL = {
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_DYLL_M10to50_MuMu.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_DYLL_M50to100_MuMu.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_DYLL_M100to200_MuMu.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_DYLL_M200to400_MuMu.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_DYLL_M400to500_MuMu.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_DYLL_M500to700_MuMu.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_DYLL_M700to800_MuMu.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_DYLL_M800to1000_MuMu.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_DYLL_M1000to1500_MuMu.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_DYLL_M1500to2000_MuMu.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_DYLL_M2000to3000_MuMu.root",
};

vector<TString> fileList_data = {
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_SingleMuon_Run2016B.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_SingleMuon_Run2016C.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_SingleMuon_Run2016D.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_SingleMuon_Run2016E.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_SingleMuon_Run2016F.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_SingleMuon_Run2016G.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_SingleMuon_Run2016Hver2.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_SingleMuon_Run2016Hver3.root",
};

vector<TString> fileList_back = {
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_ST_tW.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_ST_tbarW.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_ttbar_M0to700.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_ttbar_M700to1000.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_ttbar_M1000toInf.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_WW.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_WZ.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_ZZ.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_DYLL_M10to50_TauTau.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_DYLL_M50to100_TauTau.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_DYLL_M100to200_TauTau.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_DYLL_M200to400_TauTau.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_DYLL_M400to500_TauTau.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_DYLL_M500to700_TauTau.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_DYLL_M700to800_TauTau.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_DYLL_M800to1000_TauTau.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_DYLL_M1000to1500_TauTau.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_DYLL_M1500to2000_TauTau.root",
	"output_data/saveFile_MuMu_NoSF_NoPVz_WithDressed_DYLL_M2000to3000_TauTau.root",
};

void getUnfoldingHistograms()
{
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	gROOT->SetBatch(true);
	TFile*file_data;
	TFile*file_DYLL;
	TFile*file_back;

	TH1D*hInvMassData;
	TH1D*hInvMassDYLL;
	TH1D*hInvMassDressed;
	TH1D*hInvMassBack;
	TH2D*hInvMassMatrix;

	int nDataFiles = fileList_data.size();
	int nDYLLFiles = fileList_DYLL.size();
	int nBackFiles = fileList_back.size();

	// Get data files and place into histograms
	for(int i=0;i<nDataFiles;i++){
		file_data = new TFile(fileList_data.at(i));
		if(i==0)hInvMassData = (TH1D*)file_data->Get("histInvMassReco");
		else hInvMassData->Add((TH1D*)file_data->Get("histInvMassReco"));
	}// end loop over data files

	// Get background files and place into histograms
	for(int i=0;i<nBackFiles;i++){
		file_back = new TFile(fileList_back.at(i));
		if(i==0)hInvMassBack = (TH1D*)file_back->Get("histInvMassReco");
		else hInvMassBack->Add((TH1D*)file_back->Get("histInvMassReco"));
	}// end loop over background files

	// Get DYLL files and place into histograms
	for(int i=0;i<nDYLLFiles;i++){
		file_DYLL = new TFile(fileList_DYLL.at(i));
		if(i==0){
			hInvMassDressed = (TH1D*)file_DYLL->Get("histInvMassDressed");
			hInvMassDYLL = (TH1D*)file_DYLL->Get("histInvMassReco");
			hInvMassMatrix  = (TH2D*)file_DYLL->Get("histMatrixInvMassDressed");
		}
		else{ 
			hInvMassDressed->Add((TH1D*)file_DYLL->Get("histInvMassDressed"));
			hInvMassDYLL   ->Add((TH1D*)file_DYLL->Get("histInvMassReco"));
			hInvMassMatrix ->Add((TH2D*)file_DYLL->Get("histMatrixInvMassDressed"));
		}
	}// end loop over DYLL files
	
	TH1F*projX = (TH1F*)hInvMassMatrix->ProjectionX();
	projX->SetLineColor(kBlue);
	TH1F*projY = (TH1F*)hInvMassMatrix->ProjectionY();
	projY->SetLineColor(kRed);

	hInvMassData->SetMarkerStyle(20);
	hInvMassDressed->SetMinimum(0.1);
	hInvMassDressed->SetLineColor(kBlue);
	hInvMassBack->SetMarkerColor(kBlue);
	hInvMassBack->SetMarkerStyle(20);

	TCanvas*c1 = new TCanvas("c1","",0,0,1000,1000);
	c1->SetGrid();
	c1->SetLogy();
	c1->SetLogx();
	hInvMassDressed->Draw("hist,same");
	hInvMassData->Draw("pe,same");
	hInvMassBack->Draw("pe,same");

	TCanvas*c2 = new TCanvas("c2","",0,0,1000,1000);
	c2->SetGrid();
	c2->SetLogy();
	c2->SetLogx();
	c2->SetLogz();
	c2->SetRightMargin(0.13);
	c2->SetLeftMargin(0.13);
	hInvMassMatrix->GetXaxis()->SetMoreLogLabels();
	hInvMassMatrix->GetYaxis()->SetMoreLogLabels();
	hInvMassMatrix->GetXaxis()->SetNoExponent();
	hInvMassMatrix->GetYaxis()->SetNoExponent();
	hInvMassMatrix->GetXaxis()->SetTitle("m_{ee,true} [GeV]");
	hInvMassMatrix->GetYaxis()->SetTitle("m_{ee,obs} [GeV]");
	hInvMassMatrix->Draw("colz");
	c2->SaveAs("plots/invMassMatrixDressed.png");

	hInvMassDYLL->SetFillColor(kOrange-2);
	hInvMassDYLL->SetLineColor(kOrange+3);
	hInvMassBack->SetFillColor(kViolet+2);
	THStack*hStack = new THStack("hStack","");
	hStack->Add(hInvMassBack);
	hStack->Add(hInvMassDYLL);
	TCanvas*c3 = new TCanvas("c3","",0,0,1000,1000);
	c3->SetGrid();
	c3->SetLogy();
	c3->SetLogx();
	hStack->Draw("hist");		
	hInvMassData->Draw("pe,same");
	c3->SaveAs("plots/invMassBackDYLLData.png");


	TCanvas*c4 = new TCanvas("c4","",0,0,1000,1000);
	c4->SetGrid();
	c4->SetLogy();
	c4->SetLogx();
	c4->SetRightMargin(0.13);
	c4->SetLeftMargin(0.13);
	hInvMassDressed->SetMarkerStyle(20);
	hInvMassDressed->SetMarkerColor(kBlue);
	hInvMassDYLL->SetMarkerStyle(20);
	hInvMassDYLL->SetMarkerColor(kRed);
	hInvMassDressed->Draw("pe,same");
	hInvMassDYLL->Draw("pe,same");
	projX->Draw("hist,same");
	projY->Draw("hist,same");
	c4->SaveAs("plots/invMassMatrixProjections.png");

	hInvMassData->SetName("hInvMassData");
	hInvMassDressed->SetName("hInvMassDressed");
	hInvMassDYLL->SetName("hInvMassReco");
	hInvMassBack->SetName("hInvMassBack");
	hInvMassMatrix->SetName("hInvMassMatrix");

	TFile*save_file = new TFile("data/unfoldingHistogramsMuMu.root","recreate");
	hInvMassData->Write();
	hInvMassDYLL->Write();
	hInvMassBack->Write();
	hInvMassDressed->Write();
	hInvMassMatrix->Write();

	save_file->Close();
}
