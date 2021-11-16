
TH2D*GetResponseMatrix(TH2D*hMatrix);

vector<TString> fileList_DYLL = {
        "output_data/saveFile_EE_NoPVz_DYLL_M10to50_EE.root",
        "output_data/saveFile_EE_NoPVz_DYLL_M50to100_EE.root",
        "output_data/saveFile_EE_NoPVz_DYLL_M100to200_EE.root",
        "output_data/saveFile_EE_NoPVz_DYLL_M200to400_EE.root",
        "output_data/saveFile_EE_NoPVz_DYLL_M400to500_EE.root",
        "output_data/saveFile_EE_NoPVz_DYLL_M500to700_EE.root",
        "output_data/saveFile_EE_NoPVz_DYLL_M700to800_EE.root",
        "output_data/saveFile_EE_NoPVz_DYLL_M800to1000_EE.root",
        "output_data/saveFile_EE_NoPVz_DYLL_M1000to1500_EE.root",
        "output_data/saveFile_EE_NoPVz_DYLL_M1500to2000_EE.root",
        "output_data/saveFile_EE_NoPVz_DYLL_M2000to3000_EE.root"
};

vector<TString> fileList_data = {
	"output_data/saveFile_EE_NoPVz_crab_DoubleEG_RunB.root",
	"output_data/saveFile_EE_NoPVz_crab_DoubleEG_RunC.root",
	"output_data/saveFile_EE_NoPVz_crab_DoubleEG_RunD.root",
	"output_data/saveFile_EE_NoPVz_crab_DoubleEG_RunE.root",
	"output_data/saveFile_EE_NoPVz_crab_DoubleEG_RunF.root",
	"output_data/saveFile_EE_NoPVz_crab_DoubleEG_RunG.root",
	"output_data/saveFile_EE_NoPVz_crab_DoubleEG_RunHver2.root",
	"output_data/saveFile_EE_NoPVz_crab_DoubleEG_RunHver3.root",
};

vector<TString> fileList_back = {
	"output_data/saveFile_EE_NoPVz_ST_tW.root",
	"output_data/saveFile_EE_NoPVz_ST_tbarW.root",
	"output_data/saveFile_EE_NoPVz_ttbar_M0to700.root",
	"output_data/saveFile_EE_NoPVz_ttbar_M700to1000.root",
	"output_data/saveFile_EE_NoPVz_ttbar_M1000toInf.root",
	"output_data/saveFile_EE_NoPVz_WW.root",
	"output_data/saveFile_EE_NoPVz_WZ.root",
	"output_data/saveFile_EE_NoPVz_ZZ.root",
        "output_data/saveFile_EE_NoPVz_DYLL_M10to50_TauTau.root",
        "output_data/saveFile_EE_NoPVz_DYLL_M50to100_TauTau.root",
        "output_data/saveFile_EE_NoPVz_DYLL_M100to200_TauTau.root",
        "output_data/saveFile_EE_NoPVz_DYLL_M200to400_TauTau.root",
        "output_data/saveFile_EE_NoPVz_DYLL_M400to500_TauTau.root",
        "output_data/saveFile_EE_NoPVz_DYLL_M500to700_TauTau.root",
        "output_data/saveFile_EE_NoPVz_DYLL_M700to800_TauTau.root",
        "output_data/saveFile_EE_NoPVz_DYLL_M800to1000_TauTau.root",
        "output_data/saveFile_EE_NoPVz_DYLL_M1000to1500_TauTau.root",
        "output_data/saveFile_EE_NoPVz_DYLL_M1500to2000_TauTau.root",
        "output_data/saveFile_EE_NoPVz_DYLL_M2000to3000_TauTau.root",
//	"output_data/saveFile_EE_NoPVz_WJetsToLNu_amcatnlo_ext.root",
//	"output_data/saveFile_EE_NoPVz_WJetsToLNu_amcatnlo_ext2v5.root",
};

void getUnfoldingHistograms()
{
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	//gROOT->SetBatch(true);
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


	// Plot migration matrix
	TCanvas*c2 = new TCanvas("c2","",0,0,1000,1000);
	c2->SetGrid();
	c2->SetLogy();
	c2->SetLogx();
	c2->SetLogz();
	hInvMassMatrix->GetXaxis()->SetMoreLogLabels();
	hInvMassMatrix->GetYaxis()->SetMoreLogLabels();
	hInvMassMatrix->GetXaxis()->SetNoExponent();
	hInvMassMatrix->GetYaxis()->SetNoExponent();
	hInvMassMatrix->GetXaxis()->SetTitle("m_{ee,true} [GeV]");
	hInvMassMatrix->GetYaxis()->SetTitle("m_{ee,obs} [GeV]");
	hInvMassMatrix->Draw("colz");
	c2->SaveAs("plots/invMassMatrixDressed.png");


	// Plot Matrix projections alongside 1D distributions
	// These should be identical
	double ratioRange = 0.003;
	double upperBound = 1.0-ratioRange;
	double lowerBound = 1.0+ratioRange;
	TH1F*hRatioTrue = (TH1F*)hInvMassDressed->Clone("trueRatio");
	hRatioTrue->Divide(projX);
	hRatioTrue->SetMarkerStyle(20);
	hRatioTrue->SetMarkerColor(kBlue);
	hRatioTrue->SetMinimum(1.0-ratioRange);
	hRatioTrue->SetMaximum(1.0+ratioRange);
	TH1F*hRatioReco = (TH1F*)hInvMassDYLL->Clone("recoRatio");
	hRatioReco->Divide(projY);
	hRatioReco->SetMarkerStyle(20);
	hRatioReco->SetMarkerColor(kRed);
	hRatioReco->SetMinimum(1.0-ratioRange);
	hRatioReco->SetMaximum(1.0+ratioRange);

	TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
        legend->SetTextSize(0.02);
        legend->AddEntry(hInvMassDressed,"True Distribution");
        legend->AddEntry(hInvMassDYLL,"Observed Distribution");
        legend->AddEntry(projX,"Matrix x-projection");
        legend->AddEntry(projY,"Matrix y-projection");

	TCanvas*c4 = new TCanvas("c4","",0,0,1000,1000);
	const float padmargins = 0.03;
        const float yAxisMinimum = 0.1;
        const float yAxisMaximum = 1e7;
        TPad*pad1 = new TPad("","",0,0.3,1.0,1.0);
	pad1->SetLogx();
	pad1->SetLogy();
	pad1->SetBottomMargin(padmargins);
        pad1->SetGrid();
        pad1->SetTicks(1,1);
        pad1->Draw();
        pad1->cd();
	hInvMassDressed->SetTitle("1D Distributions vs. Matrix Projections");
	hInvMassDressed->SetLabelSize(0);
	hInvMassDressed->SetTitleSize(0);
	hInvMassDressed->SetMinimum(0.1);
	hInvMassDressed->SetMarkerStyle(20);
	hInvMassDressed->SetMarkerColor(kBlue);
	hInvMassDYLL->SetMarkerStyle(20);
	hInvMassDYLL->SetMarkerColor(kRed);
	hInvMassDYLL->SetLineColor(kRed);

	hInvMassDressed->Draw("pe,same");
	hInvMassDYLL->Draw("pe,same");
	projX->Draw("hist,same");
	projY->Draw("hist,same");
	legend->Draw("same");

	double ratioSplit = 0.18;
	c4->cd();
	TPad*pad2 = new TPad("","",0,ratioSplit,1,0.3);
	pad2->SetLogx();
	pad2->SetTopMargin(padmargins);
        pad2->SetBottomMargin(0.2);
        pad2->SetGrid();
        pad2->SetTicks(1,1);
        pad2->Draw();
        pad2->cd();
	hRatioTrue->GetYaxis()->SetLabelSize(0.06);
        hRatioTrue->GetYaxis()->SetTitleSize(0.08);
        hRatioTrue->GetYaxis()->SetTitleOffset(0.3);
        hRatioTrue->GetYaxis()->SetTitle("1D/projection");
        hRatioTrue->GetXaxis()->SetLabelSize(0);
        hRatioTrue->GetXaxis()->SetTitleSize(0);
	hRatioTrue->Draw("pe");

	c4->cd();
	TPad*pad3 = new TPad("","",0,0.05,1,ratioSplit);
	pad3->SetLogx();
	pad3->SetTopMargin(padmargins);
        pad3->SetBottomMargin(0.2);
        pad3->SetGrid();
        pad3->SetTicks(1,1);
        pad3->Draw();
        pad3->cd();
	hRatioReco->GetYaxis()->SetLabelSize(0.06);
        hRatioReco->GetYaxis()->SetTitleSize(0.08);
        hRatioReco->GetYaxis()->SetTitleOffset(0.3);
        hRatioReco->GetYaxis()->SetTitle("1D/projection");
        hRatioReco->GetXaxis()->SetLabelSize(0.1);
        hRatioReco->GetXaxis()->SetTitleSize(0.1);
        hRatioReco->GetXaxis()->SetNoExponent();
        hRatioReco->GetXaxis()->SetMoreLogLabels();
        hRatioReco->GetXaxis()->SetTitle("mass [GeV]");
	hRatioReco->Draw("pe");

	c4->SaveAs("plots/invMassMatrixProjections.png");

	hInvMassData->SetName("hInvMassData");
	hInvMassDressed->SetName("hInvMassDressed");
	hInvMassDYLL->SetName("hInvMassReco");
	hInvMassBack->SetName("hInvMassBack");
	hInvMassMatrix->SetName("hInvMassMatrix");
	TH2D*hInvMassResponse = GetResponseMatrix(hInvMassMatrix);

	TCanvas*c5=new TCanvas("c5","",0,0,1000,1000);
	c5->SetGrid();
	c5->SetLogx();
	c5->SetLogy();
	hInvMassResponse->Draw("colz");
	c5->SaveAs("plots/invMassResponseMatrix.png");

	TFile*save_file = new TFile("data/unfoldingHistogramsEE.root","recreate");
	hInvMassData->Write();
	hInvMassDYLL->Write();
	hInvMassBack->Write();
	hInvMassDressed->Write();
	hInvMassMatrix->Write();
	hInvMassResponse->Write();
	save_file->Close();
}

TH2D*GetResponseMatrix(TH2D*hMatrix)
{
	TH2D*hist = (TH2D*)hMatrix->Clone("hInvMassResponse");
	double sumReco;
	int nBinsTrue = hMatrix->GetNbinsX();
	int nBinsReco = hMatrix->GetNbinsY();
	if(nBinsTrue>=nBinsReco){
		cout << "Need more reco bins than true bins for unfolding" << endl;
		cout << "Fix this and try again" << endl;
		return hist;
	}

	for(int i=1;i<=nBinsTrue;i++){
		sumReco = 0.0;
		for(int j=1;j<=nBinsReco;j++){
			sumReco += hMatrix->GetBinContent(i,j);
		}// loop over reco bins
		for(int j=1;j<=nBinsReco;j++){
			hist->SetBinContent(i,j,hMatrix->GetBinContent(i,j)/sumReco);
		}// loop over reco bins
	}// loop over true bins
	hist->RebinY(2);

	return hist;
}
