
vector<TString> file_DYLL = {
        "output_data/saveFile_EE_NoPVz_WithDressed_DYLL_M10to50_EE_v1.root",
        "output_data/saveFile_EE_NoPVz_WithDressed_DYLL_M10to50_EE_v2.root",
        "output_data/saveFile_EE_NoPVz_WithDressed_DYLL_M10to50_EE_ext1v1.root",
        "output_data/saveFile_EE_NoPVz_WithDressed_DYLL_M50to100_EE.root",
        "output_data/saveFile_EE_NoPVz_WithDressed_DYLL_M100to200_EE.root",
        "output_data/saveFile_EE_NoPVz_WithDressed_DYLL_M200to400_EE.root",
        "output_data/saveFile_EE_NoPVz_WithDressed_DYLL_M400to500_EE.root",
        "output_data/saveFile_EE_NoPVz_WithDressed_DYLL_M500to700_EE.root",
        "output_data/saveFile_EE_NoPVz_WithDressed_DYLL_M700to800_EE.root",
        "output_data/saveFile_EE_NoPVz_WithDressed_DYLL_M800to1000_EE.root",
        "output_data/saveFile_EE_NoPVz_WithDressed_DYLL_M1000to1500_EE.root",
        "output_data/saveFile_EE_NoPVz_WithDressed_DYLL_M1500to2000_EE.root",
        "output_data/saveFile_EE_NoPVz_WithDressed_DYLL_M2000to3000_EE.root"
};

vector<TString> file_dataEE = {
	"output_data/testLoadFromKorea_crab_DoubleEG_RunB.root",
	"output_data/testLoadFromKorea_crab_DoubleEG_RunC.root",
	"output_data/testLoadFromKorea_crab_DoubleEG_RunD.root",
	"output_data/testLoadFromKorea_crab_DoubleEG_RunE.root",
	"output_data/testLoadFromKorea_crab_DoubleEG_RunF.root",
	"output_data/testLoadFromKorea_crab_DoubleEG_RunHver2.root",
	"output_data/testLoadFromKorea_crab_DoubleEG_RunHver3.root"
};

void getUnfoldingHistograms()
{
	vector<TFile*> file;
	vector<TH1D*> vInvMassHard;
	vector<TH1D*> vRapidityHard;
	vector<TH1D*> vInvMassDressed;
	vector<TH1D*> vRapidityDressed;
	vector<TH2D*> vInvMassMatrixHard;
	vector<TH2D*> vInvMassMatrixDressed;
	vector<TH2D*> vRapidityMatrixHard;
	vector<TH2D*> vRapidityMatrixDressed;

	TH1D*hInvMassHard; 
	TH1D*hRapidityHard;
	TH1D*hInvMassDressed;
	TH1D*hRapidityDressed;
	TH2D*hInvMassMatrixHard;
	TH2D*hInvMassMatrixDressed;
	TH2D*hRapidityMatrixHard;
	TH2D*hRapidityMatrixDressed;

	vector<TFile*> file_loadData;
	vector<TH1D*> vInvMassReco;
	vector<TH1D*> vRapidityReco;

	TH1D*hInvMassReco;
	TH1D*hRapidityReco;

	int nDataFiles = file_dataEE.size();
	for(int i=0;i<nDataFiles;i++){
		file_loadData.push_back(new TFile(file_dataEE.at(i)));
		if(i==0){
			hInvMassReco = (TH1D*)file_loadData.at(i)->Get("histInvMassReco");
			hRapidityReco = (TH1D*)file_loadData.at(i)->Get("histRapidityReco");
		}
		vInvMassReco.push_back((TH1D*)file_loadData.at(i)->Get("histInvMassReco"));
		vRapidityReco.push_back((TH1D*)file_loadData.at(i)->Get("histRapidityReco"));
		if(i>0){
			hInvMassReco->Add(vInvMassReco.at(i));
			hRapidityReco->Add(vRapidityReco.at(i));
		}
	}// end loop over data files

	int nFiles = file_DYLL.size();
	for(int i=0;i<nFiles;i++){
		file.push_back(new TFile(file_DYLL.at(i)));
		if(i==0){
			hInvMassHard = (TH1D*)file.at(i)->Get("histInvMassHard");
			hRapidityHard = (TH1D*)file.at(i)->Get("histRapidityHard");
			hInvMassDressed = (TH1D*)file.at(i)->Get("histInvMassDressed");
			hRapidityDressed = (TH1D*)file.at(i)->Get("histRapidityDressed");
			hInvMassMatrixHard = (TH2D*)file.at(i)->Get("histMatrixInvMassHard");
			hInvMassMatrixDressed = (TH2D*)file.at(i)->Get("histMatrixInvMassDressed");
			hRapidityMatrixHard = (TH2D*)file.at(i)->Get("histMatrixRapidityHard");
			hRapidityMatrixDressed = (TH2D*)file.at(i)->Get("histMatrixRapidityDressed");
		}
		vInvMassHard.push_back((TH1D*)file.at(i)->Get("histInvMassHard"));
		vRapidityHard.push_back((TH1D*)file.at(i)->Get("histRapidityHard"));
		vInvMassDressed.push_back((TH1D*)file.at(i)->Get("histInvMassDressed"));
		vRapidityDressed.push_back((TH1D*)file.at(i)->Get("histRapidityDressed"));
		vInvMassMatrixHard.push_back((TH2D*)file.at(i)->Get("histMatrixInvMassHard"));
		vInvMassMatrixDressed.push_back((TH2D*)file.at(i)->Get("histMatrixInvMassDressed"));
		vRapidityMatrixHard.push_back((TH2D*)file.at(i)->Get("histMatrixRapidityHard"));
		vRapidityMatrixDressed.push_back((TH2D*)file.at(i)->Get("histMatrixRapidityDressed"));

		if(i>0){
			hInvMassHard->Add(vInvMassHard.at(i));
			hRapidityHard->Add(vRapidityHard.at(i));
			hInvMassDressed->Add(vInvMassDressed.at(i));
			hRapidityDressed->Add(vRapidityDressed.at(i));
			hInvMassMatrixHard->Add(vInvMassMatrixHard.at(i));
			hInvMassMatrixDressed->Add(vInvMassMatrixDressed.at(i));
			hRapidityMatrixHard->Add(vRapidityMatrixHard.at(i));
			hRapidityMatrixDressed->Add(vRapidityMatrixDressed.at(i));
		}
	}// end loop over files

	TCanvas*c1=new TCanvas("c1","",0,0,1000,1000);
	c1->SetLogx();
	c1->SetLogy();
	c1->SetLogz();
	c1->SetGrid();
	hInvMassMatrixHard->Draw("colz");
	c1->SaveAs("plots/migrationMatrixHard.png");

	TCanvas*c2=new TCanvas("c2","",0,0,1000,1000);
	c2->SetLogx();
	c2->SetLogy();
	c2->SetLogz();
	c2->SetGrid();
	hInvMassMatrixDressed->Draw("colz");
	c2->SaveAs("plots/migrationMatrixHard.png");

	TCanvas*c3=new TCanvas("c3","",0,0,1000,1000);
	c3->SetLogx();
	c3->SetLogy();
	c3->SetGrid();
	hInvMassDressed->SetLineColor(kRed);
	hInvMassHard->SetLineColor(kBlue);
	hInvMassReco->SetMarkerStyle(20);
	hInvMassReco->SetMarkerColor(kBlack);
	hInvMassHard->Draw("hist");
	hInvMassDressed->Draw("hist,same");
	hInvMassReco->Draw("pe,same");
	c3->SaveAs("plots/testInvMassHardVsDressed.png");

	TFile*save_file = new TFile("data/unfoldingHistograms.root","recreate");
	hInvMassMatrixDressed->Write();
	hInvMassDressed->Write();
	hInvMassReco->Write();
	save_file->Close();
	
}
