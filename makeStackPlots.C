
vector<TString> file_data = {
	"output_data/crab_DoubleEG_RunB.root",           // 0
	"output_data/crab_DoubleEG_RunC.root",           // 1
	"output_data/crab_DoubleEG_RunD.root",           // 2
	"output_data/crab_DoubleEG_RunE.root",           // 3
	"output_data/crab_DoubleEG_RunF.root",           // 4
	"output_data/crab_DoubleEG_RunG.root",           // 5
	"output_data/crab_DoubleEG_RunHver2.root",       // 6
	"output_data/crab_DoubleEG_RunHver3.root",       // 7
};

vector<TString> file_diboson = {
	"output_data/WW.root",
	"output_data/WZ.root",
	"output_data/ZZ.root"
};

vector<TString> file_taus = {
	"output_data/DYLL_M1000to1500_TauTau.root",
	"output_data/DYLL_M100to200_TauTau.root",
	"output_data/DYLL_M10to50_TauTau.root",
	"output_data/DYLL_M1500to2000_TauTau.root",
	"output_data/DYLL_M2000to3000_TauTau.root",
	"output_data/DYLL_M200to400_TauTau.root",
	"output_data/DYLL_M400to500_TauTau.root",
	"output_data/DYLL_M500to700_TauTau.root",
	"output_data/DYLL_M50to100_TauTau.root",
	"output_data/DYLL_M700to800_TauTau.root",
	"output_data/DYLL_M800to1000_TauTau.root"
};

vector<TString> file_DYLL = {
	"output_data/DYLL_M1000to1500_EE.root",
	"output_data/DYLL_M100to200_EE.root",
	"output_data/DYLL_M10to50_EE.root",
	"output_data/DYLL_M1500to2000_EE.root",
	"output_data/DYLL_M2000to3000_EE.root",
	"output_data/DYLL_M200to400_EE.root",
	"output_data/DYLL_M400to500_EE.root",
	"output_data/DYLL_M500to700_EE.root",
	"output_data/DYLL_M50to100_EE.root",
	"output_data/DYLL_M700to800_EE.root",
	"output_data/DYLL_M800to1000_EE.root"
};

vector<TString> file_top = {
	"output_data/ST_tW.root",
	"output_data/ST_tbarW.root",
	"output_data/ttbar_M0to700.root",
	"output_data/ttbar_M1000toInf.root",
	"output_data/ttbar_M700to1000.root"
};

vector<TString> file_Fake = {
	"output_data/WJetsToLNu_amcatnlo.root"
};

TH1D*GetHistogram(vector<TString> filesvector);

void makeStackPlots()
{
	gROOT->SetBatch(true);
	TCanvas*c1 = new TCanvas("c1","",0,0,1000,1000);
	c1->SetGrid();
	c1->SetLogx();
	c1->SetLogy();

	TH1D*hMassDYLL = GetHistogram(file_DYLL);
	TH1D*hMassTops = GetHistogram(file_top);
	TH1D*hMassFake = GetHistogram(file_Fake);
	TH1D*hMassDiboson = GetHistogram(file_diboson);
	//TH1D*hMassTau = GetHistogram(file_taus);

	THStack*hStack_Mass = new THStack("hStack_Mass","");
	hStack_Mass->Add(hMassTops);
	hStack_Mass->Add(hMassFake);
	hStack_Mass->Add(hMassDiboson);
	//hStack_Mass->Add(hMassTau);
	hStack_Mass->Add(hMassDYLL);
	
	TH1D*hMassData = GetHistogram(file_data);
	hStack_Mass->SetMinimum(1);;
	hStack_Mass->Draw("hist");
	hMassData->Draw("pe,same");

	c1->SaveAs("plots/tempDataVsMCStack.png");
}

TH1D*GetHistogram(vector<TString> filesvector)
{
	TString histMassName = "hMassReco";
	int nFiles = filesvector.size();
	TFile*file_load[nFiles];
	TH1D*hMassArray[nFiles];
	TH1D*hMass;
	for(int i=0;i<nFiles;i++){
		file_load[i] = new TFile(filesvector.at(i));
		hMassArray[i] = (TH1D*)file_load[i]->Get(histMassName);
		if(i==0) hMass = (TH1D*)hMassArray[i]->Clone();
		else hMass->Add(hMassArray[i]);
	}

	return hMass;
}

