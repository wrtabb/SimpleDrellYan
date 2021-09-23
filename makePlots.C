vector<TString> files_data = {
	"crab_DoubleEG_RunB",
	"crab_DoubleEG_RunC",
	"crab_DoubleEG_RunD",
	"crab_DoubleEG_RunE",
	"crab_DoubleEG_RunF",
	"crab_DoubleEG_RunG",
	"crab_DoubleEG_RunHver2",
	"crab_DoubleEG_RunHver3"
};

vector<TString> files_LL = {
        "DYLL_M10to50_EE",              
        "DYLL_M50to100_EE",             
        "DYLL_M100to200_EE",            
        "DYLL_M200to400_EE",            
        "DYLL_M400to500_EE",            
        "DYLL_M500to700_EE",            
        "DYLL_M700to800_EE",            
        "DYLL_M800to1000_EE",           
        "DYLL_M1000to1500_EE",          
        "DYLL_M1500to2000_EE",          
        "DYLL_M2000to3000_EE"           
};

TH1D*GetHistogram(vector<TString> vecOfFiles);
void makePlots()
{
	gStyle->SetOptStat(0);
	gROOT->SetBatch(true);

	TH1D*hData = GetHistogram(files_data);
	hData->SetMarkerStyle(20);
	hData->SetMarkerColor(kBlack);
	hData->SetMinimum(0.01);
	TH1D*hLL = GetHistogram(files_LL);
	hLL->SetFillColor(kOrange-2);
	hLL->SetMinimum(0.01);

	TCanvas*c1 = new TCanvas("c1","",0,0,1000,1000);
	c1->SetGrid();
	c1->SetLogx();
	c1->SetLogy();
	hLL->Draw("hist");
	hData->Draw("pe,same");
	c1->SaveAs("dataVsMC_signalOnly_BadScaling.png");
}

TH1D*GetHistogram(vector<TString> vecOfFiles)
{
	const int nFiles = vecOfFiles.size();
	TFile*fileLoad[nFiles];
	TH1D*hist[nFiles];
	TH1D*histCombo;

	for(int i=0;i<nFiles;i++){
		TString histName = "hist";
		histName += vecOfFiles.at(i); 
		TString loadName = "output_data/saveFile";
		loadName += vecOfFiles.at(i);
		loadName += ".root";
		fileLoad[i] = new TFile(loadName);

		cout << "Getting histogram " << histName << " from files " << loadName << endl;
		hist[i] = (TH1D*)fileLoad[i]->Get(histName);
		if(i==0) histCombo = (TH1D*)hist[i]->Clone();
		else histCombo->Add(hist[i]);
	}

	return histCombo;
}
