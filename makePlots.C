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
//        "DYLL_M10to50_EE",              
//        "DYLL_M50to100_EE",             
//        "DYLL_M100to200_EE",            
        "DYLL_M200to400_EE",            
        "DYLL_M400to500_EE",            
        "DYLL_M500to700_EE",            
        "DYLL_M700to800_EE",            
        "DYLL_M800to1000_EE",           
        "DYLL_M1000to1500_EE",          
        "DYLL_M1500to2000_EE",          
        "DYLL_M2000to3000_EE"           
};

TH1D*GetHistogram(vector<TString> vecOfFiles,TString histTag);
void makePlots()
{
	gStyle->SetOptStat(0);
//	gROOT->SetBatch(true);

	TH1D*hDataInvMass = GetHistogram(files_data,"InvMass");
	hDataInvMass->SetMarkerStyle(20);
	hDataInvMass->SetMarkerColor(kBlack);
	hDataInvMass->SetMinimum(0.01);
	TH1D*hLLInvMass = GetHistogram(files_LL,"InvMass");
	hLLInvMass->SetFillColor(kOrange-2);
	hLLInvMass->SetMinimum(0.01);

	TH1D*hDataRapidity = GetHistogram(files_data,"Rapidity");
	hDataRapidity->SetMarkerStyle(20);
	hDataRapidity->SetMarkerColor(kBlack);
	hDataRapidity->SetMinimum(0.01);
	TH1D*hLLRapidity = GetHistogram(files_LL,"Rapidity");
	hLLRapidity->SetFillColor(kOrange-2);
	hLLRapidity->SetMinimum(0.01);

	TH1D*hDataPtLead = GetHistogram(files_data,"PtLead");
	hDataPtLead->SetMarkerStyle(20);
	hDataPtLead->SetMarkerColor(kBlack);
	hDataPtLead->SetMinimum(0.01);
	TH1D*hLLPtLead = GetHistogram(files_LL,"PtLead");
	hLLPtLead->SetFillColor(kOrange-2);
	hLLPtLead->SetMinimum(0.01);

	TH1D*hDataPtSub = GetHistogram(files_data,"PtSub");
	hDataPtSub->SetMarkerStyle(20);
	hDataPtSub->SetMarkerColor(kBlack);
	hDataPtSub->SetMinimum(0.01);
	TH1D*hLLPtSub = GetHistogram(files_LL,"PtSub");
	hLLPtSub->SetFillColor(kOrange-2);
	hLLPtSub->SetMinimum(0.01);

	TCanvas*c1 = new TCanvas("c1","",0,0,1000,1000);
	c1->SetGrid();
	c1->SetLogx();
	c1->SetLogy();
	hLLInvMass->Draw("hist");
	hDataInvMass->Draw("pe,same");
	c1->SaveAs("plots/dataVsMC_signalOnly_BadScaling_invmass.png");

	TCanvas*c2 = new TCanvas("c2","",0,0,1000,1000);
	c2->SetGrid();
	c2->SetLogy();
	hLLRapidity->Draw("hist");
	hDataRapidity->Draw("pe,same");
	c2->SaveAs("plots/dataVsMC_signalOnly_BadScaling_rapidity.png");

	TCanvas*c3 = new TCanvas("c3","",0,0,1000,1000);
	c3->SetGrid();
	c3->SetLogx();
	c3->SetLogy();
	hLLPtLead->Draw("hist");
	hDataPtLead->Draw("pe,same");
	c3->SaveAs("plots/dataVsMC_signalOnly_BadScaling_ptlead.png");

	TCanvas*c4 = new TCanvas("c4","",0,0,1000,1000);
	c4->SetGrid();
	c4->SetLogx();
	c4->SetLogy();
	hLLPtSub->Draw("hist");
	hDataPtSub->Draw("pe,same");
	c4->SaveAs("plots/dataVsMC_signalOnly_BadScaling_ptsub.png");
}

TH1D*GetHistogram(vector<TString> vecOfFiles,TString histTag)
{
	const int nFiles = vecOfFiles.size();
	TFile*fileLoad[nFiles];
	TH1D*hist[nFiles];
	TH1D*histCombo;

	for(int i=0;i<nFiles;i++){
		TString histName = "hist";
		histName += vecOfFiles.at(i); 
		histName += histTag;
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
