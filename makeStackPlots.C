
vector<TString> file_data = {
	"output_data/saveFile_EE_NoPVz_crab_DoubleEG_RunB.root",           // 0
	"output_data/saveFile_EE_NoPVz_crab_DoubleEG_RunC.root",           // 1
	"output_data/saveFile_EE_NoPVz_crab_DoubleEG_RunD.root",           // 2
	"output_data/saveFile_EE_NoPVz_crab_DoubleEG_RunE.root",           // 3
	"output_data/saveFile_EE_NoPVz_crab_DoubleEG_RunF.root",           // 4
	"output_data/saveFile_EE_NoPVz_crab_DoubleEG_RunG.root",           // 5
	"output_data/saveFile_EE_NoPVz_crab_DoubleEG_RunHver2.root",       // 6
	"output_data/saveFile_EE_NoPVz_crab_DoubleEG_RunHver3.root"       // 7
};

vector<TString> file_ew= {
	"output_data/saveFile_EE_NoPVz_WW.root",
	"output_data/saveFile_EE_NoPVz_WZ.root",
	"output_data/saveFile_EE_NoPVz_ZZ.root",
/*
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
	"output_data/saveFile_EE_NoPVz_DYLL_M2000to3000_TauTau.root"
*/
};

vector<TString> file_DYLL = {
	"output_data/saveFile_EE_NoPVz_DYLL_M10to50_EE_v1.root",
//	"output_data/saveFile_EE_NoPVz_DYLL_M10to50_EE_v2.root",
//	"output_data/saveFile_EE_NoPVz_DYLL_M10to50_EE_ext1v1.root",
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

vector<TString> file_top = {
	"output_data/saveFile_EE_NoPVz_ST_tW.root",
	"output_data/saveFile_EE_NoPVz_ST_tbarW.root",
	"output_data/saveFile_EE_NoPVz_ttbar_M0to700.root",
	"output_data/saveFile_EE_NoPVz_ttbar_M700to1000.root",
	"output_data/saveFile_EE_NoPVz_ttbar_M1000toInf.root"
};

vector<TString> file_Fake = {
	"output_data/saveFile_EE_NoPVz_WJetsToLNu_amcatnlo.root",
	"output_data/saveFile_EE_NoPVz_WJetsToLNu_amcatnlo_ext.root",
	"output_data/saveFile_EE_NoPVz_WJetsToLNu_amcatnlo_ext2v5.root"
};

vector<TString> file_QCD = {
	"output_data/saveFile_EE_NoPVz_QCDEMEnriched_Pt20to30.root",       
	"output_data/saveFile_EE_NoPVz_QCDEMEnriched_Pt30to50.root",       
	"output_data/saveFile_EE_NoPVz_QCDEMEnriched_Pt50to80.root",       
	"output_data/saveFile_EE_NoPVz_QCDEMEnriched_Pt50to80_ext1.root",  
	"output_data/saveFile_EE_NoPVz_QCDEMEnriched_Pt80to120.root",      
	"output_data/saveFile_EE_NoPVz_QCDEMEnriched_Pt80to120_ext1.root", 
	"output_data/saveFile_EE_NoPVz_QCDEMEnriched_Pt120to170.root",     
	"output_data/saveFile_EE_NoPVz_QCDEMEnriched_Pt120to170_ext1.root",
	"output_data/saveFile_EE_NoPVz_QCDEMEnriched_Pt170to300.root",     
	"output_data/saveFile_EE_NoPVz_QCDEMEnriched_Pt300toInf.root"      
};

enum Variable{
	INV_MASS,
	RAPIDITY,
	PT_LEAD,
	PT_SUB,
	ETA_LEAD,
	ETA_SUB,
	ERR
};

TH1D*GetHistogram(vector<TString> filesvector,TString variable);
vector<TString> GetPlotProperties(Variable var);
void MakePlots(Variable var);

void makeStackPlots()
{
	//gROOT->SetBatch(true);
	gStyle->SetOptStat(0);
	MakePlots(INV_MASS);
//	MakePlots(RAPIDITY);
//	MakePlots(PT_LEAD);
//	MakePlots(PT_SUB);

}

void MakePlots(Variable var)
{
	vector<TString> plotProperties = GetPlotProperties(var);

	// MC Signal
	TH1D*hDYLL = GetHistogram(file_DYLL,plotProperties.at(0));
	hDYLL->SetFillColor(kOrange-2);
	hDYLL->SetLineColor(kOrange+3);
	// Top quarks
	TH1D*hTops = GetHistogram(file_top,plotProperties.at(0));
	hTops->SetFillColor(kBlue+2);
	hTops->SetLineColor(kBlue+3);
	// W+Jets
//	TH1D*hFake = GetHistogram(file_Fake,plotProperties.at(0));
//	hFake->SetFillColor(kViolet+5);
//	hFake->SetLineColor(kViolet+3);
	// EW
	TH1D*hEW = GetHistogram(file_ew,plotProperties.at(0));
	hEW->SetFillColor(kRed+2);
	hEW->SetLineColor(kRed+4);
	// QCD
//	TH1D*hQCD = GetHistogram(file_QCD,plotProperties.at(0));
//	hQCD->SetFillColor(kGreen+2);
//	hQCD->SetLineColor(kGreen+4);

	// Total MC sum
	TString hSumName = "hSum";
	hSumName += plotProperties.at(0); 
	TH1D*hSum = (TH1D*)hDYLL->Clone(hSumName);
	hSum->Add(hTops);
//	hSum->Add(hFake);
	hSum->Add(hEW);
//	hSum->Add(hQCD);

	// Some samples can have negative bin contents due to weightings
	// Set any negative bins to zero
	int nBins = hDYLL->GetNbinsX();
	for(int i=1;i<nBins+1;i++){
		if(hDYLL->GetBinContent(i)<0) hDYLL->SetBinContent(i,0);
		if(hTops->GetBinContent(i)<0) hTops->SetBinContent(i,0);
//		if(hFake->GetBinContent(i)<0) hFake->SetBinContent(i,0);
		if(hEW->GetBinContent(i)<0) hEW->SetBinContent(i,0);
//		if(hQCD->GetBinContent(i)<0) hQCD->SetBinContent(i,0);
	}// end loop over bins

	// Place signal and background into a stack
	TString stackName = "hStack";
	stackName += plotProperties.at(0);
	THStack*hStack = new THStack(stackName,"");
//	hStack->Add(hQCD);
//	hStack->Add(hFake);
	hStack->Add(hEW);
	hStack->Add(hTops);
	hStack->Add(hDYLL);
	
	// Data
	TH1D*hData = GetHistogram(file_data,plotProperties.at(0));
	hData->SetMarkerStyle(20);
	hData->SetMarkerColor(kBlack);

	// Data over MC ratio
	TString ratioName = "ratio";
	ratioName += plotProperties.at(0);
	TH1D*hRatio = (TH1D*)hData->Clone(ratioName);
	hRatio->Divide(hSum);

	double binWidthLastBin = hDYLL->GetBinWidth(nBins);
	float x1 = hDYLL->GetBinLowEdge(0);
	float x2 = hDYLL->GetBinLowEdge(nBins)+binWidthLastBin;

	TLine*line = new TLine(x1,1,x2,1);
	line->SetLineColor(kRed);

	// Canvas
	TString canvasName = "canvas";
	canvasName += plotProperties.at(0);
	TCanvas*canvas = new TCanvas(canvasName,"",0,0,1000,1000);

	// Upper pad
	TPad*pad1 = new TPad("","",0,0.3,1,1.0);
	pad1->SetBottomMargin(0.03);
	pad1->SetGrid();
	pad1->SetLogy();
	pad1->SetTicks(1,1);
	pad1->Draw();
	pad1->cd();
	if(var==INV_MASS) pad1->SetLogx();
	hStack->Draw("hist");

	hStack->SetMinimum(0.01);
	hStack->SetMaximum(100000000);
	hStack->Draw("hist");
	hData->Draw("pe,same");

	TLegend*legend = new TLegend(0.65,0.9,0.9,0.75);
	legend->SetTextSize(0.02);
	legend->AddEntry(hData,"Data");
	legend->AddEntry(hDYLL,"#gamma^{*}/Z #rightarrow e^{-}e^{+}");
	legend->AddEntry(hTops,"t#bar{t}+tW+#bar{t}W");
	legend->AddEntry(hEW,"diboson + #gamma^{*}/Z #rightarrow #tau^{-}#tau^{+}");
//	legend->AddEntry(hFake,"W+Jets");
//	legend->AddEntry(hQCD,"QCD");
	legend->Draw("same");

	// Stack properties
	hStack->GetXaxis()->SetTitle(plotProperties.at(2));
	hStack->GetXaxis()->SetMoreLogLabels();
	hStack->GetXaxis()->SetNoExponent();
	hStack->SetTitle(plotProperties.at(1));
	hStack->SetMinimum(0.1);
	hStack->GetXaxis()->SetLabelSize(0);
	hStack->GetXaxis()->SetTitleSize(0);
	canvas->Update();
	hData->Draw("pe,same");
	legend->Draw("same");

	canvas->cd();

	// Lower pad
	TPad*pad2 = new TPad("","",0,0.05,1,0.3);	
	if(var==INV_MASS) pad2->SetLogx();
	pad2->SetTopMargin(0.03);
	pad2->SetBottomMargin(0.2);
	pad2->SetGrid();
	pad2->SetTicks(1,1);
	pad2->Draw();
	pad2->cd();

	hRatio->GetYaxis()->SetLabelSize(0.06);
	hRatio->GetYaxis()->SetTitleSize(0.08);
	hRatio->GetYaxis()->SetTitleOffset(0.3);
	hRatio->GetYaxis()->SetTitle("Data/MC");
	hRatio->GetXaxis()->SetLabelSize(0.1);
	hRatio->GetXaxis()->SetTitleSize(0.1);
	hRatio->GetXaxis()->SetNoExponent();
	hRatio->GetXaxis()->SetMoreLogLabels();
	hRatio->SetMaximum(1.3);
	hRatio->SetMinimum(0.7);
	hRatio->GetXaxis()->SetTitle(plotProperties.at(2));
	hRatio->Draw("PE");
	line->Draw("same");
	canvas->Update();

	TString saveName = "plots/dataVsMC_EE_NoPVz_WithDressed_NoFake_NoQCD_";
	saveName += plotProperties.at(0);
	saveName += ".png";
	canvas->SaveAs(saveName);
}

vector<TString> GetPlotProperties(Variable var)
{
	vector<TString> properties;
	vector<TString> errorVec = {"ERROR"};

	TString varTag;
	TString plotTitle;
	TString axisTitle;

	if(var==INV_MASS){
		varTag = "InvMass";
		plotTitle = "Dielectron invariant mass";
		axisTitle = "invariant mass [GeV]";		
	}
	else if(var==RAPIDITY){
		varTag = "Rapidity";
		plotTitle = "Dielectron rapidity";
		axisTitle = "rapidity";		
	}
	else if(var==PT_LEAD){
		varTag = "PtLead";
		plotTitle = "Leading electron p_{T}";
		axisTitle = "p_{T} [GeV]";		
	}
	else if(var==PT_SUB){
		varTag = "PtSub";
		plotTitle = "Sub-leading electron p_{T}";
		axisTitle = "p_{T} [GeV]";		
	}
	else if(var==ETA_LEAD){
		varTag = "EtaLead";
		plotTitle = "Leading electron #eta";
		axisTitle = "#eta";		
	}
	else if(var==ETA_SUB){
		varTag = "EtaSub";
		plotTitle = "Sub-leading electron #eta";
		axisTitle = "#eta";		
	}
	else{
		cout << "Variable not properly set in MakePlots()" << endl;
		return errorVec;
	}
	properties.push_back(varTag);
	properties.push_back(plotTitle);
	properties.push_back(axisTitle);

	return properties;	
}

TH1D*GetHistogram(vector<TString> filesvector,TString variable)
{
	int nFiles = filesvector.size();
	TString loadFile;
	TString loadHist;
	TFile*file_load[nFiles];
	TH1D*hMassArray[nFiles];
	TH1D*hMass;
	for(int i=0;i<nFiles;i++){
		loadFile = filesvector.at(i);
		loadHist = "hist";
		loadHist += variable;
		loadHist += "Reco";
		cout << "********************************************" << endl;
		cout << "Loading histogram: " << loadHist << endl;
		cout << "From file: " << loadFile << endl; 
		file_load[i] = new TFile(loadFile);
		if(!file_load[i]){
			cout << "ERROR loading file: " << loadFile << endl;
		}
		hMassArray[i] = (TH1D*)file_load[i]->Get(loadHist);
		if(!hMassArray[i]){
			cout << "ERROR loading histogram: " << loadHist << endl;
		}
		cout << "********************************************" << endl;
		cout << endl;
		if(i==0) hMass = (TH1D*)hMassArray[i]->Clone();
		else hMass->Add(hMassArray[i]);
	}// end loop over files

	return hMass;
}

