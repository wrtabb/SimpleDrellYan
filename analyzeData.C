#include <TMath.h>
#include <TBranch.h>
#include <TString.h>
#include <TFile.h>
#include <TH1D.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TH2F.h>
#include <TH1F.h>
#include <iostream>

TString base_directory = "root://xrootd-local.unl.edu///store/user/wtabb/DrellYan_13TeV_2016/v2p6/skims/skims_MuMu/";

vector<TString> files= {
	// Data
	"SingleMuon_Run2016B",		// 0
	"SingleMuon_Run2016C",		// 1
	"SingleMuon_Run2016D",		// 2
	"SingleMuon_Run2016E",		// 3
	"SingleMuon_Run2016F",		// 4
	"SingleMuon_Run2016G",		// 5
	"SingleMuon_Run2016Hver2",	// 6
	"SingleMuon_Run2016Hver3",	// 7

	// MC Signal
	"DYLL_M10to50_MuMu",		// 8
	"DYLL_M50to100_MuMu",		// 9
	"DYLL_M100to200_MuMu",		// 10
	"DYLL_M200to400_MuMu",		// 11
	"DYLL_M400to500_MuMu",		// 12
	"DYLL_M500to700_MuMu",		// 13
	"DYLL_M700to800_MuMu",		// 14
	"DYLL_M800to1000_MuMu",		// 15
	"DYLL_M1000to1500_MuMu",		// 16
	"DYLL_M1500to2000_MuMu",		// 17
	"DYLL_M2000to3000_MuMu",		// 18

	// Tops
	"ST_tW",			// 19
	"ST_tbarW",			// 20
	"ttbar_truncated_M0To700",	// 21
	"ttbar_M700to1000",		// 22
	"ttbar_M1000toInf",		// 23

	// EW
	"WW",				// 24
	"WZ",				// 25
	"ZZ",				// 26
	"DYLL_M10to50_TauTau",		// 27
	"DYLL_M50to100_TauTau",		// 28
	"DYLL_M100to200_TauTau",	// 29
	"DYLL_M200to400_TauTau",	// 30
	"DYLL_M400to500_TauTau",	// 31
	"DYLL_M500to700_TauTau",	// 32
	"DYLL_M700to800_TauTau",	// 33
	"DYLL_M800to1000_TauTau",	// 34
	"DYLL_M1000to1500_TauTau",	// 35
	"DYLL_M1500to2000_TauTau",	// 36
	"DYLL_M2000to3000_TauTau",	// 37

	// Fakes
	"WJetsToLNu_amcatnlo"		// 38	
};
TString treeName = "recoTree/DYTree";
vector<double> xSecVec = {
	1,1,1,1,1,1,1,1,//Data
	18610.0/3,      //DYLL_10to50 v1,v2,ext1v1 combined (NLO)
	1923.26,	//DYLL_50to100(NNLO)
	78.1258,	//DYLL_100to200(NNLO)
	2.73309,	//DYLL_200to400(NNLO) 
	0.142945,	//DYLL_400to500(NNLO)
	0.0809755,	//DYLL_500to700(NNLO)
	0.0125589,	//DYLL_700to800(NNLO)
	0.0105845,	//DYLL_800to1000(NNLO)
	0.00556507,	//DYLL_1000to1500(NNLO)
	0.000730495,	//DYLL_1500to2000(NNLO)
	0.00016844,	//DYLL_2000to3000(NNLO)
	35.85,		//ST_tW
	35.85,		//ST_tbarW
	728.74,		//ttbar_M0to700
	76.605,		//ttbar_M700to1000
	20.578,		//ttbar_M1000toInf
	118.7,		//WW
	47.13,		//WZ
	16.523,		//ZZ
	18610.0/3.0,	//10to50 (NLO)
	1923.26,	//50to100 (NNLO)
	78.1258,	//100to200 (NNLO)
	2.73309,	//200to400 (NNLO)
	0.142945,	//400to500 (NNLO)
	0.0809755,	//500to700 (NNLO)
	0.0125589,	//700to800 (NNLO)
	0.0105845,	//800to1000 (NNLO)
	0.00556507,	//1000to1500 (NNLO)
	0.000730495,	//1500to2000 (NNLO)
	0.00016844,	//2000to3000 ((NNLO)
	61526.7		//WJetsToLNu (NNLO)
};
int dataLuminosity = 35867;
const TString muonTrigger1 = "HLT_IsoMu24_v*";
const TString muonTrigger2 = "HLT_IsoTkMu24_v*";
const double etaGapLow = 1.4442;
const double etaGapHigh = 1.566;
const double etaHigh = 2.4;
const double ptLow = 17;
const double ptHigh = 28;
const float dRMinCut = 0.3;
const double ptBinHigh = 500.0;
const double etaBinLow = -2.5;
const double etaBinHigh = 2.5;
const double pi = TMath::Pi();
const int MPSIZE = 5000;
int GENnPair;//number of gen leptons per event
double GENEvt_weight;
double GENLepton_phi[MPSIZE];
double GENLepton_eta[MPSIZE];
double GENLepton_pT[MPSIZE];
double GENLepton_Px[MPSIZE];
double GENLepton_Py[MPSIZE];
double GENLepton_Pz[MPSIZE];
double GENLepton_E[MPSIZE];
int GENLepton_ID[MPSIZE];
bool GENLepton_isHardProcess[MPSIZE];
bool GENLepton_fromHardProcessFinalState[MPSIZE];
int nGenOthers;
double GenOthers_phi[MPSIZE];
double GenOthers_eta[MPSIZE];
double GenOthers_pT[MPSIZE];
double GenOthers_Px[MPSIZE];
double GenOthers_Py[MPSIZE];
double GenOthers_Pz[MPSIZE];
double GenOthers_E[MPSIZE];
int GenOthers_ID[MPSIZE];
int GenOthers_isHardProcess[MPSIZE];
int GenOthers_isPromptFinalState[MPSIZE];
int nMuon;
int Nmuons;
double Muon_pT[MPSIZE];
double Muon_Px[MPSIZE];
double Muon_Py[MPSIZE];
double Muon_Pz[MPSIZE];
double Muon_eta[MPSIZE];
double Muon_phi[MPSIZE];
int Muon_charge[MPSIZE];
double Muon_dxy[MPSIZE];
double Muon_dz[MPSIZE];
bool Muon_passTightID[MPSIZE];
double Muon_PfChargedHadronIsoR04[MPSIZE];
double Muon_PfNeutralHadronIsoR04[MPSIZE];
double Muon_PfGammaIsoR04[MPSIZE];
double Muon_PFSumPUIsoR04[MPSIZE];
double Muon_trkiso[MPSIZE];
double _prefiringweight;
double muMass = 0.1056583715;
int HLT_ntrig;
int HLT_trigType[MPSIZE];
int HLT_trigFired[MPSIZE];
std::vector<std::string> HLT_trigName;
std::vector<std::string> *pHLT_trigName = &HLT_trigName;
int nVertices;
int nPileUp;
std::vector<double> vtxTrkCkt1Pt;
std::vector<double>*pvtxTrkCkt1Pt = &vtxTrkCkt1Pt;

std::vector<double> vtxTrkCkt2Pt;
std::vector<double>*pvtxTrkCkt2Pt = &vtxTrkCkt2Pt;

std::vector<double> vtxTrkChi2;
std::vector<double>*pvtxTrkChi2 = &vtxTrkChi2;

std::vector<double> vtxTrkNdof;
std::vector<double>*pvtxTrkNdof = &vtxTrkNdof;
TBranch*b_runNum;
TBranch*b_evtNum;
TBranch*b_lumiBlock;
TBranch*b_PUweight;
TBranch*b_Nelectrons;
TBranch*b_nVertices;
TBranch*b_nPileUp;

TBranch*b__prefiringweight;
TBranch*b__prefiringweightup;
TBranch*b__prefiringweightdown;
TBranch*b_HLT_ntrig;
TBranch*b_HLT_trigType;
TBranch*b_HLT_trigFired;
TBranch*b_nMuon;
TBranch*b_Nmuons;
TBranch*b_Muon_pT;
TBranch*b_Muon_Px;
TBranch*b_Muon_Py;
TBranch*b_Muon_Pz;
TBranch*b_Muon_eta;
TBranch*b_Muon_phi;
TBranch*b_Muon_charge;
TBranch*b_Muon_dxy;
TBranch*b_Muon_dz;
TBranch*b_Muon_passTightID;
TBranch*b_Muon_PfChargedHadronIsoR04;
TBranch*b_Muon_PfNeutralHadronIsoR04;
TBranch*b_Muon_PfGammaIsoR04;
TBranch*b_Muon_PFSumPUIsoR04;
TBranch*b_Muon_trkiso;
TBranch*b_GENnPair;
TBranch*b_GENLepton_phi;
TBranch*b_GENLepton_eta;
TBranch*b_GENLepton_pT;
TBranch*b_GENLepton_Px;
TBranch*b_GENLepton_Py;
TBranch*b_GENLepton_Pz;
TBranch*b_GENLepton_E;
TBranch*b_GENLepton_mother;
TBranch*b_GENLepton_mother_pT;
TBranch*b_GENLepton_charge;
TBranch*b_GENLepton_status;
TBranch*b_GENLepton_ID;
TBranch*b_GENLepton_isPrompt;
TBranch*b_GENLepton_isPromptFinalState;
TBranch*b_GENLepton_isTauDecayProduct;
TBranch*b_GENLepton_isPromptTauDecayProduct;
TBranch*b_GENLepton_isDirectPromptTauDecayProductFinalState;
TBranch*b_GENLepton_isHardProcess;
TBranch*b_GENLepton_isLastCopy;
TBranch*b_GENLepton_isLastCopyBeforeFSR;
TBranch*b_GENLepton_isPromptDecayed;
TBranch*b_GENLepton_isDecayedLeptonHadron;
TBranch*b_GENLepton_fromHardProcessBeforeFSR;
TBranch*b_GENLepton_fromHardProcessDecayed;
TBranch*b_GENLepton_fromHardProcessFinalState;
TBranch*b_GENLepton_isMostlyLikePythia6Status3;
TBranch*b_GENEvt_weight;
TBranch*b_GENEvt_QScale;
TBranch*b_GENEvt_x1;
TBranch*b_GENEvt_x2;
TBranch*b_GENEvt_alphaQCD;
TBranch*b_GENEvt_alphaQED;
TBranch*b_nGenOthers;
TBranch*b_GenOthers_phi;
TBranch*b_GenOthers_eta;
TBranch*b_GenOthers_pT;
TBranch*b_GenOthers_Px;
TBranch*b_GenOthers_Py;
TBranch*b_GenOthers_Pz;
TBranch*b_GenOthers_E ;
TBranch*b_GenOthers_ID;
TBranch*b_GenOthers_isHardProcess;
TBranch*b_GenOthers_isPromptFinalState;
double massbins[] = {
                15,
                20,
                25,
                30,
                35,
                40,
                45,
                50,
                55,
                60,
                64,
                68,
                72,
                76,
                81,
                86,
                91,
                96,
                101,
                106,
                110,
                115,
                120,
                126,
                133,
                141,
                150,
                160,
                171,
                185,
                200,
                220,
                243,
                273,
                320,
                380,
                440,
                510,
                600,
                700,
                830,
                1000,
                1500,
                3000
};
int nMassBins = size(massbins)-1;//43;
TLorentzVector v1;
TLorentzVector v2;
void analyzeData(TString fileName)
{
	TH1::SetDefaultSumw2();

	// Get cross section for sample
	double xSec = 1.0;
	bool isFake = false;

	// Data
	if(fileName.CompareTo(files.at(0))==0) xSec = xSecVec.at(0);
	else if(fileName.CompareTo(files.at(1))==0) xSec = xSecVec.at(1);
	else if(fileName.CompareTo(files.at(2))==0) xSec = xSecVec.at(2);
	else if(fileName.CompareTo(files.at(3))==0) xSec = xSecVec.at(3);
	else if(fileName.CompareTo(files.at(4))==0) xSec = xSecVec.at(4);
	else if(fileName.CompareTo(files.at(5))==0) xSec = xSecVec.at(5);
	else if(fileName.CompareTo(files.at(6))==0) xSec = xSecVec.at(6);
	else if(fileName.CompareTo(files.at(7))==0) xSec = xSecVec.at(7);
	// MC Signal
	else if(fileName.CompareTo(files.at(8))==0) xSec = xSecVec.at(8);
	else if(fileName.CompareTo(files.at(9))==0) xSec = xSecVec.at(9);
	else if(fileName.CompareTo(files.at(10))==0) xSec = xSecVec.at(10);
	else if(fileName.CompareTo(files.at(11))==0) xSec = xSecVec.at(11);
	else if(fileName.CompareTo(files.at(12))==0) xSec = xSecVec.at(12);
	else if(fileName.CompareTo(files.at(13))==0) xSec = xSecVec.at(13);
	else if(fileName.CompareTo(files.at(14))==0) xSec = xSecVec.at(14);
	else if(fileName.CompareTo(files.at(15))==0) xSec = xSecVec.at(15);
	else if(fileName.CompareTo(files.at(16))==0) xSec = xSecVec.at(16);
	else if(fileName.CompareTo(files.at(17))==0) xSec = xSecVec.at(17);
	else if(fileName.CompareTo(files.at(18))==0) xSec = xSecVec.at(18);
	//Tops
	else if(fileName.CompareTo(files.at(19))==0) xSec = xSecVec.at(19);
	else if(fileName.CompareTo(files.at(20))==0) xSec = xSecVec.at(20);
	else if(fileName.CompareTo(files.at(21))==0) xSec = xSecVec.at(21);
	else if(fileName.CompareTo(files.at(22))==0) xSec = xSecVec.at(22);
	else if(fileName.CompareTo(files.at(23))==0) xSec = xSecVec.at(23);
	// EW
	else if(fileName.CompareTo(files.at(24))==0) xSec = xSecVec.at(24);
	else if(fileName.CompareTo(files.at(25))==0) xSec = xSecVec.at(25);
	else if(fileName.CompareTo(files.at(26))==0) xSec = xSecVec.at(26);
	else if(fileName.CompareTo(files.at(27))==0) xSec = xSecVec.at(27);
	else if(fileName.CompareTo(files.at(28))==0) xSec = xSecVec.at(28);
	else if(fileName.CompareTo(files.at(29))==0) xSec = xSecVec.at(29);
	else if(fileName.CompareTo(files.at(30))==0) xSec = xSecVec.at(30);
	else if(fileName.CompareTo(files.at(31))==0) xSec = xSecVec.at(31);
	else if(fileName.CompareTo(files.at(32))==0) xSec = xSecVec.at(32);
	else if(fileName.CompareTo(files.at(33))==0) xSec = xSecVec.at(33);
	else if(fileName.CompareTo(files.at(34))==0) xSec = xSecVec.at(34);
	else if(fileName.CompareTo(files.at(35))==0) xSec = xSecVec.at(35);
	else if(fileName.CompareTo(files.at(36))==0) xSec = xSecVec.at(36);
	else if(fileName.CompareTo(files.at(37))==0) xSec = xSecVec.at(37);
	// Fakes
	else if(fileName.CompareTo(files.at(38))==0){
		xSec = xSecVec.at(38);
		isFake = true;
	}

	TChain*chain;
	TH1D*hInvMass;
	TH1D*hRapidity;
	TH1D*hPtLead;
	TH1D*hPtSub;

	// Get histograms needed for weights and scale factors
	TFile*fRecoSF  = new TFile("data/Reco_SF.root");
	TFile*fMedIDSF = new TFile("data/MediumID_SF.root");
	TFile*fLeg2SF  = new TFile("data/Leg2_SF.root");
	TFile*fPileup  = new TFile("data/pileup.root");

	TH2F*hRecoSF  = (TH2F*)fRecoSF->Get("EGamma_SF2D");
	TH2F*hMedIDSF = (TH2F*)fMedIDSF->Get("EGamma_SF2D");
	TH2F*hLeg2SF  = (TH2F*)fLeg2SF->Get("EGamma_SF2D");
	TH1F*hPileup  = (TH1F*)fPileup->Get("hPileupRatio");

	// Load trees
	cout << "Loadig tree from file: " << fileName << endl;
	cout << "From directory: " << base_directory << endl;
	TString loadFile = base_directory+fileName;
	loadFile += ".root";
	chain = new TChain(treeName);
	chain->Add(loadFile);
	
	// Define histograms
	TString histName = "hist";
	histName += fileName;	
	TString histNameInvMass = histName+"InvMass";
	TString histNameRapidity = histName+"Rapidity";
	TString histNamePtLead = histName+"PtLead";
	TString histNamePtSub = histName+"PtSub";

	hInvMass = new TH1D(histNameInvMass,"",nMassBins,massbins);
	hRapidity = new TH1D(histNameRapidity,"",100,-2.5,2.5);
	hPtLead = new TH1D(histNamePtLead,"",100,0,500);
	hPtSub = new TH1D(histNamePtSub,"",100,0,500);

	Long64_t nEntries = chain->GetEntries();
	cout << "Loading " << fileName << endl;
	cout << nEntries << " entries loaded. " << endl;

	bool isMC = true;
	TBranch*testBranch = 
		(TBranch*)chain->GetListOfBranches()->FindObject("GENEvt_weight");
	if(!testBranch) isMC = false;

	// Define branches
	chain->SetBranchAddress("nMuon",&nMuon,&b_nMuon);
	chain->SetBranchAddress("Muon_pT",&Muon_pT,&b_Muon_pT);
	chain->SetBranchAddress("Muon_eta",&Muon_eta,&b_Muon_eta);
	chain->SetBranchAddress("Muon_phi",&Muon_phi,&b_Muon_phi);
	chain->SetBranchAddress("Muon_passTightID",&Muon_passTightID,
				&b_Muon_passTightID);
	chain->SetBranchAddress("Muon_charge",&Muon_charge,&b_Muon_charge);
	chain->SetBranchAddress("Muon_PfChargedHadronIsoR04",
				&Muon_PfChargedHadronIsoR04,
				&b_Muon_PfChargedHadronIsoR04);
	chain->SetBranchAddress("Muon_PfNeutralHadronIsoR04",
				&Muon_PfNeutralHadronIsoR04,
				&b_Muon_PfNeutralHadronIsoR04);
	chain->SetBranchAddress("Muon_PfGammaIsoR04",
				&Muon_PfGammaIsoR04,
				&b_Muon_PfGammaIsoR04);
	chain->SetBranchAddress("Muon_PFSumPUIsoR04",
				&Muon_PFSumPUIsoR04,
				&b_Muon_PFSumPUIsoR04);
	chain->SetBranchAddress("HLT_ntrig",&HLT_ntrig,&b_HLT_ntrig);
	chain->SetBranchAddress("HLT_trigType",&HLT_trigType,&b_HLT_trigType);
	chain->SetBranchAddress("HLT_trigFired",&HLT_trigFired,&b_HLT_trigFired);
	chain->SetBranchAddress("HLT_trigName",&pHLT_trigName);
	chain->SetBranchAddress("nVertices",&nVertices,&b_nVertices);
	chain->SetBranchAddress("nPileUp",&nPileUp,&b_nPileUp);
	chain->SetBranchAddress("vtxTrkCkt1Pt",&pvtxTrkCkt1Pt);
	chain->SetBranchAddress("vtxTrkCkt2Pt",&pvtxTrkCkt2Pt);
	chain->SetBranchAddress("vtxTrkChi2",&pvtxTrkChi2);
	chain->SetBranchAddress("vtxTrkNdof",&pvtxTrkNdof);

	if(isMC){
		chain->SetBranchAddress("_prefiringweight", &_prefiringweight,
					&b__prefiringweight);
		chain->SetBranchAddress("GENnPair", &GENnPair, &b_GENnPair);
		chain->SetBranchAddress("GENLepton_eta", &GENLepton_eta, 
					&b_GENLepton_eta);
		chain->SetBranchAddress("GENLepton_phi",&GENLepton_phi, 
					&b_GENLepton_phi);
		chain->SetBranchAddress("GENLepton_pT",&GENLepton_pT, 
					&b_GENLepton_pT);
		chain->SetBranchAddress("GENLepton_ID",&GENLepton_ID, 
					&b_GENLepton_ID);
		chain->SetBranchAddress("GENLepton_isHardProcess",
					&GENLepton_isHardProcess,
					&b_GENLepton_isHardProcess);
		chain->SetBranchAddress("GENLepton_fromHardProcessFinalState",
					&GENLepton_fromHardProcessFinalState,
					&b_GENLepton_fromHardProcessFinalState);
		chain->SetBranchAddress("nGenOthers",&nGenOthers,&b_nGenOthers);
		chain->SetBranchAddress("GenOthers_eta",&GenOthers_eta,
					&b_GenOthers_eta);
		chain->SetBranchAddress("GenOthers_phi",&GenOthers_phi,
					&b_GenOthers_phi);
		chain->SetBranchAddress("GenOthers_pT",&GenOthers_pT,
					&b_GenOthers_pT);
		chain->SetBranchAddress("GenOthers_ID",&GenOthers_ID,
					&b_GenOthers_ID);
		chain->SetBranchAddress("GenOthers_isHardProcess",
					&GenOthers_isHardProcess,
					&b_GenOthers_isHardProcess);
		chain->SetBranchAddress("GenOthers_isPromptFinalState",
					&GenOthers_isPromptFinalState,
					&b_GenOthers_isPromptFinalState);
		chain->SetBranchAddress("GENEvt_weight",&GENEvt_weight,
					&b_GENEvt_weight);
	}// end if monte carlo

	// Load xSec weights
	double xSecWeight = dataLuminosity*xSec/1.0;

	// Find the gen weight sum
	double sumGenWeight = 0.0;
	double genWeight;

	if(isMC){
		Long64_t localEntry;
		for(Long64_t iGen=0;iGen<nEntries;iGen++){
			localEntry = chain->LoadTree(iGen);
			b_GENEvt_weight->GetEntry(localEntry);
			genWeight = GENEvt_weight/fabs(GENEvt_weight);
			sumGenWeight += genWeight;
		}
		if(sumGenWeight<0){
			cout << "Gen weight sum < 0 for sample " << fileName << endl;
		}
	}
	// Loop over events
	for(Long64_t iEntry=0;iEntry<nEntries;iEntry++){
		chain->GetEntry(iEntry);
		if(nMuon<2) continue;

		// Check if event passes HLT cut
		TString trigName;
		int trigNameSize = pHLT_trigName->size();
		bool passHLT = false;
		for(int iHLT=0;iHLT<trigNameSize;iHLT++){
			trigName = pHLT_trigName->at(iHLT);
			if(((trigName.CompareTo(muonTrigger1)==0)||
			    (trigName.CompareTo(muonTrigger2)==0)) && 
			     HLT_trigFired[iHLT]==1){
				passHLT = true;
				break;
			} // end if trigName
		}// end loop over triggers

		if(!passHLT) continue;

		// Get Reco Muons
		double ptLead = -1000;
		double ptSub = -1000;
		double etaLead = -1000;
		double etaSub = -1000;
		double phiLead = -1000;
		double phiSub = -1000;

		int idxLead = -1;
		int idxSub = -1;

		// Find lead pT muon
		// NOTE: Need to choose two muons by smallest vertex chi2
		// Choosing highest two pT is a temporary placeholder 
		for(int iMu=0;iMu<nMuon;iMu++){
			if(!Muon_passTightID[iMu]) continue;
			// need to add pfIso/pt still
			if(Muon_pT[iMu] > ptLead){
				ptLead = Muon_pT[iMu];
				etaLead = Muon_eta[iMu];
				phiLead = Muon_phi[iMu];
				idxLead = iMu;
			}
		}// end lead pt Loop

		// Find subleading pT muon
		for(int iMu=0;iMu<nMuon;iMu++) {
			if(!Muon_passTightID[iMu]) continue;
			if(Muon_pT[iMu] > ptSub && Muon_pT[iMu] < ptLead){
				ptSub = Muon_pT[iMu];
				etaSub = Muon_eta[iMu];
				phiSub = Muon_phi[iMu];
				idxSub = iMu;
			}
		}//end sub pt loop

		// If either lead or subleading muon not defined, skip to next event
		if(idxLead<0 || idxSub<0) continue;

		// Kinematic cuts
		if(abs(etaLead)>etaGapLow && abs(etaLead)<etaGapHigh) continue;
		if(abs(etaSub)>etaGapLow && abs(etaSub)<etaGapHigh) continue;
		if(abs(etaLead)>etaHigh||abs(etaSub)>etaHigh) continue;
		if(!(ptLead>ptHigh && ptSub>ptLow)) continue;

		v1.SetPtEtaPhiM(ptLead,etaLead,phiLead,muMass);
		v2.SetPtEtaPhiM(ptSub,etaSub,phiSub,muMass);

		// dimuon invariant mass
		double invMassReco = (v1+v2).M();
		
		// dimuon rapidity
		double rapidity = (v1+v2).Rapidity();

		double sfWeight = 1.0;
		double puWeight = 1.0;
		double prefireWeight = 1.0;

		// Need Rochester correction for muons for data and MC
		if(isMC){
			// Get Scale factors
			if(ptLead>ptBinHigh) ptLead = ptBinHigh;
			if(ptSub>ptBinHigh) ptSub = ptBinHigh;
			if(etaLead>etaBinHigh) etaLead = etaBinHigh;
			if(etaSub>etaBinHigh) etaSub = etaBinHigh;
			if(etaLead<etaBinLow) etaLead = etaBinLow;
			if(etaSub<etaBinLow) etaSub = etaBinLow;

			if(!isFake){
			double sfReco1 = 
				hRecoSF->GetBinContent(hRecoSF->FindBin(etaLead,ptLead));
			double sfReco2 = 
				hRecoSF->GetBinContent(hRecoSF->FindBin(etaSub,ptSub));
			double sfID1 = 
				hMedIDSF->GetBinContent(hMedIDSF->FindBin(etaLead,ptLead));
			double sfID2 = 
				hMedIDSF->GetBinContent(hMedIDSF->FindBin(etaSub,ptSub));
			double sfHLT = 
				(hLeg2SF->GetBinContent(hLeg2SF->FindBin(etaLead,ptLead)))*
				 (hLeg2SF->GetBinContent(hLeg2SF->FindBin(etaLead,ptLead)));
			// Need to acquire muon SFs before adding this
			sfWeight = 1.0;
			}// end isFake

			// Get gen weight
			genWeight = (GENEvt_weight/fabs(GENEvt_weight))/sumGenWeight;

			// Pileup weight
			puWeight = hPileup->GetBinContent(hPileup->FindBin(nPileUp));

			// Prefire weight
			prefireWeight = _prefiringweight;
		}

		double weight = xSecWeight*genWeight*sfWeight*puWeight;
		if(!isMC) weight = 1.0;
		hInvMass->Fill(invMassReco,weight);
		hRapidity->Fill(rapidity,weight);
		hPtLead->Fill(ptLead,weight);
		hPtSub->Fill(ptSub,weight);

	}// end loop over entries
	TString saveName = "output_data/saveFile_MuMu_";
	saveName += fileName;
	saveName += ".root";
	TFile*file;
	file = new TFile(saveName,"recreate");
	hInvMass->Write();
	hRapidity->Write();
	hPtLead->Write();
	hPtSub->Write();
	file->Close();

}

