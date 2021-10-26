#include <TMath.h>
#include <TBranch.h>
#include <TString.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TH2F.h>
#include <TH1F.h>
#include <iostream>

//-----Functions-----//
double GetCrossSection(TString fileName);
bool IsSampleFake(TString fileName);
bool PassDileptonSelection(double eta1,double eta2,double pt1,double pt2,int idx1,int idx2);
vector<double> GetVariables(double eta1,double eta2,double pt1,double pt2,double phi1,
			    double phi2);
bool GetRecoLeptons(int &idxRecoLead, int &idxRecoSub);
bool GetHardLeptons(int &idxHard1,int &idxHard2);
std::vector<TLorentzVector> GetDressedLeptons(int &idxHardLead,int &idxHardSub);
void Counter(Long64_t event,Long64_t total);

TString base_directory = "root://xrootd-local.unl.edu///store/user/wtabb/DrellYan_13TeV_2016/v2p6/skims/skims_EE/";
vector<TString> files= {
	// Data
	"crab_DoubleEG_RunB",		// 0
	"crab_DoubleEG_RunC",		// 1
	"crab_DoubleEG_RunD",		// 2
	"crab_DoubleEG_RunE",		// 3
	"crab_DoubleEG_RunF",		// 4
	"crab_DoubleEG_RunG",		// 5
	"crab_DoubleEG_RunHver2",	// 6
	"crab_DoubleEG_RunHver3",	// 7

	// MC Signal
	"DYLL_M10to50_EE",		// 8
	"DYLL_M50to100_EE",		// 9
	"DYLL_M100to200_EE",		// 10
	"DYLL_M200to400_EE",		// 11
	"DYLL_M400to500_EE",		// 12
	"DYLL_M500to700_EE",		// 13
	"DYLL_M700to800_EE",		// 14
	"DYLL_M800to1000_EE",		// 15
	"DYLL_M1000to1500_EE",		// 16
	"DYLL_M1500to2000_EE",		// 17
	"DYLL_M2000to3000_EE",		// 18

	// Tops
	"ST_tW",			// 19
	"ST_tbarW",			// 20
	"ttbar_M0to700",		// 21
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
	"WJetsToLNu_amcatnlo",		// 38	
	"WJetsToLNu_amcatnlo_ext",	// 39	
	"WJetsToLNu_amcatnlo_ext2v5",	// 40	

	// QCD
	"QCDEMEnriched_Pt20to30",	// 41
	"QCDEMEnriched_Pt30to50",	// 42
	"QCDEMEnriched_Pt50to80",	// 43
	"QCDEMEnriched_Pt50to80_ext1",	// 44
	"QCDEMEnriched_Pt80to120",	// 45
	"QCDEMEnriched_Pt80to120_ext1",	// 46
	"QCDEMEnriched_Pt120to170",	// 47
	"QCDEMEnriched_Pt120to170_ext1",// 48
	"QCDEMEnriched_Pt170to300",	// 49
	"QCDEMEnriched_Pt300toInf"	// 50
};
TString treeName = "recoTree/DYTree";
vector<double> xSecVec = {
	1,1,1,1,1,1,1,1,	//Data
	18610.0/3,		//DYLL_10to50 v1,v2,ext1v1 combined (NLO)
	1923.26,		//DYLL_50to100(NNLO)
	78.1258,		//DYLL_100to200(NNLO)
	2.73309,		//DYLL_200to400(NNLO) 
	0.142945,		//DYLL_400to500(NNLO)
	0.0809755,		//DYLL_500to700(NNLO)
	0.0125589,		//DYLL_700to800(NNLO)
	0.0105845,		//DYLL_800to1000(NNLO)
	0.00556507,		//DYLL_1000to1500(NNLO)
	0.000730495,		//DYLL_1500to2000(NNLO)
	0.00016844,		//DYLL_2000to3000(NNLO)
	35.85,			//ST_tW
	35.85,			//ST_tbarW
	728.74,			//ttbar_M0to700
	76.605,			//ttbar_M700to1000
	20.578,			//ttbar_M1000toInf
	118.7,			//WW
	47.13,			//WZ
	16.523,			//ZZ
	18610.0/3.0,		//10to50 (NLO)
	1923.26,		//50to100 (NNLO)
	78.1258,		//100to200 (NNLO)
	2.73309,		//200to400 (NNLO)
	0.142945,		//400to500 (NNLO)
	0.0809755,		//500to700 (NNLO)
	0.0125589,		//700to800 (NNLO)
	0.0105845,		//800to1000 (NNLO)
	0.00556507,		//1000to1500 (NNLO)
	0.000730495,		//1500to2000 (NNLO)
	0.00016844,		//2000to3000 ((NNLO)
	61526.7,		//WJetsToLNu (NNLO)
	61526.7,		//WJetsToLNu_ext (NNLO)
	61526.7,		//WJetsToLNu_ext2v5 (NNLO)
	557600000*0.0096,	//QCDEMEnriched_Pt20to30 (LO)	
	136000000*0.073,	//QCDEMEnriched_Pt30to50 (LO)
	19800000*0.146,		//QCDEMEnriched_Pt50to80 (LO)
	19800000*0.146,		//QCDEMEnriched_Pt50to80_ext1 (LO)
	2800000*0.125,		//QCDEMEnriched_Pt80to120 (LO)
	2800000*0.125,		//QCDEMEnriched_Pt80to120_ext1 (LO)
	477000*0.132,		//QCDEMEnriched_Pt120to170 (LO)
	477000*0.132,		//QCDEMEnriched_Pt120to170_ext1 (LO)
	114000*0.165,		//QCDEMEnriched_Pt170to300 (LO)
	9000*0.15		//QCDEMEnriched_Pt300toInf (LO)
};
int dataLuminosity = 35867;
const TString electronTrigger = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*";
const double etaGapLow = 1.4442;
const double etaGapHigh = 1.566;
const double etaHigh = 2.4;
const double ptLow = 17;
const double ptHigh = 28;
const float dRMinCut = 0.3;
const double ptBinHigh = 499;
const double ptBinLow = 26;
const double etaBinLow = -2.5;
const double etaBinHigh = 2.5;
const double pi = TMath::Pi();
const int MPSIZE = 2000;
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
int Nelectrons;
double Electron_Energy[MPSIZE];  //no muon
double Electron_pT[MPSIZE];
double Electron_Px[MPSIZE];
double Electron_Py[MPSIZE];
double Electron_Pz[MPSIZE];
double Electron_eta[MPSIZE];
double Electron_phi[MPSIZE];
int Electron_charge[MPSIZE];
double Electron_etaSC[MPSIZE]; //no muon
double Electron_phiSC[MPSIZE]; //no muon
double Electron_dxy[MPSIZE];
double Electron_dz[MPSIZE];
double Electron_EnergySC[MPSIZE]; //no muon
double Electron_etSC[MPSIZE]; //no muon
bool Electron_passMediumID[MPSIZE];
double _prefiringweight;
double eMass = 0.000510998;
int HLT_ntrig;
int HLT_trigType[MPSIZE];
int HLT_trigFired[MPSIZE];
std::vector<std::string> HLT_trigName;
std::vector<std::string> *pHLT_trigName = &HLT_trigName;
int nVertices;
int nPileUp;
double PVz;
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
TBranch*b_PVz;

TBranch*b__prefiringweight;
TBranch*b__prefiringweightup;
TBranch*b__prefiringweightdown;
TBranch*b_HLT_ntrig;
TBranch*b_HLT_trigType;
TBranch*b_HLT_trigFired;
TBranch*b_Electron_Energy;
TBranch*b_Electron_pT;
TBranch*b_Electron_Px;
TBranch*b_Electron_Py;
TBranch*b_Electron_Pz;
TBranch*b_Electron_eta;
TBranch*b_Electron_phi;
TBranch*b_Electron_charge;
TBranch*b_Electron_etaSC;
TBranch*b_Electron_phiSC;
TBranch*b_Electron_dxy;
TBranch*b_Electron_dz;
TBranch*b_Electron_EnergySC;
TBranch*b_Electron_etSC;
TBranch*b_Electron_passMediumID;
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
double massbins2[] = {15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,
 50,52.5,55,57.5,60,62,64,66,68,70,72,74,76,78.5,81,83.5,86,88.5,91,93.5,96,98.5,101,103.5,
 106,108,110,112.5,115,117.5,120,123,126,129.5,133,137,141,145.5,150,155,160,165.5,171,178,
 185,192.5,200,210,220,231.5,243,258,273,296.5,320,350,380,410,440,475,510,555,600,650,700,
 765,830,915,1000,1250,1500,2250,3000};
int nMassBins2 = size(massbins2)-1;

void analyzeData(TString fileName)
{
	TH1::SetDefaultSumw2();

	// Get cross section for sample
	bool isFake = IsSampleFake(fileName);
        double xSec = GetCrossSection(fileName);

	TChain*chain;
	TH1D*hInvMassReco;
	TH1D*hRapidityReco;
	TH1D*hPtLeadReco;
	TH1D*hPtSubReco;

	TH1D*hInvMassHard;
	TH1D*hRapidityHard;
	TH1D*hPtLeadHard;
	TH1D*hPtSubHard;

	TH1D*hInvMassDressed;
	TH1D*hRapidityDressed;
	TH1D*hPtLeadDressed;
	TH1D*hPtSubDressed;

	TH2D*hMatrixInvMassHard;
	TH2D*hMatrixInvMassDressed;
	TH2D*hMatrixRapidityHard;
	TH2D*hMatrixRapidityDressed;

	TH2D*hDressedVsHard;
        TH2D*hDressedVsFSR;

	// Get histograms needed for weights and scale factors
	TFile*fRecoSF  = new TFile("data/Reco_SF.root");
	TFile*fMedIDSF = new TFile("data/MediumID_SF.root");
	TFile*fLeg2SF  = new TFile("data/Leg2_SF.root");
	TFile*fPileup  = new TFile("data/pileup.root");
	TFile*fPVzSF   = new TFile("data/PVz.root");

	TH2F*hRecoSF  = (TH2F*)fRecoSF->Get("EGamma_SF2D");
	TH2F*hMedIDSF = (TH2F*)fMedIDSF->Get("EGamma_SF2D");
	TH2F*hLeg2SF  = (TH2F*)fLeg2SF->Get("EGamma_SF2D");
	TH1F*hPileup  = (TH1F*)fPileup->Get("hPileupRatio");
	TH1D*hPVzSF   = (TH1D*)fPVzSF->Get("PVz_SF");

	// Load trees
	TString loadFile = base_directory+fileName;
	loadFile += ".root";
	chain = new TChain(treeName);
	chain->Add(loadFile);
	
	// Define Reco histograms
	TString histName = "hist";
	TString histNameInvMassReco = histName+"InvMassReco";
	TString histNameRapidityReco = histName+"RapidityReco";
	TString histNamePtLeadReco = histName+"PtLeadReco";
	TString histNamePtSubReco = histName+"PtSubReco";

	hInvMassReco = new TH1D(histNameInvMassReco,"",nMassBins,massbins);
	hRapidityReco = new TH1D(histNameRapidityReco,"",100,-2.5,2.5);
	hPtLeadReco = new TH1D(histNamePtLeadReco,"",100,0,500);
	hPtSubReco = new TH1D(histNamePtSubReco,"",100,0,500);

	// Define Hard histograms
	TString histNameInvMassHard = histName+"InvMassHard";
	TString histNameRapidityHard = histName+"RapidityHard";
	TString histNamePtLeadHard = histName+"PtLeadHard";
	TString histNamePtSubHard = histName+"PtSubHard";

	hInvMassHard = new TH1D(histNameInvMassHard,"",nMassBins,massbins);
	hRapidityHard = new TH1D(histNameRapidityHard,"",100,-2.5,2.5);
	hPtLeadHard = new TH1D(histNamePtLeadHard,"",100,0,500);
	hPtSubHard = new TH1D(histNamePtSubHard,"",100,0,500);

	// Define Dressed histograms
	TString histNameInvMassDressed = histName+"InvMassDressed";
	TString histNameRapidityDressed = histName+"RapidityDressed";
	TString histNamePtLeadDressed = histName+"PtLeadDressed";
	TString histNamePtSubDressed = histName+"PtSubDressed";

	hInvMassDressed = new TH1D(histNameInvMassDressed,"",nMassBins,massbins);
	hRapidityDressed = new TH1D(histNameRapidityDressed,"",100,-2.5,2.5);
	hPtLeadDressed = new TH1D(histNamePtLeadDressed,"",100,0,500);
	hPtSubDressed = new TH1D(histNamePtSubDressed,"",100,0,500);

	// Define migration matrices
	TString matrixName = "histMatrix";
	TString matrixNameInvMassHard = matrixName+"InvMassHard";
	TString matrixNameRapidityHard = matrixName+"RapidityHard";

	TString matrixNameInvMassDressed = matrixName+"InvMassDressed";
	TString matrixNameRapidityDressed = matrixName+"RapidityDressed";

	hMatrixInvMassHard = new 
		TH2D(matrixNameInvMassHard,"",nMassBins,massbins,nMassBins2,massbins2); 
	hMatrixInvMassDressed = new 
		TH2D(matrixNameInvMassDressed,"",nMassBins,massbins,nMassBins2,massbins2); 
	hMatrixRapidityHard = new 
		TH2D(matrixNameRapidityHard,"",100,-2.5,2.5,200,-2.5,2.5); 
	hMatrixRapidityDressed = new 
		TH2D(matrixNameRapidityDressed,"",100,-2.5,2.5,200,-2.5,2.5); 

	// Define dressed vs hard and fsr histograms
	TString DressedVsHard = "histDressedVsHard";
        TString DressedVsFSR  = "histDressedVsFSR";
        hDressedVsHard = new TH2D(DressedVsHard,"",nMassBins,massbins,nMassBins,massbins);
        hDressedVsFSR  = new TH2D(DressedVsFSR,"",nMassBins,massbins,nMassBins,massbins);;

	Long64_t nEntries = chain->GetEntries();
	cout << "Loading " << fileName << endl;
	cout << nEntries << " entries loaded. " << endl;

	bool isMC = true;
	TBranch*testBranch = 
		(TBranch*)chain->GetListOfBranches()->FindObject("GENEvt_weight");
	if(!testBranch) isMC = false;

	// Define branches
	chain->SetBranchAddress("Nelectrons",&Nelectrons,&b_Nelectrons);
	chain->SetBranchAddress("Electron_pT",&Electron_pT,&b_Electron_pT);
	chain->SetBranchAddress("Electron_Px",&Electron_Px,&b_Electron_Px);
	chain->SetBranchAddress("Electron_Py",&Electron_Py,&b_Electron_Py);
	chain->SetBranchAddress("Electron_Pz",&Electron_Pz,&b_Electron_Pz);
	chain->SetBranchAddress("Electron_Energy",&Electron_Energy,&b_Electron_Energy);
	chain->SetBranchAddress("Electron_eta",&Electron_eta,&b_Electron_eta);
	chain->SetBranchAddress("Electron_phi",&Electron_phi,&b_Electron_phi);
	chain->SetBranchAddress("Electron_passMediumID",&Electron_passMediumID,
				 &b_Electron_passMediumID);
	chain->SetBranchAddress("HLT_ntrig",&HLT_ntrig,&b_HLT_ntrig);
	chain->SetBranchAddress("HLT_trigType",&HLT_trigType,&b_HLT_trigType);
	chain->SetBranchAddress("HLT_trigFired",&HLT_trigFired,&b_HLT_trigFired);
	chain->SetBranchAddress("HLT_trigName",&pHLT_trigName);
	chain->SetBranchAddress("PVz",&PVz,&b_PVz);
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
		chain->SetBranchAddress("GenOthers_Px",&GenOthers_Px,&b_GenOthers_Px);
		chain->SetBranchAddress("GenOthers_Py",&GenOthers_Py,&b_GenOthers_Py);
		chain->SetBranchAddress("GenOthers_Pz",&GenOthers_Pz,&b_GenOthers_Pz);
		chain->SetBranchAddress("GenOthers_E",&GenOthers_E,&b_GenOthers_E);
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
	}// end isMC for gen weight sum calculation

	// Loop over events
	for(Long64_t iEntry=0;iEntry<nEntries;iEntry++){
		chain->GetEntry(iEntry);
		Counter(iEntry,nEntries);

		// Check if event passes HLT cut
		TString trigName;
		int trigNameSize = pHLT_trigName->size();
		bool passHLT = false;
		for(int iHLT=0;iHLT<trigNameSize;iHLT++){
			trigName = pHLT_trigName->at(iHLT);
			if(trigName.CompareTo(electronTrigger==0) && HLT_trigFired[iHLT]==1){
				passHLT = true;
				break;
			} // end if trigName
		}// end loop over triggers

		//-----Get Reconstructed Quantities-----//
		double invMassReco	= -1000;
		double rapidityReco	= -1000;
		double leadPtReco	= -1000;
		double subPtReco	= -1000;

		double ptRecoLead  = -1000;
		double ptRecoSub   = -1000;
		double etaRecoLead = -1000;
		double etaRecoSub  = -1000;
		double phiRecoLead = -1000;
		double phiRecoSub  = -1000;

		int idxRecoLead = -1;
		int idxRecoSub  = -1;

		// Find lead pT electron
		if(passHLT){
			bool recoLep = GetRecoLeptons(idxRecoLead,idxRecoSub);
			if(recoLep){
				ptRecoLead  = Electron_pT[idxRecoLead];
				ptRecoSub   = Electron_pT[idxRecoSub];
				etaRecoLead = Electron_eta[idxRecoLead];
				etaRecoSub  = Electron_eta[idxRecoSub];
				phiRecoLead = Electron_phi[idxRecoLead];
				phiRecoSub  = Electron_phi[idxRecoSub];
			}
		}

		// Determine if reco leptons pass selection
		bool passRecoSelection = PassDileptonSelection(etaRecoLead,etaRecoSub,
							       ptRecoLead,ptRecoSub,
							       idxRecoLead,idxRecoSub);
		// Get Reco Variables
		vector<double> recoVariables;
		recoVariables = GetVariables(etaRecoLead,etaRecoSub,ptRecoLead,ptRecoSub,
					     phiRecoLead,phiRecoSub);


		if(passRecoSelection){
			invMassReco	= recoVariables.at(0);
			rapidityReco	= recoVariables.at(1);
			leadPtReco	= recoVariables.at(2);
			subPtReco	= recoVariables.at(3);
		}

		//-----Get Hard Process Quantities-----//
		double invMassHard	= -1000;
		double rapidityHard	= -1000;
		double leadPtHard	= -1000;
		double subPtHard	= -1000;

		double ptHardLead  = -1000;
		double ptHardSub   = -1000;
		double etaHardLead = -1000;
		double etaHardSub  = -1000;
		double phiHardLead = -1000;
		double phiHardSub  = -1000;

		int idxHardLead = -1;
		int idxHardSub  = -1;

		bool hardLep = GetHardLeptons(idxHardLead,idxHardSub);
		if(hardLep){
			ptHardLead  = GENLepton_pT[idxHardLead];
			ptHardSub   = GENLepton_pT[idxHardSub];
			etaHardLead = GENLepton_eta[idxHardLead];
			etaHardSub  = GENLepton_eta[idxHardSub];
			phiHardLead = GENLepton_phi[idxHardLead];
			phiHardSub  = GENLepton_phi[idxHardSub];
		}
		bool passHardSelection = false;
		if(ptHardLead >=0 || ptHardSub >= 0 || etaHardLead > -3 ||
		   etaHardSub > -3 || phiHardLead > -100 || phiHardSub > -100){
			passHardSelection = PassDileptonSelection(etaHardLead,etaHardSub,
					   		          ptHardLead,ptHardSub,
							          idxHardLead,idxHardSub);
		}

		// Get Hard Variables
		vector<double> hardVariables;
		hardVariables = GetVariables(etaHardLead,etaHardSub,ptHardLead,ptHardSub,
					     phiHardLead,phiHardSub);


		if(passHardSelection){
			invMassHard	= hardVariables.at(0);
			rapidityHard	= hardVariables.at(1);
			leadPtHard	= hardVariables.at(2);
			subPtHard	= hardVariables.at(3);
		}

		//-----Get Dressed Quantities-----//
		double invMassDressed	= -1000;
		double rapidityDressed	= -1000;
		double leadPtDressed	= -1000;
		double subPtDressed	= -1000;

		int idxDressedLead = -1;
		int idxDressedSub  = -1;

		vector<TLorentzVector> dressedLeptons = 
			GetDressedLeptons(idxDressedLead,idxDressedSub);

		double ptDressedLead  = -1000;
		double ptDressedSub   = -1000;
		double etaDressedLead = -1000;
		double etaDressedSub  = -1000;
		double phiDressedLead = -1000;
		double phiDressedSub  = -1000;

		bool passDressedSelection = PassDileptonSelection(etaDressedLead,
								  etaDressedSub,
							          ptDressedLead,
								  ptDressedSub,
							          idxDressedLead,
								  idxDressedSub);

		ptDressedLead  = dressedLeptons.at(0).Pt(); 
		ptDressedSub   = dressedLeptons.at(1).Pt();
		etaDressedLead = dressedLeptons.at(0).Eta();
		etaDressedSub  = dressedLeptons.at(1).Eta();
		phiDressedLead = dressedLeptons.at(0).Phi();
		phiDressedSub  = dressedLeptons.at(1).Phi();

		// Get Dressed Variables
		vector<double> dressedVariables;
		dressedVariables = GetVariables(etaDressedLead,etaDressedSub,ptDressedLead,
					        ptDressedSub,phiDressedLead,phiDressedSub);


		if(passDressedSelection){
			invMassDressed	= dressedVariables.at(0);
			rapidityDressed	= dressedVariables.at(1);
			leadPtDressed	= dressedVariables.at(2);
			subPtDressed	= dressedVariables.at(3);
		}

		//-----Get Weights-----//
		double sfWeight = 1.0;
		double pvzWeight = 1.0;
		double puWeight = 1.0;
		double prefireWeight = 1.0;
		double genWeight = 1.0;

		double pt1  = ptRecoLead;
		double pt2  = ptRecoSub;
		double eta1 = etaRecoLead;
		double eta2 = etaRecoSub;

		if(isMC){
			// Get Scale factors
			if(pt1>ptBinHigh) pt1 = ptBinHigh;
			if(pt2>ptBinHigh) pt2 = ptBinHigh;
			if(pt1<ptBinLow) pt1 = ptBinLow;
			if(pt2<ptBinLow) pt2 = ptBinLow;

			if(!isFake){
			double sfReco1 = 
				hRecoSF->GetBinContent(hRecoSF->FindBin(eta1,pt1));
			double sfReco2 = 
				hRecoSF->GetBinContent(hRecoSF->FindBin(eta2,pt2));
			double sfID1 = 
				hMedIDSF->GetBinContent(hMedIDSF->FindBin(eta1,pt1));
			double sfID2 = 
				hMedIDSF->GetBinContent(hMedIDSF->FindBin(eta2,pt2));
			double sfHLT1 =hLeg2SF->GetBinContent(hLeg2SF->FindBin(eta1,pt1));
			double sfHLT2 = hLeg2SF->GetBinContent(hLeg2SF->FindBin(eta2,pt2));
			sfWeight = sfReco1*sfReco2*sfID1*sfID2*sfHLT1*sfHLT2;
			}// end isFake

			// Get gen weight
			genWeight = (GENEvt_weight/fabs(GENEvt_weight))/sumGenWeight;

			// PVz weight
//			pvzWeight = hPVzSF->GetBinContent(hPVzSF->FindBin(PVz));

			// Pileup weight
			puWeight = hPileup->GetBinContent(hPileup->FindBin(nPileUp));

			// Prefire weight
			prefireWeight = _prefiringweight;
		}// end isMC for weight calculation

		double recoWeight = 
			xSecWeight*genWeight*sfWeight*pvzWeight*puWeight*prefireWeight;
		double hardWeight = 
			xSecWeight*genWeight*pvzWeight*puWeight*prefireWeight;

		if(!isMC) recoWeight = 1.0;

		// Fill reco histograms
		hInvMassReco->Fill(invMassReco,recoWeight);
		hRapidityReco->Fill(rapidityReco,recoWeight);
		hPtLeadReco->Fill(leadPtReco,recoWeight);
		hPtSubReco->Fill(subPtReco,recoWeight);

		// Fill gen-level from hard process histograms
		hInvMassHard->Fill(invMassHard,hardWeight);
		hRapidityHard->Fill(rapidityHard,hardWeight);
		hPtLeadHard->Fill(leadPtHard,hardWeight);
		hPtSubHard->Fill(subPtHard,hardWeight);

		// Fill dressed histograms
		hInvMassDressed->Fill(invMassDressed,hardWeight);
		hRapidityDressed->Fill(rapidityDressed,hardWeight);
		hPtLeadDressed->Fill(leadPtDressed,hardWeight);
		hPtSubDressed->Fill(subPtDressed,hardWeight);

		// Fill migration matrices
		hMatrixInvMassHard->Fill(invMassHard,invMassReco,recoWeight);
		hMatrixInvMassHard->Fill(invMassHard,0.0,recoWeight*(1-sfWeight));
		hMatrixInvMassDressed->Fill(invMassDressed,invMassReco,recoWeight); 
		hMatrixInvMassDressed->Fill(invMassDressed,0.0,recoWeight*(1-sfWeight)); 
		hMatrixRapidityHard->Fill(rapidityHard,rapidityReco,recoWeight); 
		hMatrixRapidityHard->Fill(rapidityHard,0.0,recoWeight*(1-sfWeight)); 
		hMatrixRapidityDressed->Fill(rapidityDressed,rapidityReco,recoWeight);
		hMatrixRapidityDressed->Fill(rapidityDressed,rapidityReco,recoWeight*(1-sfWeight));

		// Fill dressedVs histograms
		hDressedVsHard->Fill(invMassDressed,invMassHard,hardWeight);
		//hDressedVsFSR ->Fill(invMassDressed,invMassFSR);
		
	}// end loop over entries

	// Save results to output file
	TString saveName = "output_data/saveFile_EE_NoPVz_WithDressed_";
	saveName += fileName;
	saveName += ".root";
	TFile*file;
	file = new TFile(saveName,"recreate");
	hInvMassReco->Write();
	hRapidityReco->Write();
	hPtLeadReco->Write();
	hPtSubReco->Write();
	hInvMassHard->Write();
	hRapidityHard->Write();
	hPtLeadHard->Write();
	hPtSubHard->Write();
	hInvMassDressed->Write();
	hRapidityDressed->Write();
	hPtLeadDressed->Write();
	hPtSubDressed->Write();
	hMatrixInvMassHard->Write();
	hMatrixInvMassDressed->Write();
	hMatrixRapidityHard->Write();
	hMatrixRapidityDressed->Write();	
	hDressedVsHard->Write();
	//hDressedVsFSR->Write();
	file->Close();
}


bool PassDileptonSelection(double eta1,double eta2,double pt1,double pt2,int idx1,int idx2)
{
	if(idx1<0 || idx2<0) return false;
	if(abs(eta1)>etaGapLow && abs(eta1)<etaGapHigh) return false;
	if(abs(eta2)>etaGapLow && abs(eta2)<etaGapHigh) return false;
	if(abs(eta1)>etaHigh||abs(eta2)>etaHigh) return false;
	if(pt1>pt2 && (pt1<ptHigh || pt2<ptLow)) return false;
	if(pt2>pt1 && (pt2<ptHigh || pt1<ptLow)) return false;

	return true;
}

vector<double> GetVariables(double eta1,double eta2,double pt1,double pt2,double phi1,
			    double phi2)
{
	TLorentzVector v1;
	TLorentzVector v2;
	
	v1.SetPtEtaPhiM(pt1,eta1,phi1,eMass);
	v2.SetPtEtaPhiM(pt2,eta2,phi2,eMass);

	// dielectron invariant mass
	double invMass = (v1+v2).M();
	
	// dielectron rapidity
	double rapidity = (v1+v2).Rapidity();

	double ptLead;
	double ptSub;

	if(pt1>pt2){
		ptLead = pt1;
		ptSub  = pt2;
	}
	else{
		ptLead = pt2;
		ptSub  = pt1;
	}

	vector<double> variableReturn = {
		invMass,	// 0
		rapidity,	// 1
		ptLead,		// 2
		ptSub		// 3
	};

	return variableReturn;
}

bool GetRecoLeptons(int &idxRecoLead, int &idxRecoSub)
{
	int nDileptons = 0;
	for(int iEle=0;iEle<Nelectrons;iEle++){
		if(!Electron_passMediumID[iEle]) continue;
		for(int jEle=iEle+1;jEle<Nelectrons;jEle++){
			if(!Electron_passMediumID[jEle]) continue;
			if(Electron_pT[iEle] > Electron_pT[jEle]){
				idxRecoLead = iEle;
				idxRecoSub = jEle;
			}// end if iEle is leading electron
			else{
				idxRecoLead = jEle;
				idxRecoSub = iEle;
			}// end if jEle is leading electron
			nDileptons++;
		}// end inner reco lepton loop
	}// end reco lepton Loop
	if(nDileptons==1) return true;
	else return false;
}// end GetRecoLeptons()

bool GetHardLeptons(int &idxHardLead,int &idxHardSub)
{
	int nHardLeptons = 0;
	for(int iLep=0;iLep<GENnPair;iLep++){
		for(int jLep=iLep+1;jLep<GENnPair;jLep++){
			if(!(abs(GENLepton_ID[iLep])==11 && abs(GENLepton_ID[jLep])==11))
				continue;
			if(GENLepton_ID[iLep]*GENLepton_ID[jLep]>0) continue;
			if(GENLepton_isHardProcess[iLep]==1 && 
			   GENLepton_isHardProcess[jLep]==1){
				if(GENLepton_pT[iLep] > GENLepton_pT[jLep]){
					idxHardLead = iLep;
					idxHardSub = jLep;
					nHardLeptons++;
				}// end if iLep is leading electron
				else{
					idxHardLead = jLep;
					idxHardSub = iLep;
					nHardLeptons++;
				}// end if jLep is leading electron
			}// end if hard process
		}//end inner loop over gen leptons
	}//end outer loop over gen leptons

	if(nHardLeptons==0){
		return false;
	}
	if(nHardLeptons>1){
		cout << "Too many electrons found from hard process" << endl;
		cout << "Very odd; this shouldn't happen" << endl;
		return false;
	}
	return true;
}// end GetHardLeptons()

double GetCrossSection(TString fileName)
{
        int nFiles = files.size();
        double xsec = 1.0;
        for(int i=0;i<nFiles;i++){
                if(fileName.CompareTo(files.at(i))==0) xsec = xSecVec.at(i);
        }// end loop over possible samples      

        return xsec;
}// end GetCrossSection()

bool IsSampleFake(TString fileName)
{
        for(int i=38;i<41;i++){
                if(fileName.CompareTo(files.at(i))==0) return true;
        }

        return false;
}// end IsSampleFake()

std::vector<TLorentzVector> GetDressedLeptons(int &idxDressedLead,int &idxDressedSub)
{
	TLorentzVector dressed1;
	TLorentzVector dressed2;

	// loop over electrons and select two post-fsr gen-level electrons
	for(int iLep=0;iLep<GENnPair;iLep++){
		for(int jLep=iLep+1;jLep<GENnPair;jLep++){
			
			// require that they both be electrons
			if(!(abs(GENLepton_ID[iLep])==11 && abs(GENLepton_ID[jLep])==11))
				continue;
		
			// require that they have opposite charge
			if(GENLepton_ID[iLep]*GENLepton_ID[jLep]>0) continue;

			// require that they be post-fsr
			if(GENLepton_fromHardProcessFinalState[iLep]==1 && 
			   GENLepton_fromHardProcessFinalState[jLep]==1){

				// determine which electron is lead and sub-lead
				if(GENLepton_pT[iLep] > GENLepton_pT[jLep]){
					idxDressedLead = iLep;
					idxDressedSub = jLep;
				}// end if iLep is leading electron
				else{
					idxDressedLead = jLep;
					idxDressedSub = iLep;
				}// end if jLep is leading electron
			}// end if hard process
		}//end inner loop over gen leptons
	}//end outer loop over gen leptons

	double px1 = GENLepton_Px[idxDressedLead];
	double px2 = GENLepton_Px[idxDressedSub];
	double py1 = GENLepton_Py[idxDressedLead];
	double py2 = GENLepton_Py[idxDressedSub];
	double pz1 = GENLepton_Pz[idxDressedLead];
	double pz2 = GENLepton_Pz[idxDressedSub];
	double E1 = GENLepton_E[idxDressedLead];
	double E2 = GENLepton_E[idxDressedSub];
	dressed1.SetPxPyPzE(px1,py1,pz1,E1);
	dressed2.SetPxPyPzE(px2,py2,pz2,E2);

	double eta1 = GENLepton_eta[idxDressedLead];
	double eta2 = GENLepton_eta[idxDressedSub];
	double phi1 = GENLepton_phi[idxDressedLead];
	double phi2 = GENLepton_phi[idxDressedSub];

	double dRMin = 0.1;
	TLorentzVector phoVec;
	double etaPho,phiPho;
	double etaDiff1,phiDiff1;

	double etaDiff2,phiDiff2;
	double dR1Squared,dR1;
	double dR2Squared,dR2;
	double pxPho,pyPho,pzPho,EPho;

	// Photon loop
	for(int iPho=0;iPho<nGenOthers;iPho++){
		
		// require that they be photons and are prompt final state
		if(abs(GenOthers_ID[iPho])!=22 ||
		   GenOthers_isPromptFinalState[iPho]!=1) continue;

			// location of photon
			etaPho = GenOthers_eta[iPho];
			phiPho = GenOthers_phi[iPho];

			// find location related to leading electron
			etaDiff1 = eta1-etaPho;
			phiDiff1 = phi1-phiPho;
			dR1Squared = etaDiff1*etaDiff1+phiDiff1*phiDiff1;
			dR1 = sqrt(dR1Squared);

			// find location related to subleading electron
			etaDiff2 = eta2-etaPho;
			phiDiff2 = phi2-phiPho;
			dR2Squared = etaDiff2*etaDiff2+phiDiff2*phiDiff2;
			dR2 = sqrt(dR2Squared);

			// place current photon into a lorentz vector
			pxPho = GenOthers_Px[iPho];
			pyPho = GenOthers_Py[iPho];
			pzPho = GenOthers_Pz[iPho];
			EPho  = GenOthers_E[iPho];
			phoVec.SetPxPyPzE(pxPho,pyPho,pzPho,EPho);

			// only keep photon if it is within the cone dRMin
			if(dR1>dRMin && dR2>dRMin) continue;

			// associate photon with whichever electron it is closest to 
			if(dR1<dR2){
				dressed1 += phoVec;				
			}
			else{
				dressed2 += phoVec;				
			}
	}// end loop over photons	

	vector<TLorentzVector> returnVector;
	returnVector.push_back(dressed1);
	returnVector.push_back(dressed2);

	return returnVector;
}// end GetDressedLeptons()

void Counter(Long64_t event,Long64_t total)
{
	int P = 100*(event)/(total);
	if(event%(total/100)==0)
		cout << P << "%" << endl;
	 return;
}
