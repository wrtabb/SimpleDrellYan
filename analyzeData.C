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

TString base_directory = "root://xrootd-local.unl.edu///store/user/wtabb/DrellYan_13TeV_2016/v2p6/skims/skims_EE/";
vector<TString> files= {
	"crab_DoubleEG_RunB",		// 0
	"crab_DoubleEG_RunC",		// 1
	"crab_DoubleEG_RunD",		// 2
	"crab_DoubleEG_RunE",		// 3
	"crab_DoubleEG_RunF",		// 4
	"crab_DoubleEG_RunG",		// 5
	"crab_DoubleEG_RunHver2",	// 6
	"crab_DoubleEG_RunHver3",	// 7
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
	"DYLL_M2000to3000_EE"		// 18
};
TString treeName = "recoTree/DYTree";
vector<double> xSecVec = {
	1,1,1,1,1,1,1,1,//Data
	18610.0/3,        //DYLL_10to50 v1,v2,ext1v1 combined (NLO)
	1923.26,      //DYLL_50to100(NNLO)
	78.1258,      //DYLL_100to200(NNLO)
	2.73309,      //DYLL_200to400(NNLO) 
	0.142945,     //DYLL_400to500(NNLO)
	0.0809755,    //DYLL_500to700(NNLO)
	0.0125589,    //DYLL_700to800(NNLO)
	0.0105845,    //DYLL_800to1000(NNLO)
	0.00556507,   //DYLL_1000to1500(NNLO)
	0.000730495,  //DYLL_1500to2000(NNLO)
	0.00016844    //DYLL_2000to3000(NNLO)
};
int dataLuminosity = 35867;
const TString electronTrigger = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*";
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
TLorentzVector v1;
TLorentzVector v2;
void analyzeData(TString fileName)
{
	TH1::SetDefaultSumw2();

	// Get cross section for sample
	double xSec = 1.0;
	if(fileName.CompareTo(files.at(0))==0) xSec = xSecVec.at(0);
	else if(fileName.CompareTo(files.at(1))==0) xSec = xSecVec.at(1);
	else if(fileName.CompareTo(files.at(2))==0) xSec = xSecVec.at(2);
	else if(fileName.CompareTo(files.at(3))==0) xSec = xSecVec.at(3);
	else if(fileName.CompareTo(files.at(4))==0) xSec = xSecVec.at(4);
	else if(fileName.CompareTo(files.at(5))==0) xSec = xSecVec.at(5);
	else if(fileName.CompareTo(files.at(6))==0) xSec = xSecVec.at(6);
	else if(fileName.CompareTo(files.at(7))==0) xSec = xSecVec.at(7);
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
	chain->SetBranchAddress("Nelectrons",&Nelectrons,&b_Nelectrons);
	chain->SetBranchAddress("Electron_pT",&Electron_pT,&b_Electron_pT);
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
		if(Nelectrons<2) continue;

		// Check if event passes HLT cut
		TString trigName;
		int trigNameSize = pHLT_trigName->size();
		bool passHLT = false;
		for(int iHLT=0;iHLT<trigNameSize;iHLT++){
			trigName = pHLT_trigName->at(iHLT);
			if(trigName.CompareTo(electronTrigger) && HLT_trigFired[iHLT]==1){
				passHLT = true;
				break;
			} // end if trigName
		}// end loop over triggers

		if(!passHLT) continue;

		// Get Reco Electrons
		double ptLead = -1000;
		double ptSub = -1000;
		double etaLead = -1000;
		double etaSub = -1000;
		double phiLead = -1000;
		double phiSub = -1000;

		int idxLead = -1;
		int idxSub = -1;

		// Find lead pT electron
		for(int iEle=0;iEle<Nelectrons;iEle++){
			if(!Electron_passMediumID[iEle]) continue;
			if(Electron_pT[iEle] > ptLead){
				ptLead = Electron_pT[iEle];
				etaLead = Electron_eta[iEle];
				phiLead = Electron_phi[iEle];
				idxLead = iEle;
			}
		}// end lead pt Loop

		// Find subleading pT electron
		for(int iEle=0;iEle<Nelectrons;iEle++) {
			if(!Electron_passMediumID[iEle]) continue;
			if(Electron_pT[iEle] > ptSub && Electron_pT[iEle] < ptLead){
				ptSub = Electron_pT[iEle];
				etaSub = Electron_eta[iEle];
				phiSub = Electron_phi[iEle];
				idxSub = iEle;
			}
		}//end sub pt loop

		// If either lead or subleading electron not defined, skip to next event
		if(idxLead<0 || idxSub<0) continue;

		// Kinematic cuts
		if(abs(etaLead)>etaGapLow && abs(etaLead)<etaGapHigh) continue;
		if(abs(etaSub)>etaGapLow && abs(etaSub)<etaGapHigh) continue;
		if(abs(etaLead)>etaHigh||abs(etaSub)>etaHigh) continue;
		if(!((ptLead>ptLow && ptSub>ptHigh) || 
		     (ptLead>ptHigh && ptSub>ptLow))) continue;

		v1.SetPtEtaPhiM(ptLead,etaLead,phiLead,eMass);
		v2.SetPtEtaPhiM(ptSub,etaSub,phiSub,eMass);

		// dielectron invariant mass
		double invMassReco = (v1+v2).M();
		
		// dielectron rapidity
		double rapidity = (v1+v2).Rapidity();

		double sfWeight = 1.0;
		double pvzWeight = 1.0;
		double puWeight = 1.0;
		double prefireWeight = 1.0;
		// NOTE: Fakes do not get Scale Factors; need to ensure this is the case when it is implemented
		if(isMC){
			// Get Scale factors
			if(ptLead>ptBinHigh) ptLead = ptBinHigh;
			if(ptSub>ptBinHigh) ptSub = ptBinHigh;
			if(etaLead>etaBinHigh) etaLead = etaBinHigh;
			if(etaSub>etaBinHigh) etaSub = etaBinHigh;
			if(etaLead<etaBinLow) etaLead = etaBinLow;
			if(etaSub<etaBinLow) etaSub = etaBinLow;

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
			sfWeight = sfReco1*sfReco2*sfID1*sfID2*sfHLT;

			// Get gen weight
			genWeight = (GENEvt_weight/fabs(GENEvt_weight))/sumGenWeight;

			// PVz weight
			pvzWeight = hPVzSF->GetBinContent(hPVzSF->FindBin(PVz));

			// Pileup weight
			puWeight = hPileup->GetBinContent(hPileup->FindBin(nPileUp));

			// Prefire weight
			prefireWeight = _prefiringweight;
		}
		double weight = xSecWeight*genWeight*sfWeight*pvzWeight*puWeight;
		if(!isMC) weight = 1.0;
		hInvMass->Fill(invMassReco,weight);
		hRapidity->Fill(rapidity,weight);
		hPtLead->Fill(ptLead,weight);
		hPtSub->Fill(ptSub,weight);

	}// end loop over entries
	TString saveName = "output_data/saveFile";
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

