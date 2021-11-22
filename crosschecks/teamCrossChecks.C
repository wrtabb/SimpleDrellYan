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
#include "files_DYEE_50toInf.h"

//-----Functions-----//
bool PassDileptonSelection(double eta1,double eta2,double pt1,double pt2);
vector<double> GetVariables(double eta1,double eta2,double pt1,double pt2,double phi1,
			    double phi2);
bool GetRecoLeptons(int &idxRecoLead, int &idxRecoSub);
bool GetHardLeptons(int &idxHard1,int &idxHard2);
void Counter(Long64_t event,Long64_t total);

//TString base_directory = "root://xrootd-local.unl.edu///store/user/wtabb/DrellYan_13TeV_2016/v2p6/DYLL_M50toInf/base/";
TString base_directory = "/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/wtabb/DrellYan_13TeV_2016/v2p6/DYLL_M50toInf/base/*.root";

TString treeName = "recoTree/DYTree";
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
int GENLepton_isHardProcess[MPSIZE];
int GENLepton_fromHardProcessFinalState[MPSIZE];
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

void teamCrossChecks()
{
	TH1::SetDefaultSumw2();

	TChain*chain;

	// Load trees
	chain = new TChain(treeName);
	int nFiles = files.size();
	chain->Add(base_directory);

//	for(int i=0;i<nFiles;i++){
//		TString loadFile = base_directory+files.at(i);
//		chain->Add(loadFile);
//	}
	
	// Define histograms
	TH1D*h_gen_el_pt     = new TH1D("h_gen_el_pt","",10000,0,10000);
	TH1D*h_gen_el_eta    = new TH1D("h_gen__el_eta","",200,-10,10);
	TH1D*h_gen_el_phi    = new TH1D("h_gen_el_phi","",80,-4,4);
	TH1D*h_gen_diEl_mass = new TH1D("h_gen_diEl_mass","",10000,0,10000);
	TH1D*h_gen_diEl_pt   = new TH1D("h_gen_diEl_pt","",10000,0,10000);
	TH1D*h_gen_diEl_rap  = new TH1D("h_gen_diEl_rap","",200,-10,10);

	Long64_t nEntries = chain->GetEntries();
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

	if(isMC){
		chain->SetBranchAddress("_prefiringweight", &_prefiringweight,
					&b__prefiringweight);
		chain->SetBranchAddress("GENnPair", &GENnPair, &b_GENnPair);
		chain->SetBranchAddress("GENLepton_Px",&GENLepton_Px,&b_GENLepton_Px);
		chain->SetBranchAddress("GENLepton_Py",&GENLepton_Py,&b_GENLepton_Py);
		chain->SetBranchAddress("GENLepton_Pz",&GENLepton_Pz,&b_GENLepton_Pz);
		chain->SetBranchAddress("GENLepton_E",&GENLepton_E,&b_GENLepton_E);
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

	// Find the gen weight sum
	double sumGenWeight = 0.0;
	double genWeight;
	if(isMC){
		Long64_t localEntry;
		cout << "Getting gen weights" << endl;
		for(Long64_t iGen=0;iGen<nEntries;iGen++){
			Counter(iGen,nEntries);
			localEntry = chain->LoadTree(iGen);
			b_GENEvt_weight->GetEntry(localEntry);
			genWeight = GENEvt_weight/fabs(GENEvt_weight);
			sumGenWeight += genWeight;
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
							       ptRecoLead,ptRecoSub);
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
		double diPtHard		= -1000;
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
		if(hardLep && passHLT){
			ptHardLead  = GENLepton_pT[idxHardLead];
			ptHardSub   = GENLepton_pT[idxHardSub];
			etaHardLead = GENLepton_eta[idxHardLead];
			etaHardSub  = GENLepton_eta[idxHardSub];
			phiHardLead = GENLepton_phi[idxHardLead];
			phiHardSub  = GENLepton_phi[idxHardSub];
		}
		vector<double>hardVariables = GetVariables(etaHardLead,etaHardSub,ptHardLead,ptHardSub,
                                             phiHardLead,phiHardSub);
		invMassHard     = hardVariables.at(0);
                rapidityHard    = hardVariables.at(1);
		diPtHard	= hardVariables.at(4);

		double genWeight = 1.0;
		if(isMC) genWeight = (GENEvt_weight/fabs(GENEvt_weight))/sumGenWeight;

	h_gen_el_pt	->Fill(ptHardLead,genWeight);
	h_gen_el_pt	->Fill(ptHardSub,genWeight);
	h_gen_el_eta	->Fill(etaHardLead,genWeight);
	h_gen_el_eta	->Fill(etaHardSub,genWeight);
	h_gen_el_phi	->Fill(phiHardLead,genWeight);
	h_gen_el_phi	->Fill(phiHardSub,genWeight);
	h_gen_diEl_mass	->Fill(invMassHard,genWeight);
	h_gen_diEl_pt	->Fill(diPtHard,genWeight);
	h_gen_diEl_rap	->Fill(rapidityHard,genWeight);

	}// end loop over entries

	// Save results to output file
	TString saveName = "output_data/crossChecks.root";
	TFile*file;
	file = new TFile(saveName,"recreate");
	h_gen_el_pt	->Write();
	h_gen_el_pt	->Write();
	h_gen_el_eta	->Write();
	h_gen_el_eta	->Write();
	h_gen_el_phi	->Write();
	h_gen_el_phi	->Write();
	h_gen_diEl_mass	->Write();
	h_gen_diEl_pt	->Write();
	h_gen_diEl_rap	->Write();
	file->Close();
}


bool PassDileptonSelection(double eta1,double eta2,double pt1,double pt2)
{
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
	double diPt = (v1+v2).Pt();

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
		ptSub,		// 3
		diPt
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
			if(GENLepton_fromHardProcessFinalState[iLep]==1 && 
			   GENLepton_fromHardProcessFinalState[jLep]==1){
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

void Counter(Long64_t event,Long64_t total)
{
	int P = 100*(event)/(total);
	if(event%(total/100)==0)
		cout << P << "%" << endl;
	 return;
}
