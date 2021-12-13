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
bool GetRecoLeptons(int &idxRecoLead, int &idxRecoSub);
bool GetHardLeptons(int &idxHard1,int &idxHard2);
void Counter(Long64_t event,Long64_t total);

//TString base_directory = "root://xrootd-local.unl.edu///store/user/wtabb/DrellYan_13TeV_2016/v2p6/DYLL_M50toInf/base/";
TString base_directory = "/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/wtabb/DrellYan_13TeV_2016/v2p6/DYLL_M50toInf/base/";

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
int GENLepton_ID[MPSIZE];
int GENLepton_fromHardProcessFinalState[MPSIZE];
int Nelectrons;
double Electron_pT[MPSIZE];
double Electron_eta[MPSIZE];
double Electron_phi[MPSIZE];
double Electron_etaSC[MPSIZE]; //no muon
double Electron_phiSC[MPSIZE]; //no muon
double Electron_etSC[MPSIZE]; //no muon
bool Electron_passMediumID[MPSIZE];
double eMass = 0.000510998;
int HLT_ntrig;
int HLT_trigType[MPSIZE];
int HLT_trigFired[MPSIZE];
std::vector<std::string> HLT_trigName;
std::vector<std::string> *pHLT_trigName = &HLT_trigName;

TBranch*b_Nelectrons;
TBranch*b_HLT_ntrig;
TBranch*b_HLT_trigType;
TBranch*b_HLT_trigFired;
TBranch*b_Electron_pT;
TBranch*b_Electron_eta;
TBranch*b_Electron_phi;
TBranch*b_Electron_etaSC;
TBranch*b_Electron_phiSC;
TBranch*b_Electron_etSC;
TBranch*b_Electron_passMediumID;
TBranch*b_GENnPair;
TBranch*b_GENLepton_phi;
TBranch*b_GENLepton_eta;
TBranch*b_GENLepton_pT;
TBranch*b_GENLepton_ID;
TBranch*b_GENLepton_fromHardProcessFinalState;
TBranch*b_GENEvt_weight;

void teamCrossChecks(TString fileName)
{
	TH1::SetDefaultSumw2();

	TString loadName = base_directory;
	loadName += fileName;
	loadName += ".root";
	TFile*loadFile = new TFile(loadName);
	TTree*chain=(TTree*)loadFile->Get(treeName);

	//-----Define histograms-----/
	// gen electrons
	TH1D*h_gen_el_pt     = new TH1D("h_gen_el_pt","",10000,0,10000);
	TH1D*h_gen_el_eta    = new TH1D("h_gen_el_eta","",200,-10,10);
	TH1D*h_gen_el_phi    = new TH1D("h_gen_el_phi","",80,-4,4);
	TH1D*h_gen_diEl_mass = new TH1D("h_gen_diEl_mass","",10000,0,10000);
	TH1D*h_gen_diEl_pt   = new TH1D("h_gen_diEl_pt","",10000,0,10000);
	TH1D*h_gen_diEl_rap  = new TH1D("h_gen_diEl_rap","",200,-10,10);

	// gen electrons passing acceptance
	TH1D*h_gen_acc_el_pt     = new TH1D("h_gen_acc_el_pt","",10000,0,10000);
	TH1D*h_gen_acc_el_eta    = new TH1D("h_gen_acc_el_eta","",200,-10,10);
	TH1D*h_gen_acc_el_phi    = new TH1D("h_gen_acc_el_phi","",80,-4,4);
	TH1D*h_gen_acc_diEl_mass = new TH1D("h_gen_acc_diEl_mass","",10000,0,10000);
	TH1D*h_gen_acc_diEl_pt   = new TH1D("h_gen_acc_diEl_pt","",10000,0,10000);
	TH1D*h_gen_acc_diEl_rap  = new TH1D("h_gen_acc_diEl_rap","",200,-10,10);

	// reco electrons
	TH1D*h_reco_el_pt 	= new TH1D("h_reco_el_pt","",10000,0,10000);
	TH1D*h_reco_el_eta 	= new TH1D("h_reco_el_eta","",60,-3,3);
	TH1D*h_reco_el_phi 	= new TH1D("h_reco_el_phi","",80,-4,4);
	TH1D*h_reco_el_lead_pt 	= new TH1D("h_reco_el_lead_pt","",10000,0,10000);
	TH1D*h_reco_el_lead_eta = new TH1D("h_reco_el_lead_eta","",60,-3,3);
	TH1D*h_reco_el_lead_phi = new TH1D("h_reco_el_lead_phi","",80,-4,4);
	TH1D*h_reco_el_sub_pt 	= new TH1D("h_reco_el_sub_pt","",10000,0,10000);
	TH1D*h_reco_el_sub_eta 	= new TH1D("h_reco_el_sub_eta","",60,-3,3);
	TH1D*h_reco_el_sub_phi 	= new TH1D("h_reco_el_sub_phi","",80,-4,4);
	TH1D*h_reco_diEl_mass 	= new TH1D("h_reco_diEl_mass","",10000,0,10000);
	TH1D*h_reco_diEl_pt 	= new TH1D("h_reco_diEl_pt","",10000,0,10000);
	TH1D*h_reco_diEl_rap 	= new TH1D("h_reco_diEl_rap","",60,-3,3);

	Long64_t nEntries = chain->GetEntries();
	cout << nEntries << " entries loaded. " << endl;

	bool isMC = true;
	TBranch*testBranch = 
		(TBranch*)chain->GetListOfBranches()->FindObject("GENEvt_weight");
	if(!testBranch) isMC = false;

	// Define branches
	chain->SetBranchAddress("Nelectrons",&Nelectrons,&b_Nelectrons);
	chain->SetBranchAddress("Electron_pT",&Electron_pT,&b_Electron_pT);
	chain->SetBranchAddress("Electron_eta",&Electron_eta,&b_Electron_eta);
	chain->SetBranchAddress("Electron_etaSC",&Electron_etaSC,&b_Electron_etaSC);
	chain->SetBranchAddress("Electron_phi",&Electron_phi,&b_Electron_phi);
	chain->SetBranchAddress("Electron_passMediumID",&Electron_passMediumID,
				 &b_Electron_passMediumID);
	chain->SetBranchAddress("HLT_ntrig",&HLT_ntrig,&b_HLT_ntrig);
	chain->SetBranchAddress("HLT_trigType",&HLT_trigType,&b_HLT_trigType);
	chain->SetBranchAddress("HLT_trigFired",&HLT_trigFired,&b_HLT_trigFired);
	chain->SetBranchAddress("HLT_trigName",&pHLT_trigName);

	if(isMC){
		chain->SetBranchAddress("GENnPair", &GENnPair, &b_GENnPair);
		chain->SetBranchAddress("GENLepton_eta", &GENLepton_eta, 
					&b_GENLepton_eta);
		chain->SetBranchAddress("GENLepton_phi",&GENLepton_phi, 
					&b_GENLepton_phi);
		chain->SetBranchAddress("GENLepton_pT",&GENLepton_pT, 
					&b_GENLepton_pT);
		chain->SetBranchAddress("GENLepton_ID",&GENLepton_ID, 
					&b_GENLepton_ID);
		chain->SetBranchAddress("GENLepton_fromHardProcessFinalState",
					&GENLepton_fromHardProcessFinalState,
					&b_GENLepton_fromHardProcessFinalState);
		chain->SetBranchAddress("GENEvt_weight",&GENEvt_weight,
                                        &b_GENEvt_weight);
	}// end if monte carlo

	// Find the gen weight sum

	// Loop over events
	for(Long64_t iEntry=0;iEntry<nEntries;iEntry++){
		chain->GetEntry(iEntry);
//		Counter(iEntry,nEntries);

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

		//-----Get Hard Process Quantities-----//
		double invMassHard	= -1000;
		double rapidityHard	= -1000;
		double diPtHard		= -1000;
		double leadPtHard	= -1000;
		double subPtHard	= -1000;

		double ptHard1  = -1000;
		double ptHard2   = -1000;
		double etaHard1 = -1000;
		double etaHard2  = -1000;
		double phiHard1 = -1000;
		double phiHard2  = -1000;

		int idxHard1 = -1;
		int idxHard2  = -1;

		bool hardLep = GetHardLeptons(idxHard1,idxHard2);
		if(!hardLep) continue;

		ptHard1  = GENLepton_pT[idxHard1];
		ptHard2   = GENLepton_pT[idxHard2];
		etaHard1 = GENLepton_eta[idxHard1];
		etaHard2  = GENLepton_eta[idxHard2];
		phiHard1 = GENLepton_phi[idxHard1];
		phiHard2  = GENLepton_phi[idxHard2];

		TLorentzVector v1;
		TLorentzVector v2;

		v1.SetPtEtaPhiM(ptHard1,etaHard1,phiHard1,eMass);
		v2.SetPtEtaPhiM(ptHard2,etaHard2,phiHard2,eMass);

		invMassHard     = (v1+v2).M();
                rapidityHard    = (v1+v2).Rapidity();
		diPtHard	= (v1+v2).Pt();

		//-----Get Reconstructed Quantities-----//
		double invMassReco	= -1000;
		double rapidityReco	= -1000;
		double diPtReco		= -1000;

		double ptRecoLead  = -1000;
		double ptRecoSub   = -1000;
		double etaSCRecoLead = -1000;
		double etaSCRecoSub  = -1000;
		double etaRecoLead = -1000;
		double etaRecoSub  = -1000;
		double phiRecoLead = -1000;
		double phiRecoSub  = -1000;

		int idxRecoLead = -1;
		int idxRecoSub  = -1;

		bool recoLep 	= GetRecoLeptons(idxRecoLead,idxRecoSub);
		ptRecoLead  	= Electron_pT[idxRecoLead];
		ptRecoSub   	= Electron_pT[idxRecoSub];
		etaSCRecoLead 	= Electron_etaSC[idxRecoLead];
		etaSCRecoSub  	= Electron_etaSC[idxRecoSub];
		etaRecoLead 	= Electron_eta[idxRecoLead];
		etaRecoSub  	= Electron_eta[idxRecoSub];
		phiRecoLead 	= Electron_phi[idxRecoLead];
		phiRecoSub  	= Electron_phi[idxRecoSub];

		v1.SetPtEtaPhiM(ptRecoLead,etaRecoLead,phiRecoLead,eMass);
		v2.SetPtEtaPhiM(ptRecoSub,etaRecoSub,phiRecoSub,eMass);

		invMassReco     = (v1+v2).M();
                rapidityReco    = (v1+v2).Rapidity();
		diPtReco	= (v1+v2).Pt();

		double genWeight = 1.0;
		if(isMC){
			if(GENEvt_weight<0) genWeight = -1.0;
		}// end isMC for gen weight sum calculation

		// fill gen histograms
		h_gen_el_pt	->Fill(ptHard1,genWeight);
		h_gen_el_pt	->Fill(ptHard2,genWeight);
		h_gen_el_eta	->Fill(etaHard1,genWeight);
		h_gen_el_eta	->Fill(etaHard2,genWeight);
		h_gen_el_phi	->Fill(phiHard1,genWeight);
		h_gen_el_phi	->Fill(phiHard2,genWeight);
		h_gen_diEl_mass	->Fill(invMassHard,genWeight);
		h_gen_diEl_pt	->Fill(diPtHard,genWeight);
		h_gen_diEl_rap	->Fill(rapidityHard,genWeight);

		bool passGenAcc = PassDileptonSelection(etaHard1,etaHard2,ptHard1,ptHard2);
		if(passGenAcc && passHLT){
			// fill gen histograms for leptons within acceptance
			h_gen_acc_el_pt		->Fill(ptHard1,genWeight);
			h_gen_acc_el_pt		->Fill(ptHard2,genWeight);
			h_gen_acc_el_eta	->Fill(etaHard1,genWeight);
			h_gen_acc_el_eta	->Fill(etaHard2,genWeight);
			h_gen_acc_el_phi	->Fill(phiHard1,genWeight);
			h_gen_acc_el_phi	->Fill(phiHard2,genWeight);
			h_gen_acc_diEl_mass	->Fill(invMassHard,genWeight);
			h_gen_acc_diEl_pt	->Fill(diPtHard,genWeight);
			h_gen_acc_diEl_rap	->Fill(rapidityHard,genWeight);
		}

		if(!recoLep) continue;

		bool passRecoAcc = PassDileptonSelection(etaSCRecoLead,etaSCRecoSub,
				   		         ptRecoLead,ptRecoSub);
		if(passRecoAcc && passHLT){
			// fill reco histograms
			h_reco_diEl_mass	->Fill(invMassReco,genWeight);
			h_reco_diEl_pt		->Fill(diPtReco,genWeight);
			h_reco_diEl_rap		->Fill(rapidityReco,genWeight);
			h_reco_el_pt		->Fill(ptRecoLead,genWeight);
			h_reco_el_pt		->Fill(ptRecoSub,genWeight);
			h_reco_el_eta		->Fill(etaRecoLead,genWeight);
			h_reco_el_eta		->Fill(etaRecoSub,genWeight);
			h_reco_el_phi		->Fill(phiRecoLead,genWeight);
			h_reco_el_phi		->Fill(phiRecoSub,genWeight);
			h_reco_el_lead_pt	->Fill(ptRecoLead,genWeight);
			h_reco_el_lead_eta	->Fill(etaRecoLead,genWeight);
			h_reco_el_lead_phi	->Fill(phiRecoLead,genWeight);
			h_reco_el_sub_pt	->Fill(ptRecoSub,genWeight);
			h_reco_el_sub_eta	->Fill(etaRecoSub,genWeight);
			h_reco_el_sub_phi	->Fill(phiRecoSub,genWeight);
		}
	}// end loop over entries

	// Save results to output file
	TString saveName = "output_data/";
	saveName += fileName;
	saveName += ".root";
	TFile*file;
	file = new TFile(saveName,"recreate");
	h_gen_el_pt		->Write();
	h_gen_el_pt		->Write();
	h_gen_el_eta		->Write();
	h_gen_el_eta		->Write();
	h_gen_el_phi		->Write();
	h_gen_el_phi		->Write();
	h_gen_diEl_mass		->Write();
	h_gen_diEl_pt		->Write();
	h_gen_diEl_rap		->Write();

	h_gen_acc_el_pt		->Write();
	h_gen_acc_el_pt		->Write();
	h_gen_acc_el_eta	->Write();
	h_gen_acc_el_eta	->Write();
	h_gen_acc_el_phi	->Write();
	h_gen_acc_el_phi	->Write();
	h_gen_acc_diEl_mass	->Write();
	h_gen_acc_diEl_pt	->Write();
	h_gen_acc_diEl_rap	->Write();

	h_reco_el_pt		->Write();
	h_reco_el_pt		->Write();
	h_reco_el_eta		->Write();
	h_reco_el_eta		->Write();
	h_reco_el_phi		->Write();
	h_reco_el_phi		->Write();
	h_reco_el_lead_pt	->Write();
	h_reco_el_lead_eta	->Write();
	h_reco_el_lead_phi	->Write();
	h_reco_el_sub_pt	->Write();
	h_reco_el_sub_eta	->Write();
	h_reco_el_sub_phi	->Write();

	h_reco_diEl_mass	->Write();
	h_reco_diEl_pt		->Write();
	h_reco_diEl_rap		->Write();
	file->Close();
}


bool PassDileptonSelection(double eta1,double eta2,double pt1,double pt2)
{
	if(abs(eta1)>etaGapLow && abs(eta1)<etaGapHigh) return false;
	if(abs(eta2)>etaGapLow && abs(eta2)<etaGapHigh) return false;
	if(abs(eta1)>etaHigh || abs(eta2)>etaHigh) return false;
	if(pt1>pt2 && (pt1<ptHigh || pt2<ptLow)) return false;
	if(pt2>pt1 && (pt2<ptHigh || pt1<ptLow)) return false;

	return true;
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

bool GetHardLeptons(int &idxHard1,int &idxHard2)
{
	int nHardLeptons = 0;
	for(int iLep=0;iLep<GENnPair;iLep++){
		for(int jLep=iLep+1;jLep<GENnPair;jLep++){
			if(!(abs(GENLepton_ID[iLep])==11 && abs(GENLepton_ID[jLep])==11))
				continue;
			if(GENLepton_fromHardProcessFinalState[iLep]==1 && 
			   GENLepton_fromHardProcessFinalState[jLep]==1){
				idxHard1 = iLep;
				idxHard2 = jLep;
				nHardLeptons++;
			}// end if hard process
		}//end inner loop over gen leptons
	}//end outer loop over gen leptons

	if(nHardLeptons==1) return true;
	else return false;
}// end GetHardLeptons()

void Counter(Long64_t event,Long64_t total)
{
	int P = 100*(event)/(total);
	if(event%(total/100)==0)
		cout << P << "%" << endl;
	 return;
}
