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

enum LepVariables{
	PT_LEAD,
	PT_SUB,
	ETA_LEAD,
	ETA_SUB,
	PHI_LEAD,
	PHI_SUB,
	MASS,
	RAPIDITY,
	DI_PT
};

// Functions
bool PassDileptonSelection(double eta1,double eta2,double pt1,double pt2);
double GetVar(LepVariables var,int idxLead,int idxSub);
bool GetRecoLeptons(int &idxRecoLead, int &idxRecoSub);
void Counter(Long64_t event,Long64_t total);

TString base_directory = "/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/wtabb/DrellYan_13TeV_2016/v2p6/skims/skims_MuMu/";

TString treeName = "recoTree/DYTree";
int dataLuminosity = 35867;
const TString muonTrigger1 = "HLT_IsoMu24_v*";
const TString muonTrigger2 = "HLT_IsoTkMu24_v*";
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
int nMuon;
double Muon_pT[MPSIZE];
double Muon_eta[MPSIZE];
double Muon_phi[MPSIZE];
bool Muon_passTightID[MPSIZE];
double Muon_PfChargedHadronIsoR04[MPSIZE];
double Muon_PfNeutralHadronIsoR04[MPSIZE];
double Muon_PfGammaIsoR04[MPSIZE];
double Muon_PFSumPUIsoR04[MPSIZE];
double Muon_trkiso[MPSIZE];
double muMass = 0.1056583715;
int HLT_ntrig;
int HLT_trigType[MPSIZE];
int HLT_trigFired[MPSIZE];
std::vector<std::string> HLT_trigName;
std::vector<std::string> *pHLT_trigName = &HLT_trigName;

TBranch*b_HLT_ntrig;
TBranch*b_HLT_trigType;
TBranch*b_HLT_trigFired;
TBranch*b_nMuon;
TBranch*b_Muon_pT;
TBranch*b_Muon_eta;
TBranch*b_Muon_phi;
TBranch*b_Muon_Inner_pT;
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

        TH1D*h_reco_mu_pt       = new TH1D("h_reco_mu_pt","",10000,0,10000);
        TH1D*h_reco_mu_eta      = new TH1D("h_reco_mu_eta","",60,-3,3);
        TH1D*h_reco_mu_phi      = new TH1D("h_reco_mu_phi","",80,-4,4);
        TH1D*h_reco_mu_lead_pt  = new TH1D("h_reco_mu_lead_pt","",10000,0,10000);
        TH1D*h_reco_mu_lead_eta = new TH1D("h_reco_mu_lead_eta","",60,-3,3);
        TH1D*h_reco_mu_lead_phi = new TH1D("h_reco_mu_lead_phi","",80,-4,4);
        TH1D*h_reco_mu_sub_pt   = new TH1D("h_reco_mu_sub_pt","",10000,0,10000);
        TH1D*h_reco_mu_sub_eta  = new TH1D("h_reco_mu_sub_eta","",60,-3,3);
        TH1D*h_reco_mu_sub_phi  = new TH1D("h_reco_mu_sub_phi","",80,-4,4);
        TH1D*h_reco_diMu_mass   = new TH1D("h_reco_diMu_mass","",10000,0,10000);
        TH1D*h_reco_diMu_pt     = new TH1D("h_reco_diMu_pt","",10000,0,10000);
        TH1D*h_reco_diMu_rap    = new TH1D("h_reco_diMu_rap","",60,-3,3);

	Long64_t nEntries = chain->GetEntries();
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
			if(((trigName.CompareTo(muonTrigger1)==0)||
			    (trigName.CompareTo(muonTrigger2)==0)) && 
			     HLT_trigFired[iHLT]==1){
				passHLT = true;
				break;
			} // end if trigName
		}// end loop over triggers

		if(!passHLT) continue;

		//-----Get Reconstructed Quantities-----//
		double invMassReco      = -1000;
                double rapidityReco     = -1000;
		double diPtReco         = -1000;
		double ptRecoLead  	= -1000;
		double ptRecoSub   	= -1000;
		double etaRecoLead 	= -1000;
		double etaRecoSub  	= -1000;
		double phiRecoLead 	= -1000;
		double phiRecoSub  	= -1000;

		int idxRecoLead = -1;
		int idxRecoSub  = -1;

		bool recoLep = GetRecoLeptons(idxRecoLead,idxRecoSub);

		ptRecoLead  	= GetVar(PT_LEAD,idxRecoLead,idxRecoSub);
		ptRecoSub   	= GetVar(PT_SUB,idxRecoLead,idxRecoSub);
		etaRecoLead 	= GetVar(ETA_LEAD,idxRecoLead,idxRecoSub);
		etaRecoSub  	= GetVar(ETA_SUB,idxRecoLead,idxRecoSub);
		phiRecoLead 	= GetVar(PHI_LEAD,idxRecoLead,idxRecoSub);
		phiRecoSub  	= GetVar(PHI_SUB,idxRecoLead,idxRecoSub);
		invMassReco 	= GetVar(MASS,idxRecoLead,idxRecoSub);
		rapidityReco	= GetVar(RAPIDITY,idxRecoLead,idxRecoSub);
		diPtReco	= GetVar(DI_PT,idxRecoLead,idxRecoSub);	

		if(!recoLep) continue;

		bool passRecoAcc = PassDileptonSelection(etaRecoLead,etaRecoSub,
							 ptRecoLead,ptRecoSub);
		if(!passRecoAcc) continue;

		h_reco_mu_pt            ->Fill(ptRecoLead);
		h_reco_mu_pt            ->Fill(ptRecoSub);
		h_reco_mu_eta           ->Fill(etaRecoLead);
		h_reco_mu_eta           ->Fill(etaRecoSub);
		h_reco_mu_phi           ->Fill(phiRecoLead);
		h_reco_mu_phi           ->Fill(phiRecoSub);
		h_reco_mu_lead_pt       ->Fill(ptRecoLead);
		h_reco_mu_lead_eta      ->Fill(etaRecoLead);
		h_reco_mu_lead_phi      ->Fill(phiRecoLead);
		h_reco_mu_sub_pt        ->Fill(ptRecoSub);
		h_reco_mu_sub_eta       ->Fill(etaRecoSub);
		h_reco_mu_sub_phi       ->Fill(phiRecoSub);
		h_reco_diMu_mass        ->Fill(invMassReco);
		h_reco_diMu_pt          ->Fill(diPtReco);
		h_reco_diMu_rap         ->Fill(rapidityReco);

	}// end loop over entries


	// Save results to output file
	TString saveName = "output_data/";
	saveName += fileName;
	saveName += ".root";
	TFile*file;
	file = new TFile(saveName,"recreate");
        h_reco_mu_pt            ->Write();
        h_reco_mu_pt            ->Write();
        h_reco_mu_eta           ->Write();
        h_reco_mu_eta           ->Write();
        h_reco_mu_phi           ->Write();
        h_reco_mu_phi           ->Write();
        h_reco_mu_lead_pt       ->Write();
        h_reco_mu_lead_eta      ->Write();
        h_reco_mu_lead_phi      ->Write();
        h_reco_mu_sub_pt        ->Write();
        h_reco_mu_sub_eta       ->Write();
        h_reco_mu_sub_phi       ->Write();

        h_reco_diMu_mass        ->Write();
        h_reco_diMu_pt          ->Write();
        h_reco_diMu_rap         ->Write();
	file->Close();
}// end analyze

bool PassDileptonSelection(double eta1,double eta2,double pt1,double pt2)
{
        if(abs(eta1)>etaHigh || abs(eta2)>etaHigh) return false;
        if(pt1>pt2 && (pt1<ptHigh || pt2<ptLow)) return false;
        if(pt2>pt1 && (pt2<ptHigh || pt1<ptLow)) return false;

        return true;
}// end PassDileptonSelection()

bool GetRecoLeptons(int &idxRecoLead, int &idxRecoSub)
{
	// isolation-related variables
	double chargedIso;
	double neutralIso;
	double gammaIso;
	double sumPUPt;
	double iso_dBeta;

	// kinematic variables
	double pt1,pt2;  
	double eta1,eta2;
	double phi1,phi2;

	int nDileptons = 0;

	for(int iMu=0;iMu<nMuon;iMu++){
		if(!Muon_passTightID[iMu]) continue;
		chargedIso = Muon_PfChargedHadronIsoR04[iMu];
		neutralIso = Muon_PfNeutralHadronIsoR04[iMu];
		gammaIso   = Muon_PfGammaIsoR04[iMu];
		sumPUPt    = Muon_PFSumPUIsoR04[iMu];
		pt1        = Muon_pT[iMu];
		eta1       = Muon_eta[iMu];
		phi1       = Muon_phi[iMu];
		iso_dBeta = 
			(chargedIso+max(0.0,neutralIso+gammaIso-0.5*sumPUPt))/pt1;
		if(iso_dBeta > 0.15) continue;

		for(int jMu=iMu+1;jMu<nMuon;jMu++){
			if(!Muon_passTightID[jMu]) continue;
			chargedIso = Muon_PfChargedHadronIsoR04[jMu];
			neutralIso = Muon_PfNeutralHadronIsoR04[jMu];
			gammaIso   = Muon_PfGammaIsoR04[jMu];
			sumPUPt    = Muon_PFSumPUIsoR04[jMu];
			pt2        = Muon_pT[jMu];
			eta2       = Muon_eta[jMu];
			phi2       = Muon_phi[jMu];
			iso_dBeta = 
				(chargedIso+max(0.0,neutralIso+gammaIso-0.5*sumPUPt))/pt2;
			if(iso_dBeta > 0.15) continue;

			if(pt1 > pt2){
				idxRecoLead = iMu;
				idxRecoSub  = jMu;
			}
			else{
				idxRecoLead = jMu;
				idxRecoSub  = iMu;
			} 
			nDileptons++;
		}//end inner muon loop
	}// end outer muon loop 

	if(nDileptons==1) return true; 
	else return false;
}// end GetRecoLeptons()

void Counter(Long64_t event,Long64_t total)
{
        int P = 100*(event)/(total);
        if(event%(total/100)==0) 
                cout << P << "%" << endl;
         return;
}


double GetVar(LepVariables var,int idxLead,int idxSub)
{
	double ptLead	= Muon_pT[idxLead];
	double ptSub	= Muon_pT[idxSub];
	double etaLead	= Muon_eta[idxLead];
	double etaSub	= Muon_eta[idxSub];
	double phiLead	= Muon_phi[idxLead];
	double phiSub	= Muon_phi[idxSub];

	TLorentzVector v1;
	TLorentzVector v2;
	v1.SetPtEtaPhiM(ptLead,etaLead,phiLead,muMass);
	v2.SetPtEtaPhiM(ptSub,etaSub,phiSub,muMass);

	double invMass	= (v1+v2).M();
	double rapidity	= (v1+v2).Rapidity();
	double diPt	= (v1+v2).Pt();

	if(var==PT_LEAD)	return ptLead;
	else if(var==PT_SUB)	return ptSub;
	else if(var==ETA_LEAD)	return etaLead;
	else if(var==ETA_SUB)	return etaSub;
	else if(var==PHI_LEAD)	return phiLead;
	else if(var==PHI_SUB)	return phiSub;
	else if(var==MASS)	return invMass;
	else if(var==RAPIDITY)	return rapidity;
	else if(var==DI_PT)	return diPt;

	else{
		cout << "LepVariables enum not properly defined for GetVar()" << endl;
		return -10000;
	}
}
