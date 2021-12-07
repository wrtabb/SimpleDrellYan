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

// Functions
bool PassDileptonSelection(double eta1,double eta2,double pt1,double pt2);
vector<double> GetVariables(double eta1,double eta2,double pt1,double pt2,double phi1,
                            double phi2);
bool GetRecoLeptons(int &idxRecoLead, int &idxRecoSub);
bool GetHardLeptons(int &idxHard1,int &idxHard2);
void Counter(Long64_t event,Long64_t total);

TString base_directory = "/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/wtabb/DrellYan_13TeV_2016/v2p6/DYLL_M50toInf/base/";

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
	TH1D*h_gen_mu_pt     = new TH1D("h_gen_mu_pt","",10000,0,10000);
        TH1D*h_gen_mu_eta    = new TH1D("h_gen_mu_eta","",200,-10,10);
        TH1D*h_gen_mu_phi    = new TH1D("h_gen_mu_phi","",80,-4,4);
        TH1D*h_gen_diMu_mass = new TH1D("h_gen_diMu_mass","",10000,0,10000);
        TH1D*h_gen_diMu_pt   = new TH1D("h_gen_diMu_pt","",10000,0,10000);
        TH1D*h_gen_diMu_rap  = new TH1D("h_gen_diMu_rap","",200,-10,10);

        TH1D*h_gen_acc_mu_pt     = new TH1D("h_gen_acc_mu_pt","",10000,0,10000);
        TH1D*h_gen_acc_mu_eta    = new TH1D("h_gen_acc_mu_eta","",200,-10,10);
        TH1D*h_gen_acc_mu_phi    = new TH1D("h_gen_acc_mu_phi","",80,-4,4);
        TH1D*h_gen_acc_diMu_mass = new TH1D("h_gen_acc_diMu_mass","",10000,0,10000);
        TH1D*h_gen_acc_diMu_pt   = new TH1D("h_gen_acc_diMu_pt","",10000,0,10000);
        TH1D*h_gen_acc_diMu_rap  = new TH1D("h_gen_acc_diMu_rap","",200,-10,10);

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
	}// end if monte carlo

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

		//-----Get Hard Process Quantities-----//
		double invMassHard      = -1000;
                double rapidityHard     = -1000;
		double diPtHard		= -1000;
                double leadPtHard       = -1000;
                double subPtHard        = -1000;

                double ptHardLead  = -1000;
                double ptHardSub   = -1000;
                double etaHardLead = -1000;
                double etaHardSub  = -1000;
                double phiHardLead = -1000;
                double phiHardSub  = -1000;

                int idxHardLead = -1;
                int idxHardSub  = -1;

                bool hardLep = GetHardLeptons(idxHardLead,idxHardSub);
		if(!hardLep) continue;

		ptHardLead  = GENLepton_pT[idxHardLead];
		ptHardSub   = GENLepton_pT[idxHardSub];
		etaHardLead = GENLepton_eta[idxHardLead];
		etaHardSub  = GENLepton_eta[idxHardSub];
		phiHardLead = GENLepton_phi[idxHardLead];
		phiHardSub  = GENLepton_phi[idxHardSub];

		// Get Hard Variables
		vector<double> hardVariables;
                hardVariables = GetVariables(etaHardLead,etaHardSub,ptHardLead,ptHardSub,
                                             phiHardLead,phiHardSub);

		invMassHard     = hardVariables.at(0);
		rapidityHard    = hardVariables.at(1);
		diPtHard	= hardVariables.at(4);

		//-----Get Reconstructed Quantities-----//
		double invMassReco      = -1000;
                double rapidityReco     = -1000;
                double leadPtReco       = -1000;
                double subPtReco        = -1000;
		double diPtReco         = -1000;

		double ptRecoLead  = -1000;
		double ptRecoSub   = -1000;
		double etaRecoLead = -1000;
		double etaRecoSub  = -1000;
		double phiRecoLead = -1000;
		double phiRecoSub  = -1000;

		int idxRecoLead = -1;
		int idxRecoSub = -1;

		bool recoLep = GetRecoLeptons(idxRecoLead,idxRecoSub);
		ptRecoLead  = Muon_pT[idxRecoLead]; 
		ptRecoSub   = Muon_pT[idxRecoSub];
		etaRecoLead = Muon_eta[idxRecoLead];
		etaRecoSub  = Muon_eta[idxRecoSub];
		phiRecoLead = Muon_phi[idxRecoLead];
		phiRecoSub  = Muon_phi[idxRecoSub];

		vector<double> recoVariables;
                recoVariables = GetVariables(etaRecoLead,etaRecoSub,ptRecoLead,ptRecoSub,
                                             phiRecoLead,phiRecoSub);
		invMassReco     = recoVariables.at(0);
		rapidityReco    = recoVariables.at(1);
		leadPtReco      = recoVariables.at(2);
		subPtReco       = recoVariables.at(3);
		diPtReco	= recoVariables.at(4);
	
		double genWeight = 1.0;
                if(isMC){
                        if(GENEvt_weight<0) genWeight = -1.0;
                }// end isMC for gen weight sum calculation

		h_gen_mu_pt     ->Fill(ptHardLead,genWeight);
                h_gen_mu_pt     ->Fill(ptHardSub,genWeight);
                h_gen_mu_eta    ->Fill(etaHardLead,genWeight);
                h_gen_mu_eta    ->Fill(etaHardSub,genWeight);
                h_gen_mu_phi    ->Fill(phiHardLead,genWeight);
                h_gen_mu_phi    ->Fill(phiHardSub,genWeight);
                h_gen_diMu_mass ->Fill(invMassHard,genWeight);
                h_gen_diMu_pt   ->Fill(diPtHard,genWeight);
                h_gen_diMu_rap  ->Fill(rapidityHard,genWeight);

		bool passGenAcc = PassDileptonSelection(etaHardLead,etaHardSub,
			   			        ptHardLead,ptHardSub);

		if(passGenAcc){
                        h_gen_acc_mu_pt         ->Fill(ptHardLead,genWeight);
                        h_gen_acc_mu_pt         ->Fill(ptHardSub,genWeight);
                        h_gen_acc_mu_eta        ->Fill(etaHardLead,genWeight);
                        h_gen_acc_mu_eta        ->Fill(etaHardSub,genWeight);
                        h_gen_acc_mu_phi        ->Fill(phiHardLead,genWeight);
                        h_gen_acc_mu_phi        ->Fill(phiHardSub,genWeight);
                        h_gen_acc_diMu_mass     ->Fill(invMassHard,genWeight);
                        h_gen_acc_diMu_pt       ->Fill(diPtHard,genWeight);
                        h_gen_acc_diMu_rap      ->Fill(rapidityHard,genWeight);
                }

		if(!recoLep) continue;


		bool passRecoAcc = PassDileptonSelection(etaRecoLead,etaRecoSub,
							 ptRecoLead,ptRecoSub);
		if(passRecoAcc){
			h_reco_mu_pt            ->Fill(ptRecoLead,genWeight);
			h_reco_mu_pt            ->Fill(ptRecoSub,genWeight);
			h_reco_mu_eta           ->Fill(etaRecoLead,genWeight);
			h_reco_mu_eta           ->Fill(etaRecoSub,genWeight);
			h_reco_mu_phi           ->Fill(phiRecoLead,genWeight);
			h_reco_mu_phi           ->Fill(phiRecoSub,genWeight);
			h_reco_mu_lead_pt       ->Fill(ptRecoLead,genWeight);
			h_reco_mu_lead_eta      ->Fill(etaRecoLead,genWeight);
			h_reco_mu_lead_phi      ->Fill(phiRecoLead,genWeight);
			h_reco_mu_sub_pt        ->Fill(ptRecoSub,genWeight);
			h_reco_mu_sub_eta       ->Fill(etaRecoSub,genWeight);
			h_reco_mu_sub_phi       ->Fill(phiRecoSub,genWeight);
			h_reco_diMu_mass        ->Fill(invMassReco,genWeight);
			h_reco_diMu_pt          ->Fill(diPtReco,genWeight);
			h_reco_diMu_rap         ->Fill(rapidityReco,genWeight);
		}

	}// end loop over entries


	// Save results to output file
	TString saveName = "output_data/";
	saveName += fileName;
	saveName += ".root";
	TFile*file;
	file = new TFile(saveName,"recreate");
	h_gen_mu_pt             ->Write();
        h_gen_mu_pt             ->Write();
        h_gen_mu_eta            ->Write();
        h_gen_mu_eta            ->Write();
        h_gen_mu_phi            ->Write();
        h_gen_mu_phi            ->Write();
        h_gen_diMu_mass         ->Write();
        h_gen_diMu_pt           ->Write();
        h_gen_diMu_rap          ->Write();

        h_gen_acc_mu_pt         ->Write();
        h_gen_acc_mu_pt         ->Write();
        h_gen_acc_mu_eta        ->Write();
        h_gen_acc_mu_eta        ->Write();
        h_gen_acc_mu_phi        ->Write();
        h_gen_acc_mu_phi        ->Write();
        h_gen_acc_diMu_mass     ->Write();
        h_gen_acc_diMu_pt       ->Write();
        h_gen_acc_diMu_rap      ->Write();

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

vector<double> GetVariables(double eta1,double eta2,double pt1,double pt2,double phi1,
                            double phi2)
{
        TLorentzVector v1;
        TLorentzVector v2;

        v1.SetPtEtaPhiM(pt1,eta1,phi1,muMass);
        v2.SetPtEtaPhiM(pt2,eta2,phi2,muMass);

	// Dimuon invariant mass
        double invMass = (v1+v2).M();

	// Dimuon rapidity
        double rapidity = (v1+v2).Rapidity();

	// dimuon pT
	double diPt = (v1+v2).Pt();

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
                invMass,        // 0
                rapidity,       // 1
                ptLead,         // 2
                ptSub,          // 3
		diPt		// 4
        };

        return variableReturn;
}// end GetVariables()

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

bool GetHardLeptons(int &idxHardLead,int &idxHardSub)
{
	int nDileptons = 0;
        for(int iLep=0;iLep<GENnPair;iLep++){
                for(int jLep=iLep+1;jLep<GENnPair;jLep++){
                        if(!(abs(GENLepton_ID[iLep])==13 && abs(GENLepton_ID[jLep])==13))
                                continue;
                        if(GENLepton_fromHardProcessFinalState[iLep]==1 &&
                           GENLepton_fromHardProcessFinalState[jLep]==1){
                                if(GENLepton_pT[iLep] > GENLepton_pT[jLep]){
                                        idxHardLead = iLep;
                                        idxHardSub = jLep;
                                }// end if iLep is leading electron
                                else{
                                        idxHardLead = jLep;
                                        idxHardSub = iLep;
                                }// end if jLep is leading electron
				nDileptons++;
                        }// end if hard process
                }//end inner loop over gen leptons
        }//end outer loop over gen leptons

        if(nDileptons==1) return true;
        else return false;
}// end GetHardLeptons()

void Counter(Long64_t event,Long64_t total)
{
        int P = 100*(event)/(total);
        if(event%(total/100)==0) 
                cout << P << "%" << endl;
         return;
}

