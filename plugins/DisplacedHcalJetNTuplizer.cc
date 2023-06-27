// -*- C++ -*-
// Class:      DisplacedHcalJetNTuplizer
/*
	Description: Base class for miniAOD analysis with CRAB
*/
//   Author:  Cristián Peña and Si Xie.
//   Created:  Thu, 26 March 2019 15:00:06 GMT

#include "DisplacedHcalJetNTuplizer.h"
#include <fastjet/internal/base.hh>
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

// 
// ************************************************************************************
// Constructor & Destructor 
// ************************************************************************************
// 

// ------------------------------------------------------------------------------------
DisplacedHcalJetNTuplizer::DisplacedHcalJetNTuplizer(const edm::ParameterSet& iConfig):
	// Inputs from Config
	debug(iConfig.getParameter<bool>( "debug" )),
	isData_(iConfig.getParameter<bool>( "isData" )),
	isSignal_(iConfig.getParameter<bool>( "isSignal" )),
	// Trigger
	triggerBitsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"))),
	triggerObjectsToken_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
	triggerPrescalesToken_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("triggerPrescales"))),
	triggerPathNamesFile_(iConfig.getParameter<string>( "triggerPathNamesFile" )),
	// General 
	verticesToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
	primaryVertexAssociationToken_(consumes<edm::Association<vector<reco::Vertex> > >(edm::InputTag("primaryVertexAssociation","original"))),
	primaryVertexAssociationValueMapToken_(consumes<edm::ValueMap<int> >(edm::InputTag("primaryVertexAssociation","original"))),
	// Event-Level Info
	metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("met"))),
	bsTag_(iConfig.getUntrackedParameter<edm::InputTag>("offlineBeamSpot", edm::InputTag("offlineBeamSpot"))),
	bsToken_(consumes<reco::BeamSpot>(bsTag_)),
	// Physics Objects
	electronsToken_(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
	muonsToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
	tausToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
	photonsToken_(consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
	jetsToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("pfjetsAK4"))),
	calojetsToken_(consumes<reco::CaloJetCollection>(iConfig.getParameter<edm::InputTag>("calojetsAK4"))),
	LRJetsToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("pfjetsAK8"))),
	caloLRJetsToken_(consumes<reco::CaloJetCollection>(iConfig.getParameter<edm::InputTag>("calojetsAK8"))),
	l1jetsToken_(consumes<BXVector<l1t::Jet>>(iConfig.getParameter<edm::InputTag>("l1jets"))), // GK added for L1 jets access
	// Low-Level Objects
	generalTracksToken_(consumes<std::vector<reco::Track>>(edm::InputTag("generalTracks"))),
	secondaryVerticesToken_(consumes<vector<reco::VertexCompositePtrCandidate> >(iConfig.getParameter<edm::InputTag>("secondaryVertices"))),	
	PFCandsToken_(consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
	// RecHits
	hcalRecHitsHBHEToken_(consumes<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>>( iConfig.getParameter<edm::InputTag>("hbRecHits") )),
	// Other
	electron_cutbasedID_decisions_loose_Token_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electron_cutbasedID_decisions_loose"))),
	electron_cutbasedID_decisions_medium_Token_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electron_cutbasedID_decisions_medium"))),
	electron_cutbasedID_decisions_tight_Token_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electron_cutbasedID_decisions_tight"))),	
	photon_cutbasedID_decisions_loose_Token_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photon_cutbasedID_decisions_loose"))),
	photon_cutbasedID_decisions_medium_Token_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photon_cutbasedID_decisions_medium"))),
	photon_cutbasedID_decisions_tight_Token_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photon_cutbasedID_decisions_tight"))),
	// MC
	genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
	// Geometry
	caloGeometryToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()), // GK
	castorGeometryToken_(esConsumes<PCaloGeometry, PCastorRcd>()), // GK
	gemGeoToken_(esConsumes<RecoIdealGeometry, GEMRecoGeometryRcd>()), 
	magneticFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()), // GK
	propagatorToken_(esConsumes<Propagator, TrackingComponentsRecord>(edm::ESInputTag("", "PropagatorWithMaterial"))) // GK // thePropagatorToken(esConsumes<Propagator, TrackingComponentsRecord>(edm::ESInputTag("", thePropagatorName))),

	// MC

	{

	cout<<"Option Settings:"<<endl;
	cout<<" --> isData = "<<isData_<<endl;

	//declare the TFileService for output
	edm::Service<TFileService> fs;

	// ----- Set up output tree ----- // 

	if( debug ) cout<<"Setting up output tree..."<<endl;

	output_tree = fs->make<TTree>("Events", "");
	NEvents = fs->make<TH1F>( "NEvents", " ; ; NEvents;" , 1, -0.5, 0.5 );
	if (!isData_) {
		sumWeights = fs->make<TH1D>( "sumWeights", " ; ; sumWeights;", 1, -0.5, 0.5 );

		sumWeights->Sumw2();
	}
	else {
		sumWeights = 0;
	}

	// ----- Get Triggers ----- // 

	if( debug ) cout<<"Getting triggers..."<<endl;

	triggerPathNames.push_back("1 HLT_HT200_L1SingleLLPJet_DisplacedDijet35_Inclusive1PtrkShortSig5");
	triggerPathNames.push_back("2 HLT_HT200_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5");
	triggerPathNames.push_back("3 HLT_HT170_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack");
	triggerPathNames.push_back("4 HLT_HT200_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack");
	triggerPathNames.push_back("5 HLT_HT270_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack");
	triggerPathNames.push_back("6 HLT_HT200_L1SingleLLPJet_DisplacedDijet60_DisplacedTrack");
	triggerPathNames.push_back("7 HLT_HT320_L1SingleLLPJet_DisplacedDijet60_Inclusive");
	triggerPathNames.push_back("8 HLT_HT420_L1SingleLLPJet_DisplacedDijet60_Inclusive");
	triggerPathNames.push_back("9 HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay0p5nsTrackless");
	triggerPathNames.push_back("10 HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay1nsInclusive");
	triggerPathNames.push_back("11 HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay1nsTrackless");
	triggerPathNames.push_back("12 HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay2nsInclusive");

	if( debug ) cout<<"Defining trigger hists..."<<endl;

	const int NTriggers = triggerPathNames.size();
	NEvents_HLT = fs->make<TH1F>( "NEvents_HLT", " ; ; NEvents_HLT;" , NTriggers, 0, NTriggers );

	if( debug ) cout<<"Defining trigger hists..."<<endl;

	for( int i=0; i<NTriggers; i++ ){
		NEvents_HLT->GetXaxis()->SetBinLabel(i+1,  Form("%s", triggerPathNames.at(i).c_str() ) );
		triggerPathNamesIndices[triggerPathNames.at(i)] = i;
	}

	cout<<"DisplacedHcalJetNTuplizer DONE"<<endl;

}

// ------------------------------------------------------------------------------------
DisplacedHcalJetNTuplizer::~DisplacedHcalJetNTuplizer(){
};

// 
// ************************************************************************************
// Load Event 
// ************************************************************************************
// 

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::loadEvent(const edm::Event& iEvent){ 
	// Load the miniAOD objects and reset tree variables for each event

	if( debug ) cout<<"Running DisplacedHcalJetNTuplizer::loadEvent"<<endl;	

	iEvent.getByToken(triggerBitsToken_, triggerBits);
	iEvent.getByToken(triggerObjectsToken_, triggerObjects);
	if( isData_ ) iEvent.getByToken(triggerPrescalesToken_, triggerPrescales);

	iEvent.getByToken(verticesToken_, vertices);
	iEvent.getByToken(primaryVertexAssociationToken_,primaryVertexAssociation);
	iEvent.getByToken(primaryVertexAssociationValueMapToken_,primaryVertexAssociationValueMap);

	iEvent.getByToken(metToken_, met);

	// Physics objects
	iEvent.getByToken(electronsToken_, electrons);
	iEvent.getByToken(muonsToken_, muons);
	iEvent.getByToken(tausToken_, taus);
	iEvent.getByToken(photonsToken_, photons);
	iEvent.getByToken(jetsToken_, jets);
	iEvent.getByToken(calojetsToken_, calojets);
	iEvent.getByToken(LRJetsToken_, LRJets);
	iEvent.getByToken(caloLRJetsToken_, caloLRJets);
	iEvent.getByToken(l1jetsToken_, l1jets);

	// Low-level objects
	iEvent.getByToken(generalTracksToken_,generalTracks);	
	iEvent.getByToken(secondaryVerticesToken_,secondaryVertices);	
	iEvent.getByToken(PFCandsToken_, pfCands);

	// Rechits
	iEvent.getByToken(hcalRecHitsHBHEToken_, hcalRecHitsHBHE);
	// iEvent.getByToken(ebRecHitsToken_,ebRecHits); // GK remove ebRecHits for LLP HCAL skim

	// Other
	iEvent.getByToken(electron_cutbasedID_decisions_loose_Token_, electron_cutbasedID_decisions_loose); 
	iEvent.getByToken(electron_cutbasedID_decisions_medium_Token_, electron_cutbasedID_decisions_medium); 
	iEvent.getByToken(electron_cutbasedID_decisions_tight_Token_, electron_cutbasedID_decisions_tight); 
	iEvent.getByToken(photon_cutbasedID_decisions_loose_Token_, photon_cutbasedID_decisions_loose); 
	iEvent.getByToken(photon_cutbasedID_decisions_medium_Token_, photon_cutbasedID_decisions_medium); 
	iEvent.getByToken(photon_cutbasedID_decisions_tight_Token_, photon_cutbasedID_decisions_tight); 

	iEvent.getByToken(genParticlesToken_,genParticles);

	if( debug ) cout<<"Done DisplacedHcalJetNTuplizer::loadEvent"<<endl;	

}

// 
// ************************************************************************************
// Enable Output NTuple Branches
// ************************************************************************************
// 

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::EnableBranches(){

	if( debug ) cout<<"Running DisplacedHcalJetNTuplizer::EnableBranches"<<endl;

	EnableEventInfoBranches();
	EnablePVBranches();
	EnableTriggerBranches();
	EnableMetBranches();
	// Standard Objects
	EnableElectronBranches();
	EnableMuonBranches();
	EnableTauBranches();
	EnablePhotonBranches();
	EnableJetBranches();
	// Low-Level Objects
	EnableTrackBranches();
	EnableSecondaryVerticesBranches();
	EnablePFCandidateBranches();
	// Hits
	EnableEcalRechitBranches();
	EnableHcalRechitBranches();
	// MC
	//EnablePileupBranches();
	//EnableMCBranches();
	EnableGenParticleBranches();

	if( debug ) cout<<"Done DisplacedHcalJetNTuplizer::EnableBranches"<<endl;

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::EnableEventInfoBranches(){

	output_tree->Branch( "isData", &isData );
	output_tree->Branch( "runNum", &runNumber );
	output_tree->Branch( "lumiNum", &lumiNumber );
	output_tree->Branch( "eventNum", &eventNumber );
	output_tree->Branch( "eventTime", &eventTime );

}

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::EnablePVBranches(){

	output_tree->Branch( "n_PV", &n_PV );	
	output_tree->Branch( "PV_X", &PV_X );
	output_tree->Branch( "PV_Y", &PV_Y );
	output_tree->Branch( "PV_Z", &PV_Z );

	output_tree->Branch( "n_PVCand", &n_PVCand );
	output_tree->Branch( "PVCand_X", &PVCand_X );
	output_tree->Branch( "PVCand_Y", &PVCand_Y );
	output_tree->Branch( "PVCand_Z", &PVCand_Z );
	output_tree->Branch( "PVCand_LogSumPtSq", &PVCand_LogSumPtSq );
	output_tree->Branch( "PVCand_SumPx", &PVCand_SumPx );
	output_tree->Branch( "PVCand_SumPy", &PVCand_SumPy );

	output_tree->Branch( "PVTrack_Pt", &PVTrack_Pt );
	output_tree->Branch( "PVTrack_Eta", &PVTrack_Eta );
	output_tree->Branch( "PVTrack_Phi", &PVTrack_Phi );

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::EnableTriggerBranches(){

	// FIX
	//nameHLT = new std::vector<std::string>; nameHLT->clear();
	//output_tree->Branch("HLTDecision", &triggerDecision, ("HLTDecision[" + std::to_string(NTriggersMAX) +  "]/O").c_str());
	//output_tree->Branch("HLTPrescale", &triggerHLTPrescale, ("HLTPrescale[" + std::to_string(NTriggersMAX) +  "]/I").c_str());

	//output_tree->Branch("HLT_Names", &HLT_Names);
	output_tree->Branch("HLT_Decision", &HLT_Decision);
	output_tree->Branch("HLT_Prescale", &HLT_Prescale);

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::EnableMetBranches(){

	output_tree->Branch( "met_Pt", &met_Pt );
	output_tree->Branch( "met_Phi", &met_Phi );
	output_tree->Branch( "met_SumEt", &met_SumEt );

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::EnableElectronBranches(){

	output_tree->Branch( "n_ele", &n_ele );
	output_tree->Branch( "ele_Pt", &ele_Pt );
	output_tree->Branch( "ele_Eta", &ele_Eta );
	output_tree->Branch( "ele_Phi", &ele_Phi );
	output_tree->Branch( "ele_E", &ele_E );
	output_tree->Branch( "ele_Charge", &ele_Charge );
	output_tree->Branch( "ele_d0", &ele_d0 );
	output_tree->Branch( "ele_dZ", &ele_dZ );
	output_tree->Branch( "ele_dEta", &ele_dEta );
	output_tree->Branch( "ele_dPhi", &ele_dPhi );
	output_tree->Branch( "ele_passCutBasedIDLoose", &ele_passCutBasedIDLoose );
	output_tree->Branch( "ele_passCutBasedIDMedium", &ele_passCutBasedIDMedium );
	output_tree->Branch( "ele_passCutBasedIDTight", &ele_passCutBasedIDTight );
	output_tree->Branch( "ele_EcalRechitIDs", &ele_EcalRechitIDs );
	output_tree->Branch( "ele_EcalRechitIndices", &ele_EcalRechitIndices );
	output_tree->Branch( "ele_SeedRechitID", &ele_SeedRechitID );
	output_tree->Branch( "ele_SeedRechitIndex", &ele_SeedRechitIndex );

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::EnableMuonBranches(){

	output_tree->Branch( "n_muon", &n_muon );
	output_tree->Branch( "muon_Pt", &muon_Pt );
	output_tree->Branch( "muon_Eta", &muon_Eta );
	output_tree->Branch( "muon_Phi", &muon_Phi );
	output_tree->Branch( "muon_E", &muon_E );
	output_tree->Branch( "muon_Charge", &muon_Charge );
	output_tree->Branch( "muon_d0", &muon_d0 );
	output_tree->Branch( "muon_dZ", &muon_dZ );
	output_tree->Branch( "muon_IsLoose", &muon_IsLoose );
	output_tree->Branch( "muon_IsMedium", &muon_IsMedium );
	output_tree->Branch( "muon_IsTight", &muon_IsTight );

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::EnableTauBranches(){

	output_tree->Branch( "n_tau", &n_tau );
	output_tree->Branch( "tau_E", &tau_E );
	output_tree->Branch( "tau_Pt", &tau_Pt );
	output_tree->Branch( "tau_Eta", &tau_Eta );
	output_tree->Branch( "tau_Phi", &tau_Phi );
	output_tree->Branch( "tau_IsLoose", &tau_IsLoose );
	output_tree->Branch( "tau_IsMedium", &tau_IsMedium );
	output_tree->Branch( "tau_IsTight", &tau_IsTight );

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::EnablePhotonBranches(){

	output_tree->Branch( "n_pho", &n_pho );
	output_tree->Branch( "pho_Pt", &pho_Pt );
	output_tree->Branch( "pho_Eta", &pho_Eta );
	output_tree->Branch( "pho_Phi", &pho_Phi );
	output_tree->Branch( "pho_E", &pho_E );
	output_tree->Branch( "pho_sminor", &pho_sminor );
	output_tree->Branch( "pho_smajor", &pho_smajor );
	output_tree->Branch( "pho_convType", &pho_convType );
	output_tree->Branch( "pho_convTrkZ", &pho_convTrkZ );
	output_tree->Branch( "pho_passCutBasedIDLoose", &pho_passCutBasedIDLoose );
	output_tree->Branch( "pho_passCutBasedIDMedium", &pho_passCutBasedIDMedium );
	output_tree->Branch( "pho_passCutBasedIDTight", &pho_passCutBasedIDTight );
	output_tree->Branch( "pho_SeedRechitID", &pho_SeedRechitID );
	output_tree->Branch( "pho_SeedRechitIndex", &pho_SeedRechitIndex );
	output_tree->Branch( "pho_EcalRechitIDs", &pho_EcalRechitIDs );
	output_tree->Branch( "pho_EcalRechitIndices", &pho_EcalRechitIndices );

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::EnableJetBranches(){

	// AK4 PF Jets
	output_tree->Branch( "n_jet", &n_jet );
	output_tree->Branch( "jet_Pt", &jet_Pt );
	output_tree->Branch( "jet_Eta", &jet_Eta );
	output_tree->Branch( "jet_Phi", &jet_Phi );
	output_tree->Branch( "jet_E", &jet_E );
	output_tree->Branch( "jet_Mass", &jet_Mass );
	output_tree->Branch( "jet_JetArea", &jet_JetArea );
	output_tree->Branch( "jet_ChargedHadEFrac", &jet_ChargedHadEFrac );
	output_tree->Branch( "jet_NeutralHadEFrac", &jet_NeutralHadEFrac );
	output_tree->Branch( "jet_PhoEFrac", &jet_PhoEFrac );
	output_tree->Branch( "jet_EleEFrac", &jet_EleEFrac );
	output_tree->Branch( "jet_MuonEFrac", &jet_MuonEFrac );
	output_tree->Branch( "jet_ChargedHadMult", &jet_ChargedHadMult );
	output_tree->Branch( "jet_NeutralHadMult", &jet_NeutralHadMult );
	output_tree->Branch( "jet_PhoMult", &jet_PhoMult );
	output_tree->Branch( "jet_EleMult", &jet_EleMult );
	output_tree->Branch( "jet_MuonMult", &jet_MuonMult );
	output_tree->Branch( "jet_PtAllTracks", &jet_PtAllTracks );
	output_tree->Branch( "jet_PtAllPVTracks", &jet_PtAllPVTracks );
	output_tree->Branch( "jet_NVertexTracks", &jet_NVertexTracks );
	output_tree->Branch( "jet_NSelectedTracks", &jet_NSelectedTracks );
	output_tree->Branch( "jet_NSV", &jet_NSV );
	output_tree->Branch( "jet_NSVCand", &jet_NSVCand );
	output_tree->Branch( "jet_SV_x", &jet_SV_x );
	output_tree->Branch( "jet_SV_y", &jet_SV_y );
	output_tree->Branch( "jet_SV_z", &jet_SV_z );
	output_tree->Branch( "jet_SV_NTracks", &jet_SV_NTracks );
	output_tree->Branch( "jet_SV_Mass", &jet_SV_Mass );
	output_tree->Branch( "jet_SV_DRJet", &jet_SV_DRJet );
	output_tree->Branch( "jet_FlightDist2D", &jet_FlightDist2D );
	output_tree->Branch( "jet_FlightDist2DErr", &jet_FlightDist2DErr );
	output_tree->Branch( "jet_FlightDist3D", &jet_FlightDist3D );
	output_tree->Branch( "jet_FlightDist3DErr", &jet_FlightDist3DErr );
	output_tree->Branch( "jet_PassIDLoose", &jet_PassIDLoose );
	output_tree->Branch( "jet_PassIDTight", &jet_PassIDTight );
	output_tree->Branch( "jet_PartonFlavor", &jet_PartonFlavor );
	output_tree->Branch( "jet_HadronFlavor", &jet_HadronFlavor );
	output_tree->Branch( "jet_PtAllTracks_wp", &jet_PtAllTracks_wp );
	output_tree->Branch( "jet_PtAllPVTracks_wp", &jet_PtAllPVTracks_wp );
	output_tree->Branch( "jet_MedianIP_wp", &jet_MedianIP_wp );
	output_tree->Branch( "jet_MinDeltaRAllTracks_wp", &jet_MinDeltaRAllTracks_wp );
	output_tree->Branch( "jet_MinDeltaRPVTracks_wp", &jet_MinDeltaRPVTracks_wp );
	output_tree->Branch( "jet_NTracks", &jet_NTracks );
	output_tree->Branch( "jet_TrackIndices", &jet_TrackIndices );
	output_tree->Branch( "jet_NPFCands", &jet_NPFCands );
	output_tree->Branch( "jet_PFCandIndices", &jet_PFCandIndices );
	output_tree->Branch( "jet_SeedHcalRechitID", &jet_SeedHcalRechitID );
	output_tree->Branch( "jet_SeedHcalRechitIndex", &jet_SeedHcalRechitIndex );
	output_tree->Branch( "jet_EcalRechitIDs", &jet_EcalRechitIDs );
	output_tree->Branch( "jet_EcalRechitIndices", &jet_EcalRechitIndices );
	output_tree->Branch( "jet_HcalRechitIDs", &jet_HcalRechitIDs );
	output_tree->Branch( "jet_HcalRechitIndices", &jet_HcalRechitIndices );
	//output_tree->Branch( "jet_passHLTFilter", &jet_passHLTFilter );

	// AK4 Calo Jets
	output_tree->Branch( "n_calojet", &n_calojet );
	output_tree->Branch( "calojet_Pt", &calojet_Pt );
	output_tree->Branch( "calojet_Eta", &calojet_Eta );
	output_tree->Branch( "calojet_Phi", &calojet_Phi );
	output_tree->Branch( "calojet_E", &calojet_E );
	output_tree->Branch( "calojet_Mass", &calojet_Mass );
	output_tree->Branch( "calojet_NPFCands", &calojet_NPFCands );
	output_tree->Branch( "calojet_PFCandIndices", &calojet_PFCandIndices );
	output_tree->Branch( "calojet_SeedHcalRechitID", &calojet_SeedHcalRechitID );
	output_tree->Branch( "calojet_SeedHcalRechitIndex", &calojet_SeedHcalRechitIndex );
	output_tree->Branch( "calojet_EcalRechitIDs", &calojet_EcalRechitIDs );
	output_tree->Branch( "calojet_EcalRechitIndices", &calojet_EcalRechitIndices );
	output_tree->Branch( "calojet_HcalRechitIDs", &calojet_HcalRechitIDs );
	output_tree->Branch( "calojet_HcalRechitIndices", &calojet_HcalRechitIndices );
	//output_tree->Branch( "calojet_passHLTFilter", &calojet_passHLTFilter );	

	// L1 Jets
	output_tree->Branch( "n_l1jet", &n_l1jet);
	output_tree->Branch( "l1jet_Pt", &l1jet_Pt);
	output_tree->Branch( "l1jet_Eta", &l1jet_Eta);
	output_tree->Branch( "l1jet_Phi", &l1jet_Phi);
	output_tree->Branch( "l1jet_E", &l1jet_E);	
	output_tree->Branch( "l1jet_hwQual", &l1jet_hwQual);

	// AK8 PF Jets
	output_tree->Branch( "n_LRJet", &n_LRJet );
	output_tree->Branch( "LRJet_Pt", &LRJet_Pt );
	output_tree->Branch( "LRJet_Eta", &LRJet_Eta );
	output_tree->Branch( "LRJet_Phi", &LRJet_Phi );
	output_tree->Branch( "LRJet_E", &LRJet_E );
	output_tree->Branch( "LRJet_TrimmedM", &LRJet_TrimmedM );
	output_tree->Branch( "LRJet_PrunedM", &LRJet_PrunedM );
	output_tree->Branch( "LRJet_FilteredM", &LRJet_FilteredM );
	output_tree->Branch( "LRJet_SoftDropM", &LRJet_SoftDropM );
	output_tree->Branch( "LRJet_JetArea", &LRJet_JetArea );
	output_tree->Branch( "LRJet_PtAllTracks", &LRJet_PtAllTracks );
	output_tree->Branch( "LRJet_PtAllPVTracks", &LRJet_PtAllPVTracks );
	output_tree->Branch( "LRJet_NVertexTracks", &LRJet_NVertexTracks );
	output_tree->Branch( "LRJet_NSelectedTracks", &LRJet_NSelectedTracks );
	output_tree->Branch( "LRJet_ChargedHadEFrac", &LRJet_ChargedHadEFrac );
	output_tree->Branch( "LRJet_NeutralHadEFrac", &LRJet_NeutralHadEFrac );
	output_tree->Branch( "LRJet_EleEFrac", &LRJet_EleEFrac );
	output_tree->Branch( "LRJet_PhoEFrac", &LRJet_PhoEFrac );
	output_tree->Branch( "LRJet_EleMult", &LRJet_EleMult );
	output_tree->Branch( "LRJet_PhoMult", &LRJet_PhoMult );
	output_tree->Branch( "LRJet_MuonMult", &LRJet_MuonMult );
	output_tree->Branch( "LRJet_NSV", &LRJet_NSV );
	output_tree->Branch( "LRJet_NSVCand", &LRJet_NSVCand );
	output_tree->Branch( "LRJet_SV_x", &LRJet_SV_x );
	output_tree->Branch( "LRJet_SV_y", &LRJet_SV_y );
	output_tree->Branch( "LRJet_SV_z", &LRJet_SV_z );
	output_tree->Branch( "LRJet_SV_NTracks", &LRJet_SV_NTracks );
	output_tree->Branch( "LRJet_SV_Mass", &LRJet_SV_Mass );
	output_tree->Branch( "LRJet_SV_DRJet", &LRJet_SV_DRJet );
	output_tree->Branch( "LRJet_FlightDist2D", &LRJet_FlightDist2D );
	output_tree->Branch( "LRJet_FlightDist2DErr", &LRJet_FlightDist2DErr );
	output_tree->Branch( "LRJet_FlightDist3D", &LRJet_FlightDist3D );
	output_tree->Branch( "LRJet_FlightDist3DErr", &LRJet_FlightDist3DErr );
	output_tree->Branch( "LRJet_PassIDLoose", &LRJet_PassIDLoose );
	output_tree->Branch( "LRJet_PassIDTight", &LRJet_PassIDTight );
	output_tree->Branch( "LRJet_NPFCands", &LRJet_NPFCands );
	output_tree->Branch( "LRJet_PFCandIndices", &LRJet_PFCandIndices );
	output_tree->Branch( "LRJet_SeedHcalRechitID", &LRJet_SeedHcalRechitID );
	output_tree->Branch( "LRJet_SeedHcalRechitIndex", &LRJet_SeedHcalRechitIndex );
	output_tree->Branch( "LRJet_EcalRechitIDs", &LRJet_EcalRechitIDs );
	output_tree->Branch( "LRJet_EcalRechitIndices", &LRJet_EcalRechitIndices );
	output_tree->Branch( "LRJet_HcalRechitIDs", &LRJet_HcalRechitIDs );
	output_tree->Branch( "LRJet_HcalRechitIndices", &LRJet_HcalRechitIndices );
	//output_tree->Branch( "LRJet_passHLTFilter", &LRJet_passHLTFilter );	

	// AK8 Calo Jets
	output_tree->Branch( "n_caloLRJet", &n_caloLRJet );
	output_tree->Branch( "caloLRJet_Pt", &caloLRJet_Pt );
	output_tree->Branch( "caloLRJet_Eta", &caloLRJet_Eta );
	output_tree->Branch( "caloLRJet_Phi", &caloLRJet_Phi );
	output_tree->Branch( "caloLRJet_E", &caloLRJet_E );
	output_tree->Branch( "caloLRJet_Mass", &caloLRJet_Mass );
	output_tree->Branch( "caloLRJet_NPFCands", &caloLRJet_NPFCands );
	output_tree->Branch( "caloLRJet_PFCandIndices", &caloLRJet_PFCandIndices );
	output_tree->Branch( "caloLRJet_SeedHcalRechitID", &caloLRJet_SeedHcalRechitID );
	output_tree->Branch( "caloLRJet_SeedHcalRechitIndex", &caloLRJet_SeedHcalRechitIndex );
	output_tree->Branch( "caloLRJet_EcalRechitIDs", &caloLRJet_EcalRechitIDs );
	output_tree->Branch( "caloLRJet_EcalRechitIndices", &caloLRJet_EcalRechitIndices );
	output_tree->Branch( "caloLRJet_HcalRechitIDs", &caloLRJet_HcalRechitIDs );
	output_tree->Branch( "caloLRJet_HcalRechitIndices", &caloLRJet_HcalRechitIndices );
	//output_tree->Branch( "caloLRJet_passHLTFilter", &caloLRJet_passHLTFilter );

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::EnableTrackBranches(){

	output_tree->Branch( "n_track", &n_track );
	output_tree->Branch( "track_Pt", &track_Pt );
	output_tree->Branch( "track_Eta", &track_Eta );
	output_tree->Branch( "track_Phi", &track_Phi );
	output_tree->Branch( "track_Charge", &track_Charge ); 
	output_tree->Branch( "track_bestVertexIndex", &track_bestVertexIndex );
	output_tree->Branch( "track_nMissingInnerHits", &track_nMissingInnerHits );
	output_tree->Branch( "track_nMissingOuterHits", &track_nMissingOuterHits );
	output_tree->Branch( "track_nPixelHits", &track_nPixelHits );
	output_tree->Branch( "track_nHits", &track_nHits );
	output_tree->Branch( "track_angle", &track_angle );
	output_tree->Branch( "track_dxyToBS", &track_dxyToBS );
	output_tree->Branch( "track_dxyErr", &track_dxyErr );
	output_tree->Branch( "track_dzToPV", &track_dzToPV );
	output_tree->Branch( "track_dzErr", &track_dzErr );
	output_tree->Branch( "track_chi2", &track_chi2 );
	output_tree->Branch( "track_ndof", &track_ndof );
	output_tree->Branch( "track_index", &track_index );

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::EnableSecondaryVerticesBranches(){

	output_tree->Branch( "n_secVtx", &n_secVtx );
	output_tree->Branch( "secVtx_Pt", &secVtx_Pt );
	output_tree->Branch( "secVtx_Eta", &secVtx_Eta );
	output_tree->Branch( "secVtx_Phi", &secVtx_Phi );
	output_tree->Branch( "secVtx_Charge", &secVtx_Charge );
	output_tree->Branch( "secVtx_nConstituents", &secVtx_nConstituents );
	output_tree->Branch( "secVtx_X", &secVtx_X );
	output_tree->Branch( "secVtx_Y", &secVtx_Y );
	output_tree->Branch( "secVtx_Z", &secVtx_Z );
	output_tree->Branch( "secVtx_Distance", &secVtx_Distance );
	output_tree->Branch( "secVtx_DistanceError", &secVtx_Distance );

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::EnablePFCandidateBranches(){

	PFCand_BranchesEnabled = true;

	output_tree->Branch( "n_PFCand", &n_PFCand );
	output_tree->Branch( "PFCand_Pt", &PFCand_Pt );
	output_tree->Branch( "PFCand_Eta", &PFCand_Eta );
	output_tree->Branch( "PFCand_Phi", &PFCand_Phi );
	output_tree->Branch( "PFCand_PdgId", &PFCand_PdgId );
	output_tree->Branch( "PFCand_TrackIndex", &PFCand_TrackIndex );
	output_tree->Branch( "PFCand_GeneralTrackIndex", &PFCand_GeneralTrackIndex );
	output_tree->Branch( "PFCand_PVIndex", &PFCand_PVIndex );

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::EnableEcalRechitBranches(){

	EcalRechit_BranchesEnabled = true;

	output_tree->Branch( "n_ecalRechit", &n_ecalRechit );
	output_tree->Branch( "ecalRechit_Eta", &ecalRechit_Eta );
	output_tree->Branch( "ecalRechit_Phi", &ecalRechit_Phi );
	output_tree->Branch( "ecalRechit_E", &ecalRechit_E );
	output_tree->Branch( "ecalRechit_X", &ecalRechit_X );
	output_tree->Branch( "ecalRechit_Y", &ecalRechit_Y );
	output_tree->Branch( "ecalRechit_Z", &ecalRechit_Z );	
	output_tree->Branch( "ecalRechit_time", &ecalRechit_time );
	output_tree->Branch( "ecalRechit_E_err", &ecalRechit_E_err );
	output_tree->Branch( "ecalRechit_time_err", &ecalRechit_time_err );
	output_tree->Branch( "ecalRechit_kSaturatedflag", &ecalRechit_kSaturatedflag );
	output_tree->Branch( "ecalRechit_kLeadingEdgeRecoveredflag", &ecalRechit_kLeadingEdgeRecoveredflag );
	output_tree->Branch( "ecalRechit_kPoorRecoflag", &ecalRechit_kPoorRecoflag );
	output_tree->Branch( "ecalRechit_kWeirdflag", &ecalRechit_kWeirdflag );
	output_tree->Branch( "ecalRechit_kDiWeirdflag", &ecalRechit_kDiWeirdflag );

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::EnableHcalRechitBranches(){

	HcalRechit_BranchesEnabled = true;

	output_tree->Branch( "n_hbheRechit", &n_hbheRechit );
	output_tree->Branch( "hbheRechit_Eta", &hbheRechit_Eta );
	output_tree->Branch( "hbheRechit_Phi", &hbheRechit_Phi );
	output_tree->Branch( "hbheRechit_E", &hbheRechit_E );
	output_tree->Branch( "hbheRechit_X", &hbheRechit_X );
	output_tree->Branch( "hbheRechit_Y", &hbheRechit_Y );
	output_tree->Branch( "hbheRechit_Z", &hbheRechit_Z );
	output_tree->Branch( "hbheRechit_time", &hbheRechit_time );
	output_tree->Branch( "hbheRechit_auxTDC", &hbheRechit_auxTDC );
	output_tree->Branch( "hbheRechit_timeFalling", &hbheRechit_timeFalling );
	//output_tree->Branch( "hbheRechit_cctime", &hbheRechit_cctime );
	output_tree->Branch( "hbheRechit_iEta", &hbheRechit_iEta );
	output_tree->Branch( "hbheRechit_iPhi", &hbheRechit_iPhi );
	output_tree->Branch( "hbheRechit_depth", &hbheRechit_depth );

	output_tree->Branch( "n_hoRechit", &n_hoRechit );
	output_tree->Branch( "hoRechit_Eta", &hoRechit_Eta );
	output_tree->Branch( "hoRechit_Phi", &hoRechit_Phi );
	output_tree->Branch( "hoRechit_E", &hoRechit_E );
	output_tree->Branch( "hoRechit_X", &hoRechit_X );
	output_tree->Branch( "hoRechit_Y", &hoRechit_Y );
	output_tree->Branch( "hoRechit_Z", &hoRechit_Z );
	output_tree->Branch( "hoRechit_T", &hoRechit_T );

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::EnablePileupBranches(){

	/*output_tree->Branch( "nBunchXing", &nBunchXing );
	output_tree->Branch( "BunchXing", &BunchXing );
	output_tree->Branch( "nPU", &nPU );
	output_tree->Branch( "nPUmean", &nPUmean );*/

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::EnableMCBranches(){

/*  output_tree->Branch( "nGenJets", &nGenJets );
	output_tree->Branch( "genJetE", &genJetE );
	output_tree->Branch( "genJetPt", &genJetPt );
	output_tree->Branch( "genJetEta", &genJetEta );
	output_tree->Branch( "genJetPhi", &genJetPhi );
	output_tree->Branch( "genMetPtCalo", &genMetPtCalo );
	output_tree->Branch( "genMetPhiCalo", &genMetPhiCalo );
	output_tree->Branch( "genMetPtTrue", &genMetPtTrue );
	output_tree->Branch( "genMetPhiTrue", &genMetPhiTrue );
	output_tree->Branch( "genVertexX", &genVertexX );
	output_tree->Branch( "genVertexY", &genVertexY );
	output_tree->Branch( "genVertexZ", &genVertexZ );
	output_tree->Branch( "genVertexT", &genVertexT );
	output_tree->Branch( "genWeight", &genWeight );
	output_tree->Branch( "genSignalProcessID", &genSignalProcessID );
	output_tree->Branch( "genQScale", &genQScale );
	output_tree->Branch( "genAlphaQCD", &genAlphaQCD );
	output_tree->Branch( "genAlphaQED", &genAlphaQED );
	output_tree->Branch( "lheComments", &lheComments );

	/* //scaleWeights = new std::vector<float>; scaleWeights->clear();
	//pdfWeights = new std::vector<float>; pdfWeights->clear();
	//alphasWeights = new std::vector<float>; alphasWeights->clear();
	output_tree->Branch("scaleWeights", &scaleWeights);
	output_tree->Branch("pdfWeights", &pdfWeights);
	output_tree->Branch("alphasWeights", &alphasWeights);
	*/

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::EnableGenParticleBranches(){

	output_tree->Branch( "n_gParticle", &n_gParticle );
	output_tree->Branch( "gParticle_ParentIndex", &gParticle_ParentIndex );
	output_tree->Branch( "gParticle_ParentId", &gParticle_ParentId );
	output_tree->Branch( "gParticle_Status", &gParticle_Status );
	output_tree->Branch( "gParticle_Id", &gParticle_Id );

	output_tree->Branch( "gParticle_Pt", &gParticle_Pt );
	output_tree->Branch( "gParticle_Px", &gParticle_Px );
	output_tree->Branch( "gParticle_Py", &gParticle_Py );
	output_tree->Branch( "gParticle_Pz", &gParticle_Pz );
	output_tree->Branch( "gParticle_Eta", &gParticle_Eta );
	output_tree->Branch( "gParticle_Phi", &gParticle_Phi );
	output_tree->Branch( "gParticle_E", &gParticle_E );
	output_tree->Branch( "gParticle_ProdVtx_X", &gParticle_ProdVtx_X );
	output_tree->Branch( "gParticle_ProdVtx_Y", &gParticle_ProdVtx_Y );
	output_tree->Branch( "gParticle_ProdVtx_Z", &gParticle_ProdVtx_Z );
	//output_tree->Branch( "gParticle_DecayVtx_X", &gParticle_DecayVtx_X );
	//output_tree->Branch( "gParticle_DecayVtx_Y", &gParticle_DecayVtx_Y );
	//output_tree->Branch( "gParticle_DecayVtx_Z", &gParticle_DecayVtx_Z );

	output_tree->Branch( "n_gLLP", &n_gLLP ); 
	output_tree->Branch( "gLLP_Pt", &gLLP_Pt ); 
	output_tree->Branch( "gLLP_Eta", &gLLP_Eta ); 
	output_tree->Branch( "gLLP_Phi", &gLLP_Phi ); 
	output_tree->Branch( "gLLP_E", &gLLP_E ); 
	output_tree->Branch( "gLLP_Beta", &gLLP_Beta ); 
	output_tree->Branch( "gLLP_TravelTime", &gLLP_TravelTime ); 
	output_tree->Branch( "gLLP_DecayVtx_X", &gLLP_DecayVtx_X ); 
	output_tree->Branch( "gLLP_DecayVtx_Y", &gLLP_DecayVtx_Y ); 
	output_tree->Branch( "gLLP_DecayVtx_Z", &gLLP_DecayVtx_Z ); 	
	output_tree->Branch( "gLLP_ProdVtx_X", &gLLP_ProdVtx_X ); 
	output_tree->Branch( "gLLP_ProdVtx_Y", &gLLP_ProdVtx_Y ); 
	output_tree->Branch( "gLLP_ProdVtx_Z", &gLLP_ProdVtx_Z ); 

};

// 
// ************************************************************************************
// Reset Output NTuple Branches
// ************************************************************************************
// 

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::ResetBranches(){
	// Called by the loadEvent() method

	if( debug ) cout<<"Running DisplacedHcalJetNTuplizer::ResetBranches"<<endl; 

	// Event Level Info
	ResetEventInfoBranches();
	ResetPVBranches();
	ResetTriggerBranches();
	ResetMetBranches();
	// Standard Objects
	ResetElectronBranches();
	ResetMuonBranches();
	ResetTauBranches();
	ResetPhotonBranches();
	ResetJetBranches();
	// Low-Level Objects
	ResetTrackBranches();
	ResetPFCandidateBranches();
	ResetSecondaryVerticesBranches();
	// Hits
	ResetEcalRechitBranches();
	ResetHcalRechitBranches();
	// MC
	ResetPileupBranches();
	ResetMCBranches();
	ResetGenParticleBranches();

	if( debug ) cout<<"Done DisplacedHcalJetNTuplizer::ResetBranches"<<endl;  

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::ResetEventInfoBranches(){

	isData 			= false;
	runNumber 	= 0;
	lumiNumber 	= 0;
	eventNumber = 0;
	eventTime 	= 0;

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::ResetPVBranches(){

	n_PV 		= -9999;
	PV_X 	= -9999.9;
	PV_Y  = -9999.9;
	PV_Z  = -9999.9;

	n_PVCand = 0;
	PVCand_X.clear();
	PVCand_Y.clear();
	PVCand_Z.clear();
	PVCand_LogSumPtSq.clear();
	PVCand_SumPx.clear();
	PVCand_SumPy.clear();

	n_PVTrack = 0;
	PVTrack_Pt.clear();
	PVTrack_Eta.clear();
	PVTrack_Phi.clear();

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::ResetTriggerBranches(){

	//HLT_Names.clear();
	HLT_Decision.clear();
	HLT_Prescale.clear();

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::ResetMetBranches(){

	met_Pt = -9999.9;
	met_Phi = -9999.9;
	met_SumEt = -9999.9;

}

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::ResetElectronBranches(){

	n_ele = 0;
	ele_Pt.clear();
	ele_Eta.clear();
	ele_Phi.clear();
	ele_E.clear();
	ele_Charge.clear();
	ele_d0.clear();
	ele_dZ.clear();
	ele_dEta.clear();
	ele_dPhi.clear();
	ele_passCutBasedIDLoose.clear();
	ele_passCutBasedIDMedium.clear();
	ele_passCutBasedIDTight.clear();
	ele_SeedRechitID.clear();
	ele_SeedRechitIndex.clear();
	ele_EcalRechitIDs.clear();
	ele_EcalRechitIndices.clear();
	//ele_passHLTFilter.clear();

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::ResetMuonBranches(){

	n_muon = 0; 
	muon_Pt.clear();
	muon_Eta.clear();
	muon_Phi.clear();
	muon_E.clear();
	muon_Charge.clear();
	muon_d0.clear();
	muon_dZ.clear();
	muon_IsLoose.clear();
	muon_IsMedium.clear();
	muon_IsTight.clear();
	//muon_passHLTFilter.clear();

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::ResetTauBranches(){

	n_tau = 0;
	tau_Pt.clear();
	tau_Eta.clear();
	tau_Phi.clear();
	tau_E.clear();
	tau_IsLoose.clear();
	tau_IsMedium.clear();
	tau_IsTight.clear();

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::ResetPhotonBranches(){

	n_pho = 0; 
	pho_Pt.clear();
	pho_Eta.clear();
	pho_Phi.clear();
	pho_E.clear();
	pho_sminor.clear();
	pho_smajor.clear();
	pho_convType.clear();
	pho_convTrkZ.clear();
	pho_passCutBasedIDLoose.clear();
	pho_passCutBasedIDMedium.clear();
	pho_passCutBasedIDTight.clear();
	pho_SeedRechitID.clear();
	pho_SeedRechitIndex.clear();
	pho_EcalRechitIDs.clear();
	pho_EcalRechitIndices.clear();
	//pho_passHLTFilter.clear();

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::ResetJetBranches(){

	// AK4 PF Jets 
	n_jet = 0;
	jet_Pt.clear();
	jet_Eta.clear();
	jet_Phi.clear();
	jet_E.clear();
	jet_Mass.clear();
	jet_JetArea.clear();
	jet_ChargedHadEFrac.clear();
	jet_NeutralHadEFrac.clear();
	jet_PhoEFrac.clear();
	jet_EleEFrac.clear();
	jet_MuonEFrac.clear();
	jet_ChargedHadMult.clear();
	jet_NeutralHadMult.clear();
	jet_PhoMult.clear();
	jet_EleMult.clear();
	jet_MuonMult.clear();
	jet_PtAllTracks.clear();
	jet_PtAllPVTracks.clear();
	jet_NVertexTracks.clear();
	jet_NSelectedTracks.clear();
	jet_NSV.clear();
	jet_NSVCand.clear();
	jet_SV_x.clear();
	jet_SV_y.clear();
	jet_SV_z.clear();
	jet_SV_NTracks.clear();
	jet_SV_Mass.clear();
	jet_SV_DRJet.clear();
	jet_FlightDist2D.clear();
	jet_FlightDist2DErr.clear();
	jet_FlightDist3D.clear();
	jet_FlightDist3DErr.clear();
	jet_PassIDLoose.clear();
	jet_PassIDTight.clear();
	jet_PartonFlavor.clear();
	jet_HadronFlavor.clear();
	jet_PtAllTracks_wp.clear();
	jet_PtAllPVTracks_wp.clear();
	jet_MedianIP_wp.clear();
	jet_MinDeltaRAllTracks_wp.clear();
	jet_MinDeltaRPVTracks_wp.clear();
	jet_TrackIndices.clear();
	jet_NPFCands.clear();
	jet_PFCandIndices.clear();
	jet_SeedHcalRechitID.clear();
	jet_SeedHcalRechitIndex.clear();
	jet_EcalRechitIDs.clear();
	jet_EcalRechitIndices.clear();
	jet_HcalRechitIDs.clear();
	jet_HcalRechitIndices.clear();

	// AK4 Calo Jets 
	n_calojet = 0;
	calojet_Pt.clear();
	calojet_Eta.clear();
	calojet_Phi.clear();
	calojet_E.clear();
	calojet_Mass.clear();
	calojet_NPFCands.clear();
	calojet_PFCandIndices.clear();
	calojet_SeedHcalRechitID.clear();
	calojet_SeedHcalRechitIndex.clear();
	calojet_EcalRechitIDs.clear();
	calojet_EcalRechitIndices.clear();
	calojet_HcalRechitIDs.clear();
	calojet_HcalRechitIndices.clear();

	// L1 jets
	n_l1jet = 0;
	l1jet_Pt.clear();
	l1jet_Eta.clear();
	l1jet_Phi.clear();
	l1jet_E.clear();
	l1jet_hwQual.clear();

	// AK8 PF Jets 
	n_LRJet = 0;
	LRJet_Pt.clear();
	LRJet_Eta.clear();
	LRJet_Phi.clear();
	LRJet_E.clear();
	LRJet_TrimmedM.clear();
	LRJet_PrunedM.clear();
	LRJet_FilteredM.clear();
	LRJet_SoftDropM.clear();
	LRJet_JetArea.clear();
	LRJet_PtAllTracks.clear();
	LRJet_PtAllPVTracks.clear();
	LRJet_NVertexTracks.clear();
	LRJet_NSelectedTracks.clear();
	LRJet_ChargedHadEFrac.clear();
	LRJet_NeutralHadEFrac.clear();
	LRJet_EleEFrac.clear();
	LRJet_PhoEFrac.clear();
	LRJet_EleMult.clear();
	LRJet_PhoMult.clear();
	LRJet_MuonMult.clear();
	LRJet_NSV.clear();
	LRJet_NSVCand.clear();
	LRJet_SV_x.clear();
	LRJet_SV_y.clear();
	LRJet_SV_z.clear();
	LRJet_SV_NTracks.clear();
	LRJet_SV_Mass.clear();
	LRJet_SV_DRJet.clear();
	LRJet_FlightDist2D.clear();
	LRJet_FlightDist2DErr.clear();
	LRJet_FlightDist3D.clear();
	LRJet_FlightDist3DErr.clear();
	LRJet_PassIDLoose.clear();
	LRJet_PassIDTight.clear();
	LRJet_NPFCands.clear();
	LRJet_PFCandIndices.clear();
	LRJet_SeedHcalRechitID.clear();
	LRJet_SeedHcalRechitIndex.clear();
	LRJet_EcalRechitIDs.clear();
	LRJet_EcalRechitIndices.clear();
	LRJet_HcalRechitIDs.clear();
	LRJet_HcalRechitIndices.clear();

	// AK8 Calo Jets 
	n_caloLRJet = 0;
	caloLRJet_Pt.clear();
	caloLRJet_Eta.clear();
	caloLRJet_Phi.clear();
	caloLRJet_E.clear();
	caloLRJet_Mass.clear();
	caloLRJet_NPFCands.clear();
	caloLRJet_PFCandIndices.clear();
	caloLRJet_SeedHcalRechitID.clear();
	caloLRJet_SeedHcalRechitIndex.clear();
	caloLRJet_EcalRechitIDs.clear();
	caloLRJet_EcalRechitIndices.clear();
	caloLRJet_HcalRechitIDs.clear();
	caloLRJet_HcalRechitIndices.clear();

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::ResetTrackBranches(){

	n_track = 0;
	track_Pt.clear();
	track_Eta.clear();
	track_Phi.clear();
	track_Charge.clear();
	track_bestVertexIndex.clear();
	track_nMissingInnerHits.clear();
	track_nMissingOuterHits.clear();
	track_nPixelHits.clear();
	track_nHits.clear();
	track_angle.clear();
	track_dxyToBS.clear();
	track_dxyErr.clear();
	track_dzToPV.clear();
	track_dzErr.clear();
	track_chi2.clear();
	track_ndof.clear();
	track_index.clear();

}

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::ResetSecondaryVerticesBranches(){

	n_secVtx = 0;
	secVtx_Pt.clear();
	secVtx_Eta.clear();
	secVtx_Phi.clear();
	secVtx_Charge.clear();
	secVtx_nConstituents.clear();
	secVtx_X.clear();
	secVtx_Y.clear();
	secVtx_Z.clear();
	secVtx_Distance.clear();
	secVtx_DistanceError.clear();

}

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::ResetPFCandidateBranches(){

	n_PFCand = 0;
	PFCand_PdgId.clear();
	PFCand_Pt.clear();
	PFCand_Eta.clear();
	PFCand_Phi.clear();
	PFCand_TrackIndex.clear();
	PFCand_GeneralTrackIndex.clear();
	PFCand_PVIndex.clear();

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::ResetEcalRechitBranches(){

	n_ecalRechit = 0;
	ecalRechit_Eta.clear();
	ecalRechit_Phi.clear();
	ecalRechit_E.clear();
	ecalRechit_X.clear();
	ecalRechit_Y.clear();
	ecalRechit_Z.clear();	
	ecalRechit_time.clear();
	ecalRechit_E_err.clear();
	ecalRechit_time_err.clear();
	ecalRechit_kSaturatedflag.clear();
	ecalRechit_kLeadingEdgeRecoveredflag.clear();
	ecalRechit_kPoorRecoflag.clear();
	ecalRechit_kWeirdflag.clear();
	ecalRechit_kDiWeirdflag.clear();

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::ResetHcalRechitBranches(){

	n_hbheRechit = 0;
	hbheRechit_Phi.clear();
	hbheRechit_Eta.clear();
	hbheRechit_X.clear();
	hbheRechit_Y.clear();
	hbheRechit_Z.clear();
	hbheRechit_E.clear();
	hbheRechit_time.clear();
	hbheRechit_auxTDC.clear();
	hbheRechit_timeFalling.clear();
	//hbheRechit_cctime.clear();
	hbheRechit_iEta.clear();
	hbheRechit_iPhi.clear();
	hbheRechit_depth.clear();

	n_hoRechit = 0;
	hoRechit_Phi.clear();
	hoRechit_Eta.clear();
	hoRechit_X.clear();
	hoRechit_Y.clear();
	hoRechit_Z.clear();
	hoRechit_E.clear();
	hoRechit_T.clear();

}

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::ResetPileupBranches(){

	/*nBunchXing = 0;
	BunchXing.clear();
	nPU.clear();
	nPUmean.clear();*/

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::ResetMCBranches(){

	/*
	nGenJets = 0;
	for ( int i = 0; i < OBJECTARRAYSIZE; i++ )
	{
		genJetE[i] = -999.;
		genJetPt[i] = -999.;
		genJetEta[i] = -999.;
		genJetPhi[i] = -999.;
		genJetMET[i] = -999.;
		//genJet_match_jet_index[i] = 666;
		//genJet_min_delta_r_match_jet[i] = -666.;
	}

	genMetPtCalo  = -999.;
	genMetPhiCalo = -999.;
	genMetPtTrue  = -999.;
	genMetPhiTrue = -999.;

	genVertexX = -999.;
	genVertexY = -999.;
	genVertexZ = -999.;
	genVertexT = -999.;

	genWeight = -999.;
	genSignalProcessID = 0;
	genQScale = -999.;
	genAlphaQCD = -999.;
	genAlphaQED = -999.;

	if (scaleWeights) scaleWeights->clear();
	if (pdfWeights) pdfWeights->clear();
	if (alphasWeights) alphasWeights->clear();
		*/

	return;
};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::ResetGenParticleBranches(){

	n_gParticle = 0;
	gParticle_Id.clear();
	gParticle_Status.clear();
	gParticle_ParentId.clear();
	gParticle_ParentIndex.clear();

	gParticle_Pt.clear();
	gParticle_Px.clear();
	gParticle_Py.clear();
	gParticle_Pz.clear();
	gParticle_Eta.clear();
	gParticle_Phi.clear();
	gParticle_E.clear();
	gParticle_ProdVtx_X.clear();
	gParticle_ProdVtx_Y.clear();
	gParticle_ProdVtx_Z.clear();
	gParticle_DecayVtx_X.clear();
	gParticle_DecayVtx_Y.clear();
	gParticle_DecayVtx_Z.clear();

	n_gLLP = 0;
	gLLP_Pt.clear();
	gLLP_Eta.clear();
	gLLP_Phi.clear();
	gLLP_E.clear();
	gLLP_Beta.clear();
	gLLP_TravelTime.clear();
	gLLP_ProdVtx_X.clear();
	gLLP_ProdVtx_Y.clear();
	gLLP_ProdVtx_Z.clear();
	gLLP_DecayVtx_X.clear();
	gLLP_DecayVtx_Y.clear();
	gLLP_DecayVtx_Z.clear();

	return;
};

// 
// ************************************************************************************
// Protected Classes to Execute
// ************************************************************************************
// 

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::beginJob(){
	// Method called once each job just before starting event loop

	if( debug ) cout<<"Running DisplacedHcalJetNTuplizer::beginJob"<<endl; 

	// ----- Enable Branches ----- //

	EnableBranches();

	if( debug ) cout<<"Done DisplacedHcalJetNTuplizer::beginJob"<<endl; 

}

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
	// Method called for each run:
	// read LHE header if present and determine which weights to read for pdf and alphas uncertainties based on the
	// central pdf set used
	// This is semi-hardcoded for now to work with current centrally produced samples
	// covering SUSY signal samples, LO madgraph, NLO madgraph_aMC@NLO, and NLO powheg
	// generated with nnpdf30
	// More robust selection will require some parsing of <initrwgt> block

	if( debug ) cout<<"Running DisplacedHcalJetNTuplizer::beginRun"<<endl; 
	
	// Code		

	if( debug ) cout<<"Done DisplacedHcalJetNTuplizer::beginRun"<<endl; 

}

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup) {
	// Method called for each lumi block

	if( debug ) cout<<"Running DisplacedHcalJetNTuplizer::beginLuminosityBlock"<<endl; 

	//if (useGen_) {
	//	iLumi.getByToken(genLumiHeaderToken_,genLumiHeader);
	//}

	//fill lhe comment lines with SUSY model parameter information
	string lheComments_temp = "";
	/* FIX if (genLumiHeader.isValid()) {
		lheComments_temp = genLumiHeader->configDescription();
	}*/
	cout<<"Begin Luminosity Block -- LHEComments: "<<lheComments_temp<<endl;

	if( debug ) cout<<"Done DisplacedHcalJetNTuplizer::beginLuminosityBlock"<<endl; 

}

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup) {
	if (debug) cout << "Running DisplacedHcalJetNTuplizer::endLuminosityBlock" << endl;
}

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
	if( debug ) cout<<"Running DisplacedHcalJetNTuplizer::endRun"<<endl; 
}

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
	// Method called for each event

	if( debug ) cout<<"Running DisplacedHcalJetNTuplizer::analyze"<<endl; 

	using namespace edm;

	// ----- Initialize Event ----- //

	if( debug ) cout<<"Initializing Event..."<<endl;

	loadEvent(iEvent); // Load objects and resets tree branches
	NEvents->Fill(0); // Event Count

	// ----- Reset Branches ----- //

	if( debug ) cout<<"Resetting Branches..."<<endl;

	ResetBranches();

	// ----- Getting CaloGeometry ----- //

	if( debug ) cout<<"Getting CaloGeometry..."<<endl;

	auto const& geoHandle = iSetup.getData(caloGeometryToken_);
	caloGeometry_EB = geoHandle.getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
	caloGeometry_EE = geoHandle.getSubdetectorGeometry(DetId::Ecal, EcalEndcap);
	caloGeometry_HB = geoHandle.getSubdetectorGeometry(DetId::Hcal, HcalBarrel);
	caloGeometry_HE = geoHandle.getSubdetectorGeometry(DetId::Hcal, HcalEndcap);
	caloGeometry_HO = geoHandle.getSubdetectorGeometry(DetId::Hcal, HcalOuter); 

	// ----- Fill Branches ----- //

	if( debug ) cout<<"Filling Branches..."<<endl;

	//FillBranches();

	// Event Level Info
	FillEventInfoBranches( iEvent );
	FillPVBranches( iEvent );
	FillTriggerBranches( iEvent );
	FillMetBranches( iEvent );

	// Standard Objects
	FillElectronBranches( iEvent );
	FillMuonBranches( iEvent );
	FillPhotonBranches( iEvent );
	FillJetBranches( iEvent, iSetup );

	// Low-Level Objects
	FillTrackBranches( iEvent ); //, iSetup );

	// Hits
	FillEcalRechitBranches( iEvent, iSetup ); 
	FillHcalRechitBranches( iEvent, iSetup );

	// MC
	if( !isData ){
		FillGenParticleBranches();
	}

	//fillPVAll();
	//if ( enableTriggerInfo_ ) 
	//	fillTrigger( iEvent );

	output_tree->Fill();

	if( debug ) cout<<"Done DisplacedHcalJetNTuplizer::analyze"<<endl; 

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::endJob(){
	// Method called once each job just after ending the event loop
};

// 
// ************************************************************************************
// Fill Output NTuple Branches
// ************************************************************************************
// 

// ------------------------------------------------------------------------------------
bool DisplacedHcalJetNTuplizer::FillEventInfoBranches(const edm::Event& iEvent){
	// Fill tree variables

	if( debug ) cout<<"Running DisplacedHcalJetNTuplizer::FillEventInfoBranches"<<endl; 

	isData 		= isData_;
	runNumber 	= iEvent.id().run();
	lumiNumber 	= iEvent.luminosityBlock();
	eventNumber = iEvent.id().event();
	eventTime 	= iEvent.eventAuxiliary().time().unixTime();

	if( debug ) cout<<"Done DisplacedHcalJetNTuplizer::FillEventInfoBranches"<<endl; 

	return true;

};

// ------------------------------------------------------------------------------------
bool DisplacedHcalJetNTuplizer::FillPVBranches( const edm::Event& iEvent ){

	if( debug ) cout<<"Running DisplacedHcalJetNTuplizer::FillPVBranches"<<endl; 

	// Select the primary vertex (if any) //

	n_PV 		 = 0;
	PV_global 	 = &(vertices->front());
	bool foundPV = false;

	for( unsigned int i = 0; i < vertices->size(); i++ ){
		if(vertices->at(i).isValid() && !vertices->at(i).isFake()){
			if (!foundPV){
				PV_global = &(vertices->at(i));
				foundPV = true;
			}
			n_PV++;
		}
	}

	PV_X = PV_global->x();
	PV_Y = PV_global->y();
	PV_Z = PV_global->z();

	// get rho -- FIX
	// fixedGridRhoFastjetAll = *rhoFastjetAll;

	// Get All PVs //

	n_PVCand = vertices->size(); 
	for (int ipv = 0; ipv < n_PVCand; ++ipv){
		const reco::Vertex &vtx_temp = vertices->at(ipv);
		PVCand_X.push_back( vtx_temp.x() );
		PVCand_Y.push_back( vtx_temp.y() );
		PVCand_Z.push_back( vtx_temp.z() );

		PVCand_LogSumPtSq.push_back( log(0) ); // FIX
		PVCand_SumPx.push_back( 0 ); // FIX
		PVCand_SumPy.push_back( 0 ); // FIX
	} 

	if( debug ) cout<<"Done DisplacedHcalJetNTuplizer::FillPVBranches"<<endl; 

	return true;

};


// ------------------------------------------------------------------------------------
bool DisplacedHcalJetNTuplizer::FillTriggerBranches(const edm::Event& iEvent){

	if( debug ) cout<<"Running DisplacedHcalJetNTuplizer::FillTriggerBranches"<<endl; 

	//fill trigger information
	const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

	if( debug ) cout<<"N TRIGGER BITS = "<<triggerBits->size()<<endl;

	for( uint i = 0; i < triggerBits->size(); ++i) { // TODO: CHECK OUT INDEXING

		if( debug && triggerBits->accept(i) ) cout<<"PASS TRIGGER "<<names.triggerName(i)<<endl;

		string hltPathNameReq = "HLT_";
		
		if( (names.triggerName(i)).find(hltPathNameReq) == string::npos ) continue;
		if( (names.triggerName(i)).find_last_of("_") == string::npos ) continue;

		int lastUnderscorePos = (names.triggerName(i)).find_last_of("_");
		string hltPathNameWithoutVersionNumber = (names.triggerName(i)).substr(0,lastUnderscorePos);

		//for( uint j = 0; j < triggerPathNames->size(); ++j) {
		for( auto triggerPathName: triggerPathNames ){
			if( triggerPathName == "" ) continue;
			if( hltPathNameWithoutVersionNumber == triggerPathName ){
				//HLT_Names.push_back( triggerPathName );
				HLT_Decision.push_back( triggerBits->accept(i) );
				HLT_Prescale.push_back( -1 );

				if( triggerBits->accept(i) )
					NEvents_HLT->Fill(triggerPathNamesIndices[triggerPathName]); // FIX WEIGHTS
				//if (isData_) HLT_Prescale.push_back( triggerPrescales->getPrescaleForIndex(i) );
				//else HLT_Prescale.push_back( 1 ); // TODO Need to figure out yields
			}
		}
	}

	if( debug ) cout<<"Done DisplacedHcalJetNTuplizer::FillTriggerBranches"<<endl; 

	return true;
};

// ------------------------------------------------------------------------------------
bool DisplacedHcalJetNTuplizer::FillMetBranches(const edm::Event& iEvent){
	
	if( debug ) cout<<"Running DisplacedHcalJetNTuplizer::FillMetBranches"<<endl; 

	//const reco::PFMET &Met = mets->front();
	//const pat::MET &Met = mets->front();
	auto &met_temp = met->front();
	met_Pt = met_temp.pt();
	met_Phi = met_temp.phi();
	met_SumEt = met_temp.sumEt();

	if( debug ) cout<<"Done DisplacedHcalJetNTuplizer::FillMetBranches"<<endl; 

	return true;
};

// ------------------------------------------------------------------------------------
bool DisplacedHcalJetNTuplizer::FillElectronBranches(const edm::Event& iEvent){

	if( debug ) cout<<"Running DisplacedHcalJetNTuplizer::FillElectronBranches"<<endl; 

	// ----- Loop Over Electrons ------ // 

	//for(const pat::Electron &ele : *electrons) {
	for( uint i = 0; i < electrons->size(); ++i){

		if( debug ) cout<<" -- ele idx "<<i<<endl;

		const reco::GsfElectron ele = (*electrons)[i];
		//const auto eleRef = electrons->ptrAt(i);
		reco::GsfElectronRef eleRef(electrons, i);

		if( ele.pt() < 5 ) continue;

		n_ele++; 

		// Basics //

		ele_Pt.push_back( ele.pt() );
		ele_Eta.push_back( ele.eta() );
		ele_Phi.push_back( ele.phi() );
		ele_E.push_back( ele.energy() );
		ele_Charge.push_back( ele.charge() );

		// Features //

		ele_d0.push_back( -ele.gsfTrack().get()->dxy(PV_global->position()) );
		ele_dZ.push_back( ele.gsfTrack().get()->dz(PV_global->position()) );
		ele_dEta.push_back( ele.deltaEtaSuperClusterTrackAtVtx() - ele.superCluster()->eta() + ele.superCluster()->seed()->eta() );
		ele_dPhi.push_back( ele.deltaPhiSuperClusterTrackAtVtx() );

		// ID //

		ele_passCutBasedIDLoose.push_back( (*electron_cutbasedID_decisions_loose)[eleRef] );
		ele_passCutBasedIDMedium.push_back( (*electron_cutbasedID_decisions_medium)[eleRef] );
		ele_passCutBasedIDTight.push_back( (*electron_cutbasedID_decisions_tight)[eleRef] );

		// Rechits Association //

		ele_SeedRechitID.push_back( 0 ); // FIXURGENT
		ele_SeedRechitIndex.push_back( 0 ); // FIXURGENT
		ele_EcalRechitIDs.push_back( vector<uint>{0} ); // FIXURGENT
		ele_EcalRechitIndices.push_back( vector<uint>{0} ); // FIXURGENT

		// Other //

		// FIX ele_HoverE.push_back( ele.hcalOverEcal() );

	}

	if( debug ) cout<<"Done DisplacedHcalJetNTuplizer::FillElectronBranches"<<endl; 	

	return true;

};

// ------------------------------------------------------------------------------------
bool DisplacedHcalJetNTuplizer::FillMuonBranches(const edm::Event& iEvent) {

	if( debug ) cout<<"Running DisplacedHcalJetNTuplizer::FillMuonBranches"<<endl; 	

	// ----- Loop Over Muons ----- // 

	for( const pat::Muon &muon : *muons) {

		if( debug ) cout<<" -- muon idx "<<n_muon<<endl;

		n_muon++;

		// Basics // 

		muon_Pt.push_back( muon.pt() );
		muon_Eta.push_back( muon.eta() );
		muon_Phi.push_back( muon.phi() );
		muon_E.push_back( muon.energy() );

		// Features //

		muon_Charge.push_back( muon.charge() );
		muon_d0.push_back( -muon.muonBestTrack()->dxy(PV_global->position()) );
		muon_dZ.push_back( muon.muonBestTrack()->dz(PV_global->position()) );

		// ID //

		muon_IsLoose.push_back( muon.isLooseMuon() );
		muon_IsMedium.push_back( muon.isMediumMuon() );
		muon_IsTight.push_back( muon.isTightMuon(*PV_global) );

	}

	if( debug ) cout<<"Done DisplacedHcalJetNTuplizer::FillMuonBranches"<<endl; 	

	return true;
	
};

// ------------------------------------------------------------------------------------
bool DisplacedHcalJetNTuplizer::FillTauBranches(const edm::Event& iEvent){

	if( debug ) cout<<"Running DisplacedHcalJetNTuplizer::FillTauBranches"<<endl; 	

	for ( auto &tau : *taus) {

		if (tau.pt() < 18) continue;

		if( debug ) cout<<" -- tau idx "<<n_tau<<endl;

		n_tau++;

		// Basics //

		tau_Pt.push_back( tau.pt() );
		tau_Eta.push_back( tau.eta() );
		tau_Phi.push_back( tau.phi() );
		tau_E.push_back( tau.energy() );

		// ID //

		tau_IsLoose.push_back( bool(tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")) );
		tau_IsMedium.push_back( bool(tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits")) );
		tau_IsTight.push_back( bool(tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits")) );

	} 

	if( debug ) cout<<"Done DisplacedHcalJetNTuplizer::FillTauBranches"<<endl; 		

	return true;

};

// ------------------------------------------------------------------------------------
bool DisplacedHcalJetNTuplizer::FillPhotonBranches( const edm::Event& iEvent ){

	if( debug ) cout<<"Running DisplacedHcalJetNTuplizer::FillPhotonBranches"<<endl; 		

	//noZS::EcalClusterLazyTools *lazyToolnoZS = new noZS::EcalClusterLazyTools(iEvent, iSetup, ebRecHitsToken_, eeRecHitsToken_);
	//for (const reco::Photon &pho : *photons) {
	for (uint i=0; i < photons->size() ; i++) {
		const reco::Photon pho = (*photons)[i];
		reco::PhotonRef phoRef(photons, i);

		if( debug ) cout<<" -- pho idx "<<n_tau<<endl;

		if( pho.pt() < 20 ) continue;

		n_pho++;

		// Basics //

		pho_E.push_back( pho.energy() );
		pho_Pt.push_back( pho.pt() );
		pho_Eta.push_back( pho.eta() ); //correct this for the vertex
		pho_Phi.push_back( pho.phi() ); //correct this for the vertex

		// Features //

		pho_sminor.push_back( -1 ); // FIXURGENT 
		pho_smajor.push_back( -1 ); // FIXURGENT 
		pho_convType.push_back( -1 ); // FIX
		pho_convTrkZ.push_back( -1 ); // FIX

		// ID  // 

		pho_passCutBasedIDLoose.push_back( (*photon_cutbasedID_decisions_loose)[phoRef] );
		pho_passCutBasedIDMedium.push_back( (*photon_cutbasedID_decisions_medium)[phoRef] );
		pho_passCutBasedIDTight.push_back( (*photon_cutbasedID_decisions_tight)[phoRef] );

		// Rechits Association // 
		pho_SeedRechitID.push_back( 0 ); // FIXURGENT
		pho_SeedRechitIndex.push_back( 0 ); // FIXURGENT
		pho_EcalRechitIDs.push_back( vector<uint>{0} ); // FIXURGENT
		pho_EcalRechitIndices.push_back( vector<uint>{0} ); // FIXURGENT

	}

	if( debug ) cout<<"Done DisplacedHcalJetNTuplizer::FillPhotonBranches"<<endl; 		

	return true;

};

// ------------------------------------------------------------------------------------
bool DisplacedHcalJetNTuplizer::FillJetBranches( const edm::Event& iEvent, const edm::EventSetup& iSetup ){

	if( debug ) cout<<"Running DisplacedHcalJetNTuplizer::FillJetBranches"<<endl; 		

	// FIX URGENT
	//  edm::ESHandle<CaloGeometry> geoHandle;
	//  iSetup.get<CaloGeometryRecord>().get(geoHandle);
	// const CaloSubdetectorGeometry *barrelGeometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
	// const CaloSubdetectorGeometry *endcapGeometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);
	//  const CaloSubdetectorGeometry *hbGeometry = geoHandle->getSubdetectorGeometry(DetId::Hcal, HcalBarrel);
	//  const CaloSubdetectorGeometry *heGeometry = geoHandle->getSubdetectorGeometry(DetId::Hcal, HcalEndcap);
	// const CaloSubdetectorGeometry *hoGeometry = geoHandle->getSubdetectorGeometry(DetId::Hcal, HcalOuter); 

	// ********************************************************
	// AK4 Jets
	// ********************************************************

	if( debug ) cout<<" -> AK4 PF Jets"<<endl;  

	//for (const pat::Jet &j : *jets) {
	for( auto &jet : *jets ) {

		if( jet.pt() < 10 || fabs(jet.eta()) > 1.5 ) continue;

		if( debug ) cout<<" ------ jet idx"<<n_jet<<endl;  

		n_jet++;
		// if (fabs(j.eta()) > 2.4) continue;
		//NO (KK): only fill more detailed information for jets with pT >= 20
		//if (j.pt() >= 20) {
		//TLorentzVector thisJet;
		//thisJet.SetPtEtaPhiE(jetPt[nJets], jetEta[nJets], jetPhi[nJets], jetE[nJets]);

		// ----- Basics ----- // 

		jet_E.push_back( jet.energy() );
		jet_Pt.push_back( jet.pt() );
		jet_Eta.push_back( jet.eta() );
		jet_Phi.push_back( jet.phi() );
		jet_Mass.push_back( jet.mass() );

		// ----- ID ----- //

		if( debug ) cout<<" ------ 1"<<endl;  

		jet_PassIDLoose.push_back( 0 ); // FIX passJetID(&j, 0) );
		jet_PassIDTight.push_back( 0 ); // FIX passJetID(&j, 1) );

		if( !isData_ ) {
			jet_PartonFlavor.push_back( jet.partonFlavour() );
			jet_HadronFlavor.push_back( jet.hadronFlavour() );
		} else {
			jet_PartonFlavor.push_back( -9999 );
			jet_HadronFlavor.push_back( -9999 );      
		}

		// ----- General Features ----- // 

		if( debug ) cout<<" ------ 2"<<endl;  

		jet_JetArea.push_back( jet.jetArea() );

		jet_ChargedHadEFrac.push_back( jet.chargedHadronEnergyFraction() );
		jet_NeutralHadEFrac.push_back( jet.neutralHadronEnergyFraction() );
		jet_PhoEFrac.push_back( jet.neutralEmEnergyFraction() );
		jet_EleEFrac.push_back( jet.electronEnergyFraction() );
		jet_MuonEFrac.push_back( jet.muonEnergyFraction() );
		jet_ChargedHadMult.push_back( jet.chargedHadronMultiplicity() );
		jet_NeutralHadMult.push_back( jet.neutralHadronMultiplicity() );
		jet_PhoMult.push_back( jet.photonMultiplicity() );
		jet_EleMult.push_back( jet.electronMultiplicity() );
		jet_MuonMult.push_back( jet.muonMultiplicity() );
		jet_PtAllTracks.push_back( -1. ); // FIX
		jet_PtAllPVTracks.push_back( -1. ); // FIX
		jet_NVertexTracks.push_back( -1. ); // FIX
		jet_NSelectedTracks.push_back( -1. ); // FIX

		//jet_PileupE.push_back( jet.pileup() );
		//jet_PileupId.push_back( jet.userFloat("pileupJetId:fullDiscriminant") );
		//jet_PileupIdFlag.push_back( jet.userInt("pileupJetId:fullId") ); //A bit map for loose, medium, and tight working points

		// ----- Secondary Vertex Features ----- //

		if( debug ) cout<<" ------ 3"<<endl;  

		if( jet.tagInfoLabels().size() > 0 && jet.hasTagInfo("pfSecondaryVertex") ){

			auto jet_SVCandInfo_temp = jet.tagInfoCandSecondaryVertex("pfSecondaryVertex");

			jet_NSV.push_back( jet_SVCandInfo_temp->nVertices() );
			jet_NSVCand.push_back( jet_SVCandInfo_temp->nVertexCandidates() );
			//jet.tagInfoCandSecondaryVertex("pfSecondaryVertex")->nVertexTracks()
			//jet.tagInfoCandSecondaryVertex("pfSecondaryVertex")->nSelectedTracks() 

			if( jet_SVCandInfo_temp->nVertices() > 0 ){
				jet_SV_x.push_back( jet_SVCandInfo_temp->secondaryVertex(0).vx() );
				jet_SV_y.push_back( jet_SVCandInfo_temp->secondaryVertex(0).vy() );
				jet_SV_z.push_back( jet_SVCandInfo_temp->secondaryVertex(0).vz() );
				jet_SV_NTracks.push_back( jet_SVCandInfo_temp->nVertexTracks(0) );
				jet_SV_Mass.push_back( jet_SVCandInfo_temp->secondaryVertex(0).mass() );
				jet_FlightDist2D.push_back( jet_SVCandInfo_temp->flightDistance(0, true).value() );
				jet_FlightDist2DErr.push_back( jet_SVCandInfo_temp->flightDistance(0, true).error() );
				jet_FlightDist3D.push_back( jet_SVCandInfo_temp->flightDistance(0, false).value() );
				jet_FlightDist3DErr.push_back( jet_SVCandInfo_temp->flightDistance(0, false).error() );
			} else{
				jet_SV_x.push_back( -9999.9 );
				jet_SV_y.push_back( -9999.9 );
				jet_SV_z.push_back( -9999.9 );
				jet_SV_NTracks.push_back( -9999 );
				jet_SV_Mass.push_back( -9999.9 );
				jet_FlightDist2D.push_back( -9999.9 );
				jet_FlightDist2DErr.push_back( -9999.9 );
				jet_FlightDist3D.push_back( -9999.9 );
				jet_FlightDist3DErr.push_back( -9999.9 );
			}

			if( jet_SVCandInfo_temp->taggingVariables().getList(reco::btau::vertexJetDeltaR,false).size() > 0 ){
				jet_SV_DRJet.push_back( jet_SVCandInfo_temp->taggingVariables().getList(reco::btau::vertexJetDeltaR,false)[0] );
			} else {
				jet_SV_DRJet.push_back( -9999.9 );
			}

		} else {
			jet_NSV.push_back( -9999 );
			jet_NSVCand.push_back( -9999 );
			jet_SV_x.push_back( -9999.9 );
			jet_SV_y.push_back( -9999.9 );
			jet_SV_z.push_back( -9999.9 );
			jet_SV_NTracks.push_back( -9999 );
			jet_SV_Mass.push_back( -9999.9 );
			jet_FlightDist2D.push_back( -9999.9 );
			jet_FlightDist2DErr.push_back( -9999.9 );
			jet_FlightDist3D.push_back( -9999.9 );
			jet_FlightDist3DErr.push_back( -9999.9 );
			jet_SV_DRJet.push_back( -9999.9 );

		}

		// ----- Trackless Variables ----- //

		// FIXURGENT

		// ----- Find PFCandidates Inside Jet ----- //

		if( debug ) cout<<" ------ 4"<<endl;  

		vector<uint> jet_PFCandIndices_temp; 

		for( uint i_pfcand = 0; i_pfcand < pfCands->size(); i_pfcand++ ){
			//const reco::PFCandidate *p = &(*pfCands)[i_pfcand];
			reco::PFCandidatePtr p_ptr(pfCands, i_pfcand);

			for( uint i_constituent=0; i_constituent < jet.getPFConstituents().size(); i_constituent++ ){
				if( p_ptr == jet.getPFConstituents()[i_constituent] ) jet_PFCandIndices_temp.push_back(i_pfcand);
			}

		}

		jet_NPFCands.push_back( jet.getPFConstituents().size() ); 
		jet_PFCandIndices.push_back( jet_PFCandIndices_temp ); 

		// ----- Find Tracks Inside Jet ----- //

		if( debug ) cout<<" ------ 5"<<endl;  

		vector<uint> jet_TrackIndices_temp; 

		for( uint it = 0; it < generalTracks->size(); it ++){
			reco::Track generalTrack = generalTracks->at(it);

			// DR
			if( deltaR( jet.eta(), jet.phi(), generalTrack.eta(), generalTrack.phi()) > 0.5 ) continue;
			// Energy 
			if( generalTrack.pt() < 1 ) continue;

			jet_TrackIndices_temp.push_back( it );

		}
		jet_NTracks.push_back( jet_TrackIndices_temp.size() );
		jet_TrackIndices.push_back( jet_TrackIndices_temp );

		// ----- Find Ecal Rechits Inside Jet ----- //

		if( debug ) cout<<" ------ 6"<<endl;

		// ----- Find Hcal Rechits Inside Jet ----- //

		if( debug ) cout<<" ------ 7 ( size of hcalRecHitsHBHE: "<<hcalRecHitsHBHE->size()<<" )"<<endl;  

		vector<uint> jet_HcalRechitIDs_temp;
		vector<uint> jet_HcalRechitIndices_temp;

		for( uint ih = 0; ih < hcalRecHitsHBHE->size(); ih++ ){

			if( debug ) cout<<" -------- hcalRecHitsHBHE loop: "<<ih<<endl;  

			const HBHERecHit *recHit = &(*hcalRecHitsHBHE)[ih];
			const DetId recHitId = recHit->detid();

			// Barrel
			if( recHit->detid().subdetId() != HcalBarrel ) continue;
			const auto recHitPos = caloGeometry_HB->getGeometry(recHitId)->getPosition();

			// DR
			if( deltaR( jet.eta(), jet.phi(), recHitPos.eta(), recHitPos.phi()) > 0.5) continue; 
			// Energy 
			if( recHit->energy() < 0.1 ) continue;

			jet_HcalRechitIDs_temp.push_back( 0 ); // FIX DetId );
			jet_HcalRechitIndices_temp.push_back( ih );

			//cout << "Save this Rechit: " << q << " | " << SaveThisEBRechit[q] << " : " << recHit->energy() << " " << recHitPos.eta() << " " << recHitPos.phi() << "\n";	
		
		}
		jet_HcalRechitIDs.push_back( jet_HcalRechitIDs_temp );
		jet_HcalRechitIndices.push_back( jet_HcalRechitIndices_temp );

		if( debug ) cout<<" ------ 8"<<endl;  

	} //loop over jets

	// ********************************************************
	// L1 Jets
	// ******************************************************** 

	for ( auto &l1jet : *l1jets) {

		if( l1jet.pt() < 10 || fabs(l1jet.eta()) > 1.5 ) continue;

		if( debug ) cout<<" ------ L1 jet idx"<<n_l1jet<<endl;  

		n_l1jet++;
		
		// ----- Basics ----- // 

		l1jet_E.push_back( l1jet.energy() );
		l1jet_Pt.push_back( l1jet.pt() );
		l1jet_Eta.push_back( l1jet.eta() );
		l1jet_Phi.push_back( l1jet.phi() );
		l1jet_hwQual.push_back( l1jet.hwQual() );
	}

	if( debug ) cout<<"Done DisplacedHcalJetNTuplizer::FillJetBranches"<<endl; 			

	return true;
};

// ------------------------------------------------------------------------------------
bool DisplacedHcalJetNTuplizer::FillTrackBranches( const edm::Event& iEvent ){

	if( debug ) cout<<"Running DisplacedHcalJetNTuplizer::FillJetBranches"<<endl; 			

	edm::Handle<reco::BeamSpot> beamSpot;
	iEvent.getByToken(bsToken_, beamSpot);

	for( uint it = 0; it < generalTracks->size(); it ++ ){

		// reco::Track generalTrack = generalTracks->at(iTrack);
		// //const auto& generalTrack = generalTracks->at(iTrack);
		// // TLorentzVector generalTrackVecTemp;
		// // generalTrackVecTemp.SetPtEtaPhiM(generalTrack.pt(),generalTrack.eta(),generalTrack.phi(),0);

		// Check if we should save
		bool save_track = false; 
		for( int ij = 0; ij < n_jet; ij++ ){
			for( int ijt = 0; ijt < (int)jet_TrackIndices.at(ij).size(); ijt++ ){
				if( it == jet_TrackIndices.at(ij).at(ijt) ){
					save_track = true; 
					break;
				}
			}
			if( save_track == true ) break;
		}

		if( !save_track ) continue;

		n_track++;

		//reco::TrackBaseRef tref(generalTracks, it);
		reco::Track track = generalTracks->at(it);

		// make transient track (unfolding effects of B field ?)
		//reco::TransientTrack tt(generalTrackHandle->at(iTrack),magneticFieldGeometry);
		//if(!tt.isValid()) continue;

		track_Pt.push_back( track.pt() );
		track_Eta.push_back( track.eta() );
		track_Phi.push_back( track.phi() );
		
		track_Charge.push_back( track.charge() );

		// Find the best vertex for this track
		float maxWeight = 0;
		int bestVertexIndex = -1;
		/* FIX for( int k = 0; k < (int)vertices->size();k++ ){
			if( vertices->at(k).trackWeight(track) > maxWeight ){
				maxWeight = vertices->at(k).trackWeight(track);
				bestVertexIndex = k;
			}
		}*/

		track_bestVertexIndex.push_back( bestVertexIndex );
		track_nMissingInnerHits.push_back( track.hitPattern().numberOfLostTrackerHits(reco::HitPattern::MISSING_INNER_HITS) );
		track_nMissingOuterHits.push_back( track.hitPattern().numberOfLostTrackerHits(reco::HitPattern::MISSING_OUTER_HITS) );
		track_nPixelHits.push_back( track.hitPattern().numberOfValidPixelHits() );
		track_nHits.push_back( track.hitPattern().numberOfValidHits() );
		track_dxyToBS.push_back( track.dxy(beamSpot->position()) ); 
		track_dxyErr.push_back( track.dxyError() );
		track_dzToPV.push_back( track.dz(beamSpot->position()) );
		track_dzErr.push_back( track.dzError() ); 
		track_chi2.push_back( track.chi2() ); 
		track_ndof.push_back( track.ndof() ); 
		track_index.push_back( it ); 

	} //loop over tracks

	if( debug ) cout<<"Done DisplacedHcalJetNTuplizer::FillJetBranches"<<endl; 			

	return true; 

}

// ------------------------------------------------------------------------------------
bool DisplacedHcalJetNTuplizer::FillPFCandidateBranches( const edm::Event& iEvent ){

	if( debug ) cout<<"Running DisplacedHcalJetNTuplizer::FillJetBranches"<<endl; 			


	if( debug ) cout<<"Done DisplacedHcalJetNTuplizer::FillJetBranches"<<endl; 			

	return true; 
	
}

// ------------------------------------------------------------------------------------
bool DisplacedHcalJetNTuplizer::FillSecondaryVerticesBranches( const edm::Event& iEvent ){

	if( debug ) cout<<"Running DisplacedHcalJetNTuplizer::FillJetBranches"<<endl; 			


	if( debug ) cout<<"Done DisplacedHcalJetNTuplizer::FillJetBranches"<<endl; 			

	return true; 
	
}

bool DisplacedHcalJetNTuplizer::FillEcalRechitBranches(const edm::Event& iEvent, const edm::EventSetup& iSetup){

	if( debug ) cout<<"Running DisplacedHcalJetNTuplizer::FillEcalRechitBranches"<<endl; 			

	if( debug ) cout<<"Done DisplacedHcalJetNTuplizer::FillEcalRechitBranches"<<endl; 			

	return true; 

}

bool DisplacedHcalJetNTuplizer::FillHcalRechitBranches(const edm::Event& iEvent, const edm::EventSetup& iSetup){

	if( debug ) cout<<"Running DisplacedHcalJetNTuplizer::FillHcalRechitBranches"<<endl; 			

	for( uint ih = 0; ih < hcalRecHitsHBHE->size(); ih++ ){

		if( debug ) cout<<" -- hcalrechit idx "<<ih<<endl;

		const HBHERecHit *recHit = &(*hcalRecHitsHBHE)[ih];
		const HcalDetId recHitId = recHit->detid();
		if( recHit->detid().subdetId() != HcalBarrel ) continue;

		// Check if we should save
		bool save_hit = false; 
		save_hit = true;

		/* Save all hits for now

		for( int ij = 0; ij < n_jet; ij++ ){
			
			for( int ijh = 0; ijh < (int)jet_HcalRechitIndices.at(ij).size(); ijh++ ){
				if( ih == jet_HcalRechitIndices.at(ij).at(ijh) ){
					save_hit = true; 
					break;
				}
			}

			if( save_hit == true ) continue;
		}*/

		if( !save_hit ) continue;

		if (recHit->energy() < 0.5) continue;

		const auto recHitPos = caloGeometry_HB->getGeometry(recHitId)->getPosition();

		n_hbheRechit++;

		hbheRechit_Eta.push_back( recHitPos.eta() );
		hbheRechit_Phi.push_back( recHitPos.phi() );
		hbheRechit_E.push_back( recHit->energy() );
		// Features  
		hbheRechit_time.push_back( recHit->time() );
		hbheRechit_auxTDC.push_back( recHit->auxTDC() );
		hbheRechit_timeFalling.push_back( recHit->timeFalling() );
		//hbheRechit_cctime.push_back( );
		hbheRechit_X.push_back( recHitPos.x() );
		hbheRechit_Y.push_back( recHitPos.y() );
		hbheRechit_Z.push_back( recHitPos.z() );
		hbheRechit_iEta.push_back( recHitId.ieta() );
		hbheRechit_iPhi.push_back( recHitId.iphi() );
		hbheRechit_depth.push_back( recHitId.depth() );

	} 

	if( debug ) cout<<"Done DisplacedHcalJetNTuplizer::FillHcalRechitBranches"<<endl; 			

	return true; 
}


// ------------------------------------------------------------------------------------
bool DisplacedHcalJetNTuplizer::FillPileupBranches(){

	if( debug ) cout<<"Running DisplacedHcalJetNTuplizer::FillPileup"<<endl; 

	/*for(const PileupSummaryInfo &pu : *puInfo){
		nBunchXing++;
		BunchXing.push_back( pu.getBunchCrossing() );
		nPU.push_back( pu.getPU_NumInteractions() );
		nPUmean.push_back( pu.getTrueNumInteractions() );
	}*/
	
	if( debug ) cout<<"Done DisplacedHcalJetNTuplizer::FillPileup"<<endl; 

	return true;
};

// ------------------------------------------------------------------------------------
bool DisplacedHcalJetNTuplizer::FillMCBranches() {

	return true;

};

// ------------------------------------------------------------------------------------
bool DisplacedHcalJetNTuplizer::FillGenParticleBranches(){

	if( debug ) cout<<"Running DisplacedHcalJetNTuplizer::FillGenParticleBranches"<<endl; 

	//Fills selected gen particles
	const double pt_cut = 0.0;

	vector<int> LLP_PdgIds;
	LLP_PdgIds.push_back(9000006);
	LLP_PdgIds.push_back(9000007);
	LLP_PdgIds.push_back(1023);
	LLP_PdgIds.push_back(1000023);
	LLP_PdgIds.push_back(1000025);
	LLP_PdgIds.push_back(6000113);
	LLP_PdgIds.push_back(9900012);
	LLP_PdgIds.push_back(9900014);
	LLP_PdgIds.push_back(9900016);

	//Allows easier comparison for mother finding
	std::vector<const reco::Candidate*> prunedV;

	if( debug ) cout<<"Looping over genParticles"<<endl;

	for( size_t i=0; i < genParticles->size(); i++ ){

		if( debug ) cout<<" -- genParticle idx "<<i<<endl;

		bool already_saved = false;

		// TODO: Add higgs status requirement (pythia progression) -- want to last one

		if( abs( (*genParticles)[i].pdgId() ) == 25 && (*genParticles)[i].status() != 62 ){ 
			if( debug ) cout<<"GenParticle Higgs STatus: "<<(*genParticles)[i].status()<<", isPromptFinalState():"<<(*genParticles)[i].isPromptFinalState()<<", isPromptDecayed():"<<(*genParticles)[i].isPromptDecayed()<<endl;
			continue;
		}
			

		if(
			(abs((*genParticles)[i].pdgId()) >= 1 && abs((*genParticles)[i].pdgId()) <= 6 && ( (*genParticles)[i].status() < 30 ))
			|| (abs((*genParticles)[i].pdgId()) >= 11 && abs((*genParticles)[i].pdgId()) <= 16)
			|| (abs((*genParticles)[i].pdgId()) == 21 && (*genParticles)[i].status() < 30)
			|| (abs((*genParticles)[i].pdgId()) == 22 && (*genParticles)[i].pt() > 10.0 )
			|| (abs((*genParticles)[i].pdgId()) >= 23 && abs((*genParticles)[i].pdgId()) <= 25)
			|| (abs((*genParticles)[i].pdgId()) >= 32 && abs((*genParticles)[i].pdgId()) <= 42)
			//|| (abs((*genParticles)[i].pdgId()) >= 100 && abs((*genParticles)[i].pdgId()) <= 350)
			|| (abs((*genParticles)[i].pdgId()) == 1023)
			|| (abs((*genParticles)[i].pdgId()) >= 1000001 && abs((*genParticles)[i].pdgId()) <= 1000039)
			|| (abs((*genParticles)[i].pdgId()) == 9000006 || abs((*genParticles)[i].pdgId()) == 9000007)
			|| (abs((*genParticles)[i].pdgId()) == 6000113 )
			|| (abs((*genParticles)[i].pdgId()) == 9900012 || abs((*genParticles)[i].pdgId()) == 9900014 || abs((*genParticles)[i].pdgId()) == 9900016)
			
			) {

			if( (*genParticles)[i].pt() > pt_cut ){
				prunedV.push_back(&(*genParticles)[i]);
				already_saved = true;
				continue;
			}
		}

		if( already_saved ) continue;
	

	} // End loop over genParticles

	// Total number of gen particles
	n_gParticle = prunedV.size();

	if( debug ) cout<<"Looping over prunedV"<<endl;	

	//Look for mother particle and Fill gen variables
	for( unsigned int i = 0; i < prunedV.size(); i++ ) {

		if( debug ) cout<<" -- prunedV idx "<<i<<endl;

		n_gParticle++;

		gParticle_Id.push_back( prunedV[i]->pdgId() );
		gParticle_Status.push_back( prunedV[i]->status() );

		gParticle_ParentId.push_back( prunedV[i]->mother()->pdgId() ); // prunedV[i]->mother() this returns a pointer to mother particle
		for( unsigned int j=0; j < prunedV.size(); j++ ) {
		  if (prunedV[j]->pdgId() == prunedV[i]->mother()->pdgId() && prunedV[j]->vx() == prunedV[i]->mother()->vx() && prunedV[j]->vy() == prunedV[i]->mother()->vy() && prunedV[j]->vz() == prunedV[i]->mother()->vz() && prunedV[j]->pt() == prunedV[i]->mother()->pt()) {
		    gParticle_ParentIndex.push_back(j);
		  }
		  else {
		    gParticle_ParentIndex.push_back( -1 ); // if no mother particle was found to match -- could be the case because not all gen particles are saved in prunedV
		  }
		}
		    
		gParticle_Pt.push_back( prunedV[i]->pt() );
		gParticle_Px.push_back( prunedV[i]->px() );
		gParticle_Py.push_back( prunedV[i]->py() );
		gParticle_Pz.push_back( prunedV[i]->pz() );
		gParticle_Eta.push_back( prunedV[i]->eta() );
		gParticle_Phi.push_back( prunedV[i]->phi() );
		gParticle_E.push_back( prunedV[i]->energy() );

		gParticle_ProdVtx_X.push_back( prunedV[i]->vx() );
		gParticle_ProdVtx_Y.push_back( prunedV[i]->vy() );
		gParticle_ProdVtx_Z.push_back( prunedV[i]->vz() );

		// ----- Gen LLP Info ----- //

		if( debug ) cout<<" ------ 1"<<endl;

		// Match with one of the entries in the llpIDs List
		bool matchedLLP = false;
		int matchedLLPID = 0;
		if( gParticle_Status.at(i) == 22 ) {
			for( uint d=0; d < LLP_PdgIds.size() ; d++ ){
				if( abs(gParticle_Id.at(i)) != LLP_PdgIds[d] ) continue;

				matchedLLPID = gParticle_Id.at(i);
				matchedLLP = true;
				
			}
		}

		if( debug ) cout<<" ------ 2"<<endl;

		if( !matchedLLP ) continue;

		n_gLLP++;

		gLLP_Pt.push_back( gParticle_Pt.at(i) ); 
		gLLP_Eta.push_back( gParticle_Eta.at(i) );
		gLLP_Phi.push_back( gParticle_Phi.at(i) );
		gLLP_E.push_back( gParticle_E.at(i) );

		float beta_temp = sqrt(gParticle_Px.at(i)*gParticle_Px.at(i)+gParticle_Py.at(i)*gParticle_Py.at(i)+gParticle_Pz.at(i)*gParticle_Pz.at(i))/gParticle_E.at(i);
		gLLP_Beta.push_back( beta_temp );

		// Get Decay Vertex Info 

		if( debug ) cout<<" ------ 3"<<endl;

		const reco::Candidate *llp_child = 0;
		bool found_llp_child = false;
		const reco::Candidate *tmpParticle = prunedV[i];		

		while( true ){
			if( tmpParticle->numberOfDaughters() == 0 ) break;
			llp_child = tmpParticle->daughter(0);
			if( llp_child && (llp_child->pdgId() != matchedLLPID) ){
				found_llp_child = true;
				break;
			} else {
				tmpParticle = llp_child;
			}
		}

		if( debug ) cout<<" ------ 4"<<endl;

		gLLP_ProdVtx_X.push_back( prunedV[i]->vx() );
		gLLP_ProdVtx_Y.push_back( prunedV[i]->vy() );
		gLLP_ProdVtx_Z.push_back( prunedV[i]->vz() );		

		if( found_llp_child ){
			gLLP_DecayVtx_X.push_back( llp_child->vx() ); // note this is in cm
			gLLP_DecayVtx_Y.push_back( llp_child->vy() );
			gLLP_DecayVtx_Z.push_back( llp_child->vz() );
		} else {
			gLLP_DecayVtx_X.push_back( -9999.9 );
			gLLP_DecayVtx_Y.push_back( -9999.9 );
			gLLP_DecayVtx_Z.push_back( -9999.9 );			
		}

		if( debug ) cout<<" ------ 5"<<endl;

		//cout<<i<<" "<<gLLP_DecayVtx_X.size()<<"  "<<gLLP_Beta.size()<<endl;

		float travel_time_temp = sqrt(
				  pow( gLLP_DecayVtx_X.back() - gLLP_ProdVtx_X.back(),2)
				+ pow( gLLP_DecayVtx_Y.back() - gLLP_ProdVtx_Y.back(),2)
				+ pow( gLLP_DecayVtx_Z.back() - gLLP_ProdVtx_Z.back(),2) )
				/(30. * gLLP_Beta.back() ); //1/30 is to convert cm to ns
		
		if( debug ) cout<<" ------ 6"<<endl;		

		if( !found_llp_child ) travel_time_temp = -9999.9;

		gLLP_TravelTime.push_back( travel_time_temp );		

		if( debug ) cout<<" ------ 7"<<endl;

	} // end loop over prunced v

	if( debug ) cout<<"Done DisplacedHcalJetNTuplizer::FillGenParticleBranches"<<endl; 

	return true;

};

double DisplacedHcalJetNTuplizer::deltaPhi(double phi1, double phi2)
{
	double dphi = phi1-phi2;
	while (dphi > TMath::Pi())
	{
		dphi -= TMath::TwoPi();
	}
	while (dphi <= -TMath::Pi())
	{
		dphi += TMath::TwoPi();
	}
	return dphi;
};

double DisplacedHcalJetNTuplizer::deltaR(double eta1, double phi1, double eta2, double phi2)
{
	double dphi = deltaPhi(phi1,phi2);
	double deta = eta1 - eta2;
	return sqrt( dphi*dphi + deta*deta);
};


//define this as a plug-in
DEFINE_FWK_MODULE(DisplacedHcalJetNTuplizer);