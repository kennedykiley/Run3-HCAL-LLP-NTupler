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
	rand_(0), // random seed // GK for JER
	// Trigger
	triggerBitsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"))),
	triggerObjectsToken_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
	triggerPrescalesToken_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("triggerPrescales"))),
	triggerPathNamesFile_(iConfig.getParameter<string>( "triggerPathNamesFile" )),
	//eleHLTFilterNamesFile_(iConfig.getParameter<string>( "electronHLTFilterNamesFile" )),
	//muonHLTFilterNamesFile_(iConfig.getParameter<string>( "muonHLTFilterNamesFile" )),
	//phoHLTFilterNamesFile_(iConfig.getParameter<string>( "photonHLTFilterNamesFile" )),
	//jetHLTFilterNamesFile_(iConfig.getParameter<string>( "jetHLTFilterNamesFile" )),
	// General 
	//metFilterBitsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterBits"))),
    metFilterBitsToken_(
        iConfig.getParameter<edm::InputTag>("metFilterBits").label().empty() ?
        edm::EDGetTokenT<edm::TriggerResults>() :  // default-initialized token
        mayConsume<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterBits"))
    ),
	verticesToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
	primaryVertexAssociationToken_(consumes<edm::Association<vector<reco::Vertex> > >(edm::InputTag("primaryVertexAssociation","original"))),
	primaryVertexAssociationValueMapToken_(consumes<edm::ValueMap<int> >(edm::InputTag("primaryVertexAssociation","original"))),
  	rhoFastjetAllToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetAll"))),
	// Event-Level Info
	metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("met"))),
	metPuppiToken_(consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("metPuppi"))), // GK, PUPPI MET access
	bsTag_(iConfig.getUntrackedParameter<edm::InputTag>("offlineBeamSpot", edm::InputTag("offlineBeamSpot"))),
	bsToken_(consumes<reco::BeamSpot>(bsTag_)),
	// Physics Objects
	electronsToken_(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
	muonsToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
	tausToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
	photonsToken_(consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
	jetsToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("pfjetsAK4"))),
	jetsCorrToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("pfjetsAK4_corrected"))), // GK, adding JECs
	jetsPuppiCorrToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("pfjetsAK4Puppi"))), // GK, adding JECs for PUPPI jets
	calojetsToken_(consumes<reco::CaloJetCollection>(iConfig.getParameter<edm::InputTag>("calojetsAK4"))),
	LRJetsToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("pfjetsAK8"))),
	caloLRJetsToken_(consumes<reco::CaloJetCollection>(iConfig.getParameter<edm::InputTag>("calojetsAK8"))),
	l1jetsToken_(consumes<BXVector<l1t::Jet>>(iConfig.getParameter<edm::InputTag>("l1jets"))), // GK added for L1 jets access
	// // also add JER 
	// jerResToken_CHS_(esConsumes<JME::JetResolution, JetResolutionRcd>(edm::ESInputTag("", "AK4PFchs_pt"))),
	// jerSFToken_CHS_(esConsumes<JME::JetResolutionScaleFactor, JetResolutionScaleFactorRcd>(edm::ESInputTag("", "AK4PFchs"))),
	// jerResToken_(esConsumes<JME::JetResolution, JetResolutionRcd>(edm::ESInputTag("", "AK4PFPuppi_pt"))),
	// jerSFToken_(esConsumes<JME::JetResolutionScaleFactor, JetResolutionScaleFactorRcd>(edm::ESInputTag("", "AK4PFPuppi"))),	
	// Low-Level Objects
	//tracksToken_(consumes<edm::View<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracks"))),
	generalTracksToken_(consumes<std::vector<reco::Track>>(edm::InputTag("generalTracks"))),
	secondaryVerticesToken_(consumes<vector<reco::VertexCompositePtrCandidate> >(iConfig.getParameter<edm::InputTag>("secondaryVertices"))),	
	PFCandsToken_(consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
	// RecHits
	// ebRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("ebRecHits"))),
	// hcalRecHitsHBHEToken_(consumes<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>(edm::InputTag("reducedHcalRecHits","hbhereco"))),
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
    puInfoToken_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("puInfo"))),
	// Geometry
	caloGeometryToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()), // GK
	castorGeometryToken_(esConsumes<PCaloGeometry, PCastorRcd>()), // GK
	// gemGeoToken_(esConsumes<RecoIdealGeometry, GEMRecoGeometryRcd>()), // commenting out to test 350 GeV MC
	magneticFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()), // GK
	propagatorToken_(esConsumes<Propagator, TrackingComponentsRecord>(edm::ESInputTag("", "PropagatorWithMaterial"))), // GK // thePropagatorToken(esConsumes<Propagator, TrackingComponentsRecord>(edm::ESInputTag("", thePropagatorName))),

	// Noise Filters
	Flag_HBHENoiseFilterToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_HBHENoiseFilter"))),
	Flag_HBHENoiseIsoFilterToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_HBHENoiseIsoFilter"))),
	Flag_CSCTightHaloFilterToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_CSCTightHaloFilter"))),
	Flag_CSCTightHaloTrkMuUnvetoFilterToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_CSCTightHaloTrkMuUnvetoFilter"))),
	Flag_CSCTightHalo2015FilterToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_CSCTightHalo2015Filter"))),
	Flag_globalTightHalo2016FilterToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_globalTightHalo2016Filter"))),
	Flag_globalSuperTightHalo2016FilterToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_globalSuperTightHalo2016Filter"))),
	Flag_HcalStripHaloFilterToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_HcalStripHaloFilter"))),
	Flag_hcalLaserEventFilterToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_hcalLaserEventFilter"))),
	Flag_EcalDeadCellTriggerPrimitiveFilterToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_EcalDeadCellTriggerPrimitiveFilter"))),
	Flag_EcalDeadCellBoundaryEnergyFilterToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_EcalDeadCellBoundaryEnergyFilter"))),
	Flag_ecalBadCalibFilterToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_ecalBadCalibFilter"))),
	Flag_goodVerticesToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_goodVertices"))),
	Flag_eeBadScFilterToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_eeBadScFilter"))),
	Flag_ecalLaserCorrFilterToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_ecalLaserCorrFilter"))),
	Flag_trkPOGFiltersToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_trkPOGFilters"))),
	Flag_chargedHadronTrackResolutionFilterToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_chargedHadronTrackResolutionFilter"))),
	Flag_muonBadTrackFilterToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_muonBadTrackFilter"))),
	Flag_BadChargedCandidateFilterToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_BadChargedCandidateFilter"))),
	Flag_BadPFMuonFilterToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_BadPFMuonFilter"))),
	Flag_BadChargedCandidateSummer16FilterToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_BadChargedCandidateSummer16Filter"))),
	Flag_BadPFMuonSummer16FilterToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_BadPFMuonSummer16Filter"))),
	Flag_BadPFMuonDzFilterToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_BadPFMuonDzFilter"))),
	Flag_hfNoisyHitsFilterToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_hfNoisyHitsFilter"))),
	Flag_trkPOG_manystripclus53XToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_trkPOG_manystripclus53X"))),
	Flag_trkPOG_toomanystripclus53XToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_trkPOG_toomanystripclus53X"))),
	Flag_trkPOG_logErrorTooManyClustersToken_(mayConsume<edm::HLTPathStatus>(iConfig.getParameter<edm::InputTag>("Flag_trkPOG_logErrorTooManyClusters")))
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
		//sumScaleWeights = fs->make<TH1D>("sumScaleWeights",";;sumScaleWeights;",9,-0.5,8.5);
		//sumPdfWeights = fs->make<TH1D>("sumPdfWeights",";;sumPdfWeights;",100,-0.5,99.5);
		//sumAlphasWeights = fs->make<TH1D>("sumAlphasWeights",";;sumAlphasWeights;",2,-0.5,1.5);

		sumWeights->Sumw2();
		//sumScaleWeights->Sumw2();
		//sumPdfWeights->Sumw2();
		//sumAlphasWeights->Sumw2();

		// GK, JER setup
		jerRes_ = JME::JetResolution(iConfig.getParameter<edm::FileInPath>("jer_PtResolution").fullPath());
		jerSF_ = JME::JetResolutionScaleFactor(iConfig.getParameter<edm::FileInPath>("jer_ScaleFactor").fullPath());
	}
	else {
		sumWeights = 0;
		//sumScaleWeights = 0;
		//sumPdfWeights = 0;
		//sumAlphasWeights = 0;
	}
	// JEC uncertainty, data and MC txt paths are listed in .py file
	// std::string jecUncPath = iConfig.getParameter<std::string>("jec_Uncertainty"); // for cms.string
	// jecUnc_ = std::make_unique<JetCorrectionUncertainty>(jecUncPath);
	// edm::FileInPath jecUncPath = iConfig.getParameter<edm::FileInPath>("jec_Uncertainty"); // no uncertainties for Run 3 yet
	// JetCorrectorParameters jecUncParams(jecUncPath.fullPath());
    // jecUnc_ = std::make_unique<JetCorrectionUncertainty>(jecUncParams);

	// ----- Get Triggers ----- // 

	if( debug ) cout<<"Getting triggers..."<<endl;

	/* // trigger names as of June 13, run 368822
	HLT_L1SingleLLPJet_v2
	HLT_HT200_L1SingleLLPJet_DisplacedDijet35_Inclusive1PtrkShortSig5_v5
	HLT_HT200_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5_v5
	HLT_HT240_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5_v2 
	HLT_HT280_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5_v2
	HLT_HT170_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack_v5
	HLT_HT200_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack_v5
	HLT_HT270_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack_v5
	HLT_HT200_L1SingleLLPJet_DisplacedDijet60_DisplacedTrack_v5
	HLT_HT320_L1SingleLLPJet_DisplacedDijet60_Inclusive_v5
	HLT_HT420_L1SingleLLPJet_DisplacedDijet60_Inclusive_v5
	HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay0p5nsTrackless_v5 
	HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay1nsInclusive_v5
	HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay1nsTrackless_v5
	HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay2nsInclusive_v5
	*/

// 	triggerPathNames.push_back("HLT_L1SingleLLPJet");
//   //triggerPathNames.push_back("HLT_HT200_L1SingleLLPJet_DisplacedDijet35_Inclusive1PtrkShortSig5");
// 	triggerPathNames.push_back("HLT_HT200_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5");
// 	triggerPathNames.push_back("HLT_HT240_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5");
// 	triggerPathNames.push_back("HLT_HT280_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5");
// 	triggerPathNames.push_back("HLT_HT170_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack");
// 	triggerPathNames.push_back("HLT_HT200_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack");
// 	triggerPathNames.push_back("HLT_HT270_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack");
// 	triggerPathNames.push_back("HLT_HT200_L1SingleLLPJet_DisplacedDijet60_DisplacedTrack");
// 	triggerPathNames.push_back("HLT_HT320_L1SingleLLPJet_DisplacedDijet60_Inclusive");
// 	triggerPathNames.push_back("HLT_HT420_L1SingleLLPJet_DisplacedDijet60_Inclusive");
// 	triggerPathNames.push_back("HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay0p5nsTrackless");
// 	triggerPathNames.push_back("HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay1nsInclusive");
// 	triggerPathNames.push_back("HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay1nsTrackless");
// 	triggerPathNames.push_back("HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay2nsInclusive");

	// Dec 2023 from https://github.com/cms-sw/cmssw/blob/master/HLTrigger/Configuration/python/HLTrigger_Datasets_GRun_cff.py#L1287
	triggerPathNames.push_back("HLT_CaloMET60_DTCluster50");
    triggerPathNames.push_back("HLT_CaloMET60_DTClusterNoMB1S50");
    triggerPathNames.push_back("HLT_HT170_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack");
    triggerPathNames.push_back("HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay0p5nsTrackless");
    triggerPathNames.push_back("HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay1nsInclusive");
    triggerPathNames.push_back("HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay1nsTrackless");
    triggerPathNames.push_back("HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay2nsInclusive");
    triggerPathNames.push_back("HLT_HT200_L1SingleLLPJet_DisplacedDijet35_Inclusive1PtrkShortSig5");
    triggerPathNames.push_back("HLT_HT200_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack");
    triggerPathNames.push_back("HLT_HT200_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5");
    triggerPathNames.push_back("HLT_HT200_L1SingleLLPJet_DisplacedDijet60_DisplacedTrack");
    triggerPathNames.push_back("HLT_HT240_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5");
    triggerPathNames.push_back("HLT_HT270_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack");
    triggerPathNames.push_back("HLT_HT280_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5");
    triggerPathNames.push_back("HLT_HT320_L1SingleLLPJet_DisplacedDijet60_Inclusive");
    triggerPathNames.push_back("HLT_HT320_L1SingleLLPJet_DisplacedDijet60_Inclusive");
    triggerPathNames.push_back("HLT_HT200_L1SingleLLPJet_PFJet60_NeutralHadronFrac0p7");
    triggerPathNames.push_back("HLT_HT200_L1SingleLLPJet_PFJet60_NeutralHadronFrac0p8");
    triggerPathNames.push_back("HLT_HT350_DelayedJet40_SingleDelay3nsInclusive");
    triggerPathNames.push_back("HLT_HT350_DelayedJet40_SingleDelay3p25nsInclusive");
    triggerPathNames.push_back("HLT_HT350_DelayedJet40_SingleDelay3p5nsInclusive");
    triggerPathNames.push_back("HLT_HT350");
    triggerPathNames.push_back("HLT_HT400_DisplacedDijet40_DisplacedTrack");
    triggerPathNames.push_back("HLT_HT420_L1SingleLLPJet_DisplacedDijet60_Inclusive");
    triggerPathNames.push_back("HLT_HT425");
    triggerPathNames.push_back("HLT_HT430_DelayedJet40_DoubleDelay0p5nsInclusive");
    triggerPathNames.push_back("HLT_HT430_DelayedJet40_DoubleDelay0p5nsTrackless");
    triggerPathNames.push_back("HLT_HT430_DelayedJet40_DoubleDelay0p75nsTrackless");
    triggerPathNames.push_back("HLT_HT430_DelayedJet40_DoubleDelay1nsInclusive");
    triggerPathNames.push_back("HLT_HT430_DelayedJet40_DoubleDelay1nsTrackless");
    triggerPathNames.push_back("HLT_HT430_DelayedJet40_DoubleDelay1p25nsInclusive");
    triggerPathNames.push_back("HLT_HT430_DelayedJet40_DoubleDelay1p5nsInclusive");
    triggerPathNames.push_back("HLT_HT430_DelayedJet40_SingleDelay0p5nsInclusive");
    triggerPathNames.push_back("HLT_HT430_DelayedJet40_SingleDelay0p5nsTrackless");
    triggerPathNames.push_back("HLT_HT430_DelayedJet40_SingleDelay1nsInclusive");
    triggerPathNames.push_back("HLT_HT430_DelayedJet40_SingleDelay1nsTrackless");
    triggerPathNames.push_back("HLT_HT430_DelayedJet40_SingleDelay1p25nsTrackless");
    triggerPathNames.push_back("HLT_HT430_DelayedJet40_SingleDelay1p5nsInclusive");
    triggerPathNames.push_back("HLT_HT430_DelayedJet40_SingleDelay1p5nsTrackless");
    triggerPathNames.push_back("HLT_HT430_DelayedJet40_SingleDelay2nsInclusive");
    triggerPathNames.push_back("HLT_HT430_DelayedJet40_SingleDelay2p25nsInclusive");
    triggerPathNames.push_back("HLT_HT430_DelayedJet40_SingleDelay2p5nsInclusive");
    triggerPathNames.push_back("HLT_HT550_DisplacedDijet60_Inclusive");
    triggerPathNames.push_back("HLT_L1MET_DTCluster50");
	triggerPathNames.push_back("HLT_L1MET_DTClusterNoMB1S50");
    triggerPathNames.push_back("HLT_L1Mu6HT240");
    triggerPathNames.push_back("HLT_L1SingleLLPJet");
    triggerPathNames.push_back("HLT_L1Tau_DelayedJet40_DoubleDelay0p5nsTrackless");
    triggerPathNames.push_back("HLT_L1Tau_DelayedJet40_DoubleDelay0p75nsInclusive");
    triggerPathNames.push_back("HLT_L1Tau_DelayedJet40_DoubleDelay1nsTrackless");
    triggerPathNames.push_back("HLT_L1Tau_DelayedJet40_DoubleDelay1p25nsInclusive");
    triggerPathNames.push_back("HLT_L1Tau_DelayedJet40_DoubleDelay1p25nsTrackless");
    triggerPathNames.push_back("HLT_L1Tau_DelayedJet40_DoubleDelay1p5nsInclusive");
    triggerPathNames.push_back("HLT_L1Tau_DelayedJet40_DoubleDelay1p5nsTrackless");
    triggerPathNames.push_back("HLT_L1Tau_DelayedJet40_DoubleDelay1p75nsInclusive");
    triggerPathNames.push_back("HLT_L1Tau_DelayedJet40_SingleDelay2p5nsTrackless");
    triggerPathNames.push_back("HLT_L1Tau_DelayedJet40_SingleDelay2p75nsTrackless");
    triggerPathNames.push_back("HLT_L1Tau_DelayedJet40_SingleDelay3nsTrackless");
    triggerPathNames.push_back("HLT_L1Tau_DelayedJet40_SingleDelay3p5nsInclusive");
    triggerPathNames.push_back("HLT_L1Tau_DelayedJet40_SingleDelay3p75nsInclusive");
    triggerPathNames.push_back("HLT_L1Tau_DelayedJet40_SingleDelay4nsInclusive");
    triggerPathNames.push_back("HLT_Mu6HT240_DisplacedDijet30_Inclusive1PtrkShortSig5_DisplacedLoose");
    triggerPathNames.push_back("HLT_Mu6HT240_DisplacedDijet35_Inclusive0PtrkShortSig5");
    triggerPathNames.push_back("HLT_Mu6HT240_DisplacedDijet35_Inclusive1PtrkShortSig5_DisplacedLoose");
    triggerPathNames.push_back("HLT_Mu6HT240_DisplacedDijet40_Inclusive0PtrkShortSig5");
    triggerPathNames.push_back("HLT_Mu6HT240_DisplacedDijet40_Inclusive1PtrkShortSig5_DisplacedLoose");
    triggerPathNames.push_back("HLT_Mu6HT240_DisplacedDijet45_Inclusive0PtrkShortSig5");
    triggerPathNames.push_back("HLT_Mu6HT240_DisplacedDijet50_Inclusive0PtrkShortSig5");
    triggerPathNames.push_back("HLT_PFJet200_TimeGt2p5ns");
    triggerPathNames.push_back("HLT_PFJet200_TimeLtNeg2p5ns");

        /*
	ifstream myfile (edm::FileInPath(triggerPathNamesFile_.c_str()).fullPath().c_str()) ;

	if( myfile.is_open() ){
		std::string line;
		int index;
		std::string hltpathname;

		while( myfile>>index>>hltpathname ){
			triggerPathNames.push_back( hltpathname );
                        cout<<hltpathname<<endl;
		}
		myfile.close();
	} else {
		std::cout << "ERROR!!! Could not open trigger path name file : " << edm::FileInPath(triggerPathNamesFile_.c_str()).fullPath().c_str() << "\n";
	}	
        */

	if( debug ) cout<<"Defining trigger hists..."<<endl;

	const int NTriggers = triggerPathNames.size();
	NEvents_HLT = fs->make<TH1F>( "NEvents_HLT", " ; ; NEvents_HLT;" , NTriggers, 0, NTriggers );

		bool save_hit = false; 
	if( debug ) cout<<"Defining trigger hists..."<<endl;

	for( int i=0; i<NTriggers; i++ ){
		NEvents_HLT->GetXaxis()->SetBinLabel(i+1,  Form("%s", triggerPathNames.at(i).c_str() ) );
		triggerPathNamesIndices[triggerPathNames.at(i)] = i;
	}

    // Check metFilterBitsToken initialization:
    metFilterBitsToken_isValid = false;
    if (!iConfig.getParameter<edm::InputTag>("metFilterBits").label().empty()) {
        metFilterBitsToken_ = mayConsume<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterBits"));
        metFilterBitsToken_isValid = true;
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
	iEvent.getByToken(rhoFastjetAllToken_,rhoFastjetAll);

	iEvent.getByToken(metToken_, met);
	iEvent.getByToken(metPuppiToken_, metPuppi); 

	// Physics objects
	iEvent.getByToken(electronsToken_, electrons);
	iEvent.getByToken(muonsToken_, muons);
	iEvent.getByToken(tausToken_, taus);
	iEvent.getByToken(photonsToken_, photons);
	iEvent.getByToken(jetsToken_, jets);
	iEvent.getByToken(jetsCorrToken_, jetsCorr); // GK, adding JECs
	iEvent.getByToken(jetsPuppiCorrToken_, jetsPuppiCorr);
	iEvent.getByToken(calojetsToken_, calojets);
	iEvent.getByToken(LRJetsToken_, LRJets);
	iEvent.getByToken(caloLRJetsToken_, caloLRJets);
	iEvent.getByToken(l1jetsToken_, l1jets);

	// Low-level objects
	//iEvent.getByToken(tracksToken_,tracks);
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
	iEvent.getByToken(puInfoToken_,puInfo);
	
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

	output_tree->Branch( "fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll );

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
	output_tree->Branch( "met_PUPPI_Pt", &met_PUPPI_Pt );
	output_tree->Branch( "met_PUPPI_Phi", &met_PUPPI_Phi );
	output_tree->Branch( "met_PUPPI_SumEt", &met_PUPPI_SumEt );

    output_tree->Branch( "Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter);
    output_tree->Branch( "Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter);
    output_tree->Branch( "Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter);
    output_tree->Branch( "Flag_CSCTightHaloTrkMuUnvetoFilter", &Flag_CSCTightHaloTrkMuUnvetoFilter);
    output_tree->Branch( "Flag_CSCTightHalo2015Filter", &Flag_CSCTightHalo2015Filter);
    output_tree->Branch( "Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter);
    output_tree->Branch( "Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter);
    output_tree->Branch( "Flag_HcalStripHaloFilter", &Flag_HcalStripHaloFilter);
    output_tree->Branch( "Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter);
    output_tree->Branch( "Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter);
    output_tree->Branch( "Flag_EcalDeadCellBoundaryEnergyFilter", &Flag_EcalDeadCellBoundaryEnergyFilter);
    output_tree->Branch( "Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter);
    output_tree->Branch( "Flag_goodVertices", &Flag_goodVertices);
    output_tree->Branch( "Flag_eeBadScFilter", &Flag_eeBadScFilter);
    output_tree->Branch( "Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter);
    output_tree->Branch( "Flag_trkPOGFilters", &Flag_trkPOGFilters);
    output_tree->Branch( "Flag_chargedHadronTrackResolutionFilter", &Flag_chargedHadronTrackResolutionFilter);
    output_tree->Branch( "Flag_muonBadTrackFilter", &Flag_muonBadTrackFilter);
    output_tree->Branch( "Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter);
    output_tree->Branch( "Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter);
    output_tree->Branch( "Flag_BadPFMuonDzFilter", &Flag_BadPFMuonDzFilter);
    output_tree->Branch( "Flag_hfNoisyHitsFilter", &Flag_hfNoisyHitsFilter);
    output_tree->Branch( "Flag_BadChargedCandidateSummer16Filter", &Flag_BadChargedCandidateSummer16Filter);
    output_tree->Branch( "Flag_BadPFMuonSummer16Filter", &Flag_BadPFMuonSummer16Filter);
    output_tree->Branch( "Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X);
    output_tree->Branch( "Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X);
    output_tree->Branch( "Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters);
    output_tree->Branch( "Flag_METFilters_2022_2023_PromptReco", &Flag_METFilters_2022_2023_PromptReco);

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
	output_tree->Branch( "ele_EtaSC", &ele_EtaSC);
	output_tree->Branch( "ele_passCutBasedIDLoose", &ele_passCutBasedIDLoose );
	output_tree->Branch( "ele_passCutBasedIDMedium", &ele_passCutBasedIDMedium );
	output_tree->Branch( "ele_passCutBasedIDTight", &ele_passCutBasedIDTight );
	output_tree->Branch( "ele_pileupIso", &ele_pileupIso );
	output_tree->Branch( "ele_chargedIso", &ele_chargedIso );
	output_tree->Branch( "ele_photonIso", &ele_photonIso );
	output_tree->Branch( "ele_neutralHadIso", &ele_neutralHadIso );
	output_tree->Branch( "ele_EcalRechitIDs", &ele_EcalRechitIDs );
	output_tree->Branch( "ele_EcalRechitIndices", &ele_EcalRechitIndices );
	output_tree->Branch( "ele_SeedRechitID", &ele_SeedRechitID );
	output_tree->Branch( "ele_SeedRechitIndex", &ele_SeedRechitIndex );
	//output_tree->Branch( "ele_passHLTFilter", &ele_passHLTFilter );

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
	output_tree->Branch( "muon_pileupIso", &muon_pileupIso );
	output_tree->Branch( "muon_chargedIso", &muon_chargedIso );
	output_tree->Branch( "muon_photonIso", &muon_photonIso );
	output_tree->Branch( "muon_neutralHadIso", &muon_neutralHadIso );
	output_tree->Branch( "muon_ip3dSignificance", &muon_ip3dSignificance );
	output_tree->Branch( "muon_dB", &muon_dB );
	output_tree->Branch( "muon_edB", &muon_edB );
	//output_tree->Branch( "muon_passHLTFilter", &muon_passHLTFilter );

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

	//output_tree->Branch( "pho_passHLTFilter", &pho_passHLTFilter );

};

// ------------------------------------------------------------------------------------
void DisplacedHcalJetNTuplizer::EnableJetBranches(){

	// AK4 PF Jets
	output_tree->Branch( "n_jet", &n_jet );
	output_tree->Branch( "jetRaw_Pt", &jetRaw_Pt );
	output_tree->Branch( "jetRaw_E", &jetRaw_E );
	output_tree->Branch( "jet_Pt", &jet_Pt );
	output_tree->Branch( "jet_Eta", &jet_Eta );
	output_tree->Branch( "jet_Phi", &jet_Phi );
	output_tree->Branch( "jet_E", &jet_E );
	output_tree->Branch( "jet_Mass", &jet_Mass );
	output_tree->Branch( "jet_JetArea", &jet_JetArea );
	output_tree->Branch( "jet_Pt_noJER", &jet_Pt_noJER );
	output_tree->Branch( "jet_E_noJER", &jet_E_noJER );
	output_tree->Branch( "jet_Mass_noJER", &jet_Mass_noJER );
	output_tree->Branch( "jet_Pt_JER_up", &jet_Pt_JER_up );
	output_tree->Branch( "jet_E_JER_up", &jet_E_JER_up );
	output_tree->Branch( "jet_Mass_JER_up", &jet_Mass_JER_up );
	output_tree->Branch( "jet_Pt_JER_down", &jet_Pt_JER_down );
	output_tree->Branch( "jet_E_JER_down", &jet_E_JER_down );
	output_tree->Branch( "jet_Mass_JER_down", &jet_Mass_JER_down );
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
	output_tree->Branch( "jet_PileupE", &jet_PileupE );
	output_tree->Branch( "jet_PileupId", &jet_PileupId );
	output_tree->Branch( "jet_PileupIdFlag", &jet_PileupIdFlag );
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
    output_tree->Branch( "jet_DeepCSV_prob_b", &jet_DeepCSV_prob_b ); 
    output_tree->Branch( "jet_DeepCSV_prob_c", &jet_DeepCSV_prob_c ); 
    output_tree->Branch( "jet_DeepCSV_prob_bb", &jet_DeepCSV_prob_bb ); 
    output_tree->Branch( "jet_DeepCSV_prob_udsg", &jet_DeepCSV_prob_udsg ); 
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

        output_tree->Branch( "n_BunchXing", &n_BunchXing );
        output_tree->Branch( "BunchXing", &BunchXing );
	output_tree->Branch( "nPU", &nPU );
	output_tree->Branch( "nPUmean", &nPUmean );

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

	fixedGridRhoFastjetAll 	= -9999.9;

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
	met_PUPPI_Pt = -9999.9;
	met_PUPPI_Phi = -9999.9;
	met_PUPPI_SumEt = -9999.9;

	Flag_HBHENoiseFilter = false;
	Flag_HBHENoiseIsoFilter = false;
	Flag_CSCTightHaloFilter = false;
	Flag_CSCTightHaloTrkMuUnvetoFilter = false;
	Flag_CSCTightHalo2015Filter = false;
	Flag_globalTightHalo2016Filter = false;
	Flag_globalSuperTightHalo2016Filter = false;
	Flag_HcalStripHaloFilter = false;
	Flag_hcalLaserEventFilter = false;
	Flag_EcalDeadCellTriggerPrimitiveFilter = false;
	Flag_EcalDeadCellBoundaryEnergyFilter = false;
	Flag_ecalBadCalibFilter = false;
	Flag_goodVertices = false;
	//Flag_trackingFailureFilter = false;
	Flag_eeBadScFilter = false;
	Flag_ecalLaserCorrFilter = false;
	Flag_trkPOGFilters = false;
	Flag_chargedHadronTrackResolutionFilter = false;
	Flag_muonBadTrackFilter = false;
	Flag_BadChargedCandidateFilter = false;
	Flag_BadPFMuonFilter = false;
	Flag_BadChargedCandidateSummer16Filter = false;
	Flag_BadPFMuonSummer16Filter = false;
	Flag_BadPFMuonDzFilter = false;
	Flag_hfNoisyHitsFilter = false;
	Flag_trkPOG_manystripclus53X = false;
	Flag_trkPOG_toomanystripclus53X = false;
	Flag_trkPOG_logErrorTooManyClusters = false;
	Flag_METFilters_2022_2023_PromptReco = false;

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
	ele_EtaSC.clear();
	ele_passCutBasedIDLoose.clear();
	ele_passCutBasedIDMedium.clear();
	ele_passCutBasedIDTight.clear();
	ele_pileupIso.clear();
	ele_chargedIso.clear();
	ele_photonIso.clear();
	ele_neutralHadIso.clear();
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
	muon_pileupIso.clear();
	muon_chargedIso.clear();
	muon_photonIso.clear();
	muon_neutralHadIso.clear();
	muon_ip3dSignificance.clear();
	muon_dB.clear();
	muon_edB.clear();
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
	jetRaw_Pt.clear();
	jetRaw_E.clear();
	jet_Pt_noJER.clear();
	jet_E_noJER.clear();
	jet_Mass_noJER.clear();
	jet_Pt_JER_up.clear();
	jet_E_JER_up.clear();
	jet_Mass_JER_up.clear();
	jet_Pt_JER_down.clear();
	jet_E_JER_down.clear();
	jet_Mass_JER_down.clear();
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
	jet_PileupE.clear();
	jet_PileupId.clear();
	jet_PileupIdFlag.clear();
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
    jet_DeepCSV_prob_b.clear();
    jet_DeepCSV_prob_c.clear();
    jet_DeepCSV_prob_bb.clear();
    jet_DeepCSV_prob_udsg.clear();
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

        n_BunchXing = 0;
        BunchXing.clear();
	nPU.clear();
	nPUmean.clear();

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
	//FillTauBranches( iEvent );
	FillPhotonBranches( iEvent );
	FillJetBranches( iEvent, iSetup );

	// Low-Level Objects
	FillTrackBranches( iEvent ); //, iSetup );
	//FillPFCandidateBranches( iEvent, iSetup );
	//FillSecondaryVerticesBranches( iEvent, iSetup );

	// Hits
	FillEcalRechitBranches( iEvent, iSetup ); 
	FillHcalRechitBranches( iEvent, iSetup );

	// MC
	if( !isData ){
		//FillPileupBranches();
		//FillMCBranches( iEvent, iSetup );
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

	// get rho
	fixedGridRhoFastjetAll = *rhoFastjetAll;

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

	string hltPathNameReq = "HLT_";

	//fill trigger information
	const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

	if( debug ) cout<<"N TRIGGER BITS = "<<triggerBits->size()<<endl;

	for( auto triggerPathName: triggerPathNames ){

		if( debug ) cout<<"  -- "<<triggerPathName<<endl;

		if( triggerPathName == "" ) continue;

		bool found_trigger = false;

		for( uint i = 0; i < triggerBits->size(); ++i) { // TODO: CHECK OUT INDEXING

			if( (names.triggerName(i)).find(hltPathNameReq) == string::npos ) continue;

			int lastUnderscorePos = (names.triggerName(i)).find_last_of("_");
			//if( lastUnderscorePos == string::npos ) continue;

			string hltPathNameWithoutVersionNumber = (names.triggerName(i)).substr(0,lastUnderscorePos);
			if( hltPathNameWithoutVersionNumber != triggerPathName ) continue;

			if( debug && triggerBits->accept(i) ) cout<<"PASS TRIGGER "<<names.triggerName(i)<<endl;

			// Fill Branches //

			HLT_Decision.push_back( triggerBits->accept(i) );
			HLT_Prescale.push_back( -1 );			
                        // if (isData_) HLT_Prescale.push_back( triggerPrescales->getPrescaleForIndex(i) ); // FIX
                        // else HLT_Prescale.push_back( 1 ); // TODO Need to figure out yields

			if( triggerBits->accept(i) ) NEvents_HLT->Fill(triggerPathNamesIndices[triggerPathName]); // FIX WEIGHTS

			found_trigger = true;

			break;
	
		}
		
		if( !found_trigger ){
			if( debug ) cout<<"    --> trig not found in triggerBits"<<endl;
			HLT_Decision.push_back( false );
			HLT_Prescale.push_back( -2 );
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
	auto &met_CHS_temp = met->front(); // CHS MET
	met_Pt = met_CHS_temp.pt();
	met_Phi = met_CHS_temp.phi();
	met_SumEt = met_CHS_temp.sumEt();
	auto &met_temp = metPuppi->front();
	met_PUPPI_Pt = met_temp.pt();
	met_PUPPI_Phi = met_temp.phi();
	met_PUPPI_SumEt = met_temp.sumEt();

	if( metFilterBitsToken_isValid ){ 
		iEvent.getByToken(metFilterBitsToken_, metFilterBits);
	
    	const edm::TriggerNames &metNames = iEvent.triggerNames(*metFilterBits);

	    for(unsigned int i = 0, n = metFilterBits->size(); i < n; ++i){
	        if( debug && metNames.triggerName(i).find("Flag_") != std::string::npos ){ cout<<metNames.triggerName(i).c_str()<<"  "<<metFilterBits->accept(i)<<endl; }
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_HBHENoiseFilter") == 0)                    Flag_HBHENoiseFilter = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_HBHENoiseIsoFilter") == 0)                 Flag_HBHENoiseIsoFilter = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_CSCTightHaloFilter") == 0)                 Flag_CSCTightHaloFilter = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_CSCTightHaloTrkMuUnvetoFilter") == 0)      Flag_CSCTightHaloTrkMuUnvetoFilter = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_CSCTightHalo2015Filter") == 0)             Flag_CSCTightHalo2015Filter = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_globalTightHalo2016Filter") == 0)          Flag_globalTightHalo2016Filter = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_globalSuperTightHalo2016Filter") == 0)     Flag_globalSuperTightHalo2016Filter = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_HcalStripHaloFilter") == 0)                Flag_HcalStripHaloFilter = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_hcalLaserEventFilter") == 0)               Flag_hcalLaserEventFilter = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_EcalDeadCellTriggerPrimitiveFilter") == 0) Flag_EcalDeadCellTriggerPrimitiveFilter = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_EcalDeadCellBoundaryEnergyFilter") == 0)   Flag_EcalDeadCellBoundaryEnergyFilter = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_ecalBadCalibFilter") == 0)                 Flag_ecalBadCalibFilter = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_goodVertices") == 0)                       Flag_goodVertices = metFilterBits->accept(i);
	        //if(strcmp(metNames.triggerName(i).c_str(), "Flag_trackingFailureFilter") == 0) Flag_trackingFailureFilter = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_eeBadScFilter") == 0)                      Flag_eeBadScFilter = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_ecalLaserCorrFilter") == 0)                Flag_ecalLaserCorrFilter = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_trkPOGFilters") == 0)                      Flag_trkPOGFilters = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_chargedHadronTrackResolutionFilter") == 0) Flag_chargedHadronTrackResolutionFilter = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_muonBadTrackFilter") == 0)                 Flag_muonBadTrackFilter = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_BadChargedCandidateFilter") == 0)          Flag_BadChargedCandidateFilter = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_BadPFMuonFilter") == 0)                    Flag_BadPFMuonFilter = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_BadChargedCandidateSummer16Filter") == 0)  Flag_BadChargedCandidateSummer16Filter = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_BadPFMuonSummer16Filter") == 0)            Flag_BadPFMuonSummer16Filter = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_BadPFMuonDzFilter") == 0)                  Flag_BadPFMuonDzFilter = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_hfNoisyHitsFilter") == 0)                  Flag_hfNoisyHitsFilter = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_trkPOG_manystripclus53X") == 0)            Flag_trkPOG_manystripclus53X = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_trkPOG_toomanystripclus53X") == 0)         Flag_trkPOG_toomanystripclus53X = metFilterBits->accept(i);
	        if(strcmp(metNames.triggerName(i).c_str(), "Flag_trkPOG_logErrorTooManyClusters") == 0)     Flag_trkPOG_logErrorTooManyClusters = metFilterBits->accept(i);
	    }

	} else {
		edm::Handle<edm::HLTPathStatus> Flag_HBHENoiseFilterHandle;
		edm::Handle<edm::HLTPathStatus> Flag_HBHENoiseIsoFilterHandle;
		edm::Handle<edm::HLTPathStatus> Flag_CSCTightHaloFilterHandle;
		edm::Handle<edm::HLTPathStatus> Flag_CSCTightHaloTrkMuUnvetoFilterHandle;
		edm::Handle<edm::HLTPathStatus> Flag_CSCTightHalo2015FilterHandle;
		edm::Handle<edm::HLTPathStatus> Flag_globalTightHalo2016FilterHandle;
		edm::Handle<edm::HLTPathStatus> Flag_globalSuperTightHalo2016FilterHandle;
		edm::Handle<edm::HLTPathStatus> Flag_HcalStripHaloFilterHandle;
		edm::Handle<edm::HLTPathStatus> Flag_hcalLaserEventFilterHandle;
		edm::Handle<edm::HLTPathStatus> Flag_EcalDeadCellTriggerPrimitiveFilterHandle;
		edm::Handle<edm::HLTPathStatus> Flag_EcalDeadCellBoundaryEnergyFilterHandle;
		edm::Handle<edm::HLTPathStatus> Flag_ecalBadCalibFilterHandle;
		edm::Handle<edm::HLTPathStatus> Flag_goodVerticesHandle;
		edm::Handle<edm::HLTPathStatus> Flag_eeBadScFilterHandle;
		edm::Handle<edm::HLTPathStatus> Flag_ecalLaserCorrFilterHandle;
		edm::Handle<edm::HLTPathStatus> Flag_trkPOGFiltersHandle;
		edm::Handle<edm::HLTPathStatus> Flag_chargedHadronTrackResolutionFilterHandle;
		edm::Handle<edm::HLTPathStatus> Flag_muonBadTrackFilterHandle;
		edm::Handle<edm::HLTPathStatus> Flag_BadChargedCandidateFilterHandle;
		edm::Handle<edm::HLTPathStatus> Flag_BadPFMuonFilterHandle;
		edm::Handle<edm::HLTPathStatus> Flag_BadChargedCandidateSummer16FilterHandle;
		edm::Handle<edm::HLTPathStatus> Flag_BadPFMuonSummer16FilterHandle;
		edm::Handle<edm::HLTPathStatus> Flag_BadPFMuonDzFilterHandle;
		edm::Handle<edm::HLTPathStatus> Flag_hfNoisyHitsFilterHandle;
		edm::Handle<edm::HLTPathStatus> Flag_trkPOG_manystripclus53XHandle;
		edm::Handle<edm::HLTPathStatus> Flag_trkPOG_toomanystripclus53XHandle;
		edm::Handle<edm::HLTPathStatus> Flag_trkPOG_logErrorTooManyClustersHandle;

		iEvent.getByToken(Flag_HBHENoiseFilterToken_, Flag_HBHENoiseFilterHandle);
		iEvent.getByToken(Flag_HBHENoiseIsoFilterToken_, Flag_HBHENoiseIsoFilterHandle);
		iEvent.getByToken(Flag_CSCTightHaloFilterToken_, Flag_CSCTightHaloFilterHandle);
		iEvent.getByToken(Flag_CSCTightHaloTrkMuUnvetoFilterToken_, Flag_CSCTightHaloTrkMuUnvetoFilterHandle);
		iEvent.getByToken(Flag_CSCTightHalo2015FilterToken_, Flag_CSCTightHalo2015FilterHandle);
		iEvent.getByToken(Flag_globalTightHalo2016FilterToken_, Flag_globalTightHalo2016FilterHandle);
		iEvent.getByToken(Flag_globalSuperTightHalo2016FilterToken_, Flag_globalSuperTightHalo2016FilterHandle);
		iEvent.getByToken(Flag_HcalStripHaloFilterToken_, Flag_HcalStripHaloFilterHandle);
		iEvent.getByToken(Flag_hcalLaserEventFilterToken_, Flag_hcalLaserEventFilterHandle);
		iEvent.getByToken(Flag_EcalDeadCellTriggerPrimitiveFilterToken_, Flag_EcalDeadCellTriggerPrimitiveFilterHandle);
		iEvent.getByToken(Flag_EcalDeadCellBoundaryEnergyFilterToken_, Flag_EcalDeadCellBoundaryEnergyFilterHandle);
		iEvent.getByToken(Flag_ecalBadCalibFilterToken_, Flag_ecalBadCalibFilterHandle);
		iEvent.getByToken(Flag_goodVerticesToken_, Flag_goodVerticesHandle);
		iEvent.getByToken(Flag_eeBadScFilterToken_, Flag_eeBadScFilterHandle);
		iEvent.getByToken(Flag_ecalLaserCorrFilterToken_, Flag_ecalLaserCorrFilterHandle);
		iEvent.getByToken(Flag_trkPOGFiltersToken_, Flag_trkPOGFiltersHandle);
		iEvent.getByToken(Flag_chargedHadronTrackResolutionFilterToken_, Flag_chargedHadronTrackResolutionFilterHandle);
		iEvent.getByToken(Flag_muonBadTrackFilterToken_, Flag_muonBadTrackFilterHandle);
		iEvent.getByToken(Flag_BadChargedCandidateFilterToken_, Flag_BadChargedCandidateFilterHandle);
		iEvent.getByToken(Flag_BadPFMuonFilterToken_, Flag_BadPFMuonFilterHandle);
		iEvent.getByToken(Flag_BadChargedCandidateSummer16FilterToken_, Flag_BadChargedCandidateSummer16FilterHandle);
		iEvent.getByToken(Flag_BadPFMuonSummer16FilterToken_, Flag_BadPFMuonSummer16FilterHandle);
		iEvent.getByToken(Flag_BadPFMuonDzFilterToken_, Flag_BadPFMuonDzFilterHandle);
		iEvent.getByToken(Flag_hfNoisyHitsFilterToken_, Flag_hfNoisyHitsFilterHandle);
		iEvent.getByToken(Flag_trkPOG_manystripclus53XToken_, Flag_trkPOG_manystripclus53XHandle);
		iEvent.getByToken(Flag_trkPOG_toomanystripclus53XToken_, Flag_trkPOG_toomanystripclus53XHandle);
		iEvent.getByToken(Flag_trkPOG_logErrorTooManyClustersToken_, Flag_trkPOG_logErrorTooManyClustersHandle);

		Flag_HBHENoiseFilter = Flag_HBHENoiseFilterHandle->wasrun() && Flag_HBHENoiseFilterHandle->accept();
		Flag_HBHENoiseIsoFilter = Flag_HBHENoiseIsoFilterHandle->wasrun() && Flag_HBHENoiseIsoFilterHandle->accept();
		Flag_CSCTightHaloFilter = Flag_CSCTightHaloFilterHandle->wasrun() && Flag_CSCTightHaloFilterHandle->accept();
		Flag_CSCTightHaloTrkMuUnvetoFilter = Flag_CSCTightHaloTrkMuUnvetoFilterHandle->wasrun() && Flag_CSCTightHaloTrkMuUnvetoFilterHandle->accept();
		Flag_CSCTightHalo2015Filter = Flag_CSCTightHalo2015FilterHandle->wasrun() && Flag_CSCTightHalo2015FilterHandle->accept();
		Flag_globalTightHalo2016Filter = Flag_globalTightHalo2016FilterHandle->wasrun() && Flag_globalTightHalo2016FilterHandle->accept();
		Flag_globalSuperTightHalo2016Filter = Flag_globalSuperTightHalo2016FilterHandle->wasrun() && Flag_globalSuperTightHalo2016FilterHandle->accept();
		Flag_HcalStripHaloFilter = Flag_HcalStripHaloFilterHandle->wasrun() && Flag_HcalStripHaloFilterHandle->accept();
		Flag_hcalLaserEventFilter = Flag_hcalLaserEventFilterHandle->wasrun() && Flag_hcalLaserEventFilterHandle->accept();
		Flag_EcalDeadCellTriggerPrimitiveFilter = Flag_EcalDeadCellTriggerPrimitiveFilterHandle->wasrun() && Flag_EcalDeadCellTriggerPrimitiveFilterHandle->accept();
		Flag_EcalDeadCellBoundaryEnergyFilter = Flag_EcalDeadCellBoundaryEnergyFilterHandle->wasrun() && Flag_EcalDeadCellBoundaryEnergyFilterHandle->accept();
		Flag_ecalBadCalibFilter = Flag_ecalBadCalibFilterHandle->wasrun() && Flag_ecalBadCalibFilterHandle->accept();
		Flag_goodVertices = Flag_goodVerticesHandle->wasrun() && Flag_goodVerticesHandle->accept();
		Flag_eeBadScFilter = Flag_eeBadScFilterHandle->wasrun() && Flag_eeBadScFilterHandle->accept();
		Flag_ecalLaserCorrFilter = Flag_ecalLaserCorrFilterHandle->wasrun() && Flag_ecalLaserCorrFilterHandle->accept();
		Flag_trkPOGFilters = Flag_trkPOGFiltersHandle->wasrun() && Flag_trkPOGFiltersHandle->accept();
		Flag_chargedHadronTrackResolutionFilter = Flag_chargedHadronTrackResolutionFilterHandle->wasrun() && Flag_chargedHadronTrackResolutionFilterHandle->accept();
		Flag_muonBadTrackFilter = Flag_muonBadTrackFilterHandle->wasrun() && Flag_muonBadTrackFilterHandle->accept();
		Flag_BadChargedCandidateFilter = Flag_BadChargedCandidateFilterHandle->wasrun() && Flag_BadChargedCandidateFilterHandle->accept();
		Flag_BadPFMuonFilter = Flag_BadPFMuonFilterHandle->wasrun() && Flag_BadPFMuonFilterHandle->accept();
		Flag_BadChargedCandidateSummer16Filter = Flag_BadChargedCandidateSummer16FilterHandle->wasrun() && Flag_BadChargedCandidateSummer16FilterHandle->accept();
		Flag_BadPFMuonSummer16Filter = Flag_BadPFMuonSummer16FilterHandle->wasrun() && Flag_BadPFMuonSummer16FilterHandle->accept();
		Flag_BadPFMuonDzFilter = Flag_BadPFMuonDzFilterHandle->wasrun() && Flag_BadPFMuonDzFilterHandle->accept();
		Flag_hfNoisyHitsFilter = Flag_hfNoisyHitsFilterHandle->wasrun() && Flag_hfNoisyHitsFilterHandle->accept();
		Flag_trkPOG_manystripclus53X = Flag_trkPOG_manystripclus53XHandle->wasrun() && Flag_trkPOG_manystripclus53XHandle->accept();
		Flag_trkPOG_toomanystripclus53X = Flag_trkPOG_toomanystripclus53XHandle->wasrun() && Flag_trkPOG_toomanystripclus53XHandle->accept();
		Flag_trkPOG_logErrorTooManyClusters = Flag_trkPOG_logErrorTooManyClustersHandle->wasrun() && Flag_trkPOG_logErrorTooManyClustersHandle->accept();

		if( debug ) { 
			cout<<"Flag_HBHENoiseFilter  "<<Flag_HBHENoiseFilter<<endl;
			cout<<"Flag_HBHENoiseIsoFilter  "<<Flag_HBHENoiseIsoFilter<<endl;
			cout<<"Flag_CSCTightHaloFilter  "<<Flag_CSCTightHaloFilter<<endl;
			cout<<"Flag_CSCTightHaloTrkMuUnvetoFilter  "<<Flag_CSCTightHaloTrkMuUnvetoFilter<<endl;
			cout<<"Flag_CSCTightHalo2015Filter  "<<Flag_CSCTightHalo2015Filter<<endl;
			cout<<"Flag_globalTightHalo2016Filter  "<<Flag_globalTightHalo2016Filter<<endl;
			cout<<"Flag_globalSuperTightHalo2016Filter  "<<Flag_globalSuperTightHalo2016Filter<<endl;
			cout<<"Flag_HcalStripHaloFilter  "<<Flag_HcalStripHaloFilter<<endl;
			cout<<"Flag_hcalLaserEventFilter  "<<Flag_hcalLaserEventFilter<<endl;
			cout<<"Flag_EcalDeadCellTriggerPrimitiveFilter  "<<Flag_EcalDeadCellTriggerPrimitiveFilter<<endl;
			cout<<"Flag_EcalDeadCellBoundaryEnergyFilter  "<<Flag_EcalDeadCellBoundaryEnergyFilter<<endl;
			cout<<"Flag_ecalBadCalibFilter  "<<Flag_ecalBadCalibFilter<<endl;
			cout<<"Flag_goodVertices  "<<Flag_goodVertices<<endl;
			cout<<"Flag_eeBadScFilter  "<<Flag_eeBadScFilter<<endl;
			cout<<"Flag_ecalLaserCorrFilter  "<<Flag_ecalLaserCorrFilter<<endl;
			cout<<"Flag_trkPOGFilters  "<<Flag_trkPOGFilters<<endl;
			cout<<"Flag_chargedHadronTrackResolutionFilter  "<<Flag_chargedHadronTrackResolutionFilter<<endl;
			cout<<"Flag_muonBadTrackFilter  "<<Flag_muonBadTrackFilter<<endl;
			cout<<"Flag_BadChargedCandidateFilter  "<<Flag_BadChargedCandidateFilter<<endl;
			cout<<"Flag_BadPFMuonFilter  "<<Flag_BadPFMuonFilter<<endl;
			cout<<"Flag_BadChargedCandidateSummer16Filter  "<<Flag_BadChargedCandidateSummer16Filter<<endl;
			cout<<"Flag_BadPFMuonSummer16Filter  "<<Flag_BadPFMuonSummer16Filter<<endl;
			cout<<"Flag_BadPFMuonDzFilter  "<<Flag_BadPFMuonDzFilter<<endl;
			cout<<"Flag_hfNoisyHitsFilter  "<<Flag_hfNoisyHitsFilter<<endl;
			cout<<"Flag_trkPOG_manystripclus53X  "<<Flag_trkPOG_manystripclus53X<<endl;
			cout<<"Flag_trkPOG_toomanystripclus53X  "<<Flag_trkPOG_toomanystripclus53X<<endl;
			cout<<"Flag_trkPOG_logErrorTooManyClusters  "<<Flag_trkPOG_logErrorTooManyClusters<<endl;
		}

	}
   
    // Recommendations for Run3 (2022 and 2023 Prompt Reco)
    Flag_METFilters_2022_2023_PromptReco = Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_BadPFMuonDzFilter && Flag_hfNoisyHitsFilter && Flag_eeBadScFilter && Flag_ecalBadCalibFilter;

	// metUncorrectedPt = Met.uncorPt();
	// metUncorrectedPhi = Met.uncorPhi();
	//metType1Pt = Met.pt();
	//metType1Phi = Met.phi();

	/*
	const edm::TriggerNames &metNames = iEvent.triggerNames(*metFilterBits);


	iEvent.getByToken(globalSuperTightHalo2016FilterToken_, globalSuperTightHalo2016Filter);
	iEvent.getByToken(globalTightHalo2016FilterToken_, globalTightHalo2016Filter);
	iEvent.getByToken(BadChargedCandidateFilterToken_, BadChargedCandidateFilter);
	iEvent.getByToken(BadPFMuonFilterToken_, BadPFMuonFilter);
	//iEvent.getByToken(EcalDeadCellTriggerPrimitiveFilterToken_, EcalDeadCellTriggerPrimitiveFilter);
	iEvent.getByToken(HBHENoiseFilterToken_, HBHENoiseFilter);
	iEvent.getByToken(HBHEIsoNoiseFilterToken_, HBHEIsoNoiseFilter);
	// iEvent.getByToken(ecalBadCalibFilterToken_, ecalBadCalibReducedMINIAODFilter);
	// iEvent.getByToken(eeBadScFilterToken_, eeBadScFilter);
	// iEvent.getByToken(primaryVertexFilterToken_, primaryVertexFilter);


	Flag2_globalSuperTightHalo2016Filter = *globalSuperTightHalo2016Filter;
	Flag2_globalTightHalo2016Filter = *globalTightHalo2016Filter;
	Flag2_BadChargedCandidateFilter = *BadChargedCandidateFilter;
	Flag2_BadPFMuonFilter = *BadPFMuonFilter;
	//Flag2_EcalDeadCellTriggerPrimitiveFilter = *EcalDeadCellTriggerPrimitiveFilter;
	Flag2_HBHENoiseFilter = *HBHENoiseFilter;
	Flag2_HBHEIsoNoiseFilter = *HBHEIsoNoiseFilter;
	// Flag2_ecalBadCalibFilter = *ecalBadCalibReducedMINIAODFilter;
	// Flag2_eeBadScFilter = *eeBadScFilter;
	// Flag2_goodVertices = *primaryVertexFilter;


	// *******************************************************************************
	//For Debug printout
	// *******************************************************************************
	// for (unsigned int i = 0, n = metFilterBits->size(); i < n; ++i) {
	//  std::cout << "MET Filter " << metNames.triggerName(i).c_str() << "\n";
	// }

	for(unsigned int i = 0, n = metFilterBits->size(); i < n; ++i){
		if(strcmp(metNames.triggerName(i).c_str(), "Flag_trackingFailureFilter") == 0)
			Flag_trackingFailureFilter = metFilterBits->accept(i);
		else if(strcmp(metNames.triggerName(i).c_str(), "Flag_goodVertices") == 0)
			Flag_goodVertices = metFilterBits->accept(i);
		else if(strcmp(metNames.triggerName(i).c_str(), "Flag_globalSuperTightHalo2016Filter") == 0)
			Flag_globalSuperTightHalo2016Filter = metFilterBits->accept(i);
		else if(strcmp(metNames.triggerName(i).c_str(), "Flag_globalTightHalo2016Filter") == 0)
			Flag_globalTightHalo2016Filter = metFilterBits->accept(i);
		else if(strcmp(metNames.triggerName(i).c_str(), "Flag_trkPOGFilters") == 0)
			Flag_trkPOGFilters = metFilterBits->accept(i);
		else if(strcmp(metNames.triggerName(i).c_str(), "Flag_trkPOG_logErrorTooManyClusters") == 0)
			Flag_trkPOG_logErrorTooManyClusters = metFilterBits->accept(i);
		else if(strcmp(metNames.triggerName(i).c_str(), "Flag_EcalDeadCellTriggerPrimitiveFilter") == 0)
			Flag_EcalDeadCellTriggerPrimitiveFilter = metFilterBits->accept(i);
		else if(strcmp(metNames.triggerName(i).c_str(), "Flag_EcalDeadCellBoundaryEnergyFilter") == 0)
			Flag_EcalDeadCellBoundaryEnergyFilter = metFilterBits->accept(i);
		else if(strcmp(metNames.triggerName(i).c_str(), "Flag_ecalLaserCorrFilter") == 0)
			Flag_ecalLaserCorrFilter = metFilterBits->accept(i);
		else if(strcmp(metNames.triggerName(i).c_str(), "Flag_trkPOG_manystripclus53X") == 0)
			Flag_trkPOG_manystripclus53X = metFilterBits->accept(i);
		else if(strcmp(metNames.triggerName(i).c_str(), "Flag_eeBadScFilter") == 0)
			Flag_eeBadScFilter = metFilterBits->accept(i);
		else if(strcmp(metNames.triggerName(i).c_str(), "Flag_METFilters") == 0)
			Flag_METFilters = metFilterBits->accept(i);
		else if(strcmp(metNames.triggerName(i).c_str(), "Flag_HBHENoiseFilter") == 0)
			Flag_HBHENoiseFilter = metFilterBits->accept(i);
		else if(strcmp(metNames.triggerName(i).c_str(), "Flag_HBHENoiseIsoFilter") == 0)
			Flag_HBHEIsoNoiseFilter = metFilterBits->accept(i);
		else if(strcmp(metNames.triggerName(i).c_str(), "Flag_trkPOG_toomanystripclus53X") == 0)
			Flag_trkPOG_toomanystripclus53X = metFilterBits->accept(i);
		else if(strcmp(metNames.triggerName(i).c_str(), "Flag_hcalLaserEventFilter") == 0)
			Flag_hcalLaserEventFilter = metFilterBits->accept(i);
		else if(strcmp(metNames.triggerName(i).c_str(), "Flag_BadPFMuonFilter") == 0)
			Flag_BadPFMuonFilter = metFilterBits->accept(i);
		else if(strcmp(metNames.triggerName(i).c_str(), "Flag_BadChargedCandidateFilter") == 0)
			Flag_BadChargedCandidateFilter = metFilterBits->accept(i);
		else if(strcmp(metNames.triggerName(i).c_str(), "Flag_ecalBadCalibFilter") == 0)
			Flag_ecalBadCalibFilter = metFilterBits->accept(i);
		else if(strcmp(metNames.triggerName(i).c_str(), "Flag_badMuons") == 0) {
			Flag_badGlobalMuonFilter = metFilterBits->accept(i);
			//cout << "found bad muon flag : " << "\n";
		}
		else if(strcmp(metNames.triggerName(i).c_str(), "Flag_duplicateMuons") == 0)
			Flag_duplicateMuonFilter = metFilterBits->accept(i);
	} //loop over met filters

	*/

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
		ele_EtaSC.push_back( ele.superCluster()->eta() );

		// ID //

		ele_passCutBasedIDLoose.push_back( (*electron_cutbasedID_decisions_loose)[eleRef] );
		ele_passCutBasedIDMedium.push_back( (*electron_cutbasedID_decisions_medium)[eleRef] );
		ele_passCutBasedIDTight.push_back( (*electron_cutbasedID_decisions_tight)[eleRef] );

		// Isolation //
		ele_pileupIso.push_back( ele.pfIsolationVariables().sumPUPt );
    	ele_chargedIso.push_back( ele.pfIsolationVariables().sumChargedHadronPt );
    	ele_photonIso.push_back( ele.pfIsolationVariables().sumPhotonEt );
    	ele_neutralHadIso.push_back( ele.pfIsolationVariables().sumNeutralHadronEt );

		// Rechits Association //

		ele_SeedRechitID.push_back( 0 ); // FIXURGENT
		ele_SeedRechitIndex.push_back( 0 ); // FIXURGENT
		ele_EcalRechitIDs.push_back( vector<uint>{0} ); // FIXURGENT
		ele_EcalRechitIndices.push_back( vector<uint>{0} ); // FIXURGENT

		// Other //

		// FIX ele_HoverE.push_back( ele.hcalOverEcal() );

		/*
		// Conversion Veto		
		if( beamSpot.isValid() && conversions.isValid() ){
			ele_PassConvVeto.push_back( !ConversionTools::hasMatchedConversion(ele,(*conversions),beamSpot->position()) );
		} else {
			ele_PassConvVeto.push_back( false );
			cout << "\n\nERROR!!! conversions not found!!!\n";
		}

		// 1/E - 1/P
		if( ele.ecalEnergy() == 0 ){
			ele_OneOverEminusOneOverP.push_back( 1e30 );
		} else if( !std::isfinite(ele.ecalEnergy())){
			ele_OneOverEminusOneOverP.push_back( 1e30 );
		} else {
			ele_OneOverEminusOneOverP.push_back( 1./ele.ecalEnergy()  -  ele.eSuperClusterOverP()/ele.ecalEnergy() );
		}
		*/

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

		// Isolation //
		muon_pileupIso.push_back( muon.pfIsolationR04().sumPUPt );
    	muon_chargedIso.push_back( muon.pfIsolationR04().sumChargedHadronPt );
    	muon_photonIso.push_back( muon.pfIsolationR04().sumPhotonEt );
    	muon_neutralHadIso.push_back( muon.pfIsolationR04().sumNeutralHadronEt );

		// IP significance
		muon_ip3dSignificance.push_back( muon.dB(pat::Muon::PV3D)/muon.edB(pat::Muon::PV3D) );
		muon_dB.push_back( muon.dB(pat::Muon::PV3D) );
		muon_edB.push_back( muon.edB(pat::Muon::PV3D) ); // both muon_edB and muon_dB are set to the maximum possible value of a double
		// also checked PV2D, BS2D, BS3D, and PVDZ for muon.dB and muon.edB. These are all set to the maximum possible value of a double
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
	
	uint32_t run  = iEvent.id().run();
	uint32_t lumi = iEvent.luminosityBlock();
	uint32_t evt  = iEvent.id().event();
	uint32_t seed = run + lumi + evt;  // simple combo
	rand_.SetSeed(seed);

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

	// for( auto &jet : *jets ) { // uncorrected jets
	// for( auto &jet : *jetsPuppiCorr ) { // GK, corrected PF ak4 jets PUPPI
	for( auto &jet : *jetsCorr ) { // GK, corrected PF ak4 jets CHS

		if( jet.pt() < 20 || fabs(jet.eta()) > 2.0 ) continue;

		if( debug ) cout<<" ------ jet idx"<<n_jet<<endl;  

		n_jet++;

		// ----- Basics ----- // 
		// ----- uncorrected jet quantities (no JECs!) ----- // GK
		auto rawP4 = jet.correctedP4("Uncorrected");
    	jetRaw_Pt.push_back( rawP4.pt() );
    	jetRaw_E.push_back( rawP4.energy() );
		
		// ----- JEC only ----- // 
		if (isData) { // for data, no JER needed, so save as _Pt and _E
			jet_E.push_back( jet.energy() );
			jet_Pt.push_back( jet.pt() );
			jet_Mass.push_back( jet.mass() );
		}
		if (!isData) { // for MC, specify if JER is applied or not. Save JEC + JER as _Pt and _E
			jet_E_noJER.push_back( jet.energy() );
			jet_Pt_noJER.push_back( jet.pt() );
			jet_Mass_noJER.push_back( jet.mass() );
		}
		jet_Eta.push_back( jet.eta() );
		jet_Phi.push_back( jet.phi() );

		// // JEC, jet energy scale uncertainty
		// jecUnc_->setJetPt(jet.pt());
		// jecUnc_->setJetEta(jet.eta());
		// double unc = jecUnc_->getUncertainty(true);  // true = up variation
		// double factor_jesUp = (1 + unc);
		// double factor_jesDown = (1 - unc);
		// // correct 4 vector instead of each component
		// auto jet_JES_up = jet.p4() * factor_jesUp;
		// auto jet_JES_down = jet.p4() * factor_jesDown;
		// // Save to tree
		// jet_Pt_JES_up.		push_back(jet_JES_up.pt());
		// jet_E_JES_up.		push_back(jet_JES_up.energy());
		// jet_Mass_JES_up.	push_back(jet_JES_up.mass());
		// jet_Pt_JES_down.	push_back(jet_JES_down.pt());
		// jet_E_JES_down.		push_back(jet_JES_down.energy());
		// jet_Mass_JES_down.	push_back(jet_JES_down.mass());

		// ----- JEC + JER (only for MC) ----- // GK
		if (!isData_) {
			// build parameters for JER
			JME::JetParameters params;
			params.setJetPt(jet.pt());
			params.setJetEta(jet.eta());
			params.setRho(*rhoFastjetAll);
			params.setJetArea(jet.jetArea());

			double res = jerRes_.getResolution(params); // no up/down variations, uncertainty is on the scale factors
			// variation in JER from SF (nominal, up, down)
			double sf_nom = jerSF_.getScaleFactor(params, Variation::NOMINAL);
			double sf_up  = jerSF_.getScaleFactor(params, Variation::UP);
			double sf_down= jerSF_.getScaleFactor(params, Variation::DOWN);

			auto smearJetFactor = [&](const pat::Jet &jet, double sf) {
				double smearFactor = 1.0;
				// Option A: gen-jet match (if available)
				if (jet.genJet()) {
					double dPt = jet.pt() - jet.genJet()->pt();
					smearFactor = 1.0 + (sf - 1.0) * dPt / jet.pt();
				}
				// Option B: random smearing (no gen match)
				else {
					double sigma = res * std::sqrt(std::max(sf*sf - 1, 0.));
					smearFactor = 1.0 + rand_.Gaus(0, sigma);
				}
				return std::max(0.0, smearFactor);
			};

			auto jet_JER = jet.p4() * smearJetFactor(jet, sf_nom);
			jet_Pt			.push_back(jet_JER.pt());
			jet_E			.push_back(jet_JER.energy());
			jet_Mass		.push_back(jet_JER.mass());

			auto jet_JER_up = jet.p4() * smearJetFactor(jet, sf_up);
			jet_Pt_JER_up	.push_back(jet_JER_up.pt());
			jet_E_JER_up	.push_back(jet_JER_up.energy());
			jet_Mass_JER_up	.push_back(jet_JER_up.mass());
			
			auto jet_JER_down = jet.p4() * smearJetFactor(jet, sf_down);
			jet_Pt_JER_down	.push_back(jet_JER_down.pt());
			jet_E_JER_down	.push_back(jet_JER_down.energy());
			jet_Mass_JER_down.push_back(jet_JER_down.mass());
		}

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

		jet_PileupE.push_back( jet.pileup() );
		jet_PileupId.push_back( jet.userFloat("pileupJetId:fullDiscriminant") );
		jet_PileupIdFlag.push_back( jet.userInt("pileupJetId:fullId") ); //A bit map for loose, medium, and tight working points

		// ----- BTagging  ----- //        
		auto jet_pair_discrim = jet.getPairDiscri();
		jet_DeepCSV_prob_b.push_back( jet.bDiscriminator("pfDeepCSVJetTags:probb") );
		jet_DeepCSV_prob_c.push_back( jet.bDiscriminator("pfDeepCSVJetTags:probc") );
		jet_DeepCSV_prob_bb.push_back( jet.bDiscriminator("pfDeepCSVJetTags:probbb") );
		jet_DeepCSV_prob_udsg.push_back( jet.bDiscriminator("pfDeepCSVJetTags:probudsg") );

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
		/* GK remove ebRecHits for LLP HCAL skim

		vector<uint> jet_EcalRechitIDs_temp;
		vector<uint> jet_EcalRechitIndices_temp;

		for( uint ih = 0; ih < ebRecHits->size(); ih++ ){

			const EcalRecHit *recHit = &(*ebRecHits)[ih];
			const DetId recHitId = recHit->detid();

			// Barrel
			if( recHit->detid().subdetId() != EcalBarrel ) continue;
			const auto recHitPos = caloGeometry_EB->getGeometry(recHitId)->getPosition();

			// DR
			if( deltaR( jet.eta(), jet.phi(), recHitPos.eta(), recHitPos.phi()) > 0.5 ) continue; 
			// Energy 
			if( recHit->energy() < 0.1 ) continue;

			jet_EcalRechitIDs_temp.push_back( 0 ); // FIX DetId );
			jet_EcalRechitIndices_temp.push_back( ih );
			//cout << "Save this Rechit: " << q << " | " << SaveThisEBRechit[q] << " : " << recHit->energy() << " " << recHitPos.eta() << " " << recHitPos.phi() << "\n";	
		}

		jet_EcalRechitIDs.push_back( jet_EcalRechitIDs_temp );
		jet_EcalRechitIndices.push_back( jet_EcalRechitIndices_temp );
		*/

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

	// ********************************************************
	// AK8 Jets
	// ********************************************************

	/*
	cout<<" ---- ak8 jets"<<endl;  

	for ( auto &j : *jetsAK8) {

		if (j.pt() < 50) continue;
		nLRJets++;

		//-------------------
		//Fill Jet-Level Info
		//-------------------

		LRJet_E.push_back( j.energy() );
		LRJet_Pt.push_back( j.pt() );
		LRJet_Eta.push_back( j.eta() );
		LRJet_Phi.push_back( j.phi() );
		LRJet_UncorrectedPt.push_back( j.correctedP4(0).Pt() );

		// LRJet_ICSV.push_back( j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
		// LRJet_BoostedDoubleSVAK8.push_back( j.bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags") );

		LRJet_SoftDropM.push_back( (float) j.userFloat("ak8PFJetsCHSSoftDropMass") );
		LRJet_Tau1.push_back(  (float) j.userFloat("NjettinessAK8CHS:tau1") );
		LRJet_Tau2.push_back(  (float) j.userFloat("NjettinessAK8CHS:tau2") );
		LRJet_Tau3.push_back(  (float) j.userFloat("NjettinessAK8CHS:tau3") );

		//Q: Is this correct?
		//LRJet_PassIDLoose.push_back( passJetID(&j, 0) );
		//LRJet_PassIDTight.push_back( passJetID(&j, 1) );

		LRJet_ElectronEnergyFraction.push_back( j.electronEnergyFraction() );
		LRJet_PhotonEnergyFraction.push_back( j.neutralEmEnergyFraction() );
		LRJet_ChargedHadronEnergyFraction.push_back( j.chargedHadronEnergyFraction() );
		LRJet_NeutralHadronEnergyFraction.push_back( j.neutralHadronEnergyFraction() );
		LRJet_MuonEnergyFraction.push_back( j.muonEnergyFraction() );
		LRJet_ChargedHadronMultiplicity.push_back( j.chargedHadronMultiplicity() );
		LRJet_NeutralHadronMultiplicity.push_back( j.neutralHadronMultiplicity() );
		LRJet_PhotonMultiplicity.push_back( j.photonMultiplicity() );
		LRJet_ElectronMultiplicity.push_back( j.electronMultiplicity() );
		LRJet_MuonMultiplicity.push_back( j.muonMultiplicity() );

		//---------------------------
		//Trackless variables
		//---------------------------
		TLorentzVector thisJet;
		thisJet.SetPtEtaPhiE(LRJet_Pt.back(), LRJet_Eta.back(), LRJet_Phi.back(), LRJet_E.back());

		float alphaMax(0.0),medianTheta2D(0.0),medianIP(0.0),minDeltaRAllTracks(0.0),minDeltaRPVTracks(0.0),ptAllTracks(0.0), ptAllPVTracks(0.0);
		int nTracksPV(0);
		findTrackingVariables(thisJet,iSetup,alphaMax,medianTheta2D,medianIP,nTracksPV,ptAllPVTracks,ptAllTracks, minDeltaRAllTracks, minDeltaRPVTracks);

		LRJet_AlphaMax.push_back( alphaMax );
		LRJet_BetaMax.push_back( alphaMax * ptAllTracks/(j.pt()) );
		LRJet_GammaMax.push_back( alphaMax * ptAllTracks/(j.energy()) );
		LRJet_GammaMax_EM.push_back( alphaMax * ptAllTracks/(j.energy()*(j.chargedEmEnergyFraction()+j.neutralEmEnergyFraction())) );
		LRJet_GammaMax_Hadronic.push_back( alphaMax * ptAllTracks/(j.energy()*(j.chargedHadronEnergyFraction()+j.neutralHadronEnergyFraction())) );
		LRJet_GammaMax_ET.push_back( alphaMax * ptAllTracks/j.et() );
		LRJet_MedianTheta2D.push_back( medianTheta2D );
		LRJet_MedianIP.push_back( medianIP );
		LRJet_PtAllPVTracks.push_back( ptAllPVTracks );
		LRJet_PtAllTracks.push_back( ptAllTracks );
		LRJet_MinDeltaRAllTracks.push_back( minDeltaRAllTracks );
		LRJet_MinDeltaRPVTracks.push_back( minDeltaRPVTracks );

		// cout << "FatJet: "
		//   << j.userFloat("pileupJetId:fullDiscriminant") << " "
		//   << j.userInt("pileupJetId:fullId") << " "
		//   << "\n";

		//---------------------------
		//Find PFCandidates Inside the Jet
		//---------------------------

		 FIX
		for (uint q=0; q< pfCands->size(); q++) {
			//const reco::PFCandidate *p = &(*pfCands)[q];

			reco::PFCandidatePtr p_ptr(pfCands,q);
			bool found = false;
			for (uint l=0; l < j.getPFConstituents().size(); l++) {
		if (p_ptr == j.getPFConstituents()[l]) {
			found = true;
			break;
		}
			}

			if (found) {
		SaveThisPFCandidate[q] = true;
		LRJetAllPFCandIndex[nFatJets][LRJetNPFCands[nFatJets]] = q;
		LRJetNPFCands[nFatJets]++;
			}
		}

		int LRJet_NMatchedRechits = 0;
		for (uint q=0; q<ebRecHits->size(); q++) {
			const EcalRecHit *recHit = &(*ebRecHits)[q];
			const DetId recHitId = recHit->detid();
			const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();

			//save the rechits that are within DR 0.5 of the jet axis
			if ( deltaR(LRJetEta[nFatJets], LRJetPhi[nFatJets], recHitPos.eta(), recHitPos.phi())  < 1.0
	 && recHit->energy() > 0.2
	 ) {
			SaveThisEBRechit[q] = true;
			}

			if (recHit->checkFlag(EcalRecHit::kSaturated) || recHit->checkFlag(EcalRecHit::kLeadingEdgeRecovered) || recHit->checkFlag(EcalRecHit::kPoorReco) || recHit->checkFlag(EcalRecHit::kWeird) || recHit->checkFlag(EcalRecHit::kDiWeird)) continue;
			if (recHit->timeError() < 0 || recHit->timeError() > 100) continue;
			if (abs(recHit->time()) > 12.5) continue;

			//Calculate jet timestamps
			if ( deltaR(LRJetEta[nFatJets], LRJetPhi[nFatJets], recHitPos.eta(), recHitPos.phi())  < 0.8) {
		if (recHit->energy() > Rechit_cut) {
			LRJetRechitE_Error[nFatJets] += recHit->energyError() * recHit->energyError();
			LRJetRechitE[nFatJets] += recHit->energy();
			LRJetRechitT[nFatJets] += recHit->time()*recHit->energy();
			LRJetRechitT_rms[nFatJets] += recHit->time()*recHit->time();
			LRJet_NMatchedRechits++;
		}
			}
		}

		LRJetRechitT[nFatJets] = LRJetRechitT[nFatJets]/LRJetRechitE[nFatJets];
		LRJetNRechits[nFatJets] = LRJet_NMatchedRechits;
		LRJetRechitE_Error[nFatJets] = sqrt(LRJetRechitE_Error[nFatJets]);
		LRJetRechitT_rms[nFatJets] = sqrt(LRJetRechitT_rms[nFatJets]);

		for (uint q=0; q<eeRecHits->size(); q++) {
			const EcalRecHit *recHit = &(*eeRecHits)[q];
			const DetId recHitId = recHit->detid();
			const auto recHitPos = endcapGeometry->getGeometry(recHitId)->getPosition();

			//save the rechits that are within DR 0.5 of the jet axis
			if ( deltaR(LRJetEta[nFatJets], LRJetPhi[nFatJets], recHitPos.eta(), recHitPos.phi())  < 1.0
	 && recHit->energy() > 0.2
	 ) {
			SaveThisEERechit[q] = true;
			}
		}

	 //loop over hcal hits
	for (unsigned int iHit = 0; iHit < hcalRecHitsHBHE->size(); iHit ++){
		const HBHERecHit *recHit = &(*hcalRecHitsHBHE)[iHit];
		
		double hiteta = -999;
		double hitphi = -999;
		if (recHit->energy() < 0.1) continue;
		const HcalDetId recHitId = recHit->detid();
		if (recHit->detid().subdetId() == HcalBarrel) {
			const auto recHitPos = hbGeometry->getGeometry(recHitId)->getPosition();
			hiteta = recHitPos.eta();
			hitphi = recHitPos.phi();
		} else if (recHit->detid().subdetId() == HcalEndcap) {
			const auto recHitPos = heGeometry->getGeometry(recHitId)->getPosition();
			hiteta = recHitPos.eta();
			hitphi = recHitPos.phi();
		} else {
			cout << "Error: HCAL Rechit has detId subdet = " << recHit->detid().subdetId() << "  which is not HcalBarrel or HcalEndcap. skipping it. \n";
		}
		
		if ( deltaR(LRJetEta[nFatJets], LRJetPhi[nFatJets], hiteta, hitphi)  < 1.0 ) {
			SaveThisHCALRechit[iHit] = true;
		}
		} 
	
	
		//loop over HO hits
		for (unsigned int iHit = 0; iHit < hcalRecHitsHO->size(); iHit ++){
			const HORecHit *recHit = &(*hcalRecHitsHO)[iHit];
			const DetId recHitId = recHit->detid();
			double hiteta = -999;
			double hitphi = -999;
			if (recHit->energy() < 0.1) continue;
			const auto recHitPos = hoGeometry->getGeometry(recHitId)->getPosition();
			hiteta = recHitPos.eta();
			hitphi = recHitPos.phi();

			//save all HO rechits which have more than 1.5 GeV of energy
			if (recHit->energy() > 1.5) SaveThisHORechit[iHit] = true;
			if ( deltaR(LRJetEta[nFatJets], LRJetPhi[nFatJets], hiteta, hitphi)  < 1.0
 ) {
		SaveThisHORechit[iHit] = true;
			}
		}



		//loop over tracks
		for (unsigned int iTrack = 0; iTrack < generalTracks->size(); iTrack ++){
			reco::Track generalTrack = generalTracks->at(iTrack);
			if ( deltaR(LRJetEta[nFatJets], LRJetPhi[nFatJets], generalTrack.eta(), generalTrack.phi())  < 1.0
 && generalTrack.pt() > 1
 ) {
		SaveThisTrack[iTrack] = true;
			}
		}
		

		//if (nFatJets > OBJECTARRAYSIZE) {
		//  cout << "ERROR: nFatJets exceeded maximum array size: " << OBJECTARRAYSIZE << "\n";
		//  assert(false);
		//}

	}*/

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

	/* GK remove ebRecHits for LLP HCAL skim
	for( uint i = 0; i < ebRecHits->size(); i++ ){

		if( debug ) cout<<" -- ecalrechit idx "<<i<<endl;

		const EcalRecHit *recHit = &(*ebRecHits)[i];
		const DetId recHitId     = recHit->detid();
		if( recHit->detid().subdetId() != EcalBarrel ) continue;

		// Check if we should save
		bool save_hit = false; 
		save_hit = true;

		// Save all hits for now (comment below section)

		for( int ij = 0; ij < n_jet; ij++ ){
			
			for( int ijh = 0; ijh < (int)jet_EcalRechitIndices.at(ij).size(); ijh++ ){
				if( ih == jet_EcalRechitIndices.at(ij).at(ijh) ){
					save_hit = true; 
					break;
				}
			}

			if( save_hit == true ) continue;
		}

		if( !save_hit ) continue;

		const auto recHitPos = caloGeometry_EB->getGeometry(recHitId)->getPosition();
		n_ecalRechit++;

		ecalRechit_Eta.push_back( recHitPos.eta() );
		ecalRechit_Phi.push_back( recHitPos.phi() );
		ecalRechit_E.push_back( recHit->energy() );
		// Features
		ecalRechit_X.push_back( recHitPos.x() );
		ecalRechit_Y.push_back( recHitPos.y() );
		ecalRechit_Z.push_back( recHitPos.z() );		
		ecalRechit_time.push_back( recHit->time() );
		ecalRechit_E_err.push_back( recHit->energyError() );
		ecalRechit_time_err.push_back( recHit->timeError() );
		ecalRechit_kSaturatedflag.push_back( recHit->checkFlag(EcalRecHit::kSaturated) );
		ecalRechit_kLeadingEdgeRecoveredflag.push_back( recHit->checkFlag(EcalRecHit::kLeadingEdgeRecovered) );
		ecalRechit_kPoorRecoflag.push_back( recHit->checkFlag(EcalRecHit::kPoorReco) );
		ecalRechit_kWeirdflag.push_back( recHit->checkFlag(EcalRecHit::kWeird) );
		ecalRechit_kDiWeirdflag.push_back( recHit->checkFlag(EcalRecHit::kDiWeird) );	
	} 
	*/

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

		// Check if we should save //

		if( recHit->energy() < 0.1 ) continue;

		// Only save hits around a jet
                bool save_hit = false;
		for( int ij = 0; ij < n_jet; ij++ ){
			
			for( int ijh = 0; ijh < (int)jet_HcalRechitIndices.at(ij).size(); ijh++ ){
				if( ih == jet_HcalRechitIndices.at(ij).at(ijh) ){
					save_hit = true; 
					break;
				}
			}

			if( save_hit == true ) break;
		}

		if( !save_hit ) continue;

		const auto recHitPos = caloGeometry_HB->getGeometry(recHitId)->getPosition();

		n_hbheRechit++;

		hbheRechit_Eta.push_back( recHitPos.eta() );
		hbheRechit_Phi.push_back( recHitPos.phi() );
		hbheRechit_E.push_back( recHit->energy() );
		// Features  
		hbheRechit_time.push_back( recHit->time() );
		int six_bits_mask = 0x3f;
		int ts = 3; // aux TDC contains TS 0 - 5, with ts3 as SOI
		int SOI_TDC = CaloRecHitAuxSetter::getField(recHit->auxTDC(), six_bits_mask, ts * 6);
		hbheRechit_auxTDC.push_back( SOI_TDC );
		//hbheRechit_auxTDC.push_back( recHit->auxTDC() );
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


//Needs to be called AFTER jet filler and after muon detector cluster filler
/*
bool DisplacedHcalJetNTuplizer::fillHitsTracksAndPFCands(const edm::Event& iEvent, const edm::EventSetup& iSetup){ // GK 
// added "const edm::Event& iEvent," to fix error about "Called EventSetupRecord::get without using a ESGetToken." This portion not working yet!

	// not working with CaloGeometry uncommented! 
	std::cout << "just before getting calo geometry" << std::endl;
	//  const CaloGeometry* geoHandle = &iSetup.getData(caloGeometryToken_); // compiles with this line, referenced from https://cmssdt.cern.ch/lxr/source/Calibration/HcalAlCaRecoProducers/plugins/AlCaHcalIsotrkProducer.cc
	auto const& geoHandle = iSetup.getData(caloGeometryToken_); // compiles with this line, similar to how muon system is done, and recommended from twiki
	// both lines above give an error in running though, about : CastorGeometryRecord  Exception Message: No "PCastorRcd" record found in the EventSetup. 
	// added includes of PCaloTowerRcd.h, ESGetToken of PCastorRcd, and ESConsumes. Still have error
	std::cout << "just after getting calo geometry" << std::endl;

	//  edm::ESHandle<CaloGeometry> geoHandle;
	geoHandle = iSetup.get<CaloGeometryRecord>(geoHandleToken);
	//  iSetup.get<CaloGeometryRecord>().get(geoHandle);
	const CaloSubdetectorGeometry *barrelGeometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
	const CaloSubdetectorGeometry *endcapGeometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);
	const CaloSubdetectorGeometry *hbGeometry = geoHandle->getSubdetectorGeometry(DetId::Hcal, HcalBarrel);
	const CaloSubdetectorGeometry *heGeometry = geoHandle->getSubdetectorGeometry(DetId::Hcal, HcalEndcap);
	const CaloSubdetectorGeometry *hoGeometry = geoHandle->getSubdetectorGeometry(DetId::Hcal, HcalOuter); */

	// ********************************************************
	// Save EB Rechits inside Jets and AK8 Jets
	// ********************************************************
	/*  for (uint q=0; q<ebRecHits->size(); q++) {
		const EcalRecHit *recHit = &(*ebRecHits)[q];
		const DetId recHitId = recHit->detid();
		const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
		if (SaveThisEBRechit[q]) {
			ecalRechit_Eta[nRechits] = recHitPos.eta();
			ecalRechit_Phi[nRechits] = recHitPos.phi();
			ecalRechit_E[nRechits] = recHit->energy();
			ecalRechit_T[nRechits] = recHit->time();
			ecalRechit_E_Error[nRechits] = recHit->energyError();
			ecalRechit_T_Error[nRechits] = recHit->timeError();
			ecalRechit_kSaturatedflag[nRechits] = recHit->checkFlag(EcalRecHit::kSaturated);
			ecalRechit_kLeadingEdgeRecoveredflag[nRechits] = recHit->checkFlag(EcalRecHit::kLeadingEdgeRecovered);
			ecalRechit_kPoorRecoflag[nRechits] = recHit->checkFlag(EcalRecHit::kPoorReco);
			ecalRechit_kWeirdflag[nRechits]= recHit->checkFlag(EcalRecHit::kWeird);
			ecalRechit_kDiWeirdflag[nRechits] = recHit->checkFlag(EcalRecHit::kDiWeird);
			nRechits++;

			if (nRechits > RECHITARRAYSIZE) {
	cout << "ERROR: nRechits exceeded maximum array size: " << RECHITARRAYSIZE << "\n";
	assert(false);
			}
		}
		//cout << "After Rechit: " << q << " | " << SaveThisEBRechit[q] << " : " << recHit->energy() << " " << recHitPos.eta() << " " << recHitPos.phi() << "\n";
		} */

	// ********************************************************
	// Save EE Rechits inside Jets and AK8 Jets
	// ********************************************************
	/*  for (uint q=0; q<eeRecHits->size(); q++) {
		const EcalRecHit *recHit = &(*eeRecHits)[q];
		const DetId recHitId = recHit->detid();
		const auto recHitPos = endcapGeometry->getGeometry(recHitId)->getPosition();
		if (SaveThisEERechit[q]) {
			ecalRechit_Eta[nRechits] = recHitPos.eta();
			ecalRechit_Phi[nRechits] = recHitPos.phi();
			ecalRechit_E[nRechits] = recHit->energy();
			ecalRechit_T[nRechits] = recHit->time();
			ecalRechit_E_Error[nRechits] = recHit->energyError();
			ecalRechit_T_Error[nRechits] = recHit->timeError();
			ecalRechit_kSaturatedflag[nRechits] = recHit->checkFlag(EcalRecHit::kSaturated);
			ecalRechit_kLeadingEdgeRecoveredflag[nRechits] = recHit->checkFlag(EcalRecHit::kLeadingEdgeRecovered);
			ecalRechit_kPoorRecoflag[nRechits] = recHit->checkFlag(EcalRecHit::kPoorReco);
			ecalRechit_kWeirdflag[nRechits]= recHit->checkFlag(EcalRecHit::kWeird);
			ecalRechit_kDiWeirdflag[nRechits] = recHit->checkFlag(EcalRecHit::kDiWeird);
			nRechits++;

			if (nRechits > RECHITARRAYSIZE) {
	cout << "ERROR: nRechits exceeded maximum array size: " << RECHITARRAYSIZE << "\n";
	assert(false);
			}
		}
		//cout << "After Rechit: " << q << " | " << SaveThisEBRechit[q] << " : " << recHit->energy() << " " << recHitPos.eta() << " " << recHitPos.phi() << "\n";
		} */


	// SX: Disabled this because it was too many hits and causing the output to be excessively large
	// // ********************************************************
	// // Save HCAL Rechits inside Jets and AK8 Jets
	// // ********************************************************
	/*)
	for (unsigned int iHit = 0; iHit < hcalRecHitsHBHE->size(); iHit ++){
		const HBHERecHit *recHit = &(*hcalRecHitsHBHE)[iHit];
		const HcalDetId recHitId = recHit->detid();

		if (iHit == 0) cout << "HCALREchit iHit : subdet (HB=1, HE=2) : depth ieta iphi | energy" << "\n";
		cout << "HCALREchit " << iHit << " : " << recHit->detid().subdetId() << " : " << recHitId.depth() << " " << recHitId.ieta() << " " << recHitId.iphi() << " "
	 << " | " << recHit->energy() << " "
	 << "\n";

		if (SaveThisHCALRechit[iHit]) { // if we want to do jet matching, would likely only save HCAL rechits inside a jet. This relies on geometry, as in L4356 where the HCAL rechits inside a jet are marked to save
			hbheRechit_iEta[nHBHERechits]  = recHitId.ieta();
			hbheRechit_iPhi[nHBHERechits]  = recHitId.iphi();
			hbheRechit_depth[nHBHERechits]  = recHitId.depth();
			if (recHit->detid().subdetId() == HcalBarrel) {
	  const auto recHitPos = hbGeometry->getGeometry(recHitId)->getPosition();
		hbheRechit_Phi[nHBHERechits] = recHitPos.phi();
		hbheRechit_Eta[nHBHERechits] = recHitPos.eta();
		hbheRechit_X[nHBHERechits] = recHitPos.x();
		hbheRechit_Y[nHBHERechits] = recHitPos.y();
		hbheRechit_Z[nHBHERechits] = recHitPos.z(); 
 } else if (recHit->detid().subdetId() == HcalEndcap) {
		const auto recHitPos = heGeometry->getGeometry(recHitId)->getPosition();
		hbheRechit_Phi[nHBHERechits] = recHitPos.phi();
		hbheRechit_Eta[nHBHERechits] = recHitPos.eta();
		hbheRechit_X[nHBHERechits] = recHitPos.x();
		hbheRechit_Y[nHBHERechits] = recHitPos.y();
		hbheRechit_Z[nHBHERechits] = recHitPos.z();  
	} else {
		cout << "Error: HCAL Rechit has detId subdet = " << recHit->detid().subdetId() << "  which is not HcalBarrel or HcalEndcap. skipping it. \n";
			}

			hbheRechit_E[nHBHERechits] = recHit->energy();
			// GK various HCAL rechit time variables 
			hbheRechit_time[nHBHERechits] = recHit->time(); // time() is MAHI fit time in ns
			hbheRechit_auxTDC[nHBHERechits] = recHit->auxTDC(); // contains TDC values in 6 bits in HE (2 bits in HB) for each of 5 TS (timeslice 3 = SOI). Copy of HCAL Digis. Extra bit set to distinguish case of data packed from no data
			hbheRechit_timeFalling[nHBHERechits] = recHit->timeFalling(); // HBHE rise time
			//      hbheRechit_cctime[nHBHERechits] = recHit->ccTime(); // cross correlation time, https://indico.cern.ch/event/1142347/contributions/4793749/attachments/2412936/4129323/HCAL%20timing%20update.pdf

			if (hbheRechit_E[nHORechits] < -1) {
	cout << "HCAL Hit: " << hbheRechit_Eta[nHORechits] << " " << hbheRechit_Phi[nHBHERechits] << " : " << hbheRechit_E[nHORechits] << " | " << nHORechits << "\n";
			}

			nHBHERechits++;

			if (nHBHERechits > RECHITARRAYSIZE) {
		cout << "ERROR: nHBHERechits exceeded maximum array size: " << RECHITARRAYSIZE << "\n";
		//assert(false);
		break;
			}
		}
	}

	// ********************************************************
	// Save HO Rechits inside Jets and AK8 Jets
	// ********************************************************
	for (unsigned int iHit = 0; iHit < hcalRecHitsHO->size(); iHit ++){
		const HORecHit *recHit = &(*hcalRecHitsHO)[iHit];
		if (SaveThisHORechit[iHit]) {
				 const DetId recHitId = recHit->detid();
			const auto recHitPos = hoGeometry->getGeometry(recHitId)->getPosition();
			hoRechit_Phi[nHORechits] = recHitPos.phi();
			hoRechit_Eta[nHORechits] = recHitPos.eta();
			hoRechit_X[nHORechits] = recHitPos.x();
			hoRechit_Y[nHORechits] = recHitPos.y();
			hoRechit_Z[nHORechits] = recHitPos.z();
			hoRechit_E[nHORechits] = recHit->energy();
			hoRechit_T[nHORechits] = recHit->time(); 
			nHORechits ++;

			if (nHORechits > HORECHITARRAYSIZE) {
	cout << "ERROR: nHORechits exceeded maximum array size: " << HORECHITARRAYSIZE << "\n";
	assert(false);
			}
		}
	} 

	// ********************************************************
	// Save PF Candidates inside Jets and AK8 Jets
	// ********************************************************
	for (uint q=0; q< pfCands->size(); q++) {
		if (SaveThisPFCandidate[q]) {
			const reco::PFCandidate *p = &(*pfCands)[q];

			PFCandidatePdgId[nPFCandidates] = p->pdgId();
			PFCandidatePt[nPFCandidates] = p->pt();
			PFCandidateEta[nPFCandidates] = p->eta();
			PFCandidatePhi[nPFCandidates] = p->phi();
			PFCandidateTrackIndex[nPFCandidates] = -1;
			PFCandidatePVIndex[nPFCandidates] = -1;

			//find the matching PV
			const reco::VertexRef &PVOrig = (*primaryVertexAssociation)[reco::CandidatePtr(pfCands,q)];
			if(PVOrig.isNonnull()) {

	int matchedPVIndex = -1;
	for (int ipv = 0; ipv < nPVAll; ++ipv) {
		const reco::VertexRef vtxRef(vertices,ipv);
		if (PVOrig.id() == vtxRef.id() && PVOrig.key() == vtxRef.key()) {
			matchedPVIndex = ipv;
			break;
		}
	}

	// if (matchedPVIndex>=0) {
	//   cout << "PFCandidate " << q << " : " << matchedPVIndex << " --> " << pvAllX[matchedPVIndex] << " " << pvAllY[matchedPVIndex] << " " << pvAllZ[matchedPVIndex] << " | "
	//  << PVOrig->x() << " " << PVOrig->y() << " " << PVOrig->z() << " "
	//  << " \n";
	// }

	PFCandidatePVIndex[nPFCandidates] = matchedPVIndex;
			}

			//find track ref
			if (p->trackRef().isNonnull()) {
	for (unsigned int iTrack = 0; iTrack < generalTrackHandle->size(); iTrack ++){
		reco::TrackBaseRef tref(generalTrackHandle,iTrack);
		reco::TrackBaseRef pTrackRef(p->trackRef());
		//cout << "Track " << p->trackRef().id() << " " << p->trackRef().key() << " | " << pTrackRef.id() << " " << pTrackRef.key() << " | " << iTrack << " " << tref.id() << " " << tref.key() << " \n";
		if (pTrackRef.id() == tref.id() && pTrackRef.key() == tref.key() ) {
			//cout << "found track match: " << p->trackRef()->pt() << " " << p->trackRef()->eta() << " " << p->trackRef()->phi() << " | "
			//   << tref->pt() << " " << tref->eta() << " " << tref->phi() << " "
			//   << "\n";
			SaveThisTrack[iTrack] = true;
			PFCandidateGeneralTrackIndex[nPFCandidates] = iTrack;
			break;
		}
	}
			}
			PFCandToSavedPFCandMap[q] = nPFCandidates;
			nPFCandidates++;
			if (nPFCandidates > MAX_NPFCAND) {
	cout << "ERROR: nPFCandidates exceeded maximum array size: " << MAX_NPFCAND << "\n";
	assert(false);
			}
		}
	}

	// GK
	// ********************************************************
	// Save Tracks inside Jets and AK8 Jets
	// ********************************************************
	// Magnetic field
	//  auto const& magneticField = iSetup.getData(magneticFieldToken_); // GK
	edm::ESHandle<MagneticField> magneticField;
	iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
	magneticFieldGeometry = &*magneticField; 
	std::string thePropagatorName_ = "PropagatorWithMaterial";
	iSetup.get<TrackingComponentsRecord>().get(thePropagatorName_,thePropagator_);
	StateOnTrackerBound stateOnTracker(thePropagator_.product());

	for (unsigned int iTrack = 0; iTrack < generalTrackHandle->size(); iTrack ++){


		// reco::Track generalTrack = generalTracks->at(iTrack);
		// //const auto& generalTrack = generalTracks->at(iTrack);
		// // TLorentzVector generalTrackVecTemp;
		// // generalTrackVecTemp.SetPtEtaPhiM(generalTrack.pt(),generalTrack.eta(),generalTrack.phi(),0);

		if (SaveThisTrack[iTrack]) {

			reco::TrackBaseRef tref(generalTrackHandle,iTrack);
			// make transient track (unfolding effects of B field ?)

			reco::TransientTrack tt(generalTrackHandle->at(iTrack),magneticFieldGeometry);

			if(!tt.isValid()) {
	std::cout << "Error: Transient Track not valid ("
			<< tref->pt() << " " << tref->eta() << " " << tref->phi()
			<< "). Skipping the track\n";
	continue;
	} 

			track_Pt[nTracks] = tref->pt();
			track_Eta[nTracks] = tref->eta();
			track_Phi[nTracks] = tref->phi();

			//find the best vertex for this track
			float maxWeight = 0;
			int bestVertexIndex = -1;
			for(int k = 0; k < (int)vertices->size();k++){
	if(vertices->at(k).trackWeight(tref) > maxWeight){
		maxWeight = vertices->at(k).trackWeight(tref);
		bestVertexIndex = k;
	}
			}
			track_bestVertexIndex[nTracks] = bestVertexIndex;
			track_nMissingInnerHits[nTracks] = tref->hitPattern().numberOfLostTrackerHits(reco::HitPattern::MISSING_INNER_HITS);
			track_nMissingOuterHits[nTracks] = tref->hitPattern().numberOfLostTrackerHits(reco::HitPattern::MISSING_OUTER_HITS);
			track_nPixelHits[nTracks] = tref->hitPattern().numberOfValidPixelHits();
			track_nHits[nTracks] = tref->hitPattern().numberOfValidHits();
			track_dxyToBS[nTracks] = tref->dxy(*beamSpot);
			track_dxyErr[nTracks] = tref->dxyError();
			track_dzToPV[nTracks] = tref->dz(beamSpot->position());
			track_dzErr[nTracks] = tref->dzError();
			track_charge[nTracks] = tref->charge();
			track_chi2[nTracks] = tref->chi2();
			track_ndof[nTracks] = tref->ndof();

			// ********************************************************
			// For track angle
			// ********************************************************
			// get track trajectory info
			static GetTrackTrajInfo getTrackTrajInfo;
			vector<GetTrackTrajInfo::Result> trajInfo = getTrackTrajInfo.analyze(iSetup, (*tref));
			if ( trajInfo.size() > 0 && trajInfo[0].valid) {
	// get inner tracker hit from trajectory state
	const TrajectoryStateOnSurface& tsosInnerHit = trajInfo[0].detTSOS;

	//  here's the track angle
	// find beamspot x,y coordinates
	const reco::BeamSpot& pat_beamspot = (*beamSpot);
	TVector2 bmspot(pat_beamspot.x0(),pat_beamspot.y0());
	// find track trajectory state on surface inner hit
	GlobalPoint  innerPos = tsosInnerHit.globalPosition();
	GlobalVector innerMom = tsosInnerHit.globalMomentum();

	// calculate the difference between inner hit and beamspot
	TVector2 sv(innerPos.x(),innerPos.y());
	TVector2 diff = (sv-bmspot);
	//cout<<"bs x: "<<bmspot.X()<<" y: "<<bmspot.Y()<<endl;
	//cout<<" sv x: "<<sv.X()<<" y: "<<sv.Y()<<endl;
	//cout<<" diff phi: "<<diff.Phi()<<endl;
	TVector2 momentum(innerMom.x(),innerMom.y());
	//cout<<" p x: "<<momentum.X()<<" y: "<<momentum.Y()<<endl;
	//cout<<" p phi: "<<momentum.Phi()<<endl;
	//cout<<" dPhi: "<<diff.DeltaPhi(momentum)<<endl;
	track_angle[nTracks] = fabs( diff.DeltaPhi(momentum) ) ;
			}

			TrackToSavedTrackMap[iTrack] = nTracks;
			nTracks++; //increment saved tracks

			if (nTracks > RECHITARRAYSIZE) {
	cout << "ERROR: nTracks exceeded maximum array size: " << RECHITARRAYSIZE << "\n";
	assert(false);
			}
		} //end if save this Track
	} //loop over tracks
	

	// **********************************************************************
	// Cross-reference PFCandidate Indices in jets with Saved PF Candidates
	// **********************************************************************
	if (nJets >= 0) {
		for (uint q=0; q < uint(nJets); q++) {
			if (jetNPFCands[q] >= 0) {
	for (uint l=0; l < uint(jetNPFCands[q]); l++) {
		jetPFCandIndex[q][l] = PFCandToSavedPFCandMap[jetAllPFCandIndex[q][l]];
	}
			}
		}
	}
	if (nFatJets >= 0) {
		for (uint q=0; q < uint(nFatJets); q++) {
			if (jetNPFCands[q] >= 0) {
	for (uint l=0; l < uint(fatJetNPFCands[q]); l++) {
		fatJetPFCandIndex[q][l] = PFCandToSavedPFCandMap[fatJetAllPFCandIndex[q][l]];
	}
			}
		}
	}

	// *************************************************************
	// Cross-reference TrackIndex of PF Candidate with Saved Tracks
	// *************************************************************
	if (nPFCandidates >= 0 ) {
		for (uint q=0; q < uint(nPFCandidates); q++) {
			if (PFCandidateGeneralTrackIndex[q] >= 0) {
	PFCandidateTrackIndex[q] = TrackToSavedTrackMap[PFCandidateGeneralTrackIndex[q]];
			}
		}
	}
	

	return true;

} //end function fillHitsTracksAndPFCands
*/

/*
bool DisplacedHcalJetNTuplizer::fillSecondaryVertices(){

	VertexDistance3D vdist;
	for(unsigned int i = 0; i < secondaryVertices->size(); i++) {

		nSecondaryVertices++;

		// if (SaveThisSecondaryVertex[i]) {
		const reco::VertexCompositePtrCandidate *sv = &(secondaryVertices->at(i));
		GlobalVector flightDir(sv->vertex().x() - PV_global->x(), sv->vertex().y() - PV_global->y(),sv->vertex().z() - PV_global->z());
		Measurement1D dl= vdist.distance(*PV_global,VertexState(RecoVertex::convertPos(sv->position()),RecoVertex::convertError(sv->error())));
			
		// cout << "secVertex " << i << " : " << sv->vertex().x() << " " << sv->vertex().y()  << " " << sv->vertex().z()  << " | " 
		//     << sv->pt() << " " << sv->eta() << " " << sv->phi() << " " << sv->p4().M() << " " 
		//     << sv->charge() << " "
		//     << dl.significance() << " " 
		//     << dl.value() << " "
		//     << dl.error() << " " 
		//     << sv->numberOfSourceCandidatePtrs() << " "
		//     << "\n";  
			
		secVtx_Pt.push_back( sv->pt() );
		secVtx_Eta.push_back( sv->eta() );
		secVtx_Phi.push_back( sv->phi() );
		secVtx_charge.push_back( sv->charge() );
		secVtx_nConstituents.push_back( sv->numberOfSourceCandidatePtrs() );
		secVtx_X.push_back( sv->vertex().x() );
		secVtx_Y.push_back( sv->vertex().y() );
		secVtx_Z.push_back( sv->vertex().z() );
		secVtx_Distance.push_back( dl.value() );
		secVtx_DistanceError.push_back( dl.error() );

		// }
	}
	return true;
}*/


// ------------------------------------------------------------------------------------
bool DisplacedHcalJetNTuplizer::FillPileupBranches(){

	if( debug ) cout<<"Running DisplacedHcalJetNTuplizer::FillPileup"<<endl; 

	n_BunchXing = puInfo->size();
	for(const PileupSummaryInfo &pu : *puInfo){
                BunchXing.push_back( pu.getBunchCrossing() );
		nPU.push_back( pu.getPU_NumInteractions() );
		nPUmean.push_back( pu.getTrueNumInteractions() );
	}
	
	if( debug ) cout<<"Done DisplacedHcalJetNTuplizer::FillPileup"<<endl; 

	return true;
};

// ------------------------------------------------------------------------------------
bool DisplacedHcalJetNTuplizer::FillMCBranches() {

/*
	for(const reco::GenJet &j : *genJets) {
		//std::cout << nGenJets << std::endl;
		if (j.pt() > 10) {
			genJetE[nGenJets] = j.energy();
			genJetPt[nGenJets] = j.pt();
			genJetEta[nGenJets] = j.eta();
			genJetPhi[nGenJets] = j.phi();
			genJetMET[nGenJets] = j.invisibleEnergy();
			nGenJets++;
		}
	}

	const reco::GenMET &GenMetCalo = genMetsCalo->front();
	genMetPtCalo  = GenMetCalo.pt();
	genMetPhiCalo = GenMetCalo.phi();

	const reco::GenMET &GenMetTrue = genMetsTrue->front();
	genMetPtTrue  = GenMetTrue.pt();
	genMetPhiTrue = GenMetTrue.phi();

	bool foundGenVertex = false;
	for(size_t i=0; i<genParticles->size();i++) {
		if (!foundGenVertex) {
			for (unsigned int j=0; j<(*genParticles)[i].numberOfDaughters(); ++j) {
	const reco::Candidate *dau = (*genParticles)[i].daughter(j);
	if (dau) {
		genVertexX = dau->vx();
		genVertexY = dau->vy();
		genVertexZ = dau->vz();
		if(readGenVertexTime_) genVertexT = *genParticles_t0;
		foundGenVertex = true;
		break;
	}
			}
		}
	}

	genWeight = genInfo->weight();
	genSignalProcessID = genInfo->signalProcessID();
	genQScale = genInfo->qScale();
	genAlphaQCD = genInfo->alphaQCD();
	genAlphaQED = genInfo->alphaQED();

	//get lhe weights for systematic uncertainties:
	double nomlheweight = genInfo->weights()[0];

	//fill scale variation weights
	if (genInfo->weights().size()>=10) {
		for (unsigned int iwgt=1; iwgt<10; ++iwgt) {
			//normalize to
			double wgtval = genInfo->weights()[iwgt]*genWeight/genInfo->weights()[1];
			scaleWeights->push_back(wgtval);
		}
	}

	//fill pdf variation weights
	if (firstPdfWeight>=0 && lastPdfWeight>=0 && lastPdfWeight<int(genInfo->weights().size()) && (lastPdfWeight-firstPdfWeight+1)==100) {

		//fill pdf variation weights after converting with mc2hessian transformation
		std::array<double, 100> inpdfweights;
		for (int iwgt=firstPdfWeight, ipdf=0; iwgt<=lastPdfWeight; ++iwgt, ++ipdf) {
			inpdfweights[ipdf] = genInfo->weights()[iwgt]/genInfo->weights()[firstPdfWeight-1];
		}

		std::array<double, 60> outpdfweights;
		pdfweightshelper.DoMC2Hessian(inpdfweights.data(),outpdfweights.data());

		for (unsigned int iwgt=0; iwgt<60; ++iwgt) {
			double wgtval = outpdfweights[iwgt]*genWeight;
			pdfWeights->push_back(wgtval);
		}

		//fill alpha_s variation weights
		if (firstAlphasWeight>=0 && lastAlphasWeight>=0 && lastAlphasWeight<int(genInfo->weights().size())) {
			for (int iwgt = firstAlphasWeight; iwgt<=lastAlphasWeight; ++iwgt) {
	double wgtval = genInfo->weights()[iwgt]*genWeight/nomlheweight;
	alphasWeights->push_back(wgtval);
			}
		}

	}

	//fill sum of weights histograms
	sumWeights->Fill(0.,genWeight);

	for (unsigned int iwgt=0; iwgt<scaleWeights->size(); ++iwgt) {
		sumScaleWeights->Fill(double(iwgt),(*scaleWeights)[iwgt]);
	}
	for (unsigned int iwgt=0; iwgt<pdfWeights->size(); ++iwgt) {
		sumPdfWeights->Fill(double(iwgt),(*pdfWeights)[iwgt]);
	}
	for (unsigned int iwgt=0; iwgt<alphasWeights->size(); ++iwgt) {
		sumAlphasWeights->Fill(double(iwgt),(*alphasWeights)[iwgt]);
	}

	*/
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
	
		// If particle is a daughter of a tau, then save it

		/*if((*genParticles)[i].numberOfMothers() > 0) {
			const reco::Candidate* firstMotherWithDifferentID = findFirstMotherWithDifferentID(&(*genParticles)[i]);
			if (firstMotherWithDifferentID && abs(firstMotherWithDifferentID->pdgId()) == 15 ) {
				// cout << "GenParticles " << i << " : " << (*genParticles)[i].pdgId() << " "
				//      << (*genParticles)[i].status() << " "
				//      << (*genParticles)[i].pt() << " " << (*genParticles)[i].eta() << " "
				//      << " ---> " << firstMotherWithDifferentID->pdgId() << " " << firstMotherWithDifferentID->status() << " "
				//      << firstMotherWithDifferentID->pt() << " " << firstMotherWithDifferentID->eta() << " "
				//      << firstMotherWithDifferentID->phi() << " "
				//      << "\n";
				prunedV.push_back(&(*genParticles)[i]);
				already_saved = true;
				continue;
			}
		}*/

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


