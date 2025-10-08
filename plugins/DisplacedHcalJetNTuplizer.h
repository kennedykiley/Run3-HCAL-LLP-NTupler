// -*- C++ -*-
// Class:      JetNtupler
/*
Description: Base class for miniAOD analysis with CRAB
*/
//         Author:  Kiley Kennedy, Gillian Kopp, & the LPC LLP Team
//         Created: Feb 2023

#ifndef DisplacedHcalJetNTuplizer_H
#define DisplacedHcalJetNTuplizer_H

// system include files
#include <memory>
#include <string>
#include <vector>
#include <tuple>
#include <fstream>
#include <TRandom3.h>

using namespace std;

// CMSSW framework includes
#include "FWCore/Common/interface/TriggerNames.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h"
#include "CommonTools/Egamma/interface/ConversionTools.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CondFormats/GeometryObjects/interface/PCaloGeometry.h" // GK
#include "CondFormats/GeometryObjects/interface/RecoIdealGeometry.h" // GK
#include "JetMETCorrections/Modules/interface/JetResolution.h" // GK
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h" // GK
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h" // GK

#include "DataFormats/CSCDigi/interface/CSCComparatorDigi.h"
#include "DataFormats/CSCDigi/interface/CSCComparatorDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCStripDigi.h"
#include "DataFormats/CSCDigi/interface/CSCStripDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCWireDigi.h"
#include "DataFormats/CSCDigi/interface/CSCWireDigiCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HORecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/CaloRecHitAuxSetter.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SimDataFormats/DigiSimLinks/interface/DTDigiSimLink.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/StripDigiSimLink.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

// Geometry
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/Records/interface/PCastorRcd.h" // GK
#include "Geometry/Records/interface/GEMRecoGeometryRcd.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "PhysicsTools/RecoUtils/interface/CheckHitPattern.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoVertex/ConfigurableVertexReco/interface/ConfigurableVertexReconstructor.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertError.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexTools/interface/VertexCompatibleWithBeam.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

//tracking and vertex
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/GeomPropagators/interface/StateOnTrackerBound.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

// Geometry
//CMSSW package includes
//ECAL Rechits

//ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

//LOCAL includes
//#include "cms_lpc_llp/llp_ntupler/interface/GetTrackTrajInfo.h"
//#include "cms_lpc_llp/llp_ntupler/interface/RazorPDFWeightsHelper.h"

// ------ Class declaration ------ //

// ------------------------------------------------------------------------------------
class DisplacedHcalJetNTuplizer : public edm::one::EDAnalyzer<edm::one::WatchLuminosityBlocks,edm::one::WatchRuns> {
//class DisplacedHcalJetNTuplizer : public edm::EDAnalyzer {
public:
	//analyzer constructor and destructor
	explicit DisplacedHcalJetNTuplizer(const edm::ParameterSet&);
	~DisplacedHcalJetNTuplizer();

	// ----- Load Event ----- // 

	void loadEvent(const edm::Event& iEvent); // call at the beginning of each event to get input handles from the python config

	// ----- Enable Branches ----- // 

	virtual void EnableBranches();
	// Event Level Info
	void EnableEventInfoBranches();
	void EnableTriggerBranches();
	void EnablePVBranches();
	void EnableMetBranches();
	// Standard Objects
	void EnableElectronBranches();
	void EnableMuonBranches();
	void EnableTauBranches();
	void EnablePhotonBranches();
	void EnableJetBranches();
	// Low-Level Objects
	void EnableTrackBranches();
	void EnablePFCandidateBranches();
	void EnableSecondaryVerticesBranches();
	// Hits
	void EnableEcalRechitBranches();
	void EnableHcalRechitBranches();
	// MC
	void EnablePileupBranches();
	void EnableMCBranches();
	void EnableGenParticleBranches();

	// ----- Reset Output Variables in Memory ----- // 

	void ResetBranches();
	// Event Level Info
	void ResetEventInfoBranches();
	void ResetTriggerBranches();
	void ResetPVBranches();
	void ResetMetBranches();
	// Standard Objects
	void ResetElectronBranches();
	void ResetMuonBranches();
	void ResetTauBranches();
	void ResetPhotonBranches();
	void ResetJetBranches();
	// Low-Level Objects
	void ResetTrackBranches();
	void ResetPFCandidateBranches();
	void ResetSecondaryVerticesBranches();
	// Hits
	void ResetEcalRechitBranches();
	void ResetHcalRechitBranches();
	// MC
	void ResetPileupBranches();
	void ResetMCBranches();
	void ResetGenParticleBranches();

	// ----- Fill Branches ----- //   

	//virtual bool FillBranches();
	// Event Level Info
	bool FillEventInfoBranches( const edm::Event& iEvent );
	bool FillPVBranches( const edm::Event& iEvent );
	bool FillTriggerBranches( const edm::Event& iEvent );
	bool FillMetBranches( const edm::Event& iEvent );
	// Standard Objects
	bool FillElectronBranches( const edm::Event& iEvent );
	bool FillMuonBranches( const edm::Event& iEvent );
	bool FillTauBranches( const edm::Event& iEvent );
	bool FillPhotonBranches( const edm::Event& iEvent );
	bool FillJetBranches( const edm::Event& iEvent, const edm::EventSetup& iSetup );
	// Low-Level Objects
	bool FillTrackBranches( const edm::Event& iEvent ); 
	bool FillPFCandidateBranches( const edm::Event& iEvent ); 
	bool FillSecondaryVerticesBranches( const edm::Event& iEvent ); 
	// Hits
	bool FillEcalRechitBranches( const edm::Event& iEvent, const edm::EventSetup& iSetup ); 
	bool FillHcalRechitBranches( const edm::Event& iEvent, const edm::EventSetup& iSetup );
	// MC
	bool FillPileupBranches();
	bool FillMCBranches(); // const edm::Event& iEvent, const edm::EventSetup& iSetup );
	bool FillGenParticleBranches(); // const edm::Event& iEvent, const edm::EventSetup& iSetup );

	// ------ Helper Functions ------ //

	double deltaPhi(double phi1, double phi2);
	double deltaR(double eta1, double phi1, double eta2, double phi2);

protected:
        virtual void beginJob() override;
        virtual void beginRun(const edm::Run&, const edm::EventSetup&) override;
        virtual void beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) override;
    	virtual void endRun(edm::Run const& iRun, edm::EventSetup const&) override;
        virtual void endJob() override;

	// ----- Member data ------ //

	// Control Switches

	bool debug;
	bool isData_;
	bool isSignal_;

	// ====================================================================================
	// EDM Tokens (miniAOD Inputs)
	// ====================================================================================
	// Random engine for smearing // GK
	TRandom3 rand_;
	// JER tools // GK
	JME::JetResolution jerRes_;
	JME::JetResolutionScaleFactor jerSF_;
	// JEC uncertainty
	std::unique_ptr<JetCorrectionUncertainty> jecUnc_;

	edm::EDGetTokenT<edm::TriggerResults> triggerBitsToken_;
	edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsToken_;
	edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescalesToken_;

	string triggerPathNamesFile_;
	//string eleHLTFilterNamesFile_;
	//string muonHLTFilterNamesFile_;
	//string phoHLTFilterNamesFile_;
	//string jetHLTFilterNamesFile_;

	// General Info
	edm::EDGetTokenT<edm::TriggerResults> metFilterBitsToken_;
	edm::EDGetTokenT<reco::VertexCollection> verticesToken_;
  	edm::EDGetTokenT<edm::Association<vector<reco::Vertex> > > primaryVertexAssociationToken_;
  	edm::EDGetTokenT<edm::ValueMap<int> > primaryVertexAssociationValueMapToken_;
	edm::EDGetTokenT<double> rhoFastjetAllToken_;

	// Event Level
	edm::EDGetTokenT<pat::METCollection> metToken_;
	edm::InputTag bsTag_;
	edm::EDGetTokenT<reco::BeamSpot> bsToken_;

	// Physics Objects
	edm::EDGetTokenT<reco::GsfElectronCollection> electronsToken_;
	edm::EDGetTokenT<reco::MuonCollection> muonsToken_;
	edm::EDGetTokenT<pat::TauCollection> tausToken_;
	edm::EDGetTokenT<reco::PhotonCollection> photonsToken_;
	edm::EDGetTokenT<pat::JetCollection> jetsToken_;
	edm::EDGetTokenT<pat::JetCollection> jetsCorrToken_; // GK, adding JECs
	edm::EDGetTokenT<pat::JetCollection> jetsPuppiCorrToken_; // GK, adding JECs, PUPPI
	edm::EDGetTokenT<reco::CaloJetCollection> calojetsToken_;
	edm::EDGetTokenT<pat::JetCollection> LRJetsToken_;
	edm::EDGetTokenT<reco::CaloJetCollection> caloLRJetsToken_;
	edm::EDGetTokenT<BXVector<l1t::Jet>> l1jetsToken_;

	// Low-Level Objects
	//edm::EDGetTokenT<edm::View<reco::Track>> tracksToken_;
	edm::EDGetTokenT<vector<reco::Track>> generalTracksToken_;
	edm::EDGetTokenT<vector<reco::VertexCompositePtrCandidate> > secondaryVerticesToken_;
	edm::EDGetTokenT<reco::PFCandidateCollection> PFCandsToken_;

	// RecHits
	edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> ebRecHitsToken_;
	edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>> hcalRecHitsHBHEToken_;

	// OTHER
	edm::EDGetTokenT<edm::ValueMap<bool> > electron_cutbasedID_decisions_loose_Token_;
	edm::EDGetTokenT<edm::ValueMap<bool> > electron_cutbasedID_decisions_medium_Token_;
	edm::EDGetTokenT<edm::ValueMap<bool> > electron_cutbasedID_decisions_tight_Token_;	
	edm::EDGetTokenT<edm::ValueMap<bool> > photon_cutbasedID_decisions_loose_Token_;
	edm::EDGetTokenT<edm::ValueMap<bool> > photon_cutbasedID_decisions_medium_Token_;
	edm::EDGetTokenT<edm::ValueMap<bool> > photon_cutbasedID_decisions_tight_Token_;

	// MC
	edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puInfoToken_;

	// Calo Geometry
	const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_; // GK
	const edm::ESGetToken<PCaloGeometry, PCastorRcd> castorGeometryToken_; // GK
	// const edm::ESGetToken<RecoIdealGeometry, GEMRecoGeometryRcd> gemGeoToken_;	// commenting out to test 350 GeV MC
	const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_; // GK
	const edm::ESGetToken<Propagator, TrackingComponentsRecord> propagatorToken_; // GK

	// Noise Filters
	
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_HBHENoiseFilterToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_HBHENoiseIsoFilterToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_CSCTightHaloFilterToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_CSCTightHaloTrkMuUnvetoFilterToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_CSCTightHalo2015FilterToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_globalTightHalo2016FilterToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_globalSuperTightHalo2016FilterToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_HcalStripHaloFilterToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_hcalLaserEventFilterToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_EcalDeadCellTriggerPrimitiveFilterToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_EcalDeadCellBoundaryEnergyFilterToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_ecalBadCalibFilterToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_goodVerticesToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_eeBadScFilterToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_ecalLaserCorrFilterToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_trkPOGFiltersToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_chargedHadronTrackResolutionFilterToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_muonBadTrackFilterToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_BadChargedCandidateFilterToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_BadPFMuonFilterToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_BadChargedCandidateSummer16FilterToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_BadPFMuonSummer16FilterToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_BadPFMuonDzFilterToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_hfNoisyHitsFilterToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_trkPOG_manystripclus53XToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_trkPOG_toomanystripclus53XToken_;
	edm::EDGetTokenT<edm::HLTPathStatus> Flag_trkPOG_logErrorTooManyClustersToken_;

	// ====================================================================================
	// EDM Handles & Globals (miniAOD Inputs)
	// ====================================================================================

	// Mapping of the HLT Triggers and Filters
	vector<string> triggerPathNames;
	map<string, int> triggerPathNamesIndices;
	//vector<string> eleHLTFilterNames;
	//vector<string> muonHLTFilterNames;
	//vector<string> phoHLTFilterNames;
	//vector<string> jetHLTFilterNames;

	edm::Handle<edm::TriggerResults> triggerBits;
	edm::Handle<edm::HepMCProduct> hepMC;
	edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
	edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;

	edm::Handle<pat::METCollection> met;
	edm::Handle<edm::TriggerResults> metFilterBits;

	edm::Handle<reco::VertexCollection> vertices;
	edm::Handle<edm::Association<vector<reco::Vertex>>> primaryVertexAssociation;
	edm::Handle<edm::ValueMap<int>> primaryVertexAssociationValueMap;
	edm::Handle<double> rhoFastjetAll;
	
	// Objects
	edm::Handle<reco::GsfElectronCollection> electrons;
	edm::Handle<reco::MuonCollection> muons;
	edm::Handle<pat::TauCollection> taus;
	edm::Handle<reco::PhotonCollection> photons;
	edm::Handle<pat::JetCollection> jets;
	edm::Handle<pat::JetCollection> jetsCorr; // GK, for JECs
	edm::Handle<pat::JetCollection> jetsPuppiCorr; // GK, for JECs, PUPPI
	edm::Handle<reco::CaloJetCollection> calojets;
	edm::Handle<pat::JetCollection> LRJets;
	edm::Handle<reco::CaloJetCollection> caloLRJets;
	edm::Handle<BXVector<l1t::Jet>> l1jets;

	// Low Level Objects
	edm::Handle<edm::View<reco::Track> > tracks;
	edm::Handle<std::vector<reco::Track>> generalTracks;
	edm::Handle<vector<reco::VertexCompositePtrCandidate> > secondaryVertices;
	edm::Handle<reco::PFCandidateCollection> pfCands;
 
	// Rechits
	edm::Handle<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>> hcalRecHitsHBHE;
	edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> ebRecHits;

	// Other
	edm::Handle<edm::ValueMap<bool> > electron_cutbasedID_decisions_loose;
	edm::Handle<edm::ValueMap<bool> > electron_cutbasedID_decisions_medium;
	edm::Handle<edm::ValueMap<bool> > electron_cutbasedID_decisions_tight;	
	edm::Handle<edm::ValueMap<bool> > photon_cutbasedID_decisions_loose;
	edm::Handle<edm::ValueMap<bool> > photon_cutbasedID_decisions_medium;
	edm::Handle<edm::ValueMap<bool> > photon_cutbasedID_decisions_tight;

	// MC
	edm::Handle<reco::GenParticleCollection> genParticles;
	edm::Handle<std::vector<PileupSummaryInfo> > puInfo;

	// ----- Globals & Constants ----- //

	const reco::Vertex *PV_global;
	const reco::Vertex *PV_global_GenMatch;

	edm::ESHandle<Propagator> thePropagator_;
	edm::ESHandle<TransientTrackBuilder> theBuilder_; 

	const CaloSubdetectorGeometry *caloGeometry_EB;
	const CaloSubdetectorGeometry *caloGeometry_EE;
	const CaloSubdetectorGeometry *caloGeometry_HB;
	const CaloSubdetectorGeometry *caloGeometry_HE;
	const CaloSubdetectorGeometry *caloGeometry_HO;
	const MagneticField *magneticFieldGeometry;

	// ----- Outputs ----- // 

	TTree *output_tree;
	TH1F *NEvents;
	TH1D *sumWeights;
	TH1F *NEvents_HLT;
	//TH1D *sumScaleWeights;
	//TH1D *sumPdfWeights;
	//TH1D *sumAlphasWeights;

	// ====================================================================================
	// Variables for tree
	// ====================================================================================

	// ------------------------------------------------------------------------------------
	// Event-Level 

	// ----- Event Info ----- // 

	bool isData;
	uint runNumber;
	uint lumiNumber;
	ULong64_t eventNumber;
	uint eventTime; //in seconds, since 1970

	//float fixedGridRhoAll;
	float fixedGridRhoFastjetAll;
	/* float fixedGridRhoFastjetAllCalo;
	float fixedGridRhoFastjetCentralCalo;
	float fixedGridRhoFastjetCentralChargedPileUp;
	float fixedGridRhoFastjetCentralNeutral;*/

	// ----- PV & Candidates ----- // 

	int n_PV;
	float PV_X;
	float PV_Y;
	float PV_Z;

	// (full list of primary vertices for analysis-level vtx selection)
	int n_PVCand;
	vector<float> PVCand_X;
	vector<float> PVCand_Y;
	vector<float> PVCand_Z;
	vector<float> PVCand_LogSumPtSq;
	vector<float> PVCand_SumPx;
	vector<float> PVCand_SumPy;

	// ----- PV-Tracks ----- // 
	// (list of tracks associated with primary vertex with pt>10)

	const float PVTrack_pt_cut = 1.0;

	int   n_PVTrack;
	vector<float> PVTrack_Pt;
	vector<float> PVTrack_Eta;
	vector<float> PVTrack_Phi;

	// ----- Pileup ----- // 

	int n_BunchXing;
	vector<int> BunchXing; //[MAX_NBX];
	vector<int> nPU; //[MA_NBX];
	vector<float> nPUmean; //[MAX_NBX];

	// ----- Trigger Info ----- // 

	//vector<string> HLT_Names;
	vector<bool> HLT_Decision; //[NTriggersMAX];
	vector<int> HLT_Prescale; //[NTriggersMAX];

	//pdf weight helper
	//RazorPDFWeightsHelper pdfweightshelper;

	// ----- MET ----- //

    bool metFilterBitsToken_isValid;

	float met_Pt;
	float met_Phi;
	float met_SumEt;
 
	/*float metUncorrectedPt;
	float metUncorrectedPhi;
	float metType1Pt;
	float metType1Phi;
	float metPuppiPt;
	float metPuppiPhi;
	float metCaloPt;
	float metCaloPhi;*/

	// ------------------------------------------------------------------------------------
	// Physics Objects

	// ----- Electrons ----- // 

	int n_ele;
	// Basics
	vector<float> ele_Pt;
	vector<float> ele_Eta;
	vector<float> ele_Phi;
	vector<float> ele_E;
	// Features
	vector<float> ele_Charge;
	vector<float> ele_d0;
	vector<float> ele_dZ;
	vector<float> ele_dEta;
	vector<float> ele_dPhi;
	vector<float> ele_EtaSC;
	// ID
	vector<bool> ele_passCutBasedIDLoose;
	vector<bool> ele_passCutBasedIDMedium;
	vector<bool> ele_passCutBasedIDTight;
	// Isolation
	vector<float> ele_pileupIso;
	vector<float> ele_chargedIso;
	vector<float> ele_photonIso;
	vector<float> ele_neutralHadIso;
	// Rechits Association
	vector<uint> ele_SeedRechitID;
	vector<uint> ele_SeedRechitIndex;
	vector<vector<uint>> ele_EcalRechitIDs;
	vector<vector<uint>> ele_EcalRechitIndices;
	// HLT
	//vector<vector<bool>> ele_passHLTFilter; //[OBJECTARRAYSIZE][MAX_ElectronHLTFilters];

	// ----- Muons ----- //

	int n_muon;
	// Basics
	vector<float> muon_Pt;
	vector<float> muon_Eta;
	vector<float> muon_Phi;
	vector<float> muon_E;
	// Features
	vector<int> muon_Charge;
	vector<float> muon_d0;
	vector<float> muon_dZ;
	// ID
	vector<bool> muon_IsLoose;
	vector<bool> muon_IsMedium;
	vector<bool> muon_IsTight;
	// Isolation
	vector<float> muon_pileupIso;
	vector<float> muon_chargedIso;
	vector<float> muon_photonIso;
	vector<float> muon_neutralHadIso;
	// IP significance
	vector<double> muon_ip3dSignificance;
	vector<double> muon_dB;
	vector<double> muon_edB;
	// Rechits Association (None for now)
	// HLT
	//vector<vector<bool>> muon_passHLTFilter; //[OBJECTARRAYSIZE][MAX_ElectronHLTFilters];

	// ----- Taus ----- //

	int n_tau;
	// Basics
	vector<float> tau_Pt;
	vector<float> tau_Eta;
	vector<float> tau_Phi;
	vector<float> tau_E;
	// Features (None for now)
	// ID
	vector<bool> tau_IsLoose;
	vector<bool> tau_IsMedium;
	vector<bool> tau_IsTight;
	// Rechits Association (None for now)
	// HLT (None for now)

	// ----- Photons ----- // 

	int n_pho;
	// Basics
	vector<float> pho_Pt;
	vector<float> pho_Eta;
	vector<float> pho_Phi;
	vector<float> pho_E;
	// Features 
	vector<float> pho_sminor;
	vector<float> pho_smajor;
	vector<int> pho_convType;
	vector<float> pho_convTrkZ;  
	// ID  
	vector<bool> pho_passCutBasedIDLoose;
	vector<bool> pho_passCutBasedIDMedium;
	vector<bool> pho_passCutBasedIDTight;
	// Rechits Association
	vector<uint>  pho_SeedRechitID;
	vector<uint>  pho_SeedRechitIndex;
	vector<vector<uint>> pho_EcalRechitIDs;
	vector<vector<uint>> pho_EcalRechitIndices;
	// HLT
	//vector<vector<bool>> pho_passHLTFilter; //[OBJECTARRAYSIZE][MAX_ElectronHLTFilters];

	// ----- AK4 PF Jets ----- //

	int n_jet;
	// Basics
	vector<float> jetRaw_Pt;
	vector<float> jetRaw_E;
	vector<float> jet_Puppi_Pt;
	vector<float> jet_Puppi_E;
	vector<float> jet_Puppi_Eta;
	vector<float> jet_Puppi_Phi;
	vector<float> jet_Pt;
	vector<float> jet_E;
	vector<float> jet_Eta;
	vector<float> jet_Phi;
	vector<float> jet_Mass;
	vector<float> jet_Pt_JES_up;
	vector<float> jet_E_JES_up;
	vector<float> jet_Mass_JES_up;
	vector<float> jet_Pt_JES_down;
	vector<float> jet_E_JES_down;
	vector<float> jet_Mass_JES_down;
	vector<float> jet_Pt_noJER;
	vector<float> jet_E_noJER;
	vector<float> jet_Mass_noJER;
	vector<float> jet_Pt_JER_up;
	vector<float> jet_E_JER_up;
	vector<float> jet_Mass_JER_up;
	vector<float> jet_Pt_JER_down;
	vector<float> jet_E_JER_down;
	vector<float> jet_Mass_JER_down;
	// Features 
	vector<float> jet_JetArea;  
	vector<float> jet_ChargedHadEFrac;
	vector<float> jet_NeutralHadEFrac;
	vector<float> jet_PhoEFrac; 
	vector<float> jet_EleEFrac;
	vector<float> jet_MuonEFrac; 
	vector<int>   jet_ChargedHadMult;
	vector<int>   jet_NeutralHadMult;
	vector<int>   jet_PhoMult;
	vector<int>   jet_EleMult;
	vector<int>   jet_MuonMult;
	vector<float> jet_PtAllTracks;
	vector<float> jet_PtAllPVTracks;  
	vector<int>   jet_NVertexTracks;
	vector<int>   jet_NSelectedTracks; 
	vector<int>   jet_NSV;
	vector<int>   jet_NSVCand;
	vector<float> jet_SV_x;
	vector<float> jet_SV_y;
	vector<float> jet_SV_z;
	vector<int>   jet_SV_NTracks;
	vector<float> jet_SV_Mass;
	vector<float> jet_SV_DRJet;  
	vector<float> jet_FlightDist2D;
	vector<float> jet_FlightDist2DErr;
	vector<float> jet_FlightDist3D;
	vector<float> jet_FlightDist3DErr;  
	// ID
	vector<bool>  jet_PassIDLoose;
	vector<bool>  jet_PassIDTight;
	vector<int>   jet_PartonFlavor;
	vector<int>   jet_HadronFlavor; 
	vector<float> jet_PtAllTracks_wp;
	vector<float> jet_PtAllPVTracks_wp;
	vector<float> jet_MedianIP_wp;
	vector<float> jet_MinDeltaRAllTracks_wp;
	vector<float> jet_MinDeltaRPVTracks_wp;
    vector<float> jet_DeepCSV_prob_b;
    vector<float> jet_DeepCSV_prob_c;
    vector<float> jet_DeepCSV_prob_bb;
    vector<float> jet_DeepCSV_prob_udsg;
	// Rechits Association 
	vector<uint> jet_NTracks;
	vector<vector<uint>> jet_TrackIndices;
	vector<uint> jet_NPFCands;
	vector<vector<uint>> jet_PFCandIndices;
	vector<uint> jet_SeedHcalRechitID;
	vector<uint> jet_SeedHcalRechitIndex;
	vector<vector<uint>> jet_EcalRechitIDs;
	vector<vector<uint>> jet_EcalRechitIndices;
	vector<vector<uint>> jet_HcalRechitIDs;
	vector<vector<uint>> jet_HcalRechitIndices;
	// HLT
	//vector<vector<bool>> jet_passHLTFilter; //[OBJECTARRAYSIZE][MAX_ElectronHLTFilters];

	// ----- AK4 Calo Jets ----- //

	int n_calojet;
	// Basics
	vector<float> calojet_Pt;
	vector<float> calojet_Eta;
	vector<float> calojet_Phi;
	vector<float> calojet_E;
	vector<float> calojet_Mass;
	// Features 
	// ID
	// Rechits Association 
	vector<int>  calojet_NPFCands;
	vector<vector<int>> calojet_PFCandIndices;
	vector<uint> calojet_SeedHcalRechitID;
	vector<uint> calojet_SeedHcalRechitIndex;
	vector<vector<uint>> calojet_EcalRechitIDs;
	vector<vector<uint>> calojet_EcalRechitIndices;
	vector<vector<uint>> calojet_HcalRechitIDs;
	vector<vector<uint>> calojet_HcalRechitIndices;
	// HLT
	//vector<vector<bool>> calojet_passHLTFilter; //[OBJECTARRAYSIZE][MAX_ElectronHLTFilters];

	// ----- L1 jets ----- //
	int n_l1jet;
	vector<float> l1jet_Pt;
	vector<float> l1jet_Eta;
	vector<float> l1jet_Phi;
	vector<float> l1jet_E;
	vector<float> l1jet_hwQual;

	// ----- AK8 PF Jets ----- //

	int n_LRJet;
	// Basics
	vector<float> LRJet_Pt;
	vector<float> LRJet_Eta;
	vector<float> LRJet_Phi;
	vector<float> LRJet_E;
	vector<float> LRJet_TrimmedM;
	vector<float> LRJet_PrunedM;
	vector<float> LRJet_FilteredM;
	vector<float> LRJet_SoftDropM;
	// Features
	vector<float> LRJet_JetArea;
	vector<float> LRJet_PtAllTracks;
	vector<float> LRJet_PtAllPVTracks;  
	vector<int>   LRJet_NVertexTracks;
	vector<int>   LRJet_NSelectedTracks;   
	vector<float> LRJet_ChargedHadEFrac;
	vector<float> LRJet_NeutralHadEFrac;
	vector<float> LRJet_EleEFrac;
	vector<float> LRJet_PhoEFrac; 
	vector<int>   LRJet_EleMult;
	vector<int>   LRJet_PhoMult;
	vector<int>   LRJet_MuonMult;
	vector<int>   LRJet_NSV;
	vector<int>   LRJet_NSVCand;
	vector<float> LRJet_SV_x;
	vector<float> LRJet_SV_y;
	vector<float> LRJet_SV_z;
	vector<int>   LRJet_SV_NTracks;
	vector<float> LRJet_SV_Mass;
	vector<float> LRJet_SV_DRJet;
	vector<float> LRJet_FlightDist2D;
	vector<float> LRJet_FlightDist2DErr;
	vector<float> LRJet_FlightDist3D;
	vector<float> LRJet_FlightDist3DErr;  
	// ID
	vector<bool> LRJet_PassIDLoose;
	vector<bool> LRJet_PassIDTight;
	// Rechits Association 
	vector<int>  LRJet_NPFCands;
	vector<vector<int>> LRJet_PFCandIndices;
	vector<uint> LRJet_SeedHcalRechitID;
	vector<uint> LRJet_SeedHcalRechitIndex;
	vector<vector<uint>> LRJet_EcalRechitIDs;
	vector<vector<uint>> LRJet_EcalRechitIndices;
	vector<vector<uint>> LRJet_HcalRechitIDs;
	vector<vector<uint>> LRJet_HcalRechitIndices;
	// HLT
	//vector<vector<bool>> LRJet_passHLTFilter; //[OBJECTARRAYSIZE][MAX_ElectronHLTFilters];

	// ----- AK8 Calo Jets ----- //

	int n_caloLRJet;
	// Basics
	vector<float> caloLRJet_Pt;
	vector<float> caloLRJet_Eta;
	vector<float> caloLRJet_Phi;
	vector<float> caloLRJet_E;
	vector<float> caloLRJet_Mass;
	// Features 
	// ID
	// Rechits Association 
	vector<int>  caloLRJet_NPFCands;
	vector<vector<int>> caloLRJet_PFCandIndices;
	vector<uint> caloLRJet_SeedHcalRechitID;
	vector<uint> caloLRJet_SeedHcalRechitIndex;
	vector<vector<uint>> caloLRJet_EcalRechitIDs;
	vector<vector<uint>> caloLRJet_EcalRechitIndices;
	vector<vector<uint>> caloLRJet_HcalRechitIDs;
	vector<vector<uint>> caloLRJet_HcalRechitIndices;
	// HLT
	//vector<vector<bool>> caloLRJet_passHLTFilter;

	// ------------------------------------------------------------------------------------
	// Low-Level Objects

	// ------ Tracks ----- //

	int n_track;
	// Basics
	vector<float> track_Pt;
	vector<float> track_Eta;
	vector<float> track_Phi;
	// Features
	vector<int>   track_Charge;
	vector<int>   track_bestVertexIndex;
	vector<int>   track_nMissingInnerHits;
	vector<int>   track_nMissingOuterHits;
	vector<int>   track_nPixelHits;
	vector<int>   track_nHits;
	vector<float> track_angle;
	vector<float> track_dxyToBS;
	vector<float> track_dxyErr;
	vector<float> track_dzToPV;
	vector<float> track_dzErr;
	vector<float> track_chi2;
	vector<int>   track_ndof;
	vector<int>   track_index;

	// ------ Secondary Vertices ------ //

	int n_secVtx;
	// Basics
	vector<float> secVtx_Pt;
	vector<float> secVtx_Eta;
	vector<float> secVtx_Phi;
	// Features
	vector<int>   secVtx_Charge;
	vector<int>   secVtx_nConstituents;
	vector<float> secVtx_X;
	vector<float> secVtx_Y;
	vector<float> secVtx_Z;
	vector<float> secVtx_Distance;
	vector<float> secVtx_DistanceError;

	// ------ PFCandidates ------ //

	bool PFCand_BranchesEnabled = false;

	int n_PFCand;
	// Basics
	vector<float> PFCand_Pt;
	vector<float> PFCand_Eta;
	vector<float> PFCand_Phi;
	// Features & ID
	vector<int>   PFCand_PdgId;
	vector<int>   PFCand_TrackIndex;
	vector<int>   PFCand_GeneralTrackIndex;
	vector<int>   PFCand_PVIndex;

	// ------------------------------------------------------------------------------------
	// RecHits

	bool EcalRechit_BranchesEnabled = false;
	bool HcalRechit_BranchesEnabled = false;

	const float EcalRechit_cut = 0.5;
	const float HcalRechit_cut = 0.5;

	// ------ Ecal RecHits ------ //

	int n_ecalRechit;
	// Basics
	vector<float> ecalRechit_Eta;
	vector<float> ecalRechit_Phi;
	vector<float> ecalRechit_E;
	// Features
	vector<float> ecalRechit_X;
	vector<float> ecalRechit_Y;
	vector<float> ecalRechit_Z;	
	vector<float> ecalRechit_time;
	vector<float> ecalRechit_E_err;
	vector<float> ecalRechit_time_err;
	vector<bool> ecalRechit_kSaturatedflag;
	vector<bool> ecalRechit_kLeadingEdgeRecoveredflag;
	vector<bool> ecalRechit_kPoorRecoflag;
	vector<bool> ecalRechit_kWeirdflag;
	vector<bool> ecalRechit_kDiWeirdflag;

	// ------ Hcal RecHits ------ //

	int n_hbheRechit;
	// Basics  
	vector<float> hbheRechit_Eta;
	vector<float> hbheRechit_Phi;
	vector<float> hbheRechit_E;
	// Features  
	vector<float> hbheRechit_time;
	vector<float> hbheRechit_auxTDC;
	vector<float> hbheRechit_timeFalling;
	//vector<float> hbheRechit_cctime;
	vector<float> hbheRechit_X;
	vector<float> hbheRechit_Y;
	vector<float> hbheRechit_Z;
	vector<int  > hbheRechit_iEta;
	vector<int  > hbheRechit_iPhi;
	vector<int  > hbheRechit_depth;

	// ???
	int n_hoRechit;
	vector<float> hoRechit_Eta;
	vector<float> hoRechit_Phi;
	vector<float> hoRechit_E;
	vector<float> hoRechit_T;
	vector<float> hoRechit_X;
	vector<float> hoRechit_Y;
	vector<float> hoRechit_Z;

	// ------------------------------------------------------------------------------------
	// Flags from MET/Noise Filters (some of these applied as vetos upstream but not all)
	// Run 3 list from: https://github.com/cms-sw/cmssw/blob/a31424ec6d9989e1390981466d09b137ad068318/PhysicsTools/PatAlgos/python/slimming/metFilterPaths_cff.py

    bool Flag_HBHENoiseFilter;
    bool Flag_HBHENoiseIsoFilter;
    bool Flag_CSCTightHaloFilter;
    bool Flag_CSCTightHaloTrkMuUnvetoFilter;
    bool Flag_CSCTightHalo2015Filter;
    bool Flag_globalTightHalo2016Filter;
    bool Flag_globalSuperTightHalo2016Filter;
    bool Flag_HcalStripHaloFilter;
    bool Flag_hcalLaserEventFilter;
    bool Flag_EcalDeadCellTriggerPrimitiveFilter;
    bool Flag_EcalDeadCellBoundaryEnergyFilter;
    bool Flag_ecalBadCalibFilter;
    bool Flag_goodVertices;
    //bool Flag_trackingFailureFilter;
    bool Flag_eeBadScFilter;
    bool Flag_ecalLaserCorrFilter;
    bool Flag_trkPOGFilters;
    bool Flag_chargedHadronTrackResolutionFilter;
    bool Flag_muonBadTrackFilter;
    bool Flag_BadChargedCandidateFilter;
    bool Flag_BadPFMuonFilter;
    bool Flag_BadChargedCandidateSummer16Filter;
    bool Flag_BadPFMuonSummer16Filter;
    bool Flag_BadPFMuonDzFilter;
    bool Flag_hfNoisyHitsFilter;
    bool Flag_trkPOG_manystripclus53X;
    bool Flag_trkPOG_toomanystripclus53X;
    bool Flag_trkPOG_logErrorTooManyClusters;
    bool Flag_METFilters_2022_2023_PromptReco;

	/*bool Flag_HBHENoiseFilter;
	bool Flag_HBHETightNoiseFilter;
	bool Flag_HBHEIsoNoiseFilter;
	bool Flag_badChargedCandidateFilter;
	bool Flag_badMuonFilter;
	bool Flag_badGlobalMuonFilter;
	bool Flag_duplicateMuonFilter;
	bool Flag_globalTightHalo2016Filter;
	bool Flag_globalSuperTightHalo2016Filter;
	bool Flag_hcalLaserEventFilter;
	bool Flag_EcalDeadCellTriggerPrimitiveFilter;
	bool Flag_EcalDeadCellBoundaryEnergyFilter;
	bool Flag_goodVertices;
	bool Flag_trackingFailureFilter;
	bool Flag_eeBadScFilter;
	bool Flag_ecalLaserCorrFilter;
	bool Flag_trkPOGFilters;
	bool Flag_trkPOG_manystripclus53X;
	bool Flag_trkPOG_toomanystripclus53X;
	bool Flag_trkPOG_logErrorTooManyClusters;
	bool Flag_BadPFMuonFilter;
	bool Flag_BadChargedCandidateFilter;
	bool Flag_ecalBadCalibFilter;
	bool Flag_METFilters;

	bool Flag2_globalSuperTightHalo2016Filter;
	bool Flag2_globalTightHalo2016Filter;
	bool Flag2_goodVertices;
	bool Flag2_BadChargedCandidateFilter;
	bool Flag2_BadPFMuonFilter;
	bool Flag2_EcalDeadCellTriggerPrimitiveFilter;
	bool Flag2_HBHENoiseFilter;
	bool Flag2_HBHEIsoNoiseFilter;
	bool Flag2_ecalBadCalibFilter;
	bool Flag2_eeBadScFilter;*/

	// ------------------------------------------------------------------------------------
	// MC Truth Info

	 // Gen info

 	int n_gParticle;
	vector<int>   gParticle_Id;
	vector<int>   gParticle_Status;
	vector<int>   gParticle_ParentId;
	vector<int>   gParticle_ParentIndex;
	vector<float> gParticle_Pt;
	vector<float> gParticle_Px;
	vector<float> gParticle_Py;
	vector<float> gParticle_Pz;
	vector<float> gParticle_Eta;
	vector<float> gParticle_Phi;	
	vector<float> gParticle_E;
	vector<float> gParticle_ProdVtx_X;
	vector<float> gParticle_ProdVtx_Y;
	vector<float> gParticle_ProdVtx_Z;
	vector<float> gParticle_DecayVtx_X;
	vector<float> gParticle_DecayVtx_Y;
	vector<float> gParticle_DecayVtx_Z;

 	int n_gLLP;
	vector<float> gLLP_Pt;
	vector<float> gLLP_Eta;
	vector<float> gLLP_Phi;
	vector<float> gLLP_E;
	vector<float> gLLP_Beta;
	vector<float> gLLP_TravelTime;
	vector<float> gLLP_ProdVtx_X;
	vector<float> gLLP_ProdVtx_Y;
	vector<float> gLLP_ProdVtx_Z;
	vector<float> gLLP_DecayVtx_X;
	vector<float> gLLP_DecayVtx_Y;
	vector<float> gLLP_DecayVtx_Z;

 //MC
 /*int nGenJets;
 float genJetE[OBJECTARRAYSIZE];
 float genJetPt[OBJECTARRAYSIZE];
 float genJetEta[OBJECTARRAYSIZE];
 float genJetPhi[OBJECTARRAYSIZE];
 float genJetMET[OBJECTARRAYSIZE];
 float genMetPtCalo;
 float genMetPhiCalo;
 float genMetPtTrue;
 float genMetPhiTrue;
 float genVertexX;
 float genVertexY;
 float genVertexZ;
 float genVertexT;
 float genWeight;
 unsigned int genSignalProcessID;
 float genQScale;
 float genAlphaQCD;
 float genAlphaQED;
 string lheComments;
 vector<float> *scaleWeights;
 vector<float> *pdfWeights;
 vector<float> *alphasWeights;

 int firstPdfWeight;
 int lastPdfWeight;
 int firstAlphasWeight;
 int lastAlphasWeight;

 //gen info
 int nGenParticle;
 int gParticleMotherId[GENPARTICLEARRAYSIZE];
 int gParticleMotherIndex[GENPARTICLEARRAYSIZE];
 int gParticleId[GENPARTICLEARRAYSIZE];
 int gParticleStatus[GENPARTICLEARRAYSIZE];
 float gParticleE[GENPARTICLEARRAYSIZE];
 float gParticlePt[GENPARTICLEARRAYSIZE];
 float gParticlePx[GENPARTICLEARRAYSIZE];
 float gParticlePy[GENPARTICLEARRAYSIZE];
 float gParticlePz[GENPARTICLEARRAYSIZE];
 float gParticleEta[GENPARTICLEARRAYSIZE];
 float gParticlePhi[GENPARTICLEARRAYSIZE];

 float gParticleProdVertexX[GENPARTICLEARRAYSIZE];
 float gParticleProdVertexY[GENPARTICLEARRAYSIZE];
 float gParticleProdVertexZ[GENPARTICLEARRAYSIZE];

 float gParticleDecayVertexX[GENPARTICLEARRAYSIZE];
 float gParticleDecayVertexY[GENPARTICLEARRAYSIZE];
 float gParticleDecayVertexZ[GENPARTICLEARRAYSIZE];
 float gLLP_prod_vertex_x[LLP_ARRAY_SIZE];
 float gLLP_prod_vertex_y[LLP_ARRAY_SIZE];
 float gLLP_prod_vertex_z[LLP_ARRAY_SIZE];
 float gLLP_decay_vertex_x[LLP_ARRAY_SIZE];
 float gLLP_decay_vertex_y[LLP_ARRAY_SIZE];
 float gLLP_decay_vertex_z[LLP_ARRAY_SIZE];
 float gLLP_beta[LLP_ARRAY_SIZE];
 float gLLP_travel_time[LLP_ARRAY_SIZE];
 float gLLP_pt[LLP_ARRAY_SIZE];
 float gLLP_e[LLP_ARRAY_SIZE];
 float gLLP_eta[LLP_ARRAY_SIZE];
 float gLLP_phi[LLP_ARRAY_SIZE];
 bool gLLP_csc[LLP_ARRAY_SIZE];
 bool gLLP_dt[LLP_ARRAY_SIZE];

 float photon_travel_time[LLP_DAUGHTER_ARRAY_SIZE];
 float photon_travel_time_pv[LLP_DAUGHTER_ARRAY_SIZE];

 float gen_time[LLP_DAUGHTER_ARRAY_SIZE];
 float gen_time_pv[LLP_DAUGHTER_ARRAY_SIZE];
 float gLLP_daughter_travel_time[LLP_DAUGHTER_ARRAY_SIZE];
 int   gLLP_daughter_id[LLP_DAUGHTER_ARRAY_SIZE];
 float gLLP_daughter_pt[LLP_DAUGHTER_ARRAY_SIZE];
 float gLLP_daughter_eta[LLP_DAUGHTER_ARRAY_SIZE];
 float gLLP_daughter_phi[LLP_DAUGHTER_ARRAY_SIZE];
 float gLLP_daughter_eta_ecalcorr[LLP_DAUGHTER_ARRAY_SIZE];
 float gLLP_daughter_phi_ecalcorr[LLP_DAUGHTER_ARRAY_SIZE];
 float gLLP_daughter_e[LLP_DAUGHTER_ARRAY_SIZE];
 float gLLP_daughter_mass[LLP_DAUGHTER_ARRAY_SIZE];

 //grandaughters
 float gen_time_dau[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
 float gen_time_dau_pv[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
 float photon_travel_time_dau[LLP_DAUGHTER_ARRAY_SIZE];
 float photon_travel_time_dau_pv[LLP_DAUGHTER_ARRAY_SIZE];
 float gLLP_grandaughter_travel_time[LLP_DAUGHTER_ARRAY_SIZE];

 //bool gLLP_grandaughter_EB[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
 //bool gLLP_grandaughter_ETL[LLP_GRAND_DAUGHTER_ARRAY_SIZE];

 //float gLLP_grandaughter_photon_travel_time_EB[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
 //float gLLP_grandaughter_photon_travel_time_ETL[LLP_GRAND_DAUGHTER_ARRAY_SIZE];

 //float gLLP_grandaughter_travel_time_EB[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
 //float gLLP_grandaughter_travel_time_ETL[LLP_GRAND_DAUGHTER_ARRAY_SIZE];

 //float gen_time_grandaughter_EB[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
 //float gen_time_grandaughter_ETL[LLP_GRAND_DAUGHTER_ARRAY_SIZE];

 int   gLLP_grandaughter_id[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
 float gLLP_grandaughter_pt[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
 float gLLP_grandaughter_eta[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
 float gLLP_grandaughter_phi[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
 float gLLP_grandaughter_eta_ecalcorr[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
 float gLLP_grandaughter_phi_ecalcorr[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
 float gLLP_grandaughter_e[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
 float gLLP_grandaughter_mass[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
	*/

};

#endif
