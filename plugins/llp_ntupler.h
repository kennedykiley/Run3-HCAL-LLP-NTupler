// -*- C++ -*-
// Class:      JetNtupler
/*
Description: Base class for miniAOD analysis with CRAB
*/
//         Author:  FNAL/Caltech LLP Team
//         Created:  Wed, 13 Geb 2019 15:00:06 GMT

#ifndef RAZORTUPLIZER_H
#define RAZORTUPLIZER_H

// system include files
#include <memory>
#include <string>
#include <vector>
#include <tuple>
#include <fstream>
#include <TRandom3.h>

using namespace std;

// CMSSW framework includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"


//CMSSW package includes
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Math/interface/deltaR.h"

//ECAL Rechits
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

//ECAL conditions
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"

// Geometry
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"

#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"

//ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

//LOCAL includes
#include "cms_lpc_llp/llp_ntupler/interface/EGammaMvaEleEstimatorCSA14.h"
#include "cms_lpc_llp/llp_ntupler/interface/ElectronMVAEstimatorRun2NonTrig.h"
#include "cms_lpc_llp/llp_ntupler/interface/EGammaMvaPhotonEstimator.h"
//#include "cms_lpc_llp/llp_ntupler/interface/RazorPDFWeightsHelper.h"

//------ Array Size Constants ------//
#define OBJECTARRAYSIZE 1000
#define RECHITARRAYSIZE 2000
#define GENPARTICLEARRAYSIZE 2000
#define MAX_NPV 1000
#define MAX_NTRACK 2000
#define MAX_NPU 1000
#define MAX_NBX 1000
#define LLP_ARRAY_SIZE 2
#define LLP_DAUGHTER_ARRAY_SIZE 4

//------ Class declaration ------//

class llp_ntupler : public edm::EDAnalyzer {
public:
  //analyzer constructor and destructor
  explicit llp_ntupler(const edm::ParameterSet&);
  ~llp_ntupler();

  void loadEvent(const edm::Event& iEvent); //call at the beginning of each event to get input handles from the python config

  //enable desired output variables
  virtual void setBranches();
  void enableEventInfoBranches();
  void enablePVAllBranches();
  void enablePVTracksBranches();
  void enablePileUpBranches();
  void enableMuonBranches();
  void enableElectronBranches();
  void enableTauBranches();
  void enableIsoPFCandidateBranches();
  void enablePhotonBranches();
  void enableMuonSystemBranches();
  void enableEcalRechitBranches();
  void enableJetBranches();
  void enableCaloJetBranches();
  void enableJetAK8Branches();
  void enableMetBranches();
  void enableTriggerBranches();
  void enableMCBranches();
  void enableGenParticleBranches();


  //reset output variables in memory
  virtual void resetBranches();
  void resetEventInfoBranches();
  void resetPVAllBranches();
  void resetPVTracksBranches();
  void resetPileUpBranches();
  void resetMuonBranches();
  void resetElectronBranches();
  void resetTauBranches();
  void resetIsoPFCandidateBranches();//need to implement yet
  void resetPhotonBranches();
  void resetEcalRechitBranches();//need to implement yet
  void resetJetBranches();
  void resetCaloJetBranches();
  void resetMuonSystemBranches();
  void resetJetAK8Branches();//need to implement yet
  void resetMetBranches();//need to implement yet
  void resetMCBranches();
  void resetGenParticleBranches();
  void resetTriggerBranches();

  //------ HELPER FUNCTIONS ------//
  bool passCaloJetID( const reco::CaloJet *jetCalo, int cutLevel);
  bool passJetID( const reco::PFJet *jet, int cutLevel);
  double deltaPhi(double phi1, double phi2);
  double deltaR(double eta1, double phi1, double eta2, double phi2);
  void findTrackingVariables(const TLorentzVector &jetVec,const edm::EventSetup& iSetup,float &alphaMax,float &medianTheta2D,float &medianIP, int &nTracksPV,float &ptAllPVTracks,float &ptAllTracks,float &minDeltaRAllTracks, float &minDeltaRPVTracks);
  void jet_second_moments(std::vector<double> &et,std::vector<double> &eta,std::vector<double> &phi,double &sig1,double &sig2);
  int cscChamber(double x,double y, double z);
  int cscStation(double x,double y, double z);

  //bool fill_fat_jet(const edm::EventSetup& iSetup);
  bool fillEventInfo(const edm::Event& iEvent);
  bool fillPVAll();
  bool fillPVTracks();
  bool fillPileUp();
  bool fillMuons(const edm::Event& iEvent);
  bool fillMuonSystem(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  bool fillElectrons(const edm::Event& iEvent);
  bool fillTaus();
  bool fillPhotons(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  bool fillJets(const edm::EventSetup& iSetup);
  bool fillCaloJets(const edm::EventSetup& iSetup);
  bool fillMet(const edm::Event& iEvent);
  bool fillMC();
  bool fillGenParticles();
  bool fillTrigger(const edm::Event& iEvent);
  const reco::Candidate* findFirstMotherWithDifferentID(const reco::Candidate *particle);
  const reco::Candidate* findOriginalMotherWithSameID(const reco::Candidate *particle);

protected:
  virtual void beginJob() override;
  virtual void beginRun(const edm::Run&, const edm::EventSetup&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;


  //----- Member data ------//

  // Control Switches
  bool    isData_;
  bool    useGen_;
  bool    isFastsim_;
  bool enableTriggerInfo_;
  bool enableCaloJet_;
  bool enableGenLLPInfo_;
  bool enableEcalRechits_;
  bool readGenVertexTime_;
  bool enableAK8Jets_;
  int  llpId_;
  //bool isFourJet_;
  //bool isQCD_;

  // Mapping of the HLT Triggers and Filters
  string triggerPathNamesFile_;
  string eleHLTFilterNamesFile_;
  string muonHLTFilterNamesFile_;
  string photonHLTFilterNamesFile_;
  static const int NTriggersMAX = 602;
  string triggerPathNames[NTriggersMAX];
  static const int MAX_ElectronHLTFilters = 100;
  string eleHLTFilterNames[MAX_ElectronHLTFilters];
  static const int MAX_MuonHLTFilters = 100;
  string muonHLTFilterNames[MAX_MuonHLTFilters];
  static const int MAX_PhotonHLTFilters = 100;
  string photonHLTFilterNames[MAX_PhotonHLTFilters];

  //EDM tokens for each miniAOD input object
  edm::EDGetTokenT<reco::VertexCollection> verticesToken_;
  //edm::InputTag tracksTag_;
  //edm::InputTag trackTimeTag_;
  //edm::InputTag trackTimeResoTag_;
  edm::EDGetTokenT<edm::View<reco::Track> > tracksTag_;
  edm::EDGetTokenT<edm::ValueMap<float> > trackTimeTag_;
  edm::EDGetTokenT<edm::ValueMap<float>> trackTimeResoTag_;

  edm::EDGetTokenT<CSCSegmentCollection> cscSegmentInputToken_;
  edm::EDGetTokenT<DTRecSegment4DCollection> dtSegmentInputToken_;
  edm::EDGetTokenT<RPCRecHitCollection> rpcRecHitInputToken_;

  edm::EDGetTokenT<reco::MuonCollection> muonsToken_;
  edm::EDGetTokenT<reco::GsfElectronCollection> electronsToken_;
  edm::EDGetTokenT<reco::PFTauCollection> tausToken_;
  edm::EDGetTokenT<reco::PhotonCollection> photonsToken_;
  edm::EDGetTokenT<reco::CaloJetCollection> jetsCaloToken_;
  // edm::EDGetTokenT<reco::PFJetCollection> jetsPFToken_;
  edm::EDGetTokenT<reco::PFJetCollection> jetsToken_;
  edm::EDGetTokenT<reco::PFJetCollection> jetsPuppiToken_;
  edm::EDGetTokenT<reco::PFJetCollection> jetsAK8Token_;
  edm::EDGetTokenT<reco::PFCandidateCollection> PFCandsToken_;
  edm::EDGetTokenT<reco::PFClusterCollection> PFClustersToken_;
//  edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenParticlesToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<reco::GenJetCollection> genJetsToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBitsToken_;
  edm::EDGetTokenT<edm::HepMCProduct> hepMCToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsToken_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescalesToken_;
  edm::EDGetTokenT<reco::GenMETCollection> genMetCaloToken_;
  edm::EDGetTokenT<reco::GenMETCollection> genMetTrueToken_;
  edm::EDGetTokenT<reco::PFMETCollection> metToken_;
  edm::EDGetTokenT<reco::PFMETCollection> metNoHFToken_;
  edm::EDGetTokenT<reco::PFMETCollection> metPuppiToken_;
  edm::EDGetTokenT<edm::TriggerResults> metFilterBitsToken_;
  //edm::EDGetTokenT<bool> hbheNoiseFilterToken_;
  //edm::EDGetTokenT<bool> hbheTightNoiseFilterToken_;
  //edm::EDGetTokenT<bool> hbheIsoNoiseFilterToken_;
  //edm::EDGetTokenT<bool> badChargedCandidateFilterToken_;
  //edm::EDGetTokenT<bool> badMuonFilterToken_;
//  edm::InputTag lheRunInfoTag_;
//  edm::EDGetTokenT<LHERunInfoProduct> lheRunInfoToken_;
//  edm::EDGetTokenT<LHEEventProduct> lheInfoToken_;
  edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;
  edm::EDGetTokenT<GenLumiInfoHeader> genLumiHeaderToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puInfoToken_;
//  edm::EDGetTokenT<HcalNoiseSummary> hcalNoiseInfoToken_;
  edm::EDGetTokenT<vector<reco::VertexCompositePtrCandidate> > secondaryVerticesToken_;
  edm::EDGetTokenT<double> rhoAllToken_;
  edm::EDGetTokenT<double> rhoFastjetAllToken_;
  edm::EDGetTokenT<double> rhoFastjetAllCaloToken_;
  edm::EDGetTokenT<double> rhoFastjetCentralCaloToken_;
  edm::EDGetTokenT<double> rhoFastjetCentralChargedPileUpToken_;
  edm::EDGetTokenT<double> rhoFastjetCentralNeutralToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > ebRecHitsToken_;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > eeRecHitsToken_;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > esRecHitsToken_;
  edm::EDGetTokenT<vector<reco::CaloCluster> > ebeeClustersToken_;
  edm::EDGetTokenT<vector<reco::CaloCluster> > esClustersToken_;
  edm::EDGetTokenT<vector<reco::Conversion> > conversionsToken_;
  edm::EDGetTokenT<vector<reco::Conversion> > singleLegConversionsToken_;
  edm::EDGetTokenT<vector<reco::GsfElectronCore> > gedGsfElectronCoresToken_;
  edm::EDGetTokenT<vector<reco::PhotonCore> > gedPhotonCoresToken_;
  edm::EDGetTokenT<vector<reco::Track> > generalTrackToken_;
  //  edm::EDGetTokenT<vector<reco::SuperCluster> > superClustersToken_;
  //  edm::EDGetTokenT<vector<reco::PFCandidate> > lostTracksToken_;
  edm::EDGetTokenT<float> genParticles_t0_Token_;

  edm::EDGetTokenT<edm::ValueMap<float> > mvaGeneralPurposeValuesMapToken_;
  edm::EDGetTokenT<edm::ValueMap<int> > mvaGeneralPurposeCategoriesMapToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > mvaHZZValuesMapToken_;
  edm::EDGetTokenT<edm::ValueMap<int> > mvaHZZCategoriesMapToken_;

  //EDM handles for each miniAOD input object
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<edm::HepMCProduct> hepMC;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  edm::Handle<edm::TriggerResults> metFilterBits;
  edm::Handle<reco::VertexCollection> vertices;
  edm::Handle<edm::View<reco::Track> > tracks;
  edm::Handle<edm::ValueMap<float> > times;
  edm::Handle<edm::ValueMap<float> > timeResos;
  edm::Handle<reco::PFCandidateCollection> pfCands;
  edm::Handle<reco::PFClusterCollection> pfClusters;
  edm::Handle<reco::MuonCollection> muons;
  edm::Handle<reco::GsfElectronCollection> electrons;
  edm::Handle<reco::PhotonCollection> photons;
  edm::Handle<reco::PFTauCollection> taus;
  edm::Handle<reco::CaloJetCollection> jetsCalo;
  // edm::Handle<reco::PFJetCollection> jetsPF;
  edm::Handle<reco::PFJetCollection> jets;
  edm::Handle<reco::PFJetCollection> jetsPuppi;
  edm::Handle<reco::PFJetCollection> jetsAK8;
  edm::Handle<reco::GenMETCollection> genMetsCalo;
  edm::Handle<reco::GenMETCollection> genMetsTrue;
  //edm::Handle<reco::GenMETCollection> mets;
  edm::Handle<reco::PFMETCollection> mets;
//  edm::Handle<reco::PFMETCollection> metsNoHF;
  edm::Handle<reco::PFMETCollection> metsPuppi;
//  edm::Handle<edm::View<reco::GenParticle> > prunedGenParticles;
  edm::Handle<reco::GenParticleCollection> genParticles;
  edm::Handle<reco::GenJetCollection> genJets;
//  edm::Handle<LHEEventProduct> lheInfo;
  edm::Handle<GenEventInfoProduct> genInfo;
  edm::Handle<GenLumiInfoHeader> genLumiHeader;
  edm::Handle<std::vector<PileupSummaryInfo> > puInfo;
//  edm::Handle<HcalNoiseSummary> hcalNoiseInfo;
  //edm::Handle<bool> hbheNoiseFilter;
  //edm::Handle<bool> hbheTightNoiseFilter;
  //edm::Handle<bool> hbheIsoNoiseFilter;
  //edm::Handle<bool> badChargedCandidateFilter;
  //edm::Handle<bool> badMuonFilter;
  edm::Handle<vector<reco::VertexCompositePtrCandidate> > secondaryVertices;
  edm::Handle<double> rhoAll;
  edm::Handle<double> rhoFastjetAll;
  edm::Handle<double> rhoFastjetAllCalo;
  edm::Handle<double> rhoFastjetCentralCalo;
  edm::Handle<double> rhoFastjetCentralChargedPileUp;
  edm::Handle<double> rhoFastjetCentralNeutral;
  edm::Handle<reco::BeamSpot> beamSpot;
  edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > ebRecHits;
  edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > eeRecHits;
  edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > esRecHits;
  edm::Handle<vector<reco::CaloCluster> > ebeeClusters;
  edm::Handle<vector<reco::CaloCluster> > esClusters;
  edm::Handle<vector<reco::Conversion> > conversions;
  edm::Handle<vector<reco::Conversion>> singleLegConversions;
  edm::Handle<vector<reco::GsfElectronCore> > gedGsfElectronCores;
  edm::Handle<vector<reco::PhotonCore> > gedPhotonCores;
  edm::Handle<std::vector<reco::Track>> generalTracks;
  //  edm::Handle<vector<reco::SuperCluster> > superClusters;
  //  edm::Handle<vector<reco::PFCandidate> > lostTracks;
  edm::Handle<float> genParticles_t0;
  edm::Handle<CSCSegmentCollection> cscSegments;
  edm::Handle<DTRecSegment4DCollection> dtSegments;
  edm::Handle<RPCRecHitCollection> rpcRecHits;





  //MVAs for triggering and non-triggering electron ID
  EGammaMvaEleEstimatorCSA14* myMVATrig;
  ElectronMVAEstimatorRun2NonTrig* myMVANonTrig;
  EGammaMvaPhotonEstimator* myPhotonMVA;


  const reco::Vertex *myPV;
  const reco::Vertex *myPV_GenMatch;


  //output tree
  TTree *llpTree;
  TH1F *NEvents;
  TH1D *sumWeights;
  TH1D *sumScaleWeights;
  TH1D *sumPdfWeights;
  TH1D *sumAlphasWeights;

  //------ Variables for tree ------//

  //event info
  bool isData;
  int nPV;
  int nSlimmedSecondV;
  uint runNum;
  uint lumiNum;
  uint eventNum;
  uint eventTime;//in second, since 1970

  float pvX;
  float pvY;
  float pvZ;
  float fixedGridRhoAll;
  float fixedGridRhoFastjetAll;
  float fixedGridRhoFastjetAllCalo;
  float fixedGridRhoFastjetCentralCalo;
  float fixedGridRhoFastjetCentralChargedPileUp;
  float fixedGridRhoFastjetCentralNeutral;

  //PVAll (full list of primary vertices for analysis-level vtx selection)
 int   nPVAll;
 float pvAllX[MAX_NPV];
 float pvAllY[MAX_NPV];
 float pvAllZ[MAX_NPV];
 float pvAllLogSumPtSq[MAX_NPV];
 float pvAllSumPx[MAX_NPV];
 float pvAllSumPy[MAX_NPV];

//PV-Tacks (list of tracks associated with primary vertex with pt>10)
int   nPVTracks;
float pvTrackPt[OBJECTARRAYSIZE];
float pvTrackEta[OBJECTARRAYSIZE];
float pvTrackPhi[OBJECTARRAYSIZE];

 //PU
 int nBunchXing;
 int BunchXing[MAX_NBX];
 int nPU[MAX_NBX];
 float nPUmean[MAX_NBX];

 //Muons
 int nMuons;
 float muonE[OBJECTARRAYSIZE];
 float muonPt[OBJECTARRAYSIZE];
 float muonEta[OBJECTARRAYSIZE];
 float muonPhi[OBJECTARRAYSIZE];
 int muonCharge[OBJECTARRAYSIZE];//muon charge
 bool muonIsLoose[OBJECTARRAYSIZE];
 bool muonIsMedium[OBJECTARRAYSIZE];
 bool muonIsTight[OBJECTARRAYSIZE];
 float muon_d0[OBJECTARRAYSIZE];//transverse impact paramenter
 float muon_dZ[OBJECTARRAYSIZE];//impact parameter
 float muon_ip3d[OBJECTARRAYSIZE];//3d impact paramenter
 float muon_ip3dSignificance[OBJECTARRAYSIZE];//3d impact paramenter/error
 unsigned int muonType[OBJECTARRAYSIZE];//muonTypeBit: global, tracker, standalone
 unsigned int muonQuality[OBJECTARRAYSIZE];//muonID Quality Bits
 float muon_pileupIso[OBJECTARRAYSIZE];
 float muon_chargedIso[OBJECTARRAYSIZE];
 float muon_photonIso[OBJECTARRAYSIZE];
 float muon_neutralHadIso[OBJECTARRAYSIZE];
 float muon_ptrel[OBJECTARRAYSIZE];
 float muon_chargedMiniIso[OBJECTARRAYSIZE];
 float muon_photonAndNeutralHadronMiniIso[OBJECTARRAYSIZE];
 float muon_chargedPileupMiniIso[OBJECTARRAYSIZE];
 float muon_activityMiniIsoAnnulus[OBJECTARRAYSIZE];
 bool  muon_passSingleMuTagFilter[OBJECTARRAYSIZE];
 bool  muon_passHLTFilter[OBJECTARRAYSIZE][MAX_MuonHLTFilters];
 float muon_validFractionTrackerHits[OBJECTARRAYSIZE];
 bool  muon_isGlobal[OBJECTARRAYSIZE];
 float muon_normChi2[OBJECTARRAYSIZE];
 float muon_chi2LocalPosition[OBJECTARRAYSIZE];
 float muon_kinkFinder[OBJECTARRAYSIZE];
 float muon_segmentCompatability[OBJECTARRAYSIZE];
 bool muonIsICHEPMedium[OBJECTARRAYSIZE];

 //Electrons
 int nElectrons;
 float eleE[OBJECTARRAYSIZE];
 float elePt[OBJECTARRAYSIZE];
 float eleEta[OBJECTARRAYSIZE];
 float elePhi[OBJECTARRAYSIZE];
 float eleCharge[OBJECTARRAYSIZE];
 float eleE_SC[OBJECTARRAYSIZE];
 //float SC_ElePt[OBJECTARRAYSIZE];
 float eleEta_SC[OBJECTARRAYSIZE];
 float elePhi_SC[OBJECTARRAYSIZE];
 float eleSigmaIetaIeta[OBJECTARRAYSIZE];
 float eleFull5x5SigmaIetaIeta[OBJECTARRAYSIZE];
 float eleR9[OBJECTARRAYSIZE];
 float ele_dEta[OBJECTARRAYSIZE];
 float ele_dPhi[OBJECTARRAYSIZE];
 float ele_HoverE[OBJECTARRAYSIZE];
 float ele_d0[OBJECTARRAYSIZE];
 float ele_dZ[OBJECTARRAYSIZE];
 float ele_ip3d[OBJECTARRAYSIZE];
 float ele_ip3dSignificance[OBJECTARRAYSIZE];
 float ele_pileupIso[OBJECTARRAYSIZE];
 float ele_chargedIso[OBJECTARRAYSIZE];
 float ele_photonIso[OBJECTARRAYSIZE];
 float ele_neutralHadIso[OBJECTARRAYSIZE];
 int   ele_MissHits[OBJECTARRAYSIZE];
 bool  ele_PassConvVeto[OBJECTARRAYSIZE];
 float ele_OneOverEminusOneOverP[OBJECTARRAYSIZE];
 float ele_IDMVAGeneralPurpose[OBJECTARRAYSIZE];
 int   ele_IDMVACategoryGeneralPurpose[OBJECTARRAYSIZE];
 float ele_IDMVAHZZ[OBJECTARRAYSIZE];
 int   ele_IDMVACategoryHZZ[OBJECTARRAYSIZE];
 float ele_RegressionE[OBJECTARRAYSIZE];
 float ele_CombineP4[OBJECTARRAYSIZE];
 float ele_ptrel[OBJECTARRAYSIZE];
 float ele_chargedMiniIso[OBJECTARRAYSIZE];
 float ele_photonAndNeutralHadronMiniIso[OBJECTARRAYSIZE];
 float ele_chargedPileupMiniIso[OBJECTARRAYSIZE];
 float ele_activityMiniIsoAnnulus[OBJECTARRAYSIZE];
 bool ele_passSingleEleTagFilter[OBJECTARRAYSIZE];
 bool ele_passTPOneTagFilter[OBJECTARRAYSIZE];
 bool ele_passTPTwoTagFilter[OBJECTARRAYSIZE];
 bool ele_passTPOneProbeFilter[OBJECTARRAYSIZE];
 bool ele_passTPTwoProbeFilter[OBJECTARRAYSIZE];
 bool ele_passHLTFilter[OBJECTARRAYSIZE][MAX_ElectronHLTFilters];
 vector<vector<uint> > ele_EcalRechitID;
 vector<vector<uint> > *ele_EcalRechitIndex;
 vector<uint> ele_SeedRechitID;
 vector<uint> *ele_SeedRechitIndex;

 //Taus
 int nTaus;
 float tauE[OBJECTARRAYSIZE];
 float tauPt[OBJECTARRAYSIZE];
 float tauEta[OBJECTARRAYSIZE];
 float tauPhi[OBJECTARRAYSIZE];
 bool tau_IsLoose[OBJECTARRAYSIZE];
 bool tau_IsMedium[OBJECTARRAYSIZE];
 bool tau_IsTight[OBJECTARRAYSIZE];
 bool tau_passEleVetoLoose[OBJECTARRAYSIZE];
 bool tau_passEleVetoMedium[OBJECTARRAYSIZE];
 bool tau_passEleVetoTight[OBJECTARRAYSIZE];
 bool tau_passMuVetoLoose[OBJECTARRAYSIZE];
 bool tau_passMuVetoMedium[OBJECTARRAYSIZE];
 bool tau_passMuVetoTight[OBJECTARRAYSIZE];
 UInt_t tau_ID[OBJECTARRAYSIZE];//tauID Bits
 float tau_combinedIsoDeltaBetaCorr3Hits[OBJECTARRAYSIZE];
 float tau_chargedIsoPtSum[OBJECTARRAYSIZE];
 float tau_neutralIsoPtSum[OBJECTARRAYSIZE];
 float tau_puCorrPtSum[OBJECTARRAYSIZE];
 float tau_eleVetoMVA[OBJECTARRAYSIZE];
 int tau_eleVetoCategory[OBJECTARRAYSIZE];
 float tau_muonVetoMVA[OBJECTARRAYSIZE];
 float tau_isoMVAnewDMwLT[OBJECTARRAYSIZE];
 float tau_isoMVAnewDMwoLT[OBJECTARRAYSIZE];
 float tau_leadCandPt[OBJECTARRAYSIZE];
 int tau_leadCandID[OBJECTARRAYSIZE];
 float tau_leadChargedHadrCandPt[OBJECTARRAYSIZE];
 int tau_leadChargedHadrCandID[OBJECTARRAYSIZE];

 //IsolatedChargedPFCandidates
 int nIsoPFCandidates;
 float isoPFCandidatePt[OBJECTARRAYSIZE];
 float isoPFCandidateEta[OBJECTARRAYSIZE];
 float isoPFCandidatePhi[OBJECTARRAYSIZE];
 float isoPFCandidateIso04[OBJECTARRAYSIZE];
 float isoPFCandidateD0[OBJECTARRAYSIZE];
 int   isoPFCandidatePdgId[OBJECTARRAYSIZE];

 //Photons
 int nPhotons;
 int nPhotons_overlap;
 float phoE[OBJECTARRAYSIZE];
 float phoPt[OBJECTARRAYSIZE];
 float phoEta[OBJECTARRAYSIZE];
 float phoPhi[OBJECTARRAYSIZE];
 float phoSigmaIetaIeta[OBJECTARRAYSIZE];
 float phoFull5x5SigmaIetaIeta[OBJECTARRAYSIZE];
 float phoR9[OBJECTARRAYSIZE];
 float pho_sminor[OBJECTARRAYSIZE];
 float pho_smajor[OBJECTARRAYSIZE];
 float pho_HoverE[OBJECTARRAYSIZE];
 float pho_sumChargedHadronPtAllVertices[OBJECTARRAYSIZE][MAX_NPV];
 float pho_sumChargedHadronPt[OBJECTARRAYSIZE];
 float pho_sumNeutralHadronEt[OBJECTARRAYSIZE];
 float pho_sumPhotonEt[OBJECTARRAYSIZE];
 float pho_ecalPFClusterIso[OBJECTARRAYSIZE];
 float pho_hcalPFClusterIso[OBJECTARRAYSIZE];
 float pho_trkSumPtHollowConeDR03[OBJECTARRAYSIZE];
 float pho_sumWorstVertexChargedHadronPt[OBJECTARRAYSIZE];
 float pho_pfIsoChargedHadronIso[OBJECTARRAYSIZE];
 float pho_pfIsoChargedHadronIsoWrongVtx[OBJECTARRAYSIZE];
 float pho_pfIsoNeutralHadronIso[OBJECTARRAYSIZE];
 float pho_pfIsoPhotonIso[OBJECTARRAYSIZE];
 float pho_pfIsoModFrixione[OBJECTARRAYSIZE];
 float pho_pfIsoSumPUPt[OBJECTARRAYSIZE];
 bool  pho_isConversion[OBJECTARRAYSIZE];
 bool  pho_passEleVeto[OBJECTARRAYSIZE];
 float pho_RegressionE[OBJECTARRAYSIZE];
 float pho_RegressionEUncertainty[OBJECTARRAYSIZE];
 float pho_IDMVA[OBJECTARRAYSIZE];
 float pho_superClusterEnergy[OBJECTARRAYSIZE];
 float pho_superClusterRawEnergy[OBJECTARRAYSIZE];
 float pho_superClusterEta[OBJECTARRAYSIZE];
 float pho_superClusterPhi[OBJECTARRAYSIZE];
 float pho_superClusterX[OBJECTARRAYSIZE];
 float pho_superClusterY[OBJECTARRAYSIZE];
 float pho_superClusterZ[OBJECTARRAYSIZE];
 bool pho_hasPixelSeed[OBJECTARRAYSIZE];
 bool pho_passHLTFilter[OBJECTARRAYSIZE][MAX_PhotonHLTFilters];
 int pho_convType[OBJECTARRAYSIZE];
 float pho_convTrkZ[OBJECTARRAYSIZE];
 float pho_convTrkClusZ[OBJECTARRAYSIZE];
 float pho_vtxSumPx[OBJECTARRAYSIZE][MAX_NPV];
 float pho_vtxSumPy[OBJECTARRAYSIZE][MAX_NPV];
 bool  pho_isStandardPhoton[OBJECTARRAYSIZE];
 bool  pho_seedRecHitSwitchToGain6[OBJECTARRAYSIZE];
 bool  pho_seedRecHitSwitchToGain1[OBJECTARRAYSIZE];
 bool  pho_anyRecHitSwitchToGain6[OBJECTARRAYSIZE];
 bool  pho_anyRecHitSwitchToGain1[OBJECTARRAYSIZE];
 vector<vector<uint> > pho_EcalRechitID;
 vector<vector<uint> > *pho_EcalRechitIndex;
 vector<uint>  pho_SeedRechitID;
 vector<uint>  *pho_SeedRechitIndex;
 //extra stuff
 float pho_sumChargedHadronPt_NewPV_NoTiming[OBJECTARRAYSIZE];
 float pho_sumChargedHadronPt_NewPV_Timing50_TrkVtx[OBJECTARRAYSIZE];
 float pho_sumChargedHadronPt_NewPV_Timing80_TrkVtx[OBJECTARRAYSIZE];
 float pho_sumChargedHadronPt_NewPV_Timing100_TrkVtx[OBJECTARRAYSIZE];
 float pho_sumChargedHadronPt_NewPV_Timing120_TrkVtx[OBJECTARRAYSIZE];
 float pho_sumChargedHadronPt_NewPV_Timing50_TrkPho[OBJECTARRAYSIZE];
 float pho_sumChargedHadronPt_NewPV_Timing80_TrkPho[OBJECTARRAYSIZE];
 float pho_sumChargedHadronPt_NewPV_Timing100_TrkPho[OBJECTARRAYSIZE];
 float pho_sumChargedHadronPt_NewPV_Timing120_TrkPho[OBJECTARRAYSIZE];
 float pho_superClusterSeedX[OBJECTARRAYSIZE];
float pho_superClusterSeedY[OBJECTARRAYSIZE];
float pho_superClusterSeedZ[OBJECTARRAYSIZE];
float pho_superClusterSeedT[OBJECTARRAYSIZE];
float pho_superClusterSeedE[OBJECTARRAYSIZE];
float pho_pfClusterSeedE[OBJECTARRAYSIZE];


 //Ecal RecHits
 const float Rechit_cut = 0.5;
 int nRechits;
 float ecalRechit_Eta[RECHITARRAYSIZE];
 float ecalRechit_Phi[RECHITARRAYSIZE];
 float ecalRechit_E[RECHITARRAYSIZE];
 float ecalRechit_T[RECHITARRAYSIZE];
 float ecalRechit_E_Error[RECHITARRAYSIZE];
 float ecalRechit_T_Error[RECHITARRAYSIZE];
 bool ecalRechit_kSaturatedflag[RECHITARRAYSIZE];
 bool ecalRechit_kLeadingEdgeRecoveredflag[RECHITARRAYSIZE];
 bool ecalRechit_kPoorRecoflag[RECHITARRAYSIZE];
 bool ecalRechit_kWeirdflag[RECHITARRAYSIZE];
 bool ecalRechit_kDiWeirdflag[RECHITARRAYSIZE];

  //Muon system
  int nCsc;
  float cscPhi[OBJECTARRAYSIZE];
  float cscEta[OBJECTARRAYSIZE];
  float cscX[OBJECTARRAYSIZE];
  float cscY[OBJECTARRAYSIZE];
  float cscZ[OBJECTARRAYSIZE];
  float cscDirectionX[OBJECTARRAYSIZE];
  float cscDirectionY[OBJECTARRAYSIZE];
  float cscDirectionZ[OBJECTARRAYSIZE];
  float cscNRecHits[OBJECTARRAYSIZE];
  float cscNRecHits_flag[OBJECTARRAYSIZE];
  float cscT[OBJECTARRAYSIZE];
  float cscChi2[OBJECTARRAYSIZE];

  int           nCscClusters;
  float         cscClusterJetVetoPt[OBJECTARRAYSIZE];   
  float         cscClusterJetVetoE[OBJECTARRAYSIZE];   
  float         cscClusterCaloJetVeto[OBJECTARRAYSIZE];
  float         cscClusterMuonVetoPt[OBJECTARRAYSIZE];   
  float         cscClusterMuonVetoE[OBJECTARRAYSIZE];   
  float         cscClusterX[OBJECTARRAYSIZE];   
  float         cscClusterY[OBJECTARRAYSIZE];   
  float         cscClusterZ[OBJECTARRAYSIZE];   
  float         cscClusterTime[OBJECTARRAYSIZE];   
  float         cscClusterTimeSpread[OBJECTARRAYSIZE];
  float         cscClusterGenMuonDeltaR[OBJECTARRAYSIZE];
  float         cscClusterMajorAxis[OBJECTARRAYSIZE];
  float         cscClusterMinorAxis[OBJECTARRAYSIZE];
  float         cscClusterXSpread[OBJECTARRAYSIZE];   
  float         cscClusterYSpread[OBJECTARRAYSIZE];   
  float         cscClusterZSpread[OBJECTARRAYSIZE];   
  float         cscClusterEtaPhiSpread[OBJECTARRAYSIZE];   
  float         cscClusterEtaSpread[OBJECTARRAYSIZE];   
  float         cscClusterPhiSpread[OBJECTARRAYSIZE];   
  float         cscClusterEta[OBJECTARRAYSIZE];   
  float         cscClusterPhi[OBJECTARRAYSIZE];   
  int           cscClusterSize[OBJECTARRAYSIZE];
  int           cscClusterNSegmentChamberPlus11[OBJECTARRAYSIZE];
  int           cscClusterNSegmentChamberPlus12[OBJECTARRAYSIZE];
  int           cscClusterNSegmentChamberPlus13[OBJECTARRAYSIZE];
  int           cscClusterNSegmentChamberPlus21[OBJECTARRAYSIZE];
  int           cscClusterNSegmentChamberPlus22[OBJECTARRAYSIZE];
  int           cscClusterNSegmentChamberPlus31[OBJECTARRAYSIZE];
  int           cscClusterNSegmentChamberPlus32[OBJECTARRAYSIZE];
  int           cscClusterNSegmentChamberPlus41[OBJECTARRAYSIZE];
  int           cscClusterNSegmentChamberPlus42[OBJECTARRAYSIZE];
  int           cscClusterNSegmentChamberMinus11[OBJECTARRAYSIZE];
  int           cscClusterNSegmentChamberMinus12[OBJECTARRAYSIZE];
  int           cscClusterNSegmentChamberMinus13[OBJECTARRAYSIZE];
  int           cscClusterNSegmentChamberMinus21[OBJECTARRAYSIZE];
  int           cscClusterNSegmentChamberMinus22[OBJECTARRAYSIZE];
  int           cscClusterNSegmentChamberMinus31[OBJECTARRAYSIZE];
  int           cscClusterNSegmentChamberMinus32[OBJECTARRAYSIZE];
  int           cscClusterNSegmentChamberMinus41[OBJECTARRAYSIZE];
  int           cscClusterNSegmentChamberMinus42[OBJECTARRAYSIZE];
  float         cscClusterMe11Ratio[OBJECTARRAYSIZE];
  float         cscClusterMe12Ratio[OBJECTARRAYSIZE];
  float         cscClusterVertexR[OBJECTARRAYSIZE];   
  float         cscClusterVertexZ[OBJECTARRAYSIZE];   
  int           cscClusterVertexN[OBJECTARRAYSIZE];   
  int           cscClusterVertexN1[OBJECTARRAYSIZE];   
  int           cscClusterVertexN5[OBJECTARRAYSIZE];   
  int           cscClusterVertexN10[OBJECTARRAYSIZE];   
  int           cscClusterVertexN15[OBJECTARRAYSIZE];   
  int           cscClusterVertexN20[OBJECTARRAYSIZE];   
  float         cscClusterVertexChi2[OBJECTARRAYSIZE];   
  float         cscClusterVertexDis[OBJECTARRAYSIZE];   
  float         cscClusterMaxStationRatio[OBJECTARRAYSIZE];   
  int           cscClusterMaxStation[OBJECTARRAYSIZE];   
  int           cscClusterNStation[OBJECTARRAYSIZE];
  float         cscClusterMaxChamberRatio[OBJECTARRAYSIZE];   
  int           cscClusterMaxChamber[OBJECTARRAYSIZE];   
  int           cscClusterNChamber[OBJECTARRAYSIZE];

  int nRpc;
  float rpcPhi[OBJECTARRAYSIZE];
  float rpcEta[OBJECTARRAYSIZE];
  float rpcX[OBJECTARRAYSIZE];
  float rpcY[OBJECTARRAYSIZE];
  float rpcZ[OBJECTARRAYSIZE];
  float rpcT[OBJECTARRAYSIZE];
  float rpcTError[OBJECTARRAYSIZE];

  int nDt;
  float dtPhi[OBJECTARRAYSIZE];
  float dtEta[OBJECTARRAYSIZE];
  float dtX[OBJECTARRAYSIZE];
  float dtY[OBJECTARRAYSIZE];
  float dtZ[OBJECTARRAYSIZE];
  float dtDirX[OBJECTARRAYSIZE];
  float dtDirY[OBJECTARRAYSIZE];
  float dtDirZ[OBJECTARRAYSIZE];
  float dtT[OBJECTARRAYSIZE];
  float dtTError[OBJECTARRAYSIZE];

 //AK4 Jets
 int nJets;
 float jetE[OBJECTARRAYSIZE];
 float jetEt[OBJECTARRAYSIZE];
 float jetPt[OBJECTARRAYSIZE];
 float jetEta[OBJECTARRAYSIZE];
 float jetPhi[OBJECTARRAYSIZE];
 float jetCSV[OBJECTARRAYSIZE];
 float jetCISV[OBJECTARRAYSIZE];
 float jetProbb[OBJECTARRAYSIZE];
 float jetProbc[OBJECTARRAYSIZE];
 float jetProbudsg[OBJECTARRAYSIZE];
 float jetProbbb[OBJECTARRAYSIZE];
 float jetMass[OBJECTARRAYSIZE];
 float jetJetArea[OBJECTARRAYSIZE];
 float jetPileupE[OBJECTARRAYSIZE];
 float jetPileupId[OBJECTARRAYSIZE];
 int   jetPileupIdFlag[OBJECTARRAYSIZE];
 bool  jetPassIDLoose[OBJECTARRAYSIZE];
 bool  jetPassIDTight[OBJECTARRAYSIZE];
 bool  jetPassMuFrac[OBJECTARRAYSIZE];
 bool  jetPassEleFrac[OBJECTARRAYSIZE];
 int   jetPartonFlavor[OBJECTARRAYSIZE];
 int   jetHadronFlavor[OBJECTARRAYSIZE];
 float jetChargedEMEnergyFraction[OBJECTARRAYSIZE];
 float jetNeutralEMEnergyFraction[OBJECTARRAYSIZE];
 float jetChargedHadronEnergyFraction[OBJECTARRAYSIZE];
 float jetNeutralHadronEnergyFraction[OBJECTARRAYSIZE];
 float jetMuonEnergyFraction[OBJECTARRAYSIZE];
 float jetHOEnergyFraction[OBJECTARRAYSIZE];
 float jetHFHadronEnergyFraction[OBJECTARRAYSIZE];
 float jetHFEMEnergyFraction[OBJECTARRAYSIZE];
 float jetAllMuonPt[OBJECTARRAYSIZE];
 float jetAllMuonEta[OBJECTARRAYSIZE];
 float jetAllMuonPhi[OBJECTARRAYSIZE];
 float jetAllMuonM[OBJECTARRAYSIZE];
 float jetPtWeightedDZ[OBJECTARRAYSIZE];
 int   jetNRechits[OBJECTARRAYSIZE];
 float jetRechitE[OBJECTARRAYSIZE];
 float jetRechitT[OBJECTARRAYSIZE];
 float jetRechitT_rms[OBJECTARRAYSIZE];
 float jetRechitE_Error[OBJECTARRAYSIZE];
 float jetRechitT_Error[OBJECTARRAYSIZE];

 float jetGammaMax[OBJECTARRAYSIZE];
 float jetGammaMax_ET[OBJECTARRAYSIZE];
 float jetGammaMax_EM[OBJECTARRAYSIZE];
 float jetGammaMax_Hadronic[OBJECTARRAYSIZE];
 float jetAlphaMax[OBJECTARRAYSIZE];
 float jetBetaMax[OBJECTARRAYSIZE];

 float jetPtAllTracks[OBJECTARRAYSIZE];
 float jetPtAllPVTracks[OBJECTARRAYSIZE];
 float jetMedianTheta2D[OBJECTARRAYSIZE];
 float jetMedianIP[OBJECTARRAYSIZE];
 float jetMinDeltaRAllTracks[OBJECTARRAYSIZE];
 float jetMinDeltaRPVTracks[OBJECTARRAYSIZE];

 float jet_energy_frac[OBJECTARRAYSIZE];
 float jet_sig_et1[OBJECTARRAYSIZE];
 float jet_sig_et2[OBJECTARRAYSIZE];
 bool jet_matched[OBJECTARRAYSIZE];

 //Calo Jets
 int nCaloJets;
 float calojetE[OBJECTARRAYSIZE];
 float calojetEt[OBJECTARRAYSIZE];
 float calojetPt[OBJECTARRAYSIZE];
 float calojetEta[OBJECTARRAYSIZE];
 float calojetPhi[OBJECTARRAYSIZE];

 float calojetMass[OBJECTARRAYSIZE];
 float calojetJetArea[OBJECTARRAYSIZE];
 float calojetPileupE[OBJECTARRAYSIZE];
 float calojetPileupId[OBJECTARRAYSIZE];
 int   calojetPileupIdFlag[OBJECTARRAYSIZE];
 bool  calojetPassIDLoose[OBJECTARRAYSIZE];
 bool  calojetPassIDTight[OBJECTARRAYSIZE];

 int   calojetNRechits[OBJECTARRAYSIZE];
 float calojetRechitE[OBJECTARRAYSIZE];
 float calojetRechitT[OBJECTARRAYSIZE];
 float calojetRechitT_rms[OBJECTARRAYSIZE];

 float calojetAlphaMax[OBJECTARRAYSIZE];
 float calojetBetaMax[OBJECTARRAYSIZE];
 float calojetGammaMax[OBJECTARRAYSIZE];
 float calojetGammaMax_ET[OBJECTARRAYSIZE];
 float calojetGammaMax_EM[OBJECTARRAYSIZE];
 float calojetGammaMax_Hadronic[OBJECTARRAYSIZE];
 float calojet_EMEnergyFraction[OBJECTARRAYSIZE];
 float calojet_HadronicEnergyFraction[OBJECTARRAYSIZE];

 float calojetPtAllTracks[OBJECTARRAYSIZE];
 float calojetPtAllPVTracks[OBJECTARRAYSIZE];
 float calojetMedianTheta2D[OBJECTARRAYSIZE];
 float calojetMedianIP[OBJECTARRAYSIZE];
 float calojetMinDeltaRAllTracks[OBJECTARRAYSIZE];
 float calojetMinDeltaRPVTracks[OBJECTARRAYSIZE];

 //AK8 Jets
 int nFatJets;
 float fatJetE[OBJECTARRAYSIZE];
 float fatJetPt[OBJECTARRAYSIZE];
 float fatJetEta[OBJECTARRAYSIZE];
 float fatJetPhi[OBJECTARRAYSIZE];
 float fatJetCorrectedPt[OBJECTARRAYSIZE];
 float fatJetCorrectedEta[OBJECTARRAYSIZE];
 float fatJetCorrectedPhi[OBJECTARRAYSIZE];
 float fatJetTrimmedM[OBJECTARRAYSIZE];
 float fatJetPrunedM[OBJECTARRAYSIZE];
 float fatJetFilteredM[OBJECTARRAYSIZE];
 float fatJetSoftDropM[OBJECTARRAYSIZE];
 float fatJetCorrectedSoftDropM[OBJECTARRAYSIZE];
 float fatJetUncorrectedSoftDropM[OBJECTARRAYSIZE];
 float fatJetTau1[OBJECTARRAYSIZE];
 float fatJetTau2[OBJECTARRAYSIZE];
 float fatJetTau3[OBJECTARRAYSIZE];
 float fatJetMaxSubjetCSV[OBJECTARRAYSIZE];
 bool fatJetPassIDLoose[OBJECTARRAYSIZE];
 bool fatJetPassIDTight[OBJECTARRAYSIZE];

 //MET
 float metPt;
 float metPhi;
 float sumMET;
 float UncMETdpx;
 float UncMETdpy;
 float UncMETdSumEt;

 float metEGCleanPt;
 float metEGCleanPhi;
 float metMuEGCleanPt;
 float metMuEGCleanPhi;
 float metMuEGCleanCorrPt;
 float metMuEGCleanCorrPhi;
 float metUncorrectedPt;
 float metUncorrectedPhi;
 float metType0Pt;
 float metType0Phi;
 float metType1Pt;
 float metType1Pt_raw;
 float metType1Px;
 float metType1Py;
 float metType1Eta;
 float metType1Phi_raw;
 float metType1Phi;
 float metType0Plus1Pt;
 float metType0Plus1Phi;
 float metPtRecomputed;
 float metPhiRecomputed;
 float metNoHFPt;
 float metNoHFPhi;
 float metPuppiPt;
 float metPuppiPhi;
 float metCaloPt;
 float metCaloPhi;

 float metType1PtJetResUp;
 float metType1PtJetResDown;
 float metType1PtJetEnUp;
 float metType1PtJetEnDown;
 float metType1PtMuonEnUp;
 float metType1PtMuonEnDown;
 float metType1PtElectronEnUp;
 float metType1PtElectronEnDown;
 float metType1PtTauEnUp;
 float metType1PtTauEnDown;
 float metType1PtUnclusteredEnUp;
 float metType1PtUnclusteredEnDown;
 float metType1PtPhotonEnUp;
 float metType1PtPhotonEnDown;
 float metType1PtMETUncertaintySize;
 float metType1PtJetResUpSmear;
 float metType1PtJetResDownSmear;
 float metType1PtMETFullUncertaintySize;

 float metType1PhiJetResUp;
 float metType1PhiJetResDown;
 float metType1PhiJetEnUp;
 float metType1PhiJetEnDown;
 float metType1PhiMuonEnUp;
 float metType1PhiMuonEnDown;
 float metType1PhiElectronEnUp;
 float metType1PhiElectronEnDown;
 float metType1PhiTauEnUp;
 float metType1PhiTauEnDown;
 float metType1PhiUnclusteredEnUp;
 float metType1PhiUnclusteredEnDown;
 float metType1PhiPhotonEnUp;
 float metType1PhiPhotonEnDown;
 float metType1PhiMETUncertaintySize;
 float metType1PhiJetResUpSmear;
 float metType1PhiJetResDownSmear;
 float metType1PhiMETFullUncertaintySize;


 bool Flag_HBHENoiseFilter;
 bool Flag_HBHETightNoiseFilter;
 bool Flag_HBHEIsoNoiseFilter;
 bool Flag_badChargedCandidateFilter;
 bool Flag_badMuonFilter;
 bool Flag_badGlobalMuonFilter;
 bool Flag_duplicateMuonFilter;
 bool Flag_CSCTightHaloFilter;
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


 //MC
 int nGenJets;
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

 float photon_travel_time[LLP_DAUGHTER_ARRAY_SIZE];
 float photon_travel_time_pv[LLP_DAUGHTER_ARRAY_SIZE];

 float gen_time[LLP_DAUGHTER_ARRAY_SIZE];
 float gen_time_pv[LLP_DAUGHTER_ARRAY_SIZE];
 int   gLLP_daughter_id[LLP_DAUGHTER_ARRAY_SIZE];
 float gLLP_daughter_travel_time[LLP_DAUGHTER_ARRAY_SIZE];
 float gLLP_daughter_pt[LLP_DAUGHTER_ARRAY_SIZE];
 float gLLP_daughter_eta[LLP_DAUGHTER_ARRAY_SIZE];
 float gLLP_daughter_phi[LLP_DAUGHTER_ARRAY_SIZE];
 float gLLP_daughter_eta_ecalcorr[LLP_DAUGHTER_ARRAY_SIZE];
 float gLLP_daughter_phi_ecalcorr[LLP_DAUGHTER_ARRAY_SIZE];
 float gLLP_daughter_e[LLP_DAUGHTER_ARRAY_SIZE];
 unsigned int gLLP_daughter_match_jet_index[LLP_DAUGHTER_ARRAY_SIZE];
 float gLLP_min_delta_r_match_jet[LLP_DAUGHTER_ARRAY_SIZE];
 unsigned int gLLP_daughter_match_calojet_index[LLP_DAUGHTER_ARRAY_SIZE];
 float gLLP_min_delta_r_match_calojet[LLP_DAUGHTER_ARRAY_SIZE];


 //razor variables
 float HLTMR, HLTRSQ;

 //trigger info
 std::vector<std::string>  *nameHLT;
 bool triggerDecision[NTriggersMAX];
 int  triggerHLTPrescale[NTriggersMAX];

 const float pvTrack_pt_cut = 1.0;


};

#endif
