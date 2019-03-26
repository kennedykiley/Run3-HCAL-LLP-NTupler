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

using namespace std;

// CMSSW framework includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

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

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/TauReco/interface/PFTau.h"
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

//ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

//------ Array Size Constants ------//
#define OBJECTARRAYSIZE 2000
#define RECHITARRAYSIZE 2000
#define GENPARTICLEARRAYSIZE 2000
#define MAX_NPV 600
#define MAX_NTRACK 2000
#define MAX_NPU 300
#define MAX_NBX 40
#define LLP_ARRAY_SIZE 2
#define LLP_DAUGHTER_ARRAY_SIZE 4

//------ Class declaration ------//

class JetNtupler : public edm::EDAnalyzer {
public:
  //analyzer constructor and destructor
  explicit JetNtupler(const edm::ParameterSet&);
  ~JetNtupler();

  void loadEvent(const edm::Event& iEvent); //call at the beginning of each event to get input handles from the python config
  virtual void resetBranches();

  //enable desired output variables
  virtual void setBranches();
  void reset_event_variables();
  void reset_photon_variable();
  void reset_jet_variables();
  void reset_fat_jet_variables();
  void reset_gen_llp_variable();
  void reset_gen_jet_variable();
  void reset_qcd_variables();

  //------ HELPER FUNCTIONS ------//
  bool passJetID( const reco::PFJet *jet, int cutLevel);
  double deltaPhi(double phi1, double phi2);
  double deltaR(double eta1, double phi1, double eta2, double phi2);
  void enableFatJetBranches();
  void enableMCBranches();
  void enableGenParticleBranches();
  void enableTriggerBranches();
  void enableQCDBranches();


  bool fill_fat_jet(const edm::EventSetup& iSetup);
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
  bool isData_;
  bool isFourJet_;
  bool useGen_;
  bool isFastsim_;
  bool isQCD_;
  bool enableTriggerInfo_;
  bool enableRecHitInfo_;
  bool readGenVertexTime_;

  // Mapping of the HLT Triggers and Filters
  string triggerPathNamesFile_;
  string eleHLTFilterNamesFile_;
  string muonHLTFilterNamesFile_;
  string photonHLTFilterNamesFile_;
  static const int NTriggersMAX = 601;
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

  edm::EDGetTokenT<reco::MuonCollection> muonsToken_;
  edm::EDGetTokenT<reco::GsfElectronCollection> electronsToken_;
  edm::EDGetTokenT<reco::PFTauCollection> tausToken_;
  edm::EDGetTokenT<reco::PhotonCollection> photonsToken_;
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
  //  edm::EDGetTokenT<vector<reco::SuperCluster> > superClustersToken_;
  //  edm::EDGetTokenT<vector<reco::PFCandidate> > lostTracksToken_;
  edm::EDGetTokenT<float> genParticles_t0_Token_;


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
  edm::Handle<reco::PFJetCollection> jets;
  edm::Handle<reco::PFJetCollection> jetsPuppi;
  edm::Handle<reco::PFJetCollection> jetsAK8;
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
  //  edm::Handle<vector<reco::SuperCluster> > superClusters;
  //  edm::Handle<vector<reco::PFCandidate> > lostTracks;
  edm::Handle<float> genParticles_t0;

  const reco::Vertex *myPV;
  const reco::Vertex *myPV_GenMatch;

  //output tree
  TTree *llpTree;
  TH1F *NEvents;

  //------ Variables for tree ------//


  //AK4 Jets
  int nJets;
  float jetE[OBJECTARRAYSIZE];
  float jetPt[OBJECTARRAYSIZE];
  float jetEta[OBJECTARRAYSIZE];
  float jetPhi[OBJECTARRAYSIZE];
  float jetCISV[OBJECTARRAYSIZE];
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
  float jet_charged_hadron_multiplicity[OBJECTARRAYSIZE];
  float jet_neutral_hadron_multiplicity[OBJECTARRAYSIZE];
  float jet_photon_multiplicity[OBJECTARRAYSIZE];
  float jet_electron_multiplicity[OBJECTARRAYSIZE];
  float jet_muon_multiplicity[OBJECTARRAYSIZE];
  float jet_HF_hadron_multiplicity[OBJECTARRAYSIZE];
  float jet_HF_em_multiplicity[OBJECTARRAYSIZE];
  float jet_charged_multiplicity[OBJECTARRAYSIZE];
  float jet_neutral_multiplicity[OBJECTARRAYSIZE];
  float jetMatchedGenPt[OBJECTARRAYSIZE];
  float jetMatchedGenEta[OBJECTARRAYSIZE];
  float jetMatchedGenPhi[OBJECTARRAYSIZE];
  float jetMatchedGenMass[OBJECTARRAYSIZE];
  float jetMatchedGenTime[OBJECTARRAYSIZE];
  int   jet_n_rechits[OBJECTARRAYSIZE];
  float jet_rechit_E[OBJECTARRAYSIZE];
  float jet_rechit_T[OBJECTARRAYSIZE];
  float jet_pv_rechit_T[OBJECTARRAYSIZE];

  float jet_rechit_E_Ecut4[OBJECTARRAYSIZE];
  float jet_rechit_T_Ecut4[OBJECTARRAYSIZE];
  float jet_rechit_E_Ecut3[OBJECTARRAYSIZE];
  float jet_rechit_T_Ecut3[OBJECTARRAYSIZE];
  float jet_rechit_E_Ecut2[OBJECTARRAYSIZE];
  float jet_rechit_T_Ecut2[OBJECTARRAYSIZE];
  float jet_rechit_E_Ecut1p5[OBJECTARRAYSIZE];
  float jet_rechit_T_Ecut1p5[OBJECTARRAYSIZE];
  float jet_rechit_E_Ecut1[OBJECTARRAYSIZE];
  float jet_rechit_T_Ecut1[OBJECTARRAYSIZE];
  float jet_rechit_E_Ecut0p5[OBJECTARRAYSIZE];
  float jet_rechit_T_Ecut0p5[OBJECTARRAYSIZE];
  float jet_rechits_E[OBJECTARRAYSIZE][OBJECTARRAYSIZE];
  float jet_rechits_T[OBJECTARRAYSIZE][OBJECTARRAYSIZE];

  float jet_pv_rechit_T_Ecut4[OBJECTARRAYSIZE];
  float jet_pv_rechit_T_Ecut3[OBJECTARRAYSIZE];
  float jet_pv_rechit_T_Ecut2[OBJECTARRAYSIZE];
  float jet_pv_rechit_T_Ecut1p5[OBJECTARRAYSIZE];
  float jet_pv_rechit_T_Ecut1[OBJECTARRAYSIZE];
  float jet_pv_rechit_T_Ecut0p5[OBJECTARRAYSIZE];
  float jet_pv_rechits_T[OBJECTARRAYSIZE][OBJECTARRAYSIZE];

  //AK8 Jets
  int n_fat_Jets;
  float fat_jetE[OBJECTARRAYSIZE];
  float fat_jetPt[OBJECTARRAYSIZE];
  float fat_jetEta[OBJECTARRAYSIZE];
  float fat_jetPhi[OBJECTARRAYSIZE];
  float fat_jetCISV[OBJECTARRAYSIZE];
  float fat_jetMass[OBJECTARRAYSIZE];
  float fat_jetJetArea[OBJECTARRAYSIZE];
  float fat_jetPileupE[OBJECTARRAYSIZE];
  float fat_jetPileupId[OBJECTARRAYSIZE];
  int   fat_jetPileupIdFlag[OBJECTARRAYSIZE];
  bool  fat_jetPassIDLoose[OBJECTARRAYSIZE];
  bool  fat_jetPassIDTight[OBJECTARRAYSIZE];
  bool  fat_jetPassMuFrac[OBJECTARRAYSIZE];
  bool  fat_jetPassEleFrac[OBJECTARRAYSIZE];
  int   fat_jetPartonFlavor[OBJECTARRAYSIZE];
  int   fat_jetHadronFlavor[OBJECTARRAYSIZE];
  float fat_jetChargedEMEnergyFraction[OBJECTARRAYSIZE];
  float fat_jetNeutralEMEnergyFraction[OBJECTARRAYSIZE];
  float fat_jetChargedHadronEnergyFraction[OBJECTARRAYSIZE];
  float fat_jetNeutralHadronEnergyFraction[OBJECTARRAYSIZE];
  float fat_jet_charged_hadron_multiplicity[OBJECTARRAYSIZE];
  float fat_jet_neutral_hadron_multiplicity[OBJECTARRAYSIZE];
  float fat_jet_photon_multiplicity[OBJECTARRAYSIZE];
  float fat_jet_electron_multiplicity[OBJECTARRAYSIZE];
  float fat_jet_muon_multiplicity[OBJECTARRAYSIZE];
  float fat_jet_HF_hadron_multiplicity[OBJECTARRAYSIZE];
  float fat_jet_HF_em_multiplicity[OBJECTARRAYSIZE];
  float fat_jet_charged_multiplicity[OBJECTARRAYSIZE];
  float fat_jet_neutral_multiplicity[OBJECTARRAYSIZE];
  float fat_jetMatchedGenPt[OBJECTARRAYSIZE];
  float fat_jetMatchedGenEta[OBJECTARRAYSIZE];
  float fat_jetMatchedGenPhi[OBJECTARRAYSIZE];
  float fat_jetMatchedGenMass[OBJECTARRAYSIZE];
  float fat_jetMatchedGenTime[OBJECTARRAYSIZE];
  int   fat_jet_n_rechits[OBJECTARRAYSIZE];
  float fat_jet_rechit_E[OBJECTARRAYSIZE];
  float fat_jet_rechit_T[OBJECTARRAYSIZE];
  float fat_jet_pv_rechit_T[OBJECTARRAYSIZE];

  float fat_jet_rechit_E_Ecut4[OBJECTARRAYSIZE];
  float fat_jet_rechit_T_Ecut4[OBJECTARRAYSIZE];
  float fat_jet_rechit_E_Ecut3[OBJECTARRAYSIZE];
  float fat_jet_rechit_T_Ecut3[OBJECTARRAYSIZE];
  float fat_jet_rechit_E_Ecut2[OBJECTARRAYSIZE];
  float fat_jet_rechit_T_Ecut2[OBJECTARRAYSIZE];
  float fat_jet_rechit_E_Ecut1p5[OBJECTARRAYSIZE];
  float fat_jet_rechit_T_Ecut1p5[OBJECTARRAYSIZE];
  float fat_jet_rechit_E_Ecut1[OBJECTARRAYSIZE];
  float fat_jet_rechit_T_Ecut1[OBJECTARRAYSIZE];
  float fat_jet_rechit_E_Ecut0p5[OBJECTARRAYSIZE];
  float fat_jet_rechit_T_Ecut0p5[OBJECTARRAYSIZE];
  float fat_jet_rechits_E[OBJECTARRAYSIZE][OBJECTARRAYSIZE];
  float fat_jet_rechits_T[OBJECTARRAYSIZE][OBJECTARRAYSIZE];

  float fat_jet_pv_rechit_T_Ecut4[OBJECTARRAYSIZE];
  float fat_jet_pv_rechit_T_Ecut3[OBJECTARRAYSIZE];
  float fat_jet_pv_rechit_T_Ecut2[OBJECTARRAYSIZE];
  float fat_jet_pv_rechit_T_Ecut1p5[OBJECTARRAYSIZE];
  float fat_jet_pv_rechit_T_Ecut1[OBJECTARRAYSIZE];
  float fat_jet_pv_rechit_T_Ecut0p5[OBJECTARRAYSIZE];
  float fat_jet_pv_rechits_T[OBJECTARRAYSIZE][OBJECTARRAYSIZE];





  //All Photons Match To the Jet (Take Seed RecHit as a reference)
  Int_t                   fJetNPhotons;
  Float_t                 fJetPhotonPt[OBJECTARRAYSIZE];
  Float_t                 fJetPhotonEta[OBJECTARRAYSIZE];
  Float_t                 fJetPhotonPhi[OBJECTARRAYSIZE];
  Float_t                 fJetPhotonSeedRecHitE[OBJECTARRAYSIZE];
  Float_t                 fJetPhotonSeedRecHitEta[OBJECTARRAYSIZE];
  Float_t                 fJetPhotonSeedRecHitPhi[OBJECTARRAYSIZE];
  Float_t                 fJetPhotonSeedRecHitTime[OBJECTARRAYSIZE];

  //All RecHits matched to the Photon within the Jet
  vector<float> *fJetPhotonRecHitE;
  vector<float> *fJetPhotonRecHitEta;
  vector<float> *fJetPhotonRecHitPhi;
  vector<float> *fJetPhotonRecHitTime;

  //event info
  bool isData;
  bool isFourJet;
  bool isQCD;
  uint runNum;
  uint lumiNum;
  uint eventNum;
  float pvX;
  float pvY;
  float pvZ;
  int nPV;
  float Rho;
  int nPUmean;
  int nPU;

  //MC
int nGenJets;
float genJetE[OBJECTARRAYSIZE];
float genJetME[OBJECTARRAYSIZE];
float genJetPt[OBJECTARRAYSIZE];
float genJetEta[OBJECTARRAYSIZE];
float genJetPhi[OBJECTARRAYSIZE];
float genMetPt;
float genMetPhi;
float genVertexX;
float genVertexY;
float genVertexZ;
float genVertexT;
float genWeight;
unsigned int genSignalProcessID;
float genQScale;
float genAlphaQCD;
float genAlphaQED;
unsigned int genJet_match_jet_index[OBJECTARRAYSIZE];
float genJet_min_delta_r_match_jet[OBJECTARRAYSIZE];




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


float photon_travel_time[LLP_DAUGHTER_ARRAY_SIZE];
float photon_travel_time_pv[LLP_DAUGHTER_ARRAY_SIZE];

float gen_time[LLP_DAUGHTER_ARRAY_SIZE];
float gen_time_pv[LLP_DAUGHTER_ARRAY_SIZE];
float gLLP_daughter_travel_time[LLP_DAUGHTER_ARRAY_SIZE];
float gLLP_daughter_pt[LLP_DAUGHTER_ARRAY_SIZE];
float gLLP_daughter_eta[LLP_DAUGHTER_ARRAY_SIZE];
float gLLP_daughter_phi[LLP_DAUGHTER_ARRAY_SIZE];
float gLLP_daughter_eta_ecalcorr[LLP_DAUGHTER_ARRAY_SIZE];
float gLLP_daughter_phi_ecalcorr[LLP_DAUGHTER_ARRAY_SIZE];
float gLLP_daughter_eta_hcalcorr[LLP_DAUGHTER_ARRAY_SIZE];
float gLLP_daughter_phi_hcalcorr[LLP_DAUGHTER_ARRAY_SIZE];
float gLLP_daughter_e[LLP_DAUGHTER_ARRAY_SIZE];
unsigned int gLLP_daughter_match_genJet_index[LLP_DAUGHTER_ARRAY_SIZE];
float gLLP_min_delta_r_match_genJet[LLP_DAUGHTER_ARRAY_SIZE];
unsigned int gLLP_daughter_match_jet_index_hcal[LLP_DAUGHTER_ARRAY_SIZE];
float gLLP_min_delta_r_match_jet_hcal[LLP_DAUGHTER_ARRAY_SIZE];
unsigned int gLLP_daughter_match_jet_index_loose[LLP_DAUGHTER_ARRAY_SIZE];
float gLLP_min_delta_r_match_jet_loose[LLP_DAUGHTER_ARRAY_SIZE];
unsigned int gLLP_daughter_match_jet_index_hcal_loose[LLP_DAUGHTER_ARRAY_SIZE];
float gLLP_min_delta_r_match_jet_hcal_loose[LLP_DAUGHTER_ARRAY_SIZE];
unsigned int gLLP_daughter_match_jet_index[LLP_DAUGHTER_ARRAY_SIZE];
float gLLP_min_delta_r_match_jet[LLP_DAUGHTER_ARRAY_SIZE];
float gLLP_min_delta_r_nocorr_match_jet[LLP_DAUGHTER_ARRAY_SIZE];

// QCD
int nGenQCDParticles;
float genQCD_pt[GENPARTICLEARRAYSIZE];
float genQCD_eta[GENPARTICLEARRAYSIZE];
float genQCD_phi[GENPARTICLEARRAYSIZE];
float genQCD_e[GENPARTICLEARRAYSIZE];
unsigned int genParticleQCD_match_jet_index[GENPARTICLEARRAYSIZE];
float genParticleQCD_min_delta_r_match_jet[GENPARTICLEARRAYSIZE];


//trigger info
std::vector<string>  *nameHLT;
bool triggerDecision[NTriggersMAX];
int  triggerHLTPrescale[NTriggersMAX];
};

#endif
