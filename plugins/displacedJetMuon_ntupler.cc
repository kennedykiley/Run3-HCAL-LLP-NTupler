// -*- C++ -*-
// Class:      displacedJetMuon_ntupler
/*
  Description: Base class for miniAOD analysis with CRAB
*/
//         Author:  Cristián Peña and Si Xie.
//         Created:  Thu, 26 March 2019 15:00:06 GMT

#include "displacedJetMuon_ntupler.h"
#include <fastjet/internal/base.hh>
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"
#include "cms_lpc_llp/llp_ntupler/interface/DBSCAN.h"

struct muonCscLayers{
  int id;
  int nSimHits;
  int nRecHits;
  int nWireDigis;
  int nStripDigis;
  bool Me1112;
};
struct largest_muonCscLayers
{
  inline bool operator() (const muonCscLayers& c1, const muonCscLayers& c2){return c1.nSimHits > c2.nSimHits;}
} my_largest_muonCscLayers;

//------ Constructors and destructor ------//
displacedJetMuon_ntupler::displacedJetMuon_ntupler(const edm::ParameterSet& iConfig):
  //get inputs from config file
  isData_(iConfig.getParameter<bool> ("isData")),
  useGen_(iConfig.getParameter<bool> ("useGen")),
  isFastsim_(iConfig.getParameter<bool> ("isFastsim")),
  enableTriggerInfo_(iConfig.getParameter<bool> ("enableTriggerInfo")),
  enableGenLLPInfo_(iConfig.getParameter<bool> ("enableGenLLPInfo")),
  enableEcalRechits_(iConfig.getParameter<bool> ("enableEcalRechits")),
  readGenVertexTime_(iConfig.getParameter<bool> ("readGenVertexTime")),
  readMuonDigis_(iConfig.getParameter<bool> ("readMuonDigis")),
  triggerPathNamesFile_(iConfig.getParameter<string> ("triggerPathNamesFile")),
  eleHLTFilterNamesFile_(iConfig.getParameter<string> ("eleHLTFilterNamesFile")),
  muonHLTFilterNamesFile_(iConfig.getParameter<string> ("muonHLTFilterNamesFile")),
  photonHLTFilterNamesFile_(iConfig.getParameter<string> ("photonHLTFilterNamesFile")),
  verticesToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  tracksTag_(consumes<edm::View<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracks"))),
  cscSegmentInputToken_(consumes<CSCSegmentCollection>(edm::InputTag("cscSegments"))),
  cscRechitInputToken_(consumes<CSCRecHit2DCollection>(edm::InputTag("csc2DRecHits"))),
  dtSegmentInputToken_(consumes<DTRecSegment4DCollection>(edm::InputTag("dt4DSegments"))),
  dtCosmicSegmentInputToken_(consumes<DTRecSegment4DCollection>(edm::InputTag("dt4DCosmicSegments"))),
  dtRechitInputToken_(consumes<DTRecHitCollection>(edm::InputTag("dt1DRecHits"))),
  dtCosmicRechitInputToken_(consumes<DTRecHitCollection>(edm::InputTag("dt1DCosmicRecHits"))),
  rpcRecHitInputToken_(consumes<RPCRecHitCollection>(edm::InputTag("rpcRecHits"))),
  MuonCSCSimHitsToken_(consumes<vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("MuonCSCSimHits"))),
  MuonCSCComparatorDigiToken_(consumes<MuonDigiCollection<CSCDetId,CSCComparatorDigi>>(iConfig.getParameter<edm::InputTag>("MuonCSCComparatorDigi"))),
  MuonCSCStripDigiToken_(consumes<CSCStripDigiCollection>(iConfig.getParameter<edm::InputTag>("MuonCSCStripDigi"))),
  MuonCSCWireDigiToken_(consumes<CSCWireDigiCollection>(iConfig.getParameter<edm::InputTag>("MuonCSCWireDigi"))),
  MuonCSCStripDigiSimLinksToken_(consumes<edm::DetSetVector<StripDigiSimLink>>(iConfig.getParameter<edm::InputTag>("MuonCSCStripDigiSimLinks"))),
  MuonCSCWireDigiSimLinksToken_(consumes<edm::DetSetVector<StripDigiSimLink>>(iConfig.getParameter<edm::InputTag>("MuonCSCWireDigiSimLinks"))),
  muonsToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronsToken_(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  tausToken_(consumes<reco::PFTauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  photonsToken_(consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  jetsToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  jetsPuppiToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jetsPuppi"))),
  //jetsAK8Token_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jetsAK8"))),
  jetsAK8Token_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetsAK8"))),
  PFCandsToken_(consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  PFClustersToken_(consumes<reco::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("pfClusters"))),
  //genParticlesToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"))),
  //genParticlesToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"))),
  genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  genJetsToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets"))),
  triggerBitsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"))),
  hepMCToken_(consumes<edm::HepMCProduct>(iConfig.getParameter<edm::InputTag>("hepMC"))),
  //triggerObjectsToken_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
  //triggerPrescalesToken_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("triggerPrescales"))),
  genMetCaloToken_(consumes<reco::GenMETCollection>(iConfig.getParameter<edm::InputTag>("genMetsCalo"))),
  genMetTrueToken_(consumes<reco::GenMETCollection>(iConfig.getParameter<edm::InputTag>("genMetsTrue"))),
  metToken_(consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
//  metNoHFToken_(consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("metsNoHF"))),
  metPuppiToken_(consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("metsPuppi"))),
  metFilterBitsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterBits"))),
  globalSuperTightHalo2016FilterToken_(consumes<bool>(edm::InputTag("globalSuperTightHalo2016Filter"))),
  globalTightHalo2016FilterToken_(consumes<bool>(edm::InputTag("globalTightHalo2016Filter"))),
  BadChargedCandidateFilterToken_(consumes<bool>(edm::InputTag("BadChargedCandidateFilter"))),
  BadPFMuonFilterToken_(consumes<bool>(edm::InputTag("BadPFMuonFilter"))),
  EcalDeadCellTriggerPrimitiveFilterToken_(consumes<bool>(edm::InputTag("EcalDeadCellTriggerPrimitiveFilter"))),
  HBHENoiseFilterToken_(consumes<bool>(edm::InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"))),
  HBHEIsoNoiseFilterToken_(consumes<bool>(edm::InputTag("HBHENoiseFilterResultProducer","HBHEIsoNoiseFilterResult"))),
  ecalBadCalibFilterToken_(consumes<bool>(edm::InputTag("ecalBadCalibReducedMINIAODFilter"))),
  eeBadScFilterToken_(consumes<bool>(edm::InputTag("eeBadScFilter"))),
  primaryVertexFilterToken_(consumes<bool>(edm::InputTag("primaryVertexFilter"))),

  genInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genInfo"))),
  genLumiHeaderToken_(consumes<GenLumiInfoHeader,edm::InLumi>(edm::InputTag("generator",""))),
  puInfoToken_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("puInfo"))),
  secondaryVerticesToken_(consumes<vector<reco::VertexCompositePtrCandidate> >(iConfig.getParameter<edm::InputTag>("secondaryVertices"))),
  rhoAllToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoAll"))),
  rhoFastjetAllToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetAll"))),
  rhoFastjetAllCaloToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetAllCalo"))),
  rhoFastjetCentralCaloToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetCentralCalo"))),
  rhoFastjetCentralChargedPileUpToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetCentralChargedPileUp"))),
  rhoFastjetCentralNeutralToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetCentralNeutral"))),
  beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
  ebRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("ebRecHits"))),
  eeRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("eeRecHits"))),
  esRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("esRecHits"))),
  hcalRecHitsHOToken_(consumes<edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit>>>(edm::InputTag("horeco"))),
  hcalRecHitsHBHEToken_(consumes<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>>(edm::InputTag("hbhereco"))),
  ebeeClustersToken_(consumes<vector<reco::CaloCluster> >(iConfig.getParameter<edm::InputTag>("ebeeClusters"))),
  esClustersToken_(consumes<vector<reco::CaloCluster> >(iConfig.getParameter<edm::InputTag>("esClusters"))),
  conversionsToken_(consumes<vector<reco::Conversion> >(iConfig.getParameter<edm::InputTag>("conversions"))),
  singleLegConversionsToken_(consumes<vector<reco::Conversion> >(iConfig.getParameter<edm::InputTag>("singleLegConversions"))),
  gedGsfElectronCoresToken_(consumes<vector<reco::GsfElectronCore> >(iConfig.getParameter<edm::InputTag>("gedGsfElectronCores"))),
  gedPhotonCoresToken_(consumes<vector<reco::PhotonCore> >(iConfig.getParameter<edm::InputTag>("gedPhotonCores"))),
  generalTrackToken_(consumes<std::vector<reco::Track>>(edm::InputTag("generalTracks"))),
  generalTrackHandleToken_(consumes<edm::View<reco::Track>>(edm::InputTag("generalTracks")))
{
  //declare the TFileService for output
  edm::Service<TFileService> fs;

  //set up output tree
  displacedJetMuonTree = fs->make<TTree>("llp", "selected AOD information for llp analyses");
  NEvents = fs->make<TH1F>("NEvents",";;NEvents;",1,-0.5,0.5);
  if (!isData_) {
    sumWeights = fs->make<TH1D>("sumWeights",";;sumWeights;",1,-0.5,0.5);
    sumScaleWeights = fs->make<TH1D>("sumScaleWeights",";;sumScaleWeights;",9,-0.5,8.5);
    sumPdfWeights = fs->make<TH1D>("sumPdfWeights",";;sumPdfWeights;",100,-0.5,99.5);
    sumAlphasWeights = fs->make<TH1D>("sumAlphasWeights",";;sumAlphasWeights;",2,-0.5,1.5);
    
    sumWeights->Sumw2();
    sumScaleWeights->Sumw2();
    sumPdfWeights->Sumw2();
    sumAlphasWeights->Sumw2();
  }
  else {
    sumWeights = 0;
    sumScaleWeights = 0;
    sumPdfWeights = 0;
    sumAlphasWeights = 0;
  }


  //------------------------------------------------------------------
  //Read in HLT Trigger Path List from config file
  //------------------------------------------------------------------
  for (int i = 0; i<NTriggersMAX; ++i) triggerPathNames[i] = "";
  ifstream myfile (edm::FileInPath(triggerPathNamesFile_.c_str()).fullPath().c_str()) ;
  if (myfile.is_open())
  {
    std::string line;
    int index;
    std::string hltpathname;

    while(myfile>>index>>hltpathname)
    {
      if (index < NTriggersMAX)
      {
        triggerPathNames[index] = hltpathname;
      }
    }
    myfile.close();
  }
  else
  {
    std::cout << "ERROR!!! Could not open trigger path name file : " << edm::FileInPath(triggerPathNamesFile_.c_str()).fullPath().c_str() << "\n";
  }

  // if(enableTriggerInfo_)
  // {
  //   std::cout << "\n";
  //   std::cout << "****************** Trigger Paths Defined For Razor Ntuple ******************\n";
  //   for (int i = 0; i<NTriggersMAX; ++i)
  //   {
  //     if (triggerPathNames[i] != "") std::cout << "Trigger " << i << " " << triggerPathNames[i] << "\n";
  //   }
  //   std::cout << "****************************************************************************\n";
  //   std::cout << "\n";
  // }

  if(readGenVertexTime_) genParticles_t0_Token_ = consumes<float>(iConfig.getParameter<edm::InputTag>("genParticles_t0"));

}

displacedJetMuon_ntupler::~displacedJetMuon_ntupler()
{
};

//***********************************************
//Enable output ntuple branches
//***********************************************

void displacedJetMuon_ntupler::setBranches()
{
  enableEventInfoBranches();
  enablePVAllBranches();
  enablePileUpBranches();
  enableMuonBranches();
  enableElectronBranches();
  enableTauBranches();
  // enableIsoPFCandidateBranches();
  enablePhotonBranches();
  enableMuonSystemBranches();
  enableHORechitBranches();
  enableEcalRechitBranches();
  enableJetBranches();  
  enableJetAK8Branches();
  enableMetBranches();
  enableTriggerBranches();
  enableMCBranches();
  enableGenParticleBranches();
};

void displacedJetMuon_ntupler::enableEventInfoBranches()
{
  displacedJetMuonTree->Branch("isData", &isData, "isData/O");
  displacedJetMuonTree->Branch("nPV", &nPV, "nPV/I");
  displacedJetMuonTree->Branch("runNum", &runNum, "runNum/i");
  displacedJetMuonTree->Branch("nSlimmedSecondV", &nSlimmedSecondV, "nSlimmedSecondV/i");
  displacedJetMuonTree->Branch("lumiNum", &lumiNum, "lumiNum/i");
  displacedJetMuonTree->Branch("eventNum", &eventNum, "eventNum/i");
  displacedJetMuonTree->Branch("eventTime", &eventTime, "eventTime/i");
  displacedJetMuonTree->Branch("pvX", &pvX, "pvX/F");
  displacedJetMuonTree->Branch("pvY", &pvY, "pvY/F");
  displacedJetMuonTree->Branch("pvZ", &pvZ, "pvZ/F");
  displacedJetMuonTree->Branch("fixedGridRhoAll", &fixedGridRhoAll, "fixedGridRhoAll/F");
  displacedJetMuonTree->Branch("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll, "fixedGridRhoFastjetAll/F");
  displacedJetMuonTree->Branch("fixedGridRhoFastjetAllCalo", &fixedGridRhoFastjetAllCalo, "fixedGridRhoFastjetAllCalo/F");
  displacedJetMuonTree->Branch("fixedGridRhoFastjetCentralCalo", &fixedGridRhoFastjetCentralCalo, "fixedGridRhoFastjetCentralCalo/F");
  displacedJetMuonTree->Branch("fixedGridRhoFastjetCentralChargedPileUp", &fixedGridRhoFastjetCentralChargedPileUp, "fixedGridRhoFastjetCentralChargedPileUp/F");
  displacedJetMuonTree->Branch("fixedGridRhoFastjetCentralNeutral", &fixedGridRhoFastjetCentralNeutral, "fixedGridRhoFastjetCentralNeutral/F");
}

void displacedJetMuon_ntupler::enablePVAllBranches()
{
  displacedJetMuonTree->Branch("nPVAll", &nPVAll,"nPVAll/I");
  displacedJetMuonTree->Branch("pvAllX", pvAllX,"pvAllX[nPVAll]/F");
  displacedJetMuonTree->Branch("pvAllY", pvAllY,"pvAllY[nPVAll]/F");
  displacedJetMuonTree->Branch("pvAllZ", pvAllZ,"pvAllZ[nPVAll]/F");
  displacedJetMuonTree->Branch("pvAllLogSumPtSq", pvAllLogSumPtSq,"pvAllLogSumPtSq[nPVAll]/F");
  displacedJetMuonTree->Branch("pvAllSumPx", pvAllSumPx,"pvAllSumPx[nPVAll]/F");
  displacedJetMuonTree->Branch("pvAllSumPy", pvAllSumPy,"pvAllSumPy[nPVAll]/F");
};


void displacedJetMuon_ntupler::enablePileUpBranches()
{
  displacedJetMuonTree->Branch("nBunchXing", &nBunchXing, "nBunchXing/I");
  displacedJetMuonTree->Branch("BunchXing", BunchXing, "BunchXing[nBunchXing]/I");
  displacedJetMuonTree->Branch("nPU", nPU, "nPU[nBunchXing]/I");
  displacedJetMuonTree->Branch("nPUmean", nPUmean, "nPUmean[nBunchXing]/F");
};

void displacedJetMuon_ntupler::enableMuonBranches()
{
  displacedJetMuonTree->Branch("nMuons", &nMuons,"nMuons/I");
  displacedJetMuonTree->Branch("muonE", muonE,"muonE[nMuons]/F");
  displacedJetMuonTree->Branch("muonPt", muonPt,"muonPt[nMuons]/F");
  displacedJetMuonTree->Branch("muonEta", muonEta,"muonEta[nMuons]/F");
  displacedJetMuonTree->Branch("muonPhi", muonPhi,"muonPhi[nMuons]/F");
  displacedJetMuonTree->Branch("muonCharge", muonCharge, "muonCharge[nMuons]/I");
  displacedJetMuonTree->Branch("muonIsLoose", muonIsLoose,"muonIsLoose[nMuons]/O");
  displacedJetMuonTree->Branch("muonIsMedium", muonIsMedium,"muonIsMedium[nMuons]/O");
  displacedJetMuonTree->Branch("muonIsTight", muonIsTight,"muonIsTight[nMuons]/O");
  displacedJetMuonTree->Branch("muon_d0", muon_d0, "muon_d0[nMuons]/F");
  displacedJetMuonTree->Branch("muon_dZ", muon_dZ, "muon_dZ[nMuons]/F");
  displacedJetMuonTree->Branch("muon_ip3d", muon_ip3d, "muon_ip3d[nMuons]/F");
  displacedJetMuonTree->Branch("muon_ip3dSignificance", muon_ip3dSignificance, "muon_ip3dSignificance[nMuons]/F");
  displacedJetMuonTree->Branch("muonType", muonType, "muonType[nMuons]/i");
  displacedJetMuonTree->Branch("muonQuality", muonQuality, "muonQuality[nMuons]/i");
  displacedJetMuonTree->Branch("muon_pileupIso", muon_pileupIso, "muon_pileupIso[nMuons]/F");
  displacedJetMuonTree->Branch("muon_chargedIso", muon_chargedIso, "muon_chargedIso[nMuons]/F");
  displacedJetMuonTree->Branch("muon_photonIso", muon_photonIso, "muon_photonIso[nMuons]/F");
  displacedJetMuonTree->Branch("muon_neutralHadIso", muon_neutralHadIso, "muon_neutralHadIso[nMuons]/F");
  displacedJetMuonTree->Branch("muon_ptrel", muon_ptrel, "muon_ptrel[nMuons]/F");
  displacedJetMuonTree->Branch("muon_chargedMiniIso", muon_chargedMiniIso, "muon_chargedMiniIso[nMuons]/F");
  displacedJetMuonTree->Branch("muon_photonAndNeutralHadronMiniIso", muon_photonAndNeutralHadronMiniIso, "muon_photonAndNeutralHadronMiniIso[nMuons]/F");
  displacedJetMuonTree->Branch("muon_chargedPileupMiniIso", muon_chargedPileupMiniIso, "muon_chargedPileupMiniIso[nMuons]/F");
  displacedJetMuonTree->Branch("muon_activityMiniIsoAnnulus", muon_activityMiniIsoAnnulus, "muon_activityMiniIsoAnnulus[nMuons]/F");
  displacedJetMuonTree->Branch("muon_passSingleMuTagFilter", muon_passSingleMuTagFilter, "muon_passSingleMuTagFilter[nMuons]/O");
  displacedJetMuonTree->Branch("muon_passHLTFilter", &muon_passHLTFilter, Form("muon_passHLTFilter[nMuons][%d]/O",MAX_MuonHLTFilters));
  displacedJetMuonTree->Branch("muon_validFractionTrackerHits", muon_validFractionTrackerHits, "muon_validFractionTrackerHits[nMuons]/F");
  displacedJetMuonTree->Branch("muon_isGlobal", muon_isGlobal,"muon_isGlobal[nMuons]/O");
  displacedJetMuonTree->Branch("muon_normChi2", muon_normChi2,"muon_normChi2[nMuons]/F");
  displacedJetMuonTree->Branch("muon_chi2LocalPosition", muon_chi2LocalPosition,"muon_chi2LocalPosition[nMuons]/F");
  displacedJetMuonTree->Branch("muon_kinkFinder", muon_kinkFinder,"muon_kinkFinder[nMuons]/F");
  displacedJetMuonTree->Branch("muon_segmentCompatability", muon_segmentCompatability,"muon_segmentCompatability[nMuons]/F");
  displacedJetMuonTree->Branch("muonIsICHEPMedium", muonIsICHEPMedium,"muonIsICHEPMedium[nMuons]/O");
};

void displacedJetMuon_ntupler::enableElectronBranches()
{
  displacedJetMuonTree->Branch("nElectrons", &nElectrons,"nElectrons/I");
  displacedJetMuonTree->Branch("eleE", eleE,"eleE[nElectrons]/F");
  displacedJetMuonTree->Branch("elePt", elePt,"elePt[nElectrons]/F");
  displacedJetMuonTree->Branch("eleEta", eleEta,"eleEta[nElectrons]/F");
  displacedJetMuonTree->Branch("elePhi", elePhi,"elePhi[nElectrons]/F");
  displacedJetMuonTree->Branch("eleCharge", eleCharge, "eleCharge[nElectrons]/F");
  //displacedJetMuonTree->Branch("EleE_SC", eleE_SC,"eleE_SC[nElectrons]/F");
  displacedJetMuonTree->Branch("eleEta_SC", eleEta_SC,"eleEta_SC[nElectrons]/F");
  //displacedJetMuonTree->Branch("elePhi_SC", elePhi_SC,"elePhi_SC[nElectrons]/F");
  displacedJetMuonTree->Branch("eleSigmaIetaIeta", eleSigmaIetaIeta, "eleSigmaIetaIeta[nElectrons]/F");
  displacedJetMuonTree->Branch("eleFull5x5SigmaIetaIeta", eleFull5x5SigmaIetaIeta, "eleFull5x5SigmaIetaIeta[nElectrons]/F");
  displacedJetMuonTree->Branch("eleR9", eleR9, "eleR9[nElectrons]/F");
  displacedJetMuonTree->Branch("ele_dEta", ele_dEta, "ele_dEta[nElectrons]/F");
  displacedJetMuonTree->Branch("ele_dPhi", ele_dPhi, "ele_dPhi[nElectrons]/F");
  displacedJetMuonTree->Branch("ele_HoverE", ele_HoverE, "ele_HoverE[nElectrons]/F");
  displacedJetMuonTree->Branch("ele_d0", ele_d0, "ele_d0[nElectrons]/F");
  displacedJetMuonTree->Branch("ele_dZ", ele_dZ, "ele_dZ[nElectrons]/F");
  displacedJetMuonTree->Branch("ele_ip3d", ele_ip3d, "ele_ip3d[nElectrons]/F");
  displacedJetMuonTree->Branch("ele_ip3dSignificance", ele_ip3dSignificance, "ele_ip3dSignificance[nElectrons]/F");
  displacedJetMuonTree->Branch("ele_pileupIso", ele_pileupIso, "ele_pileupIso[nElectrons]/F");
  displacedJetMuonTree->Branch("ele_chargedIso", ele_chargedIso, "ele_chargedIso[nElectrons]/F");
  displacedJetMuonTree->Branch("ele_photonIso", ele_photonIso, "ele_photonIso[nElectrons]/F");
  displacedJetMuonTree->Branch("ele_neutralHadIso", ele_neutralHadIso, "ele_neutralHadIso[nElectrons]/F");
  displacedJetMuonTree->Branch("ele_MissHits", ele_MissHits, "ele_MissHits[nElectrons]/I");
  displacedJetMuonTree->Branch("ele_PassConvVeto", ele_PassConvVeto, "ele_PassConvVeto[nElectrons]/O");
  displacedJetMuonTree->Branch("ele_OneOverEminusOneOverP", ele_OneOverEminusOneOverP, "ele_OneOverEminusOneOverP[nElectrons]/F");
  displacedJetMuonTree->Branch("ele_IDMVAGeneralPurpose", ele_IDMVAGeneralPurpose, "ele_IDMVAGeneralPurpose[nElectrons]/F");
  displacedJetMuonTree->Branch("ele_IDMVACategoryGeneralPurpose", ele_IDMVACategoryGeneralPurpose, "ele_IDMVACategoryGeneralPurpose[nElectrons]/I");
  displacedJetMuonTree->Branch("ele_IDMVAHZZ", ele_IDMVAHZZ, "ele_IDMVAHZZ[nElectrons]/F");
  displacedJetMuonTree->Branch("ele_IDMVACategoryHZZ", ele_IDMVACategoryHZZ, "ele_IDMVACategoryHZZ[nElectrons]/I");
  displacedJetMuonTree->Branch("ele_RegressionE", ele_RegressionE, "ele_RegressionE[nElectrons]/F");
  displacedJetMuonTree->Branch("ele_CombineP4", ele_CombineP4, "ele_CombineP4[nElectrons]/F");
  displacedJetMuonTree->Branch("ele_ptrel", ele_ptrel, "ele_ptrel[nElectrons]/F");
  displacedJetMuonTree->Branch("ele_chargedMiniIso", ele_chargedMiniIso, "ele_chargedMiniIso[nElectrons]/F");
  displacedJetMuonTree->Branch("ele_photonAndNeutralHadronMiniIso", ele_photonAndNeutralHadronMiniIso, "ele_photonAndNeutralHadronMiniIso[nElectrons]/F");
  displacedJetMuonTree->Branch("ele_chargedPileupMiniIso", ele_chargedPileupMiniIso, "ele_chargedPileupMiniIso[nElectrons]/F");
  displacedJetMuonTree->Branch("ele_activityMiniIsoAnnulus", ele_activityMiniIsoAnnulus, "ele_activityMiniIsoAnnulus[nElectrons]/F");
  displacedJetMuonTree->Branch("ele_passSingleEleTagFilter", ele_passSingleEleTagFilter, "ele_passSingleEleTagFilter[nElectrons]/O");
  displacedJetMuonTree->Branch("ele_passTPOneTagFilter", ele_passTPOneTagFilter, "ele_passTPOneTagFilter[nElectrons]/O");
  displacedJetMuonTree->Branch("ele_passTPTwoTagFilter", ele_passTPTwoTagFilter, "ele_passTPTwoTagFilter[nElectrons]/O");
  displacedJetMuonTree->Branch("ele_passTPOneProbeFilter", ele_passTPOneProbeFilter, "ele_passTPOneProbeFilter[nElectrons]/O");
  displacedJetMuonTree->Branch("ele_passTPTwoProbeFilter", ele_passTPTwoProbeFilter, "ele_passTPTwoProbeFilter[nElectrons]/O");
  displacedJetMuonTree->Branch("ele_passHLTFilter", &ele_passHLTFilter, Form("ele_passHLTFilter[nElectrons][%d]/O",MAX_ElectronHLTFilters));

};

void displacedJetMuon_ntupler::enableTauBranches()
{
  displacedJetMuonTree->Branch("nTaus", &nTaus,"nTaus/I");
  displacedJetMuonTree->Branch("tauE", tauE,"tauE[nTaus]/F");
  displacedJetMuonTree->Branch("tauPt", tauPt,"tauPt[nTaus]/F");
  displacedJetMuonTree->Branch("tauEta", tauEta,"tauEta[nTaus]/F");
  displacedJetMuonTree->Branch("tauPhi", tauPhi,"tauPhi[nTaus]/F");
  displacedJetMuonTree->Branch("tau_IsLoose", tau_IsLoose, "tau_IsLoose[nTaus]/O");
  displacedJetMuonTree->Branch("tau_IsMedium", tau_IsMedium, "tau_IsMedium[nTaus]/O");
  displacedJetMuonTree->Branch("tau_IsTight", tau_IsTight, "tau_IsTight[nTaus]/O");
  displacedJetMuonTree->Branch("tau_passEleVetoLoose", tau_passEleVetoLoose, "tau_passEleVetoLoose[nTaus]/O");
  displacedJetMuonTree->Branch("tau_passEleVetoMedium", tau_passEleVetoMedium, "tau_passEleVetoMedium[nTaus]/O");
  displacedJetMuonTree->Branch("tau_passEleVetoTight", tau_passEleVetoTight, "tau_passEleVetoTight[nTaus]/O");
  displacedJetMuonTree->Branch("tau_passMuVetoLoose", tau_passMuVetoLoose, "tau_passMuVetoLoose[nTaus]/O");
  displacedJetMuonTree->Branch("tau_passMuVetoMedium", tau_passMuVetoMedium, "tau_passMuVetoMedium[nTaus]/O");
  displacedJetMuonTree->Branch("tau_passMuVetoTight", tau_passMuVetoTight, "tau_passMuVetoTight[nTaus]/O");
  displacedJetMuonTree->Branch("tau_ID", tau_ID, "tau_ID[nTaus]/i");
  displacedJetMuonTree->Branch("tau_combinedIsoDeltaBetaCorr3Hits", tau_combinedIsoDeltaBetaCorr3Hits, "tau_combinedIsoDeltaBetaCorr3Hits[nTaus]/F");
  displacedJetMuonTree->Branch("tau_chargedIsoPtSum", tau_chargedIsoPtSum, "tau_chargedIsoPtSum[nTaus]/F");
  displacedJetMuonTree->Branch("tau_neutralIsoPtSum", tau_neutralIsoPtSum, "tau_neutralIsoPtSum[nTaus]/F");
  displacedJetMuonTree->Branch("tau_puCorrPtSum", tau_puCorrPtSum, "tau_puCorrPtSum[nTaus]/F");
  displacedJetMuonTree->Branch("tau_eleVetoMVA", tau_eleVetoMVA, "tau_eleVetoMVA[nTaus]/F");
  displacedJetMuonTree->Branch("tau_eleVetoCategory", tau_eleVetoCategory, "tau_eleVetoCategory[nTaus]/I");
  displacedJetMuonTree->Branch("tau_muonVetoMVA", tau_muonVetoMVA, "tau_muonVetoMVA[nTaus]/F");
  displacedJetMuonTree->Branch("tau_isoMVAnewDMwLT", tau_isoMVAnewDMwLT, "tau_isoMVAnewDMwLT[nTaus]/F");
  displacedJetMuonTree->Branch("tau_isoMVAnewDMwoLT", tau_isoMVAnewDMwoLT, "tau_isoMVAnewDMwoLT[nTaus]/F");
  displacedJetMuonTree->Branch("tau_leadCandPt", tau_leadCandPt, "tau_leadCandPt[nTaus]/F");
  displacedJetMuonTree->Branch("tau_leadCandID", tau_leadCandID, "tau_leadCandID[nTaus]/I");
  displacedJetMuonTree->Branch("tau_leadChargedHadrCandPt", tau_leadChargedHadrCandPt, "tau_leadChargedHadrCandPt[nTaus]/F");
  displacedJetMuonTree->Branch("tau_leadChargedHadrCandID", tau_leadChargedHadrCandID, "tau_leadChargedHadrCandID[nTaus]/I");
};



void displacedJetMuon_ntupler::enablePhotonBranches()
{
  displacedJetMuonTree->Branch("nPhotons", &nPhotons,"nPhotons/I");
  displacedJetMuonTree->Branch("nPhotons_overlap", &nPhotons_overlap,"nPhotons_overlap/I");
  displacedJetMuonTree->Branch("phoE", phoE,"phoE[nPhotons]/F");
  displacedJetMuonTree->Branch("phoPt", phoPt,"phoPt[nPhotons]/F");
  displacedJetMuonTree->Branch("phoEta", phoEta,"phoEta[nPhotons]/F");
  displacedJetMuonTree->Branch("phoPhi", phoPhi,"phoPhi[nPhotons]/F");
  displacedJetMuonTree->Branch("phoSigmaIetaIeta", phoSigmaIetaIeta, "phoSigmaIetaIeta[nPhotons]/F");
  displacedJetMuonTree->Branch("phoFull5x5SigmaIetaIeta", phoFull5x5SigmaIetaIeta, "phoFull5x5SigmaIetaIeta[nPhotons]/F");
  displacedJetMuonTree->Branch("phoR9", phoR9, "phoR9[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_sminor", pho_sminor, "pho_sminor[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_smajor", pho_smajor, "pho_smajor[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_HoverE", pho_HoverE, "pho_HoverE[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_sumChargedHadronPtAllVertices", &pho_sumChargedHadronPtAllVertices,Form("pho_sumChargedHadronPtAllVertices[nPhotons][%d]/F",MAX_NPV));
  displacedJetMuonTree->Branch("pho_sumChargedHadronPt", &pho_sumChargedHadronPt, "pho_sumChargedHadronPt[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_sumNeutralHadronEt", pho_sumNeutralHadronEt, "pho_sumNeutralHadronEt[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_sumPhotonEt", pho_sumPhotonEt, "pho_sumPhotonEt[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_ecalPFClusterIso", pho_ecalPFClusterIso, "pho_ecalPFClusterIso[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_hcalPFClusterIso", pho_hcalPFClusterIso, "pho_hcalPFClusterIso[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_trkSumPtHollowConeDR03", pho_trkSumPtHollowConeDR03, "pho_trkSumPtHollowConeDR03[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_sumWorstVertexChargedHadronPt", pho_sumWorstVertexChargedHadronPt, "pho_sumWorstVertexChargedHadronPt[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_pfIsoChargedHadronIso", pho_pfIsoChargedHadronIso, "pho_pfIsoChargedHadronIso[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_pfIsoChargedHadronIsoWrongVtx", pho_pfIsoChargedHadronIsoWrongVtx, "pho_pfIsoChargedHadronIsoWrongVtx[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_pfIsoNeutralHadronIso", pho_pfIsoNeutralHadronIso, "pho_pfIsoNeutralHadronIso[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_pfIsoPhotonIso", pho_pfIsoPhotonIso, "pho_pfIsoPhotonIso[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_pfIsoModFrixione", pho_pfIsoModFrixione, "pho_pfIsoModFrixione[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_pfIsoSumPUPt", pho_pfIsoSumPUPt, "pho_pfIsoSumPUPt[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_isConversion", pho_isConversion, "pho_isConversion[nPhotons]/O");
  displacedJetMuonTree->Branch("pho_passEleVeto", pho_passEleVeto, "pho_passEleVeto[nPhotons]/O");
  displacedJetMuonTree->Branch("pho_RegressionE", pho_RegressionE, "pho_RegressionE[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_RegressionEUncertainty", pho_RegressionEUncertainty, "pho_RegressionEUncertainty[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_IDMVA", pho_IDMVA, "pho_IDMVA[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_superClusterEnergy", pho_superClusterEnergy, "pho_superClusterEnergy[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_superClusterRawEnergy", pho_superClusterRawEnergy, "pho_superClusterRawEnergy[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_superClusterEta", pho_superClusterEta, "pho_superClusterEta[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_superClusterPhi", pho_superClusterPhi, "pho_superClusterPhi[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_superClusterX", pho_superClusterX, "pho_superClusterX[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_superClusterY", pho_superClusterY, "pho_superClusterY[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_superClusterZ", pho_superClusterZ, "pho_superClusterZ[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_hasPixelSeed", pho_hasPixelSeed, "pho_hasPixelSeed[nPhotons]/O");
  displacedJetMuonTree->Branch("pho_passHLTFilter", &pho_passHLTFilter, Form("pho_passHLTFilter[nPhotons][%d]/O",MAX_PhotonHLTFilters));
  displacedJetMuonTree->Branch("pho_convType", pho_convType, "pho_convType[nPhotons]/I");
  displacedJetMuonTree->Branch("pho_convTrkZ", pho_convTrkZ, "pho_convTrkZ[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_convTrkClusZ", pho_convTrkClusZ, "pho_convTrkClusZ[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_vtxSumPx", &pho_vtxSumPx,Form("pho_vtxSumPx[nPhotons][%d]/F",MAX_NPV));
  displacedJetMuonTree->Branch("pho_vtxSumPy", &pho_vtxSumPy,Form("pho_vtxSumPy[nPhotons][%d]/F",MAX_NPV));
  displacedJetMuonTree->Branch("pho_isStandardPhoton", pho_isStandardPhoton, "pho_isStandardPhoton[nPhotons]/O");
  displacedJetMuonTree->Branch("pho_seedRecHitSwitchToGain6", pho_seedRecHitSwitchToGain6, "pho_seedRecHitSwitchToGain6[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_seedRecHitSwitchToGain1", pho_seedRecHitSwitchToGain1, "pho_seedRecHitSwitchToGain1[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_anyRecHitSwitchToGain6", pho_anyRecHitSwitchToGain6, "pho_anyRecHitSwitchToGain6[nPhotons]/F");
  displacedJetMuonTree->Branch("pho_anyRecHitSwitchToGain1", pho_anyRecHitSwitchToGain1, "pho_anyRecHitSwitchToGain1[nPhotons]/F");


};

void displacedJetMuon_ntupler::enableMuonSystemBranches()
{

    displacedJetMuonTree->Branch("nCscWireDigis",&nCscWireDigis,"nCscWireDigis/I");
    displacedJetMuonTree->Branch("nCscStripDigis",&nCscStripDigis,"nCscStripDigis/I");

    displacedJetMuonTree->Branch("nCscSeg",&nCscSeg,"nCscSeg/I");
    displacedJetMuonTree->Branch("cscSegPhi",cscSegPhi,"cscSegPhi[nCscSeg]/F");
    displacedJetMuonTree->Branch("cscSegEta",cscSegEta,"cscSegEta[nCscSeg]/F");
    displacedJetMuonTree->Branch("cscSegX",cscSegX,"cscSegX[nCscSeg]/F");
    displacedJetMuonTree->Branch("cscSegY",cscSegY,"cscSegY[nCscSeg]/F");
    displacedJetMuonTree->Branch("cscSegZ",cscSegZ,"cscSegZ[nCscSeg]/F");
    displacedJetMuonTree->Branch("cscSegT",cscSegT,"cscSegT[nCscSeg]/F");
    displacedJetMuonTree->Branch("cscSegChi2",cscSegChi2,"cscSegChi2[nCscSeg]/F");
    displacedJetMuonTree->Branch("cscSegChamber",cscSegChamber,"cscSegChamber[nCscSeg]/I");
    displacedJetMuonTree->Branch("cscSegStation",cscSegStation,"cscSegStation[nCscSeg]/I");

    displacedJetMuonTree->Branch("cscSegNRecHits",cscSegNRecHits,"cscSegNRecHits[nCscSeg]/I");

    displacedJetMuonTree->Branch("ncscRechits",&ncscRechits,"ncscRechits/I");
    displacedJetMuonTree->Branch("cscRechitsPhi",cscRechitsPhi,"cscRechitsPhi[ncscRechits]/F");
    displacedJetMuonTree->Branch("cscRechitsEta",cscRechitsEta,"cscRechitsEta[ncscRechits]/F");
    displacedJetMuonTree->Branch("cscRechitsX",cscRechitsX,"cscRechitsX[ncscRechits]/F");
    displacedJetMuonTree->Branch("cscRechitsY",cscRechitsY,"cscRechitsY[ncscRechits]/F");
    displacedJetMuonTree->Branch("cscRechitsZ",cscRechitsZ,"cscRechitsZ[ncscRechits]/F");
    displacedJetMuonTree->Branch("cscRechitsE",cscRechitsE,"cscRechitsE[ncscRechits]/F");
    displacedJetMuonTree->Branch("cscRechitsTpeak",cscRechitsTpeak,"cscRechitsTpeak[ncscRechits]/F");
    displacedJetMuonTree->Branch("cscRechitsTwire",cscRechitsTwire,"cscRechitsTwire[ncscRechits]/F");
    displacedJetMuonTree->Branch("cscRechitsQuality",cscRechitsQuality,"cscRechitsQuality[ncscRechits]/I");
    displacedJetMuonTree->Branch("cscRechitsChamber",cscRechitsChamber,"cscRechitsChamber[ncscRechits]/I");
    displacedJetMuonTree->Branch("cscRechitsStation",cscRechitsStation,"cscRechitsStation[ncscRechits]/I");

    // displacedJetMuonTree->Branch("nMatchedSimHits",nMatchedSimHits,"nMatchedSimHits[nCscRecHits]/I");
    displacedJetMuonTree->Branch("cscRechitsChannels",cscRechitsChannels,"cscRechitsChannels[ncscRechits]/I");
    displacedJetMuonTree->Branch("cscRechitsNStrips",cscRechitsNStrips,"cscRechitsNStrips[ncscRechits]/i");
    displacedJetMuonTree->Branch("cscRechitsHitWire",cscRechitsHitWire,"cscRechitsHitWire[ncscRechits]/I");
    displacedJetMuonTree->Branch("cscRechitsWGroupsBX",cscRechitsWGroupsBX,"cscRechitsWGroupsBX[ncscRechits]/I");
    displacedJetMuonTree->Branch("cscRechitsNWireGroups",cscRechitsNWireGroups,"cscRechitsNWireGroups[ncscRechits]/i");
    displacedJetMuonTree->Branch("cscRechitsDetId",cscRechitsDetId,"cscRechitsDetId[ncscRechits]/I");
    // all csc clusters


    displacedJetMuonTree->Branch("nCscRechitClusters",             &nCscRechitClusters, "nCscRechitClusters/I");
    displacedJetMuonTree->Branch("cscRechitCluster_match_cscSegCluster_minDeltaR",             cscRechitCluster_match_cscSegCluster_minDeltaR,             "cscRechitCluster_match_cscSegCluster_minDeltaR[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitCluster_match_cscSegCluster_index",             cscRechitCluster_match_cscSegCluster_index,             "cscRechitCluster_match_cscSegCluster_index[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitCluster_match_gParticle_minDeltaR",             cscRechitCluster_match_gParticle_minDeltaR,             "cscRechitCluster_match_gParticle_minDeltaR[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitCluster_match_gParticle_index",             cscRechitCluster_match_gParticle_index,             "cscRechitCluster_match_gParticle_index[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitCluster_match_gParticle_id",             cscRechitCluster_match_gParticle_id,             "cscRechitCluster_match_gParticle_id[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterX",             cscRechitClusterX,             "cscRechitClusterX[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterY",             cscRechitClusterY,             "cscRechitClusterY[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterZ",             cscRechitClusterZ,             "cscRechitClusterZ[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterTime",             cscRechitClusterTime,             "cscRechitClusterTime[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterTimeSpread",             cscRechitClusterTimeSpread,             "cscRechitClusterTimeSpread[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterGenMuonDeltaR",             cscRechitClusterGenMuonDeltaR,             "cscRechitClusterGenMuonDeltaR[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterMajorAxis",             cscRechitClusterMajorAxis,             "cscRechitClusterMajorAxis[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterMinorAxis",             cscRechitClusterMinorAxis,             "cscRechitClusterMinorAxis[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterEtaPhiSpread",             cscRechitClusterEtaPhiSpread,             "cscRechitClusterEtaPhiSpread[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterPhiSpread",             cscRechitClusterPhiSpread,             "cscRechitClusterPhiSpread[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterEtaSpread",             cscRechitClusterEtaSpread,             "cscRechitClusterEtaSpread[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterXSpread",             cscRechitClusterXSpread,             "cscRechitClusterXSpread[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterYSpread",             cscRechitClusterYSpread,             "cscRechitClusterYSpread[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterZSpread",             cscRechitClusterZSpread,             "cscRechitClusterZSpread[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterPhi",             cscRechitClusterPhi,             "cscRechitClusterPhi[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterEta",             cscRechitClusterEta,             "cscRechitClusterEta[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterJetVetoPt",       cscRechitClusterJetVetoPt,       "cscRechitClusterJetVetoPt[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterJetVetoE",        cscRechitClusterJetVetoE,        "cscRechitClusterJetVetoE[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterMuonVetoPt",      cscRechitClusterMuonVetoPt,      "cscRechitClusterMuonVetoPt[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterMuonVetoE",       cscRechitClusterMuonVetoE,       "cscRechitClusterMuonVetoE[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterCaloJetVeto",     cscRechitClusterCaloJetVeto,     "cscRechitClusterCaloJetVeto[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterSize",            cscRechitClusterSize,            "cscRechitClusterSize[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterNRechitChamberPlus11",             cscRechitClusterNRechitChamberPlus11,             "cscRechitClusterNRechitChamberPlus11[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterNRechitChamberPlus12",             cscRechitClusterNRechitChamberPlus12,             "cscRechitClusterNRechitChamberPlus12[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterNRechitChamberPlus13",             cscRechitClusterNRechitChamberPlus13,             "cscRechitClusterNRechitChamberPlus13[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterNRechitChamberPlus21",             cscRechitClusterNRechitChamberPlus21,             "cscRechitClusterNRechitChamberPlus21[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterNRechitChamberPlus22",             cscRechitClusterNRechitChamberPlus22,             "cscRechitClusterNRechitChamberPlus22[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterNRechitChamberPlus31",             cscRechitClusterNRechitChamberPlus31,             "cscRechitClusterNRechitChamberPlus31[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterNRechitChamberPlus32",             cscRechitClusterNRechitChamberPlus32,             "cscRechitClusterNRechitChamberPlus32[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterNRechitChamberPlus41",             cscRechitClusterNRechitChamberPlus41,             "cscRechitClusterNRechitChamberPlus41[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterNRechitChamberPlus42",             cscRechitClusterNRechitChamberPlus42,             "cscRechitClusterNRechitChamberPlus42[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterNRechitChamberMinus11",            cscRechitClusterNRechitChamberMinus11,            "cscRechitClusterNRechitChamberMinus11[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterNRechitChamberMinus12",            cscRechitClusterNRechitChamberMinus12,            "cscRechitClusterNRechitChamberMinus12[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterNRechitChamberMinus13",            cscRechitClusterNRechitChamberMinus13,            "cscRechitClusterNRechitChamberMinus13[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterNRechitChamberMinus21",            cscRechitClusterNRechitChamberMinus21,            "cscRechitClusterNRechitChamberMinus21[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterNRechitChamberMinus22",            cscRechitClusterNRechitChamberMinus22,            "cscRechitClusterNRechitChamberMinus22[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterNRechitChamberMinus31",            cscRechitClusterNRechitChamberMinus31,            "cscRechitClusterNRechitChamberMinus31[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterNRechitChamberMinus32",            cscRechitClusterNRechitChamberMinus32,            "cscRechitClusterNRechitChamberMinus32[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterNRechitChamberMinus41",            cscRechitClusterNRechitChamberMinus41,            "cscRechitClusterNRechitChamberMinus41[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterNRechitChamberMinus42",            cscRechitClusterNRechitChamberMinus42,            "cscRechitClusterNRechitChamberMinus42[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterMe11Ratio",             cscRechitClusterMe11Ratio,             "cscRechitClusterMe11Ratio[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterMe12Ratio",             cscRechitClusterMe12Ratio,             "cscRechitClusterMe12Ratio[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterNStation",             cscRechitClusterNStation,             "cscRechitClusterNStation[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterMaxStation",             cscRechitClusterMaxStation,             "cscRechitClusterMaxStation[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterMaxStationRatio",             cscRechitClusterMaxStationRatio,             "cscRechitClusterMaxStationRatio[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterNChamber",             cscRechitClusterNChamber,             "cscRechitClusterNChamber[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterMaxChamber",             cscRechitClusterMaxChamber,             "cscRechitClusterMaxChamber[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterMaxChamberRatio",             cscRechitClusterMaxChamberRatio,             "cscRechitClusterMaxChamberRatio[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterVertexR",             cscRechitClusterVertexR,             "cscRechitClusterVertexR[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterVertexZ",             cscRechitClusterVertexZ,             "cscRechitClusterVertexZ[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterVertexDis",             cscRechitClusterVertexDis,             "cscRechitClusterVertexDis[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterVertexChi2",             cscRechitClusterVertexChi2,             "cscRechitClusterVertexChi2[nCscRechitClusters]/F");
    displacedJetMuonTree->Branch("cscRechitClusterVertexN1",             cscRechitClusterVertexN1,             "cscRechitClusterVertexN1[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterVertexN5",             cscRechitClusterVertexN5,             "cscRechitClusterVertexN5[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterVertexN10",             cscRechitClusterVertexN10,             "cscRechitClusterVertexN10[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterVertexN15",             cscRechitClusterVertexN15,             "cscRechitClusterVertexN15[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterVertexN20",             cscRechitClusterVertexN20,             "cscRechitClusterVertexN20[nCscRechitClusters]/I");
    displacedJetMuonTree->Branch("cscRechitClusterVertexN",             cscRechitClusterVertexN,             "cscRechitClusterVertexN[nCscRechitClusters]/I");

    //segment Clusters
    // all csc clusters
    displacedJetMuonTree->Branch("nCscSegClusters",             &nCscSegClusters, "nCscSegClusters/I");
    displacedJetMuonTree->Branch("cscSegCluster_match_gParticle_minDeltaR",             cscSegCluster_match_gParticle_minDeltaR,             "cscSegCluster_match_gParticle_minDeltaR[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegCluster_match_gParticle_index",             cscSegCluster_match_gParticle_index,             "cscSegCluster_match_gParticle_index[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegCluster_match_gParticle_id",             cscSegCluster_match_gParticle_id,             "cscSegCluster_match_gParticle_id[nCscSegClusters]/I");

    displacedJetMuonTree->Branch("cscSegClusterX",             cscSegClusterX,             "cscSegClusterX[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterY",             cscSegClusterY,             "cscSegClusterY[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterZ",             cscSegClusterZ,             "cscSegClusterZ[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterTime",             cscSegClusterTime,             "cscSegClusterTime[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterTimeSpread",             cscSegClusterTimeSpread,             "cscSegClusterTimeSpread[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterGenMuonDeltaR",             cscSegClusterGenMuonDeltaR,             "cscSegClusterGenMuonDeltaR[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterMajorAxis",             cscSegClusterMajorAxis,             "cscSegClusterMajorAxis[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterMinorAxis",             cscSegClusterMinorAxis,             "cscSegClusterMinorAxis[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterEtaPhiSpread",             cscSegClusterEtaPhiSpread,             "cscSegClusterEtaPhiSpread[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterPhiSpread",             cscSegClusterPhiSpread,             "cscSegClusterPhiSpread[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterEtaSpread",             cscSegClusterEtaSpread,             "cscSegClusterEtaSpread[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterXSpread",             cscSegClusterXSpread,             "cscSegClusterXSpread[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterYSpread",             cscSegClusterYSpread,             "cscSegClusterYSpread[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterZSpread",             cscSegClusterZSpread,             "cscSegClusterZSpread[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterPhi",             cscSegClusterPhi,             "cscSegClusterPhi[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterEta",             cscSegClusterEta,             "cscSegClusterEta[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterJetVetoPt",       cscSegClusterJetVetoPt,       "cscSegClusterJetVetoPt[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterJetVetoE",        cscSegClusterJetVetoE,        "cscSegClusterJetVetoE[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterMuonVetoPt",      cscSegClusterMuonVetoPt,      "cscSegClusterMuonVetoPt[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterMuonVetoE",       cscSegClusterMuonVetoE,       "cscSegClusterMuonVetoE[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterCaloJetVeto",     cscSegClusterCaloJetVeto,     "cscSegClusterCaloJetVeto[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterSize",            cscSegClusterSize,            "cscSegClusterSize[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterNSegmentChamberPlus11",             cscSegClusterNSegmentChamberPlus11,             "cscSegClusterNSegmentChamberPlus11[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterNSegmentChamberPlus12",             cscSegClusterNSegmentChamberPlus12,             "cscSegClusterNSegmentChamberPlus12[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterNSegmentChamberPlus13",             cscSegClusterNSegmentChamberPlus13,             "cscSegClusterNSegmentChamberPlus13[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterNSegmentChamberPlus21",             cscSegClusterNSegmentChamberPlus21,             "cscSegClusterNSegmentChamberPlus21[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterNSegmentChamberPlus22",             cscSegClusterNSegmentChamberPlus22,             "cscSegClusterNSegmentChamberPlus22[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterNSegmentChamberPlus31",             cscSegClusterNSegmentChamberPlus31,             "cscSegClusterNSegmentChamberPlus31[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterNSegmentChamberPlus32",             cscSegClusterNSegmentChamberPlus32,             "cscSegClusterNSegmentChamberPlus32[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterNSegmentChamberPlus41",             cscSegClusterNSegmentChamberPlus41,             "cscSegClusterNSegmentChamberPlus41[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterNSegmentChamberPlus42",             cscSegClusterNSegmentChamberPlus42,             "cscSegClusterNSegmentChamberPlus42[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterNSegmentChamberMinus11",            cscSegClusterNSegmentChamberMinus11,            "cscSegClusterNSegmentChamberMinus11[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterNSegmentChamberMinus12",            cscSegClusterNSegmentChamberMinus12,            "cscSegClusterNSegmentChamberMinus12[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterNSegmentChamberMinus13",            cscSegClusterNSegmentChamberMinus13,            "cscSegClusterNSegmentChamberMinus13[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterNSegmentChamberMinus21",            cscSegClusterNSegmentChamberMinus21,            "cscSegClusterNSegmentChamberMinus21[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterNSegmentChamberMinus22",            cscSegClusterNSegmentChamberMinus22,            "cscSegClusterNSegmentChamberMinus22[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterNSegmentChamberMinus31",            cscSegClusterNSegmentChamberMinus31,            "cscSegClusterNSegmentChamberMinus31[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterNSegmentChamberMinus32",            cscSegClusterNSegmentChamberMinus32,            "cscSegClusterNSegmentChamberMinus32[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterNSegmentChamberMinus41",            cscSegClusterNSegmentChamberMinus41,            "cscSegClusterNSegmentChamberMinus41[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterNSegmentChamberMinus42",            cscSegClusterNSegmentChamberMinus42,            "cscSegClusterNSegmentChamberMinus42[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterMe11Ratio",             cscSegClusterMe11Ratio,             "cscSegClusterMe11Ratio[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterMe12Ratio",             cscSegClusterMe12Ratio,             "cscSegClusterMe12Ratio[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterNStation",             cscSegClusterNStation,             "cscSegClusterNStation[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterMaxStation",             cscSegClusterMaxStation,             "cscSegClusterMaxStation[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterMaxStationRatio",             cscSegClusterMaxStationRatio,             "cscSegClusterMaxStationRatio[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterNChamber",             cscSegClusterNChamber,             "cscSegClusterNChamber[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterMaxChamber",             cscSegClusterMaxChamber,             "cscSegClusterMaxChamber[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterMaxChamberRatio",             cscSegClusterMaxChamberRatio,             "cscSegClusterMaxChamberRatio[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterVertexR",             cscSegClusterVertexR,             "cscSegClusterVertexR[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterVertexZ",             cscSegClusterVertexZ,             "cscSegClusterVertexZ[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterVertexDis",             cscSegClusterVertexDis,             "cscSegClusterVertexDis[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterVertexChi2",             cscSegClusterVertexChi2,             "cscSegClusterVertexChi2[nCscSegClusters]/F");
    displacedJetMuonTree->Branch("cscSegClusterVertexN1",             cscSegClusterVertexN1,             "cscSegClusterVertexN1[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterVertexN5",             cscSegClusterVertexN5,             "cscSegClusterVertexN5[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterVertexN10",             cscSegClusterVertexN10,             "cscSegClusterVertexN10[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterVertexN15",             cscSegClusterVertexN15,             "cscSegClusterVertexN15[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterVertexN20",             cscSegClusterVertexN20,             "cscSegClusterVertexN20[nCscSegClusters]/I");
    displacedJetMuonTree->Branch("cscSegClusterVertexN",             cscSegClusterVertexN,             "cscSegClusterVertexN[nCscSegClusters]/I");


    displacedJetMuonTree->Branch("nDtRechits",             &nDtRechits,             "nDtRechits/I");
    displacedJetMuonTree->Branch("dtRechitX",             dtRechitX,             "dtRechitX[nDtRechits]/F");
    displacedJetMuonTree->Branch("dtRechitY",             dtRechitY,             "dtRechitY[nDtRechits]/F");
    displacedJetMuonTree->Branch("dtRechitZ",             dtRechitZ,             "dtRechitZ[nDtRechits]/F");
    displacedJetMuonTree->Branch("dtRechitEta",             dtRechitEta,             "dtRechitEta[nDtRechits]/F");
    displacedJetMuonTree->Branch("dtRechitPhi",             dtRechitPhi,             "dtRechitPhi[nDtRechits]/F");
    displacedJetMuonTree->Branch("dtRechitTime",             dtRechitTime,             "dtRechitTime[nDtRechits]/F");
    displacedJetMuonTree->Branch("dtRechitStation",             dtRechitStation,             "dtRechitStation[nDtRechits]/I");
    displacedJetMuonTree->Branch("dtRechitWheel",             dtRechitWheel,             "dtRechitWheel[nDtRechits]/I");

    displacedJetMuonTree->Branch("nDtRechitClusters",             &nDtRechitClusters, "nDtRechitClusters/I");
    displacedJetMuonTree->Branch("dtRechitCluster_match_gParticle_minDeltaR",             dtRechitCluster_match_gParticle_minDeltaR,             "dtRechitCluster_match_gParticle_minDeltaR[nDtRechitClusters]/F");
    displacedJetMuonTree->Branch("dtRechitCluster_match_gParticle_index",             dtRechitCluster_match_gParticle_index,             "dtRechitCluster_match_gParticle_index[nDtRechitClusters]/I");
    displacedJetMuonTree->Branch("dtRechitCluster_match_gParticle_id",             dtRechitCluster_match_gParticle_id,             "dtRechitCluster_match_gParticle_id[nDtRechitClusters]/I");
    displacedJetMuonTree->Branch("dtRechitClusterX",             dtRechitClusterX,             "dtRechitClusterX[nDtRechitClusters]/F");
    displacedJetMuonTree->Branch("dtRechitClusterY",             dtRechitClusterY,             "dtRechitClusterY[nDtRechitClusters]/F");
    displacedJetMuonTree->Branch("dtRechitClusterZ",             dtRechitClusterZ,             "dtRechitClusterZ[nDtRechitClusters]/F");
    displacedJetMuonTree->Branch("dtRechitClusterTime",             dtRechitClusterTime,             "dtRechitClusterTime[nDtRechitClusters]/F");
    displacedJetMuonTree->Branch("dtRechitClusterTimeSpread",             dtRechitClusterTimeSpread,             "dtRechitClusterTimeSpread[nDtRechitClusters]/F");
    displacedJetMuonTree->Branch("dtRechitClusterGenMuonDeltaR",             dtRechitClusterGenMuonDeltaR,             "dtRechitClusterGenMuonDeltaR[nDtRechitClusters]/F");
    displacedJetMuonTree->Branch("dtRechitClusterMajorAxis",             dtRechitClusterMajorAxis,             "dtRechitClusterMajorAxis[nDtRechitClusters]/F");
    displacedJetMuonTree->Branch("dtRechitClusterMinorAxis",             dtRechitClusterMinorAxis,             "dtRechitClusterMinorAxis[nDtRechitClusters]/F");
    displacedJetMuonTree->Branch("dtRechitClusterEtaPhiSpread",             dtRechitClusterEtaPhiSpread,             "dtRechitClusterEtaPhiSpread[nDtRechitClusters]/F");
    displacedJetMuonTree->Branch("dtRechitClusterPhiSpread",             dtRechitClusterPhiSpread,             "dtRechitClusterPhiSpread[nDtRechitClusters]/F");
    displacedJetMuonTree->Branch("dtRechitClusterEtaSpread",             dtRechitClusterEtaSpread,             "dtRechitClusterEtaSpread[nDtRechitClusters]/F");
    displacedJetMuonTree->Branch("dtRechitClusterXSpread",             dtRechitClusterXSpread,             "dtRechitClusterXSpread[nDtRechitClusters]/F");
    displacedJetMuonTree->Branch("dtRechitClusterYSpread",             dtRechitClusterYSpread,             "dtRechitClusterYSpread[nDtRechitClusters]/F");
    displacedJetMuonTree->Branch("dtRechitClusterZSpread",             dtRechitClusterZSpread,             "dtRechitClusterZSpread[nDtRechitClusters]/F");
    displacedJetMuonTree->Branch("dtRechitClusterPhi",             dtRechitClusterPhi,             "dtRechitClusterPhi[nDtRechitClusters]/F");
    displacedJetMuonTree->Branch("dtRechitClusterEta",             dtRechitClusterEta,             "dtRechitClusterEta[nDtRechitClusters]/F");
    displacedJetMuonTree->Branch("dtRechitClusterJetVetoPt",       dtRechitClusterJetVetoPt,       "dtRechitClusterJetVetoPt[nDtRechitClusters]/F");
    displacedJetMuonTree->Branch("dtRechitClusterJetVetoE",        dtRechitClusterJetVetoE,        "dtRechitClusterJetVetoE[nDtRechitClusters]/F");
    displacedJetMuonTree->Branch("dtRechitClusterMuonVetoPt",      dtRechitClusterMuonVetoPt,      "dtRechitClusterMuonVetoPt[nDtRechitClusters]/F");
    displacedJetMuonTree->Branch("dtRechitClusterMuonVetoE",       dtRechitClusterMuonVetoE,       "dtRechitClusterMuonVetoE[nDtRechitClusters]/F");
    displacedJetMuonTree->Branch("dtRechitClusterCaloJetVeto",     dtRechitClusterCaloJetVeto,     "dtRechitClusterCaloJetVeto[nDtRechitClusters]/F");
    displacedJetMuonTree->Branch("dtRechitClusterSize",            dtRechitClusterSize,            "dtRechitClusterSize[nDtRechitClusters]/I");
    displacedJetMuonTree->Branch("dtRechitClusterNStation",             dtRechitClusterNStation,             "dtRechitClusterNStation[nDtRechitClusters]/I");
    displacedJetMuonTree->Branch("dtRechitClusterMaxStation",             dtRechitClusterMaxStation,             "dtRechitClusterMaxStation[nDtRechitClusters]/I");
    displacedJetMuonTree->Branch("dtRechitClusterMaxStationRatio",             dtRechitClusterMaxStationRatio,             "dtRechitClusterMaxStationRatio[nDtRechitClusters]/F");
    displacedJetMuonTree->Branch("dtRechitClusterNChamber",             dtRechitClusterNChamber,             "dtRechitClusterNChamber[nDtRechitClusters]/I");
    displacedJetMuonTree->Branch("dtRechitClusterMaxChamber",             dtRechitClusterMaxChamber,             "dtRechitClusterMaxChamber[nDtRechitClusters]/I");
    displacedJetMuonTree->Branch("dtRechitClusterMaxChamberRatio",             dtRechitClusterMaxChamberRatio,             "dtRechitClusterMaxChamberRatio[nDtRechitClusters]/F");
    displacedJetMuonTree->Branch("dtRechitClusterNSegmentStation1",            dtRechitClusterNSegmentStation1,            "dtRechitClusterNSegmentStation1[nDtRechitClusters]/I");
    displacedJetMuonTree->Branch("dtRechitClusterNSegmentStation2",            dtRechitClusterNSegmentStation2,            "dtRechitClusterNSegmentStation2[nDtRechitClusters]/I");
    displacedJetMuonTree->Branch("dtRechitClusterNSegmentStation3",            dtRechitClusterNSegmentStation3,            "dtRechitClusterNSegmentStation3[nDtRechitClusters]/I");
    displacedJetMuonTree->Branch("dtRechitClusterNSegmentStation4",            dtRechitClusterNSegmentStation4,            "dtRechitClusterNSegmentStation4[nDtRechitClusters]/I");


    displacedJetMuonTree->Branch("nRpc",&nRpc,"nRpc/I");
    displacedJetMuonTree->Branch("rpcPhi",rpcPhi,"rpcPhi[nRpc]/F");
    displacedJetMuonTree->Branch("rpcEta",rpcEta,"rpcEta[nRpc]/F");
    displacedJetMuonTree->Branch("rpcX",rpcX,"rpcX[nRpc]/F");
    displacedJetMuonTree->Branch("rpcY",rpcY,"rpcY[nRpc]/F");
    displacedJetMuonTree->Branch("rpcZ",rpcZ,"rpcZ[nRpc]/F");
    displacedJetMuonTree->Branch("rpcT",rpcT,"rpcT[nRpc]/F");
    displacedJetMuonTree->Branch("rpcBx",rpcBx,"rpcBx[nRpc]/I");
    displacedJetMuonTree->Branch("rpcTError",rpcTError,"rpcTError[nRpc]/F");


    displacedJetMuonTree->Branch("nDtSeg",&nDtSeg,"nDtSeg/I");
    displacedJetMuonTree->Branch("dtSegPhi",dtSegPhi,"dtSegPhi[nDtSeg]/F");
    displacedJetMuonTree->Branch("dtSegEta",dtSegEta,"dtSegEta[nDtSeg]/F");
    displacedJetMuonTree->Branch("dtSegX",dtSegX,"dtSegX[nDtSeg]/F");
    displacedJetMuonTree->Branch("dtSegY",dtSegY,"dtSegY[nDtSeg]/F");
    displacedJetMuonTree->Branch("dtSegZ",dtSegZ,"dtSegZ[nDtSeg]/F");
    displacedJetMuonTree->Branch("dtSegStation",dtSegStation,"dtSegStation[nDtSeg]/I");
    displacedJetMuonTree->Branch("dtSegWheel",dtSegWheel,"dtSegWheel[nDtSeg]/I");
    displacedJetMuonTree->Branch("dtSegTime",dtSegTime,"dtSegTime[nDtSeg]/F");
    displacedJetMuonTree->Branch("dtSegTimeError",dtSegTimeError,"dtSegTimeError[nDtSeg]/F");

    displacedJetMuonTree->Branch("nDtSegClusters",             &nDtSegClusters, "nDtSegClusters/I");
    displacedJetMuonTree->Branch("dtSegCluster_match_gParticle_minDeltaR",             dtSegCluster_match_gParticle_minDeltaR,             "dtSegCluster_match_gParticle_minDeltaR[nDtSegClusters]/F");
    displacedJetMuonTree->Branch("dtSegCluster_match_gParticle_index",             dtSegCluster_match_gParticle_index,             "dtSegCluster_match_gParticle_index[nDtSegClusters]/I");
    displacedJetMuonTree->Branch("dtSegCluster_match_gParticle_id",             dtSegCluster_match_gParticle_id,             "dtSegCluster_match_gParticle_id[nDtSegClusters]/I");
    displacedJetMuonTree->Branch("dtSegClusterX",             dtSegClusterX,             "dtSegClusterX[nDtSegClusters]/F");
    displacedJetMuonTree->Branch("dtSegClusterY",             dtSegClusterY,             "dtSegClusterY[nDtSegClusters]/F");
    displacedJetMuonTree->Branch("dtSegClusterZ",             dtSegClusterZ,             "dtSegClusterZ[nDtSegClusters]/F");
    displacedJetMuonTree->Branch("dtSegClusterTime",             dtSegClusterTime,             "dtSegClusterTime[nDtSegClusters]/F");
    displacedJetMuonTree->Branch("dtSegClusterTimeSpread",             dtSegClusterTimeSpread,             "dtSegClusterTimeSpread[nDtSegClusters]/F");
    displacedJetMuonTree->Branch("dtSegClusterGenMuonDeltaR",             dtSegClusterGenMuonDeltaR,             "dtSegClusterGenMuonDeltaR[nDtSegClusters]/F");
    displacedJetMuonTree->Branch("dtSegClusterMajorAxis",             dtSegClusterMajorAxis,             "dtSegClusterMajorAxis[nDtSegClusters]/F");
    displacedJetMuonTree->Branch("dtSegClusterMinorAxis",             dtSegClusterMinorAxis,             "dtSegClusterMinorAxis[nDtSegClusters]/F");
    displacedJetMuonTree->Branch("dtSegClusterEtaPhiSpread",             dtSegClusterEtaPhiSpread,             "dtSegClusterEtaPhiSpread[nDtSegClusters]/F");
    displacedJetMuonTree->Branch("dtSegClusterPhiSpread",             dtSegClusterPhiSpread,             "dtSegClusterPhiSpread[nDtSegClusters]/F");
    displacedJetMuonTree->Branch("dtSegClusterEtaSpread",             dtSegClusterEtaSpread,             "dtSegClusterEtaSpread[nDtSegClusters]/F");
    displacedJetMuonTree->Branch("dtSegClusterXSpread",             dtSegClusterXSpread,             "dtSegClusterXSpread[nDtSegClusters]/F");
    displacedJetMuonTree->Branch("dtSegClusterYSpread",             dtSegClusterYSpread,             "dtSegClusterYSpread[nDtSegClusters]/F");
    displacedJetMuonTree->Branch("dtSegClusterZSpread",             dtSegClusterZSpread,             "dtSegClusterZSpread[nDtSegClusters]/F");
    displacedJetMuonTree->Branch("dtSegClusterPhi",             dtSegClusterPhi,             "dtSegClusterPhi[nDtSegClusters]/F");
    displacedJetMuonTree->Branch("dtSegClusterEta",             dtSegClusterEta,             "dtSegClusterEta[nDtSegClusters]/F");
    displacedJetMuonTree->Branch("dtSegClusterJetVetoPt",       dtSegClusterJetVetoPt,       "dtSegClusterJetVetoPt[nDtSegClusters]/F");
    displacedJetMuonTree->Branch("dtSegClusterJetVetoE",        dtSegClusterJetVetoE,        "dtSegClusterJetVetoE[nDtSegClusters]/F");
    displacedJetMuonTree->Branch("dtSegClusterMuonVetoPt",      dtSegClusterMuonVetoPt,      "dtSegClusterMuonVetoPt[nDtSegClusters]/F");
    displacedJetMuonTree->Branch("dtSegClusterMuonVetoE",       dtSegClusterMuonVetoE,       "dtSegClusterMuonVetoE[nDtSegClusters]/F");
    displacedJetMuonTree->Branch("dtSegClusterCaloJetVeto",     dtSegClusterCaloJetVeto,     "dtSegClusterCaloJetVeto[nDtSegClusters]/F");
    displacedJetMuonTree->Branch("dtSegClusterSize",            dtSegClusterSize,            "dtSegClusterSize[nDtSegClusters]/I");
    displacedJetMuonTree->Branch("dtSegClusterNSegmentStation1",            dtSegClusterNSegmentStation1,            "dtSegClusterNSegmentStation1[nDtSegClusters]/I");
    displacedJetMuonTree->Branch("dtSegClusterNSegmentStation2",            dtSegClusterNSegmentStation2,            "dtSegClusterNSegmentStation2[nDtSegClusters]/I");
    displacedJetMuonTree->Branch("dtSegClusterNSegmentStation3",            dtSegClusterNSegmentStation3,            "dtSegClusterNSegmentStation3[nDtSegClusters]/I");
    displacedJetMuonTree->Branch("dtSegClusterNSegmentStation4",            dtSegClusterNSegmentStation4,            "dtSegClusterNSegmentStation4[nDtSegClusters]/I");

    // displacedJetMuonTree->Branch("dtSegClusterMe11Ratio",             dtSegClusterMe11Ratio,             "dtSegClusterMe11Ratio[nDtSegClusters]/F");
    // displacedJetMuonTree->Branch("dtSegClusterMe12Ratio",             dtSegClusterMe12Ratio,             "dtSegClusterMe12Ratio[nDtSegClusters]/F");
    displacedJetMuonTree->Branch("dtSegClusterNStation",             dtSegClusterNStation,             "dtSegClusterNStation[nDtSegClusters]/I");
    displacedJetMuonTree->Branch("dtSegClusterMaxStation",             dtSegClusterMaxStation,             "dtSegClusterMaxStation[nDtSegClusters]/I");
    displacedJetMuonTree->Branch("dtSegClusterMaxStationRatio",             dtSegClusterMaxStationRatio,             "dtSegClusterMaxStationRatio[nDtSegClusters]/F");
    displacedJetMuonTree->Branch("dtSegClusterNChamber",             dtSegClusterNChamber,             "dtSegClusterNChamber[nDtSegClusters]/I");
    displacedJetMuonTree->Branch("dtSegClusterMaxChamber",             dtSegClusterMaxChamber,             "dtSegClusterMaxChamber[nDtSegClusters]/I");
    displacedJetMuonTree->Branch("dtSegClusterMaxChamberRatio",             dtSegClusterMaxChamberRatio,             "dtSegClusterMaxChamberRatio[nDtSegClusters]/F");
    // displacedJetMuonTree->Branch("dtSegClusterVertexR",             dtSegClusterVertexR,             "dtSegClusterVertexR[nDtSegClusters]/F");
    // displacedJetMuonTree->Branch("dtSegClusterVertexZ",             dtSegClusterVertexZ,             "dtSegClusterVertexZ[nDtSegClusters]/F");
    // displacedJetMuonTree->Branch("dtSegClusterVertexDis",             dtSegClusterVertexDis,             "dtSegClusterVertexDis[nDtSegClusters]/F");
    // displacedJetMuonTree->Branch("dtSegClusterVertexChi2",             dtSegClusterVertexChi2,             "dtSegClusterVertexChi2[nDtSegClusters]/F");
    // displacedJetMuonTree->Branch("dtSegClusterVertexN1",             dtSegClusterVertexN1,             "dtSegClusterVertexN1[nDtSegClusters]/I");
    // displacedJetMuonTree->Branch("dtSegClusterVertexN5",             dtSegClusterVertexN5,             "dtSegClusterVertexN5[nDtSegClusters]/I");
    // displacedJetMuonTree->Branch("dtSegClusterVertexN10",             dtSegClusterVertexN10,             "dtSegClusterVertexN10[nDtSegClusters]/I");
    // displacedJetMuonTree->Branch("dtSegClusterVertexN15",             dtSegClusterVertexN15,             "dtSegClusterVertexN15[nDtSegClusters]/I");
    // displacedJetMuonTree->Branch("dtSegClusterVertexN20",             dtSegClusterVertexN20,             "dtSegClusterVertexN20[nDtSegClusters]/I");
    // displacedJetMuonTree->Branch("dtSegClusterVertexN",             dtSegClusterVertexN,             "dtSegClusterVertexN[nDtSegClusters]/I");

};
void displacedJetMuon_ntupler::enableHORechitBranches()
{
  displacedJetMuonTree->Branch("nHORechits", &nHORechits,"nHORechits/I");
  displacedJetMuonTree->Branch("hoRechit_Eta", hoRechit_Eta,"hoRechit_Eta[nHORechits]/F");
  displacedJetMuonTree->Branch("hoRechit_Phi",hoRechit_Phi, "hoRechit_Phi[nHORechits]/F");
  displacedJetMuonTree->Branch("hoRechit_E", hoRechit_E, "hoRechit_E[nHORechits]/F");
  displacedJetMuonTree->Branch("hoRechit_X", hoRechit_X, "hoRechit_X[nHORechits]/F");
  displacedJetMuonTree->Branch("hoRechit_Y", hoRechit_Y, "hoRechit_Y[nHORechits]/F");
  displacedJetMuonTree->Branch("hoRechit_Z", hoRechit_Z, "hoRechit_Z[nHORechits]/F");
  displacedJetMuonTree->Branch("hoRechit_T", hoRechit_T, "hoRechit_T[nHORechits]/F");

};

void displacedJetMuon_ntupler::enableEcalRechitBranches()
{
  displacedJetMuonTree->Branch("nRechits", &nRechits,"nRechits/I");
  displacedJetMuonTree->Branch("ecalRechit_Eta", ecalRechit_Eta,"ecalRechit_Eta[nRechits]/F");
  displacedJetMuonTree->Branch("ecalRechit_Phi",ecalRechit_Phi, "ecalRechit_Phi[nRechits]/F");
  displacedJetMuonTree->Branch("ecalRechit_E", ecalRechit_E, "ecalRechit_E[nRechits]/F");
  displacedJetMuonTree->Branch("ecalRechit_T", ecalRechit_T, "ecalRechit_T[nRechits]/F");
  displacedJetMuonTree->Branch("ecalRechit_E_Error", ecalRechit_E_Error, "ecalRechit_E_Error[nRechits]/F");
  displacedJetMuonTree->Branch("ecalRechit_T_Error", ecalRechit_T_Error, "ecalRechit_T_Error[nRechits]/F");
  displacedJetMuonTree->Branch("ecalRechit_kSaturatedflag", ecalRechit_kSaturatedflag, "ecalRechit_kSaturatedflag[nRechits]/O");
  displacedJetMuonTree->Branch("ecalRechit_kLeadingEdgeRecoveredflag", ecalRechit_kLeadingEdgeRecoveredflag, "ecalRechit_kLeadingEdgeRecoveredflag[nRechits]/O");
  displacedJetMuonTree->Branch("ecalRechit_kPoorRecoflag", ecalRechit_kPoorRecoflag, "ecalRechit_kPoorRecoflag[nRechits]/O");
  displacedJetMuonTree->Branch("ecalRechit_kWeirdflag", ecalRechit_kWeirdflag, "ecalRechit_kWeirdflag[nRechits]/O");
  displacedJetMuonTree->Branch("ecalRechit_kDiWeirdflag", ecalRechit_kDiWeirdflag, "ecalRechit_kDiWeirdflag[nRechits]/O");

  displacedJetMuonTree->Branch("nHBHERechits", &nHBHERechits,"nHBHERechits/I");
  displacedJetMuonTree->Branch("hbheRechit_Eta", hbheRechit_Eta,"hbheRechit_Eta[nHBHERechits]/F");
  displacedJetMuonTree->Branch("hbheRechit_Phi",hbheRechit_Phi, "hbheRechit_Phi[nHBHERechits]/F");
  displacedJetMuonTree->Branch("hbheRechit_E", hbheRechit_E, "hbheRechit_E[nHBHERechits]/F");
  displacedJetMuonTree->Branch("hbheRechit_X", hbheRechit_X, "hbheRechit_X[nHBHERechits]/F");
  displacedJetMuonTree->Branch("hbheRechit_Y", hbheRechit_Y, "hbheRechit_Y[nHBHERechits]/F");
  displacedJetMuonTree->Branch("hbheRechit_Z", hbheRechit_Z, "hbheRechit_Z[nHBHERechits]/F");
  displacedJetMuonTree->Branch("hbheRechit_T", hbheRechit_T, "hbheRechit_T[nHBHERechits]/F");
  displacedJetMuonTree->Branch("hbheRechit_iEta", hbheRechit_iEta,"hbheRechit_iEta[nHBHERechits]/I");
  displacedJetMuonTree->Branch("hbheRechit_iPhi", hbheRechit_iPhi,"hbheRechit_iPhi[nHBHERechits]/I");
  displacedJetMuonTree->Branch("hbheRechit_depth", hbheRechit_depth,"hbheRechit_depth[nHBHERechits]/I");

  displacedJetMuonTree->Branch("nTracks", &nTracks,"nTracks/I");
  displacedJetMuonTree->Branch("track_Pt", track_Pt,"track_Pt[nTracks]/F");
  displacedJetMuonTree->Branch("track_Eta", track_Eta,"track_Eta[nTracks]/F");
  displacedJetMuonTree->Branch("track_Phi", track_Phi,"track_Phi[nTracks]/F");
  displacedJetMuonTree->Branch("track_bestVertexIndex", track_bestVertexIndex,"track_bestVertexIndex[nTracks]/I");
  displacedJetMuonTree->Branch("track_nMissingInnerHits", track_nMissingInnerHits,"track_nMissingInnerHits[nTracks]/I");
  displacedJetMuonTree->Branch("track_nMissingOuterHits", track_nMissingOuterHits,"track_nMissingOuterHits[nTracks]/I");
  displacedJetMuonTree->Branch("track_dxyToBS", track_dxyToBS,"track_dxyToBS[nTracks]/F");
  displacedJetMuonTree->Branch("track_dxyErr", track_dxyErr,"track_dxyErr[nTracks]/F");
  displacedJetMuonTree->Branch("track_dzToPV", track_dzToPV,"track_dzToPV[nTracks]/F");
  displacedJetMuonTree->Branch("track_dzErr", track_dzErr,"track_dzErr[nTracks]/F");

};

void displacedJetMuon_ntupler::enableJetBranches()
{
  displacedJetMuonTree->Branch("nJets", &nJets,"nJets/I");
  displacedJetMuonTree->Branch("jetE", jetE,"jetE[nJets]/F");
  displacedJetMuonTree->Branch("jetPt", jetPt,"jetPt[nJets]/F");
  displacedJetMuonTree->Branch("jetEta", jetEta,"jetEta[nJets]/F");

  displacedJetMuonTree->Branch("jetPhi", jetPhi,"jetPhi[nJets]/F");
  displacedJetMuonTree->Branch("jetCSV", jetCSV,"jetCSV[nJets]/F");
  displacedJetMuonTree->Branch("jetCISV", jetCISV,"jetCISV[nJets]/F");
  displacedJetMuonTree->Branch("jetProbb", jetProbb,"jetProbb[nJets]/F");
  displacedJetMuonTree->Branch("jetProbc", jetProbc,"jetProbc[nJets]/F");
  displacedJetMuonTree->Branch("jetProbudsg", jetProbudsg,"jetProbudsg[nJets]/F");
  displacedJetMuonTree->Branch("jetProbbb", jetProbbb,"jetProbbb[nJets]/F");
  displacedJetMuonTree->Branch("jetMass", jetMass, "jetMass[nJets]/F");
  displacedJetMuonTree->Branch("jetJetArea", jetJetArea, "jetJetArea[nJets]/F");
  displacedJetMuonTree->Branch("jetPileupE", jetPileupE, "jetPileupE[nJets]/F");
  displacedJetMuonTree->Branch("jetPileupId", jetPileupId, "jetPileupId[nJets]/F");
  displacedJetMuonTree->Branch("jetPileupIdFlag", jetPileupIdFlag, "jetPileupIdFlag[nJets]/I");
  displacedJetMuonTree->Branch("jetPassIDLoose", jetPassIDLoose, "jetPassIDLoose[nJets]/O");
  displacedJetMuonTree->Branch("jetPassIDTight", jetPassIDTight, "jetPassIDTight[nJets]/O");
  displacedJetMuonTree->Branch("jetPassMuFrac", jetPassMuFrac, "jetPassMuFrac[nJets]/O");
  displacedJetMuonTree->Branch("jetPassEleFrac", jetPassEleFrac, "jetPassEleFrac[nJets]/O");
  displacedJetMuonTree->Branch("jetPartonFlavor", jetPartonFlavor, "jetPartonFlavor[nJets]/I");
  displacedJetMuonTree->Branch("jetHadronFlavor", jetHadronFlavor, "jetHadronFlavor[nJets]/I");
  displacedJetMuonTree->Branch("jetChargedEMEnergyFraction", jetChargedEMEnergyFraction, "jetChargedEMEnergyFraction[nJets]/F");
  displacedJetMuonTree->Branch("jetNeutralEMEnergyFraction", jetNeutralEMEnergyFraction, "jetNeutralEMEnergyFraction[nJets]/F");
  displacedJetMuonTree->Branch("jetChargedHadronEnergyFraction", jetChargedHadronEnergyFraction, "jetChargedHadronEnergyFraction[nJets]/F");
  displacedJetMuonTree->Branch("jetNeutralHadronEnergyFraction", jetNeutralHadronEnergyFraction, "jetNeutralHadronEnergyFraction[nJets]/F");
  displacedJetMuonTree->Branch("jetMuonEnergyFraction", jetMuonEnergyFraction, "jetMuonEnergyFraction[nJets]/F");
  displacedJetMuonTree->Branch("jetHOEnergyFraction", jetHOEnergyFraction, "jetHOEnergyFraction[nJets]/F");
  displacedJetMuonTree->Branch("jetHFHadronEnergyFraction", jetHFHadronEnergyFraction, "jetHFHadronEnergyFraction[nJets]/F");
  displacedJetMuonTree->Branch("jetHFEMEnergyFraction",jetHFEMEnergyFraction, "jetHFEMEnergyFraction[nJets]/F");
  displacedJetMuonTree->Branch("jetAllMuonPt", jetAllMuonPt,"jetAllMuonPt[nJets]/F");
  displacedJetMuonTree->Branch("jetAllMuonEta", jetAllMuonEta,"jetAllMuonEta[nJets]/F");
  displacedJetMuonTree->Branch("jetAllMuonPhi", jetAllMuonPhi,"jetAllMuonPhi[nJets]/F");
  displacedJetMuonTree->Branch("jetAllMuonM", jetAllMuonM,"jetAllMuonM[nJets]/F");
  displacedJetMuonTree->Branch("jetPtWeightedDZ", jetPtWeightedDZ,"jetPtWeightedDZ[nJets]/F");
  displacedJetMuonTree->Branch("jetNRechits", jetNRechits,"jetNRechits[nJets]/I");
  displacedJetMuonTree->Branch("jetRechitE", jetRechitE,"jetRechitE[nJets]/F");
  displacedJetMuonTree->Branch("jetRechitT", jetRechitT,"jetRechitT[nJets]/F");
  displacedJetMuonTree->Branch("jetRechitT_rms", jetRechitT_rms,"jetRechitT_rms[nJets]/F");

  displacedJetMuonTree->Branch("jetRechitE_Error", jetRechitE_Error,"jetRechitE_Error[nJets]/F");
  displacedJetMuonTree->Branch("jetRechitT_Error", jetRechitT_Error,"jetRechitT_Error[nJets]/F");
  displacedJetMuonTree->Branch("jetAlphaMax",jetAlphaMax,"jetAlphaMax[nJets]/F");
  displacedJetMuonTree->Branch("jetBetaMax",jetBetaMax,"jetBetaMax[nJets]/F");
  displacedJetMuonTree->Branch("jetGammaMax_ET",jetGammaMax_ET,"jetGammaMax_ET[nJets]/F");
  displacedJetMuonTree->Branch("jetGammaMax_EM",jetGammaMax_EM,"jetGammaMax_EM[nJets]/F");
  displacedJetMuonTree->Branch("jetGammaMax_Hadronic",jetGammaMax_Hadronic,"jetGammaMax_Hadronic[nJets]/F");
  displacedJetMuonTree->Branch("jetGammaMax",jetGammaMax,"jetGammaMax[nJets]/F");
  displacedJetMuonTree->Branch("jetPtAllTracks",jetPtAllTracks,"jetPtAllTracks[nJets]/F");
  displacedJetMuonTree->Branch("jetPtAllPVTracks",jetPtAllPVTracks,"jetPtAllPVTracks[nJets]/F");
  displacedJetMuonTree->Branch("jetMedianTheta2D",jetMedianTheta2D,"jetMedianTheta2D[nJets]/F");
  displacedJetMuonTree->Branch("jetMedianIP",jetMedianIP,"jetMedianIP[nJets]/F");
  displacedJetMuonTree->Branch("jetMinDeltaRAllTracks",jetMinDeltaRAllTracks,"jetMinDeltaRAllTracks[nJets]/F");
  displacedJetMuonTree->Branch("jetMinDeltaRPVTracks",jetMinDeltaRPVTracks,"jetMinDeltaRPVTracks[nJets]/F");
  displacedJetMuonTree->Branch("jet_sig_et1",jet_sig_et1,"jet_sig_et1[nJets]/F");
  displacedJetMuonTree->Branch("jet_sig_et2",jet_sig_et2,"jet_sig_et2[nJets]/F");
  displacedJetMuonTree->Branch("jet_energy_frac",jet_energy_frac,"jet_energy_frac[nJets]/F");

  displacedJetMuonTree->Branch("jetAlphaMax_wp",jetAlphaMax_wp,"jetAlphaMax_wp[nJets]/F");
  displacedJetMuonTree->Branch("jetBetaMax_wp",jetBetaMax_wp,"jetBetaMax_wp[nJets]/F");
  displacedJetMuonTree->Branch("jetGammaMax_ET_wp",jetGammaMax_ET_wp,"jetGammaMax_ET_wp[nJets]/F");
  displacedJetMuonTree->Branch("jetGammaMax_EM_wp",jetGammaMax_EM_wp,"jetGammaMax_EM_wp[nJets]/F");
  displacedJetMuonTree->Branch("jetGammaMax_Hadronic_wp",jetGammaMax_Hadronic_wp,"jetGammaMax_Hadronic_wp[nJets]/F");
  displacedJetMuonTree->Branch("jetGammaMax_wp",jetGammaMax_wp,"jetGammaMax_wp[nJets]/F");
  displacedJetMuonTree->Branch("jetPtAllTracks_wp",jetPtAllTracks_wp,"jetPtAllTracks_wp[nJets]/F");
  displacedJetMuonTree->Branch("jetPtAllPVTracks_wp",jetPtAllPVTracks_wp,"jetPtAllPVTracks_wp[nJets]/F");
  displacedJetMuonTree->Branch("jetMedianTheta2D_wp",jetMedianTheta2D_wp,"jetMedianTheta2D_wp[nJets]/F");
  displacedJetMuonTree->Branch("jetMedianIP_wp",jetMedianIP_wp,"jetMedianIP_wp[nJets]/F");
  displacedJetMuonTree->Branch("jetMinDeltaRAllTracks_wp",jetMinDeltaRAllTracks_wp,"jetMinDeltaRAllTracks_wp[nJets]/F");
  displacedJetMuonTree->Branch("jetMinDeltaRPVTracks_wp",jetMinDeltaRPVTracks_wp,"jetMinDeltaRPVTracks_wp[nJets]/F");


};



void displacedJetMuon_ntupler::enableJetAK8Branches()
{
  displacedJetMuonTree->Branch("nFatJets", &nFatJets,"nFatJets/i");
  displacedJetMuonTree->Branch("fatJetE", fatJetE,"fatJetE[nFatJets]/F");
  displacedJetMuonTree->Branch("fatJetPt", fatJetPt,"fatJetPt[nFatJets]/F");
  displacedJetMuonTree->Branch("fatJetEta", fatJetEta,"fatJetEta[nFatJets]/F");
  displacedJetMuonTree->Branch("fatJetPhi", fatJetPhi,"fatJetPhi[nFatJets]/F");
  displacedJetMuonTree->Branch("fatJetCorrectedPt", fatJetCorrectedPt,"fatJetCorrectedPt[nFatJets]/F");
  // displacedJetMuonTree->Branch("fatJetCorrectedEta", fatJetCorrectedEta,"fatJetCorrectedEta[nFatJets]/F");
  // displacedJetMuonTree->Branch("fatJetCorrectedPhi", fatJetCorrectedPhi,"fatJetCorrectedPhi[nFatJets]/F");
  displacedJetMuonTree->Branch("fatJetPrunedM", fatJetPrunedM,"fatJetPrunedM[nFatJets]/F");
  displacedJetMuonTree->Branch("fatJetTrimmedM", fatJetTrimmedM,"fatJetTrimmedM[nFatJets]/F");
  displacedJetMuonTree->Branch("fatJetFilteredM", fatJetFilteredM,"fatJetFilteredM[nFatJets]/F");
  displacedJetMuonTree->Branch("fatJetSoftDropM", fatJetSoftDropM,"fatJetSoftDropM[nFatJets]/F");
  displacedJetMuonTree->Branch("fatJetCorrectedSoftDropM", fatJetCorrectedSoftDropM,"fatJetCorrectedSoftDropM[nFatJets]/F");
  displacedJetMuonTree->Branch("fatJetUncorrectedSoftDropM", fatJetUncorrectedSoftDropM,"fatJetUncorrectedSoftDropM[nFatJets]/F");
  displacedJetMuonTree->Branch("fatJetTau1", fatJetTau1,"fatJetTau1[nFatJets]/F");
  displacedJetMuonTree->Branch("fatJetTau2", fatJetTau2,"fatJetTau2[nFatJets]/F");
  displacedJetMuonTree->Branch("fatJetTau3", fatJetTau3,"fatJetTau3[nFatJets]/F");
  displacedJetMuonTree->Branch("fatJetMaxSubjetCSV", fatJetMaxSubjetCSV, "fatJetMaxSubjetCSV[nFatJets]/F");
  displacedJetMuonTree->Branch("fatJetPassIDLoose", fatJetPassIDLoose,"fatJetPassIDLoose[nFatJets]/O");
  displacedJetMuonTree->Branch("fatJetPassIDTight", fatJetPassIDTight,"fatJetPassIDTight[nFatJets]/O");
};

void displacedJetMuon_ntupler::enableMetBranches()
{
  displacedJetMuonTree->Branch("metPt", &metPt, "metPt/F");
  displacedJetMuonTree->Branch("metPhi", &metPhi, "metPhi/F");
  displacedJetMuonTree->Branch("sumMET", &sumMET, "sumMET/F");
  displacedJetMuonTree->Branch("metType0Pt", &metType0Pt, "metType0Pt/F");
  displacedJetMuonTree->Branch("metType0Phi", &metType0Phi, "metType0Phi/F");
  displacedJetMuonTree->Branch("metType1Pt_raw", &metType1Pt_raw, "metType1Pt_raw/F");
  displacedJetMuonTree->Branch("metType1Pt", &metType1Pt, "metType1Pt/F");
  displacedJetMuonTree->Branch("metType1Px", &metType1Px, "metType1Px/F");
  displacedJetMuonTree->Branch("metType1Py", &metType1Py, "metType1Py/F");
  displacedJetMuonTree->Branch("metType1Eta", &metType1Eta, "metType1Eta/F");
  displacedJetMuonTree->Branch("metType1Phi", &metType1Phi, "metType1Phi/F");
  displacedJetMuonTree->Branch("metType1Phi_raw", &metType1Phi_raw, "metType1Phi_raw/F");
  displacedJetMuonTree->Branch("metType0Plus1Pt", &metType0Plus1Pt, "metType0Plus1Pt/F");
  displacedJetMuonTree->Branch("metType0Plus1Phi", &metType0Plus1Phi, "metType0Plus1Phi/F");
  displacedJetMuonTree->Branch("metNoHFPt", &metNoHFPt, "metNoHFPt/F");
  displacedJetMuonTree->Branch("metNoHFPhi", &metNoHFPhi, "metNoHFPhi/F");
  displacedJetMuonTree->Branch("metPuppiPt", &metPuppiPt, "metPuppiPt/F");
  displacedJetMuonTree->Branch("metPuppiPhi", &metPuppiPhi, "metPuppiPhi/F");
  displacedJetMuonTree->Branch("metCaloPt", &metCaloPt, "metCaloPt/F");
  displacedJetMuonTree->Branch("metCaloPhi", &metCaloPhi, "metCaloPhi/F");

  displacedJetMuonTree->Branch("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, "Flag_HBHENoiseFilter/O");
  displacedJetMuonTree->Branch("Flag_HBHETightNoiseFilter", &Flag_HBHETightNoiseFilter, "Flag_HBHETightNoiseFilter/O");
  displacedJetMuonTree->Branch("Flag_HBHEIsoNoiseFilter", &Flag_HBHEIsoNoiseFilter, "Flag_HBHEIsoNoiseFilter/O");
  displacedJetMuonTree->Branch("Flag_badChargedCandidateFilter", &Flag_badChargedCandidateFilter, "Flag_badChargedCandidateFilter/O");
  displacedJetMuonTree->Branch("Flag_badMuonFilter", &Flag_badMuonFilter, "Flag_badMuonFilter/O");
  displacedJetMuonTree->Branch("Flag_badGlobalMuonFilter", &Flag_badGlobalMuonFilter, "Flag_badGlobalMuonFilter/O");
  displacedJetMuonTree->Branch("Flag_duplicateMuonFilter", &Flag_duplicateMuonFilter, "Flag_duplicateMuonFilter/O");
  displacedJetMuonTree->Branch("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter, "Flag_globalSuperTightHalo2016Filter/O");
  displacedJetMuonTree->Branch("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter, "Flag_globalTightHalo2016Filter/O");

  displacedJetMuonTree->Branch("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, "Flag_hcalLaserEventFilter/O");
  displacedJetMuonTree->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, "Flag_EcalDeadCellTriggerPrimitiveFilter/O");
  displacedJetMuonTree->Branch("Flag_EcalDeadCellBoundaryEnergyFilter", &Flag_EcalDeadCellBoundaryEnergyFilter, "Flag_EcalDeadCellBoundaryEnergyFilter/O");
  displacedJetMuonTree->Branch("Flag_goodVertices", &Flag_goodVertices, "Flag_goodVertices/O");
  displacedJetMuonTree->Branch("Flag_trackingFailureFilter", &Flag_trackingFailureFilter, "Flag_trackingFailureFilter/O");
  displacedJetMuonTree->Branch("Flag_eeBadScFilter", &Flag_eeBadScFilter, "Flag_eeBadScFilter/O");
  displacedJetMuonTree->Branch("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, "Flag_ecalLaserCorrFilter/O");
  displacedJetMuonTree->Branch("Flag_trkPOGFilters", &Flag_trkPOGFilters, "Flag_trkPOGFilters/O");
  displacedJetMuonTree->Branch("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, "Flag_trkPOG_manystripclus53X/O");
  displacedJetMuonTree->Branch("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, "Flag_trkPOG_toomanystripclus53X/O");
  displacedJetMuonTree->Branch("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, "Flag_trkPOG_logErrorTooManyClusters/O");
  displacedJetMuonTree->Branch("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, "Flag_BadPFMuonFilter/O");
  displacedJetMuonTree->Branch("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, "Flag_BadChargedCandidateFilter/O");
  displacedJetMuonTree->Branch("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, "Flag_ecalBadCalibFilter/O");
  displacedJetMuonTree->Branch("Flag_METFilters", &Flag_METFilters, "Flag_METFilters/O");

  displacedJetMuonTree->Branch("Flag2_globalSuperTightHalo2016Filter", &Flag2_globalSuperTightHalo2016Filter, "Flag2_globalSuperTightHalo2016Filter/O");
  displacedJetMuonTree->Branch("Flag2_globalTightHalo2016Filter", &Flag2_globalTightHalo2016Filter, "Flag2_globalTightHalo2016Filter/O");
  displacedJetMuonTree->Branch("Flag2_goodVertices", &Flag2_goodVertices, "Flag2_goodVertices/O");
  displacedJetMuonTree->Branch("Flag2_BadChargedCandidateFilter", &Flag2_BadChargedCandidateFilter, "Flag2_BadChargedCandidateFilter/O");
  displacedJetMuonTree->Branch("Flag2_BadPFMuonFilter", &Flag2_BadPFMuonFilter, "Flag2_BadPFMuonFilter/O");
  displacedJetMuonTree->Branch("Flag2_EcalDeadCellTriggerPrimitiveFilter", &Flag2_EcalDeadCellTriggerPrimitiveFilter, "Flag2_EcalDeadCellTriggerPrimitiveFilter/O");
  displacedJetMuonTree->Branch("Flag2_HBHENoiseFilter", &Flag2_HBHENoiseFilter, "Flag2_HBHENoiseFilter/O");
  displacedJetMuonTree->Branch("Flag2_HBHEIsoNoiseFilter", &Flag2_HBHEIsoNoiseFilter, "Flag2_HBHEIsoNoiseFilter/O");
  displacedJetMuonTree->Branch("Flag2_ecalBadCalibFilter", &Flag2_ecalBadCalibFilter, "Flag2_ecalBadCalibFilter/O");
  displacedJetMuonTree->Branch("Flag2_eeBadScFilter", &Flag2_eeBadScFilter, "Flag2_eeBadScFilter/O");

};

void displacedJetMuon_ntupler::enableTriggerBranches()
{
  nameHLT = new std::vector<std::string>; nameHLT->clear();
  displacedJetMuonTree->Branch("HLTDecision", &triggerDecision, ("HLTDecision[" + std::to_string(NTriggersMAX) +  "]/O").c_str());
  displacedJetMuonTree->Branch("HLTPrescale", &triggerHLTPrescale, ("HLTPrescale[" + std::to_string(NTriggersMAX) +  "]/I").c_str()); 
};

void displacedJetMuon_ntupler::enableMCBranches()
{
  displacedJetMuonTree->Branch("nGenJets", &nGenJets, "nGenJets/I");
  displacedJetMuonTree->Branch("genJetE", genJetE, "genJetE[nGenJets]/F");
  displacedJetMuonTree->Branch("genJetPt", genJetPt, "genJetPt[nGenJets]/F");
  displacedJetMuonTree->Branch("genJetEta", genJetEta, "genJetEta[nGenJets]/F");
  displacedJetMuonTree->Branch("genJetPhi", genJetPhi, "genJetPhi[nGenJets]/F");
  displacedJetMuonTree->Branch("genMetPtCalo", &genMetPtCalo, "genMetPtCalo/F");
  displacedJetMuonTree->Branch("genMetPhiCalo", &genMetPhiCalo, "genMetPhiCalo/F");
  displacedJetMuonTree->Branch("genMetPtTrue", &genMetPtTrue, "genMetPtTrue/F");
  displacedJetMuonTree->Branch("genMetPhiTrue", &genMetPhiTrue, "genMetPhiTrue/F");
  displacedJetMuonTree->Branch("genVertexX", &genVertexX, "genVertexX/F");
  displacedJetMuonTree->Branch("genVertexY", &genVertexY, "genVertexY/F");
  displacedJetMuonTree->Branch("genVertexZ", &genVertexZ, "genVertexZ/F");
  displacedJetMuonTree->Branch("genVertexT", &genVertexT, "genVertexT/F");
  displacedJetMuonTree->Branch("genWeight", &genWeight, "genWeight/F");
  displacedJetMuonTree->Branch("genSignalProcessID", &genSignalProcessID, "genSignalProcessID/i");
  displacedJetMuonTree->Branch("genQScale", &genQScale, "genQScale/F");
  displacedJetMuonTree->Branch("genAlphaQCD", &genAlphaQCD, "genAlphaQCD/F");
  displacedJetMuonTree->Branch("genAlphaQED", &genAlphaQED, "genAlphaQED/F");
  scaleWeights = new std::vector<float>; scaleWeights->clear();
  pdfWeights = new std::vector<float>; pdfWeights->clear();
  alphasWeights = new std::vector<float>; alphasWeights->clear();
  displacedJetMuonTree->Branch("lheComments", "std::string",&lheComments);
  displacedJetMuonTree->Branch("scaleWeights", "std::vector<float>",&scaleWeights);
  displacedJetMuonTree->Branch("pdfWeights", "std::vector<float>",&pdfWeights);
  displacedJetMuonTree->Branch("alphasWeights", "std::vector<float>",&alphasWeights);
};

void displacedJetMuon_ntupler::enableGenParticleBranches()
{
  displacedJetMuonTree->Branch("nGenParticle", &nGenParticle, "nGenParticle/I");
  displacedJetMuonTree->Branch("gParticleMotherId", gParticleMotherId, "gParticleMotherId[nGenParticle]/I");
  displacedJetMuonTree->Branch("gParticleMotherIndex", gParticleMotherIndex, "gParticleMotherIndex[nGenParticle]/I");
  displacedJetMuonTree->Branch("gParticleId", gParticleId, "gParticleId[nGenParticle]/I");
  displacedJetMuonTree->Branch("gParticleStatus", gParticleStatus, "gParticleStatus[nGenParticle]/I");
  displacedJetMuonTree->Branch("gParticleE", gParticleE, "gParticleE[nGenParticle]/F");
  displacedJetMuonTree->Branch("gParticlePt", gParticlePt, "gParticlePt[nGenParticle]/F");
  displacedJetMuonTree->Branch("gParticleEta", gParticleEta, "gParticleEta[nGenParticle]/F");
  displacedJetMuonTree->Branch("gParticlePhi", gParticlePhi, "gParticlePhi[nGenParticle]/F");
  displacedJetMuonTree->Branch("gParticleProdVertexX", gParticleProdVertexX, "gParticleProdVertexX[nGenParticle]/F");
  displacedJetMuonTree->Branch("gParticleProdVertexY", gParticleProdVertexY, "gParticleProdVertexY[nGenParticle]/F");
  displacedJetMuonTree->Branch("gParticleProdVertexZ", gParticleProdVertexZ, "gParticleProdVertexZ[nGenParticle]/F");
  displacedJetMuonTree->Branch("gParticleDecayVertexX", gParticleDecayVertexX, "gParticleDecayVertexX[nGenParticle]/F");
  displacedJetMuonTree->Branch("gParticleDecayVertexY", gParticleDecayVertexY, "gParticleDecayVertexY[nGenParticle]/F");
  displacedJetMuonTree->Branch("gParticleDecayVertexZ", gParticleDecayVertexZ, "gParticleDecayVertexZ[nGenParticle]/F");

   // displacedJetMuonTree->Branch("gLLP_prod_vertex_x", gLLP_prod_vertex_x, "gLLP_prod_vertex_x[2]/F");
   // displacedJetMuonTree->Branch("gLLP_prod_vertex_y", gLLP_prod_vertex_y, "gLLP_prod_vertex_y[2]/F");
   // displacedJetMuonTree->Branch("gLLP_prod_vertex_z", gLLP_prod_vertex_z, "gLLP_prod_vertex_z[2]/F");
   displacedJetMuonTree->Branch("gLLP_decay_vertex_x", gLLP_decay_vertex_x, "gLLP_decay_vertex_x[2]/F");
   displacedJetMuonTree->Branch("gLLP_decay_vertex_y", gLLP_decay_vertex_y, "gLLP_decay_vertex_y[2]/F");
   displacedJetMuonTree->Branch("gLLP_decay_vertex_z", gLLP_decay_vertex_z, "gLLP_decay_vertex_z[2]/F");
   displacedJetMuonTree->Branch("gLLP_beta", gLLP_beta, "gLLP_beta[2]/F");
   displacedJetMuonTree->Branch("gLLP_e", gLLP_e, "gLLP_e[2]/F");
   displacedJetMuonTree->Branch("gLLP_pt", gLLP_pt, "gLLP_pt[2]/F");
   displacedJetMuonTree->Branch("gLLP_eta", gLLP_eta, "gLLP_eta[2]/F");
   displacedJetMuonTree->Branch("gLLP_phi", gLLP_phi, "gLLP_phi[2]/F");
   displacedJetMuonTree->Branch("gLLP_csc", gLLP_csc, "gLLP_csc[2]/O");
   displacedJetMuonTree->Branch("gLLP_dt", gLLP_dt, "gLLP_dt[2]/O");
   displacedJetMuonTree->Branch("gLLP_travel_time", gLLP_travel_time, "gLLP_travel_time[2]/F");

   displacedJetMuonTree->Branch("gLLP_daughter_id", gLLP_daughter_id, "gLLP_daughter_id[4]/I");
   displacedJetMuonTree->Branch("gLLP_daughter_pt", gLLP_daughter_pt, "gLLP_daughter_pt[4]/F");
   displacedJetMuonTree->Branch("gLLP_daughter_eta", gLLP_daughter_eta, "gLLP_daughter_eta[4]/F");
   displacedJetMuonTree->Branch("gLLP_daughter_phi", gLLP_daughter_phi, "gLLP_daughter_phi[4]/F");
   displacedJetMuonTree->Branch("gLLP_daughter_eta_ecalcorr", gLLP_daughter_eta_ecalcorr, "gLLP_daughter_eta_ecalcorr[4]/F");
   displacedJetMuonTree->Branch("gLLP_daughter_phi_ecalcorr", gLLP_daughter_phi_ecalcorr, "gLLP_daughter_phi_ecalcorr[4]/F");
   displacedJetMuonTree->Branch("gLLP_daughter_e", gLLP_daughter_e, "gLLP_daughter_e[4]/F");
   displacedJetMuonTree->Branch("gLLP_daughter_mass", gLLP_daughter_mass, "gLLP_daughter_mass[4]/F");


   displacedJetMuonTree->Branch("gLLP_grandaughter_id", gLLP_grandaughter_id, "gLLP_grandaughter_id[4]/I");
   displacedJetMuonTree->Branch("gLLP_grandaughter_pt", gLLP_grandaughter_pt, "gLLP_grandaughter_pt[4]/F");
   displacedJetMuonTree->Branch("gLLP_grandaughter_eta", gLLP_grandaughter_eta, "gLLP_grandaughter_eta[4]/F");
   displacedJetMuonTree->Branch("gLLP_grandaughter_phi", gLLP_grandaughter_phi, "gLLP_grandaughter_phi[4]/F");
   displacedJetMuonTree->Branch("gLLP_grandaughter_eta_ecalcorr", gLLP_grandaughter_eta_ecalcorr, "gLLP_grandaughter_eta_ecalcorr[4]/F");
   displacedJetMuonTree->Branch("gLLP_grandaughter_phi_ecalcorr", gLLP_grandaughter_phi_ecalcorr, "gLLP_grandaughter_phi_ecalcorr[4]/F");
   displacedJetMuonTree->Branch("gLLP_grandaughter_e", gLLP_grandaughter_e, "gLLP_grandaughter_e[4]/F");
   displacedJetMuonTree->Branch("gLLP_grandaughter_mass", gLLP_grandaughter_mass, "gLLP_grandaughter_mass[4]/F");






};

//------ Load the miniAOD objects and reset tree variables for each event ------//
void displacedJetMuon_ntupler::loadEvent(const edm::Event& iEvent)//load all miniAOD objects for the current event
{
  iEvent.getByToken(hcalRecHitsHOToken_,hcalRecHitsHO);
  iEvent.getByToken(hcalRecHitsHBHEToken_,hcalRecHitsHBHE);
  iEvent.getByToken(triggerBitsToken_, triggerBits);
  iEvent.getByToken(hepMCToken_, hepMC);
  iEvent.getByToken(metFilterBitsToken_, metFilterBits);


  iEvent.getByToken(verticesToken_, vertices);
  iEvent.getByToken(cscSegmentInputToken_,cscSegments);
  iEvent.getByToken(cscRechitInputToken_,cscRechits);
  iEvent.getByToken(dtSegmentInputToken_,dtSegments);
  iEvent.getByToken(dtCosmicSegmentInputToken_,dtCosmicSegments);
  iEvent.getByToken(dtRechitInputToken_,dtRechits);
  iEvent.getByToken(dtCosmicRechitInputToken_,dtCosmicRechits);

  iEvent.getByToken(rpcRecHitInputToken_,rpcRecHits);
  if (!isData_) {
    iEvent.getByToken(MuonCSCSimHitsToken_, MuonCSCSimHits);
    if (readMuonDigis_) {
      iEvent.getByToken(MuonCSCStripDigiSimLinksToken_, MuonCSCStripDigiSimLinks);
      iEvent.getByToken(MuonCSCWireDigiSimLinksToken_, MuonCSCWireDigiSimLinks);
    }
  }
  if (readMuonDigis_) {
    iEvent.getByToken(MuonCSCStripDigiToken_, MuonCSCStripDigi);
    iEvent.getByToken(MuonCSCWireDigiToken_, MuonCSCWireDigi);
  }

  iEvent.getByToken(tracksTag_,tracks);
  iEvent.getByToken(PFCandsToken_, pfCands);
  iEvent.getByToken(PFClustersToken_, pfClusters);
  iEvent.getByToken(muonsToken_, muons);
  iEvent.getByToken(electronsToken_, electrons);
  iEvent.getByToken(photonsToken_, photons);
  iEvent.getByToken(tausToken_, taus);
  iEvent.getByToken(jetsToken_, jets);
  iEvent.getByToken(jetsPuppiToken_, jetsPuppi);
  iEvent.getByToken(jetsAK8Token_, jetsAK8);
  iEvent.getByToken(genMetCaloToken_, genMetsCalo);
  iEvent.getByToken(genMetTrueToken_, genMetsTrue);
  iEvent.getByToken(metToken_, mets);
  iEvent.getByToken(metPuppiToken_, metsPuppi);
  iEvent.getByToken(secondaryVerticesToken_,secondaryVertices);
  iEvent.getByToken(rhoAllToken_,rhoAll);
  iEvent.getByToken(rhoFastjetAllToken_,rhoFastjetAll);
  iEvent.getByToken(rhoFastjetAllCaloToken_,rhoFastjetAllCalo);
  iEvent.getByToken(rhoFastjetCentralCaloToken_,rhoFastjetCentralCalo);
  iEvent.getByToken(rhoFastjetCentralChargedPileUpToken_,rhoFastjetCentralChargedPileUp);
  iEvent.getByToken(rhoFastjetCentralNeutralToken_,rhoFastjetCentralNeutral);
  iEvent.getByToken(beamSpotToken_,beamSpot);
  iEvent.getByToken(ebRecHitsToken_,ebRecHits);
  iEvent.getByToken(eeRecHitsToken_,eeRecHits);
  iEvent.getByToken(esRecHitsToken_,esRecHits);
  iEvent.getByToken(ebeeClustersToken_,ebeeClusters);
  iEvent.getByToken(esClustersToken_,esClusters);
  iEvent.getByToken(conversionsToken_,conversions);
  iEvent.getByToken(singleLegConversionsToken_,singleLegConversions);
  iEvent.getByToken(gedGsfElectronCoresToken_,gedGsfElectronCores);
  iEvent.getByToken(gedPhotonCoresToken_, gedPhotonCores);
  iEvent.getByToken(generalTrackToken_,generalTracks);
  iEvent.getByToken(generalTrackHandleToken_,generalTrackHandle);
  if(readGenVertexTime_) iEvent.getByToken(genParticles_t0_Token_,genParticles_t0);
  if (useGen_) {
    iEvent.getByToken(genParticlesToken_,genParticles);
    iEvent.getByToken(genJetsToken_,genJets);
    iEvent.getByToken(genInfoToken_,genInfo);
    iEvent.getByToken(puInfoToken_,puInfo);
  }


}

//called by the loadEvent() method
void displacedJetMuon_ntupler::resetBranches()
{
    //reset tree variables
    resetEventInfoBranches();
    resetPVAllBranches();
    resetPileUpBranches();
    resetMuonBranches();
    resetElectronBranches();
    resetTauBranches();
    resetPhotonBranches();
    resetJetBranches();
    resetMuonSystemBranches();
    resetMetBranches();
    resetGenParticleBranches();
    resetMCBranches();
    resetTriggerBranches();
    resetEcalRechitBranches();
    resetHORechitBranches();

};

void displacedJetMuon_ntupler::resetEventInfoBranches()
{
  //Event
  nPV = -1;
  eventNum = 0;
  eventTime = 0;
  lumiNum = 0;
  runNum = 0;
  nSlimmedSecondV = 0;
  pvX = -999.0;
  pvY = -999.0;
  pvZ = -999.0;
  fixedGridRhoAll = -999.0;
  fixedGridRhoFastjetAll = -999.0;
  fixedGridRhoFastjetAllCalo = -999.0;
  fixedGridRhoFastjetCentralCalo = -999.0;
  fixedGridRhoFastjetCentralChargedPileUp = -999.0;
  fixedGridRhoFastjetCentralNeutral = -999.0;
  return;
};

void displacedJetMuon_ntupler::resetPVAllBranches()
{
  nPVAll = 0;
  for(int i = 0; i < MAX_NPV; i++)
  {
    pvAllX[i] = -999.;
    pvAllY[i] = -999.;
    pvAllZ[i] = -999.;
    pvAllLogSumPtSq[i] = -999.;
    pvAllSumPx[i] = -999.;
    pvAllSumPy[i] = -999.;
  }
};


void displacedJetMuon_ntupler::resetPileUpBranches()
{
  nBunchXing = 0;
  for(int i = 0; i < MAX_NBX; i++)
  {
    BunchXing[i] = -999;
    nPU[i]       = -999;
    nPUmean[i]   = -999.;
  }
};

void displacedJetMuon_ntupler::resetMuonBranches()
{
  nMuons = 0;
  for(int i = 0; i < OBJECTARRAYSIZE; i++)
  {
    muonE[i] = 0.0;
    muonPt[i] = 0.0;
    muonEta[i] = 0.0;
    muonPhi[i] = 0.0;
    muonCharge[i] = -99;
    muonIsLoose[i] = false;
    muonIsMedium[i] = false;
    muonIsTight[i] = false;
    muon_d0[i] = -99.0;
    muon_dZ[i] = -99.0;
    muon_ip3d[i] = -99.0;
    muon_ip3dSignificance[i] = -99.0;
    muonType[i] = 0;
    muonQuality[i] = 0;
    muon_pileupIso[i] = -99.0;
    muon_chargedIso[i] = -99.0;
    muon_photonIso[i] = -99.0;
    muon_neutralHadIso[i] = -99.0;
    muon_ptrel[i] = -99.0;
    muon_chargedMiniIso[i] = -99.0;
    muon_photonAndNeutralHadronMiniIso[i] = -99.0;
    muon_chargedPileupMiniIso[i] = -99.0;
    muon_activityMiniIsoAnnulus[i] = -99.0;
    muon_passSingleMuTagFilter[i] = false;
    for (int q=0;q<MAX_MuonHLTFilters;q++) muon_passHLTFilter[i][q] = false;
    muon_validFractionTrackerHits[i] = -99.0;
    muon_isGlobal[i] = false;
    muon_normChi2[i] = -99.0;
    muon_chi2LocalPosition[i] = -99.0;
    muon_kinkFinder[i] = -99.0;
    muon_segmentCompatability[i] = -99.0;
    muonIsICHEPMedium[i] = false;
  }
};

void displacedJetMuon_ntupler::resetElectronBranches()
{
  nElectrons = 0;
  for(int i = 0; i < OBJECTARRAYSIZE; i++)
  {
    eleE[i] = 0.0;
    elePt[i] = 0.0;
    eleEta[i] = 0.0;
    elePhi[i] = 0.0;
    eleE_SC[i] = -99.0;
    eleEta_SC[i] = -99.0;
    elePhi_SC[i] = -99.0;
    eleSigmaIetaIeta[i] = -99.0;
    eleFull5x5SigmaIetaIeta[i] = -99.0;
    eleR9[i] = -99;
    ele_dEta[i] = -99;
    ele_dPhi[i] = -99;
    ele_HoverE[i] = -99;
    ele_d0[i] = -99;
    ele_dZ[i] = -99;
    ele_ip3d[i] = -99;
    ele_ip3dSignificance[i] = -99;
    ele_pileupIso[i] = -99.0;
    ele_chargedIso[i] = -99.0;
    ele_photonIso[i] = -99.0;
    ele_neutralHadIso[i] = -99.0;
    ele_MissHits[i] = -99;
    ele_PassConvVeto[i] = false;
    ele_OneOverEminusOneOverP[i] = -99.0;
    ele_IDMVAGeneralPurpose[i] = -99.0;
    ele_IDMVACategoryGeneralPurpose[i] = -1;
    ele_IDMVAHZZ[i] = -99.0;
    ele_IDMVACategoryHZZ[i] = -1;
    ele_RegressionE[i] = -99.0;
    ele_CombineP4[i] = -99.0;
    ele_ptrel[i] = -99.0;
    ele_chargedMiniIso[i] = -99.0;
    ele_photonAndNeutralHadronMiniIso[i] = -99.0;
    ele_chargedPileupMiniIso[i] = -99.0;
    ele_activityMiniIsoAnnulus[i] = -99.0;
    ele_passSingleEleTagFilter[i] = false;
    ele_passTPOneTagFilter[i] = false;
    ele_passTPTwoTagFilter[i] = false;
    ele_passTPOneProbeFilter[i] = false;
    ele_passTPTwoProbeFilter[i] = false;
    for (int q=0;q<MAX_ElectronHLTFilters;q++) ele_passHLTFilter[i][q] = false;
  }

};

void displacedJetMuon_ntupler::resetTauBranches()
{
  nTaus = 0;
  for(int i = 0; i < OBJECTARRAYSIZE; i++)
  {
    tauE[i] = -999.;
    tauPt[i] = -999.;
    tauEta[i] = -999.;
    tauPhi[i] = -999.;
    tau_IsLoose[i] = false;
    tau_IsMedium[i] = false;
    tau_IsTight[i] = false;
    tau_passEleVetoLoose[i] = false;
    tau_passEleVetoMedium[i] = false;
    tau_passEleVetoTight[i] = false;
    tau_passMuVetoLoose[i] = false;
    tau_passMuVetoMedium[i] = false;
    tau_passMuVetoTight[i] = false;
    tau_ID[i] = 0;//tauID Bits
    tau_combinedIsoDeltaBetaCorr3Hits[i] = -999.;;
    tau_chargedIsoPtSum[i] = -999.;;
    tau_neutralIsoPtSum[i] = -999.;;
    tau_puCorrPtSum[i] = -999.;;
    tau_eleVetoMVA[i] = -999.;;
    tau_eleVetoCategory[i] = -999.;;
    tau_muonVetoMVA[i] = -999.;;
    tau_isoMVAnewDMwLT[i] = -999.;;
    tau_isoMVAnewDMwoLT[i] = -999.;;
    tau_leadCandPt[i] = -999.;;
    tau_leadCandID[i] = -999.;;
    tau_leadChargedHadrCandPt[i] = -999.;;
    tau_leadChargedHadrCandID[i] = -999.;;
  }
};

void displacedJetMuon_ntupler::resetPhotonBranches()
{
  nPhotons = 0;
  nPhotons_overlap = 0;
  //Photon
  for(int i = 0; i < OBJECTARRAYSIZE; i++)
  {
    phoE[i] = 0.0;
    phoPt[i] = 0.0;
    phoEta[i] = 0.0;
    phoPhi[i] = 0.0;
    phoSigmaIetaIeta[i] = -99.0;
    phoFull5x5SigmaIetaIeta[i] = -99.0;
    phoR9[i] = -99.0;
    pho_sminor[i] = -99.0;
    pho_smajor[i] = -99.0;
    pho_HoverE[i] = -99.0;
    pho_sumChargedHadronPt[i] = -99.0;
    pho_sumNeutralHadronEt[i] = -99.0;
    pho_sumPhotonEt[i] = -99.0;
    pho_ecalPFClusterIso[i] = -99.0;
    pho_hcalPFClusterIso[i] = -99.0;
    pho_trkSumPtHollowConeDR03[i] = -99.0;
    pho_sumWorstVertexChargedHadronPt[i] = -99.0;
    pho_pfIsoChargedHadronIso[i] = -99.0;
    pho_pfIsoChargedHadronIsoWrongVtx[i] = -99.0;
    pho_pfIsoNeutralHadronIso[i] = -99.0;
    pho_pfIsoPhotonIso[i] = -99.0;
    pho_pfIsoModFrixione[i] = -99.0;
    pho_pfIsoSumPUPt[i] = -99.0;
    pho_isConversion[i] = false;
    pho_passEleVeto[i] = false;
    pho_RegressionE[i] = -99.0;
    pho_RegressionEUncertainty[i] = -99.0;
    pho_IDMVA[i] = -99.0;
    pho_superClusterEnergy[i] = -99.0;
    pho_superClusterRawEnergy[i] = -99.0;
    pho_superClusterEta[i]    = -99.0;
    pho_superClusterPhi[i]    = -99.0;
    pho_superClusterX[i]      = -99.0;
    pho_superClusterY[i]      = -99.0;
    pho_superClusterZ[i]      = -99.0;
    pho_hasPixelSeed[i] = false;
    for (int q=0;q<MAX_PhotonHLTFilters;q++) pho_passHLTFilter[i][q] = false;
    pho_convType[i] = -99;
    pho_convTrkZ[i] = -99.;
    pho_convTrkClusZ[i] = -99.;

    for (int ipv=0; ipv < MAX_NPV; ++ipv) {
      pho_sumChargedHadronPtAllVertices[i][ipv] = -99.0;
      pho_vtxSumPx[i][ipv] = 0.;
      pho_vtxSumPy[i][ipv] = 0.;
    }
    pho_isStandardPhoton[i] = true;
    pho_seedRecHitSwitchToGain6[i] = false;
    pho_seedRecHitSwitchToGain1[i] = false;
    pho_anyRecHitSwitchToGain6[i] = false;
    pho_anyRecHitSwitchToGain1[i] = false;
  }
  return;
};

void displacedJetMuon_ntupler::resetMuonSystemBranches()
{
  nCscSeg = 0;
  ncscRechits = 0;
  nCscSegClusters = 0;
  nCscRechitClusters = 0;
  nCscWireDigis = 0;
  nCscStripDigis = 0;
  for ( int i = 0; i < CSCRECHITARRAYSIZE; i++) {
    cscRechitsPhi[i] = 0.0;
    cscRechitsEta[i] = 0.0;
    cscRechitsX[i] = 0.0;
    cscRechitsY[i] = 0.0;
    cscRechitsZ[i] = 0.0;
    cscRechitsE[i] = 0.0;
    cscRechitsTwire[i] = 0.0;
    cscRechitsTpeak[i] = 0.0;
    cscRechitsQuality[i] = 0.0;
    cscRechitsChannels[i] = 0;
    cscRechitsNStrips[i] = 0;
    cscRechitsHitWire[i] = 0;
    cscRechitsWGroupsBX[i] = 0;
    cscRechitsNWireGroups[i] = 0;
    cscRechitsDetId[i] = 0;
    cscRechitsStation[i] = 0;
    cscRechitsChamber[i] = 0;
    cscSegPhi[i] = 0.0;
    cscSegEta[i] = 0.0;
    cscSegX[i] = 0.0;
    cscSegY[i] = 0.0;
    cscSegZ[i] = 0.0;
    cscSegDirectionX[i] = 0.0;
    cscSegDirectionY[i] = 0.0;
    cscSegDirectionZ[i] = 0.0;
    cscSegT[i] = 0.0;
    cscSegChi2[i] = 0.0;
    cscSegNRecHits[i] = 0;
    cscSegStation[i] = 0;
    cscSegChamber[i] = 0;
  }

  for ( int i = 0; i < OBJECTARRAYSIZE; i++) {
    cscSegCluster_match_gParticle_id[i] = 999;
    cscSegCluster_match_gParticle_index[i] = 999;
    cscSegCluster_match_gParticle_minDeltaR[i] = 999.;
    cscSegClusterJetVetoPt[i] = 0.0;
    cscSegClusterJetVetoE[i] = 0.0;
    cscSegClusterCaloJetVeto[i] = 0.0;
    cscSegClusterMuonVetoPt[i] = 0.0;
    cscSegClusterMuonVetoE[i] = 0.0;
    cscSegClusterX[i] = 0.0;
    cscSegClusterY[i] = 0.0;
    cscSegClusterZ[i] = 0.0;
    cscSegClusterTime[i] = 0.0;
    cscSegClusterTimeSpread[i] = 0.0;
    cscSegClusterGenMuonDeltaR[i] = 0.0;
    cscSegClusterMajorAxis[i] = 0.0;
    cscSegClusterMinorAxis[i] = 0.0;
    cscSegClusterXSpread[i] = 0.0;
    cscSegClusterYSpread[i] = 0.0;
    cscSegClusterZSpread[i] = 0.0;
    cscSegClusterEtaPhiSpread[i] = 0.0;
    cscSegClusterEtaSpread[i] = 0.0;
    cscSegClusterPhiSpread[i] = 0.0;
    cscSegClusterEta[i] = 0.0;
    cscSegClusterPhi[i] = 0.0;
    cscSegClusterSize[i] = 0;
    cscSegClusterNSegmentChamberPlus11[i] = 0;
    cscSegClusterNSegmentChamberPlus12[i] = 0;
    cscSegClusterNSegmentChamberPlus13[i] = 0;
    cscSegClusterNSegmentChamberPlus21[i] = 0;
    cscSegClusterNSegmentChamberPlus22[i] = 0;
    cscSegClusterNSegmentChamberPlus31[i] = 0;
    cscSegClusterNSegmentChamberPlus32[i] = 0;
    cscSegClusterNSegmentChamberPlus41[i] = 0;
    cscSegClusterNSegmentChamberPlus42[i] = 0;
    cscSegClusterNSegmentChamberMinus11[i] = 0;
    cscSegClusterNSegmentChamberMinus12[i] = 0;
    cscSegClusterNSegmentChamberMinus13[i] = 0;
    cscSegClusterNSegmentChamberMinus21[i] = 0;
    cscSegClusterNSegmentChamberMinus22[i] = 0;
    cscSegClusterNSegmentChamberMinus31[i] = 0;
    cscSegClusterNSegmentChamberMinus32[i] = 0;
    cscSegClusterNSegmentChamberMinus41[i] = 0;
    cscSegClusterNSegmentChamberMinus42[i] = 0;
    cscSegClusterMe11Ratio[i] = 0.0;
    cscSegClusterMe12Ratio[i] = 0.0;
    cscSegClusterVertexR[i] = 0.0;
    cscSegClusterVertexZ[i] = 0.0;
    cscSegClusterVertexN[i] = 0;
    cscSegClusterVertexN1[i] = 0;
    cscSegClusterVertexN5[i] = 0;
    cscSegClusterVertexN10[i] = 0;
    cscSegClusterVertexN15[i] = 0;
    cscSegClusterVertexN20[i] = 0;
    cscSegClusterVertexChi2[i] = 0.0;
    cscSegClusterVertexDis[i] = 0.0;
    cscSegClusterMaxStationRatio[i] = 0.0;
    cscSegClusterMaxStation[i] = 0;
    cscSegClusterNStation[i] = 0;
    cscSegClusterMaxChamberRatio[i] = 0.0;
    cscSegClusterMaxChamber[i] = 0;
    cscSegClusterNChamber[i] = 0;

    cscRechitCluster_match_cscSegCluster_index[i] = 999;
    cscRechitCluster_match_cscSegCluster_minDeltaR[i] = 999.;
    cscRechitCluster_match_gParticle_id[i] = 999;
    cscRechitCluster_match_gParticle_index[i] = 999;
    cscRechitCluster_match_gParticle_minDeltaR[i] = 999.;
    cscRechitClusterJetVetoPt[i] = 0.0;
    cscRechitClusterJetVetoE[i] = 0.0;
    cscRechitClusterCaloJetVeto[i] = 0.0;
    cscRechitClusterMuonVetoPt[i] = 0.0;
    cscRechitClusterMuonVetoE[i] = 0.0;
    cscRechitClusterX[i] = 0.0;
    cscRechitClusterY[i] = 0.0;
    cscRechitClusterZ[i] = 0.0;
    cscRechitClusterTime[i] = 0.0;
    cscRechitClusterTimeSpread[i] = 0.0;
    cscRechitClusterGenMuonDeltaR[i] = 0.0;
    cscRechitClusterMajorAxis[i] = 0.0;
    cscRechitClusterMinorAxis[i] = 0.0;
    cscRechitClusterXSpread[i] = 0.0;
    cscRechitClusterYSpread[i] = 0.0;
    cscRechitClusterZSpread[i] = 0.0;
    cscRechitClusterEtaPhiSpread[i] = 0.0;
    cscRechitClusterEtaSpread[i] = 0.0;
    cscRechitClusterPhiSpread[i] = 0.0;
    cscRechitClusterEta[i] = 0.0;
    cscRechitClusterPhi[i] = 0.0;
    cscRechitClusterSize[i] = 0;
    cscRechitClusterNRechitChamberPlus11[i] = 0;
    cscRechitClusterNRechitChamberPlus12[i] = 0;
    cscRechitClusterNRechitChamberPlus13[i] = 0;
    cscRechitClusterNRechitChamberPlus21[i] = 0;
    cscRechitClusterNRechitChamberPlus22[i] = 0;
    cscRechitClusterNRechitChamberPlus31[i] = 0;
    cscRechitClusterNRechitChamberPlus32[i] = 0;
    cscRechitClusterNRechitChamberPlus41[i] = 0;
    cscRechitClusterNRechitChamberPlus42[i] = 0;
    cscRechitClusterNRechitChamberMinus11[i] = 0;
    cscRechitClusterNRechitChamberMinus12[i] = 0;
    cscRechitClusterNRechitChamberMinus13[i] = 0;
    cscRechitClusterNRechitChamberMinus21[i] = 0;
    cscRechitClusterNRechitChamberMinus22[i] = 0;
    cscRechitClusterNRechitChamberMinus31[i] = 0;
    cscRechitClusterNRechitChamberMinus32[i] = 0;
    cscRechitClusterNRechitChamberMinus41[i] = 0;
    cscRechitClusterNRechitChamberMinus42[i] = 0;
    cscRechitClusterMe11Ratio[i] = 0.0;
    cscRechitClusterMe12Ratio[i] = 0.0;
    cscRechitClusterVertexR[i] = 0.0;
    cscRechitClusterVertexZ[i] = 0.0;
    cscRechitClusterVertexN[i] = 0;
    cscRechitClusterVertexN1[i] = 0;
    cscRechitClusterVertexN5[i] = 0;
    cscRechitClusterVertexN10[i] = 0;
    cscRechitClusterVertexN15[i] = 0;
    cscRechitClusterVertexN20[i] = 0;
    cscRechitClusterVertexChi2[i] = 0.0;
    cscRechitClusterVertexDis[i] = 0.0;
    cscRechitClusterMaxStationRatio[i] = 0.0;
    cscRechitClusterMaxStation[i] = 0;
    cscRechitClusterNStation[i] = 0;
    cscRechitClusterMaxChamberRatio[i] = 0.0;
    cscRechitClusterMaxChamber[i] = 0;
    cscRechitClusterNChamber[i] = 0;
  }
  
  nRpc = 0;
  for ( int i = 0; i < OBJECTARRAYSIZE; i++) {
    rpcPhi[i] = 0.0;
    rpcEta[i] = 0.0;
    rpcX[i] = 0.0;
    rpcY[i] = 0.0;
    rpcZ[i] = 0.0;
    rpcT[i] = 0.0;
    rpcBx[i] = 0;
    rpcTError[i] = 0.0;
  }

  nDtSeg = 0;
  nDtRechits = 0;
  nDtSegClusters = 0;
  nDtRechitClusters = 0;
  for ( int i = 0; i < CSCRECHITARRAYSIZE; i++) {
    dtRechitX[i] = 0.0;
    dtRechitY[i] = 0.0;
    dtRechitZ[i] = 0.0;
    dtRechitEta[i] = 0.0;
    dtRechitPhi[i] = 0.0;
    dtRechitTime[i] = 0.0;
    dtRechitStation[i] = 0;
    dtRechitWheel[i] = 0;

    dtSegPhi[i] = 0.0;
    dtSegEta[i] = 0.0;
    dtSegX[i] = 0.0;
    dtSegY[i] = 0.0;
    dtSegZ[i] = 0.0;
    dtSegStation[i] = 0;
    dtSegWheel[i] = 0;
    dtSegTime[i] = -9999.0;
    dtSegTimeError[i] = -9999.0;
  }
  
  for ( int i = 0; i < OBJECTARRAYSIZE; i++) {
    dtRechitCluster_match_gParticle_id[i] = 999;
    dtRechitCluster_match_gParticle_index[i] = 999;
    dtRechitCluster_match_gParticle_minDeltaR[i] = 999.;
    dtRechitClusterJetVetoPt[i] = 0.0;
    dtRechitClusterJetVetoE[i] = 0.0;
    dtRechitClusterCaloJetVeto[i] = 0.0;
    dtRechitClusterMuonVetoPt[i] = 0.0;
    dtRechitClusterMuonVetoE[i] = 0.0;
    dtRechitClusterX[i] = 0.0;
    dtRechitClusterY[i] = 0.0;
    dtRechitClusterZ[i] = 0.0;
    dtRechitClusterTime[i] = 0.0;
    dtRechitClusterTimeSpread[i] = 0.0;
    dtRechitClusterGenMuonDeltaR[i] = 0.0;
    dtRechitClusterMajorAxis[i] = 0.0;
    dtRechitClusterMinorAxis[i] = 0.0;
    dtRechitClusterXSpread[i] = 0.0;
    dtRechitClusterYSpread[i] = 0.0;
    dtRechitClusterZSpread[i] = 0.0;
    dtRechitClusterEtaPhiSpread[i] = 0.0;
    dtRechitClusterEtaSpread[i] = 0.0;
    dtRechitClusterPhiSpread[i] = 0.0;
    dtRechitClusterEta[i] = 0.0;
    dtRechitClusterPhi[i] = 0.0;
    dtRechitClusterSize[i] = 0;
    dtRechitClusterNSegmentStation1[i] = 0;
    dtRechitClusterNSegmentStation2[i] = 0;
    dtRechitClusterNSegmentStation3[i] = 0;
    dtRechitClusterNSegmentStation4[i] = 0;
    dtRechitClusterMaxStationRatio[i] = 0.0;
    dtRechitClusterMaxStation[i] = 0;
    dtRechitClusterNStation[i] = 0;
    dtRechitClusterMaxChamberRatio[i] = 0.0;
    dtRechitClusterMaxChamber[i] = 0;
    dtRechitClusterNChamber[i] = 0;    
 
    dtSegCluster_match_gParticle_id[i] = 999;
    dtSegCluster_match_gParticle_index[i] = 999;
    dtSegCluster_match_gParticle_minDeltaR[i] = 999.;
    dtSegClusterJetVetoPt[i] = 0.0;
    dtSegClusterJetVetoE[i] = 0.0;
    dtSegClusterCaloJetVeto[i] = 0.0;
    dtSegClusterMuonVetoPt[i] = 0.0;
    dtSegClusterMuonVetoE[i] = 0.0;
    dtSegClusterX[i] = 0.0;
    dtSegClusterY[i] = 0.0;
    dtSegClusterZ[i] = 0.0;
    dtSegClusterTime[i] = 0.0;
    dtSegClusterTimeSpread[i] = 0.0;
    dtSegClusterGenMuonDeltaR[i] = 0.0;
    dtSegClusterMajorAxis[i] = 0.0;
    dtSegClusterMinorAxis[i] = 0.0;
    dtSegClusterXSpread[i] = 0.0;
    dtSegClusterYSpread[i] = 0.0;
    dtSegClusterZSpread[i] = 0.0;
    dtSegClusterEtaPhiSpread[i] = 0.0;
    dtSegClusterEtaSpread[i] = 0.0;
    dtSegClusterPhiSpread[i] = 0.0;
    dtSegClusterEta[i] = 0.0;
    dtSegClusterPhi[i] = 0.0;
    dtSegClusterSize[i] = 0;
    dtSegClusterNSegmentStation1[i] = 0;
    dtSegClusterNSegmentStation2[i] = 0;
    dtSegClusterNSegmentStation3[i] = 0;
    dtSegClusterNSegmentStation4[i] = 0;
    dtSegClusterMaxStationRatio[i] = 0.0;
    dtSegClusterMaxStation[i] = 0;
    dtSegClusterNStation[i] = 0;
    dtSegClusterMaxChamberRatio[i] = 0.0;
    dtSegClusterMaxChamber[i] = 0;
    dtSegClusterNChamber[i] = 0;
  }

  return;
};
void displacedJetMuon_ntupler::resetJetBranches()
{
  nJets = 0;
  for ( int i = 0; i < OBJECTARRAYSIZE; i++)
  {
    jetE[i] = 0.0;
    jetPt[i] = 0.0;
    jetEta[i] = 0.0;
    jetPhi[i] = 0.0;
    jetCSV[i] = 0.0;
    jetCISV[i] = 0.0;
    jetMass[i] =  -99.0;
    jetJetArea[i] = -99.0;
    jetPileupE[i] = -99.0;
    jetPileupId[i] = -99.0;
    jetPileupIdFlag[i] = -1;
    jetPassIDLoose[i] = false;
    jetPassIDTight[i] = false;
    jetPassMuFrac[i] = false;
    jetPassEleFrac[i] = false;
    jetPartonFlavor[i] = 0;
    jetHadronFlavor[i] = 0;
    jetChargedEMEnergyFraction[i] = -99.0;
    jetNeutralEMEnergyFraction[i] = -99.0;
    jetChargedHadronEnergyFraction[i] = -99.0;
    jetNeutralHadronEnergyFraction[i] = -99.0;
    jetMuonEnergyFraction[i] = -99.0;
    jetHOEnergyFraction[i] = -99.0;
    jetHFHadronEnergyFraction[i] = -99.0;
    jetHFEMEnergyFraction[i] = -99.0;
    jetAllMuonPt[i] = 0.0;
    jetAllMuonEta[i] = 0.0;
    jetAllMuonPhi[i] = 0.0;
    jetAllMuonM[i] = 0.0;
    jetPtWeightedDZ[i] = 0.0;
    jetNRechits[i] = 0;
    jetRechitE[i] = 0.0;
    jetRechitT[i] = 0.0;
    jetRechitT_rms[i] = 0.0;

    jetRechitE_Error[i] = 0.0;
    jetRechitT_Error[i] = 0.0;
    jetGammaMax[i] = -99.0;
    jetGammaMax_EM[i] = -99.0;
    jetGammaMax_Hadronic[i] = -99.0;
    jetGammaMax_ET[i] = -99.0;
    jetAlphaMax[i] = -99.0;
    jetBetaMax[i] = -99.0;
    jetPtAllTracks[i] = -99.0;
    jetPtAllPVTracks[i] = -99.0;
    jetMedianTheta2D[i] = -99.0;
    jetMedianIP[i] = -99.0;
    jetMinDeltaRAllTracks[i] =-99.0;
    jetMinDeltaRPVTracks[i] = -99.0;
    jet_sig_et1[i] = -99.0;
    jet_sig_et2[i] = -99.0;
    jet_energy_frac[i] = 0.0;

    jetGammaMax_wp[i] = -99.0;
    jetGammaMax_EM_wp[i] = -99.0;
    jetGammaMax_Hadronic_wp[i] = -99.0;
    jetGammaMax_ET_wp[i] = -99.0;
    jetAlphaMax_wp[i] = -99.0;
    jetBetaMax_wp[i] = -99.0;
    jetPtAllTracks_wp[i] = -99.0;
    jetPtAllPVTracks_wp[i] = -99.0;
    jetMedianTheta2D_wp[i] = -99.0;
    jetMedianIP_wp[i] = -99.0;
    jetMinDeltaRAllTracks_wp[i] =-99.0;
    jetMinDeltaRPVTracks_wp[i] = -99.0;
  }

  nFatJets = 0;
  for ( int i = 0; i < OBJECTARRAYSIZE; i++) {
    fatJetE[i] = 0.0;
    fatJetPt[i] = 0.0;
    fatJetEta[i] = 0.0;
    fatJetPhi[i] = 0.0;
    fatJetCorrectedPt[i] = 0.0;
    fatJetCorrectedEta[i] = 0.0;
    fatJetCorrectedPhi[i] = 0.0;
    fatJetPrunedM[i] = 0.0;
    fatJetTrimmedM[i] = 0.0;
    fatJetFilteredM[i] = 0.0;
    fatJetSoftDropM[i] = 0.0;
    fatJetCorrectedSoftDropM[i] = 0.0;
    fatJetUncorrectedSoftDropM[i] = 0.0;
    fatJetTau1[i] = 0.0;
    fatJetTau2[i] = 0.0;
    fatJetTau3[i] = 0.0;     
  }

  return;
};
void displacedJetMuon_ntupler::resetHORechitBranches()
{
  nHORechits = 0;
  for ( int i = 0; i < HORECHITARRAYSIZE; i++)
  {

	hoRechit_Phi[i] = -999.;
	hoRechit_Eta[i] = -999.;
	hoRechit_X[i] = -999.;
	hoRechit_Y[i] = -999.;
	hoRechit_Z[i] = -999.;
	hoRechit_E[i] = -999.;
	hoRechit_T[i] = -999.;
  }
  return;
};
void displacedJetMuon_ntupler::resetEcalRechitBranches()
{
  nRechits = 0;
  for ( int i = 0; i < RECHITARRAYSIZE; i++) {
    ecalRechit_Eta[i] = -999.;
    ecalRechit_Phi[i] = -999.;
    ecalRechit_E[i] = -999.;
    ecalRechit_T[i] = -999.;
    ecalRechit_E_Error[i] = -999.;
    ecalRechit_T_Error[i] = -999.;
    ecalRechit_kSaturatedflag[i] = false;
    ecalRechit_kLeadingEdgeRecoveredflag[i] = false;
    ecalRechit_kPoorRecoflag[i] = false;
    ecalRechit_kWeirdflag[i] = false;
    ecalRechit_kDiWeirdflag[i] = false;
  }

  nHBHERechits = 0;
  for ( int i = 0; i < RECHITARRAYSIZE; i++) {
    hbheRechit_Phi[i] = -999.;
    hbheRechit_Eta[i] = -999.;
    hbheRechit_X[i] = -999.;
    hbheRechit_Y[i] = -999.;
    hbheRechit_Z[i] = -999.;
    hbheRechit_E[i] = -999.;
    hbheRechit_T[i] = -999.;
    hbheRechit_iEta[i] = -999;
    hbheRechit_iPhi[i] = -999;
    hbheRechit_depth[i] = -999;
  }

  nTracks = 0;
  for ( int i = 0; i < RECHITARRAYSIZE; i++) {
    track_Pt[i] = -999.;
    track_Eta[i] = -999.;
    track_Phi[i] = -999.;
    track_bestVertexIndex[i] = -1;
    track_nMissingInnerHits[i] = -1;
    track_nMissingOuterHits[i] = -1;
    track_angle[i] = -999.;
    track_dxyToBS[i] = -999.;
    track_dxyErr[i] = -999.;
    track_dzToPV[i] = -999.;
    track_dzErr[i] = -999.;
  }

  return;
};


void displacedJetMuon_ntupler::resetMetBranches()
{
  metPt = -999;
  metPhi = -999;
  sumMET = -99.0;
  UncMETdpx = -99.0;
  UncMETdpy = -99.0;
  UncMETdSumEt = -99.0;
  metType0Pt = -99.0;
  metType0Phi = -99.0;
  metType1Pt_raw = -99.0;
  metType1Pt = -99.0;
  metType1Px = -99.0;
  metType1Py = -99.0;
  metType1Eta = -99.0;
  metType1Phi = -99.0;
  metType1Phi_raw = -99.0;
  metType0Plus1Pt = -99.0;
  metType0Plus1Phi = -99.0;
  metPtRecomputed = -99.0;
  metPhiRecomputed = -99.0;
  metNoHFPt = -99.0;
  metNoHFPhi = -99.0;
  metPuppiPt = -99.0;
  metPuppiPhi = -99.0;
  metCaloPt = -999;
  metCaloPhi = -999;
  Flag_HBHENoiseFilter = false;
  Flag_HBHETightNoiseFilter = false;
  Flag_HBHEIsoNoiseFilter = false;
  Flag_badChargedCandidateFilter = false;
  Flag_badMuonFilter = false;
  Flag_badGlobalMuonFilter = false;
  Flag_duplicateMuonFilter = false;
  Flag_globalTightHalo2016Filter = false;
  Flag_globalSuperTightHalo2016Filter = false;
  Flag_hcalLaserEventFilter = false;
  Flag_EcalDeadCellTriggerPrimitiveFilter = false;
  Flag_EcalDeadCellBoundaryEnergyFilter = false;
  Flag_goodVertices = false;
  Flag_trackingFailureFilter = false;
  Flag_eeBadScFilter = false;
  Flag_ecalLaserCorrFilter = false;
  Flag_trkPOGFilters = false;
  Flag_trkPOG_manystripclus53X = false;
  Flag_trkPOG_toomanystripclus53X = false;
  Flag_trkPOG_logErrorTooManyClusters = false;
  Flag_BadPFMuonFilter = false;
  Flag_BadChargedCandidateFilter = false;
  Flag_ecalBadCalibFilter = false;
  Flag_METFilters = false;

  Flag2_globalSuperTightHalo2016Filter = false;
  Flag2_globalTightHalo2016Filter = false;
  Flag2_goodVertices = false;
  Flag2_BadChargedCandidateFilter = false;
  Flag2_BadPFMuonFilter   = false;
  Flag2_EcalDeadCellTriggerPrimitiveFilter  = false;
  Flag2_HBHENoiseFilter  = false;
  Flag2_HBHEIsoNoiseFilter  = false;
  Flag2_ecalBadCalibFilter  = false;
  Flag2_eeBadScFilter  = false;


  metType1PtJetResUp=-999.;
  metType1PtJetResDown=-999.;
  metType1PtJetEnUp=-999.;
  metType1PtJetEnDown=-999.;
  metType1PtMuonEnUp=-999.;
  metType1PtMuonEnDown=-999.;
  metType1PtElectronEnUp=-999.;
  metType1PtElectronEnDown=-999.;
  metType1PtTauEnUp=-999.;
  metType1PtTauEnDown=-999.;
  metType1PtUnclusteredEnUp=-999.;
  metType1PtUnclusteredEnDown=-999.;
  metType1PtPhotonEnUp=-999.;
  metType1PtPhotonEnDown=-999.;
  metType1PtMETUncertaintySize=-999.;
  metType1PtJetResUpSmear=-999.;
  metType1PtJetResDownSmear=-999.;
  metType1PtMETFullUncertaintySize=-999.;

  metType1PhiJetResUp=-999.;
  metType1PhiJetResDown=-999.;
  metType1PhiJetEnUp=-999.;
  metType1PhiJetEnDown=-999.;
  metType1PhiMuonEnUp=-999.;
  metType1PhiMuonEnDown=-999.;
  metType1PhiElectronEnUp=-999.;
  metType1PhiElectronEnDown=-999.;
  metType1PhiTauEnUp=-999.;
  metType1PhiTauEnDown=-999.;
  metType1PhiUnclusteredEnUp=-999.;
  metType1PhiUnclusteredEnDown=-999.;
  metType1PhiPhotonEnUp=-999.;
  metType1PhiPhotonEnDown=-999.;
  metType1PhiMETUncertaintySize=-999.;
  metType1PhiJetResUpSmear=-999.;
  metType1PhiJetResDownSmear=-999.;
  metType1PhiMETFullUncertaintySize=-999.;
};
void displacedJetMuon_ntupler::resetGenParticleBranches()
{
  for(int i = 0; i < GENPARTICLEARRAYSIZE; i++)
  {
    gParticleMotherId[i] = -99999;
    gParticleMotherIndex[i] = -99999;
    gParticleId[i] = -99999;
    gParticleStatus[i] = -99999;
    gParticleE[i] = -99999.0;
    gParticlePt[i] = -99999.0;
    gParticlePx[i] = -99999.0;
    gParticlePy[i] = -99999.0;
    gParticlePz[i] = -99999.0;
    gParticleEta[i] = -99999.0;
    gParticlePhi[i] = -99999.0;

    gParticleProdVertexX[i] = -99999.0;
    gParticleProdVertexY[i] = -99999.0;
    gParticleProdVertexZ[i] = -99999.0;
    gParticleDecayVertexX[i] = -99999.0;
    gParticleDecayVertexY[i] = -99999.0;
    gParticleDecayVertexZ[i] = -99999.0;

  }
  for ( int i = 0; i < LLP_ARRAY_SIZE; i++ )
  {
    gLLP_prod_vertex_x[i] = -666.;
    gLLP_prod_vertex_y[i] = -666.;
    gLLP_prod_vertex_z[i] = -666.;
    gLLP_decay_vertex_x[i] = -666.;
    gLLP_decay_vertex_y[i] = -666.;
    gLLP_decay_vertex_z[i] = -666.;
    gLLP_beta[i] = -666.;
    gLLP_pt[i] = -666.;
    gLLP_eta[i] = -666.;
    gLLP_e[i] = -666.;
    gLLP_phi[i] = -666.;
    gLLP_csc[i] = false;
    gLLP_dt[i] = false;
    gLLP_travel_time[i] = -666.;
  }

  for ( int i = 0; i < LLP_DAUGHTER_ARRAY_SIZE; i++ )
  {
    gLLP_daughter_id[i] = 0;
    gLLP_daughter_pt[i] = -666.;
    gLLP_daughter_eta[i] = -666.;
    gLLP_daughter_phi[i] = -666.;
    gLLP_daughter_eta_ecalcorr[i] = -666.;
    gLLP_daughter_phi_ecalcorr[i] = -666.;
    gLLP_daughter_e[i] = -666.;
    gLLP_daughter_mass[i] = -666.;

    gLLP_daughter_travel_time[i] = -666.;
    gen_time[i] = -666.;
    gen_time_pv[i] = -666.;
    photon_travel_time[i] = -666.;
    photon_travel_time_pv[i] = -666.;
}

 //grandaughters
  for ( int i = 0; i < LLP_GRAND_DAUGHTER_ARRAY_SIZE; i++ ) {
    // gLLP_grandaughter_EB[i] = false;
    // gLLP_grandaughter_ETL[i] = false;

    // gLLP_grandaughter_photon_travel_time_EB[i] = -666.;
    // gLLP_grandaughter_photon_travel_time_ETL[i] = -666.;

    // gLLP_grandaughter_travel_time_EB[i] = -666.;
    // gLLP_grandaughter_travel_time_ETL[i] = -666.;

    // gen_time_grandaughter_EB[i] = -666.;
    // gen_time_grandaughter_ETL[i] = -666.;

    gLLP_grandaughter_id[i] = 0;
    gLLP_grandaughter_pt[i] = -666.;
    gLLP_grandaughter_eta[i] = -666.;
    gLLP_grandaughter_phi[i] = -666.;
    gLLP_grandaughter_eta_ecalcorr[i] = -666.;
    gLLP_grandaughter_phi_ecalcorr[i] = -666.;
    gLLP_grandaughter_e[i] = -666.;
    gLLP_grandaughter_mass[i] = -666.;
  }

  return;
};

void displacedJetMuon_ntupler::resetMCBranches()
{
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

  scaleWeights->clear();
  pdfWeights->clear();
  alphasWeights->clear();

  return;
};

void displacedJetMuon_ntupler::resetTriggerBranches()
{
  for( int i = 0; i < NTriggersMAX; i++ )
  {
    triggerDecision[i] = false;
    triggerHLTPrescale[i] = 0;
  }
};
//------ Method called for each run ------//

void displacedJetMuon_ntupler::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {

  //read LHE header if present and determine which weights to read for pdf and alphas uncertainties based on the 
  //central pdf set used
  //This is semi-hardcoded for now to work with current centrally produced samples
  //covering SUSY signal samples, LO madgraph, NLO madgraph_aMC@NLO, and NLO powheg
  //generated with nnpdf30
  //More robust selection will require some parsing of <initrwgt> block
  
  if (!isData_) {
  
      //hardcode this for now. All SUSY signals use 5 flavor scheme, LO madgraph
      //to do it properly, we would need to parse the LHE header
      int pdfidx = 263000; 
            
      if (pdfidx == 263000) {
	//NNPDF30_lo_as_0130 (nf5) for LO madgraph samples and SUSY signals
	//pdfweightshelper.Init(100,60,edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/NNPDF30_lo_as_0130_hessian_60.csv"));
	firstPdfWeight = 11;
	lastPdfWeight = 110;
	firstAlphasWeight = -1;
	lastAlphasWeight = -1;      
      }
      else if (pdfidx == 263400) {
	//NNPdf30_lo_as_0130_nf4 for LO madgraph samples
	//pdfweightshelper.Init(100,60,edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/NNPDF30_lo_as_0130_nf_4_hessian_60.csv"));
	firstPdfWeight = 112;
	lastPdfWeight = 211;
	firstAlphasWeight = -1;
	lastAlphasWeight = -1;            
      }
      else if (pdfidx == 260000 || pdfidx == -1) {
	//NNPdf30_nlo_as_0118 (nf5) for NLO powheg samples
	//(work around apparent bug in current powheg samples by catching "-1" as well)
	//pdfweightshelper.Init(100,60,edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/NNPDF30_nlo_as_0118_hessian_60.csv"));
	firstPdfWeight = 10;
	lastPdfWeight = 109;
	firstAlphasWeight = 110;
	lastAlphasWeight = 111; 
      }
      else if (pdfidx == 292200) {
	//NNPdf30_nlo_as_0118 (nf5) with built-in alphas variations for NLO aMC@NLO samples
	//pdfweightshelper.Init(100,60,edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/NNPDF30_nlo_as_0118_hessian_60.csv"));
	firstPdfWeight = 10;
	lastPdfWeight = 109;
	firstAlphasWeight = 110;
	lastAlphasWeight = 111; 
      }   
      else if (pdfidx == 292000) {
	//NNPdf30_nlo_as_0118_nf4 with built-in alphas variations for NLO aMC@NLO samples
	//pdfweightshelper.Init(100,60,edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/NNPDF30_nlo_as_0118_nf_4_hessian_60.csv"));
	firstPdfWeight = 10;
	lastPdfWeight = 109;
	firstAlphasWeight = 110;
	lastAlphasWeight = 111;
      }
      else {
	firstPdfWeight = -1;
	lastPdfWeight = -1;
	firstAlphasWeight = -1;
	lastAlphasWeight = -1;
      }    
  } else {
    firstPdfWeight = -1;
    lastPdfWeight = -1;
    firstAlphasWeight = -1;
    lastAlphasWeight = -1;
  }      
  


}


//------ Method called for each lumi block ------//
void displacedJetMuon_ntupler::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
  if (useGen_) {  
    iLumi.getByToken(genLumiHeaderToken_,genLumiHeader);
  }
   
  //fill lhe comment lines with SUSY model parameter information
  lheComments = "";
  if (genLumiHeader.isValid()) {
    lheComments = genLumiHeader->configDescription();    
  }    
  cout << "LHEComments: " << lheComments << "\n";
  
}


//------ Method called for each event ------//
//------ Method called for each event ------//

void displacedJetMuon_ntupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;

  //initialize
  loadEvent(iEvent); //loads objects and resets tree branches
  NEvents->Fill(0); //increment event count

  //resetting output tree branches

  resetBranches();
  fillEventInfo(iEvent);
  fillPVAll();  
  fillMuons(iEvent);
  fillElectrons(iEvent);
  fillPhotons(iEvent, iSetup);  
  fillJets(iSetup);
  fillMet(iEvent);
  
  if (!isData) {
    fillPileUp();
    fillMC();
    fillGenParticles();
  }
  //fillHOSystem(iEvent,iSetup);
  fillMuonSystem(iEvent, iSetup);
  if ( enableTriggerInfo_ ) fillTrigger( iEvent );
  displacedJetMuonTree->Fill();

};


//------ Method called once each job just before starting event loop ------//
void displacedJetMuon_ntupler::beginJob()
{
  setBranches();
}

//------ Method called once each job just after ending the event loop ------//
void displacedJetMuon_ntupler::endJob(){};



//------ Methods to fill tree variables ------//
bool displacedJetMuon_ntupler::fillEventInfo(const edm::Event& iEvent)
{
  isData = isData_;
  runNum = iEvent.id().run();
  lumiNum = iEvent.luminosityBlock();
  eventNum = iEvent.id().event();
  eventTime = iEvent.eventAuxiliary().time().unixTime();

  //number of slimmedSecondaryVertices
  nSlimmedSecondV = secondaryVertices->size();

  //select the primary vertex, if any
  nPV = 0;
  myPV = &(vertices->front());
  bool foundPV = false;
  for(unsigned int i = 0; i < vertices->size(); i++)
  {
    if(vertices->at(i).isValid() && !vertices->at(i).isFake())
    {
      if (!foundPV)
      {
        myPV = &(vertices->at(i));
        foundPV = true;
      }
      nPV++;
    }
  }

  pvX = myPV->x();
  pvY = myPV->y();
  pvZ = myPV->z();

  //get rho
  fixedGridRhoFastjetAll = *rhoFastjetAll;

  return true;
};


bool displacedJetMuon_ntupler::fillHOSystem(const edm::Event& iEvent, const edm::EventSetup& iSetup){
    edm::ESHandle<CaloGeometry> geoHandle;
    iSetup.get<CaloGeometryRecord>().get(geoHandle);
    const CaloSubdetectorGeometry *hoGeometry = geoHandle->getSubdetectorGeometry(DetId::Hcal, HcalOuter);
    for (HORecHitCollection::const_iterator recHit=hcalRecHitsHO->begin(); recHit!=hcalRecHitsHO->end(); recHit++) {
	if (recHit->energy() < 1.5) continue;
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
    }
    return true;
}
bool displacedJetMuon_ntupler::fillMuonSystem(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::ESHandle<CSCGeometry> cscG;
  edm::ESHandle<DTGeometry> dtG;
  edm::ESHandle<RPCGeometry> rpcG;

  iSetup.get<MuonGeometryRecord>().get(cscG);
  iSetup.get<MuonGeometryRecord>().get(dtG);
  iSetup.get<MuonGeometryRecord>().get(rpcG);

  //*****************
  //** DIGIS
  //*****************
  if (readMuonDigis_) {
    nCscStripDigis = 0;
    CSCStripDigiCollection::DigiRangeIterator stripDetIt;
    for (stripDetIt = MuonCSCStripDigi->begin(); stripDetIt != MuonCSCStripDigi->end(); stripDetIt++){
      // const CSCDetId &id = (*stripDetIt).first;
      // int tempDetId = CSCDetId::rawIdMaker(CSCDetId::endcap(id), CSCDetId::station(id), CSCDetId::ring(id), CSCDetId::chamber(id), CSCDetId::layer(id));

      const CSCStripDigiCollection::Range &range = (*stripDetIt).second;
      for (CSCStripDigiCollection::const_iterator digiIt = range.first; digiIt != range.second; ++digiIt) {
        //ADC count is 8 time stamps
        std::vector<int> myADCVals = digiIt->getADCCounts();
	bool thisStripFired = false;
	float thisPedestal = 0.5 * (float)(myADCVals[0] + myADCVals[1]);
	float threshold = STRIP_DIGI_THRESHOLD;
	float diff = 0.;
	for (unsigned int iCount = 0; iCount < myADCVals.size(); iCount++) {
	  diff = (float)myADCVals[iCount] - thisPedestal;
	  if (diff > threshold) {
	    thisStripFired = true;
	  }
	}
	if (thisStripFired) nCscStripDigis++;
      }
    }
    nCscWireDigis = 0;
    CSCWireDigiCollection::DigiRangeIterator wireDetIt;
    for (wireDetIt = MuonCSCWireDigi->begin(); wireDetIt != MuonCSCWireDigi->end(); wireDetIt++){
      const CSCWireDigiCollection::Range &range = (*wireDetIt).second;
      for (CSCWireDigiCollection::const_iterator digiIt = range.first; digiIt != range.second; ++digiIt) {

        nCscWireDigis++;
      }

    }
  }







  //*****************
  //** CSC Segments
  //*****************
  vector<Point> points; //vector defined for DBScan Clustering algorithm
  for (const CSCSegment cscSegment : *cscSegments) {

    CSCDetId id  = (CSCDetId)(cscSegment).cscDetId();

    int endcap = CSCDetId::endcap(id) == 1 ? 1 : -1;
    LocalPoint segPos = (cscSegment).localPosition();
    LocalVector segDirection = (cscSegment).localDirection();
    const CSCChamber* cscchamber = cscG->chamber(id);
    if (cscchamber) {
      GlobalPoint globalPosition = cscchamber->toGlobal(segPos);
      GlobalVector globalDirection = cscchamber->toGlobal(segDirection);

      cscSegNRecHits[nCscSeg] = cscSegment.nRecHits();
      cscSegX[nCscSeg] = globalPosition.x();
      cscSegY[nCscSeg] = globalPosition.y();
      cscSegZ[nCscSeg] = globalPosition.z();
      cscSegPhi[nCscSeg] =  globalPosition.phi();
      cscSegEta[nCscSeg] =globalPosition.eta();
      cscSegDirectionX[nCscSeg] = globalDirection.x();
      cscSegDirectionY[nCscSeg] = globalDirection.y();
      cscSegDirectionZ[nCscSeg] = globalDirection.z();
      cscSegT[nCscSeg] = cscSegment.time();
      cscSegChi2[nCscSeg] = cscSegment.chi2();
      cscSegStation[nCscSeg] = endcap*CSCDetId::station(id);
      cscSegChamber[nCscSeg] = endcap* (CSCDetId::station(id)*10 + CSCDetId::ring(id));
      //https://github.com/cms-sw/cmssw/blob/master/DataFormats/MuonDetId/interface/CSCDetId.h#L173-L177
      // chamber 14 is also in ME11
      if (CSCDetId::ring(id) == 4) cscRechitsChamber[ncscRechits] = endcap * (CSCDetId::station(id)*10 + 1);

      // cscSegRing[nCscSeg] = CSCDetId::ring(id);


      // for dbscan
      Point p;
      p.phi = globalPosition.phi();
      p.eta = globalPosition.eta();
      p.x = globalPosition.x();
      p.y = globalPosition.y();
      p.z = globalPosition.z();
      p.t = cscSegment.time();
      p.dirX = globalDirection.x();
      p.dirY = globalDirection.y();
      p.dirZ = globalDirection.z();
      p.station = cscSegStation[nCscSeg];
      // p.ring = CSCDetId::ring(id);
	    p.chamber = cscSegChamber[nCscSeg];
      p.clusterID = UNCLASSIFIED;
      points.push_back(p);
      nCscSeg++;
    }
  }


    //Do DBSCAN Clustering
    int min_point_seg = 10;  //minimum number of segments to call it a cluster
    float epsilon_seg = 0.2; //cluster radius parameter
    DBSCAN ds_seg(min_point_seg, epsilon_seg, points);
    ds_seg.run();
    ds_seg.result();
    ds_seg.clusterMoments();
    // ds_seg.vertexing();
    ds_seg.sort_clusters();

    //Save cluster information
    for ( auto &tmp : ds_seg.CscCluster ) {
      cscSegClusterX[nCscSegClusters] =tmp.x;
      cscSegClusterY[nCscSegClusters] =tmp.y;
      cscSegClusterZ[nCscSegClusters] =tmp.z;
      cscSegClusterTime[nCscSegClusters] = tmp.t;
      cscSegClusterEta[nCscSegClusters] =tmp.eta;
      cscSegClusterPhi[nCscSegClusters] = tmp.phi;
      cscSegClusterMajorAxis[nCscSegClusters] =tmp.MajorAxis;
      cscSegClusterMinorAxis[nCscSegClusters] =tmp.MinorAxis;
      cscSegClusterXSpread[nCscSegClusters] =tmp.XSpread;
      cscSegClusterYSpread[nCscSegClusters] =tmp.YSpread;
      cscSegClusterZSpread[nCscSegClusters] =tmp.ZSpread;
      cscSegClusterEtaPhiSpread[nCscSegClusters] =tmp.EtaPhiSpread;
      cscSegClusterEtaSpread[nCscSegClusters] =tmp.EtaSpread;
      cscSegClusterPhiSpread[nCscSegClusters] = tmp.PhiSpread;
      cscSegClusterTimeSpread[nCscSegClusters] = tmp.TSpread;
      cscSegClusterSize[nCscSegClusters] = tmp.nCscSegments;

      cscSegClusterNSegmentChamberPlus11[nCscSegClusters] = tmp.nCscSegmentChamberPlus11;
      cscSegClusterNSegmentChamberPlus12[nCscSegClusters] = tmp.nCscSegmentChamberPlus12;
      cscSegClusterNSegmentChamberPlus13[nCscSegClusters] = tmp.nCscSegmentChamberPlus13;
      cscSegClusterNSegmentChamberPlus21[nCscSegClusters] = tmp.nCscSegmentChamberPlus21;
      cscSegClusterNSegmentChamberPlus22[nCscSegClusters] = tmp.nCscSegmentChamberPlus22;
      cscSegClusterNSegmentChamberPlus31[nCscSegClusters] = tmp.nCscSegmentChamberPlus31;
      cscSegClusterNSegmentChamberPlus32[nCscSegClusters] = tmp.nCscSegmentChamberPlus32;
      cscSegClusterNSegmentChamberPlus41[nCscSegClusters] = tmp.nCscSegmentChamberPlus41;
      cscSegClusterNSegmentChamberPlus42[nCscSegClusters] = tmp.nCscSegmentChamberPlus42;
      cscSegClusterNSegmentChamberMinus11[nCscSegClusters] = tmp.nCscSegmentChamberMinus11;
      cscSegClusterNSegmentChamberMinus12[nCscSegClusters] = tmp.nCscSegmentChamberMinus12;
      cscSegClusterNSegmentChamberMinus13[nCscSegClusters] = tmp.nCscSegmentChamberMinus13;
      cscSegClusterNSegmentChamberMinus21[nCscSegClusters] = tmp.nCscSegmentChamberMinus21;
      cscSegClusterNSegmentChamberMinus22[nCscSegClusters] = tmp.nCscSegmentChamberMinus22;
      cscSegClusterNSegmentChamberMinus31[nCscSegClusters] = tmp.nCscSegmentChamberMinus31;
      cscSegClusterNSegmentChamberMinus32[nCscSegClusters] = tmp.nCscSegmentChamberMinus32;
      cscSegClusterNSegmentChamberMinus41[nCscSegClusters] = tmp.nCscSegmentChamberMinus41;
      cscSegClusterNSegmentChamberMinus42[nCscSegClusters] = tmp.nCscSegmentChamberMinus42;
      cscSegClusterMaxChamber[nCscSegClusters] = tmp.maxChamber;
      cscSegClusterMaxChamberRatio[nCscSegClusters] = 1.0*tmp.maxChamberSegment/tmp.nCscSegments;
      cscSegClusterNChamber[nCscSegClusters] = tmp.nChamber;
      cscSegClusterMaxStation[nCscSegClusters] = tmp.maxStation;
      cscSegClusterMaxStationRatio[nCscSegClusters] = 1.0*tmp.maxStationSegment/tmp.nCscSegments;
      cscSegClusterNStation[nCscSegClusters] = tmp.nStation;

      cscSegClusterMe11Ratio[nCscSegClusters] = tmp.Me11Ratio;
      cscSegClusterMe12Ratio[nCscSegClusters] = tmp.Me12Ratio;
      // cscSegClusterVertexR[nCscSegClusters] = tmp.vertex_r;
      // cscSegClusterVertexZ[nCscSegClusters] = tmp.vertex_z;
      // cscSegClusterVertexChi2[nCscSegClusters] = tmp.vertex_chi2;
      // cscSegClusterVertexDis[nCscSegClusters] = tmp.vertex_dis;
      // cscSegClusterVertexN[nCscSegClusters] = tmp.vertex_n;
      // cscSegClusterVertexN1[nCscSegClusters] = tmp.vertex_n1;
      // cscSegClusterVertexN5[nCscSegClusters] = tmp.vertex_n5;
      // cscSegClusterVertexN15[nCscSegClusters] = tmp.vertex_n15;
      // cscSegClusterVertexN20[nCscSegClusters] = tmp.vertex_n20;
      // cscSegClusterVertexN10[nCscSegClusters] = tmp.vertex_n10;

      //Jet veto/ muon veto
      cscSegClusterJetVetoPt[nCscSegClusters] = 0.0;
      cscSegClusterJetVetoE[nCscSegClusters] = 0.0;
      cscSegClusterCaloJetVeto[nCscSegClusters] = 0.0;
      cscSegClusterMuonVetoPt[nCscSegClusters] = 0.0;
      cscSegClusterMuonVetoE[nCscSegClusters] = 0.0;

      for (const reco::PFJet &j : *jets) {
        //if (j.pt() < 10) continue;
        if (fabs(j.eta())>3.0) continue;
        if (deltaR(tmp.eta, tmp.phi, j.eta(),j.phi()) < 0.4 && j.pt() > cscSegClusterJetVetoPt[nCscSegClusters] ) {
          cscSegClusterJetVetoPt[nCscSegClusters]  = j.pt();
        }
        if (deltaR(tmp.eta, tmp.phi, j.eta(),j.phi()) < 0.4 && j.energy() > cscSegClusterJetVetoE[nCscSegClusters] ) {
          cscSegClusterJetVetoE[nCscSegClusters]  = j.energy();
        }
      }
      for(const pat::Muon &mu : *muons) {
	//if (mu.pt() < 20.0) continue;
        if (fabs(mu.eta()) > 3.0) continue;
        if (deltaR(tmp.eta, tmp.phi, mu.eta(), mu.phi()) < 0.4 && mu.pt() > cscSegClusterMuonVetoPt[nCscSegClusters]) {
          cscSegClusterMuonVetoPt[nCscSegClusters] = mu.pt();
        }
        if (deltaR(tmp.eta, tmp.phi, mu.eta(), mu.phi()) < 0.4 && mu.energy() > cscSegClusterMuonVetoE[nCscSegClusters]) {
          cscSegClusterMuonVetoE[nCscSegClusters] = mu.energy();
        }
      }

      if (!isData) {
	//match to genparticles
	float min_deltaR = 15.;
	int index = 999;
	for(int j = 0; j < nGenParticle; j++)
	  {
	    if (abs(gParticleId[j]) >= 100 && abs(gParticleId[j]) <=350) continue;

	    double current_delta_r = deltaR(cscSegClusterEta[nCscSegClusters], cscSegClusterPhi[nCscSegClusters], gParticleEta[j], gParticlePhi[j]);
	    if (current_delta_r < min_deltaR)
	      {
		min_deltaR = current_delta_r;
		index = j;
	      }
	  }
	if (min_deltaR < 0.4)
	  {
	    cscSegCluster_match_gParticle_minDeltaR[nCscSegClusters] = min_deltaR;
	    cscSegCluster_match_gParticle_index[nCscSegClusters] = index;
	    cscSegCluster_match_gParticle_id[nCscSegClusters] = gParticleId[index];
	  }
      }

      nCscSegClusters++;
    }




  //************************************************************************************************************
  //************************************************************************************************************
  //** CSC RECHITS
  //************************************************************************************************************
  //************************************************************************************************************



  //cout << "Number of rec hits: "<<cscRechits->size()<<endl;
  points.clear();
  for (const CSCRecHit2D cscRechit : *cscRechits) {
    LocalPoint  cscRecHitLocalPosition       = cscRechit.localPosition();
    // LocalError  segmentLocalDirectionError = iDT->localDirectionError();
    CSCDetId cscdetid = cscRechit.cscDetId();
    cscRechitsDetId[ncscRechits] = CSCDetId::rawIdMaker(CSCDetId::endcap(cscdetid), CSCDetId::station(cscdetid), CSCDetId::ring(cscdetid), CSCDetId::chamber(cscdetid), CSCDetId::layer(cscdetid));
    int endcap = CSCDetId::endcap(cscdetid) == 1 ? 1 : -1;
    const CSCChamber* cscchamber = cscG->chamber(cscdetid);
    if (cscchamber) {
      GlobalPoint globalPosition = cscchamber->toGlobal(cscRecHitLocalPosition);
      cscRechitsX[ncscRechits] = globalPosition.x();
      cscRechitsY[ncscRechits] = globalPosition.y();
      cscRechitsZ[ncscRechits] = globalPosition.z();
      cscRechitsPhi[ncscRechits] = globalPosition.phi();
      cscRechitsEta[ncscRechits] = globalPosition.eta();
      cscRechitsE[ncscRechits] = cscRechit.energyDepositedInLayer();//not saved
      cscRechitsTpeak[ncscRechits] = cscRechit.tpeak();
      cscRechitsTwire[ncscRechits] = cscRechit.wireTime();
      cscRechitsQuality[ncscRechits] = cscRechit.quality();
      cscRechitsChamber[ncscRechits] = endcap * (CSCDetId::station(cscdetid)*10 + CSCDetId::ring(cscdetid));
      if (CSCDetId::ring(cscdetid) == 4) cscRechitsChamber[ncscRechits] = endcap * (CSCDetId::station(cscdetid)*10 + 1);
      // cscRechitsRing[ncscRechits] = CSCDetId::ring(cscdetid);

      cscRechitsStation[ncscRechits] = endcap *CSCDetId::station(cscdetid);
      // cscRechitsChannels[ncscRechits] = cscRechit.channels();
      cscRechitsNStrips[ncscRechits] = cscRechit.nStrips();
      cscRechitsHitWire[ncscRechits] = cscRechit.hitWire();
      cscRechitsWGroupsBX[ncscRechits] = cscRechit.wgroupsBX();
      cscRechitsNWireGroups[ncscRechits] = cscRechit.nWireGroups();

      Point p;
	    p.phi = cscRechitsPhi[ncscRechits];
	    p.eta = cscRechitsEta[ncscRechits];
	    p.x = cscRechitsX[ncscRechits];
	    p.y = cscRechitsY[ncscRechits];
	    p.z = cscRechitsZ[ncscRechits];
	    p.t = cscRechitsTpeak[ncscRechits];
	    p.station = cscRechitsStation[ncscRechits];
	    p.chamber = cscRechitsChamber[ncscRechits];
	    p.clusterID = UNCLASSIFIED;
	    points.push_back(p);

      ncscRechits++;
    }
  }
  //Do DBSCAN Clustering
  int min_point = 50;  //minimum number of segments to call it a cluster
  float epsilon = 0.2; //cluster radius parameter
  DBSCAN ds(min_point, epsilon, points);
  ds.run();
  ds.result();
  ds.clusterMoments();
  // ds.vertexing();
  ds.sort_clusters();
  //Save cluster information
  for ( auto &tmp : ds.CscCluster ) {
    cscRechitClusterX[nCscRechitClusters] =tmp.x;
    cscRechitClusterY[nCscRechitClusters] =tmp.y;
    cscRechitClusterZ[nCscRechitClusters] =tmp.z;
    cscRechitClusterTime[nCscRechitClusters] = tmp.t;
    cscRechitClusterEta[nCscRechitClusters] =tmp.eta;
    cscRechitClusterPhi[nCscRechitClusters] = tmp.phi;
    cscRechitClusterMajorAxis[nCscRechitClusters] =tmp.MajorAxis;
    cscRechitClusterMinorAxis[nCscRechitClusters] =tmp.MinorAxis;
    cscRechitClusterXSpread[nCscRechitClusters] =tmp.XSpread;
    cscRechitClusterYSpread[nCscRechitClusters] =tmp.YSpread;
    cscRechitClusterZSpread[nCscRechitClusters] =tmp.ZSpread;
    cscRechitClusterEtaPhiSpread[nCscRechitClusters] =tmp.EtaPhiSpread;
    cscRechitClusterEtaSpread[nCscRechitClusters] =tmp.EtaSpread;
    cscRechitClusterPhiSpread[nCscRechitClusters] = tmp.PhiSpread;
    cscRechitClusterTimeSpread[nCscRechitClusters] = tmp.TSpread;
    cscRechitClusterSize[nCscRechitClusters] = tmp.nCscSegments;

    cscRechitClusterNRechitChamberPlus11[nCscRechitClusters] = tmp.nCscSegmentChamberPlus11;
    cscRechitClusterNRechitChamberPlus12[nCscRechitClusters] = tmp.nCscSegmentChamberPlus12;
    cscRechitClusterNRechitChamberPlus13[nCscRechitClusters] = tmp.nCscSegmentChamberPlus13;
    cscRechitClusterNRechitChamberPlus21[nCscRechitClusters] = tmp.nCscSegmentChamberPlus21;
    cscRechitClusterNRechitChamberPlus22[nCscRechitClusters] = tmp.nCscSegmentChamberPlus22;
    cscRechitClusterNRechitChamberPlus31[nCscRechitClusters] = tmp.nCscSegmentChamberPlus31;
    cscRechitClusterNRechitChamberPlus32[nCscRechitClusters] = tmp.nCscSegmentChamberPlus32;
    cscRechitClusterNRechitChamberPlus41[nCscRechitClusters] = tmp.nCscSegmentChamberPlus41;
    cscRechitClusterNRechitChamberPlus42[nCscRechitClusters] = tmp.nCscSegmentChamberPlus42;
    cscRechitClusterNRechitChamberMinus11[nCscRechitClusters] = tmp.nCscSegmentChamberMinus11;
    cscRechitClusterNRechitChamberMinus12[nCscRechitClusters] = tmp.nCscSegmentChamberMinus12;
    cscRechitClusterNRechitChamberMinus13[nCscRechitClusters] = tmp.nCscSegmentChamberMinus13;
    cscRechitClusterNRechitChamberMinus21[nCscRechitClusters] = tmp.nCscSegmentChamberMinus21;
    cscRechitClusterNRechitChamberMinus22[nCscRechitClusters] = tmp.nCscSegmentChamberMinus22;
    cscRechitClusterNRechitChamberMinus31[nCscRechitClusters] = tmp.nCscSegmentChamberMinus31;
    cscRechitClusterNRechitChamberMinus32[nCscRechitClusters] = tmp.nCscSegmentChamberMinus32;
    cscRechitClusterNRechitChamberMinus41[nCscRechitClusters] = tmp.nCscSegmentChamberMinus41;
    cscRechitClusterNRechitChamberMinus42[nCscRechitClusters] = tmp.nCscSegmentChamberMinus42;
    cscRechitClusterMaxChamber[nCscRechitClusters] = tmp.maxChamber;
    cscRechitClusterMaxChamberRatio[nCscRechitClusters] = 1.0*tmp.maxChamberSegment/tmp.nCscSegments;
    cscRechitClusterNChamber[nCscRechitClusters] = tmp.nChamber;
    cscRechitClusterMaxStation[nCscRechitClusters] = tmp.maxStation;
    cscRechitClusterMaxStationRatio[nCscRechitClusters] = 1.0*tmp.maxStationSegment/tmp.nCscSegments;
    cscRechitClusterNStation[nCscRechitClusters] = tmp.nStation;

    cscRechitClusterMe11Ratio[nCscRechitClusters] = tmp.Me11Ratio;
    cscRechitClusterMe12Ratio[nCscRechitClusters] = tmp.Me12Ratio;
    // cscRechitClusterVertexR[nCscRechitClusters] = tmp.vertex_r;
    // cscRechitClusterVertexZ[nCscRechitClusters] = tmp.vertex_z;
    // cscRechitClusterVertexChi2[nCscRechitClusters] = tmp.vertex_chi2;
    // cscRechitClusterVertexDis[nCscRechitClusters] = tmp.vertex_dis;
    // cscRechitClusterVertexN[nCscRechitClusters] = tmp.vertex_n;
    // cscRechitClusterVertexN1[nCscRechitClusters] = tmp.vertex_n1;
    // cscRechitClusterVertexN5[nCscRechitClusters] = tmp.vertex_n5;
    // cscRechitClusterVertexN15[nCscRechitClusters] = tmp.vertex_n15;
    // cscRechitClusterVertexN20[nCscRechitClusters] = tmp.vertex_n20;
    // cscRechitClusterVertexN10[nCscRechitClusters] = tmp.vertex_n10;

    //Jet veto/ muon veto
    cscRechitClusterJetVetoPt[nCscRechitClusters] = 0.0;
    cscRechitClusterJetVetoE[nCscRechitClusters] = 0.0;
    cscRechitClusterCaloJetVeto[nCscRechitClusters] = 0.0;
    cscRechitClusterMuonVetoPt[nCscRechitClusters] = 0.0;
    cscRechitClusterMuonVetoE[nCscRechitClusters] = 0.0;
    for (const reco::PFJet &j : *jets) {
      //if (j.pt() < 10) continue;
      if (fabs(j.eta())>3.0) continue;
      if (deltaR(tmp.eta, tmp.phi, j.eta(),j.phi()) < 0.4 && j.pt() > cscRechitClusterJetVetoPt[nCscRechitClusters] ) {
        cscRechitClusterJetVetoPt[nCscRechitClusters]  = j.pt();
      }
      if (deltaR(tmp.eta, tmp.phi, j.eta(),j.phi()) < 0.4 && j.energy() > cscRechitClusterJetVetoE[nCscRechitClusters] ) {
        cscRechitClusterJetVetoE[nCscRechitClusters]  = j.energy();
      }
    }
    for(const pat::Muon &mu : *muons) {
    //if (mu.pt() < 20.0) continue;
      if (fabs(mu.eta()) > 3.0) continue;
      if (deltaR(tmp.eta, tmp.phi, mu.eta(), mu.phi()) < 0.4 && mu.pt() > cscRechitClusterMuonVetoPt[nCscRechitClusters]) {
        cscRechitClusterMuonVetoPt[nCscRechitClusters] = mu.pt();
      }
      if (deltaR(tmp.eta, tmp.phi, mu.eta(), mu.phi()) < 0.4 && mu.energy() > cscRechitClusterMuonVetoE[nCscRechitClusters]) {
        cscRechitClusterMuonVetoE[nCscRechitClusters] = mu.energy();
      }
    }

    //match to segment clusters
    float min_deltaR = 15.;
    int index = 999;

    for(int j = 0; j < nCscSegClusters; j++)
    {
      double current_delta_r = deltaR(cscRechitClusterEta[nCscRechitClusters], cscRechitClusterPhi[nCscRechitClusters], cscSegClusterEta[j], cscSegClusterPhi[j]);
      if (current_delta_r < min_deltaR)
      {
        min_deltaR = current_delta_r;
        index = j;
      }
    }
    if (min_deltaR < 0.4)
    {
      cscRechitCluster_match_cscSegCluster_minDeltaR[nCscRechitClusters] = min_deltaR;
      cscRechitCluster_match_cscSegCluster_index[nCscRechitClusters] = index;
    }
    min_deltaR = 15.;
    index = 999;


    if (!isData) {
      //match to genparticles
      for(int j = 0; j < nGenParticle; j++)
	{
	  if (abs(gParticleId[j]) >= 100 && abs(gParticleId[j]) <=350) continue;
	  double current_delta_r = deltaR(cscRechitClusterEta[nCscRechitClusters], cscRechitClusterPhi[nCscRechitClusters], gParticleEta[j], gParticlePhi[j]);
	  if (current_delta_r < min_deltaR)
	    {
	      min_deltaR = current_delta_r;
	      index = j;
	    }
	}
      if (min_deltaR < 0.4)
	{
	  cscRechitCluster_match_gParticle_minDeltaR[nCscRechitClusters] = min_deltaR;
	  cscRechitCluster_match_gParticle_index[nCscRechitClusters] = index;
	  cscRechitCluster_match_gParticle_id[nCscRechitClusters] = gParticleId[index];
	}
    }

    nCscRechitClusters++;
  }



  //************************************************************************************************************
  //************************************************************************************************************
  //** DT SEGMENTS
  //************************************************************************************************************
  //************************************************************************************************************

  // Segments
  points.clear();
  for(DTRecSegment4D dtSegment : *dtSegments){
    const DTRecSegment4D dtSegment_copy = dtSegment;
    const DTChamberRecSegment2D* phiSeg = dtSegment_copy.phiSegment();

    LocalPoint  segmentLocalPosition       = dtSegment.localPosition();
    // LocalVector segmentLocalDirection      = dtSegment.localDirection();
    // LocalError  segmentLocalPositionError  = iDT->localPositionError();
    // LocalError  segmentLocalDirectionError = iDT->localDirectionError();
    DetId geoid = dtSegment.geographicalId();
    DTChamberId dtdetid = DTChamberId(geoid);
    const DTChamber * dtchamber = dtG->chamber(dtdetid);
    if (dtchamber) {
        GlobalPoint globalPosition = dtchamber->toGlobal(segmentLocalPosition);
        // GlobalVector globalDirection = dtchamber->toGlobal(segmentLocalDirection);

        dtSegPhi[nDtSeg] = globalPosition.phi();
        dtSegEta[nDtSeg] = globalPosition.eta();
        dtSegX[nDtSeg] = globalPosition.x();
        dtSegY[nDtSeg] = globalPosition.y();
        dtSegZ[nDtSeg] = globalPosition.z();
        // dtSegDirX[nDtSeg] = globalDirection.x();
        // dtSegDirY[nDtSeg] = globalDirection.y();
        // dtSegDirZ[nDtSeg] = globalDirection.z();
        dtSegStation[nDtSeg] = dtdetid.station();
        dtSegWheel[nDtSeg] = dtdetid.wheel();
        if (phiSeg)
        {
          if(phiSeg->ist0Valid())
            {
              dtSegTime[nDtSeg] = phiSeg->t0();
              dtSegTimeError[nDtSeg] = -1;
            }
        }
        Point p;
        p.phi = dtSegPhi[nDtSeg];
        p.eta = dtSegEta[nDtSeg];
        p.x = dtSegX[nDtSeg];
        p.y = dtSegY[nDtSeg];
        p.z = dtSegZ[nDtSeg];
        p.t = dtSegTime[nDtSeg];
        p.station = dtSegStation[nDtSeg];
        p.chamber = dtSegWheel[nDtSeg];
        p.clusterID = UNCLASSIFIED;
        points.push_back(p);

        nDtSeg++;
      }

    }
    //Do DBSCAN Clustering
    int min_point_dtseg = 10;  //minimum number of segments to call it a cluster
    float epsilon_dtseg = 0.2; //cluster radius parameter
    DBSCAN ds_dtseg(min_point_dtseg, epsilon_dtseg, points);
    ds_dtseg.run();
    ds_dtseg.result();
    ds_dtseg.clusterMoments();
    // ds_dtseg.vertexing();
    ds_dtseg.sort_clusters();

    //Save cluster information
    for ( auto &tmp : ds_dtseg.CscCluster ) {
      dtSegClusterX[nDtSegClusters] =tmp.x;
      dtSegClusterY[nDtSegClusters] =tmp.y;
      dtSegClusterZ[nDtSegClusters] =tmp.z;
      dtSegClusterTime[nDtSegClusters] = tmp.t;
      dtSegClusterEta[nDtSegClusters] =tmp.eta;
      dtSegClusterPhi[nDtSegClusters] = tmp.phi;
      dtSegClusterMajorAxis[nDtSegClusters] =tmp.MajorAxis;
      dtSegClusterMinorAxis[nDtSegClusters] =tmp.MinorAxis;
      dtSegClusterXSpread[nDtSegClusters] =tmp.XSpread;
      dtSegClusterYSpread[nDtSegClusters] =tmp.YSpread;
      dtSegClusterZSpread[nDtSegClusters] =tmp.ZSpread;
      dtSegClusterEtaPhiSpread[nDtSegClusters] =tmp.EtaPhiSpread;
      dtSegClusterEtaSpread[nDtSegClusters] =tmp.EtaSpread;
      dtSegClusterPhiSpread[nDtSegClusters] = tmp.PhiSpread;
      dtSegClusterTimeSpread[nDtSegClusters] = tmp.TSpread;
      dtSegClusterSize[nDtSegClusters] = tmp.nCscSegments;

      dtSegClusterNSegmentStation1[nDtSegClusters] = tmp.nDtSegmentStation1;
      dtSegClusterNSegmentStation2[nDtSegClusters] = tmp.nDtSegmentStation2;
      dtSegClusterNSegmentStation3[nDtSegClusters] = tmp.nDtSegmentStation3;
      dtSegClusterNSegmentStation4[nDtSegClusters] = tmp.nDtSegmentStation4;

      dtSegClusterMaxChamber[nDtSegClusters] = tmp.maxChamber;
      dtSegClusterMaxChamberRatio[nDtSegClusters] = 1.0*tmp.maxChamberSegment/tmp.nCscSegments;
      dtSegClusterNChamber[nDtSegClusters] = tmp.nChamber;
      dtSegClusterMaxStation[nDtSegClusters] = tmp.maxStation;
      dtSegClusterMaxStationRatio[nDtSegClusters] = 1.0*tmp.maxStationSegment/tmp.nCscSegments;
      dtSegClusterNStation[nDtSegClusters] = tmp.nStation;

      // dtSegClusterMe11Ratio[nDtSegClusters] = tmp.Me11Ratio;
      // dtSegClusterMe12Ratio[nDtSegClusters] = tmp.Me12Ratio;
      // dtSegClusterVertexR[nDtSegClusters] = tmp.vertex_r;
      // dtSegClusterVertexZ[nDtSegClusters] = tmp.vertex_z;
      // dtSegClusterVertexChi2[nDtSegClusters] = tmp.vertex_chi2;
      // dtSegClusterVertexDis[nDtSegClusters] = tmp.vertex_dis;
      // dtSegClusterVertexN[nDtSegClusters] = tmp.vertex_n;
      // dtSegClusterVertexN1[nDtSegClusters] = tmp.vertex_n1;
      // dtSegClusterVertexN5[nDtSegClusters] = tmp.vertex_n5;
      // dtSegClusterVertexN15[nDtSegClusters] = tmp.vertex_n15;
      // dtSegClusterVertexN20[nDtSegClusters] = tmp.vertex_n20;
      // dtSegClusterVertexN10[nDtSegClusters] = tmp.vertex_n10;

      //Jet veto/ muon veto
      dtSegClusterJetVetoPt[nDtSegClusters] = 0.0;
      dtSegClusterJetVetoE[nDtSegClusters] = 0.0;
      dtSegClusterCaloJetVeto[nDtSegClusters] = 0.0;
      dtSegClusterMuonVetoPt[nDtSegClusters] = 0.0;
      dtSegClusterMuonVetoE[nDtSegClusters] = 0.0;

      for (const reco::PFJet &j : *jets) {
        //if (j.pt() < 10) continue;
        if (fabs(j.eta())>3.0) continue;
        if (deltaR(tmp.eta, tmp.phi, j.eta(),j.phi()) < 0.4 && j.pt() > dtSegClusterJetVetoPt[nDtSegClusters] ) {
          dtSegClusterJetVetoPt[nDtSegClusters]  = j.pt();
        }
        if (deltaR(tmp.eta, tmp.phi, j.eta(),j.phi()) < 0.4 && j.energy() > dtSegClusterJetVetoE[nDtSegClusters] ) {
          dtSegClusterJetVetoE[nDtSegClusters]  = j.energy();
        }
      }
      for(const pat::Muon &mu : *muons) {
  //if (mu.pt() < 20.0) continue;
        if (fabs(mu.eta()) > 3.0) continue;
        if (deltaR(tmp.eta, tmp.phi, mu.eta(), mu.phi()) < 0.4 && mu.pt() > dtSegClusterMuonVetoPt[nDtSegClusters]) {
          dtSegClusterMuonVetoPt[nDtSegClusters] = mu.pt();
        }
        if (deltaR(tmp.eta, tmp.phi, mu.eta(), mu.phi()) < 0.4 && mu.energy() > dtSegClusterMuonVetoE[nDtSegClusters]) {
          dtSegClusterMuonVetoE[nDtSegClusters] = mu.energy();
        }
      }

      if (!isData) {
	//match to genparticles
	float min_deltaR = 15.;
	int index = 999;
	for(int j = 0; j < nGenParticle; j++)
	  {
	    if (abs(gParticleId[j]) >= 100 && abs(gParticleId[j]) <=350) continue;

	    double current_delta_r = deltaR(dtSegClusterEta[nDtSegClusters], dtSegClusterPhi[nDtSegClusters], gParticleEta[j], gParticlePhi[j]);
	    if (current_delta_r < min_deltaR)
	      {
		min_deltaR = current_delta_r;
		index = j;
	      }
	  }
	if (min_deltaR < 0.4)
	  {
	    dtSegCluster_match_gParticle_minDeltaR[nDtSegClusters] = min_deltaR;
	    dtSegCluster_match_gParticle_index[nDtSegClusters] = index;
	    dtSegCluster_match_gParticle_id[nDtSegClusters] = gParticleId[index];
	  }
      }

      nDtSegClusters++;
    }



  //************************************************************************************************************
  //************************************************************************************************************
  //** DT RECHITS
  //************************************************************************************************************
  //************************************************************************************************************


  //cout<<"number of dt rechits: " <<dtRechits->size()<<endl;
  points.clear();
  for(DTRecHit1DPair dtRechit: *dtRechits){
    LocalPoint  localPosition       = dtRechit.localPosition();
    DetId geoid = dtRechit.geographicalId();
    DTChamberId dtdetid = DTChamberId(geoid);
    const DTChamber * dtchamber = dtG->chamber(dtdetid);
    if (dtchamber) {
      GlobalPoint globalPosition = dtchamber->toGlobal(localPosition);
      dtRechitPhi[nDtRechits] = globalPosition.phi();
      dtRechitEta[nDtRechits] = globalPosition.eta();
      dtRechitX[nDtRechits] = globalPosition.x();
      dtRechitY[nDtRechits] = globalPosition.y();
      dtRechitZ[nDtRechits] = globalPosition.z();
      dtRechitTime[nDtRechits] = dtRechit.digiTime();
      dtRechitStation[nDtRechits] = dtdetid.station();
      dtRechitWheel[nDtRechits] = dtdetid.wheel();
      Point p;
      p.phi = dtRechitPhi[nDtRechits];
      p.eta = dtRechitEta[nDtRechits];
      p.x = dtRechitX[nDtRechits];
      p.y = dtRechitY[nDtRechits];
      p.z = dtRechitZ[nDtRechits];
      p.t = dtRechitTime[nDtRechits];
      p.station = dtRechitStation[nDtRechits];
      p.chamber = dtRechitWheel[nDtRechits];
      p.clusterID = UNCLASSIFIED;
      points.push_back(p);

      nDtRechits++;
    }
  }
  //Do DBSCAN Clustering
  int min_point_dtrechit = 50;  //minimum number of segments to call it a cluster
  float epsilon_dtrechit = 0.2; //cluster radius parameter
  DBSCAN ds_dtRechit(min_point_dtrechit, epsilon_dtrechit, points);
  ds_dtRechit.run();
  ds_dtRechit.result();
  ds_dtRechit.clusterMoments();
  // ds_dtseg.vertexing();
  ds_dtRechit.sort_clusters();

  //Save cluster information
  for ( auto &tmp : ds_dtRechit.CscCluster ) {
    dtRechitClusterX[nDtRechitClusters] =tmp.x;
    dtRechitClusterY[nDtRechitClusters] =tmp.y;
    dtRechitClusterZ[nDtRechitClusters] =tmp.z;
    dtRechitClusterTime[nDtRechitClusters] = tmp.t;
    dtRechitClusterEta[nDtRechitClusters] =tmp.eta;
    dtRechitClusterPhi[nDtRechitClusters] = tmp.phi;
    dtRechitClusterMajorAxis[nDtRechitClusters] =tmp.MajorAxis;
    dtRechitClusterMinorAxis[nDtRechitClusters] =tmp.MinorAxis;
    dtRechitClusterXSpread[nDtRechitClusters] =tmp.XSpread;
    dtRechitClusterYSpread[nDtRechitClusters] =tmp.YSpread;
    dtRechitClusterZSpread[nDtRechitClusters] =tmp.ZSpread;
    dtRechitClusterEtaPhiSpread[nDtRechitClusters] =tmp.EtaPhiSpread;
    dtRechitClusterEtaSpread[nDtRechitClusters] =tmp.EtaSpread;
    dtRechitClusterPhiSpread[nDtRechitClusters] = tmp.PhiSpread;
    dtRechitClusterTimeSpread[nDtRechitClusters] = tmp.TSpread;
    dtRechitClusterSize[nDtRechitClusters] = tmp.nCscSegments;

    dtRechitClusterNSegmentStation1[nDtRechitClusters] = tmp.nDtSegmentStation1;
    dtRechitClusterNSegmentStation2[nDtRechitClusters] = tmp.nDtSegmentStation2;
    dtRechitClusterNSegmentStation3[nDtRechitClusters] = tmp.nDtSegmentStation3;
    dtRechitClusterNSegmentStation4[nDtRechitClusters] = tmp.nDtSegmentStation4;

    dtRechitClusterMaxChamber[nDtRechitClusters] = tmp.maxChamber;
    dtRechitClusterMaxChamberRatio[nDtRechitClusters] = 1.0*tmp.maxChamberSegment/tmp.nCscSegments;
    dtRechitClusterNChamber[nDtRechitClusters] = tmp.nChamber;
    dtRechitClusterMaxStation[nDtRechitClusters] = tmp.maxStation;
    dtRechitClusterMaxStationRatio[nDtRechitClusters] = 1.0*tmp.maxStationSegment/tmp.nCscSegments;
    dtRechitClusterNStation[nDtRechitClusters] = tmp.nStation;


    //Jet veto/ muon veto
    dtRechitClusterJetVetoPt[nDtRechitClusters] = 0.0;
    dtRechitClusterJetVetoE[nDtRechitClusters] = 0.0;
    dtRechitClusterCaloJetVeto[nDtRechitClusters] = 0.0;
    dtRechitClusterMuonVetoPt[nDtRechitClusters] = 0.0;
    dtRechitClusterMuonVetoE[nDtRechitClusters] = 0.0;

    for (const reco::PFJet &j : *jets) {
      //if (j.pt() < 10) continue;
      if (fabs(j.eta())>3.0) continue;
      if (deltaR(tmp.eta, tmp.phi, j.eta(),j.phi()) < 0.4 && j.pt() > dtRechitClusterJetVetoPt[nDtRechitClusters] ) {
        dtRechitClusterJetVetoPt[nDtRechitClusters]  = j.pt();
      }
      if (deltaR(tmp.eta, tmp.phi, j.eta(),j.phi()) < 0.4 && j.energy() > dtRechitClusterJetVetoE[nDtRechitClusters] ) {
        dtRechitClusterJetVetoE[nDtRechitClusters]  = j.energy();
      }
    }
    for(const pat::Muon &mu : *muons) {
//if (mu.pt() < 20.0) continue;
      if (fabs(mu.eta()) > 3.0) continue;
      if (deltaR(tmp.eta, tmp.phi, mu.eta(), mu.phi()) < 0.4 && mu.pt() > dtRechitClusterMuonVetoPt[nDtRechitClusters]) {
        dtRechitClusterMuonVetoPt[nDtRechitClusters] = mu.pt();
      }
      if (deltaR(tmp.eta, tmp.phi, mu.eta(), mu.phi()) < 0.4 && mu.energy() > dtRechitClusterMuonVetoE[nDtRechitClusters]) {
        dtRechitClusterMuonVetoE[nDtRechitClusters] = mu.energy();
      }
    }

    if (!isData) {
      //match to genparticles
      float min_deltaR = 15.;
      int index = 999;
      for(int j = 0; j < nGenParticle; j++)
	{
	  if (abs(gParticleId[j]) >= 100 && abs(gParticleId[j]) <=350) continue;

	  double current_delta_r = deltaR(dtRechitClusterEta[nDtRechitClusters], dtRechitClusterPhi[nDtRechitClusters], gParticleEta[j], gParticlePhi[j]);
	  if (current_delta_r < min_deltaR)
	    {
	      min_deltaR = current_delta_r;
	      index = j;
	    }
	}
      if (min_deltaR < 0.4)
	{
	  dtRechitCluster_match_gParticle_minDeltaR[nDtRechitClusters] = min_deltaR;
	  dtRechitCluster_match_gParticle_index[nDtRechitClusters] = index;
	  dtRechitCluster_match_gParticle_id[nDtRechitClusters] = gParticleId[index];
	}
    }

    nDtRechitClusters++;
  }

  for (const RPCRecHit rpcRecHit : *rpcRecHits){
    LocalPoint  rpcRecHitLocalPosition       = rpcRecHit.localPosition();
    // LocalError  segmentLocalDirectionError = iDT->localDirectionError();
    DetId geoid = rpcRecHit.geographicalId();
    RPCDetId rpcdetid = RPCDetId(geoid);
    const RPCChamber * rpcchamber = rpcG->chamber(rpcdetid);
    if (rpcchamber) {
      GlobalPoint globalPosition = rpcchamber->toGlobal(rpcRecHitLocalPosition);
      rpcX[nRpc] = globalPosition.x();
      rpcY[nRpc] = globalPosition.y();
      rpcZ[nRpc] = globalPosition.z();
      rpcPhi[nRpc] = globalPosition.phi();
      rpcEta[nRpc] = globalPosition.eta();
      rpcT[nRpc] = rpcRecHit.time();
      rpcBx[nRpc] = rpcRecHit.BunchX();
      rpcTError[nRpc] = rpcRecHit.timeError();
      nRpc++;
    }
  }
  
  return true;
}
bool displacedJetMuon_ntupler::fillGenParticles(){
  std::vector<const reco::Candidate*> prunedV;//Allows easier comparison for mother finding
  //Fills selected gen particles
  const double pt_cut = 0.0;
  // int llp_id = 6000113;

  vector<int> llpIDs;
  llpIDs.push_back(9000006);
  llpIDs.push_back(9000007);
  llpIDs.push_back(1023);
  llpIDs.push_back(1000023);
  llpIDs.push_back(1000025);

  for(size_t i=0; i<genParticles->size();i++)
  {
    if(
       (abs((*genParticles)[i].pdgId()) >= 1 && abs((*genParticles)[i].pdgId()) <= 6 && ( (*genParticles)[i].status() < 30 ))
       || (abs((*genParticles)[i].pdgId()) >= 11 && abs((*genParticles)[i].pdgId()) <= 16)
       || (abs((*genParticles)[i].pdgId()) == 21 && (*genParticles)[i].status() < 30)
       || (abs((*genParticles)[i].pdgId()) == 22 && (*genParticles)[i].pt() > 10.0 )
       || (abs((*genParticles)[i].pdgId()) >= 23 && abs((*genParticles)[i].pdgId()) <= 25)
       || (abs((*genParticles)[i].pdgId()) >= 32 && abs((*genParticles)[i].pdgId()) <= 42)
       || (abs((*genParticles)[i].pdgId()) >= 100 && abs((*genParticles)[i].pdgId()) <= 350)
       || (abs((*genParticles)[i].pdgId()) == 1023)
       || (abs((*genParticles)[i].pdgId()) >= 1000001 && abs((*genParticles)[i].pdgId()) <= 1000039)
       || (abs((*genParticles)[i].pdgId()) == 9000006 || abs((*genParticles)[i].pdgId()) == 9000007)
	)
       {
         if ((*genParticles)[i].pt()>pt_cut){
           prunedV.push_back(&(*genParticles)[i]);
         }
       }
  }
  //Total number of gen particles
  nGenParticle = prunedV.size();
  bool _found_first_llp = false;
  //Look for mother particle and Fill gen variables
  for(unsigned int i = 0; i < prunedV.size(); i++) {
    gParticleId[i] = prunedV[i]->pdgId();
    gParticleStatus[i] = prunedV[i]->status();
    gParticleE[i] = prunedV[i]->energy();
    gParticlePt[i] = prunedV[i]->pt();
    gParticlePx[i] = prunedV[i]->px();
    gParticlePy[i] = prunedV[i]->py();
    gParticlePz[i] = prunedV[i]->pz();
    gParticleEta[i] = prunedV[i]->eta();
    gParticlePhi[i] = prunedV[i]->phi();
    gParticleProdVertexX[i] = prunedV[i]->vx();
    gParticleProdVertexY[i] = prunedV[i]->vy();
    gParticleProdVertexZ[i] = prunedV[i]->vz();
    gParticleMotherId[i] = 0;
    gParticleMotherIndex[i] = -1;
   
    if(prunedV[i]->numberOfMothers() > 0) {
      //find the ID of the first mother that has a different ID than the particle itself
      const reco::Candidate* firstMotherWithDifferentID = findFirstMotherWithDifferentID(prunedV[i]);
      if (firstMotherWithDifferentID)
      {
        gParticleMotherId[i] = firstMotherWithDifferentID->pdgId();
        //gParticleDecayVertexX[i] = firstMotherWithDifferentID->vx();
        //gParticleDecayVertexY[i] = firstMotherWithDifferentID->vy();
        //gParticleDecayVertexZ[i] = firstMotherWithDifferentID->vz();
      }

      //find the mother and keep going up the mother chain if the ID's are the same
      const reco::Candidate* originalMotherWithSameID = findOriginalMotherWithSameID(prunedV[i]);
      for(unsigned int j = 0; j < prunedV.size(); j++) {
        if(prunedV[j] == originalMotherWithSameID) {
          gParticleMotherIndex[i] = j;
          break;
        }
      }
    } else {
      gParticleMotherIndex[i] = -1;
    }

    //---------------------------------------
    //Find LLPs production and decay vertices
    //---------------------------------------
    if (enableGenLLPInfo_) {

      //match with one of the entries in the llpIDs List
      bool matchedLLP = false;
      int matchedLLPID = 0;
      if (gParticleStatus[i] == 22) {
	for (uint d=0 ; d < llpIDs.size() ; d++) {
	  if ( abs(gParticleId[i]) == llpIDs[d] ) {
	    matchedLLPID = gParticleId[i];
	    matchedLLP = true;	    
	  }
	}
      }

      if ( matchedLLP ) {
        if (!_found_first_llp) {	 
          gLLP_prod_vertex_x[0] = prunedV[i]->vx();
          gLLP_prod_vertex_y[0] = prunedV[i]->vy();
          gLLP_prod_vertex_z[0] = prunedV[i]->vz();
        }
        else {    	  
          gLLP_prod_vertex_x[1] = prunedV[i]->vx();
          gLLP_prod_vertex_y[1] = prunedV[i]->vy();
          gLLP_prod_vertex_z[1] = prunedV[i]->vz();
        }

        const reco::Candidate *dau = 0;
        bool foundDaughter = false;
        bool noDaughter = false;
        const reco::Candidate *tmpParticle = prunedV[i];

        while (!foundDaughter && !noDaughter) {
          if (tmpParticle->numberOfDaughters() > 0) {
            dau = tmpParticle->daughter(0);
            if (dau && (dau->pdgId() != matchedLLPID)) {
              foundDaughter = true;
            } else {
              tmpParticle = dau;
            }
          } else {
            noDaughter = true;
          }
        }

        if (foundDaughter) {
	  gParticleDecayVertexX[i] = dau->vx();
	  gParticleDecayVertexY[i] = dau->vy();
	  gParticleDecayVertexZ[i] = dau->vz();

          if (!_found_first_llp) {
	    _found_first_llp = true;
	    gLLP_decay_vertex_x[0] = dau->vx();
            gLLP_decay_vertex_y[0] = dau->vy();
            gLLP_decay_vertex_z[0] = dau->vz();
            gLLP_pt[0] = sqrt(gParticlePx[i]*gParticlePx[i]+gParticlePy[i]*gParticlePy[i]);
            gLLP_e[0] = gParticleE[i];
            gLLP_eta[0] = gParticleEta[i];
            gLLP_phi[0] = gParticlePhi[i];
            gLLP_beta[0] = sqrt(gParticlePx[i]*gParticlePx[i]+gParticlePy[i]*gParticlePy[i]+gParticlePz[i]*gParticlePz[i])/gParticleE[i];
            gLLP_travel_time[0] = sqrt(pow(gLLP_decay_vertex_x[0]-gLLP_prod_vertex_x[0],2)
                                    +pow(gLLP_decay_vertex_y[0]-gLLP_prod_vertex_y[0],2)
                                    +pow(gLLP_decay_vertex_z[0]-gLLP_prod_vertex_z[0],2))/(30. * gLLP_beta[0]);//1/30 is to convert cm to ns

            double radius = sqrt( pow(gLLP_decay_vertex_x[0],2) + pow(gLLP_decay_vertex_y[0],2) );
	    if (abs(gLLP_eta[0]) < 2.4 && abs(gLLP_eta[0]) > 0.9
               && abs(gLLP_decay_vertex_z[0])<1100 && abs(gLLP_decay_vertex_z[0])>568
               && radius < 695.5) gLLP_csc[0] = true;
           if (radius < 740 && radius > 400 && abs(gLLP_decay_vertex_z[0])< 650 ) gLLP_dt[0] = true;
            double ecal_radius = 129.0;
            // double hcal_radius = 179.0;
            double EB_z = 268.36447217; // 129*sinh(1.479)
            double EE_z = 298.5; //where Ecal Endcap starts in z direction

            for (unsigned int id = 0; id < tmpParticle->numberOfDaughters(); id++ ) {
	      //std::cout << "====================" << std::endl;
	      //std::cout << " -> "<< tmpParticle->daughter(id)->pdgId() << std::endl;
              if( id > 1 ) break; //only allow 2 daughters
              TLorentzVector tmp;
              tmp.SetPxPyPzE(tmpParticle->daughter(id)->px(), tmpParticle->daughter(id)->py(), tmpParticle->daughter(id)->pz(), tmpParticle->daughter(id)->energy());
              if(tmp.Pt()<pt_cut) continue;
	      gLLP_daughter_id[id] = tmpParticle->daughter(id)->pdgId();
	      gLLP_daughter_pt[id] = tmp.Pt();
              gLLP_daughter_eta[id] = tmp.Eta();
              gLLP_daughter_phi[id] = tmp.Phi();
              gLLP_daughter_e[id]  = tmp.E();
	      gLLP_daughter_mass[id]  = tmp.M();

              // double gLLP_daughter_travel_time_hcal= (1./30.)*(hcal_radius-radius)/(tmp.Pt()/tmp.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns

              gLLP_daughter_travel_time[id] = (1./30.)*(ecal_radius-radius)/(tmp.Pt()/tmp.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns
            //Calculate dt from generation point to ECAL face
              double x_ecal = gLLP_decay_vertex_x[0] + 30. * (tmp.Px()/tmp.E())*gLLP_daughter_travel_time[id];
              double y_ecal = gLLP_decay_vertex_y[0] + 30. * (tmp.Py()/tmp.E())*gLLP_daughter_travel_time[id];
              double z_ecal = gLLP_decay_vertex_z[0] + 30. * (tmp.Pz()/tmp.E())*gLLP_daughter_travel_time[id];
              // double x_hcal = gLLP_decay_vertex_x[0] + 30. * (tmp.Px()/tmp.E())*gLLP_daughter_travel_time_hcal;
              // double y_hcal = gLLP_decay_vertex_y[0] + 30. * (tmp.Py()/tmp.E())*gLLP_daughter_travel_time_hcal;
              // double z_hcal = gLLP_decay_vertex_z[0] + 30. * (tmp.Pz()/tmp.E())*gLLP_daughter_travel_time_hcal;

              if( fabs(z_ecal) < EB_z && radius <= ecal_radius &&  fabs(gLLP_decay_vertex_z[0]) < EE_z) {
                photon_travel_time[id] = (1./30) * sqrt(pow(ecal_radius,2)+pow(z_ecal,2));
                photon_travel_time_pv[id] = (1./30) * sqrt(pow(x_ecal-genVertexX,2) + pow(y_ecal-genVertexY,2) + pow(z_ecal-genVertexZ,2));
                gen_time_pv[id] =  gLLP_travel_time[0] + gLLP_daughter_travel_time[id] - photon_travel_time_pv[id] + genVertexT;
                gen_time[id] = gLLP_travel_time[0] + gLLP_daughter_travel_time[id] - photon_travel_time[id] + genVertexT;

              } else {
                gLLP_daughter_travel_time[id] = -666;
                gen_time_pv[id] = -666.;
                gen_time[id] = -666.;
                photon_travel_time[id] = -666.;
                photon_travel_time_pv[id] = -666.;
              }

          // Correction of eta and phi based on ecal points
              double phi = atan((y_ecal-genVertexY)/(x_ecal-genVertexX));
              if  (x_ecal < 0.0){
                phi = TMath::Pi() + phi;
              }
              phi = deltaPhi(phi,0.0);
              double theta = atan(sqrt(pow(x_ecal-genVertexX,2)+pow(y_ecal-genVertexY,2))/abs(z_ecal-genVertexZ));
              double eta = -1.0*TMath::Sign(1.0, z_ecal-genVertexZ)*log(tan(theta/2));
              gLLP_daughter_eta_ecalcorr[id] = eta;
              gLLP_daughter_phi_ecalcorr[id] = phi;
             
                       
	      //*********************************
	      //grandaughters
	      //*********************************
	      if(tmpParticle->daughter(id)->numberOfDaughters() > 0) {
		const reco::Candidate *tmpDauParticle = tmpParticle->daughter(id);
		for (unsigned int index = 0; index < tmpDauParticle->numberOfDaughters(); index++ ) {
		  if( index > 1 ) break; //only allow two granddaughters
					   //std::cout <<"index " << index << " -> "<< tmpDauParticle->daughter(index)->pdgId() << std::endl;
		  TLorentzVector tmpdau;
		  tmpdau.SetPxPyPzE(tmpDauParticle->daughter(index)->px(), tmpDauParticle->daughter(index)->py(), tmpDauParticle->daughter(index)->pz(), tmpDauParticle->daughter(index)->energy());
		  if(tmpdau.Pt()<pt_cut) continue;
		  gLLP_grandaughter_id[index] = tmpDauParticle->daughter(index)->pdgId();
		  gLLP_grandaughter_pt[index] = tmpdau.Pt();
		  gLLP_grandaughter_eta[index] = tmpdau.Eta();
		  gLLP_grandaughter_phi[index] = tmpdau.Phi();
		  gLLP_grandaughter_e[index]  = tmpdau.E();
		  gLLP_grandaughter_mass[index]  = tmpdau.M();


		  // Correction of eta and phi based on ecal points
		  double phi = atan((y_ecal-genVertexY)/(x_ecal-genVertexX));
		  if  (x_ecal < 0.0) {
		    phi = TMath::Pi() + phi;
		  }
		  phi = deltaPhi(phi,0.0);
		  double theta = atan(sqrt(pow(x_ecal-genVertexX,2)+pow(y_ecal-genVertexY,2))/abs(z_ecal-genVertexZ));
		  double eta = -1.0*TMath::Sign(1.0, z_ecal-genVertexZ)*log(tan(theta/2));
		  gLLP_grandaughter_eta_ecalcorr[index] = eta;
		  gLLP_grandaughter_phi_ecalcorr[index] = phi;
             		
		}
	      }//loop over all gLLP granddaughters 	 

            }//loop over all gLLP daughters

          } //end if first llp
          else {
            gLLP_decay_vertex_x[1] = dau->vx();
            gLLP_decay_vertex_y[1] = dau->vy();
            gLLP_decay_vertex_z[1] = dau->vz();
            gLLP_pt[1] = sqrt(gParticlePx[i]*gParticlePx[i]+gParticlePy[i]*gParticlePy[i]);
            gLLP_e[1] = gParticleE[i];
            gLLP_eta[1] = gParticleEta[i];
            gLLP_phi[1] = gParticlePhi[i];
            gLLP_beta[1] = sqrt(gParticlePx[i]*gParticlePx[i]+gParticlePy[i]*gParticlePy[i]+gParticlePz[i]*gParticlePz[i])/gParticleE[i];
            gLLP_travel_time[1] = sqrt(pow(gLLP_decay_vertex_x[1]-gLLP_prod_vertex_x[1],2)
                    +pow(gLLP_decay_vertex_y[1]-gLLP_prod_vertex_y[1],2)
                    +pow(gLLP_decay_vertex_z[1]-gLLP_prod_vertex_z[1],2))/(30. * gLLP_beta[1]);//1/30 is to convert cm to ns
            double radius = sqrt( pow(gLLP_decay_vertex_x[1],2) + pow(gLLP_decay_vertex_y[1],2) );
            if (abs(gLLP_eta[1]) < 2.4 && abs(gLLP_eta[1]) > 0.9
		&& abs(gLLP_decay_vertex_z[1])<1100 && abs(gLLP_decay_vertex_z[1])>568
		&& radius < 695.5) gLLP_csc[1] = true;
	    if (radius < 740 && radius > 400 && abs(gLLP_decay_vertex_z[1])< 650 ) gLLP_dt[1] = true;

            double ecal_radius = 129.0;
            // double hcal_radius = 179.0;
            double EB_z = 268.36447217; // 129*sinh(1.479)
            double EE_z = 298.5; //where Ecal Endcap starts in z direction
            //--------------------------------------------------
            //Second two LLP daughters belong to LLP->pdgID()=36
            //--------------------------------------------------
            for (unsigned int id = 0; id < tmpParticle->numberOfDaughters(); id++ )
            {
              //std::cout << " -> "<< tmpParticle->daughter(id)->pdgId() << std::endl;
              if( id > 1 ) break;
              TLorentzVector tmp;
              tmp.SetPxPyPzE(tmpParticle->daughter(id)->px(), tmpParticle->daughter(id)->py(), tmpParticle->daughter(id)->pz(), tmpParticle->daughter(id)->energy());
              if(tmp.Pt()<pt_cut) continue;
	      gLLP_daughter_id[id+2] = tmpParticle->daughter(id)->pdgId();
	      gLLP_daughter_pt[id+2] = tmp.Pt();
              gLLP_daughter_eta[id+2] = tmp.Eta();
              gLLP_daughter_phi[id+2] = tmp.Phi();
              gLLP_daughter_e[id+2]  = tmp.E();
	      gLLP_daughter_mass[id+2]  = tmp.M();
	      //gLLP_daughter_travel_time[id+2] = (1./30.)*(ecal_radius-radius)/(tmp.Pt()/tmp.E()) - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns
              // double gLLP_daughter_travel_time_hcal = (1./30.)*(hcal_radius-radius)/(tmp.Pt()/tmp.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns
              gLLP_daughter_travel_time[id+2] = (1./30.)*(ecal_radius-radius)/(tmp.Pt()/tmp.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns

              //Calculate dt from generation point to ECAL face
              double x_ecal = gLLP_decay_vertex_x[1] + 30. * (tmp.Px()/tmp.E())*gLLP_daughter_travel_time[id+2];
              double y_ecal = gLLP_decay_vertex_y[1] + 30. * (tmp.Py()/tmp.E())*gLLP_daughter_travel_time[id+2];
              double z_ecal = gLLP_decay_vertex_z[1] + 30. * (tmp.Pz()/tmp.E())*gLLP_daughter_travel_time[id+2];
              //
              // double x_hcal = gLLP_decay_vertex_x[1] + 30. * (tmp.Px()/tmp.E())*gLLP_daughter_travel_time_hcal;
              // double y_hcal = gLLP_decay_vertex_y[1] + 30. * (tmp.Py()/tmp.E())*gLLP_daughter_travel_time_hcal;
              // double z_hcal = gLLP_decay_vertex_z[1] + 30. * (tmp.Pz()/tmp.E())*gLLP_daughter_travel_time_hcal;
              if( fabs(z_ecal) < EB_z && radius <= ecal_radius && fabs(gLLP_decay_vertex_z[1]) < EE_z)
              // if( fabs(z_ecal) < 10 && radius <= 0.1)
              {
                photon_travel_time[id+2] = (1./30) * sqrt(pow(ecal_radius,2)+pow(z_ecal,2));
                photon_travel_time_pv[id+2] = (1./30) * sqrt(pow(x_ecal-genVertexX,2) + pow(y_ecal-genVertexY,2) + pow(z_ecal-genVertexZ,2));
                gen_time_pv[id+2] =  gLLP_travel_time[1] + gLLP_daughter_travel_time[id+2] - photon_travel_time_pv[id+2] + genVertexT;
                gen_time[id+2] = gLLP_travel_time[1] + gLLP_daughter_travel_time[id+2] - photon_travel_time[id+2] + genVertexT;

              } else {
                gLLP_daughter_travel_time[id+2] = -666;
                gen_time_pv[id+2] = -666.;
                gen_time[id+2] = -666.;
                photon_travel_time[id+2] = -666.;
                photon_travel_time_pv[id+2] = -666.;
              }


              //Corrections for angles based on ECAL radius
              double phi = atan((y_ecal-genVertexY)/(x_ecal-genVertexX));
              if  (x_ecal < 0.0){
                phi = TMath::Pi() + phi;
              }
              phi = deltaPhi(phi,0.0);
              double theta = atan(sqrt(pow(x_ecal-genVertexX,2)+pow(y_ecal-genVertexY,2))/abs(z_ecal-genVertexZ));
              double eta = -1.0*TMath::Sign(1.0,z_ecal-genVertexZ)*log(tan(theta/2));
              gLLP_daughter_eta_ecalcorr[id+2] = eta;
              gLLP_daughter_phi_ecalcorr[id+2] = phi;                                       

	      //grandaughters
	      if(tmpParticle->daughter(id)->numberOfDaughters() > 0) {
		const reco::Candidate *tmpDauParticle = tmpParticle->daughter(id);
		for (unsigned int index = 0; index < tmpDauParticle->numberOfDaughters(); index++ ) {
		  if( index > 1 ) break;
		  //std::cout <<"index " << index << " -> "<< tmpDauParticle->daughter(index)->pdgId() << std::endl;
		  TLorentzVector tmpdau;
		  tmpdau.SetPxPyPzE(tmpDauParticle->daughter(index)->px(), tmpDauParticle->daughter(index)->py(), tmpDauParticle->daughter(index)->pz(), tmpDauParticle->daughter(index)->energy());
		  if(tmpdau.Pt()<pt_cut) continue;
		  gLLP_grandaughter_id[index+2] = tmpDauParticle->daughter(index)->pdgId();
		  gLLP_grandaughter_pt[index+2] = tmpdau.Pt();
		  gLLP_grandaughter_eta[index+2] = tmpdau.Eta();
		  gLLP_grandaughter_phi[index+2] = tmpdau.Phi();
		  gLLP_grandaughter_e[index+2]  = tmpdau.E();
		  gLLP_grandaughter_mass[index+2]  = tmpdau.M();




		  // Correction of eta and phi based on ecal points
		  double phi = atan((y_ecal-genVertexY)/(x_ecal-genVertexX));
		  if  (x_ecal < 0.0){
		    phi = TMath::Pi() + phi;
		  }
		  phi = deltaPhi(phi,0.0);
		  double theta = atan(sqrt(pow(x_ecal-genVertexX,2)+pow(y_ecal-genVertexY,2))/abs(z_ecal-genVertexZ));
		  double eta = -1.0*TMath::Sign(1.0, z_ecal-genVertexZ)*log(tan(theta/2));
		  gLLP_grandaughter_eta_ecalcorr[index+2] = eta;
		  gLLP_grandaughter_phi_ecalcorr[index+2] = phi;		   		 
		    
		} //loop of all gLLP granddaughters 	      
	      } //if daughter has granddaughters 	    
            }//for daughters loop
          }//if particle ID = 36
        }//if found daughters
      } //if LLP matched
    } //end if enableGenLLPInfo
  }// for loop of genParticles
  return true;
};


bool displacedJetMuon_ntupler::fillJets(const edm::EventSetup& iSetup)
{

  edm::ESHandle<CaloGeometry> geoHandle;
  iSetup.get<CaloGeometryRecord>().get(geoHandle);
  const CaloSubdetectorGeometry *barrelGeometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
  const CaloSubdetectorGeometry *endcapGeometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);
  const CaloSubdetectorGeometry *hbGeometry = geoHandle->getSubdetectorGeometry(DetId::Hcal, HcalBarrel);
  const CaloSubdetectorGeometry *heGeometry = geoHandle->getSubdetectorGeometry(DetId::Hcal, HcalEndcap);
  const CaloSubdetectorGeometry *hoGeometry = geoHandle->getSubdetectorGeometry(DetId::Hcal, HcalOuter);

  vector<bool> SaveThisEBRechit; SaveThisEBRechit.clear();
  vector<bool> SaveThisEERechit; SaveThisEERechit.clear();
  vector<bool> SaveThisHCALRechit; SaveThisHCALRechit.clear();
  vector<bool> SaveThisHORechit; SaveThisHORechit.clear();
  vector<bool> SaveThisTrack; SaveThisTrack.clear();
  //reset these save flags
  for (uint q=0; q<ebRecHits->size(); q++) {
    SaveThisEBRechit.push_back(false);    
  }
  for (uint q=0; q<eeRecHits->size(); q++) {
    SaveThisEERechit.push_back(false);    
  }
  for (unsigned int iTrack = 0; iTrack < generalTracks->size(); iTrack ++){
    SaveThisTrack.push_back(false);  
  }
  for (uint q=0; q<hcalRecHitsHBHE->size(); q++) {
    SaveThisHCALRechit.push_back(false);    
  }
  for (uint q=0; q<hcalRecHitsHO->size(); q++) {
    SaveThisHORechit.push_back(false);    
  }
	
  for (const reco::PFJet &j : *jets) {
    if (j.pt() < 20) continue;
    //if (fabs(j.eta()) > 2.4) continue;
    //-------------------
    //Fill Jet-Level Info
    //-------------------
    jetE[nJets] = j.energy();
    jetPt[nJets] = j.pt();
    jetEta[nJets] = j.eta();
    jetPhi[nJets] = j.phi();
    jetMass[nJets] = j.mass();

    TLorentzVector thisJet;
    thisJet.SetPtEtaPhiE(jetPt[nJets], jetEta[nJets], jetPhi[nJets], jetE[nJets]);
    //jetCISV = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

    jetJetArea[nJets] = j.jetArea();
    jetPileupE[nJets] = j.pileup();

    jetPileupIdFlag[nJets] = 0;
    jetPassIDLoose[nJets] = passJetID(&j, 0);
    jetPassIDTight[nJets] = passJetID(&j, 1);
    jetPassMuFrac[nJets]  = ( j.muonEnergyFraction() < 0.80 );
    jetPassEleFrac[nJets]  = ( j.electronEnergyFraction() < 0.90 );


    // if (useGen_) {
    //   jetPartonFlavor = j.partonFlavour();
    //   jetHadronFlavor = j.hadronFlavour();
    // }

    jetChargedEMEnergyFraction[nJets] = j.chargedEmEnergyFraction();
    jetNeutralEMEnergyFraction[nJets] = j.neutralEmEnergyFraction();
    jetChargedHadronEnergyFraction[nJets] = j.chargedHadronEnergyFraction();
    jetNeutralHadronEnergyFraction[nJets] = j.neutralHadronEnergyFraction();
    // jet_charged_hadron_multiplicity[nJets] = j.chargedHadronMultiplicity();
    // jet_neutral_hadron_multiplicity[nJets] = j.neutralHadronMultiplicity();
    //jet_photon_multiplicity[nJets] = j.photonMultiplicity();
    //jet_electron_multiplicity[nJets] = j.electronMultiplicity();
    //jet_muon_multiplicity[nJets] = j.muonMultiplicity();
    //jet_HF_hadron_multiplicity[nJets] = j.HFHadronMultiplicity();
    //jet_HF_em_multiplicity[nJets] = j.HFEMMultiplicity();
    // jet_charged_multiplicity[nJets] = j.chargedMultiplicity();
    // jet_neutral_multiplicity[nJets] = j.neutralMultiplicity();

    //---------------------------
    //Trackless variables
    //---------------------------
    float alphaMax(0.0),medianTheta2D(0.0),medianIP(0.0),minDeltaRAllTracks(0.0),minDeltaRPVTracks(0.0),ptAllTracks(0.0), ptAllPVTracks(0.0);
    int nTracksPV(0);

    findTrackingVariables(thisJet,iSetup,alphaMax,medianTheta2D,medianIP,nTracksPV,ptAllPVTracks,ptAllTracks, minDeltaRAllTracks, minDeltaRPVTracks);

    float alphaMax_wp(0.0),medianTheta2D_wp(0.0),medianIP_wp(0.0),minDeltaRAllTracks_wp(0.0),minDeltaRPVTracks_wp(0.0),ptAllTracks_wp(0.0), ptAllPVTracks_wp(0.0);
    int nTracksPV_wp(0);
    findTrackingVariablesWithoutPropagator(thisJet,iSetup,alphaMax_wp,medianTheta2D_wp,medianIP_wp,nTracksPV_wp,ptAllPVTracks_wp,ptAllTracks_wp, minDeltaRAllTracks_wp, minDeltaRPVTracks_wp);

    //jetCISV = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    jetAlphaMax[nJets] = alphaMax;
    jetBetaMax[nJets] = alphaMax * ptAllTracks/(j.pt());
    jetGammaMax[nJets] = alphaMax * ptAllTracks/(j.energy());
    jetGammaMax_EM[nJets] = alphaMax * ptAllTracks/(j.energy()*(j.chargedEmEnergyFraction()+j.neutralEmEnergyFraction()));
    jetGammaMax_Hadronic[nJets] = alphaMax * ptAllTracks/(j.energy()*(j.chargedHadronEnergyFraction()+j.neutralHadronEnergyFraction()));
    jetGammaMax_ET[nJets] = alphaMax * ptAllTracks/j.et();
    jetMedianTheta2D[nJets] = medianTheta2D;
    jetMedianIP[nJets] = medianIP;
    jetPtAllPVTracks[nJets] = ptAllPVTracks;
    jetPtAllTracks[nJets] = ptAllTracks;
    jetMinDeltaRAllTracks[nJets] = minDeltaRAllTracks;
    jetMinDeltaRPVTracks[nJets] = minDeltaRPVTracks;

    jetAlphaMax_wp[nJets] = alphaMax_wp;
    jetBetaMax_wp[nJets] = alphaMax_wp * ptAllTracks_wp/(j.pt());
    jetGammaMax_wp[nJets] = alphaMax_wp * ptAllTracks_wp/(j.energy());
    jetGammaMax_EM_wp[nJets] = alphaMax_wp * ptAllTracks_wp/(j.energy()*(j.chargedEmEnergyFraction()+j.neutralEmEnergyFraction()));
    jetGammaMax_Hadronic_wp[nJets] = alphaMax_wp * ptAllTracks_wp/(j.energy()*(j.chargedHadronEnergyFraction()+j.neutralHadronEnergyFraction()));
    jetGammaMax_ET_wp[nJets] = alphaMax_wp * ptAllTracks_wp/j.et();
    jetMedianTheta2D_wp[nJets] = medianTheta2D_wp;
    jetMedianIP_wp[nJets] = medianIP_wp;
    jetPtAllPVTracks_wp[nJets] = ptAllPVTracks_wp;
    jetPtAllTracks_wp[nJets] = ptAllTracks_wp;
    jetMinDeltaRAllTracks_wp[nJets] = minDeltaRAllTracks_wp;
    jetMinDeltaRPVTracks_wp[nJets] = minDeltaRPVTracks_wp;

     //---------------------------
    //find photons inside the jet
    //---------------------------
    /*
    for (const reco::Photon &pho : *photons) {
      //cout << "Nphoton: " << fJetNPhotons << "\n";

      if (!(deltaR(pho.eta(), pho.phi() , j.eta(), j.phi()) < 0.5)) continue;


      fJetPhotonPt[fJetNPhotons]  = pho.pt();
      fJetPhotonEta[fJetNPhotons] = pho.eta(); //correct this for the vertex
      fJetPhotonPhi[fJetNPhotons] = pho.phi(); //correct this for the vertex

      fJetPhotonSeedRecHitE[fJetNPhotons]      = pho.superCluster()->seed()->x();
      fJetPhotonSeedRecHitEta[fJetNPhotons]      = pho.superCluster()->seed()->y();
      fJetPhotonSeedRecHitPhi[fJetNPhotons]      = pho.superCluster()->seed()->z();
      fJetPhotonSeedRecHitTime[fJetNPhotons]      = pho.superCluster()->seed()->energy();

      // //get time coordinate for the seed
      // for (const reco::PFCluster &pfcluster : *pfClusters) {
      // 	if(pfcluster.seed() == pho.superCluster()->seed()->seed())
      // 	  {
      // 	    pho_superClusterSeedT[fJetNPhotons] = pfcluster.time();
      // 	    pho_pfClusterSeedE[fJetNPhotons]      = pfcluster.energy();
      // 	  }
      // }

      //-------------------------------
      //fill all rechits inside photons
      //-------------------------------

      fJetNPhotons++;

    }
    */

    //---------------------------
    //Find RecHits Inside the Jet
    //---------------------------
    // geometry (from ECAL ELF)

    //double ecal_radius = 129.0;
    int n_matched_rechits = 0;
    std::vector<double> rechitphi;
    std::vector<double> rechiteta;
    std::vector<double> rechitet;
    std::vector<double> rechitt;
    
    for (uint q=0; q<ebRecHits->size(); q++) {
    // for (EcalRecHitCollection::const_iterator recHit = ebRecHits->begin(); recHit != ebRecHits->end(); ++recHit) {
      const EcalRecHit *recHit = &(*ebRecHits)[q];
      const DetId recHitId = recHit->detid();
      const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
      
      //save the rechits that are within DR 0.5 of the jet axis
      if ( jetPt[nJets] > 30 && 
	   deltaR(jetEta[nJets], jetPhi[nJets], recHitPos.eta(), recHitPos.phi())  < 0.5
	   && recHit->energy() > 0.2	  
	  ) {
	SaveThisEBRechit[q] = true;
	//cout << "Save this Rechit: " << q << " | " << SaveThisEBRechit[q] << " : " << recHit->energy() << " " << recHitPos.eta() << " " << recHitPos.phi() << "\n";
      }

      if (recHit->checkFlag(EcalRecHit::kSaturated) || recHit->checkFlag(EcalRecHit::kLeadingEdgeRecovered) || recHit->checkFlag(EcalRecHit::kPoorReco) || recHit->checkFlag(EcalRecHit::kWeird) || recHit->checkFlag(EcalRecHit::kDiWeird)) continue;
      if (recHit->timeError() < 0 || recHit->timeError() > 100) continue;
      if (abs(recHit->time()) > 12.5) continue;

      //Calculate jet timestamps
      if ( jetPt[nJets] > 30 && deltaR(jetEta[nJets], jetPhi[nJets], recHitPos.eta(), recHitPos.phi())  < 0.4) {
        //double rechit_x = ecal_radius * cos(recHitPos.phi());
        //double rechit_y = ecal_radius * sin(recHitPos.phi());
        //double rechit_z = ecal_radius * sinh(recHitPos.eta());
        //double photon_pv_travel_time = (1./30) * sqrt(pow(pvX-rechit_x,2)+pow(pvY-rechit_y,2)+pow(pvZ-rechit_z,2));
        if ( deltaR(jetEta[nJets], jetPhi[nJets], recHitPos.eta(), recHitPos.phi())  < 0.15) {
          if (recHit->energy() > Rechit_cut) jet_energy_frac[nJets] += recHit->energy();
        }

        if (recHit->energy() > Rechit_cut) {
	  // jetRechitT_Error[nJets] += 0.0;
          jetRechitE_Error[nJets] += recHit->energyError() * recHit->energyError();
          jetRechitE[nJets] += recHit->energy();
          jetRechitT[nJets] += recHit->time()*recHit->energy();
          jetRechitT_rms[nJets] += recHit->time()*recHit->time();
          rechitphi.push_back(recHitPos.phi());
          rechiteta.push_back(recHitPos.eta());
          rechitet.push_back(recHit->energy()/cosh(recHitPos.eta()));
          rechitt.push_back(recHit->time());
          n_matched_rechits++;
        }
      }
    }//loop over EB rechits

    for (uint q=0; q<eeRecHits->size(); q++) {
      const EcalRecHit *recHit = &(*eeRecHits)[q];
      const DetId recHitId = recHit->detid();
      const auto recHitPos = endcapGeometry->getGeometry(recHitId)->getPosition();
      
      //save the rechits that are within DR 0.5 of the jet axis
      if ( jetPt[nJets] > 30 && 
	   deltaR(jetEta[nJets], jetPhi[nJets], recHitPos.eta(), recHitPos.phi())  < 0.5
	   && recHit->energy() > 0.2	  
	  ) {
	SaveThisEERechit[q] = true;
	//cout << "Save this Rechit: " << q << " | " << SaveThisEERechit[q] << " : " << recHit->energy() << " " << recHitPos.eta() << " " << recHitPos.phi() << "\n";
      }
    }//loop over EE rechits

    jetRechitT[nJets] = jetRechitT[nJets]/jetRechitE[nJets];
    jetNRechits[nJets] = n_matched_rechits;
    jetRechitE_Error[nJets] = sqrt(jetRechitE_Error[nJets]);
    jetRechitT_rms[nJets] = sqrt(jetRechitT_rms[nJets]);
    double sig1(0.0),sig2(0.0);
    jet_second_moments(rechitet,rechiteta,rechitphi,sig1,sig2);
    jet_sig_et1[nJets] = sig1;
    jet_sig_et2[nJets] = sig2;
    jet_energy_frac[nJets] = jet_energy_frac[nJets]/jetRechitE[nJets];

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

      if ( jetPt[nJets] > 30 && 
	   deltaR(jetEta[nJets], jetPhi[nJets], hiteta, hitphi)  < 0.5	   
	   ) {
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
      if ( jetPt[nJets] > 30 && 
	   deltaR(jetEta[nJets], jetPhi[nJets], hiteta, hitphi)  < 0.5	   
	   ) {
	SaveThisHORechit[iHit] = true;
      }
    }

    //loop over tracks
    for (unsigned int iTrack = 0; iTrack < generalTracks->size(); iTrack ++){
      reco::Track generalTrack = generalTracks->at(iTrack);      
      if ( jetPt[nJets] > 30 && 
	   deltaR(jetEta[nJets], jetPhi[nJets], generalTrack.eta(), generalTrack.phi())  < 0.5
	   && generalTrack.pt() > 1	   
	  ) {
	SaveThisTrack[iTrack] = true;
      }
    } 


    nJets++;
  } //loop over jets





  //********************************************************
  // AK8 Jets
  //********************************************************
  for (const pat::Jet &j : *jetsAK8) {  
  //for (const reco::PFJet &j : *jetsAK8) {
    if (j.pt() < 50) continue;
    //-------------------
    //Fill Jet-Level Info
    //-------------------
  
    fatJetE[nFatJets] = j.correctedP4(0).E();
    fatJetPt[nFatJets] = j.correctedP4(0).Pt();
    fatJetEta[nFatJets] = j.correctedP4(0).Eta();
    fatJetPhi[nFatJets] = j.correctedP4(0).Phi();    
    fatJetSoftDropM[nFatJets] = (float) j.userFloat("ak8PFJetsCHSSoftDropMass");
    fatJetTau1[nFatJets] =  (float) j.userFloat("NjettinessAK8CHS:tau1");
    fatJetTau2[nFatJets] =  (float) j.userFloat("NjettinessAK8CHS:tau2");
    fatJetTau3[nFatJets] =  (float) j.userFloat("NjettinessAK8CHS:tau3");



    for (uint q=0; q<ebRecHits->size(); q++) {
      const EcalRecHit *recHit = &(*ebRecHits)[q];
      const DetId recHitId = recHit->detid();
      const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
      
      //save the rechits that are within DR 0.5 of the jet axis
      if ( deltaR(fatJetEta[nFatJets], fatJetPhi[nFatJets], recHitPos.eta(), recHitPos.phi())  < 1.0
    	   && recHit->energy() > 0.2
    	   && recHit->timeError() > 0.0 
    	   && recHit->timeError() < 100.0
    	   ) {
    	SaveThisEBRechit[q] = true;
      }
    }

    for (uint q=0; q<eeRecHits->size(); q++) {
      const EcalRecHit *recHit = &(*eeRecHits)[q];
      const DetId recHitId = recHit->detid();
      const auto recHitPos = endcapGeometry->getGeometry(recHitId)->getPosition();
      
      //save the rechits that are within DR 0.5 of the jet axis
      if ( deltaR(fatJetEta[nFatJets], fatJetPhi[nFatJets], recHitPos.eta(), recHitPos.phi())  < 1.0
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
      
      if ( deltaR(fatJetEta[nFatJets], fatJetPhi[nFatJets], hiteta, hitphi)  < 1.0 ) {
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
      if ( deltaR(fatJetEta[nFatJets], fatJetPhi[nFatJets], hiteta, hitphi)  < 1.0	   
	   ) {
	SaveThisHORechit[iHit] = true;
      }
    }



    //loop over tracks
    for (unsigned int iTrack = 0; iTrack < generalTracks->size(); iTrack ++){
      reco::Track generalTrack = generalTracks->at(iTrack);      
      if ( deltaR(fatJetEta[nFatJets], fatJetPhi[nFatJets], generalTrack.eta(), generalTrack.phi())  < 1.0
	   && generalTrack.pt() > 1	   
	   ) {
	SaveThisTrack[iTrack] = true;
      }
    }
    
    
    
    nFatJets++;
  }

  //********************************************************
  // Save EB Rechits inside Jets and AK8 Jets
  //********************************************************
  for (uint q=0; q<ebRecHits->size(); q++) { 
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
    } 
    //cout << "After Rechit: " << q << " | " << SaveThisEBRechit[q] << " : " << recHit->energy() << " " << recHitPos.eta() << " " << recHitPos.phi() << "\n";    
  }
  
  //********************************************************
  // Save EE Rechits inside Jets and AK8 Jets
  //********************************************************
  for (uint q=0; q<eeRecHits->size(); q++) { 
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
    }
    //cout << "After Rechit: " << q << " | " << SaveThisEBRechit[q] << " : " << recHit->energy() << " " << recHitPos.eta() << " " << recHitPos.phi() << "\n";    
  }

  //********************************************************
  // Save Tracks inside Jets and AK8 Jets
  //********************************************************
  // Magnetic field
  edm::ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  magneticField_ = &*magneticField;
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
      reco::TransientTrack tt(generalTrackHandle->at(iTrack),magneticField_);
    
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

      track_dxyToBS[nTracks] = tref->dxy(*beamSpot);
      track_dxyErr[nTracks] = tref->dxyError();
      track_dzToPV[nTracks] = tref->dz(beamSpot->position());
      track_dzErr[nTracks] = tref->dzError();

      //********************************************************
      // For track angle
      //********************************************************
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



      nTracks++; //increment saved tracks
    } //end if save this Track    
  } //loop over tracks

  //********************************************************
  // Save HCAL Rechits inside Jets and AK8 Jets
  //********************************************************
  for (unsigned int iHit = 0; iHit < hcalRecHitsHBHE->size(); iHit ++){
    const HBHERecHit *recHit = &(*hcalRecHitsHBHE)[iHit];
    const HcalDetId recHitId = recHit->detid();

    //cout << "HCALREchit " << iHit << " : " << recHit->detid().subdetId() << " : " << recHitId.depth() << " " << recHitId.ieta() << " " << recHitId.iphi() << " "
	 // << " | " << recHit->energy() << " " 
	 // << "\n";
    hbheRechit_iEta[nHBHERechits]  = recHitId.ieta();
    hbheRechit_iPhi[nHBHERechits]  = recHitId.iphi();
    hbheRechit_depth[nHBHERechits]  = recHitId.depth();

    if (SaveThisHCALRechit[iHit]) {
      if (recHit->energy() < 0.1) continue;      
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
    
      hbheRechit_E[nHORechits] = recHit->energy();
      hbheRechit_T[nHORechits] = recHit->time();
    
      nHBHERechits++;
    }
  }
  
  //********************************************************
  // Save HO Rechits inside Jets and AK8 Jets
  //********************************************************
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
    }
  }

  return true;
};

bool displacedJetMuon_ntupler::fillMuons(const edm::Event& iEvent)
{
  for(const pat::Muon &mu : *muons) {
    muonE[nMuons] = mu.energy();
    muonPt[nMuons] = mu.pt();
    muonEta[nMuons] = mu.eta();
    muonPhi[nMuons] = mu.phi();
    muonCharge[nMuons] = mu.charge();
    muonIsLoose[nMuons] = mu.isLooseMuon();
    muonIsMedium[nMuons] = mu.isMediumMuon();
    muonIsTight[nMuons] = mu.isTightMuon(*myPV);
    muon_d0[nMuons] = -mu.muonBestTrack()->dxy(myPV->position());
    muon_dZ[nMuons] = mu.muonBestTrack()->dz(myPV->position());
    muon_ip3d[nMuons] = mu.dB(pat::Muon::PV3D);
    muon_ip3dSignificance[nMuons] = mu.dB(pat::Muon::PV3D)/mu.edB(pat::Muon::PV3D);
    muonType[nMuons] = mu.isMuon() + mu.isGlobalMuon() + mu.isTrackerMuon() + mu.isStandAloneMuon()
      + mu.isCaloMuon() + mu.isPFMuon() + mu.isRPCMuon();
    muonQuality[nMuons] =
      muon::isGoodMuon(mu,muon::All)
    + muon::isGoodMuon(mu,muon::AllGlobalMuons)
    + muon::isGoodMuon(mu,muon::AllStandAloneMuons)
    + muon::isGoodMuon(mu,muon::AllTrackerMuons)
    + muon::isGoodMuon(mu,muon::TrackerMuonArbitrated)
    + muon::isGoodMuon(mu,muon::AllArbitrated)
    + muon::isGoodMuon(mu,muon::GlobalMuonPromptTight)
    + muon::isGoodMuon(mu,muon::TMLastStationLoose)
    + muon::isGoodMuon(mu,muon::TMLastStationTight)
    + muon::isGoodMuon(mu,muon::TM2DCompatibilityLoose)
    + muon::isGoodMuon(mu,muon::TM2DCompatibilityTight)
    + muon::isGoodMuon(mu,muon::TMOneStationLoose)
    + muon::isGoodMuon(mu,muon::TMOneStationTight)
    + muon::isGoodMuon(mu,muon::TMLastStationOptimizedLowPtLoose)
    + muon::isGoodMuon(mu,muon::TMLastStationOptimizedLowPtTight)
    + muon::isGoodMuon(mu,muon::GMTkChiCompatibility)
    + muon::isGoodMuon(mu,muon::GMStaChiCompatibility)
    + muon::isGoodMuon(mu,muon::GMTkKinkTight)
    + muon::isGoodMuon(mu,muon::TMLastStationAngLoose)
    + muon::isGoodMuon(mu,muon::TMLastStationAngTight)
    + muon::isGoodMuon(mu,muon::TMOneStationAngLoose)
    + muon::isGoodMuon(mu,muon::TMOneStationAngTight)
    + muon::isGoodMuon(mu,muon::TMLastStationOptimizedBarrelLowPtLoose)
    + muon::isGoodMuon(mu,muon::TMLastStationOptimizedBarrelLowPtTight)
    + muon::isGoodMuon(mu,muon::RPCMuLoose);
    muon_pileupIso[nMuons] = mu.pfIsolationR04().sumPUPt;
    muon_chargedIso[nMuons] = mu.pfIsolationR04().sumChargedHadronPt;
    muon_photonIso[nMuons] = mu.pfIsolationR04().sumPhotonEt;
    muon_neutralHadIso[nMuons] = mu.pfIsolationR04().sumNeutralHadronEt;
    //muon_ptrel[nMuons] = getLeptonPtRel( jets, &mu );
    //tuple<double,double,double> PFMiniIso = getPFMiniIsolation(packedPFCands, dynamic_cast<const reco::Candidate *>(&mu), 0.05, 0.2, 10., false, false);
    //muon_chargedMiniIso[nMuons] = std::get<0>(PFMiniIso);
    //muon_photonAndNeutralHadronMiniIso[nMuons] = std::get<1>(PFMiniIso);
    //muon_chargedPileupMiniIso[nMuons] = std::get<2>(PFMiniIso);
    //muon_activityMiniIsoAnnulus[nMuons] = ActivityPFMiniIsolationAnnulus( packedPFCands, dynamic_cast<const reco::Candidate *>(&mu), 0.4, 0.05, 0.2, 10.);
    muon_validFractionTrackerHits[nMuons] = (mu.innerTrack().isNonnull() ? mu.track()->validFraction() : -99.0);
    muon_isGlobal[nMuons] = muon::isGoodMuon(mu,muon::AllGlobalMuons);
    muon_normChi2[nMuons] = ( muon::isGoodMuon(mu,muon::AllGlobalMuons) ? mu.globalTrack()->normalizedChi2() : -99.0);
    muon_chi2LocalPosition[nMuons] = mu.combinedQuality().chi2LocalPosition;
    muon_kinkFinder[nMuons] = mu.combinedQuality().trkKink;
    muon_segmentCompatability[nMuons] = muon::segmentCompatibility(mu);
    bool isGoodGlobal = mu.isGlobalMuon() && mu.globalTrack()->normalizedChi2() < 3 && mu.combinedQuality().chi2LocalPosition < 12 && mu.combinedQuality().trkKink < 20;
    muonIsICHEPMedium[nMuons] = muon::isLooseMuon(mu) && muon_validFractionTrackerHits[nMuons] > 0.49 && muon::segmentCompatibility(mu) > (isGoodGlobal ? 0.303 : 0.451);
    //-----------------------
    //Trigger Object Matching
    //-----------------------
    bool passTagMuonFilter = false;
    /*for (pat::TriggerObjectStandAlone trigObject : *triggerObjects) {
      std::cout << "muon debug 18 " << nMuons << std::endl;
      if (deltaR(trigObject.eta(), trigObject.phi(),mu.eta(),mu.phi()) > 0.3) continue;
      trigObject.unpackFilterLabels(iEvent, *triggerBits);
      std::cout << "muon debug 19 " << nMuons << std::endl;
      //check single muon filters
      if ( trigObject.hasFilterLabel("hltL3fL1sMu25L1f0Tkf27QL3trkIsoFiltered0p09") ||
    	   trigObject.hasFilterLabel("hltL3fL1sMu20Eta2p1L1f0Tkf24QL3trkIsoFiltered0p09") ||
    	   trigObject.hasFilterLabel("hltL3fL1sMu16Eta2p1L1f0Tkf20QL3trkIsoFiltered0p09") ||
    	   trigObject.hasFilterLabel("hltL3fL1sMu16L1f0Tkf20QL3trkIsoFiltered0p09") ||
    	   trigObject.hasFilterLabel("hltL3crIsoL1sMu25L1f0L2f10QL3f27QL3trkIsoFiltered0p09") ||
    	   trigObject.hasFilterLabel("hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09") ||
    	   trigObject.hasFilterLabel("hltL3crIsoL1sMu16Eta2p1L1f0L2f10QL3f20QL3trkIsoFiltered0p09") ||
    	   trigObject.hasFilterLabel("hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09")
    	   ) passTagMuonFilter = true;
         std::cout << "muon debug 20 " << nMuons << std::endl;
      //check all filters
      for ( int q=0; q<MAX_MuonHLTFilters;q++) {
    	if (trigObject.hasFilterLabel(muonHLTFilterNames[q].c_str())) muon_passHLTFilter[nMuons][q] = true;
      std::cout << "muon debug 21 " << nMuons << std::endl;
      }

    }*/

    muon_passSingleMuTagFilter[nMuons] = passTagMuonFilter;
    nMuons++;
  }

  return true;
};


bool displacedJetMuon_ntupler::passJetID( const reco::PFJet *jet, int cutLevel) {
  bool result = false;

  double NHF = jet->neutralHadronEnergyFraction();
  double NEMF = jet->neutralEmEnergyFraction();
  int NumConst = jet->chargedMultiplicity() + jet->neutralMultiplicity() ;
  double CHF = jet->chargedHadronEnergyFraction();
  double MUF = jet->muonEnergyFraction();
  double CEMF = jet->chargedEmEnergyFraction();
  int NumNeutralParticles =jet->neutralMultiplicity();
  int CHM = jet->chargedMultiplicity();

  //Loose
  if (cutLevel == 0) {
    if ( fabs(jet->eta()) <= 2.4) {
      if ( NHF  < 0.99 && NEMF < 0.99 && NumConst > 1
	   && CHF > 0 && CHM > 0 && CEMF < 0.99 ) result = true;
    } else if( fabs(jet->eta()) <= 3.0)  {
      if ( NHF  < 0.99 && NEMF < 0.99 && NumConst > 1 ) result = true;
    } else {
      if ( NEMF < 0.90 && NumNeutralParticles > 10 ) result = true;
    }
  }

  //Tight
  else if (cutLevel == 1) {
    if ( fabs(jet->eta()) <= 2.4) {
      if ( NHF  < 0.90 && NEMF < 0.90 && NumConst > 1
	   && CHF > 0 && CHM > 0 && CEMF < 0.99 ) result = true;
    } else if( fabs(jet->eta()) <= 3.0)  {
      if ( NHF  < 0.90 && NEMF < 0.90 && NumConst > 1 ) result = true;
    } else {
      if ( NEMF < 0.90 && NumNeutralParticles > 10 ) result = true;
    }
  }

  //Tight Lep Veto
  else if (cutLevel == 2) {
    if ( fabs(jet->eta()) <= 2.4) {
      if ( NHF  < 0.90 && NEMF < 0.90 && NumConst > 1
	   && CHF > 0 && CHM > 0 && CEMF < 0.99 && MUF < 0.8 ) result = true;
    } else if( fabs(jet->eta()) <= 3.0)  {
      if ( NHF  < 0.90 && NEMF < 0.90 && NumConst > 1 ) result = true;
    } else {
      if ( NEMF < 0.90 && NumNeutralParticles > 10 ) result = true;
    }
  }

  return result;
}//passJetID PFJet

void displacedJetMuon_ntupler::findTrackingVariablesWithoutPropagator(const TLorentzVector &jetVec,const edm::EventSetup& iSetup,float &alphaMax,float &medianTheta2D,float &medianIP, int &nTracksPV,float &ptAllPVTracks,float &ptAllTracks,float &minDeltaRAllTracks, float &minDeltaRPVTracks)
{
  int nTracksAll = 0;
  //Displaced jet stuff
  double ptPVTracksMax = 0.;
  minDeltaRAllTracks = 15;
  minDeltaRPVTracks = 15;
  reco::Vertex primaryVertex = vertices->at(0);
  std::vector<double> theta2Ds;
  std::vector<double> IP2Ds;
 
  for (unsigned int iTrack = 0; iTrack < generalTracks->size(); iTrack ++){
  	reco::Track generalTrack = generalTracks->at(iTrack);
  	TLorentzVector generalTrackVecTemp;
  	generalTrackVecTemp.SetPtEtaPhiM(generalTrack.pt(),generalTrack.eta(),generalTrack.phi(),0);

  	if (generalTrack.pt() > 1) {
	    if (minDeltaRAllTracks > generalTrackVecTemp.DeltaR(jetVec)) {
	      minDeltaRAllTracks =  generalTrackVecTemp.DeltaR(jetVec);
	    }
	    if (generalTrackVecTemp.DeltaR(jetVec) < 0.4){
    		nTracksAll ++;
    		//tot pt for alpha
    		ptAllTracks += generalTrack.pt();

    		// theta 2d
    		// ROOT::Math::XYZPoint innerPos = generalTrack.innerPosition();
    		// ROOT::Math::XYZPoint vertexPos = primaryVertex.position();
    		// ROOT::Math::XYZVector deltaPos = innerPos - vertexPos;
    		// ROOT::Math::XYZVector momentum = generalTrack.innerMomentum();
    		// double mag2DeltaPos = TMath::Sqrt((deltaPos.x()*deltaPos.x()) + (deltaPos.y()*deltaPos.y()));
    		// double mag2Mom = TMath::Sqrt((momentum.x()*momentum.x()) + (momentum.y()*momentum.y()));
    		// double theta2D = TMath::ACos((deltaPos.x()*momentum.x()+deltaPos.y()*momentum.y())/(mag2Mom*mag2DeltaPos));
    		// theta2Ds.push_back(theta2D);

    		//IP sig
    		edm::ESHandle<TransientTrackBuilder> theB;
    		iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
    		reco::TransientTrack transTrack = theB->build(generalTrack);
    		TrajectoryStateClosestToBeamLine traj = transTrack.stateAtBeamLine();
    		Measurement1D meas = traj.transverseImpactParameter();
    		std::pair<bool, Measurement1D> ip2d = IPTools::absoluteTransverseImpactParameter(transTrack,primaryVertex);
    		IP2Ds.push_back(ip2d.second.value()/ip2d.second.error());
      }
    }
  }
  if (ptAllTracks > 0.9){
	//No matched jets
  	for (auto vertex = vertices->begin(); vertex != vertices->end(); vertex++){
      double ptPVTracks = 0.;
      int nTracksPVTemp = 0;
      if(!vertex->isValid())continue;
      if (vertex->isFake())continue;
      for(auto pvTrack=vertex->tracks_begin(); pvTrack!=vertex->tracks_end(); pvTrack++){
    	TLorentzVector pvTrackVecTemp;
    	pvTrackVecTemp.SetPtEtaPhiM((*pvTrack)->pt(),(*pvTrack)->eta(),(*pvTrack)->phi(),0);
  	//If pv track associated with jet add pt to ptPVTracks
    	if ((*pvTrack)->pt() > 1) {
    	    if (minDeltaRPVTracks > pvTrackVecTemp.DeltaR(jetVec))
    	    {
    		     minDeltaRPVTracks =  pvTrackVecTemp.DeltaR(jetVec);
    	    }
    	    if (pvTrackVecTemp.DeltaR(jetVec) < 0.4){
        		ptPVTracks += (*pvTrack)->pt();
        		ptAllPVTracks += (*pvTrack)->pt();
        		nTracksPVTemp++;
    	    }
    		}
      }
      if (ptPVTracks > ptPVTracksMax) {
      	ptPVTracksMax = ptPVTracks;
      	nTracksPV = nTracksPVTemp;
      }
      alphaMax = ptPVTracksMax/ptAllTracks;
  	}
  }
  std::sort(IP2Ds.begin(),IP2Ds.end());
  if (IP2Ds.size() > 0){
	   medianIP = IP2Ds[IP2Ds.size()/2];
  }
  std::sort(theta2Ds.begin(),theta2Ds.end());
  if (theta2Ds.size() > 0){
    medianTheta2D = theta2Ds[theta2Ds.size()/2];
  }
};

void displacedJetMuon_ntupler::findTrackingVariables(const TLorentzVector &jetVec,const edm::EventSetup& iSetup,float &alphaMax,float &medianTheta2D,float &medianIP, int &nTracksPV,float &ptAllPVTracks,float &ptAllTracks,float &minDeltaRAllTracks, float &minDeltaRPVTracks)
{
  int nTracksAll = 0;
  //Displaced jet stuff
  double ptPVTracksMax = 0.;
  minDeltaRAllTracks = 15;
  minDeltaRPVTracks = 15;
  reco::Vertex primaryVertex = vertices->at(0);
  std::vector<double> theta2Ds;
  std::vector<double> IP2Ds;

  // propagator
  //edm::ESHandle<Propagator> thePropagator_;
  edm::ESTransientHandle<Propagator> thePropagator_;
  std::string thePropagatorName_ = "PropagatorWithMaterial";
  iSetup.get<TrackingComponentsRecord>().get(thePropagatorName_,thePropagator_);
  StateOnTrackerBound stateOnTracker(thePropagator_.product());
  
  const MagneticField* magneticField_;
  //edm::ESHandle<MagneticField> magneticField;
  edm::ESTransientHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField); 
  magneticField_ = &*magneticField; 
  //std::cout << "B " << magneticField_ << " tracks size " << tracks->size() << std::endl;
  
  for (unsigned int iTrack = 0; iTrack < generalTracks->size(); iTrack ++){

    //track propagation
    FreeTrajectoryState fts = trajectoryStateTransform::initialFreeState (tracks->at(iTrack),magneticField_); 
    TrajectoryStateOnSurface outer = stateOnTracker(fts); 
    if(!outer.isValid()) continue; 
    GlobalPoint outerPos = outer.globalPosition();
    
  	TLorentzVector generalTrackVecTemp;
  	generalTrackVecTemp.SetPtEtaPhiM((generalTracks->at(iTrack)).pt(), outerPos.eta(), outerPos.phi(), 0);

  	if ((generalTracks->at(iTrack)).pt() > 1) {
	    if (minDeltaRAllTracks > generalTrackVecTemp.DeltaR(jetVec))
	    {
		    minDeltaRAllTracks =  generalTrackVecTemp.DeltaR(jetVec);
	    }
	    if (generalTrackVecTemp.DeltaR(jetVec) < 0.4){
    		nTracksAll ++;
    		//tot pt for alpha
    		ptAllTracks += (generalTracks->at(iTrack)).pt();

      }
    }
  }
  if (ptAllTracks > 0.9){
      //No matched jets
      for (auto vertex = vertices->begin(); vertex != vertices->end(); vertex++){
      double ptPVTracks = 0.;
      int nTracksPVTemp = 0;
      if(!vertex->isValid())continue;
      if (vertex->isFake())continue;
      for(auto pvTrack=vertex->tracks_begin(); pvTrack!=vertex->tracks_end(); pvTrack++){
    
        //track propagation
        FreeTrajectoryState ftspv = trajectoryStateTransform::initialFreeState (**pvTrack, magneticField_); 
        TrajectoryStateOnSurface outerpv = stateOnTracker(ftspv); 
        if(!outerpv.isValid()) continue; 
        GlobalPoint outerpvPos = outerpv.globalPosition();
        
    	TLorentzVector pvTrackVecTemp;
    	pvTrackVecTemp.SetPtEtaPhiM((*pvTrack)->pt(),outerpvPos.eta(),outerpvPos.phi(),0);
  	//If pv track associated with jet add pt to ptPVTracks
    	if ((*pvTrack)->pt() > 1) {
    	    if (minDeltaRPVTracks > pvTrackVecTemp.DeltaR(jetVec))
    	    {
    		     minDeltaRPVTracks =  pvTrackVecTemp.DeltaR(jetVec);
    	    }
    	    if (pvTrackVecTemp.DeltaR(jetVec) < 0.4){
        		ptPVTracks += (*pvTrack)->pt();
        		ptAllPVTracks += (*pvTrack)->pt();
        		nTracksPVTemp++;
    	    }
    		}
      }
      if (ptPVTracks > ptPVTracksMax) {
      	ptPVTracksMax = ptPVTracks;
      	nTracksPV = nTracksPVTemp;
      }
      alphaMax = ptPVTracksMax/ptAllTracks;
  	}
  }
/*
  for (int iTrack = 0; iTrack < nTracks; iTrack ++){
  	TLorentzVector generalTrackVecTemp;
  	generalTrackVecTemp.SetPtEtaPhiM(TrackPt[iTrack], TrackEta[iTrack], TrackPhi[iTrack], 0);

  	if (TrackPt[iTrack] > 1) {
	    if (minDeltaRAllTracks > generalTrackVecTemp.DeltaR(jetVec))
	    {
		    minDeltaRAllTracks =  generalTrackVecTemp.DeltaR(jetVec);
	    }
	    if (generalTrackVecTemp.DeltaR(jetVec) < 0.4){
    		nTracksAll ++;
    		//tot pt for alpha
    		ptAllTracks += TrackPt[iTrack];

      	     }
         }
  }

  if (ptAllTracks > 0.9){
    //No matched jets
    for (int ipvTrack = 0; ipvTrack < npvTracks; ipvTrack++){
      double ptPVTracks = 0.;
      int nTracksPVTemp = 0;
      TLorentzVector pvTrackVecTemp;
      pvTrackVecTemp.SetPtEtaPhiM(PVTrackPt[ipvTrack], PVTrackEta[ipvTrack], PVTrackPhi[ipvTrack], 0);
      if (PVTrackPt[ipvTrack] > 1) {
    	    if (minDeltaRPVTracks > pvTrackVecTemp.DeltaR(jetVec))
    	    {
    		     minDeltaRPVTracks =  pvTrackVecTemp.DeltaR(jetVec);
    	    }
    	    if (pvTrackVecTemp.DeltaR(jetVec) < 0.4){
        		ptPVTracks += PVTrackPt[ipvTrack];
        		ptAllPVTracks += PVTrackPt[ipvTrack];
        		nTracksPVTemp++;
    	    }
    	}
      if (ptPVTracks > ptPVTracksMax) {
      	ptPVTracksMax = ptPVTracks;
      	nTracksPV = nTracksPVTemp;
      }
      alphaMax = ptPVTracksMax/ptAllTracks;
    }
  }
*/
  std::sort(IP2Ds.begin(),IP2Ds.end());
  if (IP2Ds.size() > 0){
	   medianIP = IP2Ds[IP2Ds.size()/2];
  }
  std::sort(theta2Ds.begin(),theta2Ds.end());
  if (theta2Ds.size() > 0){
    medianTheta2D = theta2Ds[theta2Ds.size()/2];
  }
};
void displacedJetMuon_ntupler::jet_second_moments(std::vector<double> &et,std::vector<double> &eta,std::vector<double> &phi,double &sig1,double &sig2)
{
  double mean_eta = 0.0;
  double mean_phi = 0.0;
  double et_squared = 0.0;
  for(unsigned int i = 0;i < eta.size();i++)
  {
    mean_eta += et[i]*et[i]*eta[i];
    mean_phi += et[i]*et[i]*phi[i];
    et_squared += et[i]*et[i];
  }
  mean_eta = mean_eta/et_squared;
  mean_phi = mean_phi/et_squared;

  double m11(0.0),m22(0.0),m12(0.0);
  for(unsigned int i = 0;i < eta.size();i++)
  {
    m11 += et[i]*et[i]*(eta[i]-mean_eta)*(eta[i]-mean_eta);
    m22 += et[i]*et[i]*(phi[i]-mean_phi)*(phi[i]-mean_phi);
    m12 += et[i]*et[i]*(phi[i]-mean_phi)*(eta[i]-mean_eta);
  }
  double a = (m11+m22)/2;
  double b = 0.5*sqrt(pow(m11+m22,2)-4*(m11*m22-m12*m12));
  sig1 = sqrt(abs(a+b)/et_squared);
  sig2 = sqrt(abs(a-b)/et_squared);

};

const reco::Candidate* displacedJetMuon_ntupler::findFirstMotherWithDifferentID(const reco::Candidate *particle){

  if( particle == 0 ){
    printf("ERROR! null candidate pointer, this should never happen\n");
    return 0;
  }

  // Is this the first parent with a different ID? If yes, return, otherwise
  // go deeper into recursion
  if (particle->numberOfMothers() > 0 && particle->pdgId() != 0) {
    if (particle->pdgId() == particle->mother(0)->pdgId()
	&& particle->mother(0)->status() != 11  // prevent infinite loop for sherpa documentation gluons
	) {
      return findFirstMotherWithDifferentID(particle->mother(0));
    } else {
      return particle->mother(0);
    }
  }

  return 0;
};

const reco::Candidate* displacedJetMuon_ntupler::findOriginalMotherWithSameID(const reco::Candidate *particle){

  if( particle == 0 ){
    printf("ERROR! null candidate pointer, this should never happen\n");
    return 0;
  }

  // Is there another parent with the same ID? If yes, go deeper into recursion
  if (particle->numberOfMothers() > 0 && particle->pdgId() != 0) {
    if (particle->mother(0)->numberOfMothers() == 0 ||
	particle->mother(0)->status() == 11 ||  // prevent infinite loop for sherpa documentation gluons
	(particle->mother(0)->numberOfMothers() > 0 && particle->mother(0)->mother(0)->pdgId() != particle->mother(0)->pdgId())
	) {
      return particle->mother(0);
    } else {
      return findOriginalMotherWithSameID(particle->mother(0));
    }
  }

  return 0;
}

bool displacedJetMuon_ntupler::fillMet(const edm::Event& iEvent)
{
  const reco::PFMET &Met = mets->front();

  sumMET = Met.sumEt();
  metType0Pt = 0;
  metType0Phi = 0;
  metType1Pt_raw = Met.pt();
  metType1Pt = Met.pt();
  metType1Px = Met.px();
  metType1Py = Met.py();
  metType1Eta = Met.eta();
  metType1Phi_raw = Met.phi();
  metType1Phi = Met.phi();
  metType0Plus1Pt = 0;
  metType0Plus1Phi = 0;

  //
  iEvent.getByToken(globalSuperTightHalo2016FilterToken_, globalSuperTightHalo2016Filter);
  iEvent.getByToken(globalTightHalo2016FilterToken_, globalTightHalo2016Filter);
  iEvent.getByToken(BadChargedCandidateFilterToken_, BadChargedCandidateFilter);
  iEvent.getByToken(BadPFMuonFilterToken_, BadPFMuonFilter);
  iEvent.getByToken(EcalDeadCellTriggerPrimitiveFilterToken_, EcalDeadCellTriggerPrimitiveFilter);
  iEvent.getByToken(HBHENoiseFilterToken_, HBHENoiseFilter);
  iEvent.getByToken(HBHEIsoNoiseFilterToken_, HBHEIsoNoiseFilter);
  iEvent.getByToken(ecalBadCalibFilterToken_, ecalBadCalibReducedMINIAODFilter);
  iEvent.getByToken(eeBadScFilterToken_, eeBadScFilter);
  // iEvent.getByToken(primaryVertexFilterToken_, primaryVertexFilter);


  Flag2_globalSuperTightHalo2016Filter = *globalSuperTightHalo2016Filter;
  Flag2_globalTightHalo2016Filter = *globalTightHalo2016Filter;
  Flag2_BadChargedCandidateFilter = *BadChargedCandidateFilter;
  Flag2_BadPFMuonFilter = *BadPFMuonFilter;
  Flag2_EcalDeadCellTriggerPrimitiveFilter = *EcalDeadCellTriggerPrimitiveFilter;
  Flag2_HBHENoiseFilter = *HBHENoiseFilter;
  Flag2_HBHEIsoNoiseFilter = *HBHEIsoNoiseFilter;
  Flag2_ecalBadCalibFilter = *ecalBadCalibReducedMINIAODFilter;
  Flag2_eeBadScFilter = *eeBadScFilter;
  // Flag2_goodVertices = *primaryVertexFilter;

  const edm::TriggerNames &metNames = iEvent.triggerNames(*metFilterBits);

  //*******************************************************************************
  //For Debug printout
  //*******************************************************************************
  // for (unsigned int i = 0, n = metFilterBits->size(); i < n; ++i) {
  // 	std::cout << "MET Filter " << metNames.triggerName(i).c_str() << "\n";
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


  return true;
};

bool displacedJetMuon_ntupler::fillTrigger(const edm::Event& iEvent)
{

  //fill trigger information
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  // std::cout << "\n === TRIGGER PATHS === " << std::endl;
  //------------------------------------------------------------------
  //Option to save all HLT path names in the ntuple per event
  //Expensive option in terms of ntuple size
  //------------------------------------------------------------------
  nameHLT->clear();
  //cout<<triggerBits->size()<<endl;
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i)
  {
    string hltPathNameReq = "HLT_";
    //if (triggerBits->accept(i))
    if ((names.triggerName(i)).find(hltPathNameReq) != string::npos) nameHLT->push_back(names.triggerName(i));
    /*
    std::cout << "Trigger " << names.triggerName(i) <<
    ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
    ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")
    << std::endl;
    if ((names.triggerName(i)).find(hltPathNameReq) != string::npos && triggerBits->accept(i)) std::cout << "Trigger " << names.triggerName(i) <<
    ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")
    << std::endl;
    */
  }
  //std::cout << "n triggers: " <<  nameHLT->size() << std::endl;
  //std::cout << "====================" << std::endl;
  //for ( unsigned int i = 0; i < nameHLT->size(); i++ )
  //{
  //  std::cout << i << " -> " << nameHLT->at(i) << std::endl;
  //}
  //------------------------------------------------------------------
  // Save trigger decisions in array of booleans
  //------------------------------------------------------------------

  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i)
  {
    string hltPathNameReq = "HLT_";
    if ((names.triggerName(i)).find(hltPathNameReq) == string::npos) continue;
    if ((names.triggerName(i)).find_last_of("_") == string::npos) continue;
    int lastUnderscorePos = (names.triggerName(i)).find_last_of("_");
    string hltPathNameWithoutVersionNumber = (names.triggerName(i)).substr(0,lastUnderscorePos);

    for (unsigned int j = 0; j < NTriggersMAX; ++j)
    {
      if (triggerPathNames[j] == "") continue;
      if (hltPathNameWithoutVersionNumber == triggerPathNames[j])
      {
        triggerDecision[j] = triggerBits->accept(i);
        //triggerHLTPrescale[j] = triggerPrescales->getPrescaleForIndex(i);
      }
    }
  }
  //------------------------------------------------------------------
  // Print Trigger Objects
  //------------------------------------------------------------------
/*
  for (pat::TriggerObjectStandAlone trigObject : *triggerObjects)
  {
    //cout << "triggerObj: " << trigObject.pt() << " " << trigObject.eta() << " " << trigObject.phi() << "\n";
    //bool foundRazor = false;
    //Need to unpack the filter labels before checking
    trigObject.unpackFilterLabels(iEvent, *triggerBits);
    for(int j=0; j<int(trigObject.filterLabels().size());j++)
    {
      //if ((trigObject.filterLabels())[j] == "hltRsqMR200Rsq0p0196MR100Calo") foundRazor = true;
      // trigObject.unpackPathNames(names);
      // cout << "filter: " << (trigObject.pathNames())[j] << " " << (trigObject.filterLabels())[j] << "\n";
      //cout << "filter: " << (trigObject.filterLabels())[j] << "\n";
    }
  }
*/
//define this as a plug-in
  return true;
};

bool displacedJetMuon_ntupler::fillMC() {
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
  
  return true;
};



bool displacedJetMuon_ntupler::fillElectrons(const edm::Event& iEvent)
{

  for(const pat::Electron &ele : *electrons) {
    if(ele.pt() < 5) continue;
    eleE[nElectrons] = ele.energy();
    elePt[nElectrons] = ele.pt();
    eleEta[nElectrons] = ele.eta();
    elePhi[nElectrons] = ele.phi();
    eleCharge[nElectrons] = ele.charge();

    eleE_SC[nElectrons] = ele.superCluster()->energy();
    eleEta_SC[nElectrons] = ele.superCluster()->eta();
    elePhi_SC[nElectrons] = ele.superCluster()->phi();

    eleSigmaIetaIeta[nElectrons] = ele.sigmaIetaIeta();
    eleFull5x5SigmaIetaIeta[nElectrons] = ele.full5x5_sigmaIetaIeta();
    eleR9[nElectrons] = ele.r9();
    ele_dEta[nElectrons] = ele.deltaEtaSuperClusterTrackAtVtx() - ele.superCluster()->eta() + ele.superCluster()->seed()->eta();

    ele_dPhi[nElectrons] = ele.deltaPhiSuperClusterTrackAtVtx();
    ele_HoverE[nElectrons] = ele.hcalOverEcal();
    ele_d0[nElectrons] = -ele.gsfTrack().get()->dxy(myPV->position());
    ele_dZ[nElectrons] = ele.gsfTrack().get()->dz(myPV->position());

    //ele_ip3d[nElectrons] = ((edm::Ptr<pat::Electron>)(ele))->dB(pat::Electron::PV3D);
    //ele_ip3dSignificance[nElectrons] = ((edm::Ptr<pat::Electron>)(ele))->dB(pat::Electron::PV3D)/((edm::Ptr<pat::Electron>)(ele))->edB(pat::Electron::PV3D);
    ele_pileupIso[nElectrons] = ele.pfIsolationVariables().sumPUPt;
    ele_chargedIso[nElectrons] = ele.pfIsolationVariables().sumChargedHadronPt;
    ele_photonIso[nElectrons] = ele.pfIsolationVariables().sumPhotonEt;
    ele_neutralHadIso[nElectrons] = ele.pfIsolationVariables().sumNeutralHadronEt;
    ele_MissHits[nElectrons] = ele.gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS);


    //---------------
    //Conversion Veto
    //---------------

    ele_PassConvVeto[nElectrons] = false;
    if( beamSpot.isValid() && conversions.isValid() )
    {
      ele_PassConvVeto[nElectrons] = !ConversionTools::hasMatchedConversion(ele,conversions,beamSpot->position());
    } else {
      cout << "\n\nERROR!!! conversions not found!!!\n";
    }
    // 1/E - 1/P
    if( ele.ecalEnergy() == 0 ){
      ele_OneOverEminusOneOverP[nElectrons] = 1e30;
    } else if( !std::isfinite(ele.ecalEnergy())){
      ele_OneOverEminusOneOverP[nElectrons] = 1e30;
    } else {
    ele_OneOverEminusOneOverP[nElectrons] = 1./ele.ecalEnergy()  -  ele.eSuperClusterOverP()/ele.ecalEnergy();
    }
    //----------------------
    //ID MVA
    //----------------------
    // ele_IDMVAGeneralPurpose[nElectrons] = (*mvaGeneralPurposeValues)[ele];
    // ele_IDMVACategoryGeneralPurpose[nElectrons] = (*mvaGeneralPurposeCategories)[ele];
    // ele_IDMVAHZZ[nElectrons] = (*mvaHZZValues)[ele];
    // ele_IDMVACategoryHZZ[nElectrons] = (*mvaHZZCategories)[ele];

    // ele_RegressionE[nElectrons] = ((edm::Ptr<pat::Electron>)(ele))->ecalRegressionEnergy();
    // ele_CombineP4[nElectrons]   = ((edm::Ptr<pat::Electron>)(ele))->ecalTrackRegressionEnergy();
    //
    // ele_ptrel[nElectrons]   = getLeptonPtRel( jets, &(*ele) );
    // tuple<double,double,double> PFMiniIso = getPFMiniIsolation(packedPFCands, dynamic_cast<const reco::Candidate *>(&(*ele)), 0.05, 0.2, 10., false, false);
    // ele_chargedMiniIso[nElectrons] = std::get<0>(PFMiniIso);
    // ele_photonAndNeutralHadronMiniIso[nElectrons] = std::get<1>(PFMiniIso);
    // ele_chargedPileupMiniIso[nElectrons] = std::get<2>(PFMiniIso);
    // ele_activityMiniIsoAnnulus[nElectrons] = ActivityPFMiniIsolationAnnulus( packedPFCands, dynamic_cast<const reco::Candidate *>(&(*ele)), 0.4, 0.05, 0.2, 10.);

    //-----------------------
    //Trigger Object Matching
    //-----------------------
    bool passSingleEleTagFilter = false;
    bool passTPOneTagFilter= false;
    bool passTPTwoTagFilter= false;
    bool passTPOneProbeFilter= false;
    bool passTPTwoProbeFilter= false;
    /*
    for (pat::TriggerObjectStandAlone trigObject : *triggerObjects)
    {
      if (deltaR(trigObject.eta(), trigObject.phi(),ele.eta(),ele.phi()) > 0.3) continue;
      trigObject.unpackFilterLabels(iEvent, *triggerBits);

      //check Single ele filters
      if (trigObject.hasFilterLabel("hltEle23WPLooseGsfTrackIsoFilter")  ||
    	  trigObject.hasFilterLabel("hltEle27WPLooseGsfTrackIsoFilter")  ||
    	  trigObject.hasFilterLabel("hltEle27WPTightGsfTrackIsoFilter")  ||
    	  trigObject.hasFilterLabel("hltEle32WPLooseGsfTrackIsoFilter")  ||
    	  trigObject.hasFilterLabel("hltEle32WPTightGsfTrackIsoFilter")
    	  )
        {
          passSingleEleTagFilter = true;
        }
      std::cout << "debug ele 7 " << nElectrons << std::endl;
      //check Tag and Probe Filters
      if (trigObject.hasFilterLabel("hltEle25WP60Ele8TrackIsoFilter")) passTPOneTagFilter = true;
      if (trigObject.hasFilterLabel("hltEle25WP60SC4TrackIsoFilter")) passTPTwoTagFilter = true;
      if (trigObject.hasFilterLabel("hltEle25WP60Ele8Mass55Filter")) passTPOneProbeFilter = true;
      if (trigObject.hasFilterLabel("hltEle25WP60SC4Mass55Filter")) passTPTwoProbeFilter = true;
      std::cout << "debug ele 8 " << nElectrons << std::endl;
      //check all filters
      for ( int q=0; q<MAX_ElectronHLTFilters;q++)
      {
        if (trigObject.hasFilterLabel(eleHLTFilterNames[q].c_str())) ele_passHLTFilter[nElectrons][q] = true;
      }
    }
    */
    //std::cout << "debug ele 9 " << nElectrons << std::endl;
    ele_passSingleEleTagFilter[nElectrons] = passSingleEleTagFilter;
    ele_passTPOneTagFilter[nElectrons] = passTPOneTagFilter;
    ele_passTPTwoTagFilter[nElectrons] = passTPTwoTagFilter;
    ele_passTPOneProbeFilter[nElectrons] = passTPOneProbeFilter;
    ele_passTPTwoProbeFilter[nElectrons] = passTPTwoProbeFilter;


    //std::cout << "debug ele 10 " << nElectrons << std::endl;
    nElectrons++;
  }

  return true;
};

bool displacedJetMuon_ntupler::fillPileUp()
{
  for(const PileupSummaryInfo &pu : *puInfo)
  {
    BunchXing[nBunchXing] = pu.getBunchCrossing();
    nPU[nBunchXing] = pu.getPU_NumInteractions();
    nPUmean[nBunchXing] = pu.getTrueNumInteractions();
    nBunchXing++;
  }
  return true;
};


bool displacedJetMuon_ntupler::fillPhotons(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  noZS::EcalClusterLazyTools *lazyToolnoZS = new noZS::EcalClusterLazyTools(iEvent, iSetup, ebRecHitsToken_, eeRecHitsToken_);
  for (const reco::Photon &pho : *photons) {
    if (pho.pt() < 20) continue;
    std::vector<float> vCov = lazyToolnoZS->localCovariances( *(pho.superCluster()->seed()) );
    //-------------------------------------------------
    //default photon 4-mometum already vertex corrected
    //-------------------------------------------------
     phoE[nPhotons]   = pho.energy();
    phoPt[nPhotons]  = pho.pt();
    phoEta[nPhotons] = pho.eta(); //correct this for the vertex
    phoPhi[nPhotons] = pho.phi(); //correct this for the vertex
  
    phoSigmaIetaIeta[nPhotons] = pho.sigmaIetaIeta();
    phoFull5x5SigmaIetaIeta[nPhotons] = pho.full5x5_sigmaIetaIeta();
    phoR9[nPhotons] = pho.full5x5_r9();
    pho_HoverE[nPhotons] = pho.hadTowOverEm();
    pho_isConversion[nPhotons] = pho.hasConversionTracks();

    //------------------------------------------
    // Fill default miniAOD isolation quantities
    //------------------------------------------
    pho_pfIsoChargedHadronIso[nPhotons] = pho.chargedHadronIso();
    pho_pfIsoChargedHadronIsoWrongVtx[nPhotons] = pho.chargedHadronIsoWrongVtx();
    pho_pfIsoNeutralHadronIso[nPhotons] = pho.neutralHadronIso();
    pho_pfIsoPhotonIso[nPhotons] = pho.photonIso();
    pho_pfIsoModFrixione[nPhotons] = pho.getPflowIsolationVariables().modFrixione;
    pho_pfIsoSumPUPt[nPhotons] = pho.sumPUPt();
    //-----------------------
    // super cluster position and time
    //-----------------------
    pho_superClusterEnergy[nPhotons] = pho.superCluster()->energy();
    pho_superClusterRawEnergy[nPhotons] = pho.superCluster()->rawEnergy();
    pho_superClusterEta[nPhotons]    = pho.superCluster()->eta();
    pho_superClusterPhi[nPhotons]    = pho.superCluster()->phi();
    pho_superClusterX[nPhotons]      = pho.superCluster()->x();
    pho_superClusterY[nPhotons]      = pho.superCluster()->y();
    pho_superClusterZ[nPhotons]      = pho.superCluster()->z();
    pho_hasPixelSeed[nPhotons]       = pho.hasPixelSeed();

    pho_superClusterSeedX[nPhotons]      = pho.superCluster()->seed()->x();
    pho_superClusterSeedY[nPhotons]      = pho.superCluster()->seed()->y();
    pho_superClusterSeedZ[nPhotons]      = pho.superCluster()->seed()->z();

    pho_superClusterSeedE[nPhotons]      = pho.superCluster()->seed()->energy();

  
    //------------------------------------------------------
    //Compute PF isolation
    //absolute uncorrected isolations with footprint removal
    //------------------------------------------------------
    const float coneSizeDR = 0.3;
    const float dxyMax = 0.1;
    const float dzMax = 0.2;
    float chargedIsoSumAllVertices[MAX_NPV];
    for (int q=0;q<MAX_NPV;++q) chargedIsoSumAllVertices[q] = 0.0;
    float chargedIsoSum = 0;
    float chargedIsoSum_NewPV_NoTiming = 0;
    float chargedIsoSum_NewPV_Timing50_TrkVtx = 0;
    float chargedIsoSum_NewPV_Timing80_TrkVtx = 0;
    float chargedIsoSum_NewPV_Timing100_TrkVtx = 0;
    float chargedIsoSum_NewPV_Timing120_TrkVtx = 0;
    float chargedIsoSum_NewPV_Timing50_TrkPho = 0;
    float chargedIsoSum_NewPV_Timing80_TrkPho = 0;
    float chargedIsoSum_NewPV_Timing100_TrkPho = 0;
    float chargedIsoSum_NewPV_Timing120_TrkPho = 0;

    float neutralHadronIsoSum = 0;
    float photonIsoSum = 0;
    // First, find photon direction with respect to the good PV
    math::XYZVector photon_directionWrtVtx(pho.superCluster()->x() - myPV->x(),pho.superCluster()->y() - myPV->y(),pho.superCluster()->z() - myPV->z());
    //math::XYZVector photon_directionWrtVtx_GenMatch(pho.superCluster()->x() - myPV_GenMatch->x(),pho.superCluster()->y() - myPV_GenMatch->y(),pho.superCluster()->z() - myPV_GenMatch->z());

    // old PV, Loop over all PF candidates
    for (const reco::PFCandidate &candidate : *pfCands)
    {
      // Check if this candidate is within the isolation cone
      float dR=deltaR(photon_directionWrtVtx.Eta(),photon_directionWrtVtx.Phi(),
      candidate.eta(), candidate.phi());
      if( dR > coneSizeDR ) continue;

      // Check if this candidate is not in the footprint

      //bool inFootprint = false;
      //for (auto itr : pho.associatedPackedPFCandidates()) {
      //if ( &(*itr) == &candidate) {
      //inFootprint = true;
      //  }
      //}
      //if( inFootprint ) continue;

      // Find candidate type
      reco::PFCandidate::ParticleType thisCandidateType = reco::PFCandidate::X;

      // the neutral hadrons and charged hadrons can be of pdgId types
      // only 130 (K0L) and +-211 (pi+-) in packed candidates
      const int pdgId = candidate.pdgId();
      if( pdgId == 22 )
      thisCandidateType = reco::PFCandidate::gamma;
      else if( abs(pdgId) == 130) // PDG ID for K0L
      thisCandidateType = reco::PFCandidate::h0;
      else if( abs(pdgId) == 211) // PDG ID for pi+-
      thisCandidateType = reco::PFCandidate::h;


      // Increment the appropriate isolation sum
      if( thisCandidateType == reco::PFCandidate::h ){
        // for charged hadrons, additionally check consistency
        // with the PV
        float dxy = -999, dz = -999;

        //For the primary vertex
        dz = candidate.trackRef()->dz(myPV->position());
        dxy =candidate.trackRef()->dxy(myPV->position());
        if (fabs(dz) <= dzMax && fabs(dxy) <= dxyMax) {
          chargedIsoSum += candidate.pt();
        }

        //loop over all vertices
        for(int q = 0; q < nPVAll; q++){
          if(!(vertices->at(q).isValid() && !vertices->at(q).isFake())) continue;

          dz = candidate.trackRef()->dz(vertices->at(q).position());
          dxy =candidate.trackRef()->dxy(vertices->at(q).position());
          if (fabs(dz) > dzMax) continue;
          if(fabs(dxy) > dxyMax) continue;
          // The candidate is eligible, increment the isolation
          chargedIsoSumAllVertices[q] += candidate.pt();
        }
      }
      if( thisCandidateType == reco::PFCandidate::h0 )
      neutralHadronIsoSum += candidate.pt();
      if( thisCandidateType == reco::PFCandidate::gamma )
      photonIsoSum += candidate.pt();
    }


    //fill the proper variables
    for(int q = 0; q < nPVAll; q++) {
      pho_sumChargedHadronPtAllVertices[nPhotons][q] = chargedIsoSumAllVertices[q];
    }
    pho_sumChargedHadronPt[nPhotons] = chargedIsoSum;

    pho_sumChargedHadronPt_NewPV_NoTiming[nPhotons] = chargedIsoSum_NewPV_NoTiming;

    pho_sumChargedHadronPt_NewPV_Timing50_TrkVtx[nPhotons] = chargedIsoSum_NewPV_Timing50_TrkVtx;
    pho_sumChargedHadronPt_NewPV_Timing80_TrkVtx[nPhotons] = chargedIsoSum_NewPV_Timing80_TrkVtx;
    pho_sumChargedHadronPt_NewPV_Timing100_TrkVtx[nPhotons] = chargedIsoSum_NewPV_Timing100_TrkVtx;
    pho_sumChargedHadronPt_NewPV_Timing120_TrkVtx[nPhotons] = chargedIsoSum_NewPV_Timing120_TrkVtx;

    pho_sumChargedHadronPt_NewPV_Timing50_TrkPho[nPhotons] = chargedIsoSum_NewPV_Timing50_TrkPho;
    pho_sumChargedHadronPt_NewPV_Timing80_TrkPho[nPhotons] = chargedIsoSum_NewPV_Timing80_TrkPho;
    pho_sumChargedHadronPt_NewPV_Timing100_TrkPho[nPhotons] = chargedIsoSum_NewPV_Timing100_TrkPho;
    pho_sumChargedHadronPt_NewPV_Timing120_TrkPho[nPhotons] = chargedIsoSum_NewPV_Timing120_TrkPho;

    pho_sumNeutralHadronEt[nPhotons] = neutralHadronIsoSum;
    pho_sumPhotonEt[nPhotons] = photonIsoSum;


    //-------------------------------------------------
    //Compute Worst Isolation Looping over all vertices
    //-------------------------------------------------
    const double ptMin = 0.0;
    const float dRvetoBarrel = 0.0;
    const float dRvetoEndcap = 0.0;
    float dRveto = 0;
    if (pho.isEB()) dRveto = dRvetoBarrel;
    else dRveto = dRvetoEndcap;

    float worstIsolation = 999;
    std::vector<float> allIsolations;
    for(unsigned int ivtx=0; ivtx<vertices->size(); ++ivtx) {

      // Shift the photon according to the vertex
      reco::VertexRef vtx(vertices, ivtx);
      math::XYZVector photon_directionWrtVtx(pho.superCluster()->x() - vtx->x(),
      pho.superCluster()->y() - vtx->y(),
      pho.superCluster()->z() - vtx->z());

      float sum = 0;
      // Loop over all PF candidates
      for (const reco::PFCandidate &candidate : *pfCands) {

        //require that PFCandidate is a charged hadron
        const int pdgId = candidate.pdgId();
        if( abs(pdgId) != 211) continue;

        if (candidate.pt() < ptMin)
        continue;

        float dxy = -999, dz = -999;
        dz = candidate.trackRef()->dz(myPV->position());
        dxy =candidate.trackRef()->dxy(myPV->position());
        if( fabs(dxy) > dxyMax) continue;
        if ( fabs(dz) > dzMax) continue;

        float dR = deltaR(photon_directionWrtVtx.Eta(), photon_directionWrtVtx.Phi(),
        candidate.eta(),      candidate.phi());
        if(dR > coneSizeDR || dR < dRveto) continue;

        sum += candidate.pt();
      }

      allIsolations.push_back(sum);
    }

    if( allIsolations.size()>0 )
    worstIsolation = * std::max_element( allIsolations.begin(), allIsolations.end() );

    pho_sumWorstVertexChargedHadronPt[nPhotons] = worstIsolation;

    //-----------------------
    //Photon ID MVA variable
    //-----------------------
    //pho_IDMVA[nPhotons] = myPhotonMVA->mvaValue( pho,  *rhoAll, photonIsoSum, chargedIsoSum, worstIsolation,lazyToolnoZS, false);

    //---------------------
    //Regression
    //---------------------
    pho_RegressionE[nPhotons]            = pho.getCorrectedEnergy( pho.getCandidateP4type() );
    pho_RegressionEUncertainty[nPhotons] = pho.getCorrectedEnergyError( pho.getCandidateP4type() );


    //conversion matching for beamspot pointing
    const reco::Conversion *convmatch = 0;
    double drmin = std::numeric_limits<double>::max();
    //double leg conversions
    for (const reco::Conversion &conv : *conversions) {
      if (conv.refittedPairMomentum().rho()<10.) continue;
      if (!conv.conversionVertex().isValid()) continue;
      if (TMath::Prob(conv.conversionVertex().chi2(),  conv.conversionVertex().ndof())<1e-6) continue;

      math::XYZVector mom(conv.refittedPairMomentum());
      math::XYZPoint scpos(pho.superCluster()->position());
      math::XYZPoint cvtx(conv.conversionVertex().position());
      math::XYZVector cscvector = scpos - cvtx;

      double dr = reco::deltaR(mom,cscvector);

      if (dr<drmin && dr<0.1) {
        drmin = dr;
        convmatch = &conv;
      }
    }
    if (!convmatch) {
      drmin = std::numeric_limits<double>::max();
      //single leg conversions
      for (const reco::Conversion &conv : *singleLegConversions) {
        math::XYZVector mom(conv.tracksPin()[0]);
        math::XYZPoint scpos(pho.superCluster()->position());
        math::XYZPoint cvtx(conv.conversionVertex().position());
        math::XYZVector cscvector = scpos - cvtx;

        double dr = reco::deltaR(mom,cscvector);

        if (dr<drmin && dr<0.1) {
          drmin = dr;
          convmatch = &conv;
        }
      }
    }

    //matched conversion, compute conversion type
    //and extrapolation to beamline
    //FIXME Both of these additional two requirements are inconsistent and make the conversion
    //selection depend on poorly defined criteria, but we keep them for sync purposes
    //if (convmatch && pho.hasConversionTracks() && conversions->size()>0) {
    if (convmatch){// && pho.hasConversionTracks() && conversions->size()>0) {
      int ntracks = convmatch->nTracks();

      math::XYZVector mom(ntracks==2 ? convmatch->refittedPairMomentum() : convmatch->tracksPin()[0]);
      math::XYZPoint scpos(pho.superCluster()->position());
      math::XYZPoint cvtx(convmatch->conversionVertex().position());
      math::XYZVector cscvector = scpos - cvtx;

      double z = cvtx.z();
      double rho = cvtx.rho();

      int legtype = ntracks==2 ? 0 : 1;
      int dettype = pho.isEB() ? 0 : 1;
      int postype =0;

      if (pho.isEB()) {
        if (rho<15.) {
          postype = 0;
        }
        else if (rho>=15. && rho<60.) {
          postype = 1;
        }
        else {
          postype = 2;
        }
      }
      else {
        if (std::abs(z) < 50.) {
          postype = 0;
        }
        else if (std::abs(z) >= 50. && std::abs(z) < 100.) {
          postype = 1;
        }
        else {
          postype = 2;
        }
      }

      pho_convType[nPhotons] = legtype + 2*dettype + 4*postype;
      pho_convTrkZ[nPhotons] = cvtx.z() - ((cvtx.x()-beamSpot->x0())*mom.x()+(cvtx.y()-beamSpot->y0())*mom.y())/mom.rho() * mom.z()/mom.rho();
      pho_convTrkClusZ[nPhotons] = cvtx.z() - ((cvtx.x()-beamSpot->x0())*cscvector.x()+(cvtx.y()-beamSpot->y0())*cscvector.y())/cscvector.rho() * cscvector.z()/cscvector.rho();
    }

    nPhotons++;
  }

 
  delete lazyToolnoZS;
  return true;

};

double displacedJetMuon_ntupler::deltaPhi(double phi1, double phi2)
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

double displacedJetMuon_ntupler::deltaR(double eta1, double phi1, double eta2, double phi2)
{
  double dphi = deltaPhi(phi1,phi2);
  double deta = eta1 - eta2;
  return sqrt( dphi*dphi + deta*deta);
};

bool displacedJetMuon_ntupler::fillPVAll()
{
  nPVAll = std::min(int(vertices->size()),int(MAX_NPV));
  for (int ipv = 0; ipv < nPVAll; ++ipv)
  {
    const reco::Vertex &vtx = vertices->at(ipv);
    pvAllX[ipv] = vtx.x();
    pvAllY[ipv] = vtx.y();
    pvAllZ[ipv] = vtx.z();
  }

  double pvAllSumPtSqD[MAX_NPV];
  double pvAllSumPxD[MAX_NPV];
  double pvAllSumPyD[MAX_NPV];

  for (int ipv=0; ipv<nPVAll; ++ipv)
  {
    pvAllSumPtSqD[ipv] = 0.;
    pvAllSumPxD[ipv]   = 0.;
    pvAllSumPyD[ipv]   = 0.;
  }

  for (int ipv=0; ipv<nPVAll; ++ipv) {
    pvAllLogSumPtSq[ipv] = log(pvAllSumPtSqD[ipv]);
    pvAllSumPx[ipv] = pvAllSumPxD[ipv];
    pvAllSumPy[ipv] = pvAllSumPyD[ipv];
  }

  return true;
};

//define this as a plug-in
DEFINE_FWK_MODULE(displacedJetMuon_ntupler);
