// -*- C++ -*-
// Class:      displacedJetMuon_rechit_studies
/*
  Description: Base class for miniAOD analysis with CRAB
*/
//         Author:  Cristián Peña and Si Xie.
//         Created:  Thu, 26 March 2019 15:00:06 GMT

#include "displacedJetMuon_rechit_studies.h"
#include <fastjet/internal/base.hh>
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"

struct muonCscLayer{
  int id;
  int nSimHits;
  int nRecHits;
  int nWireDigis;
  int nStripDigis;
  bool Me1112;
  int match_gParticle_index;
  int strip_digi_occupancy;
  int wire_digi_occupancy;
};
struct largest_muonCscLayer
{
  inline bool operator() (const muonCscLayer& c1, const muonCscLayer& c2){return c1.nSimHits > c2.nSimHits;}
} my_largest_muonCscLayer_;

//------ Constructors and destructor ------//
displacedJetMuon_rechit_studies::displacedJetMuon_rechit_studies(const edm::ParameterSet& iConfig):
  //get inputs from config file
  isData_(iConfig.getParameter<bool> ("isData")),
  //isFourJet_(iConfig.getParameter<bool> ("isFourJet")),
  useGen_(iConfig.getParameter<bool> ("useGen")),
  isFastsim_(iConfig.getParameter<bool> ("isFastsim")),
  //isQCD_(iConfig.getParameter<bool> ("isQCD")),
  enableTriggerInfo_(iConfig.getParameter<bool> ("enableTriggerInfo")),
  enableCaloJet_(iConfig.getParameter<bool> ("enableCaloJet")),
  enableGenLLPInfo_(iConfig.getParameter<bool> ("enableGenLLPInfo")),
  enableEcalRechits_(iConfig.getParameter<bool> ("enableEcalRechits")),
  readGenVertexTime_(iConfig.getParameter<bool> ("readGenVertexTime")),
  triggerPathNamesFile_(iConfig.getParameter<string> ("triggerPathNamesFile")),
  eleHLTFilterNamesFile_(iConfig.getParameter<string> ("eleHLTFilterNamesFile")),
  muonHLTFilterNamesFile_(iConfig.getParameter<string> ("muonHLTFilterNamesFile")),
  photonHLTFilterNamesFile_(iConfig.getParameter<string> ("photonHLTFilterNamesFile")),
  verticesToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  tracksTag_(consumes<edm::View<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracks"))),
  cscSegmentInputToken_(consumes<CSCSegmentCollection>(edm::InputTag("cscSegments"))),
  cscRechitInputToken_(consumes<CSCRecHit2DCollection>(edm::InputTag("csc2DRecHits"))),
  dtSegmentInputToken_(consumes<DTRecSegment4DCollection>(edm::InputTag("dt4DSegments"))),
  dtRechitInputToken_(consumes<DTRecHitCollection>(edm::InputTag("dt1DRecHits"))),
  dtCosmicSegmentInputToken_(consumes<DTRecSegment4DCollection>(edm::InputTag("dt4DCosmicSegments"))),
  rpcRecHitInputToken_(consumes<RPCRecHitCollection>(edm::InputTag("rpcRecHits"))),
  MuonDTSimHitsToken_(consumes<vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("MuonDTSimHits"))),
  MuonCSCSimHitsToken_(consumes<vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("MuonCSCSimHits"))),
  MuonCSCComparatorDigiToken_(consumes<MuonDigiCollection<CSCDetId,CSCComparatorDigi>>(iConfig.getParameter<edm::InputTag>("MuonCSCComparatorDigi"))),
  MuonCSCStripDigiToken_(consumes<MuonDigiCollection<CSCDetId,CSCStripDigi>>(iConfig.getParameter<edm::InputTag>("MuonCSCStripDigi"))),
  MuonCSCWireDigiToken_(consumes<MuonDigiCollection<CSCDetId,CSCWireDigi>>(iConfig.getParameter<edm::InputTag>("MuonCSCWireDigi"))),
  MuonCSCStripDigiSimLinksToken_(consumes<edm::DetSetVector<StripDigiSimLink>>(iConfig.getParameter<edm::InputTag>("MuonCSCStripDigiSimLinks"))),
  MuonCSCWireDigiSimLinksToken_(consumes<edm::DetSetVector<StripDigiSimLink>>(iConfig.getParameter<edm::InputTag>("MuonCSCWireDigiSimLinks"))),
  MuonDTDigiToken_(consumes<MuonDigiCollection<DTLayerId,DTDigi>>(iConfig.getParameter<edm::InputTag>("MuonDTDigi"))),

  muonsToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronsToken_(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  tausToken_(consumes<reco::PFTauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  photonsToken_(consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  jetsCaloToken_(consumes<reco::CaloJetCollection>(iConfig.getParameter<edm::InputTag>("jetsCalo"))),
  // jetsPFToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jetsPF"))),
  jetsToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  jetsPuppiToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jetsPuppi"))),
  jetsAK8Token_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jetsAK8"))),
  PFCandsToken_(consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  PFClustersToken_(consumes<reco::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("pfClusters"))),
  //genParticlesToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"))),
  //genParticlesToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"))),
  genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  genJetsToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets"))),
  // triggerBitsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"))),
  hepMCToken_(consumes<edm::HepMCProduct>(iConfig.getParameter<edm::InputTag>("hepMC"))),
  //triggerObjectsToken_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
  //triggerPrescalesToken_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("triggerPrescales"))),
  genMetCaloToken_(consumes<reco::GenMETCollection>(iConfig.getParameter<edm::InputTag>("genMetsCalo"))),
  genMetTrueToken_(consumes<reco::GenMETCollection>(iConfig.getParameter<edm::InputTag>("genMetsTrue"))),
  metToken_(consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
//  metNoHFToken_(consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("metsNoHF"))),
  metPuppiToken_(consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("metsPuppi"))),
  // metFilterBitsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterBits"))),
  //hbheNoiseFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("hbheNoiseFilter"))),
  //hbheTightNoiseFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("hbheTightNoiseFilter"))),
  //hbheIsoNoiseFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("hbheIsoNoiseFilter"))),
  //badChargedCandidateFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("BadChargedCandidateFilter"))),
  //badMuonFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("BadMuonFilter"))),
//  lheRunInfoTag_(iConfig.getParameter<edm::InputTag>("lheInfo")),
//  lheRunInfoToken_(consumes<LHERunInfoProduct,edm::InRun>(lheRunInfoTag_)),
//  lheInfoToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheInfo"))),
  genInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genInfo"))),
  genLumiHeaderToken_(consumes<GenLumiInfoHeader,edm::InLumi>(edm::InputTag("generator",""))),
  puInfoToken_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("puInfo"))),
  //hcalNoiseInfoToken_(consumes<HcalNoiseSummary>(iConfig.getParameter<edm::InputTag>("hcalNoiseInfo"))),
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
  ebeeClustersToken_(consumes<vector<reco::CaloCluster> >(iConfig.getParameter<edm::InputTag>("ebeeClusters"))),
  esClustersToken_(consumes<vector<reco::CaloCluster> >(iConfig.getParameter<edm::InputTag>("esClusters"))),
  conversionsToken_(consumes<vector<reco::Conversion> >(iConfig.getParameter<edm::InputTag>("conversions"))),
  singleLegConversionsToken_(consumes<vector<reco::Conversion> >(iConfig.getParameter<edm::InputTag>("singleLegConversions"))),
  gedGsfElectronCoresToken_(consumes<vector<reco::GsfElectronCore> >(iConfig.getParameter<edm::InputTag>("gedGsfElectronCores"))),
  gedPhotonCoresToken_(consumes<vector<reco::PhotonCore> >(iConfig.getParameter<edm::InputTag>("gedPhotonCores"))),
  generalTrackToken_(consumes<std::vector<reco::Track>>(edm::InputTag("generalTracks")))
  //superClustersToken_(consumes<vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("superClusters"))),
  //  lostTracksToken_(consumes<vector<reco::PFCandidate> >(iConfig.getParameter<edm::InputTag>("lostTracks")))
  // mvaGeneralPurposeValuesMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaGeneralPurposeValuesMap"))),
  // mvaGeneralPurposeCategoriesMapToken_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaGeneralPurposeCategoriesMap"))),
  // mvaHZZValuesMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaHZZValuesMap"))),
  // mvaHZZCategoriesMapToken_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaHZZCategoriesMap")))
{
  //declare the TFileService for output
  edm::Service<TFileService> fs;

  //set up output tree
  displacedJetMuonTree = fs->make<TTree>("llp", "selected AOD information for llp analyses");
  //displacedJetMuonTree = new TTree("Jets", "selected AOD information");
  NEvents = fs->make<TH1F>("NEvents",";;NEvents;",1,-0.5,0.5);

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

  if(enableTriggerInfo_)
  {
    std::cout << "\n";
    std::cout << "****************** Trigger Paths Defined For Razor Ntuple ******************\n";
    for (int i = 0; i<NTriggersMAX; ++i)
    {
      if (triggerPathNames[i] != "") std::cout << "Trigger " << i << " " << triggerPathNames[i] << "\n";
    }
    std::cout << "****************************************************************************\n";
    std::cout << "\n";
  }
  if(readGenVertexTime_) genParticles_t0_Token_ = consumes<float>(iConfig.getParameter<edm::InputTag>("genParticles_t0"));

}

displacedJetMuon_rechit_studies::~displacedJetMuon_rechit_studies()
{
};

//***********************************************
//Enable output ntuple branches
//***********************************************

void displacedJetMuon_rechit_studies::setBranches()
{
  enableEventInfoBranches();
  // enablePVAllBranches();
  // enablePVTracksBranches();
  // enablePileUpBranches();
  // enableMuonBranches();
  // enableElectronBranches();
  // enableTauBranches();
  // enableIsoPFCandidateBranches();
  // enablePhotonBranches();
  enableMuonSystemBranches();
  // enableEcalRechitBranches();
  // enableJetBranches();
  // enableCaloJetBranches();
  // enableJetAK8Branches();
  // enableMetBranches();
  // enableTriggerBranches();
  // enableMCBranches();
  enableGenParticleBranches();
};

void displacedJetMuon_rechit_studies::enableEventInfoBranches()
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

void displacedJetMuon_rechit_studies::enablePVAllBranches()
{
  displacedJetMuonTree->Branch("nPVAll", &nPVAll,"nPVAll/I");
  displacedJetMuonTree->Branch("pvAllX", pvAllX,"pvAllX[nPVAll]/F");
  displacedJetMuonTree->Branch("pvAllY", pvAllY,"pvAllY[nPVAll]/F");
  displacedJetMuonTree->Branch("pvAllZ", pvAllZ,"pvAllZ[nPVAll]/F");
  displacedJetMuonTree->Branch("pvAllLogSumPtSq", pvAllLogSumPtSq,"pvAllLogSumPtSq[nPVAll]/F");
  displacedJetMuonTree->Branch("pvAllSumPx", pvAllSumPx,"pvAllSumPx[nPVAll]/F");
  displacedJetMuonTree->Branch("pvAllSumPy", pvAllSumPy,"pvAllSumPy[nPVAll]/F");
};

void displacedJetMuon_rechit_studies::enablePVTracksBranches()
{
  displacedJetMuonTree->Branch("nPVTracks", &nPVTracks,"nPVTracks/I");
  displacedJetMuonTree->Branch("pvTrackPt", pvTrackPt,"pvTrackPt[nPVTracks]/F");
  displacedJetMuonTree->Branch("pvTrackEta", pvTrackEta,"pvTrackEta[nPVTracks]/F");
  displacedJetMuonTree->Branch("pvTrackPhi", pvTrackPhi,"pvTrackPhi[nPVTracks]/F");
};

void displacedJetMuon_rechit_studies::enablePileUpBranches()
{
  displacedJetMuonTree->Branch("nBunchXing", &nBunchXing, "nBunchXing/I");
  displacedJetMuonTree->Branch("BunchXing", BunchXing, "BunchXing[nBunchXing]/I");
  displacedJetMuonTree->Branch("nPU", nPU, "nPU[nBunchXing]/I");
  displacedJetMuonTree->Branch("nPUmean", nPUmean, "nPUmean[nBunchXing]/F");
};

void displacedJetMuon_rechit_studies::enableMuonBranches()
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

void displacedJetMuon_rechit_studies::enableElectronBranches()
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

void displacedJetMuon_rechit_studies::enableTauBranches()
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

void displacedJetMuon_rechit_studies::enableIsoPFCandidateBranches()
{
  displacedJetMuonTree->Branch("nIsoPFCandidates", &nIsoPFCandidates, "nIsoPFCandidates/i");
  displacedJetMuonTree->Branch("isoPFCandidatePt", isoPFCandidatePt, "isoPFCandidatePt[nIsoPFCandidates]/F");
  displacedJetMuonTree->Branch("isoPFCandidateEta", isoPFCandidateEta, "isoPFCandidateEta[nIsoPFCandidates]/F");
  displacedJetMuonTree->Branch("isoPFCandidatePhi", isoPFCandidatePhi, "isoPFCandidatePhi[nIsoPFCandidates]/F");
  displacedJetMuonTree->Branch("isoPFCandidateIso04", isoPFCandidateIso04, "isoPFCandidateIso04[nIsoPFCandidates]/F");
  displacedJetMuonTree->Branch("isoPFCandidateD0", isoPFCandidateD0, "isoPFCandidateD0[nIsoPFCandidates]/F");
  displacedJetMuonTree->Branch("isoPFCandidatePdgId", isoPFCandidatePdgId, "isoPFCandidatePdgId[nIsoPFCandidates]/I");
};

void displacedJetMuon_rechit_studies::enablePhotonBranches()
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

void displacedJetMuon_rechit_studies::enableMuonSystemBranches()
{


    // displacedJetMuonTree->Branch("nCscSeg",&nCscSeg,"nCscSeg/I");
    // displacedJetMuonTree->Branch("cscSegPhi",cscSegPhi,"cscSegPhi[nCscSeg]/F");
    // displacedJetMuonTree->Branch("cscSegEta",cscSegEta,"cscSegEta[nCscSeg]/F");
    // displacedJetMuonTree->Branch("cscSegX",cscSegX,"cscSegX[nCscSeg]/F");
    // displacedJetMuonTree->Branch("cscSegY",cscSegY,"cscSegY[nCscSeg]/F");
    // displacedJetMuonTree->Branch("cscSegZ",cscSegZ,"cscSegZ[nCscSeg]/F");
    // displacedJetMuonTree->Branch("cscSegT",cscSegT,"cscSegT[nCscSeg]/F");
    // displacedJetMuonTree->Branch("cscSegChi2",cscSegChi2,"cscSegChi2[nCscSeg]/F");
    // displacedJetMuonTree->Branch("cscSegNRecHits",cscSegNRecHits,"cscSegNRecHits[nCscSeg]/I");

    displacedJetMuonTree->Branch("nCscStripDigis",&nCscStripDigis,"nCscStripDigis/I");
    displacedJetMuonTree->Branch("cscStripDigiStripN",cscStripDigiStripN,"cscStripDigiStripN[nCscStripDigis]/I");
    displacedJetMuonTree->Branch("cscStripDigiDetId",cscStripDigiDetId,"cscStripDigiDetId[nCscStripDigis]/I");
    displacedJetMuonTree->Branch("cscStripDigiEndcap",cscStripDigiEndcap,"cscStripDigiEndcap[nCscStripDigis]/I");
    displacedJetMuonTree->Branch("cscStripDigiStation",cscStripDigiStation,"cscStripDigiStation[nCscStripDigis]/I");
    displacedJetMuonTree->Branch("cscStripDigiRing",cscStripDigiRing,"cscStripDigiRing[nCscStripDigis]/I");
    displacedJetMuonTree->Branch("cscStripDigiChamber",cscStripDigiChamber,"cscStripDigiChamber[nCscStripDigis]/I");
    displacedJetMuonTree->Branch("cscStripDigiLayer",cscStripDigiLayer,"cscStripDigiLayer[nCscStripDigis]/I");

    displacedJetMuonTree->Branch("cscStripDigiADCCounts","std::vector<vector<int>>(nCscStripDigis)",&cscStripDigiADCCounts);

    displacedJetMuonTree->Branch("nCscWireDigis",&nCscWireDigis,"nCscWireDigis/I");
    displacedJetMuonTree->Branch("cscWireDigiWireGroupN",cscWireDigiWireGroupN,"cscWireDigiWireGroupN[nCscWireDigis]/I");
    displacedJetMuonTree->Branch("cscWireDigiWireGroupBx",cscWireDigiWireGroupBx,"cscWireDigiWireGroupBx[nCscWireDigis]/I");
    displacedJetMuonTree->Branch("cscWireDigiWireTimeBin",cscWireDigiWireTimeBin,"cscWireDigiWireTimeBin[nCscWireDigis]/I");
    displacedJetMuonTree->Branch("cscWireDigiWireGroupN",cscWireDigiWireGroupN,"cscWireDigiWireGroupN[nCscWireDigis]/I");
    displacedJetMuonTree->Branch("cscWireDigiWireTimeBinsOn","std::vector<vector<int>>(cscWireDigiWireTimeBinsOn)",&cscWireDigiWireTimeBinsOn);

    displacedJetMuonTree->Branch("cscWireDigiDetId",cscWireDigiDetId,"cscWireDigiDetId[nCscWireDigis]/I");
    displacedJetMuonTree->Branch("cscWireDigiEndcap",cscWireDigiEndcap,"cscWireDigiEndcap[nCscWireDigis]/I");
    displacedJetMuonTree->Branch("cscWireDigiStation",cscWireDigiStation,"cscWireDigiStation[nCscWireDigis]/I");
    displacedJetMuonTree->Branch("cscWireDigiRing",cscWireDigiRing,"cscWireDigiRing[nCscWireDigis]/I");
    displacedJetMuonTree->Branch("cscWireDigiChamber",cscWireDigiChamber,"cscWireDigiChamber[nCscWireDigis]/I");
    displacedJetMuonTree->Branch("cscWireDigiLayer",cscWireDigiLayer,"cscWireDigiLayer[nCscWireDigis]/I");
    // displacedJetMuonTree->Branch("nCscWireDigiSimLink",&nCscWireDigiSimLink,"nCscWireDigiSimLink/I");
    // displacedJetMuonTree->Branch("cscWireDigiSimLinkNWires",&cscWireDigiSimLinkNWires,"cscWireDigiSimLinkNWires/I");
    // displacedJetMuonTree->Branch("cscWireDigiSimLinkDetId",cscWireDigiSimLinkDetId,"cscWireDigiSimLinkDetId[nCscWireDigiSimLink]/I");
    // displacedJetMuonTree->Branch("cscWireDigiSimLinkChannel",cscWireDigiSimLinkChannel,"cscWireDigiSimLinkChannel[nCscWireDigiSimLink]/I");
    // displacedJetMuonTree->Branch("cscWireDigiSimLinkSimTrackId",cscWireDigiSimLinkSimTrackId,"cscWireDigiSimLinkSimTrackId[nCscWireDigiSimLink]/I");
    // displacedJetMuonTree->Branch("cscWireDigiSimLinkCFposition",cscWireDigiSimLinkCFposition,"cscWireDigiSimLinkCFposition[nCscWireDigiSimLink]/I");
    // displacedJetMuonTree->Branch("cscWireDigiSimLinkFraction",cscWireDigiSimLinkFraction,"cscWireDigiSimLinkFraction[nCscWireDigiSimLink]/F");
    // displacedJetMuonTree->Branch("cscWireDigiSimLinkOccupancy",cscWireDigiSimLinkOccupancy,"cscWireDigiSimLinkOccupancy[nCscWireDigiSimLink]/F");
    //
    // displacedJetMuonTree->Branch("nCscStripDigiSimLink",&nCscStripDigiSimLink,"nCscStripDigiSimLink/I");
    // displacedJetMuonTree->Branch("cscStripDigiSimLinkNStrips",&cscStripDigiSimLinkNStrips,"cscStripDigiSimLinkNStrips/I");
    // displacedJetMuonTree->Branch("cscStripDigiSimLinkDetId",cscStripDigiSimLinkDetId,"cscStripDigiSimLinkDetId[nCscStripDigiSimLink]/I");
    // displacedJetMuonTree->Branch("cscStripDigiSimLinkChannel",cscStripDigiSimLinkChannel,"cscStripDigiSimLinkChannel[nCscStripDigiSimLink]/I");
    // displacedJetMuonTree->Branch("cscStripDigiSimLinkSimTrackId",cscStripDigiSimLinkSimTrackId,"cscStripDigiSimLinkSimTrackId[nCscStripDigiSimLink]/I");
    // displacedJetMuonTree->Branch("cscStripDigiSimLinkCFposition",cscStripDigiSimLinkCFposition,"cscStripDigiSimLinkCFposition[nCscStripDigiSimLink]/I");
    // displacedJetMuonTree->Branch("cscStripDigiSimLinkFraction",cscStripDigiSimLinkFraction,"cscStripDigiSimLinkFraction[nCscStripDigiSimLink]/F");
    // displacedJetMuonTree->Branch("cscStripDigiSimLinkOccupancy",cscStripDigiSimLinkOccupancy,"cscStripDigiSimLinkOccupancy[nCscStripDigiSimLink]/F");
    //
    // displacedJetMuonTree->Branch("nCscDetLayer",&nCscDetLayer,"nCscDetLayer/i");
    // displacedJetMuonTree->Branch("cscDetLayer_nCscRecHits",cscDetLayer_nCscRecHits,"cscDetLayer_nCscRecHits[nCscDetLayer]/I");
    // displacedJetMuonTree->Branch("cscDetLayer_nCscSimHits",cscDetLayer_nCscSimHits,"cscDetLayer_nCscSimHits[nCscDetLayer]/I");
    // displacedJetMuonTree->Branch("cscDetLayer_nWireDigis",cscDetLayer_nWireDigis,"cscDetLayer_nWireDigis[nCscDetLayer]/I");
    // displacedJetMuonTree->Branch("cscDetLayer_nStripDigis",cscDetLayer_nStripDigis,"cscDetLayer_nStripDigis[nCscDetLayer]/I");
    // displacedJetMuonTree->Branch("cscDetLayer_inME1112",cscDetLayer_inME1112,"cscDetLayer_inME1112[nCscDetLayer]/O");
    // displacedJetMuonTree->Branch("cscDetLayer",cscDetLayer,"cscDetLayer[nCscDetLayer]/I");
    // displacedJetMuonTree->Branch("cscDetLayer_wireDigiOccupancy",cscDetLayer_wireDigiOccupancy,"cscDetLayer_wireDigiOccupancy[nCscDetLayer]/I");
    // displacedJetMuonTree->Branch("cscDetLayer_stripDigiOccupancy",cscDetLayer_stripDigiOccupancy,"cscDetLayer_stripDigiOccupancy[nCscDetLayer]/I");
    // displacedJetMuonTree->Branch("cscDetLayer_match_gParticle_index",cscDetLayer_match_gParticle_index,"cscDetLayer_match_gParticle_index[nCscDetLayer]/I");

    displacedJetMuonTree->Branch("nCscRechits",&nCscRechits,"nCscRechits/I");
    displacedJetMuonTree->Branch("cscRechitsPhi",cscRechitsPhi,"cscRechitsPhi[nCscRechits]/F");
    displacedJetMuonTree->Branch("cscRechitsEta",cscRechitsEta,"cscRechitsEta[nCscRechits]/F");
    displacedJetMuonTree->Branch("cscRechitsX",cscRechitsX,"cscRechitsX[nCscRechits]/F");
    displacedJetMuonTree->Branch("cscRechitsY",cscRechitsY,"cscRechitsY[nCscRechits]/F");
    displacedJetMuonTree->Branch("cscRechitsZ",cscRechitsZ,"cscRechitsZ[nCscRechits]/F");
    displacedJetMuonTree->Branch("cscRechitsE",cscRechitsE,"cscRechitsE[nCscRechits]/F");
    displacedJetMuonTree->Branch("cscRechitsTpeak",cscRechitsTpeak,"cscRechitsTpeak[nCscRechits]/F");
    displacedJetMuonTree->Branch("cscRechitsTwire",cscRechitsTwire,"cscRechitsTwire[nCscRechits]/F");
    displacedJetMuonTree->Branch("cscRechitsQuality",cscRechitsQuality,"cscRechitsQuality[nCscRechits]/I");
    displacedJetMuonTree->Branch("nMatchedSimHits",nMatchedSimHits,"nMatchedSimHits[nCscRechits]/I");
    displacedJetMuonTree->Branch("cscRechitsNStrips",cscRechitsNStrips,"cscRechitsNStrips[nCscRechits]/i");
    displacedJetMuonTree->Branch("cscRechitsHitWire",cscRechitsHitWire,"cscRechitsHitWire[nCscRechits]/I");
    displacedJetMuonTree->Branch("cscRechitsWGroupsBX",cscRechitsWGroupsBX,"cscRechitsWGroupsBX[nCscRechits]/I");
    displacedJetMuonTree->Branch("cscRechitsNWireGroups",cscRechitsNWireGroups,"cscRechitsNWireGroups[nCscRechits]/i");
    displacedJetMuonTree->Branch("cscRechitsNTimeBins",cscRechitsNTimeBins,"cscRechitsNTimeBins[nCscRechits]/i");
    displacedJetMuonTree->Branch("cscRechitsBadStrip",cscRechitsBadStrip,"cscRechitsBadStrip[nCscRechits]/I");
    displacedJetMuonTree->Branch("cscRechitsBadWireGroup",cscRechitsBadWireGroup,"cscRechitsBadWireGroup[nCscRechits]/I");

    displacedJetMuonTree->Branch("cscRechitsDetId",cscRechitsDetId,"cscRechitsDetId[nCscRechits]/I");
    displacedJetMuonTree->Branch("cscRechitsEndcap",cscRechitsEndcap,"cscRechitsEndcap[nCscRechits]/I");
    displacedJetMuonTree->Branch("cscRechitsStation",cscRechitsStation,"cscRechitsStation[nCscRechits]/I");
    displacedJetMuonTree->Branch("cscRechitsChamber",cscRechitsChamber,"cscRechitsChamber[nCscRechits]/I");
    displacedJetMuonTree->Branch("cscRechitsRing",cscRechitsRing,"cscRechitsRing[nCscRechits]/I");
    displacedJetMuonTree->Branch("cscRechitsLayer",cscRechitsLayer,"cscRechitsLayer[nCscRechits]/I");


    displacedJetMuonTree->Branch("nCscSimHits",&nCscSimHits,"nCscSimHits/I");
    displacedJetMuonTree->Branch("cscSimHitsTOF",cscSimHitsTOF,"cscSimHitsTOF[nCscSimHits]/F");
    displacedJetMuonTree->Branch("cscSimHitsEnergyLoss",cscSimHitsEnergyLoss,"cscSimHitsEnergyLoss[nCscSimHits]/F");
    displacedJetMuonTree->Branch("cscSimHitsParticleType",cscSimHitsParticleType,"cscSimHitsParticleType[nCscSimHits]/I");
    displacedJetMuonTree->Branch("cscSimHitsPabs",cscSimHitsPabs,"cscSimHitsPabs[nCscSimHits]/F");
    displacedJetMuonTree->Branch("cscSimHitsProcessType",cscSimHitsProcessType,"cscSimHitsProcessType[nCscSimHits]/I");
    displacedJetMuonTree->Branch("cscSimHitsX",cscSimHitsX,"cscSimHitsX[nCscSimHits]/F");
    displacedJetMuonTree->Branch("cscSimHitsY",cscSimHitsY,"cscSimHitsY[nCscSimHits]/F");
    displacedJetMuonTree->Branch("cscSimHitsZ",cscSimHitsZ,"cscSimHitsZ[nCscSimHits]/F");
    displacedJetMuonTree->Branch("cscSimHitsPhi",cscSimHitsPhi,"cscSimHitsPhi[nCscSimHits]/F");
    displacedJetMuonTree->Branch("cscSimHitsEta",cscSimHitsEta,"cscSimHitsEta[nCscSimHits]/F");
    displacedJetMuonTree->Branch("cscSimHitsDetId",cscSimHitsDetId,"cscSimHitsDetId[nCscSimHits]/I");
    displacedJetMuonTree->Branch("cscSimHitsEndcap",cscSimHitsEndcap,"cscSimHitsEndcap[nCscSimHits]/I");
    displacedJetMuonTree->Branch("cscSimHitsChamber",cscSimHitsChamber,"cscSimHitsChamber[nCscSimHits]/I");
    displacedJetMuonTree->Branch("cscSimHitsStation",cscSimHitsStation,"cscSimHitsStation[nCscSimHits]/I");
    displacedJetMuonTree->Branch("cscSimHitsRing",cscSimHitsRing,"cscSimHitsRing[nCscSimHits]/I");
    displacedJetMuonTree->Branch("cscSimHitsLayer",cscSimHitsLayer,"cscSimHitsLayer[nCscSimHits]/I");

    // Barrel

    // displacedJetMuonTree->Branch("nDtSeg",&nDtSeg,"nDtSeg/I");
    // displacedJetMuonTree->Branch("dtSegPhi",dtSegPhi,"dtSegPhi[nDtSeg]/F");
    // displacedJetMuonTree->Branch("dtSegEta",dtSegEta,"dtSegEta[nDtSeg]/F");
    // displacedJetMuonTree->Branch("dtSegX",dtSegX,"dtSegX[nDtSeg]/F");
    // displacedJetMuonTree->Branch("dtSegY",dtSegY,"dtSegY[nDtSeg]/F");
    // displacedJetMuonTree->Branch("dtSegZ",dtSegZ,"dtSegZ[nDtSeg]/F");
    // displacedJetMuonTree->Branch("dtSegStation",dtSegStation,"dtSegStation[nDtSeg]/I");
    // displacedJetMuonTree->Branch("dtSegWheel",dtSegWheel,"dtSegWheel[nDtSeg]/I");
    // displacedJetMuonTree->Branch("dtSegTime",dtSegTime,"dtSegTime[nDtSeg]/F");
    // displacedJetMuonTree->Branch("dtSegTimeError",dtSegTimeError,"dtSegTimeError[nDtSeg]/F");

    displacedJetMuonTree->Branch("nDtRechits",             &nDtRechits,             "nDtRechits/I");
    displacedJetMuonTree->Branch("dtRechitX",             dtRechitX,             "dtRechitX[nDtRechits]/F");
    displacedJetMuonTree->Branch("dtRechitY",             dtRechitY,             "dtRechitY[nDtRechits]/F");
    displacedJetMuonTree->Branch("dtRechitZ",             dtRechitZ,             "dtRechitZ[nDtRechits]/F");
    displacedJetMuonTree->Branch("dtRechitEta",             dtRechitEta,             "dtRechitEta[nDtRechits]/F");
    displacedJetMuonTree->Branch("dtRechitPhi",             dtRechitPhi,             "dtRechitPhi[nDtRechits]/F");
    displacedJetMuonTree->Branch("dtRechitTime",             dtRechitTime,             "dtRechitTime[nDtRechits]/F");
    displacedJetMuonTree->Branch("dtRechitWireId",             dtRechitWireId,             "dtRechitWireId[nDtRechits]/I");
    displacedJetMuonTree->Branch("dtRechitStation",             dtRechitStation,             "dtRechitStation[nDtRechits]/I");
    displacedJetMuonTree->Branch("dtRechitWheel",             dtRechitWheel,             "dtRechitWheel[nDtRechits]/I");
    displacedJetMuonTree->Branch("dtRechitSector",             dtRechitSector,             "dtRechitSector[nDtRechits]/I");
    displacedJetMuonTree->Branch("dtRechitSuperLayer",             dtRechitSuperLayer,             "dtRechitSuperLayer[nDtRechits]/I");
    displacedJetMuonTree->Branch("dtRechitLayer",             dtRechitLayer,             "dtRechitLayer[nDtRechits]/I");

    displacedJetMuonTree->Branch("nDtSimHits",&nDtSimHits,"nDtSimHits/I");
    displacedJetMuonTree->Branch("dtSimHitsTOF",dtSimHitsTOF,"dtSimHitsTOF[nDtSimHits]/F");
    displacedJetMuonTree->Branch("dtSimHitsEnergyLoss",dtSimHitsEnergyLoss,"dtSimHitsEnergyLoss[nDtSimHits]/F");
    displacedJetMuonTree->Branch("dtSimHitsParticleType",dtSimHitsParticleType,"dtSimHitsParticleType[nDtSimHits]/I");
    displacedJetMuonTree->Branch("dtSimHitsPabs",dtSimHitsPabs,"dtSimHitsPabs[nDtSimHits]/F");
    displacedJetMuonTree->Branch("dtSimHitsProcessType",dtSimHitsProcessType,"dtSimHitsProcessType[nDtSimHits]/I");
    displacedJetMuonTree->Branch("dtSimHitsX",dtSimHitsX,"dtSimHitsX[nDtSimHits]/F");
    displacedJetMuonTree->Branch("dtSimHitsY",dtSimHitsY,"dtSimHitsY[nDtSimHits]/F");
    displacedJetMuonTree->Branch("dtSimHitsZ",dtSimHitsZ,"dtSimHitsZ[nDtSimHits]/F");
    displacedJetMuonTree->Branch("dtSimHitsPhi",dtSimHitsPhi,"dtSimHitsPhi[nDtSimHits]/F");
    displacedJetMuonTree->Branch("dtSimHitsEta",dtSimHitsEta,"dtSimHitsEta[nDtSimHits]/F");
    displacedJetMuonTree->Branch("dtSimHitsDetId",dtSimHitsDetId,"dtSimHitsDetId[nDtSimHits]/I");
    displacedJetMuonTree->Branch("dtSimHitsStation",             dtSimHitsStation,             "dtSimHitsStation[nDtSimHits]/I");
    displacedJetMuonTree->Branch("dtSimHitsWheel",             dtSimHitsWheel,             "dtSimHitsWheel[nDtSimHits]/I");
    displacedJetMuonTree->Branch("dtSimHitsSector",             dtSimHitsSector,             "dtSimHitsSector[nDtSimHits]/I");
    displacedJetMuonTree->Branch("dtSimHitsSuperLayer",             dtSimHitsSuperLayer,             "dtSimHitsSuperLayer[nDtSimHits]/I");
    displacedJetMuonTree->Branch("dtSimHitsLayer",             dtSimHitsLayer,             "dtSimHitsLayer[nDtSimHits]/I");

    displacedJetMuonTree->Branch("nDtDigis",                &nDtDigis,                    "nDtDigis/I");
    displacedJetMuonTree->Branch("dtDigiTime",              dtDigiTime,                 "dtDigiTime[nDtDigis]/F");
    displacedJetMuonTree->Branch("dtDigiWireId",             dtDigiWireId,             "dtDigiWireId[nDtDigis]/I");
    displacedJetMuonTree->Branch("dtDigiNumber",             dtDigiNumber,             "dtDigiNumber[nDtDigis]/I");
    displacedJetMuonTree->Branch("dtDigiCountsTDC",             dtDigiCountsTDC,             "dtDigiCountsTDC[nDtDigis]/I");
    displacedJetMuonTree->Branch("dtDigiStation",             dtDigiStation,             "dtDigiStation[nDtDigis]/I");
    displacedJetMuonTree->Branch("dtDigiSector",             dtDigiSector,             "dtDigiSector[nDtDigis]/I");
    displacedJetMuonTree->Branch("dtDigiWheel",             dtDigiWheel,             "dtDigiWheel[nDtDigis]/I");
    displacedJetMuonTree->Branch("dtDigiLayer",             dtDigiLayer,             "dtDigiLayer[nDtDigis]/I");
    displacedJetMuonTree->Branch("dtDigiSuperLayer",             dtDigiSuperLayer,             "dtDigiSuperLayer[nDtDigis]/I");


};

void displacedJetMuon_rechit_studies::enableEcalRechitBranches()
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

};

void displacedJetMuon_rechit_studies::enableJetBranches()
{
  displacedJetMuonTree->Branch("nJets", &nJets,"nJets/I");
  displacedJetMuonTree->Branch("jetE", jetE,"jetE[nJets]/F");
  displacedJetMuonTree->Branch("jetPt", jetPt,"jetPt[nJets]/F");
  displacedJetMuonTree->Branch("jetEta", jetEta,"jetEta[nJets]/F");
  displacedJetMuonTree->Branch("jetEt", jetEt,"jetEt[nJets]/F");

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
  displacedJetMuonTree->Branch("jet_matched",jet_matched,"jet_matched[nJets]/O");//matched to gen particles



};

void displacedJetMuon_rechit_studies::enableCaloJetBranches()
{
  displacedJetMuonTree->Branch("nCaloJets", &nCaloJets,"nCaloJets/I");
  displacedJetMuonTree->Branch("calojetE", calojetE,"calojetE[nCaloJets]/F");
  displacedJetMuonTree->Branch("calojetEt", calojetEt,"calojetEt[nCaloJets]/F");
  displacedJetMuonTree->Branch("calojetPt", calojetPt,"calojetPt[nCaloJets]/F");
  displacedJetMuonTree->Branch("calojetEta", calojetEta,"calojetEta[nCaloJets]/F");
  displacedJetMuonTree->Branch("calojetPhi", calojetPhi,"calojetPhi[nCaloJets]/F");
  displacedJetMuonTree->Branch("calojetMass", calojetMass, "calojetMass[nCaloJets]/F");
  displacedJetMuonTree->Branch("calojetJetArea", calojetJetArea, "calojetJetArea[nCaloJets]/F");
  displacedJetMuonTree->Branch("calojetPileupE", calojetPileupE, "calojetPileupE[nCaloJets]/F");
  displacedJetMuonTree->Branch("calojetPileupId", calojetPileupId, "calojetPileupId[nCaloJets]/F");
  displacedJetMuonTree->Branch("calojetPileupIdFlag", calojetPileupIdFlag, "calojetPileupIdFlag[nCaloJets]/I");
  displacedJetMuonTree->Branch("calojetPassIDLoose", calojetPassIDLoose, "calojetPassIDLoose[nCaloJets]/O");
  displacedJetMuonTree->Branch("calojetPassIDTight", calojetPassIDTight, "calojetPassIDTight[nCaloJets]/O");
  displacedJetMuonTree->Branch("calojetNRechits", calojetNRechits,"calojetNRechits[nCaloJets]/I");
  displacedJetMuonTree->Branch("calojetRechitE", calojetRechitE,"calojetRechitE[nCaloJets]/F");
  displacedJetMuonTree->Branch("calojetRechitT", calojetRechitT,"calojetRechitT[nCaloJets]/F");
  displacedJetMuonTree->Branch("calojetRechitT_rms", calojetRechitT_rms,"calojetRechitT_rms[nCaloJets]/F");
  displacedJetMuonTree->Branch("calojetAlphaMax",calojetAlphaMax,"calojetAlphaMax[nCaloJets]/F");
  displacedJetMuonTree->Branch("calojetBetaMax",calojetBetaMax,"calojetBetaMax[nCaloJets]/F");
  displacedJetMuonTree->Branch("calojetGammaMax_ET",calojetGammaMax_ET,"calojetGammaMax_ET[nCaloJets]/F");
  displacedJetMuonTree->Branch("calojetGammaMax_EM",calojetGammaMax_EM,"calojetGammaMax_EM[nCaloJets]/F");
  displacedJetMuonTree->Branch("calojetGammaMax_Hadronic",calojetGammaMax_Hadronic,"calojetGammaMax_Hadronic[nCaloJets]/F");
  displacedJetMuonTree->Branch("calojetGammaMax",calojetGammaMax,"calojetGammaMax[nCaloJets]/F");
  displacedJetMuonTree->Branch("calojet_HadronicEnergyFraction", calojet_HadronicEnergyFraction,"calojet_HadronicEnergyFraction[nCaloJets]/F");
  displacedJetMuonTree->Branch("calojet_EMEnergyFraction", calojet_EMEnergyFraction,"calojet_EMEnergyFraction[nCaloJets]/F");
  displacedJetMuonTree->Branch("calojetPtAllTracks",calojetPtAllTracks,"calojetPtAllTracks[nCaloJets]/F");
  displacedJetMuonTree->Branch("calojetPtAllPVTracks",calojetPtAllPVTracks,"calojetPtAllPVTracks[nCaloJets]/F");
  displacedJetMuonTree->Branch("calojetMedianTheta2D",calojetMedianTheta2D,"calojetMedianTheta2D[nCaloJets]/F");
  displacedJetMuonTree->Branch("calojetMedianIP",calojetMedianIP,"calojetMedianIP[nCaloJets]/F");
  displacedJetMuonTree->Branch("calojetMinDeltaRAllTracks",calojetMinDeltaRAllTracks,"calojetMinDeltaRAllTracks[nCaloJets]/F");
  displacedJetMuonTree->Branch("calojetMinDeltaRPVTracks",calojetMinDeltaRPVTracks,"calojetMinDeltaRPVTracks[nCaloJets]/F");
};


void displacedJetMuon_rechit_studies::enableJetAK8Branches()
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

void displacedJetMuon_rechit_studies::enableMetBranches()
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
  displacedJetMuonTree->Branch("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, "Flag_CSCTightHaloFilter/O");
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
};

void displacedJetMuon_rechit_studies::enableTriggerBranches()
{
  nameHLT = new std::vector<std::string>; nameHLT->clear();
  displacedJetMuonTree->Branch("HLTDecision", &triggerDecision, ("HLTDecision[" + std::to_string(NTriggersMAX) +  "]/O").c_str());
  displacedJetMuonTree->Branch("HLTPrescale", &triggerHLTPrescale, ("HLTPrescale[" + std::to_string(NTriggersMAX) +  "]/I").c_str());
  //displacedJetMuonTree->Branch("HLTMR", &HLTMR, "HLTMR/F");
  //displacedJetMuonTree->Branch("HLTRSQ", &HLTRSQ, "HLTRSQ/F");
};

void displacedJetMuon_rechit_studies::enableMCBranches()
{
  displacedJetMuonTree->Branch("nGenJets", &nGenJets, "nGenJets/I");
  displacedJetMuonTree->Branch("genJetE", genJetE, "genJetE[nGenJets]/F");
  displacedJetMuonTree->Branch("genJetPt", genJetPt, "genJetPt[nGenJets]/F");
  displacedJetMuonTree->Branch("genJetEta", genJetEta, "genJetEta[nGenJets]/F");
  displacedJetMuonTree->Branch("genJetPhi", genJetPhi, "genJetPhi[nGenJets]/F");
  // displacedJetMuonTree->Branch("genMetPtCalo", &genMetPtCalo, "genMetPtCalo/F");
  // displacedJetMuonTree->Branch("genMetPhiCalo", &genMetPhiCalo, "genMetPhiCalo/F");
  // displacedJetMuonTree->Branch("genMetPtTrue", &genMetPtTrue, "genMetPtTrue/F");
  // displacedJetMuonTree->Branch("genMetPhiTrue", &genMetPhiTrue, "genMetPhiTrue/F");
  // displacedJetMuonTree->Branch("genVertexX", &genVertexX, "genVertexX/F");
  // displacedJetMuonTree->Branch("genVertexY", &genVertexY, "genVertexY/F");
  // displacedJetMuonTree->Branch("genVertexZ", &genVertexZ, "genVertexZ/F");
  // displacedJetMuonTree->Branch("genVertexT", &genVertexT, "genVertexT/F");
  displacedJetMuonTree->Branch("genWeight", &genWeight, "genWeight/F");
  displacedJetMuonTree->Branch("genSignalProcessID", &genSignalProcessID, "genSignalProcessID/i");
  // displacedJetMuonTree->Branch("genQScale", &genQScale, "genQScale/F");
  // displacedJetMuonTree->Branch("genAlphaQCD", &genAlphaQCD, "genAlphaQCD/F");
  // displacedJetMuonTree->Branch("genAlphaQED", &genAlphaQED, "genAlphaQED/F");
  // scaleWeights = new std::vector<float>; scaleWeights->clear();
  // pdfWeights = new std::vector<float>; pdfWeights->clear();
  // alphasWeights = new std::vector<float>; alphasWeights->clear();
  // if (isFastsim_) {
  //   displacedJetMuonTree->Branch("lheComments", "std::string",&lheComments);
  // }
  // displacedJetMuonTree->Branch("scaleWeights", "std::vector<float>",&scaleWeights);
  // displacedJetMuonTree->Branch("pdfWeights", "std::vector<float>",&pdfWeights);
  // displacedJetMuonTree->Branch("alphasWeights", "std::vector<float>",&alphasWeights);
};

void displacedJetMuon_rechit_studies::enableGenParticleBranches()
{
  displacedJetMuonTree->Branch("nGenParticle", &nGenParticle, "nGenParticle/I");
  displacedJetMuonTree->Branch("gParticleMotherId", gParticleMotherId, "gParticleMotherId[nGenParticle]/I");
  displacedJetMuonTree->Branch("gParticleMotherIndex", gParticleMotherIndex, "gParticleMotherIndex[nGenParticle]/I");
  displacedJetMuonTree->Branch("gParticleId", gParticleId, "gParticleId[nGenParticle]/I");
  displacedJetMuonTree->Branch("gParticleStatus", gParticleStatus, "gParticleStatus[nGenParticle]/I");
  displacedJetMuonTree->Branch("gParticleE", gParticleE, "gParticleE[nGenParticle]/F");
  displacedJetMuonTree->Branch("gParticlePt", gParticlePt, "gParticlePt[nGenParticle]/F");
  // displacedJetMuonTree->Branch("gParticlePx", gParticlePx, "gParticlePx[nGenParticle]/F");
  // displacedJetMuonTree->Branch("gParticlePy", gParticlePy, "gParticlePy[nGenParticle]/F");
  // displacedJetMuonTree->Branch("gParticlePz", gParticlePz, "gParticlePz[nGenParticle]/F");
  displacedJetMuonTree->Branch("gParticleEta", gParticleEta, "gParticleEta[nGenParticle]/F");
  displacedJetMuonTree->Branch("gParticlePhi", gParticlePhi, "gParticlePhi[nGenParticle]/F");
  displacedJetMuonTree->Branch("gParticleProdVertexX", gParticleProdVertexX, "gParticleProdVertexX[nGenParticle]/F");
  displacedJetMuonTree->Branch("gParticleProdVertexY", gParticleProdVertexY, "gParticleProdVertexY[nGenParticle]/F");
  displacedJetMuonTree->Branch("gParticleProdVertexZ", gParticleProdVertexZ, "gParticleProdVertexZ[nGenParticle]/F");

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

   displacedJetMuonTree->Branch("gLLP_travel_time", gLLP_travel_time, "gLLP_travel_time[2]/F");
   displacedJetMuonTree->Branch("gLLP_gamma", gLLP_gamma, "gLLP_gamma[2]/F");
   displacedJetMuonTree->Branch("gLLP_ctau", gLLP_ctau, "gLLP_ctau[2]/F");
   displacedJetMuonTree->Branch("gLLP_Eem", gLLP_Eem, "gLLP_Eem[2]/F");
   displacedJetMuonTree->Branch("gLLP_Ehad", gLLP_Ehad, "gLLP_Ehad[2]/F");
 //   displacedJetMuonTree->Branch("gLLP_daughter_travel_time", gLLP_daughter_travel_time, "gLLP_daughter_travel_time[4]/F");
 //   displacedJetMuonTree->Branch("gLLP_daughter_pt", gLLP_daughter_pt, "gLLP_daughter_pt[4]/F");
 //   displacedJetMuonTree->Branch("gLLP_daughter_eta", gLLP_daughter_eta, "gLLP_daughter_eta[4]/F");
 //   displacedJetMuonTree->Branch("gLLP_daughter_phi", gLLP_daughter_phi, "gLLP_daughter_phi[4]/F");
 //   displacedJetMuonTree->Branch("gLLP_daughter_eta_ecalcorr", gLLP_daughter_eta_ecalcorr, "gLLP_daughter_eta_ecalcorr[4]/F");
 //   displacedJetMuonTree->Branch("gLLP_daughter_phi_ecalcorr", gLLP_daughter_phi_ecalcorr, "gLLP_daughter_phi_ecalcorr[4]/F");
 //
 //   displacedJetMuonTree->Branch("gLLP_daughter_e", gLLP_daughter_e, "gLLP_daughter_e[4]/F");
 //   displacedJetMuonTree->Branch("photon_travel_time", photon_travel_time, "photon_travel_time[4]/F");
 //   displacedJetMuonTree->Branch("gen_time", gen_time, "gen_time[4]/F");
 //   displacedJetMuonTree->Branch("gen_time_pv", gen_time_pv, "gen_time_pv[4]/F");
 //
 //   displacedJetMuonTree->Branch("gLLP_min_delta_r_match_calojet", gLLP_min_delta_r_match_calojet, "gLLP_min_delta_r_match_calojet[4]/F");
 // displacedJetMuonTree->Branch("gLLP_daughter_match_calojet_index", gLLP_daughter_match_calojet_index, "gLLP_daughter_match_calojet_index[4]/i");
 //   displacedJetMuonTree->Branch("gLLP_daughter_match_jet_index", gLLP_daughter_match_jet_index, "gLLP_daughter_match_jet_index[4]/i");
 //   displacedJetMuonTree->Branch("gLLP_min_delta_r_match_jet", gLLP_min_delta_r_match_jet, "gLLP_min_delta_r_match_jet[4]/F");

};

//------ Load the miniAOD objects and reset tree variables for each event ------//
void displacedJetMuon_rechit_studies::loadEvent(const edm::Event& iEvent)//load all miniAOD objects for the current event
{
  // iEvent.getByToken(triggerBitsToken_, triggerBits);
  iEvent.getByToken(hepMCToken_, hepMC);
  // iEvent.getByToken(metFilterBitsToken_, metFilterBits);
  iEvent.getByToken(verticesToken_, vertices);
  iEvent.getByToken(cscSegmentInputToken_,cscSegments);
  iEvent.getByToken(cscRechitInputToken_,cscRechits);
  iEvent.getByToken(dtSegmentInputToken_,dtSegments);
  iEvent.getByToken(dtRechitInputToken_,dtRechits);
  iEvent.getByToken(dtCosmicSegmentInputToken_,dtCosmicSegments);
  iEvent.getByToken(rpcRecHitInputToken_,rpcRecHits);

  iEvent.getByToken(MuonCSCStripDigiToken_, MuonCSCStripDigi);
  iEvent.getByToken(MuonCSCWireDigiToken_, MuonCSCWireDigi);
  iEvent.getByToken(MuonDTDigiToken_, MuonDTDigi);

  if (!isData)
  {
    iEvent.getByToken(MuonCSCStripDigiSimLinksToken_, MuonCSCStripDigiSimLinks);
    iEvent.getByToken(MuonCSCWireDigiSimLinksToken_, MuonCSCWireDigiSimLinks);
  }

  iEvent.getByToken(tracksTag_,tracks);
  iEvent.getByToken(PFCandsToken_, pfCands);
  iEvent.getByToken(PFClustersToken_, pfClusters);
  iEvent.getByToken(muonsToken_, muons);
  iEvent.getByToken(electronsToken_, electrons);
  iEvent.getByToken(photonsToken_, photons);
  iEvent.getByToken(tausToken_, taus);
  iEvent.getByToken(jetsCaloToken_, jetsCalo);
  // iEvent.getByToken(jetsPFToken_, jetsPF);
  iEvent.getByToken(jetsToken_, jets);
  iEvent.getByToken(jetsPuppiToken_, jetsPuppi);
  iEvent.getByToken(jetsAK8Token_, jetsAK8);
  iEvent.getByToken(genMetCaloToken_, genMetsCalo);
  iEvent.getByToken(genMetTrueToken_, genMetsTrue);
  iEvent.getByToken(metToken_, mets);
  //iEvent.getByToken(metNoHFToken_, metsNoHF);
  iEvent.getByToken(metPuppiToken_, metsPuppi);
//  iEvent.getByToken(hcalNoiseInfoToken_,hcalNoiseInfo);
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
//  iEvent.getByToken(superClustersToken_,superClusters);
//  iEvent.getByToken(lostTracksToken_,lostTracks);
//  iEvent.getByToken(hbheNoiseFilterToken_, hbheNoiseFilter);
//  iEvent.getByToken(hbheTightNoiseFilterToken_, hbheTightNoiseFilter);
//  iEvent.getByToken(hbheIsoNoiseFilterToken_, hbheIsoNoiseFilter);
  //iEvent.getByToken(badChargedCandidateFilterToken_, badChargedCandidateFilter);
  //iEvent.getByToken(badMuonFilterToken_, badMuonFilter);
  if(readGenVertexTime_) iEvent.getByToken(genParticles_t0_Token_,genParticles_t0);
  if (useGen_) {
    iEvent.getByToken(genParticlesToken_,genParticles);
    iEvent.getByToken(genJetsToken_,genJets);
    iEvent.getByToken(MuonDTSimHitsToken_, MuonDTSimHits);
    iEvent.getByToken(MuonCSCSimHitsToken_, MuonCSCSimHits);
    //for Spring16 fastsim, this has been changed and removed
//    if (!isFastsim_) iEvent.getByToken(lheInfoToken_, lheInfo);

    iEvent.getByToken(genInfoToken_,genInfo);
    iEvent.getByToken(puInfoToken_,puInfo);
  }


}

//called by the loadEvent() method
void displacedJetMuon_rechit_studies::resetBranches()
{
    //reset tree variables
    resetEventInfoBranches();
    resetPVAllBranches();
    resetPVTracksBranches();
    resetPileUpBranches();
    resetMuonBranches();
    resetElectronBranches();
    resetTauBranches();
    resetPhotonBranches();
    resetJetBranches();
    resetCaloJetBranches();
    resetMuonSystemBranches();
    resetMetBranches();
    resetGenParticleBranches();
    resetMCBranches();
    resetTriggerBranches();
    resetEcalRechitBranches();

};

void displacedJetMuon_rechit_studies::resetEventInfoBranches()
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

void displacedJetMuon_rechit_studies::resetPVAllBranches()
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

void displacedJetMuon_rechit_studies::resetPVTracksBranches()
{
  nPVTracks = 0;
  for(int i = 0; i < OBJECTARRAYSIZE; i++)
  {
    pvTrackPt[i]  = -999.;
    pvTrackEta[i] = -999.;
    pvTrackPhi[i] = -999.;
  }
};

void displacedJetMuon_rechit_studies::resetPileUpBranches()
{
  nBunchXing = 0;
  for(int i = 0; i < MAX_NBX; i++)
  {
    BunchXing[i] = -999;
    nPU[i]       = -999;
    nPUmean[i]   = -999.;
  }
};

void displacedJetMuon_rechit_studies::resetMuonBranches()
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

void displacedJetMuon_rechit_studies::resetElectronBranches()
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

void displacedJetMuon_rechit_studies::resetTauBranches()
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

void displacedJetMuon_rechit_studies::resetPhotonBranches()
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

void displacedJetMuon_rechit_studies::resetMuonSystemBranches()
{
    nCscSeg = 0;
    nCscRechits = 0;
    nCscSimHits = 0;
    nCscStripDigis = 0;
    nCscWireDigis = 0;
    nDtDigis = 0;
    // nCscDetLayer = 0;
    cscStripDigiADCCounts.clear();
    cscWireDigiWireTimeBinsOn.clear();
    // nCscWireDigiSimLink = 0;
    // nCscStripDigiSimLink = 0;
    // cscWireDigiSimLinkNWires = 0;
    // cscStripDigiSimLinkNStrips = 0;
    for ( int i = 0; i < OBJECTARRAYSIZE; i++)
    {
      // cscSegPhi[i] = 0.0;
      // cscSegEta[i] = 0.0;
      // cscSegX[i] = 0.0;
      // cscSegY[i] = 0.0;
      // cscSegZ[i] = 0.0;
      // cscSegDirectionX[i] = 0.0;
      // cscSegDirectionY[i] = 0.0;
      // cscSegDirectionZ[i] = 0.0;
      // cscSegT[i] = 0.0;
      // cscSegChi2[i] = 0.0;
      // cscSegNRechits[i] = 0;

      // cscDetLayer_nCscRechits[i] = 0;
      // cscDetLayer_nCscSimHits[i] = 0;
      // cscDetLayer_nWireDigis[i] = 0;
      // cscDetLayer_nStripDigis[i] = 0;
      // cscDetLayer_wireDigiOccupancy[i] = 0;
      // cscDetLayer_stripDigiOccupancy[i] = 0;
      // cscDetLayer_match_gParticle_index[i] = 0;
      // cscDetLayer_inME1112[i] = false;
      // cscDetLayer[i] = 0;
      cscRechitsPhi[i] = 0.0;
      cscRechitsEta[i] = 0.0;
      cscRechitsX[i] = 0.0;
      cscRechitsY[i] = 0.0;
      cscRechitsZ[i] = 0.0;
      cscRechitsE[i] = 0.0;
      cscRechitsTwire[i] = 0.0;
      cscRechitsTpeak[i] = 0.0;
      cscRechitsQuality[i] = 0.0;
      nMatchedSimHits[i] = 0;
      cscRechitsNStrips[i] = 0;
      cscRechitsHitWire[i] = 0;
      cscRechitsWGroupsBX[i] = 0;
      cscRechitsNWireGroups[i] = 0;
      cscRechitsBadStrip[i] = 0;
      cscRechitsBadWireGroup[i] = 0;
      cscRechitsNTimeBins[i] = 0;
      cscRechitsDetId[i] = 0;
      cscRechitsEndcap[i] = 0;
      cscRechitsStation[i] = 0;
      cscRechitsRing[i] = 0;
      cscRechitsChamber[i] = 0;
      cscRechitsLayer[i] = 0;

      cscSimHitsTOF[i] = 0.0;
      cscSimHitsEnergyLoss[i] = 0.0;
      cscSimHitsParticleType[i] = -999;
      cscSimHitsPabs[i] = 0.0;
      cscSimHitsProcessType[i] = -999;
      cscSimHitsX[i] = 0.0;
      cscSimHitsY[i] = 0.0;
      cscSimHitsZ[i] = 0.0;
      cscSimHitsPhi[i] = 0.0;
      cscSimHitsEta[i] = 0.0;
      // cscSimHits_match_cscRecHits_index[i] = -999;
      // cscSimHits_match_cscRecHits_minDeltaR[i] = -999.;
      // cscSimHits_match_gParticle_index[i] = -999;
      // cscSimHits_match_gParticle_minDeltaR[i] = -999.;
      cscSimHitsDetId[i] = 0;
      cscSimHitsEndcap[i] = 0;
      cscSimHitsStation[i] = 0;
      cscSimHitsChamber[i] = 0;
      cscSimHitsRing[i] = 0;
      cscSimHitsLayer[i] = 0;

      cscWireDigiDetId[i] = -999;
      cscWireDigiEndcap[i] = -999;
      cscWireDigiStation[i] = -999;
      cscWireDigiRing[i] = -999;
      cscWireDigiChamber[i] = -999;
      cscWireDigiLayer[i] = -999;
      cscWireDigiWireGroupN[i] = -999;
      cscWireDigiWireGroupBx[i] = -999;
      cscWireDigiWireTimeBin[i] = -999;



      cscStripDigiDetId[i] = -999;
      cscStripDigiStripN[i] = -999;
      cscStripDigiEndcap[i] = -999;
      cscStripDigiStation[i] = -999;
      cscStripDigiRing[i] = -999;
      cscStripDigiChamber[i] = -999;
      cscStripDigiLayer[i] = -999;

      // cscWireDigiSimLinkChannel[i] = -999;
      // cscWireDigiSimLinkCFposition[i] = -999;
      // cscWireDigiSimLinkFraction[i] = -999.;
      // cscWireDigiSimLinkDetId[i] = -999;
      // cscWireDigiSimLinkSimTrackId[i] = -999;
      // cscWireDigiSimLinkOccupancy[i] = -999;
      //
      // cscStripDigiSimLinkChannel[i] = -999;
      // cscStripDigiSimLinkCFposition[i] = -999;
      // cscStripDigiSimLinkFraction[i] = -999.;
      // cscStripDigiSimLinkDetId[i] = -999;
      // cscStripDigiSimLinkSimTrackId[i] = -999;
      // cscStripDigiSimLinkOccupancy[i] = -999;


    }

    nDtRechits = 0;
    for ( int i = 0; i < OBJECTARRAYSIZE; i++)
    {
      dtRechitX[i] = 0.0;
      dtRechitY[i] = 0.0;
      dtRechitZ[i] = 0.0;
      dtRechitEta[i] = 0.0;
      dtRechitPhi[i] = 0.0;
      dtRechitTime[i] = 0.0;
      dtRechitWireId[i] = 0;
      dtRechitStation[i] = 0;
      dtRechitWheel[i] = 0;
      dtRechitSector[i] = 0;
      dtRechitSuperLayer[i] = 0;
      dtRechitLayer[i] = 0;

      dtSimHitsTOF[i] = 0.0;
      dtSimHitsEnergyLoss[i] = 0.0;
      dtSimHitsParticleType[i] = -999;
      dtSimHitsPabs[i] = 0.0;
      dtSimHitsProcessType[i] = -999;
      dtSimHitsX[i] = 0.0;
      dtSimHitsY[i] = 0.0;
      dtSimHitsZ[i] = 0.0;
      dtSimHitsPhi[i] = 0.0;
      dtSimHitsEta[i] = 0.0;
      dtSimHitsDetId[i] = 0;
      dtSimHitsStation[i] = 0;
      dtSimHitsWheel[i] = 0;
      dtSimHitsSector[i] = 0;
      dtSimHitsSuperLayer[i] = 0;
      dtSimHitsLayer[i] = 0;



      dtDigiTime[i] = 0.0;
      dtDigiWireId[i] = 0;
      dtDigiNumber[i] = 0;
      dtDigiCountsTDC[i] = 0;
      dtDigiStation[i] = 0;
      dtDigiSector[i] = 0;
      dtDigiWheel[i] = 0;
      dtDigiLayer[i] = 0;
      dtDigiSuperLayer[i] = 0;
    }

    return;
};
void displacedJetMuon_rechit_studies::resetJetBranches()
{
  nJets = 0;
  for ( int i = 0; i < OBJECTARRAYSIZE; i++)
  {
    jetE[i] = 0.0;
    jetEt[i] = 0.0;
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
    jet_matched[i] = false;

  }
  return;
};
void displacedJetMuon_rechit_studies::resetEcalRechitBranches()
{
  nRechits = 0;
  for ( int i = 0; i < RECHITARRAYSIZE; i++)
  {

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
  return;
};
void displacedJetMuon_rechit_studies::resetCaloJetBranches()
{
  nCaloJets = 0;
  for ( int i = 0; i < OBJECTARRAYSIZE; i++)
  {
    calojetE[i] = 0.0;
    calojetEt[i] = 0.0;
    calojetPt[i] = 0.0;
    calojetEta[i] = 0.0;
    calojetPhi[i] = 0.0;
    calojetMass[i] =  -99.0;
    calojetJetArea[i] = -99.0;
    calojetPileupE[i] = -99.0;
    calojetPileupId[i] = -99.0;
    calojetPileupIdFlag[i] = -1;
    calojetPassIDLoose[i] = false;
    calojetPassIDTight[i] = false;

    calojetNRechits[i] = 0;
    calojetRechitE[i] = 0.0;
    calojetRechitT[i] = 0.0;
    calojetRechitT_rms[i] = 0.0;

    calojetAlphaMax[i] = -99.0;
    calojetBetaMax[i] = -99.0;
    calojetGammaMax[i] = -99.0;
    calojetGammaMax_EM[i] = -99.0;
    calojetGammaMax_Hadronic[i] = -99.0;
    calojetGammaMax_ET[i] = -99.0;
    calojetPtAllTracks[i] = -99.0;
    calojetPtAllPVTracks[i] = -99.0;
    calojetMedianTheta2D[i] = -99.0;
    calojetMedianIP[i] = -99.0;
    calojetMinDeltaRAllTracks[i] =-99.0;
    calojetMinDeltaRPVTracks[i] = -99.0;

    calojet_HadronicEnergyFraction[i] = -666.;
    calojet_EMEnergyFraction[i] = -666.;
  }
  return;
};


void displacedJetMuon_rechit_studies::resetMetBranches()
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
  Flag_CSCTightHaloFilter = false;
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
void displacedJetMuon_rechit_studies::resetGenParticleBranches()
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
    gLLP_travel_time[i] = -666.;
    gLLP_gamma[i] = -666.;
    gLLP_ctau[i] = -666.;
    gLLP_Eem[i] = 0.0;
    gLLP_Ehad[i] = 0.0;
  }

  for ( int i = 0; i < LLP_DAUGHTER_ARRAY_SIZE; i++ )
  {
    gLLP_daughter_pt[i] = -666.;
    gLLP_daughter_eta[i] = -666.;
    gLLP_daughter_phi[i] = -666.;
    gLLP_daughter_eta_ecalcorr[i] = -666.;
    gLLP_daughter_phi_ecalcorr[i] = -666.;
    gLLP_daughter_e[i] = -666.;
    gLLP_daughter_travel_time[i] = -666.;
    gen_time[i] = -666.;
    gen_time_pv[i] = -666.;
    photon_travel_time[i] = -666.;
    photon_travel_time_pv[i] = -666.;
    gLLP_daughter_match_calojet_index[i] = 666;
    gLLP_daughter_match_jet_index[i] = 666;
    gLLP_min_delta_r_match_calojet[i] = -666.;
    gLLP_min_delta_r_match_jet[i] = -666.;



}
  return;
};

void displacedJetMuon_rechit_studies::resetMCBranches()
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

  return;
};

void displacedJetMuon_rechit_studies::resetTriggerBranches()
{
  for( int i = 0; i < NTriggersMAX; i++ )
  {
    triggerDecision[i] = false;
    triggerHLTPrescale[i] = 0;
  }
};
//------ Method called for each run ------//

void displacedJetMuon_rechit_studies::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {


}


//------ Method called for each lumi block ------//
void displacedJetMuon_rechit_studies::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {

}


//------ Method called for each event ------//
//------ Method called for each event ------//

void displacedJetMuon_rechit_studies::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;

  //initialize
  loadEvent(iEvent); //loads objects and resets tree branches
  NEvents->Fill(0); //increment event count
  //resetting output tree branches

  resetBranches();

  fillEventInfo(iEvent);
  // fillJets(iSetup);
  // fillMuons(iEvent);
  // fillElectrons(iEvent);
  // fillMet(iEvent);
  // fillCaloJets(iSetup);
  if (!isData) {
    // fillPileUp();
    // fillMC();
    fillGenParticles();
  }
  fillMuonSystem(iEvent, iSetup);



  // if ( enableTriggerInfo_ ) fillTrigger( iEvent );

  displacedJetMuonTree->Fill();

};


//------ Method called once each job just before starting event loop ------//
void displacedJetMuon_rechit_studies::beginJob()
{
  setBranches();
}

//------ Method called once each job just after ending the event loop ------//
void displacedJetMuon_rechit_studies::endJob(){};



//------ Methods to fill tree variables ------//
bool displacedJetMuon_rechit_studies::fillEventInfo(const edm::Event& iEvent)
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


bool displacedJetMuon_rechit_studies::fillMuonSystem(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::ESHandle<CSCGeometry> cscG;
  edm::ESHandle<DTGeometry> dtG;
  edm::ESHandle<RPCGeometry> rpcG;

  iSetup.get<MuonGeometryRecord>().get(cscG);
  iSetup.get<MuonGeometryRecord>().get(dtG);
  iSetup.get<MuonGeometryRecord>().get(rpcG);
  //*****************
  //** RECHITS
  //*****************
  for (const CSCRecHit2D cscRechit : *cscRechits) {
    LocalPoint  cscRecHitLocalPosition       = cscRechit.localPosition();
    CSCDetId cscdetid = cscRechit.cscDetId();
    cscRechitsDetId[nCscRechits] = CSCDetId::rawIdMaker(CSCDetId::endcap(cscdetid), CSCDetId::station(cscdetid), CSCDetId::ring(cscdetid), CSCDetId::chamber(cscdetid), CSCDetId::layer(cscdetid));
    const CSCChamber* cscchamber = cscG->chamber(cscdetid);
    if (cscchamber) {
      GlobalPoint globalPosition = cscchamber->toGlobal(cscRecHitLocalPosition);
      cscRechitsX[nCscRechits] = globalPosition.x();
      cscRechitsY[nCscRechits] = globalPosition.y();
      cscRechitsZ[nCscRechits] = globalPosition.z();
      cscRechitsPhi[nCscRechits] = globalPosition.phi();
      cscRechitsEta[nCscRechits] = globalPosition.eta();
      cscRechitsE[nCscRechits] = cscRechit.energyDepositedInLayer();//not saved
      cscRechitsTpeak[nCscRechits] = cscRechit.tpeak();
      cscRechitsTwire[nCscRechits] = cscRechit.wireTime();
      cscRechitsQuality[nCscRechits] = cscRechit.quality();
      cscRechitsNStrips[nCscRechits] = cscRechit.nStrips();
      cscRechitsHitWire[nCscRechits] = cscRechit.hitWire();
      cscRechitsWGroupsBX[nCscRechits] = cscRechit.wgroupsBX();
      cscRechitsNWireGroups[nCscRechits] = cscRechit.nWireGroups();
      cscRechitsNTimeBins[nCscRechits] = cscRechit.nTimeBins();
      cscRechitsBadStrip[nCscRechits] = cscRechit.badStrip();
      cscRechitsBadWireGroup[nCscRechits] = cscRechit.badWireGroup();
      cscRechitsEndcap[nCscRechits] = cscdetid.endcap();
      cscRechitsStation[nCscRechits] = cscdetid.station();
      cscRechitsChamber[nCscRechits] = cscdetid.chamber();
      cscRechitsRing[nCscRechits] = cscdetid.ring();
      cscRechitsLayer[nCscRechits] = cscdetid.layer();



      nCscRechits++;
    }
  }



  for(DTRecHit1DPair dtRechit: *dtRechits){
    LocalPoint  localPosition       = dtRechit.localPosition();
    DetId geoid = dtRechit.geographicalId();
    DTLayerId dtdetid = DTLayerId(geoid);
    const DTLayer* layer = dtG->layer(DTLayerId(geoid));

    if (layer) {

      GlobalPoint globalPosition = layer->toGlobal(localPosition);

  	  dtRechitPhi[nDtRechits] = globalPosition.phi();
  	  dtRechitEta[nDtRechits] = globalPosition.eta();
  	  dtRechitX[nDtRechits] = globalPosition.x();
  	  dtRechitY[nDtRechits] = globalPosition.y();
  	  dtRechitZ[nDtRechits] = globalPosition.z();
      dtRechitTime[nDtRechits] = dtRechit.digiTime();
      dtRechitWireId[nDtRechits] = dtRechit.wireId();
  	  dtRechitStation[nDtRechits] = dtdetid.station();
  	  dtRechitWheel[nDtRechits] = dtdetid.wheel();
      dtRechitSector[nDtRechits] = dtdetid.sector();
      dtRechitSuperLayer[nDtRechits] = dtdetid.superlayer();
      dtRechitLayer[nDtRechits] = dtdetid.layer();
      nDtRechits++;
      if (nDtRechits > OBJECTARRAYSIZE) {
        cout << "ERROR: nDtRechits exceeded maximum array size: " << OBJECTARRAYSIZE << "\n";
        assert(false);
      }
  	}
  }
  //*****************
  //** DIGIS
  //*****************

  nCscStripDigis = 0;
  CSCStripDigiCollection::DigiRangeIterator stripDetIt;
  for (stripDetIt = MuonCSCStripDigi->begin(); stripDetIt != MuonCSCStripDigi->end(); stripDetIt++){
     const CSCDetId &id = (*stripDetIt).first;
     int tempDetId = CSCDetId::rawIdMaker(CSCDetId::endcap(id), CSCDetId::station(id), CSCDetId::ring(id), CSCDetId::chamber(id), CSCDetId::layer(id));

      const CSCStripDigiCollection::Range &range = (*stripDetIt).second;
      for (CSCStripDigiCollection::const_iterator digiIt = range.first; digiIt != range.second; ++digiIt) {
        //ADC count is 8 time stamps
        std::vector<int> myADCVals = digiIt->getADCCounts();
         const CSCDetId &id = (*stripDetIt).first;
         int tempDetId = CSCDetId::rawIdMaker(CSCDetId::endcap(id), CSCDetId::station(id), CSCDetId::ring(id), CSCDetId::chamber(id), CSCDetId::layer(id));
         cscStripDigiDetId[nCscStripDigis] = tempDetId;
         cscStripDigiEndcap[nCscStripDigis] = id.endcap();
         cscStripDigiStation[nCscStripDigis] = id.station();
         cscStripDigiRing[nCscStripDigis] = id.ring();
         cscStripDigiChamber[nCscStripDigis] = id.chamber();
         cscStripDigiLayer[nCscStripDigis] = id.layer();
         cscStripDigiStripN[nCscStripDigis] = (*digiIt).getStrip(); //count from 1
         cscStripDigiADCCounts.push_back(myADCVals);
         nCscStripDigis++;
      }
  }
  //https://github.com/cms-sw/cmssw/blob/93f364196736af19841a8d581013e8bd46f3972d/DataFormats/CSCDigi/interface/CSCWireDigi.h
  nCscWireDigis = 0;
  CSCWireDigiCollection::DigiRangeIterator wireDetIt;
  for (wireDetIt = MuonCSCWireDigi->begin(); wireDetIt != MuonCSCWireDigi->end(); wireDetIt++){
      const CSCWireDigiCollection::Range &range = (*wireDetIt).second;
      for (CSCWireDigiCollection::const_iterator digiIt = range.first; digiIt != range.second; ++digiIt) {
        const CSCDetId &id = (*wireDetIt).first;
        int tempDetId = CSCDetId::rawIdMaker(CSCDetId::endcap(id), CSCDetId::station(id), CSCDetId::ring(id), CSCDetId::chamber(id), CSCDetId::layer(id));
        cscWireDigiDetId[nCscWireDigis] = tempDetId;
        cscWireDigiEndcap[nCscWireDigis] = id.endcap();
        cscWireDigiStation[nCscWireDigis] = id.station();
        cscWireDigiRing[nCscWireDigis] = id.ring();
        cscWireDigiChamber[nCscWireDigis] = id.chamber();
        cscWireDigiLayer[nCscWireDigis] = id.layer();
        cscWireDigiWireGroupN[nCscWireDigis] = (*digiIt).getWireGroup(); //count from 1
        cscWireDigiWireGroupBx[nCscWireDigis] = (*digiIt).getWireGroupBX(); /// return BX assigned for the wire group (16 upper bits from the wire group number)
        cscWireDigiWireTimeBin[nCscWireDigis] = (*digiIt).getTimeBin();   /// return first tbin ON number
        std::vector<int> tbins = (*digiIt).getTimeBinsOn();
        cscWireDigiWireTimeBinsOn.push_back(tbins);
        nCscWireDigis++;
      }
  }
  //https://github.com/cms-sw/cmssw/blob/6d2f66057131baacc2fcbdd203588c41c885b42c/DataFormats/DTDigi/interface/DTDigi.h
  cout<<"here"<<endl;
  DTDigiCollection::DigiRangeIterator detUnitIt;
  for (detUnitIt = MuonDTDigi->begin(); detUnitIt != MuonDTDigi->end(); ++detUnitIt) {
    const DTLayerId &id = (*detUnitIt).first;
    const DTDigiCollection::Range &range = (*detUnitIt).second;
    cout<<"inside"<<endl;
    // Loop over the digis of this DetUnit
    for (DTDigiCollection::const_iterator digiIt = range.first; digiIt != range.second; ++digiIt) {
      cout<<"inside2"<<endl;

      dtDigiTime[nDtDigis]  = (*digiIt).time(); /// Get time in ns
      dtDigiWireId[nDtDigis]  = (*digiIt).wire();
      dtDigiNumber[nDtDigis]  = (*digiIt).number(); /// Identifies different digis within the same cell
      dtDigiCountsTDC[nDtDigis]  = (*digiIt).countsTDC(); /// Get raw TDC count
      dtDigiLayer[nDtDigis]  = id.layer(); /// Get the TDC base (counts per BX)
      dtDigiSuperLayer[nDtDigis]  =id.superlayer(); /// Get the TDC base (counts per BX)
      dtDigiStation[nDtDigis] = id.station();
      dtDigiWheel[nDtDigis] = id.wheel();
      nDtDigis++;
    }
  }
  //*****************
  //** DIGI SIM LINKS
  //*****************
  /*
  edm::DetSetVector<StripDigiSimLink>::const_iterator wireLinkIt;
  //loop through detid
  nCscWireDigiSimLink=0;
  for (wireLinkIt = MuonCSCWireDigiSimLinks->begin(); wireLinkIt != MuonCSCWireDigiSimLinks->end(); wireLinkIt++){
    edm::DetSetVector<StripDigiSimLink>::detset::const_iterator i;
    for (  i = (*wireLinkIt).data.begin(); i != (*wireLinkIt).data.end(); i++){
      const CSCDetId &id = (*wireLinkIt).detId();
      // cout<<"detid, trackid, channel, CFposition, fraction: "<<id<<","<<(*wireLinkIt).detId()<<", "<<(*i).SimTrackId()<<", "<<(*i).channel()<<", "<<(*i).CFposition()<<", "<<(*i).fraction()<<endl;
      cscWireDigiSimLinkDetId[nCscWireDigiSimLink] = (*wireLinkIt).detId();
      cscWireDigiSimLinkSimTrackId[nCscWireDigiSimLink] = (*i).SimTrackId();
      cscWireDigiSimLinkChannel[nCscWireDigiSimLink] = (*i).channel();
      cscWireDigiSimLinkCFposition[nCscWireDigiSimLink] = (*i).CFposition();
      cscWireDigiSimLinkFraction[nCscWireDigiSimLink] = (*i).fraction();
      nCscWireDigiSimLink++;
    }
  }
  //count occupancy
  int index = 0;
  int occupancy = 1;
  int previous_detid = 0;
  int previous_channel = 0;
  cscWireDigiSimLinkNWires = 0;
  for(int i = 0; i < nCscWireDigiSimLink; i++)
  {
    if (previous_channel!=cscWireDigiSimLinkChannel[i] || previous_detid!=cscWireDigiSimLinkDetId[i]){
      for(int j = i;j>=index;j--)
      {
        cscWireDigiSimLinkOccupancy[j] = occupancy;
      }
      occupancy = 0;
      cscWireDigiSimLinkNWires++;
      previous_channel = cscWireDigiSimLinkChannel[i];
      previous_detid = cscWireDigiSimLinkDetId[i];
      index = i;
    }
    occupancy++;
    if (i == nCscWireDigiSimLink-1)
    {
      for(int j = i;j>index;j--)
      {
        cscWireDigiSimLinkOccupancy[j] = occupancy;
      }
      // cscWireDigiSimLinkNWires++;
    }
  }


  // strip digi simlink
  edm::DetSetVector<StripDigiSimLink>::const_iterator stripLinkIt;
  //loop through detid
  nCscStripDigiSimLink=0;
  for (stripLinkIt = MuonCSCStripDigiSimLinks->begin(); stripLinkIt != MuonCSCStripDigiSimLinks->end(); stripLinkIt++){
    edm::DetSetVector<StripDigiSimLink>::detset::const_iterator i;
    for (  i = (*stripLinkIt).data.begin(); i != (*stripLinkIt).data.end(); i++){
      const CSCDetId &id = (*stripLinkIt).detId();
      int tempDetId = CSCDetId::rawIdMaker(CSCDetId::endcap(id), CSCDetId::station(id), CSCDetId::ring(id), CSCDetId::chamber(id), CSCDetId::layer(id));
      // cout<<"detid, trackid, channel, CFposition, fraction: "<<tempDetId<<","<<(*wireLinkIt).detId()<<", "<<(*i).SimTrackId()<<", "<<(*i).channel()<<", "<<(*i).CFposition()<<", "<<(*i).fraction()<<endl;
      cscStripDigiSimLinkDetId[nCscStripDigiSimLink] = (*stripLinkIt).detId();
      cscStripDigiSimLinkSimTrackId[nCscStripDigiSimLink] = (*i).SimTrackId();
      cscStripDigiSimLinkChannel[nCscStripDigiSimLink] = (*i).channel();
      cscStripDigiSimLinkCFposition[nCscStripDigiSimLink] = (*i).CFposition();
      cscStripDigiSimLinkFraction[nCscStripDigiSimLink] = (*i).fraction();
      nCscStripDigiSimLink++;
    }

  }
  //count occupancy
  index = 0;
  occupancy = 1;
  previous_detid = 0;
  previous_channel = 0;
  cscStripDigiSimLinkNStrips = 0;
  for(int i = 0; i < nCscStripDigiSimLink; i++)
  {
    if (previous_channel!=cscStripDigiSimLinkChannel[i] || previous_detid!=cscStripDigiSimLinkDetId[i]){
      for(int j = i;j>=index;j--)
      {
        cscStripDigiSimLinkOccupancy[j] = occupancy;
      }
      occupancy = 0;
      cscStripDigiSimLinkNStrips++;
      previous_channel = cscStripDigiSimLinkChannel[i];
      previous_detid = cscStripDigiSimLinkDetId[i];
      index = i;
    }
    occupancy++;
    if (i == nCscStripDigiSimLink-1)
    {
      for(int j = i;j>index;j--)
      {
        cscStripDigiSimLinkOccupancy[j] = occupancy;
      }
      // cscStripDigiSimLinkNStrips++;
    }
  }
  */


  //*****************
  //** SIM HITS
  //*****************
  // https://github.com/cms-sw/cmssw/blob/master/SimDataFormats/TrackingHit/interface/PSimHit.h
  // vector<int> detLayers;
  // vector<bool> layer_in_ME1112;
  if (!isData)
  {
    for(size_t i=0; i<MuonDTSimHits->size();i++)
    {
      dtSimHitsTOF[nDtSimHits] = (*MuonDTSimHits)[i].timeOfFlight();
      dtSimHitsEnergyLoss[nDtSimHits] = (*MuonDTSimHits)[i].energyLoss();
      dtSimHitsParticleType[nDtSimHits] = (*MuonDTSimHits)[i].particleType();
      dtSimHitsPabs[nDtSimHits] = (*MuonDTSimHits)[i].pabs();
      dtSimHitsProcessType[nDtSimHits] = (*MuonDTSimHits)[i].processType();

      Local3DPoint dtSimHitLocalPosition = (*MuonDTSimHits)[i].localPosition();
      // DTDetId dtdetid = (DTDetId)(*MuonDTSimHits)[i].detUnitId();
      dtSimHitsDetId[nDtSimHits] = (*MuonDTSimHits)[i].detUnitId();
      const DTLayer* dtlayer = dtG->layer(DTLayerId(dtSimHitsDetId[nDtSimHits]));
      GlobalPoint globalPosition = dtlayer->toGlobal(dtSimHitLocalPosition);



      dtSimHitsX[nDtSimHits] = globalPosition.x();
      dtSimHitsY[nDtSimHits] = globalPosition.y();
      dtSimHitsZ[nDtSimHits] = globalPosition.z();
      dtSimHitsPhi[nDtSimHits] = globalPosition.phi();
      dtSimHitsEta[nDtSimHits] = globalPosition.eta();

      nDtSimHits++;
    }

    for(size_t i=0; i<MuonCSCSimHits->size();i++)
    {
      cscSimHitsTOF[nCscSimHits] = (*MuonCSCSimHits)[i].timeOfFlight();
      cscSimHitsEnergyLoss[nCscSimHits] = (*MuonCSCSimHits)[i].energyLoss();
      cscSimHitsParticleType[nCscSimHits] = (*MuonCSCSimHits)[i].particleType();
      cscSimHitsPabs[nCscSimHits] = (*MuonCSCSimHits)[i].pabs();
      cscSimHitsProcessType[nCscSimHits] = (*MuonCSCSimHits)[i].processType();

      Local3DPoint cscSimHitLocalPosition = (*MuonCSCSimHits)[i].localPosition();
      CSCDetId cscdetid = (CSCDetId)(*MuonCSCSimHits)[i].detUnitId();
      cscSimHitsDetId[nCscSimHits] = CSCDetId::rawIdMaker(CSCDetId::endcap(cscdetid), CSCDetId::station(cscdetid), CSCDetId::ring(cscdetid), CSCDetId::chamber(cscdetid), CSCDetId::layer(cscdetid));

      // const CSCChamber* cscchamber = cscG->chamber(cscdetid);
      // GlobalPoint globalPosition = cscchamber->toGlobal(cscSimHitLocalPosition);

      const CSCLayer* csclayer = cscG->layer(cscdetid);
      GlobalPoint globalPosition = csclayer->toGlobal(cscSimHitLocalPosition);
      cscSimHitsX[nCscSimHits] = globalPosition.x();
      cscSimHitsY[nCscSimHits] = globalPosition.y();
      cscSimHitsZ[nCscSimHits] = globalPosition.z();
      cscSimHitsPhi[nCscSimHits] = globalPosition.phi();
      cscSimHitsEta[nCscSimHits] = globalPosition.eta();
      cscSimHitsEndcap[nCscSimHits] = cscdetid.endcap();
      cscSimHitsStation[nCscSimHits] = cscdetid.station();
      cscSimHitsRing[nCscSimHits] = cscdetid.ring();
      cscSimHitsChamber[nCscSimHits] =cscdetid.chamber();
      cscSimHitsLayer[nCscSimHits] = cscdetid.layer();

      nCscSimHits++;
    }
  }


  // count number of simhits per layer
  /*
  std::vector<muonCscLayer> cscLayers;
  nCscDetLayer = detLayers.size();
  for(unsigned int i = 0; i < nCscDetLayer; i++)
  {
    muonCscLayer tmpLayer;
    tmpLayer.match_gParticle_index = 999;
    tmpLayer.wire_digi_occupancy = 0;
    tmpLayer.strip_digi_occupancy = 0;
    tmpLayer.nStripDigis = 0;
    tmpLayer.nWireDigis = 0;
    tmpLayer.nSimHits = 0;
    tmpLayer.nRecHits = 0;
    for (int j = 0; j < nCscRecHits; j++)
    {
      if (detLayers[i] == cscRecHitsDetId[j])
      tmpLayer.nRecHits++;
    }
    for (int j = 0; j < nCscSimHits; j++)
    {
      if (detLayers[i] == cscSimHitsDetId[j])
      {
        tmpLayer.match_gParticle_index = cscSimHits_match_gParticle_index[j];
        tmpLayer.nSimHits++;
      }

    }
    //count number of strips digis in layer
    for (int j = 0; j < nCscStripDigis; j++){
      if (detLayers[i] == cscStripDigiDetId[j])tmpLayer.nStripDigis++;

    }
    //count number of wire digis in layer
    for (int j = 0; j < nCscWireDigis; j++){
      if (detLayers[i] == cscWireDigiDetId[j])tmpLayer.nWireDigis++;
    }
    for (int j = 0; j < nCscStripDigiSimLink; j++)
    {
      if (detLayers[i] == cscStripDigiSimLinkDetId[j])tmpLayer.strip_digi_occupancy = cscStripDigiSimLinkOccupancy[j];
    }
    for (int j = 0; j < nCscWireDigiSimLink; j++)
    {
      if (detLayers[i] == cscWireDigiSimLinkDetId[j])tmpLayer.wire_digi_occupancy = cscWireDigiSimLinkOccupancy[j];
    }
    tmpLayer.Me1112 = layer_in_ME1112[i];
    tmpLayer.id = detLayers[i];
    cscLayers.push_back(tmpLayer);
  }

  sort(cscLayers.begin(), cscLayers.end(), my_largest_muonCscLayer_);

  for(unsigned int i = 0; i < cscLayers.size();i++)
  {

    cscDetLayer[i] = cscLayers[i].id;
    cscDetLayer_nCscRecHits[i]  = cscLayers[i].nRecHits;
    cscDetLayer_nCscSimHits[i]  = cscLayers[i].nSimHits;
    cscDetLayer_nWireDigis[i]  = cscLayers[i].nWireDigis;
    cscDetLayer_nStripDigis[i]  = cscLayers[i].nStripDigis;
    cscDetLayer_inME1112[i]  = cscLayers[i].Me1112;
    cscDetLayer_stripDigiOccupancy[i] = cscLayers[i].strip_digi_occupancy;
    cscDetLayer_wireDigiOccupancy[i] = cscLayers[i].wire_digi_occupancy;
    cscDetLayer_match_gParticle_index[i] = cscLayers[i].match_gParticle_index;

  }
  */



    return true;
}

bool displacedJetMuon_rechit_studies::fillGenParticles(){
  std::vector<const reco::Candidate*> prunedV;//Allows easier comparison for mother finding
  //Fills selected gen particles
  //double pt_cut = isFourJet ? 20.:20.;//this needs to be done downstream
  const double pt_cut = 0.0;
  int llp_id = 9000006;

  for(size_t i=0; i<genParticles->size();i++)
  {
    if(
       (abs((*genParticles)[i].pdgId()) >= 1 && abs((*genParticles)[i].pdgId()) <= 6 && ( (*genParticles)[i].status() < 30 ))
       || (abs((*genParticles)[i].pdgId()) >= 11 && abs((*genParticles)[i].pdgId()) <= 16)
       || (abs((*genParticles)[i].pdgId()) == 21 && (*genParticles)[i].status() < 30)
       || (abs((*genParticles)[i].pdgId()) >= 23 && abs((*genParticles)[i].pdgId()) <= 25 && ( (*genParticles)[i].status() < 30))
       || (abs((*genParticles)[i].pdgId()) >= 32 && abs((*genParticles)[i].pdgId()) <= 42)
       || (abs((*genParticles)[i].pdgId()) >= 32 && abs((*genParticles)[i].pdgId()) <= 42)

       // || (abs((*genParticles)[i].pdgId()) >= 1000001 && abs((*genParticles)[i].pdgId()) <= 1000039)
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

  //Look for mother particle and Fill gen variables
  for(unsigned int i = 0; i < prunedV.size(); i++)
  {
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



    if(prunedV[i]->numberOfMothers() > 0)
    {
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
      for(unsigned int j = 0; j < prunedV.size(); j++)
      {
        if(prunedV[j] == originalMotherWithSameID)
        {
          gParticleMotherIndex[i] = j;
          break;
        }
      }
    }
    else
    {
      gParticleMotherIndex[i] = -1;
    }

    //---------------------------------------
    //Find LLPs production and decay vertices
    //---------------------------------------
    if (enableGenLLPInfo_)
    {
      if ( abs(gParticleId[i]) == llp_id  && gParticleStatus[i] == 22 )
      {
        if (gParticleId[i] == llp_id)
        {
          gLLP_prod_vertex_x[0] = prunedV[i]->vx();
          gLLP_prod_vertex_y[0] = prunedV[i]->vy();
          gLLP_prod_vertex_z[0] = prunedV[i]->vz();
        }
        else if (gParticleId[i] == -1*llp_id)
        {
          gLLP_prod_vertex_x[1] = prunedV[i]->vx();
          gLLP_prod_vertex_y[1] = prunedV[i]->vy();
          gLLP_prod_vertex_z[1] = prunedV[i]->vz();
        }

        const reco::Candidate *dau = 0;
        bool foundDaughter = false;
        bool noDaughter = false;
        const reco::Candidate *tmpParticle = prunedV[i];

        while (!foundDaughter && !noDaughter)
        {
          if (tmpParticle->numberOfDaughters() > 0)
          {
            dau = tmpParticle->daughter(0);
            if (dau && (dau->pdgId() != llp_id && dau->pdgId() != llp_id+1))
            {
              foundDaughter = true;
            } else
            {
              tmpParticle = dau;
            }
          }
          else
          {
            noDaughter = true;
          }
        }

        if (foundDaughter)
        {

          if (gParticleId[i] == llp_id)
          {
            gLLP_decay_vertex_x[0] = dau->vx();
            gLLP_decay_vertex_y[0] = dau->vy();
            gLLP_decay_vertex_z[0] = dau->vz();
            gLLP_pt[0] = sqrt(gParticlePx[i]*gParticlePx[i]+gParticlePy[i]*gParticlePy[i]);
            gLLP_e[0] = gParticleE[i];
            gLLP_eta[0] = gParticleEta[i];
            gLLP_phi[0] = gParticlePhi[i];
            gLLP_beta[0] = sqrt(gParticlePx[i]*gParticlePx[i]+gParticlePy[i]*gParticlePy[i]+gParticlePz[i]*gParticlePz[i])/gParticleE[i];

            gLLP_gamma[0] = 1./sqrt(1-gLLP_beta[0]*gLLP_beta[0]);
            float gLLP_decay_vertex = sqrt(pow(gLLP_decay_vertex_x[0], 2) + pow(gLLP_decay_vertex_y[0], 2) + pow(gLLP_decay_vertex_z[0],2));
            gLLP_ctau[0] = gLLP_decay_vertex/(gLLP_beta[0] * gLLP_gamma[0]);

            for(size_t j=0; j<genParticles->size();j++)
            {
               int daughterPdg = (*genParticles)[j].pdgId();
               if ((*genParticles)[j].status() != 1)continue;

               if (abs(daughterPdg)==12 || abs(daughterPdg)==14 || abs(daughterPdg)==16 || abs(daughterPdg)==13)continue; // ignore neutrinos and muons
               if (abs(daughterPdg) > 1000000) continue;//ignore BSM particles
               const reco::Candidate *tempCandidate = &(*genParticles)[j];
               while(tempCandidate->numberOfMothers() > 0 && tempCandidate->pdgId() != gParticleId[i])
               {
                 tempCandidate  = tempCandidate->mother(0);
                 // if (tempCandidate->PID == pdgCode)break;
               }
               if (tempCandidate->numberOfMothers() == 0) continue;
               cout<<daughterPdg<<","<<(*genParticles)[i].energy()<<","<<tempCandidate->mother(0)->pdgId()<<endl;
               if (abs(daughterPdg)==11 || abs(daughterPdg)==22 || abs(daughterPdg)==111)gLLP_Eem[0] += (*genParticles)[j].energy();
               else gLLP_Ehad[0] += (*genParticles)[j].energy();
            }

            gLLP_travel_time[0] = sqrt(pow(gLLP_decay_vertex_x[0]-gLLP_prod_vertex_x[0],2)
                                    +pow(gLLP_decay_vertex_y[0]-gLLP_prod_vertex_y[0],2)
                                    +pow(gLLP_decay_vertex_z[0]-gLLP_prod_vertex_z[0],2))/(30. * gLLP_beta[0]);//1/30 is to convert cm to ns

            double radius = sqrt( pow(gLLP_decay_vertex_x[0],2) + pow(gLLP_decay_vertex_y[0],2) );
            if (abs(gLLP_eta[0]) < 2.4 && abs(gLLP_eta[0]) > 0.9
               && abs(gLLP_decay_vertex_z[0])<1100 && abs(gLLP_decay_vertex_z[0])>568
               && radius < 695.5) gLLP_csc[0] = true;
            double ecal_radius = 129.0;
            // double hcal_radius = 179.0;
            double EB_z = 268.36447217; // 129*sinh(1.479)
            double EE_z = 298.5; //where Ecal Endcap starts in z direction

            for (unsigned int id = 0; id < tmpParticle->numberOfDaughters(); id++ )
            {
            //std::cout << "====================" << std::endl;
            //std::cout << " -> "<< tmpParticle->daughter(id)->pdgId() << std::endl;
              if( id > 1 ) break;
              TLorentzVector tmp;
              tmp.SetPxPyPzE(tmpParticle->daughter(id)->px(), tmpParticle->daughter(id)->py(), tmpParticle->daughter(id)->pz(), tmpParticle->daughter(id)->energy());
              if(tmp.Pt()<pt_cut) continue;
              gLLP_daughter_pt[id] = tmp.Pt();
              gLLP_daughter_eta[id] = tmp.Eta();
              gLLP_daughter_phi[id] = tmp.Phi();
              gLLP_daughter_e[id]  = tmp.E();

              // double gLLP_daughter_travel_time_hcal= (1./30.)*(hcal_radius-radius)/(tmp.Pt()/tmp.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns

              gLLP_daughter_travel_time[id] = (1./30.)*(ecal_radius-radius)/(tmp.Pt()/tmp.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns
            //Calculate dt from generation point to ECAL face
              double x_ecal = gLLP_decay_vertex_x[0] + 30. * (tmp.Px()/tmp.E())*gLLP_daughter_travel_time[id];
              double y_ecal = gLLP_decay_vertex_y[0] + 30. * (tmp.Py()/tmp.E())*gLLP_daughter_travel_time[id];
              double z_ecal = gLLP_decay_vertex_z[0] + 30. * (tmp.Pz()/tmp.E())*gLLP_daughter_travel_time[id];
              // double x_hcal = gLLP_decay_vertex_x[0] + 30. * (tmp.Px()/tmp.E())*gLLP_daughter_travel_time_hcal;
              // double y_hcal = gLLP_decay_vertex_y[0] + 30. * (tmp.Py()/tmp.E())*gLLP_daughter_travel_time_hcal;
              // double z_hcal = gLLP_decay_vertex_z[0] + 30. * (tmp.Pz()/tmp.E())*gLLP_daughter_travel_time_hcal;

              if( fabs(z_ecal) < EB_z && radius <= ecal_radius &&  fabs(gLLP_decay_vertex_z[0]) < EE_z)
              {
                photon_travel_time[id] = (1./30) * sqrt(pow(ecal_radius,2)+pow(z_ecal,2));
                photon_travel_time_pv[id] = (1./30) * sqrt(pow(x_ecal-genVertexX,2) + pow(y_ecal-genVertexY,2) + pow(z_ecal-genVertexZ,2));
                gen_time_pv[id] =  gLLP_travel_time[0] + gLLP_daughter_travel_time[id] - photon_travel_time_pv[id] + genVertexT;
                gen_time[id] = gLLP_travel_time[0] + gLLP_daughter_travel_time[id] - photon_travel_time[id] + genVertexT;

              }
              else
              {
                gLLP_daughter_travel_time[id] = -666;
                gen_time_pv[id] = -666.;
                gen_time[id] = -666.;
                photon_travel_time[id] = -666.;
                photon_travel_time_pv[id] = -666.;
              }
              double min_delta_r = 666.;
              double min_delta_r_calo = 666.;
              unsigned int match_jet_index = 666;
              unsigned int match_calojet_index = 666;


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
              for ( int i_jet = 0; i_jet < nCaloJets; i_jet++ )
              {
                double current_delta_r = deltaR(gLLP_daughter_eta_ecalcorr[id], gLLP_daughter_phi_ecalcorr[id], calojetEta[i_jet], calojetPhi[i_jet]);
                if ( current_delta_r < min_delta_r_calo )
                {
                  // min_delta_r_nocorr = deltaR(gLLP_daughter_eta[id], gLLP_daughter_phi[id], jetEta[i_jet], jetPhi[i_jet]);
                  min_delta_r_calo = current_delta_r;
                  match_calojet_index = i_jet;
                }
                  // std::cout << i_jet << " min dR = " << min_delta_r_calo << std::endl;
              }
              for ( int i_jet = 0; i_jet < nJets; i_jet++ )
              {
                double current_delta_r = deltaR(gLLP_daughter_eta_ecalcorr[id], gLLP_daughter_phi_ecalcorr[id], jetEta[i_jet], jetPhi[i_jet]);
                if ( current_delta_r < min_delta_r )
                {
                  min_delta_r = current_delta_r;
                  match_jet_index = i_jet;
            //std::cout << i_jet << " min dR = " << min_delta_r << std::endl;
                }
              }//end matching to jets using ECAL radius
              // std::cout<<"min delta r "<<min_delta_r<<std::endl;

              if ( min_delta_r < 0.3 )
              {
                gLLP_daughter_match_jet_index[id] = match_jet_index;
                gLLP_min_delta_r_match_jet[id] = min_delta_r;
                jet_matched[match_jet_index] = true;
                // std::cout<<"gllp daughted match jet index "<<gLLP_daughter_match_jet_index[id]<<std::endl;

                // std::cout<<"matched "<<jet_matched[match_jet_index]<<std::endl;
              }

              if (min_delta_r_calo < 0.3)
              {
                gLLP_daughter_match_calojet_index[id] = match_calojet_index;
                gLLP_min_delta_r_match_calojet[id] = min_delta_r_calo;

              }


            }
          }
          else if (gParticleId[i] == -1*llp_id)
          {
            gLLP_decay_vertex_x[1] = dau->vx();
            gLLP_decay_vertex_y[1] = dau->vy();
            gLLP_decay_vertex_z[1] = dau->vz();
            gLLP_pt[1] = sqrt(gParticlePx[i]*gParticlePx[i]+gParticlePy[i]*gParticlePy[i]);
            gLLP_e[1] = gParticleE[i];
            gLLP_eta[1] = gParticleEta[i];
            gLLP_phi[1] = gParticlePhi[i];
            gLLP_beta[1] = sqrt(gParticlePx[i]*gParticlePx[i]+gParticlePy[i]*gParticlePy[i]+gParticlePz[i]*gParticlePz[i])/gParticleE[i];
            gLLP_gamma[1] = 1./sqrt(1-gLLP_beta[1]*gLLP_beta[1]);
            float gLLP_decay_vertex = sqrt(pow(gLLP_decay_vertex_x[1], 2) + pow(gLLP_decay_vertex_y[1], 2) + pow(gLLP_decay_vertex_z[1],2));
            gLLP_ctau[1] = gLLP_decay_vertex/(gLLP_beta[1] * gLLP_gamma[1]);

            for(size_t j=0; j<genParticles->size();j++)
            {
               int daughterPdg = (*genParticles)[j].pdgId();
               if ((*genParticles)[j].status() != 1)continue;

               if (abs(daughterPdg)==12 || abs(daughterPdg)==14 || abs(daughterPdg)==16 || abs(daughterPdg)==13)continue; // ignore neutrinos and muons
               if (abs(daughterPdg) > 1000000) continue;//ignore BSM particles
               const reco::Candidate *tempCandidate = &(*genParticles)[j];
               while(tempCandidate->numberOfMothers() > 0 && tempCandidate->pdgId() != gParticleId[i])
               {
                 tempCandidate  = tempCandidate->mother(0);
                 // if (tempCandidate->PID == pdgCode)break;
               }
               if (tempCandidate->numberOfMothers() == 0) continue;
               if (abs(daughterPdg)==11 || abs(daughterPdg)==22 || abs(daughterPdg)==111)gLLP_Eem[1] += (*genParticles)[j].energy();
               else gLLP_Ehad[1] += (*genParticles)[j].energy();
            }

            gLLP_travel_time[1] = sqrt(pow(gLLP_decay_vertex_x[1]-gLLP_prod_vertex_x[1],2)
                    +pow(gLLP_decay_vertex_y[1]-gLLP_prod_vertex_y[1],2)
                    +pow(gLLP_decay_vertex_z[1]-gLLP_prod_vertex_z[1],2))/(30. * gLLP_beta[1]);//1/30 is to convert cm to ns
            double radius = sqrt( pow(gLLP_decay_vertex_x[1],2) + pow(gLLP_decay_vertex_y[1],2) );
            if (abs(gLLP_eta[1]) < 2.4 && abs(gLLP_eta[1]) > 0.9
               && abs(gLLP_decay_vertex_z[1])<1100 && abs(gLLP_decay_vertex_z[1])>568
               && radius < 695.5) gLLP_csc[1] = true;
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
              gLLP_daughter_pt[id+2] = tmp.Pt();
              gLLP_daughter_eta[id+2] = tmp.Eta();
              gLLP_daughter_phi[id+2] = tmp.Phi();
              gLLP_daughter_e[id+2]  = tmp.E();
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

              }
              else
              {
                gLLP_daughter_travel_time[id+2] = -666;
                gen_time_pv[id+2] = -666.;
                gen_time[id+2] = -666.;
                photon_travel_time[id+2] = -666.;
                photon_travel_time_pv[id+2] = -666.;
              }

              double min_delta_r = 666.;
              double min_delta_r_calo = 666.;
              unsigned int match_jet_index = 666;
              unsigned int match_calojet_index = 666;

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
              for ( int i_jet = 0; i_jet < nCaloJets; i_jet++ )
              {
                double current_delta_r = deltaR(gLLP_daughter_eta_ecalcorr[id+2], gLLP_daughter_phi_ecalcorr[id+2], calojetEta[i_jet], calojetPhi[i_jet]);
                if ( current_delta_r < min_delta_r_calo )
                {
                  // min_delta_r_nocorr = deltaR(gLLP_daughter_eta[id], gLLP_daughter_phi[id], jetEta[i_jet], jetPhi[i_jet]);
                  min_delta_r_calo = current_delta_r;
                  match_calojet_index = i_jet;
                 }
                      //std::cout << i_jet << " min dR = " << min_delta_r << std::endl;
             }// end matching to calojets
              for ( int i_jet = 0; i_jet < nJets; i_jet++ )
              {
                double current_delta_r = deltaR(gLLP_daughter_eta_ecalcorr[id+2], gLLP_daughter_phi_ecalcorr[id+2] , jetEta[i_jet], jetPhi[i_jet]);
                if ( current_delta_r < min_delta_r )
                {
                  min_delta_r = current_delta_r;
                  match_jet_index = i_jet;
                }
              }//end matching to jets ecal

              if ( min_delta_r < 0.3 )
              {
                gLLP_daughter_match_jet_index[id+2] = match_jet_index;
                gLLP_min_delta_r_match_jet[id+2] = min_delta_r;
                jet_matched[match_jet_index] = true;
                // std::cout<<"gllp daughted match jet index "<<gLLP_daughter_match_jet_index[id+2]<<std::endl;

                // std::cout<<"matched "<<jet_matched[match_jet_index]<<std::endl;
              }
              if ( min_delta_r_calo < 0.3 )
              {
                gLLP_daughter_match_calojet_index[id+2] = match_calojet_index;
                gLLP_min_delta_r_match_calojet[id+2] = min_delta_r_calo;
              }

            }//for daughters loop
          }//if particle ID = 36
        }//if found daughters
      }



    }




  }// for loop of genParticles
  return true;
};


bool displacedJetMuon_rechit_studies::fillJets(const edm::EventSetup& iSetup)
{

  for (const reco::PFJet &j : *jets)
  {
    //if (j.pt() < 15) continue;
    //if (fabs(j.eta()) > 2.4) continue;
    //-------------------
    //Fill Jet-Level Info
    //-------------------
    jetE[nJets] = j.energy();
    jetPt[nJets] = j.pt();
    jetEta[nJets] = j.eta();
    jetPhi[nJets] = j.phi();
    jetMass[nJets] = j.mass();
    jetEt[nJets] = j.et();

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

    edm::ESHandle<CaloGeometry> geoHandle;
    iSetup.get<CaloGeometryRecord>().get(geoHandle);
    const CaloSubdetectorGeometry *barrelGeometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
    //const CaloSubdetectorGeometry *endcapGeometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);
    //double ecal_radius = 129.0;
    int n_matched_rechits = 0;
    std::vector<double> rechitphi;
    std::vector<double> rechiteta;
    std::vector<double> rechitet;
    std::vector<double> rechitt;
    for (EcalRecHitCollection::const_iterator recHit = ebRecHits->begin(); recHit != ebRecHits->end(); ++recHit)
    {
      const DetId recHitId = recHit->detid();
      const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();

      if (enableEcalRechits_ && recHit->energy() > 0.5 && recHit->timeError() > 0.0 && recHit->timeError() < 100.0)
      {
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
        nRechits ++;

      }

      if (recHit->checkFlag(EcalRecHit::kSaturated) || recHit->checkFlag(EcalRecHit::kLeadingEdgeRecovered) || recHit->checkFlag(EcalRecHit::kPoorReco) || recHit->checkFlag(EcalRecHit::kWeird) || recHit->checkFlag(EcalRecHit::kDiWeird)) continue;
      if (recHit->timeError() < 0 || recHit->timeError() > 100) continue;
      if (abs(recHit->time()) > 12.5) continue;

      if ( deltaR(jetEta[nJets], jetPhi[nJets], recHitPos.eta(), recHitPos.phi())  < 0.4)
      {
        //double rechit_x = ecal_radius * cos(recHitPos.phi());
        //double rechit_y = ecal_radius * sin(recHitPos.phi());
        //double rechit_z = ecal_radius * sinh(recHitPos.eta());
        //double photon_pv_travel_time = (1./30) * sqrt(pow(pvX-rechit_x,2)+pow(pvY-rechit_y,2)+pow(pvZ-rechit_z,2));
        if ( deltaR(jetEta[nJets], jetPhi[nJets], recHitPos.eta(), recHitPos.phi())  < 0.15)
        {
          if (recHit->energy() > Rechit_cut) jet_energy_frac[nJets] += recHit->energy();
        }

        if (recHit->energy() > Rechit_cut)
        {
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
    }

    jetRechitT[nJets] = jetRechitT[nJets]/jetRechitE[nJets];
    jetNRechits[nJets] = n_matched_rechits;
    jetRechitE_Error[nJets] = sqrt(jetRechitE_Error[nJets]);
    jetRechitT_rms[nJets] = sqrt(jetRechitT_rms[nJets]);
    double sig1(0.0),sig2(0.0);
    jet_second_moments(rechitet,rechiteta,rechitphi,sig1,sig2);
    jet_sig_et1[nJets] = sig1;
    jet_sig_et2[nJets] = sig2;
    jet_energy_frac[nJets] = jet_energy_frac[nJets]/jetRechitE[nJets];

    nJets++;
  } //loop over jets

  return true;
};

bool displacedJetMuon_rechit_studies::fillMuons(const edm::Event& iEvent)
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

bool displacedJetMuon_rechit_studies::passCaloJetID( const reco::CaloJet *jetCalo, int cutLevel) {
  bool result = false;

  return result;
}//passJetID CaloJet

bool displacedJetMuon_rechit_studies::passJetID( const reco::PFJet *jet, int cutLevel) {
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

void displacedJetMuon_rechit_studies::findTrackingVariables(const TLorentzVector &jetVec,const edm::EventSetup& iSetup,float &alphaMax,float &medianTheta2D,float &medianIP, int &nTracksPV,float &ptAllPVTracks,float &ptAllTracks,float &minDeltaRAllTracks, float &minDeltaRPVTracks)
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
	    if (minDeltaRAllTracks > generalTrackVecTemp.DeltaR(jetVec))
	    {
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

void displacedJetMuon_rechit_studies::jet_second_moments(std::vector<double> &et,std::vector<double> &eta,std::vector<double> &phi,double &sig1,double &sig2)
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

const reco::Candidate* displacedJetMuon_rechit_studies::findFirstMotherWithDifferentID(const reco::Candidate *particle){

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

const reco::Candidate* displacedJetMuon_rechit_studies::findOriginalMotherWithSameID(const reco::Candidate *particle){

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

bool displacedJetMuon_rechit_studies::fillMet(const edm::Event& iEvent)
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

  return true;
};

bool displacedJetMuon_rechit_studies::fillTrigger(const edm::Event& iEvent)
{

  //fill trigger information
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  // std::cout << "\n === TRIGGER PATHS === " << std::endl;
  //------------------------------------------------------------------
  //Option to save all HLT path names in the ntuple per event
  //Expensive option in terms of ntuple size
  //------------------------------------------------------------------
  nameHLT->clear();
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

bool displacedJetMuon_rechit_studies::fillMC()
{
  for(const reco::GenJet &j : *genJets)
  {
    //std::cout << nGenJets << std::endl;
    genJetE[nGenJets] = j.energy();
    genJetPt[nGenJets] = j.pt();
    genJetEta[nGenJets] = j.eta();
    genJetPhi[nGenJets] = j.phi();
    genJetMET[nGenJets] = j.invisibleEnergy();
    nGenJets++;
  }

  const reco::GenMET &GenMetCalo = genMetsCalo->front();
  genMetPtCalo  = GenMetCalo.pt();
  genMetPhiCalo = GenMetCalo.phi();

  const reco::GenMET &GenMetTrue = genMetsTrue->front();
  genMetPtTrue  = GenMetTrue.pt();
  genMetPhiTrue = GenMetTrue.phi();

  bool foundGenVertex = false;
  for(size_t i=0; i<genParticles->size();i++)
  {
    if (!foundGenVertex)
    {
      for (unsigned int j=0; j<(*genParticles)[i].numberOfDaughters(); ++j)
      {
        const reco::Candidate *dau = (*genParticles)[i].daughter(j);
        if (dau)
        {
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

  return true;
};

bool displacedJetMuon_rechit_studies::fillCaloJets(const edm::EventSetup& iSetup)
{
  for (const reco::CaloJet &j : *jetsCalo) {

    //if (j.pt() < 20) continue;
    //if (fabs(j.eta()) > 2.4) continue;

    //-------------------
    //Fill Jet-Level Info
    //-------------------
    calojetE[nCaloJets] = j.energy();
    calojetEt[nCaloJets] = j.et();
    calojetPt[nCaloJets] = j.pt();
    calojetEta[nCaloJets] = j.eta();
    calojetPhi[nCaloJets] = j.phi();
    calojetMass[nCaloJets] = j.mass();
    calojet_HadronicEnergyFraction[nCaloJets] = j.energyFractionHadronic();
    calojet_EMEnergyFraction[nCaloJets] = j.emEnergyFraction();

    TLorentzVector thisJet;
    thisJet.SetPtEtaPhiE(calojetPt[nCaloJets], calojetEta[nCaloJets], calojetPhi[nCaloJets], calojetE[nCaloJets]);
    //calojetCISV = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    float alphaMax(0.0),medianTheta2D(0.0),medianIP(0.0),minDeltaRAllTracks(0.0),minDeltaRPVTracks(0.0),ptAllTracks(0.0), ptAllPVTracks(0.0);
    int nTracksPV(0);
    findTrackingVariables(thisJet,iSetup,alphaMax,medianTheta2D,medianIP,nTracksPV,ptAllPVTracks,ptAllTracks, minDeltaRAllTracks, minDeltaRPVTracks);
    //jetCISV = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    calojetAlphaMax[nCaloJets] = alphaMax;
    calojetBetaMax[nCaloJets] = alphaMax * ptAllTracks/j.pt();
    calojetGammaMax[nCaloJets] = alphaMax * ptAllTracks/(j.energy());
    calojetGammaMax_EM[nCaloJets] = alphaMax * ptAllTracks/(j.energy()*j.emEnergyFraction());
    calojetGammaMax_Hadronic[nCaloJets] =  alphaMax * ptAllTracks/(j.energy()*j.energyFractionHadronic());
    calojetGammaMax_ET[nCaloJets] = alphaMax * ptAllTracks/j.et();

    calojetMedianTheta2D[nCaloJets] = medianTheta2D;
    calojetMedianIP[nCaloJets] = medianIP;
    calojetPtAllPVTracks[nCaloJets] = ptAllPVTracks;
    calojetPtAllTracks[nCaloJets] = ptAllTracks;
    calojetMinDeltaRAllTracks[nCaloJets] = minDeltaRAllTracks;
    calojetMinDeltaRPVTracks[nCaloJets] = minDeltaRPVTracks;

    calojetJetArea[nCaloJets] = j.jetArea();
    calojetPileupE[nCaloJets] = j.pileup();

    calojetPileupIdFlag[nCaloJets] = 0;
    calojetPassIDLoose[nCaloJets] = passCaloJetID(&j, 0);
    calojetPassIDTight[nCaloJets] = passCaloJetID(&j, 1);


    //---------------------------
    //Find RecHits Inside the Jet
    //---------------------------
    // geometry (from ECAL ELF)

    edm::ESHandle<CaloGeometry> geoHandle;
    iSetup.get<CaloGeometryRecord>().get(geoHandle);
    const CaloSubdetectorGeometry *barrelGeometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
    //const CaloSubdetectorGeometry *endcapGeometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);
    //double ecal_radius = 129.0;
    int n_matched_rechits = 0;
    for (EcalRecHitCollection::const_iterator recHit = ebRecHits->begin(); recHit != ebRecHits->end(); ++recHit)
    {
      if (recHit->checkFlag(EcalRecHit::kSaturated) || recHit->checkFlag(EcalRecHit::kLeadingEdgeRecovered) || recHit->checkFlag(EcalRecHit::kPoorReco) || recHit->checkFlag(EcalRecHit::kWeird) || recHit->checkFlag(EcalRecHit::kDiWeird)) continue;
      if (recHit->timeError() < 0 || recHit->timeError() > 100) continue;
      if ( recHit->checkFlag(0) )
      {
        const DetId recHitId = recHit->detid();
        const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
        if ( deltaR(calojetEta[nCaloJets], calojetPhi[nCaloJets], recHitPos.eta(), recHitPos.phi())  < 0.4)
        {
          //double rechit_x = ecal_radius * cos(recHitPos.phi());
          //double rechit_y = ecal_radius * sin(recHitPos.phi());
          //double rechit_z = ecal_radius * sinh(recHitPos.eta());
          //double photon_pv_travel_time = (1./30) * sqrt(pow(pvX-rechit_x,2)+pow(pvY-rechit_y,2)+pow(pvZ-rechit_z,2));

          if (recHit->energy() > Rechit_cut)
          {
            calojetRechitE[nCaloJets] += recHit->energy();
            calojetRechitT[nCaloJets] += recHit->time()*recHit->energy();
            calojetRechitT_rms[nCaloJets] += recHit->time()*recHit->time();

            n_matched_rechits++;
          }

        }
      }
    }
    //cout << "Last Nphoton: " << fJetNPhotons << "\n";
    //std::cout << "n: " << n_matched_rechits << std::endl;
    if (calojetRechitE[nCaloJets] > 0.0 ){
      calojetNRechits[nCaloJets] = n_matched_rechits;
      calojetRechitT[nCaloJets] = calojetRechitT[nCaloJets]/calojetRechitE[nCaloJets];
      calojetRechitT_rms[nCaloJets] = sqrt(calojetRechitT_rms[nCaloJets]);
      nCaloJets++;
    }
  } //loop over calojets

  return true;
};


bool displacedJetMuon_rechit_studies::fillElectrons(const edm::Event& iEvent)
{

  // Get MVA values and categories (optional)
  // edm::Handle<edm::ValueMap<float> > mvaGeneralPurposeValues;
  // edm::Handle<edm::ValueMap<int> > mvaGeneralPurposeCategories;
  // edm::Handle<edm::ValueMap<float> > mvaHZZValues;
  // edm::Handle<edm::ValueMap<int> > mvaHZZCategories;
  // iEvent.getByToken(mvaGeneralPurposeValuesMapToken_,mvaGeneralPurposeValues);
  // iEvent.getByToken(mvaGeneralPurposeCategoriesMapToken_,mvaGeneralPurposeCategories);
  // iEvent.getByToken(mvaHZZValuesMapToken_,mvaHZZValues);
  // iEvent.getByToken(mvaHZZCategoriesMapToken_,mvaHZZCategories);

  // setupEgammaPostRecoSeq


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

bool displacedJetMuon_rechit_studies::fillPileUp()
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

//define this as a plug-in
DEFINE_FWK_MODULE(displacedJetMuon_rechit_studies);
