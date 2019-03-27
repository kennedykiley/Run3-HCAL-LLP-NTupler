// -*- C++ -*-
// Class:      llp_ntupler
/*
  Description: Base class for miniAOD analysis with CRAB
*/
//         Author:  Cristián Peña and Si Xie.
//         Created:  Thu, 26 March 2019 15:00:06 GMT

#include "llp_ntupler.h"
//------ Constructors and destructor ------//
llp_ntupler::llp_ntupler(const edm::ParameterSet& iConfig):
  //get inputs from config file
  isData_(iConfig.getParameter<bool> ("isData")),
  //isFourJet_(iConfig.getParameter<bool> ("isFourJet")),
  useGen_(iConfig.getParameter<bool> ("useGen")),
  isFastsim_(iConfig.getParameter<bool> ("isFastsim")),
  //isQCD_(iConfig.getParameter<bool> ("isQCD")),
  enableTriggerInfo_(iConfig.getParameter<bool> ("enableTriggerInfo")),
  //enableRecHitInfo_(iConfig.getParameter<bool> ("enableRecHitInfo")),
  readGenVertexTime_(iConfig.getParameter<bool> ("readGenVertexTime")),
  triggerPathNamesFile_(iConfig.getParameter<string> ("triggerPathNamesFile")),
  eleHLTFilterNamesFile_(iConfig.getParameter<string> ("eleHLTFilterNamesFile")),
  muonHLTFilterNamesFile_(iConfig.getParameter<string> ("muonHLTFilterNamesFile")),
  photonHLTFilterNamesFile_(iConfig.getParameter<string> ("photonHLTFilterNamesFile")),
  verticesToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  tracksTag_(consumes<edm::View<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracks"))),
  muonsToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronsToken_(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  tausToken_(consumes<reco::PFTauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  photonsToken_(consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  jetsToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  jetsPuppiToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jetsPuppi"))),
  jetsAK8Token_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jetsAK8"))),
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
  gedPhotonCoresToken_(consumes<vector<reco::PhotonCore> >(iConfig.getParameter<edm::InputTag>("gedPhotonCores")))
  //superClustersToken_(consumes<vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("superClusters"))),
  //  lostTracksToken_(consumes<vector<reco::PFCandidate> >(iConfig.getParameter<edm::InputTag>("lostTracks")))
{
  //declare the TFileService for output
  edm::Service<TFileService> fs;

  //set up output tree
  llpTree = fs->make<TTree>("llp", "selected AOD information for llp analyses");
  //llpTree = new TTree("Jets", "selected AOD information");
  NEvents = fs->make<TH1F>("NEvents",";;NEvents;",1,-0.5,0.5);
  //*****************************************************************************************
  //Read in HLT Trigger Path List from config file
  //*****************************************************************************************
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
  /*
  fJetPhotonRecHitEta = new std::vector<float>; fJetPhotonRecHitEta->clear();
  fJetPhotonRecHitPhi = new std::vector<float>; fJetPhotonRecHitPhi->clear();
  fJetPhotonRecHitE = new std::vector<float>; fJetPhotonRecHitE->clear();
  fJetPhotonRecHitTime = new std::vector<float>; fJetPhotonRecHitTime->clear();
*/
}

llp_ntupler::~llp_ntupler()
{
};

//***********************************************
//Enable output ntuple branches
//***********************************************

void llp_ntupler::setBranches()
{
  enableEventInfoBranches();
};
void llp_ntupler::enableEventInfoBranches(){
  llpTree->Branch("isData", &isData, "isData/O");
  llpTree->Branch("nPV", &nPV, "nPV/I");
  llpTree->Branch("runNum", &runNum, "runNum/i");
  llpTree->Branch("nSlimmedSecondV", &nSlimmedSecondV, "nSlimmedSecondV/i");
  llpTree->Branch("lumiNum", &lumiNum, "lumiNum/i");
  llpTree->Branch("eventNum", &eventNum, "eventNum/i");
  llpTree->Branch("eventTime", &eventTime, "eventTime/i");
  llpTree->Branch("pvX", &pvX, "pvX/F");
  llpTree->Branch("pvY", &pvY, "pvY/F");
  llpTree->Branch("pvZ", &pvZ, "pvZ/F");
  llpTree->Branch("fixedGridRhoAll", &fixedGridRhoAll, "fixedGridRhoAll/F");
  llpTree->Branch("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll, "fixedGridRhoFastjetAll/F");
  llpTree->Branch("fixedGridRhoFastjetAllCalo", &fixedGridRhoFastjetAllCalo, "fixedGridRhoFastjetAllCalo/F");
  llpTree->Branch("fixedGridRhoFastjetCentralCalo", &fixedGridRhoFastjetCentralCalo, "fixedGridRhoFastjetCentralCalo/F");
  llpTree->Branch("fixedGridRhoFastjetCentralChargedPileUp", &fixedGridRhoFastjetCentralChargedPileUp, "fixedGridRhoFastjetCentralChargedPileUp/F");
  llpTree->Branch("fixedGridRhoFastjetCentralNeutral", &fixedGridRhoFastjetCentralNeutral, "fixedGridRhoFastjetCentralNeutral/F");
}

void llp_ntupler::enablePVAllBranches()
{
  llpTree->Branch("nPVAll", &nPVAll,"nPVAll/I");
  llpTree->Branch("pvAllX", pvAllX,"pvAllX[nPVAll]/F");
  llpTree->Branch("pvAllY", pvAllY,"pvAllY[nPVAll]/F");
  llpTree->Branch("pvAllZ", pvAllZ,"pvAllZ[nPVAll]/F");
  llpTree->Branch("pvAllLogSumPtSq", pvAllLogSumPtSq,"pvAllLogSumPtSq[nPVAll]/F");
  llpTree->Branch("pvAllSumPx", pvAllSumPx,"pvAllSumPx[nPVAll]/F");
  llpTree->Branch("pvAllSumPy", pvAllSumPy,"pvAllSumPy[nPVAll]/F");
};

void llp_ntupler::enablePileUpBranches()
{
  llpTree->Branch("nBunchXing", &nBunchXing, "nBunchXing/I");
  llpTree->Branch("BunchXing", BunchXing, "BunchXing[nBunchXing]/I");
  llpTree->Branch("nPU", nPU, "nPU[nBunchXing]/I");
  llpTree->Branch("nPUmean", nPUmean, "nPUmean[nBunchXing]/F");
};

void llp_ntupler::enableMuonBranches()
{
  llpTree->Branch("nMuons", &nMuons,"nMuons/I");
  llpTree->Branch("muonE", muonE,"muonE[nMuons]/F");
  llpTree->Branch("muonPt", muonPt,"muonPt[nMuons]/F");
  llpTree->Branch("muonEta", muonEta,"muonEta[nMuons]/F");
  llpTree->Branch("muonPhi", muonPhi,"muonPhi[nMuons]/F");
  llpTree->Branch("muonCharge", muonCharge, "muonCharge[nMuons]/I");
  llpTree->Branch("muonIsLoose", muonIsLoose,"muonIsLoose[nMuons]/O");
  llpTree->Branch("muonIsMedium", muonIsMedium,"muonIsMedium[nMuons]/O");
  llpTree->Branch("muonIsTight", muonIsTight,"muonIsTight[nMuons]/O");
  llpTree->Branch("muon_d0", muon_d0, "muon_d0[nMuons]/F");
  llpTree->Branch("muon_dZ", muon_dZ, "muon_dZ[nMuons]/F");
  llpTree->Branch("muon_ip3d", muon_ip3d, "muon_ip3d[nMuons]/F");
  llpTree->Branch("muon_ip3dSignificance", muon_ip3dSignificance, "muon_ip3dSignificance[nMuons]/F");
  llpTree->Branch("muonType", muonType, "muonType[nMuons]/i");
  llpTree->Branch("muonQuality", muonQuality, "muonQuality[nMuons]/i");
  llpTree->Branch("muon_pileupIso", muon_pileupIso, "muon_pileupIso[nMuons]/F");
  llpTree->Branch("muon_chargedIso", muon_chargedIso, "muon_chargedIso[nMuons]/F");
  llpTree->Branch("muon_photonIso", muon_photonIso, "muon_photonIso[nMuons]/F");
  llpTree->Branch("muon_neutralHadIso", muon_neutralHadIso, "muon_neutralHadIso[nMuons]/F");
  llpTree->Branch("muon_ptrel", muon_ptrel, "muon_ptrel[nMuons]/F");
  llpTree->Branch("muon_chargedMiniIso", muon_chargedMiniIso, "muon_chargedMiniIso[nMuons]/F");
  llpTree->Branch("muon_photonAndNeutralHadronMiniIso", muon_photonAndNeutralHadronMiniIso, "muon_photonAndNeutralHadronMiniIso[nMuons]/F");
  llpTree->Branch("muon_chargedPileupMiniIso", muon_chargedPileupMiniIso, "muon_chargedPileupMiniIso[nMuons]/F");
  llpTree->Branch("muon_activityMiniIsoAnnulus", muon_activityMiniIsoAnnulus, "muon_activityMiniIsoAnnulus[nMuons]/F");
  llpTree->Branch("muon_passSingleMuTagFilter", muon_passSingleMuTagFilter, "muon_passSingleMuTagFilter[nMuons]/O");
  llpTree->Branch("muon_passHLTFilter", &muon_passHLTFilter, Form("muon_passHLTFilter[nMuons][%d]/O",MAX_MuonHLTFilters));
  llpTree->Branch("muon_validFractionTrackerHits", muon_validFractionTrackerHits, "muon_validFractionTrackerHits[nMuons]/F");
  llpTree->Branch("muon_isGlobal", muon_isGlobal,"muon_isGlobal[nMuons]/O");
  llpTree->Branch("muon_normChi2", muon_normChi2,"muon_normChi2[nMuons]/F");
  llpTree->Branch("muon_chi2LocalPosition", muon_chi2LocalPosition,"muon_chi2LocalPosition[nMuons]/F");
  llpTree->Branch("muon_kinkFinder", muon_kinkFinder,"muon_kinkFinder[nMuons]/F");
  llpTree->Branch("muon_segmentCompatability", muon_segmentCompatability,"muon_segmentCompatability[nMuons]/F");
  llpTree->Branch("muonIsICHEPMedium", muonIsICHEPMedium,"muonIsICHEPMedium[nMuons]/O");
};

void llp_ntupler::enableElectronBranches()
{
  llpTree->Branch("nElectrons", &nElectrons,"nElectrons/I");
  llpTree->Branch("eleE", eleE,"eleE[nElectrons]/F");
  llpTree->Branch("elePt", elePt,"elePt[nElectrons]/F");
  llpTree->Branch("eleEta", eleEta,"eleEta[nElectrons]/F");
  llpTree->Branch("elePhi", elePhi,"elePhi[nElectrons]/F");
  llpTree->Branch("eleCharge", eleCharge, "eleCharge[nElectrons]/F");
  //llpTree->Branch("EleE_SC", eleE_SC,"eleE_SC[nElectrons]/F");
  llpTree->Branch("eleEta_SC", eleEta_SC,"eleEta_SC[nElectrons]/F");
  //llpTree->Branch("elePhi_SC", elePhi_SC,"elePhi_SC[nElectrons]/F");
  llpTree->Branch("eleSigmaIetaIeta", eleSigmaIetaIeta, "eleSigmaIetaIeta[nElectrons]/F");
  llpTree->Branch("eleFull5x5SigmaIetaIeta", eleFull5x5SigmaIetaIeta, "eleFull5x5SigmaIetaIeta[nElectrons]/F");
  llpTree->Branch("eleR9", eleR9, "eleR9[nElectrons]/F");
  llpTree->Branch("ele_dEta", ele_dEta, "ele_dEta[nElectrons]/F");
  llpTree->Branch("ele_dPhi", ele_dPhi, "ele_dPhi[nElectrons]/F");
  llpTree->Branch("ele_HoverE", ele_HoverE, "ele_HoverE[nElectrons]/F");
  llpTree->Branch("ele_d0", ele_d0, "ele_d0[nElectrons]/F");
  llpTree->Branch("ele_dZ", ele_dZ, "ele_dZ[nElectrons]/F");
  llpTree->Branch("ele_ip3d", ele_ip3d, "ele_ip3d[nElectrons]/F");
  llpTree->Branch("ele_ip3dSignificance", ele_ip3dSignificance, "ele_ip3dSignificance[nElectrons]/F");
  llpTree->Branch("ele_pileupIso", ele_pileupIso, "ele_pileupIso[nElectrons]/F");
  llpTree->Branch("ele_chargedIso", ele_chargedIso, "ele_chargedIso[nElectrons]/F");
  llpTree->Branch("ele_photonIso", ele_photonIso, "ele_photonIso[nElectrons]/F");
  llpTree->Branch("ele_neutralHadIso", ele_neutralHadIso, "ele_neutralHadIso[nElectrons]/F");
  llpTree->Branch("ele_MissHits", ele_MissHits, "ele_MissHits[nElectrons]/I");
  llpTree->Branch("ele_PassConvVeto", ele_PassConvVeto, "ele_PassConvVeto[nElectrons]/O");
  llpTree->Branch("ele_OneOverEminusOneOverP", ele_OneOverEminusOneOverP, "ele_OneOverEminusOneOverP[nElectrons]/F");
  llpTree->Branch("ele_IDMVAGeneralPurpose", ele_IDMVAGeneralPurpose, "ele_IDMVAGeneralPurpose[nElectrons]/F");
  llpTree->Branch("ele_IDMVACategoryGeneralPurpose", ele_IDMVACategoryGeneralPurpose, "ele_IDMVACategoryGeneralPurpose[nElectrons]/I");
  llpTree->Branch("ele_IDMVAHZZ", ele_IDMVAHZZ, "ele_IDMVAHZZ[nElectrons]/F");
  llpTree->Branch("ele_IDMVACategoryHZZ", ele_IDMVACategoryHZZ, "ele_IDMVACategoryHZZ[nElectrons]/I");
  llpTree->Branch("ele_RegressionE", ele_RegressionE, "ele_RegressionE[nElectrons]/F");
  llpTree->Branch("ele_CombineP4", ele_CombineP4, "ele_CombineP4[nElectrons]/F");
  llpTree->Branch("ele_ptrel", ele_ptrel, "ele_ptrel[nElectrons]/F");
  llpTree->Branch("ele_chargedMiniIso", ele_chargedMiniIso, "ele_chargedMiniIso[nElectrons]/F");
  llpTree->Branch("ele_photonAndNeutralHadronMiniIso", ele_photonAndNeutralHadronMiniIso, "ele_photonAndNeutralHadronMiniIso[nElectrons]/F");
  llpTree->Branch("ele_chargedPileupMiniIso", ele_chargedPileupMiniIso, "ele_chargedPileupMiniIso[nElectrons]/F");
  llpTree->Branch("ele_activityMiniIsoAnnulus", ele_activityMiniIsoAnnulus, "ele_activityMiniIsoAnnulus[nElectrons]/F");
  llpTree->Branch("ele_passSingleEleTagFilter", ele_passSingleEleTagFilter, "ele_passSingleEleTagFilter[nElectrons]/O");
  llpTree->Branch("ele_passTPOneTagFilter", ele_passTPOneTagFilter, "ele_passTPOneTagFilter[nElectrons]/O");
  llpTree->Branch("ele_passTPTwoTagFilter", ele_passTPTwoTagFilter, "ele_passTPTwoTagFilter[nElectrons]/O");
  llpTree->Branch("ele_passTPOneProbeFilter", ele_passTPOneProbeFilter, "ele_passTPOneProbeFilter[nElectrons]/O");
  llpTree->Branch("ele_passTPTwoProbeFilter", ele_passTPTwoProbeFilter, "ele_passTPTwoProbeFilter[nElectrons]/O");
  llpTree->Branch("ele_passHLTFilter", &ele_passHLTFilter, Form("ele_passHLTFilter[nElectrons][%d]/O",MAX_ElectronHLTFilters));
  if (enableEcalRechits_) {
    ele_EcalRechitIndex = new std::vector<std::vector<uint> >; ele_EcalRechitIndex->clear();
    llpTree->Branch("ele_EcalRechitIndex", "std::vector<std::vector<uint> >",&ele_EcalRechitIndex);
    ele_SeedRechitIndex = new std::vector<uint>; ele_SeedRechitIndex->clear();
    llpTree->Branch("ele_SeedRechitIndex", "std::vector<uint>",&ele_SeedRechitIndex);
  }
};

void llp_ntupler::enableTauBranches()
{
  llpTree->Branch("nTaus", &nTaus,"nTaus/I");
  llpTree->Branch("tauE", tauE,"tauE[nTaus]/F");
  llpTree->Branch("tauPt", tauPt,"tauPt[nTaus]/F");
  llpTree->Branch("tauEta", tauEta,"tauEta[nTaus]/F");
  llpTree->Branch("tauPhi", tauPhi,"tauPhi[nTaus]/F");
  llpTree->Branch("tau_IsLoose", tau_IsLoose, "tau_IsLoose[nTaus]/O");
  llpTree->Branch("tau_IsMedium", tau_IsMedium, "tau_IsMedium[nTaus]/O");
  llpTree->Branch("tau_IsTight", tau_IsTight, "tau_IsTight[nTaus]/O");
  llpTree->Branch("tau_passEleVetoLoose", tau_passEleVetoLoose, "tau_passEleVetoLoose[nTaus]/O");
  llpTree->Branch("tau_passEleVetoMedium", tau_passEleVetoMedium, "tau_passEleVetoMedium[nTaus]/O");
  llpTree->Branch("tau_passEleVetoTight", tau_passEleVetoTight, "tau_passEleVetoTight[nTaus]/O");
  llpTree->Branch("tau_passMuVetoLoose", tau_passMuVetoLoose, "tau_passMuVetoLoose[nTaus]/O");
  llpTree->Branch("tau_passMuVetoMedium", tau_passMuVetoMedium, "tau_passMuVetoMedium[nTaus]/O");
  llpTree->Branch("tau_passMuVetoTight", tau_passMuVetoTight, "tau_passMuVetoTight[nTaus]/O");
  llpTree->Branch("tau_ID", tau_ID, "tau_ID[nTaus]/i");
  llpTree->Branch("tau_combinedIsoDeltaBetaCorr3Hits", tau_combinedIsoDeltaBetaCorr3Hits, "tau_combinedIsoDeltaBetaCorr3Hits[nTaus]/F");
  llpTree->Branch("tau_chargedIsoPtSum", tau_chargedIsoPtSum, "tau_chargedIsoPtSum[nTaus]/F");
  llpTree->Branch("tau_neutralIsoPtSum", tau_neutralIsoPtSum, "tau_neutralIsoPtSum[nTaus]/F");
  llpTree->Branch("tau_puCorrPtSum", tau_puCorrPtSum, "tau_puCorrPtSum[nTaus]/F");
  llpTree->Branch("tau_eleVetoMVA", tau_eleVetoMVA, "tau_eleVetoMVA[nTaus]/F");
  llpTree->Branch("tau_eleVetoCategory", tau_eleVetoCategory, "tau_eleVetoCategory[nTaus]/I");
  llpTree->Branch("tau_muonVetoMVA", tau_muonVetoMVA, "tau_muonVetoMVA[nTaus]/F");
  llpTree->Branch("tau_isoMVAnewDMwLT", tau_isoMVAnewDMwLT, "tau_isoMVAnewDMwLT[nTaus]/F");
  llpTree->Branch("tau_isoMVAnewDMwoLT", tau_isoMVAnewDMwoLT, "tau_isoMVAnewDMwoLT[nTaus]/F");
  llpTree->Branch("tau_leadCandPt", tau_leadCandPt, "tau_leadCandPt[nTaus]/F");
  llpTree->Branch("tau_leadCandID", tau_leadCandID, "tau_leadCandID[nTaus]/I");
  llpTree->Branch("tau_leadChargedHadrCandPt", tau_leadChargedHadrCandPt, "tau_leadChargedHadrCandPt[nTaus]/F");
  llpTree->Branch("tau_leadChargedHadrCandID", tau_leadChargedHadrCandID, "tau_leadChargedHadrCandID[nTaus]/I");
};

void llp_ntupler::enableIsoPFCandidateBranches()
{
  llpTree->Branch("nIsoPFCandidates", &nIsoPFCandidates, "nIsoPFCandidates/i");
  llpTree->Branch("isoPFCandidatePt", isoPFCandidatePt, "isoPFCandidatePt[nIsoPFCandidates]/F");
  llpTree->Branch("isoPFCandidateEta", isoPFCandidateEta, "isoPFCandidateEta[nIsoPFCandidates]/F");
  llpTree->Branch("isoPFCandidatePhi", isoPFCandidatePhi, "isoPFCandidatePhi[nIsoPFCandidates]/F");
  llpTree->Branch("isoPFCandidateIso04", isoPFCandidateIso04, "isoPFCandidateIso04[nIsoPFCandidates]/F");
  llpTree->Branch("isoPFCandidateD0", isoPFCandidateD0, "isoPFCandidateD0[nIsoPFCandidates]/F");
  llpTree->Branch("isoPFCandidatePdgId", isoPFCandidatePdgId, "isoPFCandidatePdgId[nIsoPFCandidates]/I");
};

void llp_ntupler::enablePhotonBranches()
{
  llpTree->Branch("nPhotons", &nPhotons,"nPhotons/I");
  llpTree->Branch("nPhotons_overlap", &nPhotons_overlap,"nPhotons_overlap/I");
  llpTree->Branch("phoE", phoE,"phoE[nPhotons]/F");
  llpTree->Branch("phoPt", phoPt,"phoPt[nPhotons]/F");
  llpTree->Branch("phoEta", phoEta,"phoEta[nPhotons]/F");
  llpTree->Branch("phoPhi", phoPhi,"phoPhi[nPhotons]/F");
  llpTree->Branch("phoSigmaIetaIeta", phoSigmaIetaIeta, "phoSigmaIetaIeta[nPhotons]/F");
  llpTree->Branch("phoFull5x5SigmaIetaIeta", phoFull5x5SigmaIetaIeta, "phoFull5x5SigmaIetaIeta[nPhotons]/F");
  llpTree->Branch("phoR9", phoR9, "phoR9[nPhotons]/F");
  llpTree->Branch("pho_sminor", pho_sminor, "pho_sminor[nPhotons]/F");
  llpTree->Branch("pho_smajor", pho_smajor, "pho_smajor[nPhotons]/F");
  llpTree->Branch("pho_HoverE", pho_HoverE, "pho_HoverE[nPhotons]/F");
  llpTree->Branch("pho_sumChargedHadronPtAllVertices", &pho_sumChargedHadronPtAllVertices,Form("pho_sumChargedHadronPtAllVertices[nPhotons][%d]/F",MAX_NPV));
  llpTree->Branch("pho_sumChargedHadronPt", &pho_sumChargedHadronPt, "pho_sumChargedHadronPt[nPhotons]/F");
  llpTree->Branch("pho_sumNeutralHadronEt", pho_sumNeutralHadronEt, "pho_sumNeutralHadronEt[nPhotons]/F");
  llpTree->Branch("pho_sumPhotonEt", pho_sumPhotonEt, "pho_sumPhotonEt[nPhotons]/F");
  llpTree->Branch("pho_ecalPFClusterIso", pho_ecalPFClusterIso, "pho_ecalPFClusterIso[nPhotons]/F");
  llpTree->Branch("pho_hcalPFClusterIso", pho_hcalPFClusterIso, "pho_hcalPFClusterIso[nPhotons]/F");
  llpTree->Branch("pho_trkSumPtHollowConeDR03", pho_trkSumPtHollowConeDR03, "pho_trkSumPtHollowConeDR03[nPhotons]/F");
  llpTree->Branch("pho_sumWorstVertexChargedHadronPt", pho_sumWorstVertexChargedHadronPt, "pho_sumWorstVertexChargedHadronPt[nPhotons]/F");
  llpTree->Branch("pho_pfIsoChargedHadronIso", pho_pfIsoChargedHadronIso, "pho_pfIsoChargedHadronIso[nPhotons]/F");
  llpTree->Branch("pho_pfIsoChargedHadronIsoWrongVtx", pho_pfIsoChargedHadronIsoWrongVtx, "pho_pfIsoChargedHadronIsoWrongVtx[nPhotons]/F");
  llpTree->Branch("pho_pfIsoNeutralHadronIso", pho_pfIsoNeutralHadronIso, "pho_pfIsoNeutralHadronIso[nPhotons]/F");
  llpTree->Branch("pho_pfIsoPhotonIso", pho_pfIsoPhotonIso, "pho_pfIsoPhotonIso[nPhotons]/F");
  llpTree->Branch("pho_pfIsoModFrixione", pho_pfIsoModFrixione, "pho_pfIsoModFrixione[nPhotons]/F");
  llpTree->Branch("pho_pfIsoSumPUPt", pho_pfIsoSumPUPt, "pho_pfIsoSumPUPt[nPhotons]/F");
  llpTree->Branch("pho_isConversion", pho_isConversion, "pho_isConversion[nPhotons]/O");
  llpTree->Branch("pho_passEleVeto", pho_passEleVeto, "pho_passEleVeto[nPhotons]/O");
  llpTree->Branch("pho_RegressionE", pho_RegressionE, "pho_RegressionE[nPhotons]/F");
  llpTree->Branch("pho_RegressionEUncertainty", pho_RegressionEUncertainty, "pho_RegressionEUncertainty[nPhotons]/F");
  llpTree->Branch("pho_IDMVA", pho_IDMVA, "pho_IDMVA[nPhotons]/F");
  llpTree->Branch("pho_superClusterEnergy", pho_superClusterEnergy, "pho_superClusterEnergy[nPhotons]/F");
  llpTree->Branch("pho_superClusterRawEnergy", pho_superClusterRawEnergy, "pho_superClusterRawEnergy[nPhotons]/F");
  llpTree->Branch("pho_superClusterEta", pho_superClusterEta, "pho_superClusterEta[nPhotons]/F");
  llpTree->Branch("pho_superClusterPhi", pho_superClusterPhi, "pho_superClusterPhi[nPhotons]/F");
  llpTree->Branch("pho_superClusterX", pho_superClusterX, "pho_superClusterX[nPhotons]/F");
  llpTree->Branch("pho_superClusterY", pho_superClusterY, "pho_superClusterY[nPhotons]/F");
  llpTree->Branch("pho_superClusterZ", pho_superClusterZ, "pho_superClusterZ[nPhotons]/F");
  llpTree->Branch("pho_hasPixelSeed", pho_hasPixelSeed, "pho_hasPixelSeed[nPhotons]/O");
  llpTree->Branch("pho_passHLTFilter", &pho_passHLTFilter, Form("pho_passHLTFilter[nPhotons][%d]/O",MAX_PhotonHLTFilters));
  llpTree->Branch("pho_convType", pho_convType, "pho_convType[nPhotons]/I");
  llpTree->Branch("pho_convTrkZ", pho_convTrkZ, "pho_convTrkZ[nPhotons]/F");
  llpTree->Branch("pho_convTrkClusZ", pho_convTrkClusZ, "pho_convTrkClusZ[nPhotons]/F");
  llpTree->Branch("pho_vtxSumPx", &pho_vtxSumPx,Form("pho_vtxSumPx[nPhotons][%d]/F",MAX_NPV));
  llpTree->Branch("pho_vtxSumPy", &pho_vtxSumPy,Form("pho_vtxSumPy[nPhotons][%d]/F",MAX_NPV));
  llpTree->Branch("pho_isStandardPhoton", pho_isStandardPhoton, "pho_isStandardPhoton[nPhotons]/O");
  llpTree->Branch("pho_seedRecHitSwitchToGain6", pho_seedRecHitSwitchToGain6, "pho_seedRecHitSwitchToGain6[nPhotons]/F");
  llpTree->Branch("pho_seedRecHitSwitchToGain1", pho_seedRecHitSwitchToGain1, "pho_seedRecHitSwitchToGain1[nPhotons]/F");
  llpTree->Branch("pho_anyRecHitSwitchToGain6", pho_anyRecHitSwitchToGain6, "pho_anyRecHitSwitchToGain6[nPhotons]/F");
  llpTree->Branch("pho_anyRecHitSwitchToGain1", pho_anyRecHitSwitchToGain1, "pho_anyRecHitSwitchToGain1[nPhotons]/F");
  if (enableEcalRechits_) {
    pho_EcalRechitIndex = new std::vector<std::vector<uint> >; pho_EcalRechitIndex->clear();
    llpTree->Branch("pho_EcalRechitIndex", "std::vector<std::vector<uint> >",&pho_EcalRechitIndex);
    pho_SeedRechitIndex = new std::vector<uint>; pho_SeedRechitIndex->clear();
    llpTree->Branch("pho_SeedRechitIndex", "std::vector<uint>",&pho_SeedRechitIndex);
  }

};

void llp_ntupler::enableEcalRechitBranches()
{
  ecalRechit_Eta = new std::vector<float>; ecalRechit_Eta->clear();
  ecalRechit_Phi = new std::vector<float>; ecalRechit_Phi->clear();
  ecalRechit_X = new std::vector<float>; ecalRechit_X->clear();
  ecalRechit_Y = new std::vector<float>; ecalRechit_Y->clear();
  ecalRechit_Z = new std::vector<float>; ecalRechit_Z->clear();
  ecalRechit_E = new std::vector<float>; ecalRechit_E->clear();
  ecalRechit_T = new std::vector<float>; ecalRechit_T->clear();
  ecalRechit_ID = new std::vector<uint>; ecalRechit_ID->clear();
  ecalRechit_FlagOOT = new std::vector<bool>; ecalRechit_FlagOOT->clear();
  ecalRechit_GainSwitch1 = new std::vector<bool>; ecalRechit_GainSwitch1->clear();
  ecalRechit_GainSwitch6 = new std::vector<bool>; ecalRechit_GainSwitch6->clear();
  ecalRechit_transpCorr = new std::vector<float>; ecalRechit_transpCorr->clear();

  llpTree->Branch("ecalRechit_Eta", "std::vector<float>",&ecalRechit_Eta);
  llpTree->Branch("ecalRechit_Phi", "std::vector<float>",&ecalRechit_Phi);
  llpTree->Branch("ecalRechit_X", "std::vector<float>",&ecalRechit_X);
  llpTree->Branch("ecalRechit_Y", "std::vector<float>",&ecalRechit_Y);
  llpTree->Branch("ecalRechit_Z", "std::vector<float>",&ecalRechit_Z);
  llpTree->Branch("ecalRechit_E", "std::vector<float>",&ecalRechit_E);
  llpTree->Branch("ecalRechit_T", "std::vector<float>",&ecalRechit_T);
  llpTree->Branch("ecalRechit_ID", "std::vector<uint>",&ecalRechit_ID);
  llpTree->Branch("ecalRechit_FlagOOT", "std::vector<bool>",&ecalRechit_FlagOOT);
  llpTree->Branch("ecalRechit_GainSwitch1", "std::vector<bool>",&ecalRechit_GainSwitch1);
  llpTree->Branch("ecalRechit_GainSwitch6", "std::vector<bool>",&ecalRechit_GainSwitch6);
  llpTree->Branch("ecalRechit_transpCorr", "std::vector<float>",&ecalRechit_transpCorr);

};

void llp_ntupler::enableJetBranches()
{
  llpTree->Branch("nJets", &nJets,"nJets/I");
  llpTree->Branch("jetE", jetE,"jetE[nJets]/F");
  llpTree->Branch("jetPt", jetPt,"jetPt[nJets]/F");
  llpTree->Branch("jetEta", jetEta,"jetEta[nJets]/F");
  llpTree->Branch("jetPhi", jetPhi,"jetPhi[nJets]/F");
  llpTree->Branch("jetCSV", jetCSV,"jetCSV[nJets]/F");
  llpTree->Branch("jetCISV", jetCISV,"jetCISV[nJets]/F");
  llpTree->Branch("jetProbb", jetProbb,"jetProbb[nJets]/F");
  llpTree->Branch("jetProbc", jetProbc,"jetProbc[nJets]/F");
  llpTree->Branch("jetProbudsg", jetProbudsg,"jetProbudsg[nJets]/F");
  llpTree->Branch("jetProbbb", jetProbbb,"jetProbbb[nJets]/F");
  llpTree->Branch("jetMass", jetMass, "jetMass[nJets]/F");
  llpTree->Branch("jetJetArea", jetJetArea, "jetJetArea[nJets]/F");
  llpTree->Branch("jetPileupE", jetPileupE, "jetPileupE[nJets]/F");
  llpTree->Branch("jetPileupId", jetPileupId, "jetPileupId[nJets]/F");
  llpTree->Branch("jetPileupIdFlag", jetPileupIdFlag, "jetPileupIdFlag[nJets]/I");
  llpTree->Branch("jetPassIDLoose", jetPassIDLoose, "jetPassIDLoose[nJets]/O");
  llpTree->Branch("jetPassIDTight", jetPassIDTight, "jetPassIDTight[nJets]/O");
  llpTree->Branch("jetPassMuFrac", jetPassMuFrac, "jetPassMuFrac[nJets]/O");
  llpTree->Branch("jetPassEleFrac", jetPassEleFrac, "jetPassEleFrac[nJets]/O");
  llpTree->Branch("jetPartonFlavor", jetPartonFlavor, "jetPartonFlavor[nJets]/I");
  llpTree->Branch("jetHadronFlavor", jetHadronFlavor, "jetHadronFlavor[nJets]/I");
  llpTree->Branch("jetChargedEMEnergyFraction", jetChargedEMEnergyFraction, "jetChargedEMEnergyFraction[nJets]/F");
  llpTree->Branch("jetNeutralEMEnergyFraction", jetNeutralEMEnergyFraction, "jetNeutralEMEnergyFraction[nJets]/F");
  llpTree->Branch("jetChargedHadronEnergyFraction", jetChargedHadronEnergyFraction, "jetChargedHadronEnergyFraction[nJets]/F");
  llpTree->Branch("jetNeutralHadronEnergyFraction", jetNeutralHadronEnergyFraction, "jetNeutralHadronEnergyFraction[nJets]/F");
  llpTree->Branch("jetMuonEnergyFraction", jetMuonEnergyFraction, "jetMuonEnergyFraction[nJets]/F");
  llpTree->Branch("jetHOEnergyFraction", jetHOEnergyFraction, "jetHOEnergyFraction[nJets]/F");
  llpTree->Branch("jetHFHadronEnergyFraction", jetHFHadronEnergyFraction, "jetHFHadronEnergyFraction[nJets]/F");
  llpTree->Branch("jetHFEMEnergyFraction",jetHFEMEnergyFraction, "jetHFEMEnergyFraction[nJets]/F");
  llpTree->Branch("jetAllMuonPt", jetAllMuonPt,"jetAllMuonPt[nJets]/F");
  llpTree->Branch("jetAllMuonEta", jetAllMuonEta,"jetAllMuonEta[nJets]/F");
  llpTree->Branch("jetAllMuonPhi", jetAllMuonPhi,"jetAllMuonPhi[nJets]/F");
  llpTree->Branch("jetAllMuonM", jetAllMuonM,"jetAllMuonM[nJets]/F");
  llpTree->Branch("jetPtWeightedDZ", jetPtWeightedDZ,"jetPtWeightedDZ[nJets]/F");
};

void llp_ntupler::enableJetAK8Branches()
{
  llpTree->Branch("nFatJets", &nFatJets,"nFatJets/i");
  llpTree->Branch("fatJetE", fatJetE,"fatJetE[nFatJets]/F");
  llpTree->Branch("fatJetPt", fatJetPt,"fatJetPt[nFatJets]/F");
  llpTree->Branch("fatJetEta", fatJetEta,"fatJetEta[nFatJets]/F");
  llpTree->Branch("fatJetPhi", fatJetPhi,"fatJetPhi[nFatJets]/F");
  llpTree->Branch("fatJetCorrectedPt", fatJetCorrectedPt,"fatJetCorrectedPt[nFatJets]/F");
  // llpTree->Branch("fatJetCorrectedEta", fatJetCorrectedEta,"fatJetCorrectedEta[nFatJets]/F");
  // llpTree->Branch("fatJetCorrectedPhi", fatJetCorrectedPhi,"fatJetCorrectedPhi[nFatJets]/F");
  llpTree->Branch("fatJetPrunedM", fatJetPrunedM,"fatJetPrunedM[nFatJets]/F");
  llpTree->Branch("fatJetTrimmedM", fatJetTrimmedM,"fatJetTrimmedM[nFatJets]/F");
  llpTree->Branch("fatJetFilteredM", fatJetFilteredM,"fatJetFilteredM[nFatJets]/F");
  llpTree->Branch("fatJetSoftDropM", fatJetSoftDropM,"fatJetSoftDropM[nFatJets]/F");
  llpTree->Branch("fatJetCorrectedSoftDropM", fatJetCorrectedSoftDropM,"fatJetCorrectedSoftDropM[nFatJets]/F");
  llpTree->Branch("fatJetUncorrectedSoftDropM", fatJetUncorrectedSoftDropM,"fatJetUncorrectedSoftDropM[nFatJets]/F");
  llpTree->Branch("fatJetTau1", fatJetTau1,"fatJetTau1[nFatJets]/F");
  llpTree->Branch("fatJetTau2", fatJetTau2,"fatJetTau2[nFatJets]/F");
  llpTree->Branch("fatJetTau3", fatJetTau3,"fatJetTau3[nFatJets]/F");
  llpTree->Branch("fatJetMaxSubjetCSV", fatJetMaxSubjetCSV, "fatJetMaxSubjetCSV[nFatJets]/F");
  llpTree->Branch("fatJetPassIDLoose", fatJetPassIDLoose,"fatJetPassIDLoose[nFatJets]/O");
  llpTree->Branch("fatJetPassIDTight", fatJetPassIDTight,"fatJetPassIDTight[nFatJets]/O");
};

void llp_ntupler::enableMetBranches()
{
  llpTree->Branch("metPt", &metPt, "metPt/F");
  llpTree->Branch("metPhi", &metPhi, "metPhi/F");
  llpTree->Branch("sumMET", &sumMET, "sumMET/F");
  llpTree->Branch("metType0Pt", &metType0Pt, "metType0Pt/F");
  llpTree->Branch("metType0Phi", &metType0Phi, "metType0Phi/F");
  llpTree->Branch("metType1Pt_raw", &metType1Pt_raw, "metType1Pt_raw/F");
  llpTree->Branch("metType1Pt", &metType1Pt, "metType1Pt/F");
  llpTree->Branch("metType1Px", &metType1Px, "metType1Px/F");
  llpTree->Branch("metType1Py", &metType1Py, "metType1Py/F");
  llpTree->Branch("metType1Eta", &metType1Eta, "metType1Eta/F");
  llpTree->Branch("metType1Phi", &metType1Phi, "metType1Phi/F");
  llpTree->Branch("metType1Phi_raw", &metType1Phi_raw, "metType1Phi_raw/F");
  llpTree->Branch("metType0Plus1Pt", &metType0Plus1Pt, "metType0Plus1Pt/F");
  llpTree->Branch("metType0Plus1Phi", &metType0Plus1Phi, "metType0Plus1Phi/F");
  llpTree->Branch("metNoHFPt", &metNoHFPt, "metNoHFPt/F");
  llpTree->Branch("metNoHFPhi", &metNoHFPhi, "metNoHFPhi/F");
  llpTree->Branch("metPuppiPt", &metPuppiPt, "metPuppiPt/F");
  llpTree->Branch("metPuppiPhi", &metPuppiPhi, "metPuppiPhi/F");
  llpTree->Branch("metCaloPt", &metCaloPt, "metCaloPt/F");
  llpTree->Branch("metCaloPhi", &metCaloPhi, "metCaloPhi/F");

  llpTree->Branch("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, "Flag_HBHENoiseFilter/O");
  llpTree->Branch("Flag_HBHETightNoiseFilter", &Flag_HBHETightNoiseFilter, "Flag_HBHETightNoiseFilter/O");
  llpTree->Branch("Flag_HBHEIsoNoiseFilter", &Flag_HBHEIsoNoiseFilter, "Flag_HBHEIsoNoiseFilter/O");
  llpTree->Branch("Flag_badChargedCandidateFilter", &Flag_badChargedCandidateFilter, "Flag_badChargedCandidateFilter/O");
  llpTree->Branch("Flag_badMuonFilter", &Flag_badMuonFilter, "Flag_badMuonFilter/O");
  llpTree->Branch("Flag_badGlobalMuonFilter", &Flag_badGlobalMuonFilter, "Flag_badGlobalMuonFilter/O");
  llpTree->Branch("Flag_duplicateMuonFilter", &Flag_duplicateMuonFilter, "Flag_duplicateMuonFilter/O");
  llpTree->Branch("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, "Flag_CSCTightHaloFilter/O");
  llpTree->Branch("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, "Flag_hcalLaserEventFilter/O");
  llpTree->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, "Flag_EcalDeadCellTriggerPrimitiveFilter/O");
  llpTree->Branch("Flag_EcalDeadCellBoundaryEnergyFilter", &Flag_EcalDeadCellBoundaryEnergyFilter, "Flag_EcalDeadCellBoundaryEnergyFilter/O");
  llpTree->Branch("Flag_goodVertices", &Flag_goodVertices, "Flag_goodVertices/O");
  llpTree->Branch("Flag_trackingFailureFilter", &Flag_trackingFailureFilter, "Flag_trackingFailureFilter/O");
  llpTree->Branch("Flag_eeBadScFilter", &Flag_eeBadScFilter, "Flag_eeBadScFilter/O");
  llpTree->Branch("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, "Flag_ecalLaserCorrFilter/O");
  llpTree->Branch("Flag_trkPOGFilters", &Flag_trkPOGFilters, "Flag_trkPOGFilters/O");
  llpTree->Branch("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, "Flag_trkPOG_manystripclus53X/O");
  llpTree->Branch("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, "Flag_trkPOG_toomanystripclus53X/O");
  llpTree->Branch("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, "Flag_trkPOG_logErrorTooManyClusters/O");
  llpTree->Branch("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, "Flag_BadPFMuonFilter/O");
  llpTree->Branch("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, "Flag_BadChargedCandidateFilter/O");
  llpTree->Branch("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, "Flag_ecalBadCalibFilter/O");
  llpTree->Branch("Flag_METFilters", &Flag_METFilters, "Flag_METFilters/O");
};

void llp_ntupler::enableTriggerBranches()
{
  nameHLT = new std::vector<std::string>; nameHLT->clear();
  llpTree->Branch("HLTDecision", &triggerDecision, ("HLTDecision[" + std::to_string(NTriggersMAX) +  "]/O").c_str());
  llpTree->Branch("HLTPrescale", &triggerHLTPrescale, ("HLTPrescale[" + std::to_string(NTriggersMAX) +  "]/I").c_str());
  llpTree->Branch("HLTMR", &HLTMR, "HLTMR/F");
  llpTree->Branch("HLTRSQ", &HLTRSQ, "HLTRSQ/F");
};

void llp_ntupler::enableMCBranches()
{
  llpTree->Branch("nGenJets", &nGenJets, "nGenJets/I");
  llpTree->Branch("genJetE", genJetE, "genJetE[nGenJets]/F");
  llpTree->Branch("genJetPt", genJetPt, "genJetPt[nGenJets]/F");
  llpTree->Branch("genJetEta", genJetEta, "genJetEta[nGenJets]/F");
  llpTree->Branch("genJetPhi", genJetPhi, "genJetPhi[nGenJets]/F");
  llpTree->Branch("genMetPtCalo", &genMetPtCalo, "genMetPtCalo/F");
  llpTree->Branch("genMetPhiCalo", &genMetPhiCalo, "genMetPhiCalo/F");
  llpTree->Branch("genMetPtTrue", &genMetPtTrue, "genMetPtTrue/F");
  llpTree->Branch("genMetPhiTrue", &genMetPhiTrue, "genMetPhiTrue/F");
  llpTree->Branch("genVertexX", &genVertexX, "genVertexX/F");
  llpTree->Branch("genVertexY", &genVertexY, "genVertexY/F");
  llpTree->Branch("genVertexZ", &genVertexZ, "genVertexZ/F");
  llpTree->Branch("genVertexT", &genVertexT, "genVertexT/F");
  llpTree->Branch("genWeight", &genWeight, "genWeight/F");
  llpTree->Branch("genSignalProcessID", &genSignalProcessID, "genSignalProcessID/i");
  llpTree->Branch("genQScale", &genQScale, "genQScale/F");
  llpTree->Branch("genAlphaQCD", &genAlphaQCD, "genAlphaQCD/F");
  llpTree->Branch("genAlphaQED", &genAlphaQED, "genAlphaQED/F");
  scaleWeights = new std::vector<float>; scaleWeights->clear();
  pdfWeights = new std::vector<float>; pdfWeights->clear();
  alphasWeights = new std::vector<float>; alphasWeights->clear();
  if (isFastsim_) {
    llpTree->Branch("lheComments", "std::string",&lheComments);
  }
  llpTree->Branch("scaleWeights", "std::vector<float>",&scaleWeights);
  llpTree->Branch("pdfWeights", "std::vector<float>",&pdfWeights);
  llpTree->Branch("alphasWeights", "std::vector<float>",&alphasWeights);
};

void llp_ntupler::enableGenParticleBranches()
{
  llpTree->Branch("nGenParticle", &nGenParticle, "nGenParticle/I");
  llpTree->Branch("gParticleMotherId", gParticleMotherId, "gParticleMotherId[nGenParticle]/I");
  llpTree->Branch("gParticleMotherIndex", gParticleMotherIndex, "gParticleMotherIndex[nGenParticle]/I");
  llpTree->Branch("gParticleId", gParticleId, "gParticleId[nGenParticle]/I");
  llpTree->Branch("gParticleStatus", gParticleStatus, "gParticleStatus[nGenParticle]/I");
  llpTree->Branch("gParticleE", gParticleE, "gParticleE[nGenParticle]/F");
  llpTree->Branch("gParticlePt", gParticlePt, "gParticlePt[nGenParticle]/F");
  llpTree->Branch("gParticlePx", gParticlePx, "gParticlePx[nGenParticle]/F");
  llpTree->Branch("gParticlePy", gParticlePy, "gParticlePy[nGenParticle]/F");
  llpTree->Branch("gParticlePz", gParticlePz, "gParticlePz[nGenParticle]/F");
  llpTree->Branch("gParticleEta", gParticleEta, "gParticleEta[nGenParticle]/F");
  llpTree->Branch("gParticlePhi", gParticlePhi, "gParticlePhi[nGenParticle]/F");
  llpTree->Branch("gParticleDecayVertexX", gParticleDecayVertexX, "gParticleDecayVertexX[nGenParticle]/F");
  llpTree->Branch("gParticleDecayVertexY", gParticleDecayVertexY, "gParticleDecayVertexY[nGenParticle]/F");
  llpTree->Branch("gParticleDecayVertexZ", gParticleDecayVertexZ, "gParticleDecayVertexZ[nGenParticle]/F");
};

//------ Enable the desired set of branches ------//

/*
void llp_ntupler::setBranches(){

  llpTree->Branch("isData", &isData, "isData/O");
  llpTree->Branch("isFourJet", &isFourJet, "isFourJet/O");
  llpTree->Branch("isQCD", &isQCD, "isQCD/O");
  llpTree->Branch("runNum", &runNum, "runNum/i");
  llpTree->Branch("lumiNum", &lumiNum, "lumiNum/i");
  llpTree->Branch("eventNum", &eventNum, "eventNum/i");
  llpTree->Branch("pvX", &pvX, "pvX/F");
  llpTree->Branch("pvY", &pvY, "pvY/F");
  llpTree->Branch("pvZ", &pvZ, "pvZ/F");
  llpTree->Branch("nPV", &nPV, "nPV/I");
  llpTree->Branch("Rho", &Rho, "Rho/F");
  llpTree->Branch("nPU", &nPU, "nPU/I");
  llpTree->Branch("nPUmean", &nPUmean, "nPUmean/F");

  llpTree->Branch("nJets", &nJets,"nJets/I");
  llpTree->Branch("jetE", jetE,"jetE[nJets]/F");
  llpTree->Branch("jetPt", jetPt,"jetPt[nJets]/F");
  llpTree->Branch("jetEta", jetEta,"jetEta[nJets]/F");
  llpTree->Branch("jetPhi", jetPhi,"jetPhi[nJets]/F");
  llpTree->Branch("jetCISV", jetCISV,"jetCISV[nJets]/F");
  llpTree->Branch("jetMass", jetMass, "jetMass[nJets]/F");
  llpTree->Branch("jetJetArea", jetJetArea, "jetJetArea[nJets]/F");
  llpTree->Branch("jetPileupE", jetPileupE, "jetPileupE[nJets]/F");
  llpTree->Branch("jetPileupId", jetPileupId, "jetPileupId[nJets]/F");
  llpTree->Branch("jetPileupIdFlag", jetPileupIdFlag, "jetPileupIdFlag[nJets]/I");
  llpTree->Branch("jetPassIDLoose", jetPassIDLoose, "jetPassIDLoose[nJets]/O");
  llpTree->Branch("jetPassIDTight", jetPassIDTight, "jetPassIDTight[nJets]/O");
  llpTree->Branch("jetPassMuFrac", jetPassMuFrac, "jetPassMuFrac[nJets]/O");
  llpTree->Branch("jetPassEleFrac", jetPassEleFrac, "jetPassEleFrac[nJets]/O");
  llpTree->Branch("jetPartonFlavor", jetPartonFlavor, "jetPartonFlavor[nJets]/I");
  llpTree->Branch("jetHadronFlavor", jetHadronFlavor, "jetHadronFlavor[nJets]/I");
  llpTree->Branch("jetChargedEMEnergyFraction", jetChargedEMEnergyFraction, "jetChargedEMEnergyFraction[nJets]/F");
  llpTree->Branch("jetNeutralEMEnergyFraction", jetNeutralEMEnergyFraction, "jetNeutralEMEnergyFraction[nJets]/F");
  llpTree->Branch("jetChargedHadronEnergyFraction", jetChargedHadronEnergyFraction, "jetChargedHadronEnergyFraction[nJets]/F");
  llpTree->Branch("jetNeutralHadronEnergyFraction", jetNeutralHadronEnergyFraction, "jetNeutralHadronEnergyFraction[nJets]/F");
  llpTree->Branch("jet_charged_hadron_multiplicity", jet_charged_hadron_multiplicity, "jet_charged_hadron_multiplicity[nJets]/F");
  llpTree->Branch("jet_neutral_hadron_multiplicity", jet_neutral_hadron_multiplicity, "jet_neutral_hadron_multiplicity[nJets]/F");
  llpTree->Branch("jet_photon_multiplicity", jet_photon_multiplicity, "jet_photon_multiplicity[nJets]/F");
  llpTree->Branch("jet_electron_multiplicity", jet_electron_multiplicity, "jet_electron_multiplicity[nJets]/F");
  llpTree->Branch("jet_muon_multiplicity", jet_muon_multiplicity, "jet_muon_multiplicity[nJets]/F");
  llpTree->Branch("jet_HF_hadron_multiplicity", jet_HF_hadron_multiplicity, "jet_HF_hadron_multiplicity[nJets]/F");
  llpTree->Branch("jet_HF_em_multiplicity", jet_HF_em_multiplicity, "jet_HF_em_multiplicity[nJets]/F");
  llpTree->Branch("jet_charged_multiplicity", jet_charged_multiplicity, "jet_charged_multiplicity[nJets]/F");
  llpTree->Branch("jet_neutral_multiplicity", jet_neutral_multiplicity, "jet_neutral_multiplicity[nJets]/F");
  llpTree->Branch("jetMatchedGenPt", jetMatchedGenPt,"jetMatchedGenPt[nJets]/F");
  llpTree->Branch("jetMatchedGenEta", jetMatchedGenEta,"jetMatchedGenEta[nJets]/F");
  llpTree->Branch("jetMatchedGenPhi", jetMatchedGenPhi,"jetMatchedGenPhi[nJets]/F");
  llpTree->Branch("jetMatchedGenMass", jetMatchedGenMass, "jetMatchedGenMass[nJets]/F");




  if( enableRecHitInfo_ )
  {

    llpTree->Branch("jet_n_rechits", jet_n_rechits, "jet_n_rechits[nJets]/I");
    llpTree->Branch("jet_rechit_E_Ecut3", jet_rechit_E_Ecut3, "jet_rechit_E_Ecut3[nJets]/F");
    llpTree->Branch("jet_rechit_T_Ecut3", jet_rechit_T_Ecut3, "jet_rechit_T_Ecut3[nJets]/F");
    llpTree->Branch("jet_rechit_E_Ecut4", jet_rechit_E_Ecut4, "jet_rechit_E_Ecut4[nJets]/F");
    llpTree->Branch("jet_rechit_T_Ecut4", jet_rechit_T_Ecut4, "jet_rechit_T_Ecut4[nJets]/F");
    llpTree->Branch("jet_rechit_E_Ecut2", jet_rechit_E_Ecut2, "jet_rechit_E_Ecut2[nJets]/F");
    llpTree->Branch("jet_rechit_T_Ecut2", jet_rechit_T_Ecut2, "jet_rechit_T_Ecut2[nJets]/F");
    llpTree->Branch("jet_rechit_E_Ecut1p5", jet_rechit_E_Ecut1p5, "jet_rechit_E_Ecut1p5[nJets]/F");
    llpTree->Branch("jet_rechit_T_Ecut1p5", jet_rechit_T_Ecut1p5, "jet_rechit_T_Ecut1p5[nJets]/F");
    llpTree->Branch("jet_rechit_E_Ecut1", jet_rechit_E_Ecut1, "jet_rechit_E_Ecut1[nJets]/F");
    llpTree->Branch("jet_rechit_T_Ecut1", jet_rechit_T_Ecut1, "jet_rechit_T_Ecut1[nJets]/F");
    llpTree->Branch("jet_rechit_E_Ecut0p5", jet_rechit_E_Ecut0p5, "jet_rechit_E_Ecut0p5[nJets]/F");
    llpTree->Branch("jet_rechit_T_Ecut0p5", jet_rechit_T_Ecut0p5, "jet_rechit_T_Ecut0p5[nJets]/F");
    llpTree->Branch("jet_rechit_E", jet_rechit_E, "jet_rechit_E[nJets]/F");
    llpTree->Branch("jet_rechit_T", jet_rechit_T, "jet_rechit_T[nJets]/F");


    llpTree->Branch("jet_pv_rechit_T_Ecut3", jet_pv_rechit_T_Ecut3, "jet_rechit_T_Ecut3[nJets]/F");
    llpTree->Branch("jet_pv_rechit_T_Ecut4", jet_pv_rechit_T_Ecut4, "jet_rechit_T_Ecut4[nJets]/F");
    llpTree->Branch("jet_pv_rechit_T_Ecut2", jet_pv_rechit_T_Ecut2, "jet_rechit_T_Ecut2[nJets]/F");
    llpTree->Branch("jet_pv_rechit_T_Ecut1p5", jet_pv_rechit_T_Ecut1p5, "jet_rechit_T_Ecut1p5[nJets]/F");
    llpTree->Branch("jet_pv_rechit_T_Ecut1", jet_pv_rechit_T_Ecut1, "jet_rechit_T_Ecut1[nJets]/F");
    llpTree->Branch("jet_pv_rechit_T_Ecut0p5", jet_pv_rechit_T_Ecut0p5, "jet_rechit_T_Ecut0p5[nJets]/F");
    llpTree->Branch("jet_pv_rechit_T", jet_pv_rechit_T, "jet_rechit_T[nJets]/F");

    llpTree->Branch("jet_rechits_E", jet_rechits_E, "jet_rechits_E[nJets][1000]/F");
    llpTree->Branch("jet_rechits_T", jet_rechits_T, "jet_rechits_T[nJets][1000]/F");
    llpTree->Branch("jet_pv_rechits_T", jet_pv_rechits_T, "jet_rechits_T[nJets][1000]/F");
  }

  llpTree->Branch("nPhotons", &fJetNPhotons,"nPhotons/I");
  llpTree->Branch("phoPt", fJetPhotonPt,"phoPt[nPhotons]/F");
  llpTree->Branch("phoEta", fJetPhotonEta,"phoEta[nPhotons]/F");
  llpTree->Branch("phoPhi", fJetPhotonPhi,"phoPhi[nPhotons]/F");
  llpTree->Branch("phoSeedRecHitEta", fJetPhotonSeedRecHitEta, "phoSeedRecHitEta[nPhotons]/F");
  llpTree->Branch("phoSeedRecHitPhi", fJetPhotonSeedRecHitPhi, "phoSeedRecHitPhi[nPhotons]/F");
  llpTree->Branch("phoSeedRecHitE", fJetPhotonSeedRecHitE, "phoSeedRecHitE[nPhotons]/F");
  llpTree->Branch("phoSeedRecHitT", fJetPhotonSeedRecHitTime, "phoSeedRecHitT[nPhotons]/F");

  // llpTree->Branch("fJetPhotonRecHitEta", "std::vector<float>",&fJetPhotonRecHitEta);
  // llpTree->Branch("fJetPhotonRecHitPhi", "std::vector<float>",&fJetPhotonRecHitPhi);
  // llpTree->Branch("fJetPhotonRecHitE", "std::vector<float>",&fJetPhotonRecHitE);
  // llpTree->Branch("fJetPhotonRecHitTime", "std::vector<float>",&fJetPhotonRecHitTime);

  cout << "BRANCHES\n";
  enableFatJetBranches();
  enableMCBranches();
  enableGenParticleBranches();
  if (enableTriggerInfo_) enableTriggerBranches();
  if (isQCD_)enableQCDBranches();
};

void llp_ntupler::enableFatJetBranches()
{
  llpTree->Branch("n_fat_Jets", &n_fat_Jets,"n_fat_Jets/I");
  llpTree->Branch("fat_jetE", fat_jetE,"fat_jetE[n_fat_Jets]/F");
  llpTree->Branch("fat_jetPt", fat_jetPt,"fat_jetPt[n_fat_Jets]/F");
  llpTree->Branch("fat_jetEta", fat_jetEta,"fat_jetEta[n_fat_Jets]/F");
  llpTree->Branch("fat_jetPhi", fat_jetPhi,"fat_jetPhi[n_fat_Jets]/F");
  llpTree->Branch("fat_jetCISV", fat_jetCISV,"fat_jetCISV[n_fat_Jets]/F");
  llpTree->Branch("fat_jetMass", fat_jetMass, "fat_jetMass[n_fat_Jets]/F");
  llpTree->Branch("fat_jetJetArea", fat_jetJetArea, "fat_jetJetArea[n_fat_Jets]/F");
  llpTree->Branch("fat_jetPileupE", fat_jetPileupE, "fat_jetPileupE[n_fat_Jets]/F");
  llpTree->Branch("fat_jetPileupId", fat_jetPileupId, "fat_jetPileupId[n_fat_Jets]/F");
  llpTree->Branch("fat_jetPileupIdFlag", fat_jetPileupIdFlag, "fat_jetPileupIdFlag[n_fat_Jets]/I");
  llpTree->Branch("fat_jetPassIDLoose", fat_jetPassIDLoose, "fat_jetPassIDLoose[n_fat_Jets]/O");
  llpTree->Branch("fat_jetPassIDTight", fat_jetPassIDTight, "fat_jetPassIDTight[n_fat_Jets]/O");
  llpTree->Branch("fat_jetPassMuFrac", fat_jetPassMuFrac, "fat_jetPassMuFrac[n_fat_Jets]/O");
  llpTree->Branch("fat_jetPassEleFrac", fat_jetPassEleFrac, "fat_jetPassEleFrac[n_fat_Jets]/O");
  llpTree->Branch("fat_jetPartonFlavor", fat_jetPartonFlavor, "fat_jetPartonFlavor[n_fat_Jets]/I");
  llpTree->Branch("fat_jetHadronFlavor", fat_jetHadronFlavor, "fat_jetHadronFlavor[n_fat_Jets]/I");
  llpTree->Branch("fat_jetChargedEMEnergyFraction", fat_jetChargedEMEnergyFraction, "fat_jetChargedEMEnergyFraction[n_fat_Jets]/F");
  llpTree->Branch("fat_jetNeutralEMEnergyFraction", fat_jetNeutralEMEnergyFraction, "fat_jetNeutralEMEnergyFraction[n_fat_Jets]/F");
  llpTree->Branch("fat_jetChargedHadronEnergyFraction", fat_jetChargedHadronEnergyFraction, "fat_jetChargedHadronEnergyFraction[n_fat_Jets]/F");
  llpTree->Branch("fat_jetNeutralHadronEnergyFraction", fat_jetNeutralHadronEnergyFraction, "fat_jetNeutralHadronEnergyFraction[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_charged_hadron_multiplicity", fat_jet_charged_hadron_multiplicity, "fat_jet_charged_hadron_multiplicity[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_neutral_hadron_multiplicity", fat_jet_neutral_hadron_multiplicity, "fat_jet_neutral_hadron_multiplicity[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_photon_multiplicity", fat_jet_photon_multiplicity, "fat_jet_photon_multiplicity[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_electron_multiplicity", fat_jet_electron_multiplicity, "fat_jet_electron_multiplicity[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_muon_multiplicity", fat_jet_muon_multiplicity, "fat_jet_muon_multiplicity[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_HF_hadron_multiplicity", fat_jet_HF_hadron_multiplicity, "fat_jet_HF_hadron_multiplicity[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_HF_em_multiplicity", fat_jet_HF_em_multiplicity, "fat_jet_HF_em_multiplicity[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_charged_multiplicity", fat_jet_charged_multiplicity, "fat_jet_charged_multiplicity[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_neutral_multiplicity", fat_jet_neutral_multiplicity, "fat_jet_neutral_multiplicity[n_fat_Jets]/F");
  llpTree->Branch("fat_jetMatchedGenPt", fat_jetMatchedGenPt,"fat_jetMatchedGenPt[n_fat_Jets]/F");
  llpTree->Branch("fat_jetMatchedGenEta", fat_jetMatchedGenEta,"fat_jetMatchedGenEta[n_fat_Jets]/F");
  llpTree->Branch("fat_jetMatchedGenPhi", fat_jetMatchedGenPhi,"fat_jetMatchedGenPhi[n_fat_Jets]/F");
  llpTree->Branch("fat_jetMatchedGenMass", fat_jetMatchedGenMass, "fat_jetMatchedGenMass[n_fat_Jets]/F");

  llpTree->Branch("fat_jet_rechit_E_Ecut3", fat_jet_rechit_E_Ecut3, "fat_jet_rechit_E_Ecut3[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_T_Ecut3", fat_jet_rechit_T_Ecut3, "fat_jet_rechit_T_Ecut3[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_E_Ecut4", fat_jet_rechit_E_Ecut4, "fat_jet_rechit_E_Ecut4[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_T_Ecut4", fat_jet_rechit_T_Ecut4, "fat_jet_rechit_T_Ecut4[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_E_Ecut2", fat_jet_rechit_E_Ecut2, "fat_jet_rechit_E_Ecut2[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_T_Ecut2", fat_jet_rechit_T_Ecut2, "fat_jet_rechit_T_Ecut2[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_E_Ecut1p5", fat_jet_rechit_E_Ecut1p5, "fat_jet_rechit_E_Ecut1p5[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_T_Ecut1p5", fat_jet_rechit_T_Ecut1p5, "fat_jet_rechit_T_Ecut1p5[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_E_Ecut1", fat_jet_rechit_E_Ecut1, "fat_jet_rechit_E_Ecut1[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_T_Ecut1", fat_jet_rechit_T_Ecut1, "fat_jet_rechit_T_Ecut1[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_E_Ecut0p5", fat_jet_rechit_E_Ecut0p5, "fat_jet_rechit_E_Ecut0p5[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_T_Ecut0p5", fat_jet_rechit_T_Ecut0p5, "fat_jet_rechit_T_Ecut0p5[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_E", fat_jet_rechit_E, "fat_jet_rechit_E[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_rechit_T", fat_jet_rechit_T, "fat_jet_rechit_T[n_fat_Jets]/F");

  llpTree->Branch("fat_jet_pv_rechit_T_Ecut3", fat_jet_pv_rechit_T_Ecut3, "fat_jet_rechit_T_Ecut3[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_pv_rechit_T_Ecut4", fat_jet_pv_rechit_T_Ecut4, "fat_jet_rechit_T_Ecut4[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_pv_rechit_T_Ecut2", fat_jet_pv_rechit_T_Ecut2, "fat_jet_rechit_T_Ecut2[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_pv_rechit_T_Ecut1p5", fat_jet_pv_rechit_T_Ecut1p5, "fat_jet_rechit_T_Ecut1p5[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_pv_rechit_T_Ecut1", fat_jet_pv_rechit_T_Ecut1, "fat_jet_rechit_T_Ecut1[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_pv_rechit_T_Ecut0p5", fat_jet_pv_rechit_T_Ecut0p5, "fat_jet_rechit_T_Ecut0p5[n_fat_Jets]/F");
  llpTree->Branch("fat_jet_pv_rechit_T", fat_jet_pv_rechit_T, "fat_jet_rechit_T[n_fat_Jets]/F");

  if( enableRecHitInfo_ )
  {
    llpTree->Branch("fat_jet_n_rechits", fat_jet_n_rechits, "fat_jet_n_rechits[n_fat_Jets]/I");
    llpTree->Branch("fat_jet_rechits_E", fat_jet_rechits_E, "fat_jet_rechits_E[n_fat_Jets][1000]/F");
    llpTree->Branch("fat_jet_rechits_T", fat_jet_rechits_T, "fat_jet_rechits_T[n_fat_Jets][1000]/F");
    llpTree->Branch("fat_jet_pv_rechits_T", fat_jet_pv_rechits_T, "fat_jet_rechits_T[n_fat_Jets][1000]/F");
  }


  return;
};

void llp_ntupler::enableMCBranches(){
  llpTree->Branch("nGenJets", &nGenJets, "nGenJets/I");
  llpTree->Branch("genJetE", genJetE, "genJetE[nGenJets]/F");
  llpTree->Branch("genJetPt", genJetPt, "genJetPt[nGenJets]/F");
  llpTree->Branch("genJetEta", genJetEta, "genJetEta[nGenJets]/F");
  llpTree->Branch("genJetPhi", genJetPhi, "genJetPhi[nGenJets]/F");
  llpTree->Branch("genJetME", genJetME, "genJetME[nGenJets]/F");
  llpTree->Branch("genMetPtCalo", &genMetPtCalo, "genMetPtCalo/F");
  llpTree->Branch("genMetPhiCalo", &genMetPhiCalo, "genMetPhiCalo/F");
  llpTree->Branch("genMetPtTrue", &genMetPtTrue, "genMetPtTrue/F");
  llpTree->Branch("genMetPhiTrue", &genMetPhiTrue, "genMetPhiTrue/F");
  llpTree->Branch("genVertexX", &genVertexX, "genVertexX/F");
  llpTree->Branch("genVertexY", &genVertexY, "genVertexY/F");
  llpTree->Branch("genVertexZ", &genVertexZ, "genVertexZ/F");
  llpTree->Branch("genVertexT", &genVertexT, "genVertexT/F");
  llpTree->Branch("genWeight", &genWeight, "genWeight/F");
  llpTree->Branch("genSignalProcessID", &genSignalProcessID, "genSignalProcessID/i");
  llpTree->Branch("genQScale", &genQScale, "genQScale/F");
  llpTree->Branch("genAlphaQCD", &genAlphaQCD, "genAlphaQCD/F");
  llpTree->Branch("genAlphaQED", &genAlphaQED, "genAlphaQED/F");
  llpTree->Branch("genJet_match_jet_index", &genJet_match_jet_index, "genJet_match_jet_index[nGenJets]/i");
  llpTree->Branch("genJet_min_delta_r_match_jet", &genJet_min_delta_r_match_jet, "genJet_min_delta_r_match_jet[nGenJets]/F");

  scaleWeights = new std::vector<float>; scaleWeights->clear();
  pdfWeights = new std::vector<float>; pdfWeights->clear();
  alphasWeights = new std::vector<float>; alphasWeights->clear();
  if (isFastsim_) {
    llpTree->Branch("lheComments", "std::string",&lheComments);
  }
  llpTree->Branch("scaleWeights", "std::vector<float>",&scaleWeights);
  llpTree->Branch("pdfWeights", "std::vector<float>",&pdfWeights);
  llpTree->Branch("alphasWeights", "std::vector<float>",&alphasWeights);

};
void llp_ntupler::enableQCDBranches()
{
  //QCD BRANCHES
  llpTree->Branch("nGenQCDParticles", &nGenQCDParticles, "nGenQCDParticles/I");
  llpTree->Branch("genQCD_e", genQCD_e, "genQCD_e[nGenQCDParticles]/F");
  llpTree->Branch("genQCD_pt", genQCD_pt, "genQCD_pt[nGenQCDParticles]/F");
  llpTree->Branch("genQCD_eta", genQCD_eta, "genQCD_eta[nGenQCDParticles]/F");
  llpTree->Branch("genQCD_phi", genQCD_phi, "genQCD_phi[nGenQCDParticles]/F");
  llpTree->Branch("genParticleQCD_match_jet_index", &genParticleQCD_match_jet_index, "genParticleQCD_match_jet_index[nGenQCDParticles]/i");
  llpTree->Branch("genParticleQCD_min_delta_r_match_jet", &genParticleQCD_min_delta_r_match_jet, "genParticleQCD_min_delta_r_match_jet[nGenQCDParticles]/F");
};
void llp_ntupler::enableTriggerBranches()
{
  nameHLT = new std::vector<std::string>; nameHLT->clear();
  llpTree->Branch("HLTDecision", &triggerDecision, ("HLTDecision[" + std::to_string(NTriggersMAX) +  "]/O").c_str());
  //llpTree->Branch("HLTPrescale", &triggerHLTPrescale, ("HLTPrescale[" + std::to_string(NTriggersMAX) +  "]/I").c_str());
  //llpTree->Branch("HLTMR", &HLTMR, "HLTMR/F");
  //llpTree->Branch("HLTRSQ", &HLTRSQ, "HLTRSQ/F");
};

void llp_ntupler::enableGenParticleBranches(){
  llpTree->Branch("gLLP_prod_vertex_x", gLLP_prod_vertex_x, "gLLP_prod_vertex_x[2]/F");
  llpTree->Branch("gLLP_prod_vertex_y", gLLP_prod_vertex_y, "gLLP_prod_vertex_y[2]/F");
  llpTree->Branch("gLLP_prod_vertex_z", gLLP_prod_vertex_z, "gLLP_prod_vertex_z[2]/F");
  llpTree->Branch("gLLP_decay_vertex_x", gLLP_decay_vertex_x, "gLLP_decay_vertex_x[2]/F");
  llpTree->Branch("gLLP_decay_vertex_y", gLLP_decay_vertex_y, "gLLP_decay_vertex_y[2]/F");
  llpTree->Branch("gLLP_decay_vertex_z", gLLP_decay_vertex_z, "gLLP_decay_vertex_z[2]/F");
  llpTree->Branch("gLLP_beta", gLLP_beta, "gLLP_beta[2]/F");
  llpTree->Branch("gLLP_travel_time", gLLP_travel_time, "gLLP_travel_time[2]/F");

  llpTree->Branch("gLLP_daughter_travel_time", gLLP_daughter_travel_time, "gLLP_daughter_travel_time[4]/F");
  llpTree->Branch("gLLP_daughter_pt", gLLP_daughter_pt, "gLLP_daughter_pt[4]/F");
  llpTree->Branch("gLLP_daughter_eta", gLLP_daughter_eta, "gLLP_daughter_eta[4]/F");
  llpTree->Branch("gLLP_daughter_phi", gLLP_daughter_phi, "gLLP_daughter_phi[4]/F");
  llpTree->Branch("gLLP_daughter_eta_ecalcorr", gLLP_daughter_eta_ecalcorr, "gLLP_daughter_eta_ecalcorr[4]/F");
  llpTree->Branch("gLLP_daughter_phi_ecalcorr", gLLP_daughter_phi_ecalcorr, "gLLP_daughter_phi_ecalcorr[4]/F");
  llpTree->Branch("gLLP_daughter_eta_hcalcorr", gLLP_daughter_eta_hcalcorr, "gLLP_daughter_eta_hcalcorr[4]/F");
  llpTree->Branch("gLLP_daughter_phi_hcalcorr", gLLP_daughter_phi_hcalcorr, "gLLP_daughter_phi_hcalcorr[4]/F");
  llpTree->Branch("gLLP_daughter_e", gLLP_daughter_e, "gLLP_daughter_e[4]/F");
  llpTree->Branch("photon_travel_time", photon_travel_time, "photon_travel_time[4]/F");
  llpTree->Branch("gen_time", gen_time, "gen_time[4]/F");
  llpTree->Branch("gen_time_pv", gen_time_pv, "gen_time_pv[4]/F");


  llpTree->Branch("gLLP_daughter_match_genJet_index", gLLP_daughter_match_genJet_index, "gLLP_daughter_match_genJet_index[4]/i");
  llpTree->Branch("gLLP_min_delta_r_match_genJet", gLLP_min_delta_r_match_genJet, "gLLP_min_delta_r_match_genJet[4]/F");
  llpTree->Branch("gLLP_daughter_match_jet_index_hcal", gLLP_daughter_match_jet_index_hcal, "gLLP_daughter_match_jet_index_hcal[4]/i");
  llpTree->Branch("gLLP_min_delta_r_match_jet_hcal", gLLP_min_delta_r_match_jet_hcal, "gLLP_min_delta_r_match_jet_hcal[4]/F");
  llpTree->Branch("gLLP_daughter_match_jet_index_hcal_loose", gLLP_daughter_match_jet_index_hcal_loose, "gLLP_daughter_match_jet_index_hcal_loose[4]/i");
  llpTree->Branch("gLLP_min_delta_r_match_jet_hcal_loose", gLLP_min_delta_r_match_jet_hcal_loose, "gLLP_min_delta_r_match_jet_hcal_loose[4]/F");
  llpTree->Branch("gLLP_daughter_match_jet_index_loose", gLLP_daughter_match_jet_index_loose, "gLLP_daughter_match_jet_index_loose[4]/i");
  llpTree->Branch("gLLP_min_delta_r_match_jet_loose", gLLP_min_delta_r_match_jet_loose, "gLLP_min_delta_r_match_jet_loose[4]/F");
  llpTree->Branch("gLLP_daughter_match_jet_index", gLLP_daughter_match_jet_index, "gLLP_daughter_match_jet_index[4]/i");
  llpTree->Branch("gLLP_min_delta_r_match_jet", gLLP_min_delta_r_match_jet, "gLLP_min_delta_r_match_jet[4]/F");
  llpTree->Branch("gLLP_min_delta_r_nocorr_match_jet", gLLP_min_delta_r_nocorr_match_jet, "gLLP_min_delta_r_nocorr_match_jet[4]/F");


  llpTree->Branch("nGenParticle", &nGenParticle, "nGenParticle/I");
  llpTree->Branch("gParticleMotherId", gParticleMotherId, "gParticleMotherId[nGenParticle]/I");
  llpTree->Branch("gParticleMotherIndex", gParticleMotherIndex, "gParticleMotherIndex[nGenParticle]/I");
  llpTree->Branch("gParticleId", gParticleId, "gParticleId[nGenParticle]/I");
  llpTree->Branch("gParticleStatus", gParticleStatus, "gParticleStatus[nGenParticle]/I");
  llpTree->Branch("gParticleE", gParticleE, "gParticleE[nGenParticle]/F");
  llpTree->Branch("gParticlePt", gParticlePt, "gParticlePt[nGenParticle]/F");
  llpTree->Branch("gParticlePx", gParticlePx, "gParticlePx[nGenParticle]/F");
  llpTree->Branch("gParticlePy", gParticlePy, "gParticlePy[nGenParticle]/F");
  llpTree->Branch("gParticlePz", gParticlePz, "gParticlePz[nGenParticle]/F");
  llpTree->Branch("gParticleEta", gParticleEta, "gParticleEta[nGenParticle]/F");
  llpTree->Branch("gParticlePhi", gParticlePhi, "gParticlePhi[nGenParticle]/F");
  llpTree->Branch("gParticleDecayVertexX", gParticleDecayVertexX, "gParticleDecayVertexX[nGenParticle]/F");
  llpTree->Branch("gParticleDecayVertexY", gParticleDecayVertexY, "gParticleDecayVertexY[nGenParticle]/F");
  llpTree->Branch("gParticleDecayVertexZ", gParticleDecayVertexZ, "gParticleDecayVertexZ[nGenParticle]/F");
}

*/

//------ Load the miniAOD objects and reset tree variables for each event ------//
void llp_ntupler::loadEvent(const edm::Event& iEvent){//load all miniAOD objects for the current event
  iEvent.getByToken(triggerBitsToken_, triggerBits);
  iEvent.getByToken(hepMCToken_, hepMC);
  iEvent.getByToken(triggerBitsToken_, triggerBits);
  iEvent.getByToken(metFilterBitsToken_, metFilterBits);
  iEvent.getByToken(verticesToken_, vertices);
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
//  iEvent.getByToken(superClustersToken_,superClusters);
//  iEvent.getByToken(lostTracksToken_,lostTracks);
//  iEvent.getByToken(hbheNoiseFilterToken_, hbheNoiseFilter);
//  iEvent.getByToken(hbheTightNoiseFilterToken_, hbheTightNoiseFilter);
//  iEvent.getByToken(hbheIsoNoiseFilterToken_, hbheIsoNoiseFilter);
  //iEvent.getByToken(badChargedCandidateFilterToken_, badChargedCandidateFilter);
  //iEvent.getByToken(badMuonFilterToken_, badMuonFilter);
  if(readGenVertexTime_) iEvent.getByToken(genParticles_t0_Token_,genParticles_t0);
  if (useGen_) {
//    iEvent.getByToken(genParticlesToken_,genParticles);
    iEvent.getByToken(genParticlesToken_,genParticles);
    iEvent.getByToken(genJetsToken_,genJets);

    //for Spring16 fastsim, this has been changed and removed
//    if (!isFastsim_) iEvent.getByToken(lheInfoToken_, lheInfo);

    iEvent.getByToken(genInfoToken_,genInfo);
    iEvent.getByToken(puInfoToken_,puInfo);
  }


}

//called by the loadEvent() method
void llp_ntupler::resetBranches(){
    //reset tree variables
    resetEventInfoBranches();
    resetPhotonBranches();
    resetJetBranches();
    resetGenParticleBranches();
    resetMCBranches();
    //reset_qcd_variables();
}

void llp_ntupler::resetEventInfoBranches()
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

void llp_ntupler::resetPhotonBranches()
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

void llp_ntupler::resetJetBranches()
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
  }
  return;
};

void llp_ntupler::resetGenParticleBranches()
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

    gParticleDecayVertexX[i] = -99999.0;
    gParticleDecayVertexY[i] = -99999.0;
    gParticleDecayVertexZ[i] = -99999.0;

  }
  return;
};

void llp_ntupler::resetMCBranches()
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

/*
void llp_ntupler::reset_qcd_variables()
{
  nGenQCDParticles = 0;
  for (int i = 0; i < GENPARTICLEARRAYSIZE; i++ )
  {
    genQCD_pt[i] = -666;
    genQCD_e[i] = -666;
    genQCD_eta[i] = -666;
    genQCD_phi[i] = -666;
    genParticleQCD_match_jet_index[i] = 666;
    genParticleQCD_min_delta_r_match_jet[i] = -666.;
  }
  return;
};
*/
//------ Methods to fill tree variables ------//

//------ Method called for each run ------//

void llp_ntupler::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {


}


//------ Method called for each lumi block ------//
void llp_ntupler::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {

}


//------ Method called for each event ------//

void llp_ntupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;

  //initialize
  loadEvent(iEvent); //loads objects and resets tree branches
  NEvents->Fill(0); //increment event count
  //resetting output tree branches
  resetBranches();
  //*************************************
  //Fill Event-Level Info
  //*************************************

  //store basic event info
  isData = isData_;
  //isFourJet = isFourJet_;
  //isQCD = isQCD_;
  runNum = iEvent.id().run();
  lumiNum = iEvent.luminosityBlock();
  eventNum = iEvent.id().event();

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
  //Rho = *rhoFastjetAll;

  //Fill Pileup info
  /*
  if (!isData)
  {
    for(const PileupSummaryInfo &pu : *puInfo)
    {
      if ( pu.getBunchCrossing() == 0)
      {
        nPU = pu.getPU_NumInteractions();
        nPUmean = pu.getTrueNumInteractions();
      }
    }
  }
*/
  int i_jet = 0;
  for (const reco::PFJet &j : *jets)
  {
    //resetBranches();
    if (j.pt() < 20) continue;
    if (fabs(j.eta()) > 2.4) continue;
    //*************************************
    //Fill Jet-Level Info
    //*************************************
    jetE[i_jet] = j.energy();
    jetPt[i_jet] = j.pt();
    jetEta[i_jet] = j.eta();
    jetPhi[i_jet] = j.phi();
    jetMass[i_jet] = j.mass();

    TLorentzVector thisJet;
    thisJet.SetPtEtaPhiE(jetPt[i_jet], jetEta[i_jet], jetPhi[i_jet], jetE[i_jet]);
    //jetCISV = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

    jetJetArea[i_jet] = j.jetArea();
    jetPileupE[i_jet] = j.pileup();

    jetPileupIdFlag[i_jet] = 0;
    jetPassIDLoose[i_jet] = passJetID(&j, 0);
    jetPassIDTight[i_jet] = passJetID(&j, 1);
    jetPassMuFrac[i_jet]  = ( j.muonEnergyFraction() < 0.80 );
    jetPassEleFrac[i_jet]  = ( j.electronEnergyFraction() < 0.90 );


    // if (useGen_) {
    //   jetPartonFlavor = j.partonFlavour();
    //   jetHadronFlavor = j.hadronFlavour();
    // }

    jetChargedEMEnergyFraction[i_jet] = j.chargedEmEnergyFraction();
    jetNeutralEMEnergyFraction[i_jet] = j.neutralEmEnergyFraction();
    jetChargedHadronEnergyFraction[i_jet] = j.chargedHadronEnergyFraction();
    jetNeutralHadronEnergyFraction[i_jet] = j.neutralHadronEnergyFraction();
    //jet_charged_hadron_multiplicity[i_jet] = j.chargedHadronMultiplicity();
    //jet_neutral_hadron_multiplicity[i_jet] = j.neutralHadronMultiplicity();
    //jet_photon_multiplicity[i_jet] = j.photonMultiplicity();
    //jet_electron_multiplicity[i_jet] = j.electronMultiplicity();
    //jet_muon_multiplicity[i_jet] = j.muonMultiplicity();
    //jet_HF_hadron_multiplicity[i_jet] = j.HFHadronMultiplicity();
    //jet_HF_em_multiplicity[i_jet] = j.HFEMMultiplicity();
    //jet_charged_multiplicity[i_jet] = j.chargedMultiplicity();
    //jet_neutral_multiplicity[i_jet] = j.neutralMultiplicity();


    //*************************************
    //find photons inside the jet
    //*************************************
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

      //*************************************
      //fill all rechits inside photons
      //*************************************

      fJetNPhotons++;

    }
    */

    //***************************
    //Find RecHits Inside the Jet
    //***************************
    // geometry (from ECAL ELF)

    edm::ESHandle<CaloGeometry> geoHandle;
    iSetup.get<CaloGeometryRecord>().get(geoHandle);
    const CaloSubdetectorGeometry *barrelGeometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
    //const CaloSubdetectorGeometry *endcapGeometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);
    double ecal_radius = 129.0;
    int n_matched_rechits = 0;
    for (EcalRecHitCollection::const_iterator recHit = ebRecHits->begin(); recHit != ebRecHits->end(); ++recHit)
    {
      if ( recHit->checkFlag(0) )
      {
        const DetId recHitId = recHit->detid();
        const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
        if ( deltaR(jetEta[i_jet], jetPhi[i_jet], recHitPos.eta(), recHitPos.phi())  < 0.4)
        {
          double rechit_x = ecal_radius * cos(recHitPos.phi());
    	    double rechit_y = ecal_radius * sin(recHitPos.phi());
    	    double rechit_z = ecal_radius * sinh(recHitPos.eta());
    	    double photon_pv_travel_time = (1./30) * sqrt(pow(pvX-rechit_x,2)+pow(pvY-rechit_y,2)+pow(pvZ-rechit_z,2));

          if (recHit->energy() > 1.0)
    	    {
    		    jetRechitE[i_jet] += recHit->energy();
    		    jetRechitTime[i_jet] += recHit->time()*recHit->energy();
          }
          n_matched_rechits++;
        }
      }
    }
    //cout << "Last Nphoton: " << fJetNPhotons << "\n";
    //std::cout << "n: " << n_matched_rechits << std::endl;
    jetNRechits[i_jet] = n_matched_rechits;
    jetRechitTime[i_jet] = jetRechitTime[i_jet]/jetRechitE[i_jet];
    nJets++;
    i_jet++;

  } //loop over jets

  //MC AND GEN LEVEL INFO
  fillMC();
  fillGenParticles();
  //fill_fat_jet( iSetup );
  /*if(readGenVertexTime_)
  {
    genVertexT = *genParticles_t0; //std::cout << genVertexT << std::endl;
  }
  */
  //fillMC();

  if ( enableTriggerInfo_ ) fillTrigger( iEvent );
  llpTree->Fill();
}

//------ Method called once each job just before starting event loop ------//
void llp_ntupler::beginJob(){
  setBranches();
}

//------ Method called once each job just after ending the event loop ------//
void llp_ntupler::endJob(){
}

/*
bool llp_ntupler::fill_fat_jet(const edm::EventSetup& iSetup)
{
  int i_fat_jet = 0;
  for (const reco::PFJet &j : *jetsAK8)
  {
    //resetBranches();
    if (j.pt() < 20) continue;
    if (fabs(j.eta()) > 2.4) continue;
    //*************************************
    //Fill Jet-Level Info
    //*************************************
    fat_jetE[i_fat_jet] = j.energy();
    fat_jetPt[i_fat_jet] = j.pt();
    fat_jetEta[i_fat_jet] = j.eta();
    fat_jetPhi[i_fat_jet] = j.phi();
    fat_jetMass[i_fat_jet] = j.mass();

    TLorentzVector thisJet;
    thisJet.SetPtEtaPhiE(fat_jetPt[i_fat_jet], fat_jetEta[i_fat_jet], fat_jetPhi[i_fat_jet], fat_jetE[i_fat_jet]);
    //jetCISV = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

    fat_jetJetArea[i_fat_jet] = j.jetArea();
    fat_jetPileupE[i_fat_jet] = j.pileup();

    fat_jetPileupIdFlag[i_fat_jet] = 0;
    fat_jetPassIDLoose[i_fat_jet] = passJetID(&j, 0);
    fat_jetPassIDTight[i_fat_jet] = passJetID(&j, 1);
    fat_jetPassMuFrac[i_fat_jet]  = ( j.muonEnergyFraction() < 0.80 );
    fat_jetPassEleFrac[i_fat_jet]  = ( j.electronEnergyFraction() < 0.90 );


    // if (useGen_) {
    //   fat_jetPartonFlavor = j.partonFlavour();
    //   fat_jetHadronFlavor = j.hadronFlavour();
    // }

    fat_jetChargedEMEnergyFraction[i_fat_jet] = j.chargedEmEnergyFraction();
    fat_jetNeutralEMEnergyFraction[i_fat_jet] = j.neutralEmEnergyFraction();
    fat_jetChargedHadronEnergyFraction[i_fat_jet] = j.chargedHadronEnergyFraction();
    fat_jetNeutralHadronEnergyFraction[i_fat_jet] = j.neutralHadronEnergyFraction();
    fat_jet_charged_hadron_multiplicity[i_fat_jet] = j.chargedHadronMultiplicity();
    fat_jet_neutral_hadron_multiplicity[i_fat_jet] = j.neutralHadronMultiplicity();
    fat_jet_photon_multiplicity[i_fat_jet] = j.photonMultiplicity();
    fat_jet_electron_multiplicity[i_fat_jet] = j.electronMultiplicity();
    fat_jet_muon_multiplicity[i_fat_jet] = j.muonMultiplicity();
    fat_jet_HF_hadron_multiplicity[i_fat_jet] = j.HFHadronMultiplicity();
    fat_jet_HF_em_multiplicity[i_fat_jet] = j.HFEMMultiplicity();
    fat_jet_charged_multiplicity[i_fat_jet] = j.chargedMultiplicity();
    fat_jet_neutral_multiplicity[i_fat_jet] = j.neutralMultiplicity();

    //***************************
    //Find RecHits Inside the Jet
    //***************************
    // geometry (from ECAL ELF)

    edm::ESHandle<CaloGeometry> geoHandle;
    iSetup.get<CaloGeometryRecord>().get(geoHandle);
    const CaloSubdetectorGeometry *barrelGeometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
    //const CaloSubdetectorGeometry *endcapGeometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);
    double ecal_radius = 129.0;
    int n_matched_rechits = 0;
    for (EcalRecHitCollection::const_iterator recHit = ebRecHits->begin(); recHit != ebRecHits->end(); ++recHit)
    {
      if ( recHit->checkFlag(0) )
      {
        const DetId recHitId = recHit->detid();
        const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
        if ( deltaR(fat_jetEta[i_fat_jet], fat_jetPhi[i_fat_jet], recHitPos.eta(), recHitPos.phi())  < 0.4)
        {
          fat_jet_rechit_E[i_fat_jet] += recHit->energy();
          fat_jet_rechit_T[i_fat_jet] += recHit->time()*recHit->energy();
          fat_jet_rechits_E[i_fat_jet][n_matched_rechits] = recHit->energy();
          fat_jet_rechits_T[i_fat_jet][n_matched_rechits] = recHit->time();
          double rechit_x = ecal_radius * cos(recHitPos.phi());
          double rechit_y = ecal_radius * sin(recHitPos.phi());
          double rechit_z = ecal_radius * sinh(recHitPos.eta());
          double photon_pv_travel_time = (1./30) * sqrt(pow(pvX-rechit_x,2)+pow(pvY-rechit_y,2)+pow(pvZ-rechit_z,2));
          fat_jet_pv_rechits_T[i_fat_jet][n_matched_rechits] = recHit->time()+(1./30)*ecal_radius*cosh(recHitPos.eta()) - photon_pv_travel_time;
          fat_jet_pv_rechit_T[i_fat_jet] += recHit->energy()*fat_jet_pv_rechits_T[i_fat_jet][n_matched_rechits];
          // std::cout << fat_jet_pv_rechits_T[i_fat_jet][n_matched_rechits] << fat_jet_rechits_T[i_fat_jet][n_matched_rechits] << std::endl;
          if (recHit->energy() > 0.5)
          {
            fat_jet_rechit_E_Ecut0p5[i_fat_jet] += recHit->energy();
            fat_jet_rechit_T_Ecut0p5[i_fat_jet] += recHit->time()*recHit->energy();
            fat_jet_pv_rechit_T_Ecut0p5[i_fat_jet] += fat_jet_pv_rechits_T[i_fat_jet][n_matched_rechits] *recHit->energy();
          }
          if (recHit->energy() > 1.0)
          {
            fat_jet_rechit_E_Ecut1[i_fat_jet] += recHit->energy();
            fat_jet_rechit_T_Ecut1[i_fat_jet] += recHit->time()*recHit->energy();
            fat_jet_pv_rechit_T_Ecut1[i_fat_jet] += fat_jet_pv_rechits_T[i_fat_jet][n_matched_rechits] *recHit->energy();
            // std::cout << "rechit time, with pv"<<fat_jet_rechit_T_Ecut1[i_fat_jet]<< fat_jet_pv_rechit_T_Ecut1[i_fat_jet]<< std::endl;
            // std::cout << "rechits with pv, without" <<fat_jet_pv_rechits_T[i_fat_jet][n_matched_rechits] << fat_jet_rechits_T[i_fat_jet][n_matched_rechits] << std::endl;
            // std::cout << "rechit energy and time"<<recHit->energy()<< recHit->time()<< std::endl;


          }
          if (recHit->energy() > 1.5)
          {
            fat_jet_rechit_E_Ecut1p5[i_fat_jet] += recHit->energy();
            fat_jet_rechit_T_Ecut1p5[i_fat_jet] += recHit->time()*recHit->energy();
            fat_jet_pv_rechit_T_Ecut1p5[i_fat_jet] += fat_jet_pv_rechits_T[i_fat_jet][n_matched_rechits] *recHit->energy();

          }
          if (recHit->energy() > 2.0)
          {
            fat_jet_rechit_E_Ecut2[i_fat_jet] += recHit->energy();
            fat_jet_rechit_T_Ecut2[i_fat_jet] += recHit->time()*recHit->energy();
            fat_jet_pv_rechit_T_Ecut2[i_fat_jet] += fat_jet_pv_rechits_T[i_fat_jet][n_matched_rechits] *recHit->energy();

          }
          if (recHit->energy() > 3.0)
          {
            fat_jet_rechit_E_Ecut3[i_fat_jet] += recHit->energy();
            fat_jet_rechit_T_Ecut3[i_fat_jet] += recHit->time()*recHit->energy();
            fat_jet_pv_rechit_T_Ecut3[i_fat_jet] += fat_jet_pv_rechits_T[i_fat_jet][n_matched_rechits] *recHit->energy();

          }

          if (recHit->energy() > 4.0)
          {
            fat_jet_rechit_E_Ecut4[i_fat_jet] += recHit->energy();
            fat_jet_rechit_T_Ecut4[i_fat_jet] += recHit->time()*recHit->energy();
            fat_jet_pv_rechit_T_Ecut4[i_fat_jet] += fat_jet_pv_rechits_T[i_fat_jet][n_matched_rechits] *recHit->energy();

          }
          n_matched_rechits++;
        }
      }
    }
    //cout << "Last Nphoton: " << fJetNPhotons << "\n";
    //std::cout << "n: " << n_matched_rechits << std::endl;
    fat_jet_n_rechits[i_fat_jet] = n_matched_rechits;
    fat_jet_rechit_T[i_fat_jet] = fat_jet_rechit_T[i_fat_jet]/fat_jet_rechit_E[i_fat_jet];
    fat_jet_rechit_T_Ecut4[i_fat_jet] = fat_jet_rechit_T_Ecut4[i_fat_jet]/fat_jet_rechit_E_Ecut4[i_fat_jet];
    fat_jet_rechit_T_Ecut3[i_fat_jet] = fat_jet_rechit_T_Ecut3[i_fat_jet]/fat_jet_rechit_E_Ecut3[i_fat_jet];
    fat_jet_rechit_T_Ecut2[i_fat_jet] = fat_jet_rechit_T_Ecut2[i_fat_jet]/fat_jet_rechit_E_Ecut2[i_fat_jet];
    fat_jet_rechit_T_Ecut1p5[i_fat_jet] = fat_jet_rechit_T_Ecut1p5[i_fat_jet]/fat_jet_rechit_E_Ecut1p5[i_fat_jet];
    fat_jet_rechit_T_Ecut1[i_fat_jet] =  fat_jet_rechit_T_Ecut1[i_fat_jet]/fat_jet_rechit_E_Ecut1[i_fat_jet];
    fat_jet_rechit_T_Ecut0p5[i_fat_jet] = fat_jet_rechit_T_Ecut0p5[i_fat_jet]/fat_jet_rechit_E_Ecut0p5[i_fat_jet]; //incrementing fat_jet counter
    fat_jet_pv_rechit_T[i_fat_jet] = fat_jet_pv_rechit_T[i_fat_jet]/fat_jet_rechit_E[i_fat_jet];
    fat_jet_pv_rechit_T_Ecut4[i_fat_jet] = fat_jet_pv_rechit_T_Ecut4[i_fat_jet]/fat_jet_rechit_E_Ecut4[i_fat_jet];
    fat_jet_pv_rechit_T_Ecut3[i_fat_jet] = fat_jet_pv_rechit_T_Ecut3[i_fat_jet]/fat_jet_rechit_E_Ecut3[i_fat_jet];
    fat_jet_pv_rechit_T_Ecut2[i_fat_jet] =  fat_jet_pv_rechit_T_Ecut2[i_fat_jet]/fat_jet_rechit_E_Ecut2[i_fat_jet];
    fat_jet_pv_rechit_T_Ecut1p5[i_fat_jet] = fat_jet_pv_rechit_T_Ecut1p5[i_fat_jet]/fat_jet_rechit_E_Ecut1p5[i_fat_jet];
    fat_jet_pv_rechit_T_Ecut1[i_fat_jet] = fat_jet_pv_rechit_T_Ecut1[i_fat_jet]/fat_jet_rechit_E_Ecut1[i_fat_jet];
    fat_jet_pv_rechit_T_Ecut0p5[i_fat_jet] = fat_jet_pv_rechit_T_Ecut0p5[i_fat_jet]/fat_jet_rechit_E_Ecut0p5[i_fat_jet]; //incrementing fat_jet counter
    n_fat_Jets++;
    i_fat_jet++;

  } //loop over jets
  return true;
};
*/

bool llp_ntupler::passJetID( const reco::PFJet *jet, int cutLevel) {
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
}

double llp_ntupler::deltaPhi(double phi1, double phi2)
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

double llp_ntupler::deltaR(double eta1, double phi1, double eta2, double phi2)
{
double dphi = deltaPhi(phi1,phi2);
double deta = eta1 - eta2;
return sqrt( dphi*dphi + deta*deta);
};

bool llp_ntupler::fillMC()
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
  /*
  for ( int i_genJet = 0; i_genJet < nGenJets; i_genJet++ )
  {

    unsigned int match_jet_index = 666;
    double min_delta_r = 666.;

    for (int i_jet = 0; i_jet < nJets; i_jet++)
    {

      double current_delta_r = deltaR(genJetEta[i_genJet],genJetPhi[i_genJet] , jetEta[i_jet], jetPhi[i_jet]);

      if ( current_delta_r < min_delta_r )
      {
        min_delta_r = current_delta_r;
        match_jet_index = i_jet;
      }
     }//end matching to jets
     if ( min_delta_r < 0.3 )
     {
       genJet_match_jet_index[i_genJet] = match_jet_index;
       genJet_min_delta_r_match_jet[i_genJet] = min_delta_r;
     }
   }


*/


    /*
    if (isFastsim_) {

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
    } else {

      if (lheInfo.isValid() && lheInfo->weights().size()>0) {

	double nomlheweight = lheInfo->weights()[0].wgt;

	//fill scale variation weights
	if (lheInfo->weights().size()>=9) {
	  for (unsigned int iwgt=0; iwgt<9; ++iwgt) {
	    double wgtval = lheInfo->weights()[iwgt].wgt*genWeight/nomlheweight;
	    scaleWeights->push_back(wgtval);
	  }
	}

	//fill pdf variation weights
	if (firstPdfWeight>=0 && lastPdfWeight>=0 && lastPdfWeight<int(lheInfo->weights().size()) && (lastPdfWeight-firstPdfWeight+1)==100) {

	  //fill pdf variation weights after converting with mc2hessian transformation
	  std::array<double, 100> inpdfweights;
	  for (int iwgt=firstPdfWeight, ipdf=0; iwgt<=lastPdfWeight; ++iwgt, ++ipdf) {
	    inpdfweights[ipdf] = lheInfo->weights()[iwgt].wgt/nomlheweight;
	  }

	  std::array<double, 60> outpdfweights;
	  pdfweightshelper.DoMC2Hessian(inpdfweights.data(),outpdfweights.data());

	  for (unsigned int iwgt=0; iwgt<60; ++iwgt) {
	    double wgtval = outpdfweights[iwgt]*genWeight;
	    pdfWeights->push_back(wgtval);
	  }

	  //fill alpha_s variation weights
	  if (firstAlphasWeight>=0 && lastAlphasWeight>=0 && lastAlphasWeight<int(lheInfo->weights().size())) {
	    for (int iwgt = firstAlphasWeight; iwgt<=lastAlphasWeight; ++iwgt) {
	      double wgtval = lheInfo->weights()[iwgt].wgt*genWeight/nomlheweight;
	      alphasWeights->push_back(wgtval);
	    }
	  }
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

bool llp_ntupler::fillGenParticles()
{
  return true;
};
/*
bool llp_ntupler::fillGenParticles(){
  std::vector<const reco::Candidate*> prunedV;//Allows easier comparison for mother finding
  //Fills selected gen particles
  //double pt_cut = isFourJet ? 20.:20.;//this needs to be done downstream
  const double pt_cut = 0.0;
  int llp_id = isFourJet ? 35:9000006;

  for(size_t i=0; i<genParticles->size();i++)
  {
    if( (abs((*genParticles)[i].pdgId()) >= 1 && abs((*genParticles)[i].pdgId()) <= 6 && ( (*genParticles)[i].status() < 30 ))
       || (abs((*genParticles)[i].pdgId()) >= 11 && abs((*genParticles)[i].pdgId()) <= 16)
       || (abs((*genParticles)[i].pdgId()) == 21 && (*genParticles)[i].status() < 30)
       || (abs((*genParticles)[i].pdgId()) >= 22 && abs((*genParticles)[i].pdgId()) <= 25 && ( (*genParticles)[i].status() < 30))
       || (abs((*genParticles)[i].pdgId()) >= 32 && abs((*genParticles)[i].pdgId()) <= 42)
       || (abs((*genParticles)[i].pdgId()) >= 1000001 && abs((*genParticles)[i].pdgId()) <= 1000039)
       || (abs((*genParticles)[i].pdgId()) == 9000006 || abs((*genParticles)[i].pdgId()) == 9000007))
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
    gParticleMotherId[i] = 0;
    gParticleMotherIndex[i] = -1;

    //For Neutralinos we try to find the decay vertex locaton.
    //Algorithm: Find the first daughter particle that is not a neutralino,
    //and call that the daughter. get the creation vertex of that daughter.
    if ( (gParticleId[i] == 1000022 && gParticleStatus[i] == 22) )
    {
      const reco::Candidate *dau = 0;
      bool foundDaughter = false;
      bool noDaughter = false;
      const reco::Candidate *tmpParticle = prunedV[i];

      while (!foundDaughter && !noDaughter)
      {
        if (tmpParticle->numberOfDaughters() > 0)
        {
          dau = tmpParticle->daughter(0);
          if (dau && dau->pdgId() != 1000022){
            foundDaughter = true;
          }
          else{
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
        gParticleDecayVertexX[i] = dau->vx();
        gParticleDecayVertexY[i] = dau->vy();
        gParticleDecayVertexZ[i] = dau->vz();
      }
    }


    if(prunedV[i]->numberOfMothers() > 0)
    {
      //find the ID of the first mother that has a different ID than the particle itself
      const reco::Candidate* firstMotherWithDifferentID = findFirstMotherWithDifferentID(prunedV[i]);
      if (firstMotherWithDifferentID)
      {
        gParticleMotherId[i] = firstMotherWithDifferentID->pdgId();
        gParticleDecayVertexX[i] = firstMotherWithDifferentID->vx();
        gParticleDecayVertexY[i] = firstMotherWithDifferentID->vy();
        gParticleDecayVertexZ[i] = firstMotherWithDifferentID->vz();
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

    //***************************************
    //Find LLPs production and decay vertices
    //***************************************
    if ( (gParticleId[i] == llp_id || gParticleId[i] == llp_id+1) && gParticleStatus[i] == 22 )
    {
      if (gParticleId[i] == llp_id)
      {
        gLLP_prod_vertex_x[0] = prunedV[i]->vx();
        gLLP_prod_vertex_y[0] = prunedV[i]->vy();
        gLLP_prod_vertex_z[0] = prunedV[i]->vz();
      }
      else if (gParticleId[i] == llp_id+1)
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
          gLLP_beta[0] = sqrt(gParticlePx[i]*gParticlePx[i]+gParticlePy[i]*gParticlePy[i]+gParticlePz[i]*gParticlePz[i])/gParticleE[i];
          gLLP_travel_time[0] = sqrt(pow(gLLP_decay_vertex_x[0]-gLLP_prod_vertex_x[0],2)
                                  +pow(gLLP_decay_vertex_y[0]-gLLP_prod_vertex_y[0],2)
                                  +pow(gLLP_decay_vertex_z[0]-gLLP_prod_vertex_z[0],2))/(30. * gLLP_beta[0]);//1/30 is to convert cm to ns
          double radius = sqrt( pow(gLLP_decay_vertex_x[0],2) + pow(gLLP_decay_vertex_y[0],2) );
          double ecal_radius = 129.0;
          double hcal_radius = 179.0;


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

            double gLLP_daughter_travel_time_hcal= (1./30.)*(hcal_radius-radius)/(tmp.Pt()/tmp.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns

            gLLP_daughter_travel_time[id] = (1./30.)*(ecal_radius-radius)/(tmp.Pt()/tmp.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns
          //Calculate dt from generation point to ECAL face
            double x_ecal = gLLP_decay_vertex_x[0] + 30. * (tmp.Px()/tmp.E())*gLLP_daughter_travel_time[id];
            double y_ecal = gLLP_decay_vertex_y[0] + 30. * (tmp.Py()/tmp.E())*gLLP_daughter_travel_time[id];
            double z_ecal = gLLP_decay_vertex_z[0] + 30. * (tmp.Pz()/tmp.E())*gLLP_daughter_travel_time[id];
            double x_hcal = gLLP_decay_vertex_x[0] + 30. * (tmp.Px()/tmp.E())*gLLP_daughter_travel_time_hcal;
            double y_hcal = gLLP_decay_vertex_y[0] + 30. * (tmp.Py()/tmp.E())*gLLP_daughter_travel_time_hcal;
            double z_hcal = gLLP_decay_vertex_z[0] + 30. * (tmp.Pz()/tmp.E())*gLLP_daughter_travel_time_hcal;

      //if( fabs(z_ecal) < 10 && radius <= 1)
            if( fabs(z_ecal) < 271.6561246934 && radius <= ecal_radius)
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
    	      double min_delta_r_nocorr = 666.;
            double min_delta_r_hcal = 666.;
    	      unsigned int match_jet_index = 666;
            unsigned int match_jet_index_hcal = 666;

    	      double genJet_min_delta_r = 666.;
    	      unsigned int match_genJet_index = 666;

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

        // Correction of eta and phi based on hcal points
            phi = atan((y_hcal-genVertexY)/(x_hcal-genVertexX));
            if  (x_hcal < 0.0){
              phi = TMath::Pi()  + phi;
            }
            phi = deltaPhi(phi,0.0);
            theta = atan(sqrt(pow(x_hcal-genVertexX,2)+pow(y_hcal-genVertexY,2))/abs(z_hcal-genVertexZ));
            eta = -1.0*TMath::Sign(1.0, z_hcal-genVertexZ)*log(tan(theta/2));
            gLLP_daughter_eta_hcalcorr[id] = eta;
            gLLP_daughter_phi_hcalcorr[id] = phi;


    	      for ( int i_jet = 0; i_jet < nGenJets; i_jet++)
    	      {
    		       double genJet_current_delta_r = deltaR(gLLP_daughter_eta[id], gLLP_daughter_phi[id],  genJetEta[i_jet], genJetPhi[i_jet]);
    	        //std::cout << i_jet << " current dR = " << genJet_current_delta_r << eta<<phi<<theta<<tan(theta/2.0)<<log(tan(theta/2.0))<<std::endl;
            	if ( genJet_current_delta_r < genJet_min_delta_r )
            	{
                genJet_min_delta_r = genJet_current_delta_r;
            		match_genJet_index = i_jet;
            		  //std::cout << i_jet << " min dR = " << genJet_min_delta_r << std::endl;
            	}
    	      }//end matching to genJets
    	      for ( int i_jet = 0; i_jet < nJets; i_jet++ )
    	      {
  		        double current_delta_r = deltaR(gLLP_daughter_eta_ecalcorr[id], gLLP_daughter_phi_ecalcorr[id], jetEta[i_jet], jetPhi[i_jet]);
  	          if ( current_delta_r < min_delta_r )
    	        {
    	  	      min_delta_r_nocorr = deltaR(gLLP_daughter_eta[id], gLLP_daughter_phi[id], jetEta[i_jet], jetPhi[i_jet]);
    		        min_delta_r = current_delta_r;
    		        match_jet_index = i_jet;
    		  //std::cout << i_jet << " min dR = " << min_delta_r << std::endl;
    	        }
    	      }//end matching to jets using ECAL radius
            for ( int i_jet = 0; i_jet < nJets; i_jet++ )
    	      {
  		        double current_delta_r = deltaR(gLLP_daughter_eta_hcalcorr[id], gLLP_daughter_phi_hcalcorr[id], jetEta[i_jet], jetPhi[i_jet]);
  	          if ( current_delta_r < min_delta_r_hcal )
    	        {
    		        min_delta_r_hcal = current_delta_r;
    		        match_jet_index_hcal = i_jet;
    	        }
    	      }//end matching to jets using HCAL radius
    	      if ( min_delta_r < 0.3 )
    	      {
    	        gLLP_daughter_match_jet_index[id] = match_jet_index;
    	        gLLP_min_delta_r_match_jet[id] = min_delta_r;
    	        gLLP_min_delta_r_nocorr_match_jet[id] = min_delta_r_nocorr;
    	      }
            if ( min_delta_r < 0.45 )
            {
              gLLP_daughter_match_jet_index_loose[id] = match_jet_index;
              gLLP_min_delta_r_match_jet_loose[id] = min_delta_r;
            }
            if ( min_delta_r_hcal < 0.3 )
    	      {
    	        gLLP_daughter_match_jet_index_hcal[id] = match_jet_index_hcal;
    	        gLLP_min_delta_r_match_jet_hcal[id] = min_delta_r_hcal;
    	      }
            if ( min_delta_r_hcal < 0.45 )
    	      {
    	        gLLP_daughter_match_jet_index_hcal_loose[id] = match_jet_index_hcal;
    	        gLLP_min_delta_r_match_jet_hcal_loose[id] = min_delta_r_hcal;
    	      }
    	      if ( genJet_min_delta_r < 0.3 )
    	      {
    	        gLLP_daughter_match_genJet_index[id] = match_genJet_index;
    	        gLLP_min_delta_r_match_genJet[id] = genJet_min_delta_r;
    	        //std::cout << "min dR = " << min_delta_r << " matched to jet index " << match_jet_index << std::endl;
    	      }

          }
        }
    	  else if (gParticleId[i] == llp_id+1)
    	  {
    	    gLLP_decay_vertex_x[1] = dau->vx();
    	    gLLP_decay_vertex_y[1] = dau->vy();
    	    gLLP_decay_vertex_z[1] = dau->vz();
    	    gLLP_beta[1] = sqrt(gParticlePx[i]*gParticlePx[i]+gParticlePy[i]*gParticlePy[i]+gParticlePz[i]*gParticlePz[i])/gParticleE[i];
    	    gLLP_travel_time[1] = sqrt(pow(gLLP_decay_vertex_x[1]-gLLP_prod_vertex_x[1],2)
    				      +pow(gLLP_decay_vertex_y[1]-gLLP_prod_vertex_y[1],2)
    				      +pow(gLLP_decay_vertex_z[1]-gLLP_prod_vertex_z[1],2))/(30. * gLLP_beta[1]);//1/30 is to convert cm to ns
    	    double radius = sqrt( pow(gLLP_decay_vertex_x[1],2) + pow(gLLP_decay_vertex_y[1],2) );
    	    double ecal_radius = 129.0;
          double hcal_radius = 179.0;
    	    //***************
    	    //Second two LLP daughters belong to LLP->pdgID()=36
          //************
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
            double gLLP_daughter_travel_time_hcal = (1./30.)*(hcal_radius-radius)/(tmp.Pt()/tmp.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns
            gLLP_daughter_travel_time[id+2] = (1./30.)*(ecal_radius-radius)/(tmp.Pt()/tmp.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns

    	      //Calculate dt from generation point to ECAL face
    	      double x_ecal = gLLP_decay_vertex_x[1] + 30. * (tmp.Px()/tmp.E())*gLLP_daughter_travel_time[id+2];
    	      double y_ecal = gLLP_decay_vertex_y[1] + 30. * (tmp.Py()/tmp.E())*gLLP_daughter_travel_time[id+2];
    	      double z_ecal = gLLP_decay_vertex_z[1] + 30. * (tmp.Pz()/tmp.E())*gLLP_daughter_travel_time[id+2];

            double x_hcal = gLLP_decay_vertex_x[1] + 30. * (tmp.Px()/tmp.E())*gLLP_daughter_travel_time_hcal;
    	      double y_hcal = gLLP_decay_vertex_y[1] + 30. * (tmp.Py()/tmp.E())*gLLP_daughter_travel_time_hcal;
    	      double z_hcal = gLLP_decay_vertex_z[1] + 30. * (tmp.Pz()/tmp.E())*gLLP_daughter_travel_time_hcal;
    	      if( fabs(z_ecal) < 271.6561246934 && radius <= ecal_radius)
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
    	      double genJet_min_delta_r = 666.;
            unsigned int match_genJet_index = 666;
    	      double min_delta_r = 666.;
            double min_delta_r_hcal = 666.;
    	      double min_delta_r_nocorr = 666.;
    	      unsigned int match_jet_index = 666;
            unsigned int match_jet_index_hcal = 666;

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

            //Corrections for angles based on HCAL radius
            phi = atan((y_hcal-genVertexY)/(x_hcal-genVertexX));
            if  (x_hcal < 0.0){
              phi = TMath::Pi() + phi;
            }
            phi = deltaPhi(phi,0.0);
            theta = atan(sqrt(pow(x_hcal-genVertexX,2)+pow(y_hcal-genVertexY,2))/abs(z_hcal-genVertexZ));
            eta = -1.0*TMath::Sign(1.0,z_ecal-genVertexZ)*log(tan(theta/2));
            gLLP_daughter_eta_hcalcorr[id+2] = eta;
            gLLP_daughter_phi_hcalcorr[id+2] = phi;
    	      for ( int i_jet = 0; i_jet < nGenJets; i_jet++ )
    	      {
      		    double genJet_current_delta_r = deltaR(gLLP_daughter_eta[id+2], gLLP_daughter_phi[id+2],genJetEta[i_jet], genJetPhi[i_jet]);
      	      if ( genJet_current_delta_r < genJet_min_delta_r )
      	      {
    	  	      genJet_min_delta_r = genJet_current_delta_r;
    		        match_genJet_index = i_jet;
      		  //std::cout << i_jet << " min dR = " << min_delta_r << std::endl;
      	      }
    	      }//end matching to genJets
    	      for ( int i_jet = 0; i_jet < nJets; i_jet++ )
    	      {
              double current_delta_r = deltaR(gLLP_daughter_eta_ecalcorr[id+2], gLLP_daughter_phi_ecalcorr[id+2] , jetEta[i_jet], jetPhi[i_jet]);
          		if ( current_delta_r < min_delta_r )
          		{
          		  min_delta_r_nocorr = deltaR(gLLP_daughter_eta[id+2], gLLP_daughter_phi[id+2], jetEta[i_jet], jetPhi[i_jet]);
          		  min_delta_r = current_delta_r;
          		  match_jet_index = i_jet;
          		}
    	      }//end matching to jets ecal
            for ( int i_jet = 0; i_jet < nJets; i_jet++ )
    	      {
              double current_delta_r = deltaR(gLLP_daughter_eta_hcalcorr[id+2], gLLP_daughter_phi_hcalcorr[id+2], jetEta[i_jet], jetPhi[i_jet]);
          		if ( current_delta_r < min_delta_r_hcal )
          		{
          		  min_delta_r_hcal = current_delta_r;
          		  match_jet_index_hcal = i_jet;
          		}
    	      }//end matching to jets hcal
            if ( min_delta_r < 0.3 )
    	      {
    	        gLLP_daughter_match_jet_index[id+2] = match_jet_index;
    	        gLLP_min_delta_r_match_jet[id+2] = min_delta_r;
    	        gLLP_min_delta_r_nocorr_match_jet[id] = min_delta_r_nocorr;
    	      }
            if ( min_delta_r < 0.45 )
            {
              gLLP_daughter_match_jet_index_loose[id+2] = match_jet_index;
              gLLP_min_delta_r_match_jet_loose[id+2] = min_delta_r;
            }
            if ( min_delta_r_hcal < 0.3 )
    	      {
    	        gLLP_daughter_match_jet_index_hcal[id+2] = match_jet_index_hcal;
    	        gLLP_min_delta_r_match_jet_hcal[id+2] = min_delta_r_hcal;
    	      }
            if ( min_delta_r_hcal < 0.3 )
    	      {
    	        gLLP_daughter_match_jet_index_hcal_loose[id+2] = match_jet_index_hcal;
    	        gLLP_min_delta_r_match_jet_hcal_loose[id+2] = min_delta_r_hcal;
    	      }
    	      if ( genJet_min_delta_r < 0.3 )
    	      {
    	        gLLP_daughter_match_genJet_index[id+2] = match_genJet_index;
    	        gLLP_min_delta_r_match_genJet[id+2] = genJet_min_delta_r;
    	        //std::cout << "min dR = " << min_delta_r << " matched to jet index " << match_jet_index << std::endl;
    	      }
    	    }//for daughters loop
    	  }//if particle ID = 36
    	}//if found daughters
    }
    //******************************
    //QCD Matching
    //******************************
    if (isQCD_) {
      if (abs(gParticleId[i])  <= 6 || abs(gParticleId[i]) == 21)
      {
	      if (gParticleStatus[i] == 23)
	      {
          const reco::Candidate *tmpParticle = prunedV[i];
          TLorentzVector tmp;
          tmp.SetPxPyPzE(tmpParticle->px(), tmpParticle->py(), tmpParticle->pz(), tmpParticle->energy());
          genQCD_pt[nGenQCDParticles] = tmp.Pt();
          genQCD_eta[nGenQCDParticles] = tmp.Eta();
          genQCD_phi[nGenQCDParticles] = tmp.Phi();
          genQCD_e[nGenQCDParticles]  = tmp.E();

          double min_delta_r = 666.;
          unsigned int match_jet_index = 666;

          for ( int i_jet = 0; i_jet < nJets; i_jet++ )
          {
            double current_delta_r = deltaR(genQCD_eta[nGenQCDParticles], genQCD_phi[nGenQCDParticles], jetEta[i_jet], jetPhi[i_jet]);
            if ( current_delta_r < min_delta_r )
            {
              min_delta_r = current_delta_r;
              match_jet_index = i_jet;
            }
          }//end matching to jets hcal
          if ( min_delta_r < 0.3 )
          {
            genParticleQCD_match_jet_index[nGenQCDParticles] = match_jet_index;
            genParticleQCD_min_delta_r_match_jet[nGenQCDParticles] = min_delta_r;
          }
          nGenQCDParticles ++;
        }
      }





    }
  }// for loop of genParticles
  return true;
};

*/

bool llp_ntupler::fillTrigger(const edm::Event& iEvent)
{

  //fill trigger information
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  // std::cout << "\n === TRIGGER PATHS === " << std::endl;
  //********************************************************************
  //Option to save all HLT path names in the ntuple per event
  //Expensive option in terms of ntuple size
  //********************************************************************
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
  //********************************************************************
  // Save trigger decisions in array of booleans
  //********************************************************************

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

  //********************************************************************
  // Print Trigger Objects
  //********************************************************************
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
DEFINE_FWK_MODULE(llp_ntupler);
