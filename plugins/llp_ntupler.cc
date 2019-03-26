// -*- C++ -*-
// Class:      llp_ntupler
/*
  Description: Base class for miniAOD analysis with CRAB
*/
//         Author:  Caltech razor team
//         Created:  Thu, 17 Jul 2014 15:00:06 GMT

#include "llp_ntupler.h"
//------ Constructors and destructor ------//
llp_ntupler::llp_ntupler(const edm::ParameterSet& iConfig):
  //get inputs from config file
  isData_(iConfig.getParameter<bool> ("isData")),
  isFourJet_(iConfig.getParameter<bool> ("isFourJet")),
  useGen_(iConfig.getParameter<bool> ("useGen")),
  isFastsim_(iConfig.getParameter<bool> ("isFastsim")),
  isQCD_(iConfig.getParameter<bool> ("isQCD")),
  enableTriggerInfo_(iConfig.getParameter<bool> ("enableTriggerInfo")),
  enableRecHitInfo_(iConfig.getParameter<bool> ("enableRecHitInfo")),
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

//------ Enable the desired set of branches ------//
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
  llpTree->Branch("genMetPt", &genMetPt, "genMetPt/F");
  llpTree->Branch("genMetPhi", &genMetPhi, "genMetPhi/F");
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

  /*scaleWeights = new std::vector<float>; scaleWeights->clear();
  pdfWeights = new std::vector<float>; pdfWeights->clear();
  alphasWeights = new std::vector<float>; alphasWeights->clear();
  if (isFastsim_) {
    llpTree->Branch("lheComments", "std::string",&lheComments);
  }
  llpTree->Branch("scaleWeights", "std::vector<float>",&scaleWeights);
  llpTree->Branch("pdfWeights", "std::vector<float>",&pdfWeights);
  llpTree->Branch("alphasWeights", "std::vector<float>",&alphasWeights);
  */
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
    reset_event_variables();
    reset_photon_variable();
    reset_jet_variables();
    reset_gen_llp_variable();
    reset_gen_jet_variable();
    reset_qcd_variables();
}

void llp_ntupler::reset_event_variables()
{
  eventNum = 0;
  lumiNum = 0;
  runNum = 0;
  pvX = -99.0;
  pvY = -99.0;
  pvZ = -99.0;
  nPV = -1;
  Rho = -99.0;
  nPUmean = -1;
  nPU = -1;
  return;
};

void llp_ntupler::reset_photon_variable()
{
  fJetNPhotons = 0;
  for (int i=0; i< OBJECTARRAYSIZE; i++) {
    fJetPhotonPt[i] = 0.0;
    fJetPhotonEta[i] = 0.0;
    fJetPhotonPhi[i] = 0.0;
    fJetPhotonSeedRecHitE[i]      = -99.0;
    fJetPhotonSeedRecHitEta[i]      = -99.0;
    fJetPhotonSeedRecHitPhi[i]      = -99.0;
    fJetPhotonSeedRecHitTime[i]      = -99.0;
  }
  return;
};

void llp_ntupler::reset_jet_variables()
{
  nJets = 0;
  for ( int i = 0; i < OBJECTARRAYSIZE; i++)
  {
    jetE[i] = 0.0;
    jetPt[i] = 0.0;
    jetEta[i] = 0.0;
    jetPhi[i] = 0.0;
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
    jet_charged_hadron_multiplicity[i] = -99;
    jet_neutral_hadron_multiplicity[i] = -99;
    jet_photon_multiplicity[i] = -99;
    jet_electron_multiplicity[i] = -99;
    jet_muon_multiplicity[i] = -99;
    jet_HF_hadron_multiplicity[i] = -99;
    jet_HF_em_multiplicity[i] = -99;
    jet_charged_multiplicity[i] = -99;
    jet_neutral_multiplicity[i] = -99;
    jetMatchedGenPt[i] = 0.0;
    jetMatchedGenEta[i] = 0.0;
    jetMatchedGenPhi[i] = 0.0;
    jetMatchedGenMass[i] = 0.0;
    jetMatchedGenTime[i] = 0.0;
    jet_n_rechits[i] = 0;
    jet_rechit_E[i] = 0.0;
    jet_rechit_T[i] = 0.0;
    jet_rechit_E_Ecut3[i] = 0.0; //energy with a 2 GeV cut
    jet_rechit_T_Ecut3[i] = 0.0;

    jet_rechit_E_Ecut4[i] = 0.0; //energy with a 2 GeV cut
    jet_rechit_T_Ecut4[i] = 0.0;
    jet_rechit_E_Ecut2[i] = 0.0; //energy with a 2 GeV cut
    jet_rechit_T_Ecut2[i] = 0.0;
    jet_rechit_E_Ecut1p5[i] = 0.0; //energy with a 2 GeV cut
    jet_rechit_T_Ecut1p5[i] = 0.0;
    jet_rechit_E_Ecut1[i] = 0.0; //energy with a 2 GeV cut
    jet_rechit_T_Ecut1[i] = 0.0;
    jet_rechit_E_Ecut0p5[i] = 0.0; //energy with a 2 GeV cut
    jet_rechit_T_Ecut0p5[i] = 0.0;


    jet_pv_rechit_T[i] = 0.0;
    jet_pv_rechit_T_Ecut4[i] = 0.0;
    jet_pv_rechit_T_Ecut2[i] = 0.0;
    jet_pv_rechit_T_Ecut1p5[i] = 0.0;
    jet_pv_rechit_T_Ecut1[i] = 0.0;
    jet_pv_rechit_T_Ecut0p5[i] = 0.0;

    for(int j =0; j < OBJECTARRAYSIZE;j++)
    {
      jet_rechits_E[i][j] = -666.;
      jet_rechits_T[i][j] = -666.;
      jet_pv_rechits_T[i][j] = -666.;
    }
  }
  return;
};

void llp_ntupler::reset_gen_llp_variable()
{
  for ( int i = 0; i < LLP_ARRAY_SIZE; i++ )
  {
    gLLP_prod_vertex_x[i] = -666.;
    gLLP_prod_vertex_y[i] = -666.;
    gLLP_prod_vertex_z[i] = -666.;
    gLLP_decay_vertex_x[i] = -666.;
    gLLP_decay_vertex_y[i] = -666.;
    gLLP_decay_vertex_z[i] = -666.;
    gLLP_beta[i] = -666.;
    gLLP_travel_time[i] = -666.;
  }

  for ( int i = 0; i < LLP_DAUGHTER_ARRAY_SIZE; i++ )
  {
    gLLP_daughter_pt[i] = -666.;
    gLLP_daughter_eta[i] = -666.;
    gLLP_daughter_phi[i] = -666.;
    gLLP_daughter_eta_ecalcorr[i] = -666.;
    gLLP_daughter_phi_ecalcorr[i] = -666.;
    gLLP_daughter_eta_hcalcorr[i] = -666.;
    gLLP_daughter_phi_hcalcorr[i] = -666.;
    gLLP_daughter_e[i] = -666.;
    gLLP_daughter_travel_time[i] = -666.;
    gen_time[i] = -666.;
    gen_time_pv[i] = -666.;
    photon_travel_time[i] = -666.;
    photon_travel_time_pv[i] = -666.;
    gLLP_daughter_match_jet_index[i] = 666;
    gLLP_daughter_match_jet_index_hcal[i] = 666;
    gLLP_daughter_match_jet_index_hcal_loose[i] = 666;
    gLLP_daughter_match_jet_index_loose[i] = 666;
    gLLP_min_delta_r_match_jet[i] = -666.;
    gLLP_min_delta_r_match_jet_hcal[i] = -666.;
    gLLP_min_delta_r_match_jet_loose[i] = -666.;
    gLLP_min_delta_r_match_jet_hcal_loose[i] = -666.;
    gLLP_min_delta_r_nocorr_match_jet[i] = -666.;
    gLLP_daughter_match_genJet_index[i] = 666;
    gLLP_min_delta_r_match_genJet[i] = -666.;

  }
  return;
};

void llp_ntupler::reset_gen_jet_variable()
{
  nGenJets = 0;
  for ( int i = 0; i < OBJECTARRAYSIZE; i++ )
  {
    genJetE[i] = -666.;
    genJetPt[i] = -666.;
    genJetEta[i] = -666.;
    genJetPhi[i] = -666.;
    genJetME[i] = -666.;
    genJet_match_jet_index[i] = 666;
    genJet_min_delta_r_match_jet[i] = -666.;
  }
  return;
};
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
  isFourJet = isFourJet_;
  isQCD = isQCD_;
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
  Rho = *rhoFastjetAll;

  //Fill Pileup info
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
    jet_charged_hadron_multiplicity[i_jet] = j.chargedHadronMultiplicity();
    jet_neutral_hadron_multiplicity[i_jet] = j.neutralHadronMultiplicity();
    jet_photon_multiplicity[i_jet] = j.photonMultiplicity();
    jet_electron_multiplicity[i_jet] = j.electronMultiplicity();
    jet_muon_multiplicity[i_jet] = j.muonMultiplicity();
    jet_HF_hadron_multiplicity[i_jet] = j.HFHadronMultiplicity();
    jet_HF_em_multiplicity[i_jet] = j.HFEMMultiplicity();
    jet_charged_multiplicity[i_jet] = j.chargedMultiplicity();
    jet_neutral_multiplicity[i_jet] = j.neutralMultiplicity();


    //*************************************
    //find photons inside the jet
    //*************************************
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
          jet_rechit_E[i_jet] += recHit->energy();
          jet_rechit_T[i_jet] += recHit->time()*recHit->energy();
          jet_rechits_E[i_jet][n_matched_rechits] = recHit->energy();
    	    jet_rechits_T[i_jet][n_matched_rechits] = recHit->time();
    	    double rechit_x = ecal_radius * cos(recHitPos.phi());
    	    double rechit_y = ecal_radius * sin(recHitPos.phi());
    	    double rechit_z = ecal_radius * sinh(recHitPos.eta());
    	    double photon_pv_travel_time = (1./30) * sqrt(pow(pvX-rechit_x,2)+pow(pvY-rechit_y,2)+pow(pvZ-rechit_z,2));
          jet_pv_rechits_T[i_jet][n_matched_rechits] = recHit->time()+(1./30)*ecal_radius*cosh(recHitPos.eta()) - photon_pv_travel_time;
    	    jet_pv_rechit_T[i_jet] += recHit->energy()*jet_pv_rechits_T[i_jet][n_matched_rechits];
          // std::cout << jet_pv_rechits_T[i_jet][n_matched_rechits] << jet_rechits_T[i_jet][n_matched_rechits] << std::endl;
          if (recHit->energy() > 0.5)
    	    {
        		jet_rechit_E_Ecut0p5[i_jet] += recHit->energy();
        		jet_rechit_T_Ecut0p5[i_jet] += recHit->time()*recHit->energy();
            jet_pv_rechit_T_Ecut0p5[i_jet] += jet_pv_rechits_T[i_jet][n_matched_rechits] *recHit->energy();
    	    }
          if (recHit->energy() > 1.0)
    	    {
    		    jet_rechit_E_Ecut1[i_jet] += recHit->energy();
    		    jet_rechit_T_Ecut1[i_jet] += recHit->time()*recHit->energy();
            jet_pv_rechit_T_Ecut1[i_jet] += jet_pv_rechits_T[i_jet][n_matched_rechits] *recHit->energy();
            // std::cout << "rechit time, with pv"<<jet_rechit_T_Ecut1[i_jet]<< jet_pv_rechit_T_Ecut1[i_jet]<< std::endl;
            // std::cout << "rechits with pv, without" <<jet_pv_rechits_T[i_jet][n_matched_rechits] << jet_rechits_T[i_jet][n_matched_rechits] << std::endl;
            // std::cout << "rechit energy and time"<<recHit->energy()<< recHit->time()<< std::endl;


    	    }
          if (recHit->energy() > 1.5)
    	    {
    		    jet_rechit_E_Ecut1p5[i_jet] += recHit->energy();
    		    jet_rechit_T_Ecut1p5[i_jet] += recHit->time()*recHit->energy();
            jet_pv_rechit_T_Ecut1p5[i_jet] += jet_pv_rechits_T[i_jet][n_matched_rechits] *recHit->energy();

    	    }
          if (recHit->energy() > 2.0)
    	    {
    		    jet_rechit_E_Ecut2[i_jet] += recHit->energy();
    		    jet_rechit_T_Ecut2[i_jet] += recHit->time()*recHit->energy();
            jet_pv_rechit_T_Ecut2[i_jet] += jet_pv_rechits_T[i_jet][n_matched_rechits] *recHit->energy();

    	    }
    	    if (recHit->energy() > 3.0)
          {
            jet_rechit_E_Ecut3[i_jet] += recHit->energy();
            jet_rechit_T_Ecut3[i_jet] += recHit->time()*recHit->energy();
            jet_pv_rechit_T_Ecut3[i_jet] += jet_pv_rechits_T[i_jet][n_matched_rechits] *recHit->energy();

          }

    	    if (recHit->energy() > 4.0)
          {
            jet_rechit_E_Ecut4[i_jet] += recHit->energy();
            jet_rechit_T_Ecut4[i_jet] += recHit->time()*recHit->energy();
            jet_pv_rechit_T_Ecut4[i_jet] += jet_pv_rechits_T[i_jet][n_matched_rechits] *recHit->energy();

          }
    	    n_matched_rechits++;
        }
      }
    }
    //cout << "Last Nphoton: " << fJetNPhotons << "\n";
    //std::cout << "n: " << n_matched_rechits << std::endl;
    jet_n_rechits[i_jet] = n_matched_rechits;
    jet_rechit_T[i_jet] = jet_rechit_T[i_jet]/jet_rechit_E[i_jet];
    jet_rechit_T_Ecut4[i_jet] = jet_rechit_T_Ecut4[i_jet]/jet_rechit_E_Ecut4[i_jet];
    jet_rechit_T_Ecut3[i_jet] = jet_rechit_T_Ecut3[i_jet]/jet_rechit_E_Ecut3[i_jet];
    jet_rechit_T_Ecut2[i_jet] = jet_rechit_T_Ecut2[i_jet]/jet_rechit_E_Ecut2[i_jet];
    jet_rechit_T_Ecut1p5[i_jet] = jet_rechit_T_Ecut1p5[i_jet]/jet_rechit_E_Ecut1p5[i_jet];
    jet_rechit_T_Ecut1[i_jet] =  jet_rechit_T_Ecut1[i_jet]/jet_rechit_E_Ecut1[i_jet];
    jet_rechit_T_Ecut0p5[i_jet] = jet_rechit_T_Ecut0p5[i_jet]/jet_rechit_E_Ecut0p5[i_jet]; //incrementing jet counter
    jet_pv_rechit_T[i_jet] = jet_pv_rechit_T[i_jet]/jet_rechit_E[i_jet];
    jet_pv_rechit_T_Ecut4[i_jet] = jet_pv_rechit_T_Ecut4[i_jet]/jet_rechit_E_Ecut4[i_jet];
    jet_pv_rechit_T_Ecut3[i_jet] = jet_pv_rechit_T_Ecut3[i_jet]/jet_rechit_E_Ecut3[i_jet];
    jet_pv_rechit_T_Ecut2[i_jet] =  jet_pv_rechit_T_Ecut2[i_jet]/jet_rechit_E_Ecut2[i_jet];
    jet_pv_rechit_T_Ecut1p5[i_jet] = jet_pv_rechit_T_Ecut1p5[i_jet]/jet_rechit_E_Ecut1p5[i_jet];
    jet_pv_rechit_T_Ecut1[i_jet] = jet_pv_rechit_T_Ecut1[i_jet]/jet_rechit_E_Ecut1[i_jet];
    jet_pv_rechit_T_Ecut0p5[i_jet] = jet_pv_rechit_T_Ecut0p5[i_jet]/jet_rechit_E_Ecut0p5[i_jet]; //incrementing jet counter
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
    genJetME[nGenJets] = j.invisibleEnergy();
    nGenJets++;
  }


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
    	    /*
    	    Second two LLP daughters belong to LLP->pdgID()=36
          */
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
