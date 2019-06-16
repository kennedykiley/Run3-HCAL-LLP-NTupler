import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
#------ Setup ------#

#initialize the process
process = cms.Process("LLPNtupler")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("Configuration.EventContent.EventContent_cff")

#load input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'file:/mnt/hadoop/store/user/christiw/RunIISummer16_withISR/ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_MC_prod/ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh125_mx50_pl1000_ev100000/crab_CMSSW_8_0_21_ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_ISR_mh125_mx50_pl1000_ev100000_AOD_CaltechT2_v1/190417_225744/0000/ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_step2_1.root',
    'file:/mnt/hadoop/store/user/christiw/RunIISummer16_withISR/ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_MC_prod/ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh125_mx50_pl1000_ev100000/crab_CMSSW_8_0_21_ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_ISR_mh125_mx50_pl1000_ev100000_AOD_CaltechT2_v1/190417_225744/0000/ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_step2_2.root',
    'file:/mnt/hadoop/store/user/christiw/RunIISummer16_withISR/ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_MC_prod/ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh125_mx50_pl1000_ev100000/crab_CMSSW_8_0_21_ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_ISR_mh125_mx50_pl1000_ev100000_AOD_CaltechT2_v1/190417_225744/0000/ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_step2_3.root',
    'file:/mnt/hadoop/store/user/christiw/RunIISummer16_withISR/ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_MC_prod/ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh125_mx50_pl1000_ev100000/crab_CMSSW_8_0_21_ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_ISR_mh125_mx50_pl1000_ev100000_AOD_CaltechT2_v1/190417_225744/0000/ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_step2_4.root',
    'file:/mnt/hadoop/store/user/christiw/RunIISummer16_withISR/ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_MC_prod/ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh125_mx50_pl1000_ev100000/crab_CMSSW_8_0_21_ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_ISR_mh125_mx50_pl1000_ev100000_AOD_CaltechT2_v1/190417_225744/0000/ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_step2_5.root',
    'file:/mnt/hadoop/store/user/christiw/RunIISummer16_withISR/ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_MC_prod/ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh125_mx50_pl1000_ev100000/crab_CMSSW_8_0_21_ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_ISR_mh125_mx50_pl1000_ev100000_AOD_CaltechT2_v1/190417_225744/0000/ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_step2_6.root',
    'file:/mnt/hadoop/store/user/christiw/RunIISummer16_withISR/ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_MC_prod/ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh125_mx50_pl1000_ev100000/crab_CMSSW_8_0_21_ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_ISR_mh125_mx50_pl1000_ev100000_AOD_CaltechT2_v1/190417_225744/0000/ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_step2_7.root',
    'file:/mnt/hadoop/store/user/christiw/RunIISummer16_withISR/ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_MC_prod/ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh125_mx50_pl1000_ev100000/crab_CMSSW_8_0_21_ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_ISR_mh125_mx50_pl1000_ev100000_AOD_CaltechT2_v1/190417_225744/0000/ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_step2_8.root',
    'file:/mnt/hadoop/store/user/christiw/RunIISummer16_withISR/ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_MC_prod/ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh125_mx50_pl1000_ev100000/crab_CMSSW_8_0_21_ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_ISR_mh125_mx50_pl1000_ev100000_AOD_CaltechT2_v1/190417_225744/0000/ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_step2_9.root',
    'file:/mnt/hadoop/store/user/christiw/RunIISummer16_withISR/ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_MC_prod/ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_withISR_mh125_mx50_pl1000_ev100000/crab_CMSSW_8_0_21_ppTohToSS1SS2_SS1Tobb_SS2Tobb_vh_ISR_mh125_mx50_pl1000_ev100000_AOD_CaltechT2_v1/190417_225744/0000/ppTohToSS1SS2_SS1Tobb_SS2Tobb_ggh_withISR_step2_10.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

#TFileService for output
process.TFileService = cms.Service("TFileService",
	fileName = cms.string('llp_ntupler.root'),
    closeFileFast = cms.untracked.bool(True)
)

#load run conditions
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

#------ Declare the correct global tag ------#


process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_v3'

#------ If we add any inputs beyond standard miniAOD event content, import them here ------#

process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')

#------ Analyzer ------#

#list input collections
process.ntuples = cms.EDAnalyzer('llp_ntupler',
    isData = cms.bool(False),
    useGen = cms.bool(True),
    isFastsim = cms.bool(False),
    enableTriggerInfo = cms.bool(True),
    enableEcalRechits = cms.bool(False),
    enableCaloJet = cms.bool(True),
    readGenVertexTime = cms.bool(True),#need to be false for displaced samples
    genParticles_t0 = cms.InputTag("genParticles", "t0", ""),
    triggerPathNamesFile = cms.string("cms_lpc_llp/llp_ntupler/data/trigger_names_llp_v1.dat"),
    eleHLTFilterNamesFile = cms.string("SUSYBSMAnalysis/RazorTuplizer/data/RazorElectronHLTFilterNames.dat"),
    muonHLTFilterNamesFile = cms.string("SUSYBSMAnalysis/RazorTuplizer/data/RazorMuonHLTFilterNames.dat"),
    photonHLTFilterNamesFile = cms.string("SUSYBSMAnalysis/RazorTuplizer/data/RazorPhotonHLTFilterNames.dat"),

    #vertices = cms.InputTag("offlinePrimaryVerticesWithBS"),  # for non-timing case
    vertices = cms.InputTag("offlinePrimaryVertices", "", "RECO"),
    muons = cms.InputTag("muons"),
    electrons = cms.InputTag("gedGsfElectrons"),
    taus = cms.InputTag("hpsPFTauProducer"),
    photons = cms.InputTag("gedPhotons"),
    jetsCalo = cms.InputTag("ak4CaloJets","","RECO"),
    jetsPF = cms.InputTag("ak4PFJets"),
    jets = cms.InputTag("ak4PFJetsCHS"),
    jetsPuppi = cms.InputTag("ak4PFJets"),
    jetsAK8 = cms.InputTag("ak8PFJetsCHS"),
    mets = cms.InputTag("pfMet"),
    #metsNoHF = cms.InputTag("pfMet30"),
    metsPuppi = cms.InputTag("pfMet"),
    pfCands = cms.InputTag("particleFlow","","RECO"),

    #packedPfCands = cms.InputTag("packedPFCandidates"),

    genParticles = cms.InputTag("genParticles"),

    #packedGenParticles = cms.InputTag("packedGenParticles"),
    #prunedGenParticles = cms.InputTag("prunedGenParticles"),
    genMetsCalo = cms.InputTag("genMetCalo"),
    genMetsTrue = cms.InputTag("genMetTrue"),
    genJets = cms.InputTag("ak4GenJets"),

    triggerBits = cms.InputTag("TriggerResults","","HLT"),
    #triggerBits = cms.InputTag("TriggerResults","","RECO"),
    hepMC = cms.InputTag("generatorSmeared", "", "SIM"),

    #triggerPrescales = cms.InputTag("patTrigger"),
    #triggerObjects = cms.InputTag("selectedPatTrigger"),

    metFilterBits = cms.InputTag("TriggerResults", "", "RECO"),

    #hbheNoiseFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"),
    #hbheTightNoiseFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun2Tight"),
    #hbheIsoNoiseFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHEIsoNoiseFilterResult"),

    #BadChargedCandidateFilter = cms.InputTag("BadChargedCandidateFilter",""),
    #BadMuonFilter = cms.InputTag("BadPFMuonFilter",""),

    #lheInfo = cms.InputTag("externalLHEProducer", "", ""),
    genInfo = cms.InputTag("generator", "", "SIM"),

    tracks = cms.InputTag("generalTracks", "", "RECO"),
    #trackTime = cms.InputTag("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModel"),
    #trackTimeReso = cms.InputTag("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModelResolution"),

    puInfo = cms.InputTag("addPileupInfo", "", "HLT"), #uncomment if no pre-mixing
    #puInfo = cms.InputTag("mixData", "", "HLT"), #uncomment for samples with pre-mixed pileup
    #hcalNoiseInfo = cms.InputTag("hcalnoise", "", "RECO"),

    #secondaryVertices = cms.InputTag("inclusiveSecondaryVertices", "", "RECO"),
    secondaryVertices = cms.InputTag("inclusiveCandidateSecondaryVertices","", "RECO"),

    rhoAll = cms.InputTag("fixedGridRhoAll", "", "RECO"),

    rhoFastjetAll = cms.InputTag("fixedGridRhoFastjetAll", "", "RECO"),
    rhoFastjetAllCalo = cms.InputTag("fixedGridRhoFastjetAllCalo", "", "RECO"),
    rhoFastjetCentralCalo = cms.InputTag("fixedGridRhoFastjetCentralCalo", "", "RECO"),
    rhoFastjetCentralChargedPileUp = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp", "", "RECO"),
    rhoFastjetCentralNeutral = cms.InputTag("fixedGridRhoFastjetCentralNeutral", "", "RECO"),

    beamSpot = cms.InputTag("offlineBeamSpot", "", "RECO"),
    pfClusters = cms.InputTag("particleFlowClusterECAL","","RECO"),
    ebRecHits = cms.InputTag("reducedEcalRecHitsEB", "","RECO"),
    #ebRecHits = cms.InputTag("EcalRecHit", "reducedEcalRecHitsEB", "RECO"),
    eeRecHits  = cms.InputTag("reducedEcalRecHitsEE", "","RECO"),
    esRecHits = cms.InputTag("reducedEcalRecHitsES", "","RECO"),
    #ebeeClusters = cms.InputTag("reducedEgamma", "reducedEBEEClusters", "RECO"),
    ebeeClusters = cms.InputTag("particleFlowEGamma", "EBEEClusters", "RECO"),
    esClusters = cms.InputTag("particleFlowEGamma", "ESClusters", "RECO"),
    #conversions = cms.InputTag("reducedEgamma", "reducedConversions", "RECO"),
    conversions = cms.InputTag("allConversions", "", "RECO"),

    #singleLegConversions = cms.InputTag("reducedEgamma", "reducedSingleLegConversions", "RECO"),
    singleLegConversions = cms.InputTag("particleFlowEGamma", "", "RECO"),

    gedGsfElectronCores = cms.InputTag("gedGsfElectronCores", "", "RECO"),
    gedPhotonCores = cms.InputTag("gedPhotonCore", "", "RECO"),
    #superClusters = cms.InputTag("reducedEgamma", "reducedSuperClusters", "RECO"),

    #lostTracks = cms.InputTag("lostTracks", "", "RECO")
)

#run
process.p = cms.Path( process.ntuples)
