import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# ------ Setup ------ #

#initialize the process
from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process("DisplacedHcalJetNTuplizer", Run3) # line added to fix PCastorRcd error
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load("cms_lpc_llp.Run3-HCAL-LLP-NTupler.metFilters_cff_2022")

# Fix GEM Error
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag

# ----- Parse Arguments ----- #

options = VarParsing.VarParsing()

options.register('isData',
    False, # default value
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.bool,
    "is Data"
)

options.register('isSignal',
    True, # default value
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.bool,
    "is Signal"
)

options.register('skipEvents',
    0, # default value
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.int,
    "Number of events to skip"
)

options.register('processEvents',
    -1, # default value
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.int,
    "Number of events to process"
)

options.register('inputFiles',
    "input.root",
    VarParsing.VarParsing.multiplicity.list,
    VarParsing.VarParsing.varType.string,
    "Input files"
)

options.register('outputFile',
    "output.root",
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,
    "Output file"
)

options.register('debug',
    False, # default value
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.bool,
    "debug mode"
)

options.parseArguments()

print(" ")
print("Using options:")
print(f" isData        ={options.isData}")
print(f" isSignal      ={options.isSignal}")
print(f" skipEvents    ={options.skipEvents}")
print(f" processEvents ={options.processEvents}")
print(f" inputFiles    ={options.inputFiles}")
print(f" outputFile    ={options.outputFile}")
print(f" debug         ={options.debug}")
print(" ")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.processEvents) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
)

# Load input files

inputFiles = options.inputFiles
if len(options.inputFiles) == 1 and options.inputFiles[0][-4:] == ".txt": 
    print( "Reading input files from", options.inputFiles, ":" )

    file = open(options.inputFiles[0], 'r')
    lines = file.readlines()

    unpackedFilelist = []
    for line in lines:
        line_temp = line.replace(" ","")
        if line_temp[0] == "#": continue
        unpackedFilelist.append( line.strip() )
        print( "   ->", line.strip() )

    inputFiles = unpackedFilelist

process.source = cms.Source( "PoolSource",
    fileNames  = cms.untracked.vstring( inputFiles ),
    skipEvents = cms.untracked.uint32(options.skipEvents)
)

# TFileService for output
process.TFileService = cms.Service( "TFileService",
    fileName = cms.string(options.outputFile),
    closeFileFast = cms.untracked.bool(True)
)

# ----- Load Run Conditions ----- #

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('Configuration.Geometry.GeometryIdeal_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')


# ------ Declare the correct global tag ------ #

if options.isData: process.GlobalTag.globaltag = '124X_dataRun3_Prompt_v4'
else:              process.GlobalTag = GlobalTag(process.GlobalTag,'auto:run3_mc_FULL','')

# ------ If we add any inputs beyond standard event content, import them here ------ #

process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')

process.output = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    eventAutoFlushCompressedSize = cms.untracked.int32(-900),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string('miniAOD-prod_PAT.root'),
    outputCommands = cms.untracked.vstring('keep *'),
)
  
# ------ Analyzer ------ #

# For AOD Track variables

process.load("RecoTracker.TkNavigation.NavigationSchoolESProducer_cfi")
process.MaterialPropagator = cms.ESProducer('PropagatorWithMaterialESProducer',
    ComponentName = cms.string('PropagatorWithMaterial'),
    Mass = cms.double(0.105),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('alongMomentum'),
    SimpleMagneticField = cms.string(''),
    ptMin = cms.double(-1.0),
    useRungeKutta = cms.bool(False)
)

process.TransientTrackBuilderESProducer = cms.ESProducer('TransientTrackBuilderESProducer',
    ComponentName = cms.string('TransientTrackBuilder')
)


# Main Analyzer

process.DisplacedHcalJets = cms.EDAnalyzer('DisplacedHcalJetNTuplizer',
    debug =  cms.bool( options.debug ),
    isData = cms.bool( options.isData ),
    isSignal = cms.bool( options.isSignal ),
    useGen = cms.bool(False),
    isRECO = cms.bool(True),                                
    isRAW = cms.bool(False),                                
    isBParkAOD = cms.bool(False),
    isFastsim = cms.bool(False),
    readMuonDigis = cms.bool(True),
    enableTriggerInfo = cms.bool(True),
    enableEcalRechits = cms.bool(False),
    enableCaloJet = cms.bool(True),
    enableGenLLPInfo = cms.bool(True),
    readGenVertexTime = cms.bool(False),#need to be false for displaced samples
    genParticles_t0 = cms.InputTag("genParticles", "t0", ""),
    triggerPathNamesFile = cms.string("cms_lpc_llp/Run3-HCAL-LLP-NTupler/data/HLTPathsLLPJetsHCAL.dat"),
    #triggerPathNamesFile = cms.FileInPath("/afs/cern.ch/work/k/kikenned/LLPNTupler/Run3-HCAL-LLP-NTupler/data/HLTPathsLLPJetsHCAL.dat"), #"../data/HLTPathsLLPJetsHCAL.dat"),
    #eleHLTFilterNamesFile = cms.string("SUSYBSMAnalysis/RazorTuplizer/data/RazorElectronHLTFilterNames.dat"),
    #muonHLTFilterNamesFile = cms.string("cms_lpc_llp/llp_ntupler/data/MuonHLTFilterNames.dat"),
    #photonHLTFilterNamesFile = cms.string("SUSYBSMAnalysis/RazorTuplizer/data/RazorPhotonHLTFilterNames.dat"),

    #vertices = cms.InputTag("offlinePrimaryVerticesWithBS"),  # for non-timing case
    vertices = cms.InputTag("offlinePrimaryVertices", "", "RECO"),

    electrons = cms.InputTag("gedGsfElectrons"),
    muons = cms.InputTag("muons"),
    taus = cms.InputTag("selectedPatTaus"),
    photons = cms.InputTag("gedPhotons"),
    pfjetsAK4 = cms.InputTag("selectedPatJets"),
    calojetsAK4 = cms.InputTag("ak4CaloJets","","RECO"),
    pfjetsAK8 = cms.InputTag("selectedPatJetsAK8PFCHS"),
    calojetsAK8 = cms.InputTag("ak8CaloJets","","RECO"),
    #jetsPF = cms.InputTag("ak4PFJets"),
    #jets = cms.InputTag("ak4PFJetsCHS"),
    #jets = cms.InputTag("selectedPatJets"),
    #jets = cms.InputTag("ak4PFJetsPuppi"),
    #jetsPuppi = cms.InputTag("ak4PFJets"),
    #jetsAK8 = cms.InputTag("ak8PFJetsCHS"),
    #jetsAK8 = cms.InputTag("selectedPatJetsAK8PFCHS"),

    #jetsAK8 = cms.InputTag("ak8PFJetsPuppi"),

    #mets = cms.InputTag("slimmedMETs"),
    met = cms.InputTag("patMETs"),
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

    triggerPrescales = cms.InputTag("patTrigger"),
    triggerObjects = cms.InputTag("selectedPatTrigger"),

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

    #puInfo = cms.InputTag("addPileupInfo", "", "HLT"), #uncomment if no pre-mixing
    puInfo = cms.InputTag("mixData", "", "HLT"), #uncomment for samples with pre-mixed pileup
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
    #hbRecHits = cms.InputTag("reducedHcalRecHits", "hbhereco","RECO"),
    hbRecHits = cms.InputTag("hbhereco", "","RECO"),
    #ebRecHits = cms.InputTag("EcalRecHit", "reducedEcalRecHitsEB", "RECO"),
    ebRecHits = cms.InputTag("ecalRecHit", "EcalRecHitsEB", "RECO"), 
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
    generalTracks = cms.InputTag("generalTracks", "", "RECO"),
    #superClusters = cms.InputTag("reducedEgamma", "reducedSuperClusters", "RECO"),

    #lostTracks = cms.InputTag("lostTracks", "", "RECO")

    electron_cutbasedID_decisions_veto = cms.InputTag("egmGsfElectronIDs", "cutBasedElectronID-Fall17-94X-V2-veto", ""),
    electron_cutbasedID_decisions_loose = cms.InputTag("egmGsfElectronIDs", "cutBasedElectronID-Fall17-94X-V2-loose", ""),
    electron_cutbasedID_decisions_medium = cms.InputTag("egmGsfElectronIDs", "cutBasedElectronID-Fall17-94X-V2-medium", ""),
    electron_cutbasedID_decisions_tight = cms.InputTag("egmGsfElectronIDs", "cutBasedElectronID-Fall17-94X-V2-tight", ""),
    electron_mvaIsoID_decisions_wp80 = cms.InputTag("egmGsfElectronIDs", "mvaEleID-Fall17-iso-V2-wp80", ""),
    electron_mvaIsoID_decisions_wp90 = cms.InputTag("egmGsfElectronIDs", "mvaEleID-Fall17-iso-V2-wp90", ""),
    electron_mvaIsoID_decisions_wpHZZ = cms.InputTag("egmGsfElectronIDs", "mvaEleID-Fall17-iso-V2-wpHZZ", ""),
    electron_mvaIsoID_decisions_wpLoose = cms.InputTag("egmGsfElectronIDs", "mvaEleID-Fall17-iso-V2-wpLoose", ""),
    electron_mvaNoIsoID_decisions_wp80 = cms.InputTag("egmGsfElectronIDs", "mvaEleID-Fall17-noIso-V2-wp80", ""),
    electron_mvaNoIsoID_decisions_wp90 = cms.InputTag("egmGsfElectronIDs", "mvaEleID-Fall17-noIso-V2-wp90", ""),
    electron_mvaNoIsoID_decisions_wpLoose = cms.InputTag("egmGsfElectronIDs", "mvaEleID-Fall17-noIso-V2-wpLoose", ""),
    photon_cutbasedID_decisions_loose = cms.InputTag("egmPhotonIDs", "cutBasedPhotonID-Fall17-94X-V2-loose", ""),
    photon_cutbasedID_decisions_medium = cms.InputTag("egmPhotonIDs", "cutBasedPhotonID-Fall17-94X-V2-medium", ""),
    photon_cutbasedID_decisions_tight = cms.InputTag("egmPhotonIDs", "cutBasedPhotonID-Fall17-94X-V2-tight", ""),
    photon_mvaID_decisions_wp80 = cms.InputTag("egmPhotonIDs", "mvaPhoID-RunIIFall17-v2-wp80", ""),
    photon_mvaID_decisions_wp90 = cms.InputTag("egmPhotonIDs", "mvaPhoID-RunIIFall17-v2-wp90", ""),
)

# ----- Add Additional Info ----- #

# Add jettiness for AK8 jets

process.load('RecoJets.JetProducers.nJettinessAdder_cfi')
process.NjettinessAK8CHS = process.Njettiness.clone()

# Add object ids

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
electron_id_config = cms.PSet(electron_ids = cms.vstring([                   
                    'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
                    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
                    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff',
                    'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff',
                    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff', 
                    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff',
                    ]))  
photon_id_config = cms.PSet(photon_ids = cms.vstring([                   
            'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring16_V2p2_cff',
            'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff',
            "RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V2_cff",
            "RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff"            
                    ]))  


                 
switchOnVIDElectronIdProducer(process,DataFormat.AOD)
switchOnVIDPhotonIdProducer(process,DataFormat.AOD) 
process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag("gedGsfElectrons")
process.electronMVAValueMapProducer.src = cms.InputTag("gedGsfElectrons")
process.photonMVAValueMapProducer.src = cms.InputTag("gedPhotons")
#process.electronRegressionValueMapProducer.src = cms.InputTag('reducedEgamma','gedGsfElectrons')
for idmod in electron_id_config.electron_ids.value():
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
for idmod in photon_id_config.photon_ids.value():
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

process.load("CommonTools.RecoAlgos.sortedPFPrimaryVertices_cfi")
process.primaryVertexAssociation = process.sortedPFPrimaryVertices.clone(
    qualityForPrimary = cms.int32(2),
    produceSortedVertices = cms.bool(False),
    producePileUpCollection  = cms.bool(False),  
    produceNoPileUpCollection = cms.bool(False)
    )

#PAT Stuff
process.load('PhysicsTools.PatAlgos.producersLayer1.tauProducer_cff')
process.load('PhysicsTools.PatAlgos.producersLayer1.jetProducer_cff')
process.load('PhysicsTools.PatAlgos.producersLayer1.electronProducer_cff')
process.load('PhysicsTools.PatAlgos.producersLayer1.photonProducer_cff')
process.load('PhysicsTools.PatAlgos.producersLayer1.ootPhotonProducer_cff')
process.load('PhysicsTools.PatAlgos.producersLayer1.metProducer_cff')
process.load('PhysicsTools.PatAlgos.producersLayer1.muonProducer_cff')
process.makePatJetsTask.add(process.pfImpactParameterTagInfos, 
                            process.pfSecondaryVertexTagInfos,
                            process.pfInclusiveSecondaryVertexFinderTagInfos)

process.patCandidatesTask = cms.Task(
    process.makePatElectronsTask,
    process.makePatMuonsTask,
    #process.makePatTausTask,
    process.makePatPhotonsTask,
    process.makePatOOTPhotonsTask,
    process.makePatJetsTask,
    process.makePatMETsTask
    )
process.patCandidates = cms.Sequence(process.patCandidatesTask)


process.load('PhysicsTools.PatAlgos.selectionLayer1.tauSelector_cfi')
process.load('PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi')
process.load('PhysicsTools.PatAlgos.selectionLayer1.electronSelector_cfi')
process.load('PhysicsTools.PatAlgos.selectionLayer1.muonSelector_cfi')
process.load('PhysicsTools.PatAlgos.selectionLayer1.photonSelector_cfi')
process.load('PhysicsTools.PatAlgos.selectionLayer1.ootPhotonSelector_cff')
process.selectedPatCandidatesTask = cms.Task(
    process.selectedPatElectrons,
    process.selectedPatMuons,
    #process.selectedPatTaus,
    process.selectedPatPhotons,
    process.selectedPatOOTPhotons,
    process.selectedPatJets,
 )
process.selectedPatCandidates = cms.Sequence(process.selectedPatCandidatesTask)


process.load('PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi')
process.patTrigger.onlyStandAlone = cms.bool(False)
process.patTrigger.packTriggerLabels = cms.bool(False)
process.patTrigger.packTriggerPathNames = cms.bool(False)
process.patTrigger.packTriggerPrescales = cms.bool(False) #True)

process.load('PhysicsTools.PatAlgos.slimming.selectedPatTrigger_cfi')
process.load('PhysicsTools.PatAlgos.slimming.slimmedPatTrigger_cfi')

process.patTask = cms.Task(
    process.patCandidatesTask,
    process.selectedPatCandidatesTask,
    #process.patTrigger,
    #process.selectedPatTrigger,
    #process.slimmedPatTrigger
)

# ----- Define Execution Paths ----- #

process.outputPath = cms.EndPath(process.output)
process.p = cms.Path(process.primaryVertexAssociation * process.egmGsfElectronIDSequence * process.egmPhotonIDSequence * process.NjettinessAK8CHS * process.metFilters * process.DisplacedHcalJets )
#process.p = cms.Path( process.ntuples )
process.schedule = cms.Schedule(process.p )

#Define Jet Tool Box Stuff
#listBtagDiscriminatorsAK4 = [ 
#                'pfJetProbabilityBJetTags',
#                'pfCombinedInclusiveSecondaryVertexV2BJetTags',
#                'pfCombinedMVAV2BJetTags',
#                'pfCombinedCvsLJetTags',
#                'pfCombinedCvsBJetTags',
#                ]
#from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
#jetToolbox( process, 'ak8', 'ak8JetSubs', "out", PUMethod='CHS', bTagDiscriminators=listBtagDiscriminatorsAK4, addSoftDrop=True, addNsub=True, addNsubSubjets=True, miniAOD=False )   ### For example


#Add PAT tasks for jet Toolbox to execution schedule
process.schedule.associate(process.patTask)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


#miniAOD_customize stuff
# process.patTaus.isoDeposits = cms.PSet()
# process.patTaus.addGenMatch = cms.bool(False)
# process.patTaus.embedGenMatch = cms.bool(False)
# process.patTaus.addGenJetMatch   = cms.bool(False)
# process.patTaus.embedGenJetMatch = cms.bool(False)
# process.patTaus.genParticleMatch = ''
# process.patTaus.genJetMatch      = ''
# process.selectedPatTaus.cut = cms.string("pt > 18. && tauID('decayModeFindingNewDMs')> 0.5")
# process.selectedPatJets.cut = cms.string("pt > 10")

## PU JetID
process.load("RecoJets.JetProducers.PileupJetID_cfi")
process.patTask.add(process.pileUpJetIDTask)
process.patJets.userData.userFloats.src = [ cms.InputTag("pileupJetId:fullDiscriminant"), ]
process.patJets.userData.userInts.src = [ cms.InputTag("pileupJetId:fullId"), ]

process.patJets.discriminatorSources = cms.VInputTag(
    cms.InputTag("pfJetBProbabilityBJetTags"),
    cms.InputTag("pfJetProbabilityBJetTags"),
    cms.InputTag("pfTrackCountingHighEffBJetTags"),
    cms.InputTag("pfSimpleSecondaryVertexHighEffBJetTags"),
    cms.InputTag("pfSimpleInclusiveSecondaryVertexHighEffBJetTags"),
    cms.InputTag("pfCombinedSecondaryVertexV2BJetTags"),
    cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
    cms.InputTag("softPFMuonBJetTags"),
    cms.InputTag("softPFElectronBJetTags"),
    cms.InputTag("pfCombinedMVAV2BJetTags"),   
    )
process.patJets.addTagInfos     = cms.bool(True)
process.patJets.tagInfoSources  = cms.VInputTag( 'pfImpactParameterTagInfos'
                                                 ,'pfSecondaryVertexTagInfos'
                                                 ,'pfInclusiveSecondaryVertexFinderTagInfos')
process.patJets.addGenPartonMatch   = cms.bool(False)
process.patJets.embedGenPartonMatch = cms.bool(False)
process.patJets.genPartonMatch      = ''
process.patJets.addGenJetMatch      = cms.bool(False)
process.patJets.embedGenJetMatch    = cms.bool(False)
process.patJets.genJetMatch         = ''
process.patJets.getJetMCFlavour    = cms.bool(False)
process.patJets.addJetFlavourInfo  = cms.bool(False)
process.patJets.JetPartonMapSource = ''
process.patJets.JetFlavourInfoSource = ''

# process.patJetsAK8PFCHS.addGenPartonMatch   = cms.bool(False)
# process.patJetsAK8PFCHS.embedGenPartonMatch = cms.bool(False)
# process.patJetsAK8PFCHS.genPartonMatch      = ''
# process.patJetsAK8PFCHS.addGenJetMatch      = cms.bool(False)
# process.patJetsAK8PFCHS.embedGenJetMatch    = cms.bool(False)
# process.patJetsAK8PFCHS.genJetMatch         = ''
# process.patJetsAK8PFCHS.getJetMCFlavour    = cms.bool(False)
# process.patJetsAK8PFCHS.addJetFlavourInfo  = cms.bool(False)
# process.patJetsAK8PFCHS.JetPartonMapSource = ''
# process.patJetsAK8PFCHS.JetFlavourInfoSource = ''
process.patMETs.addGenMET           = False
process.patMETs.genMETSource        = ''

