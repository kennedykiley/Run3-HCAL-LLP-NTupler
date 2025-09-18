import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# ------ Setup ------ #

#initialize the process
from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process("DisplacedHcalJetNTuplizer", Run3) # line added to fix PCastorRcd error
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")

# MET Filter Recommendations: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Run_3_2022_and_2023_data_and_MC
process.load('cms_lpc_llp.Run3-HCAL-LLP-NTupler.metFilters_Run3_cff')

# Fix GEM Error
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag

# ----- Parse Arguments ----- #

options = VarParsing.VarParsing()

options.register('isData',
    False, # default value # isData wrapper
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.bool,
    "is Data"
)

options.register('isSignal',
    True, # default value # isSignal wrapper
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.bool,
    "is Signal"
)

options.register('recoFromRAW',
    False,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.bool,
    "run reconstruction from RAW (only for *L1SingleLLPJet* triggers)"
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

# ----- Load Input Files ----- #

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
    skipEvents = cms.untracked.uint32(options.skipEvents),
    inputCommands = cms.untracked.vstring(
        "keep *",
        "drop *_gtStage2Digis_CICADAScore_*" # Run over 2024 data
    )
)

# TFileService for output
process.TFileService = cms.Service( "TFileService",
    fileName = cms.string(options.outputFile),
    closeFileFast = cms.untracked.bool(True)
)

# ----- Apply Golden Json ----- #

# Actually do this in CRAB job...
# from /eos/user/c/cmsdqm/www/CAF/certification/Collisions*

from FWCore.PythonUtilities.LumiList import LumiList
#import os

if False: #options.isData and 'CRAB_JOB_ID' not in os.environ::
    goldenjson = None
    if   "Run2022" in inputFiles[0]: goldenjson = LumiList(filename="../data/certification/Cert_Collisions2022_355100_362760_Golden.json")
    elif "Run2023" in inputFiles[0]: goldenjson = LumiList(filename="../data/certification/Cert_Collisions2023_366442_370790_Golden.json")
    elif "Run2024" in inputFiles[0]: goldenjson = LumiList(filename="../data/certification/Cert_Collisions2024_378981_386951_Golden.json")
    process.source.lumisToProcess = goldenjson.getVLuminosityBlockRange()

# ----- Load Run Conditions ----- #

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('Configuration.Geometry.GeometryIdeal_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')

# ------ Declare the correct global tag ------ #

if options.isData: process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_v15', '') 
else:              process.GlobalTag = GlobalTag(process.GlobalTag,'130X_mcRun3_2023_realistic_v14','')
# referenced from here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/GTsRun3 
# else:              process.GlobalTag = GlobalTag(process.GlobalTag,'auto:run3_mc_FULL','')

# ------ Declare Output ------ #

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
 
# ------ Reco things from RAW (needed for 2022) ------ #

from FWCore.ParameterSet.Modules import _Module
class ModuleCollector:
    def __init__(self, exclude_names=None):
        self.modules = []
        self.exclude_names = set(exclude_names) if exclude_names else set()
    def enter(self, obj):
        if isinstance(obj, _Module) and obj.label() not in self.exclude_names:
            self.modules.append(obj)
        return True
    def leave(self, obj):
        return True

def flattenSequence(seq, exclude_names=None):
    collector = ModuleCollector(exclude_names)
    seq.visit(collector)
    return collector.modules

if options.recoFromRAW: # Works for 13_2_0

    # import of standard configurations
    process.load('Configuration.StandardSequences.Services_cff')
    process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
    process.load('FWCore.MessageService.MessageLogger_cfi')
    process.load('Configuration.EventContent.EventContent_cff')
    process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
    process.load('Configuration.StandardSequences.MagneticField_cff')
    process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
    process.load('Configuration.StandardSequences.L1Reco_cff')
    process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
    process.load('Configuration.StandardSequences.EndOfProcess_cff')
    process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

    # additional
    process.load('RecoMET.METFilters.metFilters_cff')

    process.options = cms.untracked.PSet(
        FailPath = cms.untracked.vstring(),
        IgnoreCompletely = cms.untracked.vstring(),
        Rethrow = cms.untracked.vstring(),
        SkipEvent = cms.untracked.vstring(),
        accelerators = cms.untracked.vstring('*'),
        allowUnscheduled = cms.obsolete.untracked.bool,
        canDeleteEarly = cms.untracked.vstring(),
        deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
        dumpOptions = cms.untracked.bool(False),
        emptyRunLumiMode = cms.obsolete.untracked.string,
        eventSetup = cms.untracked.PSet(
            forceNumberOfConcurrentIOVs = cms.untracked.PSet(
                allowAnyLabel_=cms.required.untracked.uint32
            ),
            numberOfConcurrentIOVs = cms.untracked.uint32(0)
        ),
        fileMode = cms.untracked.string('FULLMERGE'),
        forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
        holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
        #makeTriggerResults = cms.obsolete.untracked.bool,
        modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
        numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
        numberOfConcurrentRuns = cms.untracked.uint32(1),
        numberOfStreams = cms.untracked.uint32(0),
        numberOfThreads = cms.untracked.uint32(1),
        printDependencies = cms.untracked.bool(False),
        sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
        throwIfIllegalParameter = cms.untracked.bool(True),
        wantSummary = cms.untracked.bool(False)
    )    

    # Production Info
    process.configurationMetadata = cms.untracked.PSet(
        annotation = cms.untracked.string('reco_from_raw_temp nevts:100'),
        name = cms.untracked.string('Applications'),
        version = cms.untracked.string('$Revision: 1.19 $')
    )

    process.hltFilter = cms.EDFilter("HLTHighLevel",
        TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),  # "HLT" = process name
        HLTPaths = cms.vstring("HLT_*L1SingleLLPJet*"),
        andOr = cms.bool(True),
        throw = cms.bool(False)
    )

    process.Flag_BadChargedCandidateFilter = cms.Path( process.hltFilter * process.BadChargedCandidateFilter )
    process.Flag_BadChargedCandidateSummer16Filter = cms.Path( process.hltFilter * process.BadChargedCandidateSummer16Filter )
    process.Flag_BadPFMuonDzFilter = cms.Path( process.hltFilter * process.BadPFMuonDzFilter )
    process.Flag_BadPFMuonFilter = cms.Path( process.hltFilter * process.BadPFMuonFilter )
    process.Flag_BadPFMuonSummer16Filter = cms.Path( process.hltFilter * process.BadPFMuonSummer16Filter )
    process.Flag_CSCTightHalo2015Filter = cms.Path( process.hltFilter * process.CSCTightHalo2015Filter )
    process.Flag_CSCTightHaloFilter = cms.Path( process.hltFilter * process.CSCTightHaloFilter )
    process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path( process.hltFilter * process.CSCTightHaloTrkMuUnvetoFilter )
    process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path( process.hltFilter * process.EcalDeadCellBoundaryEnergyFilter )
    process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path( process.hltFilter * process.EcalDeadCellTriggerPrimitiveFilter )
    process.Flag_HBHENoiseFilter = cms.Path( process.hltFilter * process.HBHENoiseFilterResultProducer+process.HBHENoiseFilter )
    process.Flag_HBHENoiseIsoFilter = cms.Path( process.hltFilter * process.HBHENoiseFilterResultProducer+process.HBHENoiseIsoFilter )
    process.Flag_HcalStripHaloFilter = cms.Path( process.hltFilter * process.HcalStripHaloFilter )
    process.Flag_chargedHadronTrackResolutionFilter = cms.Path( process.hltFilter * process.chargedHadronTrackResolutionFilter )
    process.Flag_ecalBadCalibFilter = cms.Path()
    process.Flag_ecalLaserCorrFilter = cms.Path( process.hltFilter * process.ecalLaserCorrFilter )
    process.Flag_eeBadScFilter = cms.Path( process.hltFilter * process.eeBadScFilter )
    process.Flag_globalSuperTightHalo2016Filter = cms.Path( process.hltFilter * process.globalSuperTightHalo2016Filter )
    process.Flag_globalTightHalo2016Filter = cms.Path( process.hltFilter * process.globalTightHalo2016Filter )
    process.Flag_goodVertices = cms.Path(process.hltFilter * process.primaryVertexFilter )
    process.Flag_hcalLaserEventFilter = cms.Path(process.hltFilter * process.hcalLaserEventFilter )
    process.Flag_hfNoisyHitsFilter = cms.Path(process.hltFilter * process.hfNoisyHitsFilter )
    process.Flag_muonBadTrackFilter = cms.Path(process.hltFilter * process.muonBadTrackFilter )
    process.Flag_trackingFailureFilter = cms.Path(process.hltFilter * process.goodVertices+process.trackingFailureFilter )
    process.Flag_trkPOGFilters = cms.Path(process.hltFilter * process.trkPOGFilters )
    process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(process.hltFilter * ~process.logErrorTooManyClusters )
    process.Flag_trkPOG_manystripclus53X = cms.Path(process.hltFilter * ~process.manystripclus53X )
    process.Flag_trkPOG_toomanystripclus53X = cms.Path(process.hltFilter * ~process.toomanystripclus53X )

# ------ Custom Additions ------ #

# For AOD Track variables
"""
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
"""

# ----- Jet Energy Corrections ----- # GK
# JEC from sqlite file # following here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JecSqliteFile
# from CondCore.CondDB.CondDB_cfi import CondDB
# CondDBJECFile = CondDB.clone(connect = cms.string( 'sqlite:Summer23Prompt23_RunCv123_V3_DATA.db' ) ) # file must be local! 
## Error: Tag "JetCorrectorParametersCollection_Run2_Run3_DATA_AK4PFchs_offline_v3" has not been found in the database. from IOVProxy::load 
# process.jec = cms.ESSource('PoolDBESSource',
#     CondDBJECFile,
#     toGet = cms.VPSet(
#         cms.PSet(
#             record = cms.string('JetCorrectionsRecord'),
#             tag    = cms.string('JetCorrectorParametersCollection_Run2_Run3_DATA_AK4PFchs_offline_v3'), # from https://cms-conddb.cern.ch/cmsDbBrowser/list/Prod/gts/140X_dataRun3_v17
#             label  = cms.untracked.string('AK4PFchs')
#         ),
#     )
# )
# process.es_prefer_jec = cms.ESPrefer('PoolDBESSource', 'jec')

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

# Define correction levels
if options.isData:
    jecLevels = ['L1FastJet','L2Relative','L3Absolute','L2L3Residual']
else:
    jecLevels = ['L1FastJet','L2Relative','L3Absolute']

# Update PAT jets with JEC
updateJetCollection(
    process,
    jetSource = cms.InputTag('selectedPatJets'),
    labelName = 'UpdatedJEC',
    jetCorrections = ('AK4PFchs', cms.vstring(jecLevels), 'None'),
    pvSource = cms.InputTag('offlinePrimaryVertices') # specifically point to the PV collection, in RECO need the offline PVs (instead of offlineSlimmedPrimaryVertices used at miniAOD)
    # svSource = cms.InputTag('slimmedSecondaryVertices'),
)
# now we have a new collection, selectedUpdatedPatJetsUpdatedJEC = corrected jets. As before, selectedPatJets = still the uncorrected PAT jets.

# ------ Analyzer ------ #

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
    vertices = cms.InputTag("offlinePrimaryVertices"), #"", "RECO"),

    electrons = cms.InputTag("gedGsfElectrons"),
    muons = cms.InputTag("muons"),
    taus = cms.InputTag("selectedPatTaus"),
    photons = cms.InputTag("gedPhotons"),
    pfjetsAK4 = cms.InputTag("selectedPatJets"),
    calojetsAK4 = cms.InputTag("ak4CaloJets"), #,"","RECO"),
    pfjetsAK8 = cms.InputTag("selectedPatJetsAK8PFCHS"),
    calojetsAK8 = cms.InputTag("ak8CaloJets"), #,"","RECO"),
    l1jets = cms.InputTag("gtStage2Digis","Jet"), #,"RECO"), # GK, added for L1 jets access
    pfjetsAK4_corrected = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC"), # GK, adding JECs # Point analyzer to corrected jets
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
    pfCands = cms.InputTag("particleFlow"), #,"","RECO"),

    #packedPfCands = cms.InputTag("packedPFCandidates"),

    genParticles = cms.InputTag("genParticles"),

    puInfo = cms.InputTag("addPileupInfo", "", "HLT"), #uncomment if no pre-mixing
    #puInfo = cms.InputTag("mixData", "", "HLT"), #uncomment for samples with pre-mixed pileup

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

    metFilterBits = cms.InputTag("TriggerResults",  "", "RECO"), 

    #hbheNoiseFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"),
    #hbheTightNoiseFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun2Tight"),
    #hbheIsoNoiseFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHEIsoNoiseFilterResult"),

    #BadChargedCandidateFilter = cms.InputTag("BadChargedCandidateFilter",""),
    #BadMuonFilter = cms.InputTag("BadPFMuonFilter",""),

    #lheInfo = cms.InputTag("externalLHEProducer", "", ""),
    genInfo = cms.InputTag("generator", "", "SIM"),

    tracks = cms.InputTag("generalTracks"), #, "", "RECO"),
    #trackTime = cms.InputTag("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModel"),
    #trackTimeReso = cms.InputTag("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModelResolution"),

    #hcalNoiseInfo = cms.InputTag("hcalnoise", "", "RECO"),

    #secondaryVertices = cms.InputTag("inclusiveSecondaryVertices", "", "RECO"),
    secondaryVertices = cms.InputTag("inclusiveCandidateSecondaryVertices"), #,"", "RECO"),

    rhoAll = cms.InputTag("fixedGridRhoAll"), #, "", "RECO"),

    rhoFastjetAll = cms.InputTag("fixedGridRhoFastjetAll"), # "", "RECO"),
    rhoFastjetAllCalo = cms.InputTag("fixedGridRhoFastjetAllCalo", "", "RECO"),
    rhoFastjetCentralCalo = cms.InputTag("fixedGridRhoFastjetCentralCalo", "", "RECO"),
    rhoFastjetCentralChargedPileUp = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp", "", "RECO"),
    rhoFastjetCentralNeutral = cms.InputTag("fixedGridRhoFastjetCentralNeutral", "", "RECO"),

    beamSpot = cms.InputTag("offlineBeamSpot"), #, "", "RECO"),
    pfClusters = cms.InputTag("particleFlowClusterECAL"), #,"","RECO"),
    #hbRecHits = cms.InputTag("reducedHcalRecHits", "hbhereco","RECO"),
    hbRecHits = cms.InputTag("hbhereco"), #, "","RECO"),
    #ebRecHits = cms.InputTag("EcalRecHit", "reducedEcalRecHitsEB", "RECO"),
    #ebRecHits = cms.InputTag("ecalRecHit", "EcalRecHitsEB", "RECO"), # GK, errors with HCAL LLP skim, as with below
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

    Flag_HBHENoiseFilter = cms.InputTag("Flag_HBHENoiseFilter"),
    Flag_HBHENoiseIsoFilter = cms.InputTag("Flag_HBHENoiseIsoFilter"),
    Flag_CSCTightHaloFilter = cms.InputTag("Flag_CSCTightHaloFilter"),
    Flag_CSCTightHaloTrkMuUnvetoFilter = cms.InputTag("Flag_CSCTightHaloTrkMuUnvetoFilter"),
    Flag_CSCTightHalo2015Filter = cms.InputTag("Flag_CSCTightHalo2015Filter"),
    Flag_globalTightHalo2016Filter = cms.InputTag("Flag_globalTightHalo2016Filter"),
    Flag_globalSuperTightHalo2016Filter = cms.InputTag("Flag_globalSuperTightHalo2016Filter"),
    Flag_HcalStripHaloFilter = cms.InputTag("Flag_HcalStripHaloFilter"),
    Flag_hcalLaserEventFilter = cms.InputTag("Flag_hcalLaserEventFilter"),
    Flag_EcalDeadCellTriggerPrimitiveFilter = cms.InputTag("Flag_EcalDeadCellTriggerPrimitiveFilter"),
    Flag_EcalDeadCellBoundaryEnergyFilter = cms.InputTag("Flag_EcalDeadCellBoundaryEnergyFilter"),
    Flag_ecalBadCalibFilter = cms.InputTag("Flag_ecalBadCalibFilter"),
    Flag_goodVertices = cms.InputTag("Flag_goodVertices"),
    Flag_eeBadScFilter = cms.InputTag("Flag_eeBadScFilter"),
    Flag_ecalLaserCorrFilter = cms.InputTag("Flag_ecalLaserCorrFilter"),
    Flag_trkPOGFilters = cms.InputTag("Flag_trkPOGFilters"),
    Flag_chargedHadronTrackResolutionFilter = cms.InputTag("Flag_chargedHadronTrackResolutionFilter"),
    Flag_muonBadTrackFilter = cms.InputTag("Flag_muonBadTrackFilter"),
    Flag_BadChargedCandidateFilter = cms.InputTag("Flag_BadChargedCandidateFilter"),
    Flag_BadPFMuonFilter = cms.InputTag("Flag_BadPFMuonFilter"),
    Flag_BadChargedCandidateSummer16Filter = cms.InputTag("Flag_BadChargedCandidateSummer16Filter"),
    Flag_BadPFMuonSummer16Filter = cms.InputTag("Flag_BadPFMuonSummer16Filter"),
    Flag_BadPFMuonDzFilter = cms.InputTag("Flag_BadPFMuonDzFilter"),
    Flag_hfNoisyHitsFilter = cms.InputTag("Flag_hfNoisyHitsFilter"),
    Flag_trkPOG_manystripclus53X = cms.InputTag("Flag_trkPOG_manystripclus53X"),
    Flag_trkPOG_toomanystripclus53X = cms.InputTag("Flag_trkPOG_toomanystripclus53X"),
    Flag_trkPOG_logErrorTooManyClusters = cms.InputTag("Flag_trkPOG_logErrorTooManyClusters"),
)

# ----- Add Additional Info ----- #

# ----- Jet Energy Resolution (only MC), and JEC uncertainties (data and MC) ----- # GK
if not options.isData: 
    process.DisplacedHcalJets.jer_PtResolution = cms.FileInPath("cms_lpc_llp/Run3-HCAL-LLP-NTupler/python/Summer23Prompt23_RunCv1234_JRV1_MC_PtResolution_AK4PFchs.txt")
    process.DisplacedHcalJets.jer_ScaleFactor  = cms.FileInPath("cms_lpc_llp/Run3-HCAL-LLP-NTupler/python/Summer23Prompt23_RunCv1234_JRV1_MC_SF_AK4PFchs.txt")
    # from https://github.com/cms-jet/JRDatabase/blob/master/textFiles/Summer23Prompt23_RunCv1234_JRV1_MC/
    # TODO need to handle different eras 
    # based on https://cms-jerc.web.cern.ch/Recommendations/#2023_1
    process.DisplacedHcalJets.jec_Uncertainty = cms.FileInPath("cms_lpc_llp/Run3-HCAL-LLP-NTupler/python/Summer23Prompt23_V3_MC_Uncertainty_AK4PFchs.txt")
    # from https://github.com/cms-jet/JECDatabase/blob/master/textFiles/Summer23Prompt23_V1_MC/Summer23Prompt23_V1_MC_Uncertainty_AK4PFPuppi.txt
    # V3 -> V2 -> V1 -> Puppi
if options.isData:
    process.DisplacedHcalJets.jec_Uncertainty = cms.FileInPath("cms_lpc_llp/Run3-HCAL-LLP-NTupler/python/Summer23Prompt23_RunCv4_V3_DATA_Uncertainty_AK4PFchs.txt")
    # same link as for MC, different folder, and DATA -> MC
    # Test FileInPath (doesn't like symlinks) and cms.string (still complains "JetCorrectorParameters: No definitions found!!!")
    # if     "Run2023C" in inputFiles[0] and ("PromptReco-v1" in inputFiles[0] or "PromptReco-v2" in inputFiles[0] or "PromptReco-v3" in inputFiles[0]): 
    #     process.DisplacedHcalJets.jec_Uncertainty = cms.string("cms_lpc_llp/Run3-HCAL-LLP-NTupler/data/JEC_JER/JECDatabase/textFiles/Summer23Prompt23_RunCv123_V3_DATA/Summer23Prompt23_RunCv123_V3_DATA_Uncertainty_AK4PFchs.txt")
    #     print("using Cv123 JEC files!")
    # elif   "Run2023C" in inputFiles[0] and "PromptReco-v4" in inputFiles[0]: 
    #     process.DisplacedHcalJets.jec_Uncertainty = cms.string("cms_lpc_llp/Run3x-HCAL-LLP-NTupler/data/JEC_JER/JECDatabase/textFiles/Summer23Prompt23_RunCv4_V3_DATA/Summer23Prompt23_RunCv4_V3_DATA_Uncertainty_AK4PFchs.txt")
    #     print("using Cv4 JEC files!")
    # elif   "Run2023D" in inputFiles[0]: process.DisplacedHcalJets.jec_Uncertainty = cms.string("cms_lpc_llp/Run3-HCAL-LLP-NTupler/data/JEC_JER/JECDatabase/textFiles/Summer23BPixPrompt23_RunD_V3_DATA/Summer23BPixPrompt23_RunD_V3_DATA_Uncertainty_AK4PFchs.txt")
    # elif   "Run2023D" in inputFiles[0]: process.DisplacedHcalJets.jec_Uncertainty = cms.string("cms_lpc_llp/Run3-HCAL-LLP-NTupler/data/JEC_JER/JECDatabase/textFiles/Summer23BPixPrompt23_RunD_V2xHFscale_DATA/Summer23BPixPrompt23_RunD_V2xHFscale_DATA_Uncertainty_AK4PFchs.txt")
    # TODO determine which one of Run D we should use
        
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
process.primaryVertexAssociationLocal = process.sortedPFPrimaryVertices.clone(
    qualityForPrimary = cms.int32(2),
    produceSortedVertices = cms.bool(False),
    producePileUpCollection  = cms.bool(False),  
    produceNoPileUpCollection = cms.bool(False)
)

# ----- PAT Stuff ----- #

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

# ----- miniAOD_customize stuff ----- #

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


#from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
#updateJetCollection(
#    process,
#    jetSource = cms.InputTag('updatedPatJetsTransientCorrected'), #slimmedJets'),
#    pvSource = cms.InputTag('offlinePrimaryVertices'),
#    svSource = cms.InputTag('slimmedSecondaryVertices'),
#    jetCorrections = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
#    btagDiscriminators = [
#        'pfParticleTransformerAK4JetTags:probb',
#    ],
#    postfix = 'SlimmedDeepFlavour',
#)

process.load('PhysicsTools.PatAlgos.recoLayer0.bTagging_cff')
process.patJets.discriminatorSources = cms.VInputTag(
    #cms.InputTag('pfParticleTransformerAK4JetTags:probb'),# 'probb'),
    cms.InputTag("pfDeepCSVJetTags:probb"),
    cms.InputTag("pfDeepCSVJetTags:probc"),
    cms.InputTag("pfDeepCSVJetTags:probudsg"),
    cms.InputTag("pfDeepCSVJetTags:probbb"),
    #cms.InputTag("pfParticleTransformerAK4JetTags:probb"),
    #cms.InputTag("pfJetBProbabilityBJetTags"),
    #cms.InputTag("pfJetProbabilityBJetTags"),
    #cms.InputTag("pfTrackCountingHighEffBJetTags"),
    # cms.InputTag("pfSimpleSecondaryVertexHighEffBJetTags"), # GK, errors with HCAL LLP skim, as with below 7
    # cms.InputTag("pfSimpleInclusiveSecondaryVertexHighEffBJetTags"),
    # cms.InputTag("pfCombinedSecondaryVertexV2BJetTags"),
    # cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
    # cms.InputTag("softPFMuonBJetTags"),
    # cms.InputTag("softPFElectronBJetTags"),
    # cms.InputTag("pfCombinedMVAV2BJetTags"),   
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

# ----- Define Execution Paths ----- #

process.outputPath = cms.EndPath(process.output)

if options.recoFromRAW:
    # Update metFilterBits, which are passed via the trigger object
    process.DisplacedHcalJets.metFilterBits = cms.InputTag("")

    # Schedule noise filters
    process.schedule = cms.Schedule(process.Flag_HBHENoiseFilter,process.Flag_HBHENoiseIsoFilter,process.Flag_CSCTightHaloFilter,process.Flag_CSCTightHaloTrkMuUnvetoFilter,process.Flag_CSCTightHalo2015Filter,process.Flag_globalTightHalo2016Filter,process.Flag_globalSuperTightHalo2016Filter,process.Flag_HcalStripHaloFilter,process.Flag_hcalLaserEventFilter,process.Flag_EcalDeadCellTriggerPrimitiveFilter,process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_ecalBadCalibFilter,process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,process.Flag_trkPOGFilters,process.Flag_chargedHadronTrackResolutionFilter,process.Flag_muonBadTrackFilter,process.Flag_BadChargedCandidateFilter,process.Flag_BadPFMuonFilter,process.Flag_BadPFMuonDzFilter,process.Flag_hfNoisyHitsFilter,process.Flag_BadChargedCandidateSummer16Filter,process.Flag_BadPFMuonSummer16Filter,process.Flag_trkPOG_manystripclus53X,process.Flag_trkPOG_toomanystripclus53X,process.Flag_trkPOG_logErrorTooManyClusters)
    process.RECO     = cms.Path( process.hltFilter * process.RawToDigi * process.gtStage2Digis  * process.reconstruction ) #* process.metFiltersAll )
    process.p        = cms.Path( process.hltFilter * process.primaryVertexAssociationLocal * process.egmGsfElectronIDSequence * process.egmPhotonIDSequence * process.NjettinessAK8CHS * process.metFiltersRecommended * process.DisplacedHcalJets )

    if False: # set to true when you want to dump all reco objects
        process.dumpEverything = cms.EDAnalyzer("EventContentAnalyzer")
        process.p = cms.Path(process.dumpEverything)

    process.schedule.insert(0, process.RECO)
    process.schedule.append(process.p)

elif False:
    process.DisplacedHcalJets.hbRecHits = cms.InputTag("")
    process.DisplacedHcalJets.pfCands   = cms.InputTag("") #patPackedCandidates") #packedPFCandidates")
    process.DisplacedHcalJets.vertices  = cms.InputTag("offlineSlimmedPrimaryVertices")
    process.DisplacedHcalJets.met       = cms.InputTag("slimmedMETs")
    process.DisplacedHcalJets.electrons = cms.InputTag("slimmedElectrons")
    process.DisplacedHcalJets.muons     = cms.InputTag("slimmedMuons")
    process.DisplacedHcalJets.photons   = cms.InputTag("slimmedPhotons")
    process.DisplacedHcalJets.pfjetsAK4 = cms.InputTag("slimmedJets")
    process.DisplacedHcalJets.l1jet     = cms.InputTag("caloStage2Digis")

    #process.metFiltersRecommended.primaryVertexFilter.src = cms.InputTag("offlineSlimmedPrimaryVertices")

    process.p = cms.Path( process.metFiltersRecommended_MINIAOD * process.DisplacedHcalJets )
    process.schedule = cms.Schedule( process.p )

else:
    process.p = cms.Path( process.primaryVertexAssociationLocal * process.egmGsfElectronIDSequence * process.egmPhotonIDSequence * process.NjettinessAK8CHS * process.metFiltersRecommended * process.DisplacedHcalJets )
    process.schedule = cms.Schedule( process.p )

#Add PAT tasks for jet Toolbox to execution schedule
#if True:
process.schedule.associate(process.patTask)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
