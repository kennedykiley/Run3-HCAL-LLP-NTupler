# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: AOD_2023D --data --eventcontent AOD --datatier AOD --conditions 124X_dataRun3_Prompt_v4 --step RAW2DIGI,L1Reco,RECO --era Run3_2023 --geometry DB:Extended --customise Configuration/DataProcessing/RecoTLR.customisePrompt --filein file:/afs/cern.ch/work/k/kikenned/Run3-HCAL-LLP-NTupler/Run3-HCAL-LLP-NTupler/run/StudyFilesRAW/DisplacedJet_Run2023D-v1_RAW_369927_files/8b60c1c5-1a0e-4d7a-b9b3-9d197b9eac81.root --fileout file:AOD_2023D.root -n 5 --no_exec --python_filename=AOD_2023D_cfg.py --outputCommands keep HBHERecHitsSorted_hbhereco__RECO
import FWCore.ParameterSet.Config as cms

#from Configuration.Eras.Era_Run3_2023_cff import Run3_2023
from Configuration.Eras.Era_Run3_cff import Run3

#TODO: Other ERAS

# ---------------------------------------------------------------------------------------
# PARSE INPUTS 
# ---------------------------------------------------------------------------------------

import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing()

options.register('isData',
    True, # default value # isData wrapper
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.bool,
    "is Data"
)

options.register('isSignal',
    False, # default value # isSignal wrapper
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

options.register('tagJEC',
    "",
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,
    "String to help select JEC and JER tags"
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
print(f" recoFromRAW   ={options.recoFromRAW}")
print(f" skipEvents    ={options.skipEvents}")
print(f" processEvents ={options.processEvents}")
print(f" inputFiles    ={options.inputFiles}")
print(f" outputFile    ={options.outputFile}")
print(f" debug         ={options.debug}")
print(" ")

# ----- Load Input Files ----- #

inputFiles = options.inputFiles
i = 0
if len(options.inputFiles) == 1 and options.inputFiles[0][-4:] == ".txt": 
    print( "Reading input files from", options.inputFiles, ":" )

    file = open(options.inputFiles[0], 'r')
    lines = file.readlines()

    unpackedFilelist = []
    for line in lines:
        line_temp = line.replace(" ","")
        if line_temp[0] == "#": continue
        unpackedFilelist.append( line.strip() )
        if i < 10: print( "   ->", line.strip() )
        i += 1

    inputFiles = unpackedFilelist

if i >= 10: print( "Loaded", i, "files...")

# ----- Parse Inputs ----- #

if options.isData:
    mapping = {
        ("Run2022C", ""):               ("Summer22_22Sep2023_RunCD_V3_DATA",  "Summer22_22Sep2023_JRV1_DATA", ""),
        ("Run2022D", ""):               ("Summer22_22Sep2023_RunCD_V3_DATA",  "Summer22_22Sep2023_JRV1_DATA", ""),
        ("Run2022E", ""):               ("Summer22EE_22Sep2023_RunE_V3_DATA", "Summer22EE_22Sep2023_JRV1_DATA", ""),
        ("Run2022F", ""):               ("Summer22EE_22Sep2023_RunF_V3_DATA", "Summer22EE_22Sep2023_JRV1_DATA", ""),
        ("Run2022G", ""):               ("Summer22EE_22Sep2023_RunG_V3_DATA", "Summer22EE_22Sep2023_JRV1_DATA", ""),
        ("Run2023C", "PromptReco-v1"):  ("Summer23Prompt23_RunCv123_V3_DATA", "Summer23Prompt23_RunCv1234_JRV1_DATA", ""),
        ("Run2023C", "PromptReco-v2"):  ("Summer23Prompt23_RunCv123_V3_DATA", "Summer23Prompt23_RunCv1234_JRV1_DATA", ""),
        ("Run2023C", "PromptReco-v3"):  ("Summer23Prompt23_RunCv123_V3_DATA", "Summer23Prompt23_RunCv1234_JRV1_DATA", ""),
        ("Run2023C", "PromptReco-v4"):  ("Summer23Prompt23_RunCv4_V3_DATA",   "Summer23Prompt23_RunCv1234_JRV1_DATA", ""),
        ("Run2023D", ""):               ("Summer23BPixPrompt23_RunD_V3_DATA", "Summer23BPixPrompt23_RunD_JRV1_DATA", ""),
    }
else:
    mapping = { # TODO make sure this agrees with MC naming scheme
        ("HToSSTo4B", "23BPix"):         ("Summer23BPixPrompt23_V3_MC", "Summer23BPixPrompt23_RunD_JRV1_MC", "2023_Summer23BPix"),
        ("HToSSTo4B", "2023Prompt_"):    ("Summer23Prompt23_V3_MC", "Summer23Prompt23_RunCv1234_JRV1_MC", "2023_Summer23"), 
        ("HToSSTo4B", "2022EE"):         ("Summer22EE_22Sep2023_V3_MC", "Summer22EE_22Sep2023_JRV1_MC", "2022_Summer22EE"),
        ("HToSSTo4B", "2022_"):          ("Summer22_22Sep2023_V3_MC", "Summer22_22Sep2023_JRV1_MC", "2022_Summer22"),
        # ("HToSSTo4B", "2023BPixPrompt"): "Summer23BPixPrompt23_V3_MC",
        # ("WJetsToLNu", "preEE"):         ("Summer22_22Sep2023_V3_MC", "Summer22_22Sep2023_JRV1_MC"),
    }

tag_name         = None
JER_tag_name     = None
BTag_SF_tag_name = None
era_name         = None

for (run, reco), (name, JERname, BTAGname) in mapping.items():
    if run in options.tagJEC and reco in options.tagJEC:
        tag_name         = name
        JER_tag_name     = JERname
        BTag_SF_tag_name = BTAGname
        era_name         = run + "_" + reco
        break
    elif options.tagJEC == "" and run in inputFiles[0] and reco in inputFiles[0]:
        tag_name         = name
        JER_tag_name     = JERname
        BTag_SF_tag_name = BTAGname
        era_name         = run + "_" + reco
        break
if tag_name is None:
    raise RuntimeError("No matching JEC tag found for input file " + inputFiles[0])
if JER_tag_name is None:
    raise RuntimeError("No matching JER tag found for input file " + inputFiles[0])
if BTag_SF_tag_name is None:
    raise RuntimeError("No matching BTAG tag found for input file " + inputFiles[0])

if options.isData: 
    era_name = era_name.split("_")[0]
else:
    if "23BPix" in era_name:        era_name = "2023postBPix"
    elif "2023Prompt_" in era_name: era_name = "2023preBPix"
    elif "2022EE" in era_name:      era_name = "2022postEE"
    elif "2022_" in era_name:       era_name = "2022preEE"

# ---------------------------------------------------------------------------------------
# SET UP PROCESS
# ---------------------------------------------------------------------------------------

#process = cms.Process('RECO',Run3_2023) #TODOFIX

process = cms.Process('DisplacedHcalJetNTuplizer',Run3) #TODOFIX

#if recoFromRAW:

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

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.processEvents),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# ----- Input Source ----- #

#print( inputFiles )

process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/k/kikenned/Run3-HCAL-LLP-NTupler/Run3-HCAL-LLP-NTupler/run/StudyFilesRAW/DisplacedJet_Run2023D-v1_RAW_369927_files/8b60c1c5-1a0e-4d7a-b9b3-9d197b9eac81.root'),
    fileNames = cms.untracked.vstring(inputFiles),
    skipEvents = cms.untracked.uint32(options.skipEvents),
    secondaryFileNames = cms.untracked.vstring()
    #inputCommands = cms.untracked.vstring(
    #    "keep *",
    #    "drop *_gtStage2Digis_CICADAScore_*" # Run over 2024 data
    #),
)

if options.isData and True:  #and 'CRAB_JOB_ID' not in os.environ:: # Golden JSON passed in crab script instead, this is not used! 
    from FWCore.PythonUtilities.LumiList import LumiList
    goldenjson = None
    if   "Run2022" in era_name: goldenjson = LumiList(filename="/eos/user/c/cmsdqm/www/CAF/certification/Collisions22/Cert_Collisions2022_355100_362760_Golden.json")
    elif "Run2023" in era_name: goldenjson = LumiList(filename="/eos/user/c/cmsdqm/www/CAF/certification/Collisions23/Cert_Collisions2023_366442_370790_Golden.json")
    elif "Run2024" in era_name: goldenjson = LumiList(filename="/eos/user/c/cmsdqm/www/CAF/certification/Collisions24/Cert_Collisions2024_378981_386951_Golden.json")
    process.source.lumisToProcess = goldenjson.getVLuminosityBlockRange()

# ----- TFileService for output ----- #

process.TFileService = cms.Service( "TFileService",
    fileName = cms.string(options.outputFile),
    closeFileFast = cms.untracked.bool(True)
)

# ----- Process Options ----- #

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
    makeTriggerResults = cms.obsolete.untracked.bool,
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

# Meta Data
#process.configurationMetadata = cms.untracked.PSet(
#    annotation = cms.untracked.string('AOD_2023D nevts:5'),
#    name = cms.untracked.string('Applications'),
#    version = cms.untracked.string('$Revision: 1.19 $')
#)

# ---------------------------------------------------------------------------------------
# RECO COMPONENT
# ---------------------------------------------------------------------------------------

# ----- Global Tag ----- #
# MC: https://twiki.cern.ch/twiki/bin/view/CMSPublic/GTsRun3#Global_Tags_used_in_official_MC
# Data: https://twiki.cern.ch/twiki/bin/view/CMSPublic/GTsRun3#Global_Tags_used_in_official_AN1

from Configuration.AlCa.GlobalTag import GlobalTag

global_tags_MC = {
    "2022preEE":    "124X_mcRun3_2022_realistic_v12",
    "2022postEE":   "124X_mcRun3_2022_realistic_postEE_v1",
    "2023preBPix":  "130X_mcRun3_2023_realistic_v14",
    "2023postBPix": "130X_mcRun3_2023_realistic_v14",
}
# 140X_dataRun3_v17

if options.isData: 
    if "2022" in era_name: process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_Prompt_v4', '') 
    else:                  process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_v15', '') 
else:                      process.GlobalTag = GlobalTag(process.GlobalTag, global_tags_MC[era_name], '') 

# ----- HLT Filter ----- #

process.hltFilter = cms.EDFilter("HLTHighLevel",
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),  # "HLT" = process name
    HLTPaths = cms.vstring("HLT_*L1SingleLLPJet*"),
    andOr = cms.bool(True),
    throw = cms.bool(False)
)

# ----- Actual Reconstruction!!!! ----- #

process.raw2digi_step = cms.Path( process.hltFilter * process.RawToDigi)
process.L1Reco_step = cms.Path( process.hltFilter * process.L1Reco)
process.reconstruction_step = cms.Path(process.hltFilter * process.reconstruction)

# ---------------------------------------------------------------------------------------
# PAT Stuff
# ---------------------------------------------------------------------------------------

# ------ Declare Output ------ #

# TFileService for output
#process.TFileService = cms.Service( "TFileService",
#    fileName = cms.string(options.outputFile), #"file:ntuple_2023D_passRECO_fromRAW.root"),
#    closeFileFast = cms.untracked.bool(True)
#)

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

# ----- PAT Candidates Task ----- #

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

# ----- Selected PAT Candidates Task ----- #

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

# ----- PAT Trigger ----- #

process.load('PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi')
process.patTrigger.onlyStandAlone = cms.bool(False)
process.patTrigger.packTriggerLabels = cms.bool(False)
process.patTrigger.packTriggerPathNames = cms.bool(False)
process.patTrigger.packTriggerPrescales = cms.bool(False) #True)

process.load('PhysicsTools.PatAlgos.slimming.selectedPatTrigger_cfi')
process.load('PhysicsTools.PatAlgos.slimming.slimmedPatTrigger_cfi')

# ----- PAT Task ----- #

process.patTask = cms.Task(
    process.patCandidatesTask,
    process.selectedPatCandidatesTask,
    #process.patTrigger,
    #process.selectedPatTrigger,
    #process.slimmedPatTrigger
)

# ----- AK8 Jets ----- #

# Add jettiness for AK8 jets

process.load('RecoJets.JetProducers.nJettinessAdder_cfi')
process.NjettinessAK8CHS = process.Njettiness.clone()

# ----- E/Gamma Object IDs ----- #

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

# ----- PAT Jets (More) ----- #

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
process.pileupJetId.jets = cms.InputTag("ak4PFJetsCHS") # Point pileupJetId to the initial jets 
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

# ---------------------------------------------------------------------------------------
# JEC and JER Stuff
# ---------------------------------------------------------------------------------------

# ----- Jet Energy Corrections and Jet Energy Resolution for MC ----- # GK
# JEC from sqlite file # following here: 
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JecSqliteFile
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyResolution#Accessing_factors_from_Global_Ta
# Tag from https://cms-conddb.cern.ch/cmsDbBrowser/list/Prod/gts/140X_dataRun3_v17, and use "conddb --db <.db file> listTags" to confirm name from .db file
# This is to use a .db file to over-ride the conditions from the GT 
# https://cms-jerc.web.cern.ch/Recommendations/#2022

from CondCore.CondDB.CondDB_cfi import CondDB

# ---------- JEC -----------
JEC_file_path = 'sqlite_file:JEC_JER/JECDatabase/SQLiteFiles/' + tag_name + '.db'
if options.tagJEC == "": # this is if processing locally, since with CRAB tagJEC is filled
    JEC_file_path = 'sqlite_file:../data/JEC_JER/JECDatabase/SQLiteFiles/' + tag_name + '.db'
CondDBJECFile = CondDB.clone(connect = cms.string(JEC_file_path))
CollectionName = 'JetCorrectorParametersCollection_' + tag_name + '_AK4PFchs'
PuppiCollectionName = 'JetCorrectorParametersCollection_' + tag_name + '_AK4PFPuppi'
process.jec = cms.ESSource('PoolDBESSource',
    CondDBJECFile,
    toGet = cms.VPSet(
        # CHS jets first
        cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string(CollectionName), 
            label  = cms.untracked.string('AK4PFchs')
        ),
        # PUPPI jets
        cms.PSet(
            record = cms.string("JetCorrectionsRecord"),
            tag    = cms.string(PuppiCollectionName),
            label  = cms.untracked.string("AK4PFPuppi")
        )
    )
)
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource', 'jec')

# # ---------- JER ----------- # actually read in from .txt files below 
# from JetMETCorrections.Modules.JetResolutionESProducer_cfi import *
# JER_file_path = 'sqlite:../data/JEC_JER/JRDatabase/SQLiteFiles/' + JER_tag_name + '.db'
# CondDBJERFile = CondDB.clone(connect = cms.string(JER_file_path))
# CollectionName = 'JR_' + JER_tag_name + '_SF_AK4PFchs' 
# PuppiCollectionName = 'JR_' + JER_tag_name + '_SF_AK4PFPuppi'
# CollectionName_Pt = 'JR_' + JER_tag_name + '_PtResolution_AK4PFchs' 
# PuppiCollectionName_Pt = 'JR_' + JER_tag_name + '_PtResolution_AK4PFPuppi'
# process.jer = cms.ESSource('PoolDBESSource',
#     CondDBJERFile,
#     toGet = cms.VPSet(
#         # CHS jets
#         cms.PSet(
#             record = cms.string('JetResolutionRcd'),
#             tag    = cms.string(CollectionName_Pt),
#             label  = cms.untracked.string('AK4PFchs_pt')
#         ),
#         cms.PSet(
#             record = cms.string('JetResolutionScaleFactorRcd'),
#             tag    = cms.string(CollectionName),
#             label  = cms.untracked.string('AK4PFchs')
#         ),
#         # PUPPI jets
#         cms.PSet(
#             record = cms.string('JetResolutionRcd'),
#             tag    = cms.string(PuppiCollectionName_Pt),
#             label  = cms.untracked.string('AK4PFPuppi_pt')
#         ),
#         cms.PSet(
#             record = cms.string('JetResolutionScaleFactorRcd'),
#             tag    = cms.string(PuppiCollectionName),
#             label  = cms.untracked.string('AK4PFPuppi')
#         )
#     )
# )
# process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')

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
from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import patJets
# GK PUPPI
from PhysicsTools.PatAlgos.JetCorrFactorsProducer_cfi import JetCorrFactorsProducer
process.patJetCorrFactorsPuppi = JetCorrFactorsProducer.clone(
    src = cms.InputTag("ak4PFJetsPuppi"),  # input RECO jets
    levels = jecLevels,  # JEC levels
    payload = 'AK4PFPuppi'  # must match the label in PoolDBESSource
)
# Step 1: make PAT jets from RECO PUPPI jets. Step 2: add correction factors, done inside here otherwise had issues with "updateJetCollection" about gen particles
# --- PUPPI PAT jets (no JEC, no gen, safe for RECO/AOD) --- #
process.patJetsPuppi = patJets.clone(
    jetSource = cms.InputTag("ak4PFJetsPuppi"),
    addBTagInfo         = cms.bool(True),
    addDiscriminators   = cms.bool(False),
    # discriminatorSources = cms.VInputTag(
    #     cms.InputTag("pfDeepCSVJetTags:probb"),
    #     cms.InputTag("pfDeepCSVJetTags:probc"),
    #     cms.InputTag("pfDeepCSVJetTags:probudsg"),
    #     cms.InputTag("pfDeepCSVJetTags:probbb")
    # ),
    addJetCorrFactors   = cms.bool(True),
    jetCorrFactorsSource = cms.VInputTag(
        cms.InputTag('patJetCorrFactorsPuppi')  # JEC for PUPPI jets
    ),
    addAssociatedTracks = cms.bool(False),
    addGenPartonMatch   = cms.bool(False),
    addGenJetMatch      = cms.bool(False),
    getJetMCFlavour     = cms.bool(False),
    addJetFlavourInfo   = cms.bool(False),
    addJetCharge        = cms.bool(False),
    JetFlavourInfoSource = cms.InputTag(""), # doesn't seem to exist at this level
    JetPartonMapSource  = cms.InputTag(""),
    trackAssociationSource = cms.InputTag(""),
    embedGenPartonMatch = cms.bool(False),
    embedGenJetMatch    = cms.bool(False),
    embedPFCandidates   = cms.bool(False),
    embedCaloTowers     = cms.bool(False),
    addTagInfos         = cms.bool(True)
)


# ------ BTAG Scale Factors and Uncertainties ------ #


# ---------------------------------------------------------------------------------------
# Attach the NTupler
# ---------------------------------------------------------------------------------------

# ------ Import MET Filter ------ #

process.load('cms_lpc_llp.Run3-HCAL-LLP-NTupler.metFilters_Run3_cff')

# ------ Analyzer ------ #

process.DisplacedHcalJets = cms.EDAnalyzer('DisplacedHcalJetNTuplizer',
    debug =  cms.bool( options.debug ),
    isData = cms.bool( options.isData ),
    isSignal = cms.bool( options.isSignal ),
    era    = cms.string( era_name ),
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
    triggerPathNamesFile = cms.string("cms_lpc_llp/Run3-HCAL-LLP-NTupler/data/HLTPathsLLPJetsHCAL.dat"), # note this is not used! hard-coded in .cc file instead
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
    # GK PUPPI
    # jets = cms.InputTag("ak4PFJetsPuppi"),
    pfjetsAK4Puppi = cms.InputTag("patJetsPuppi"), # selectedUpdatedPatJetsPuppiUpdatedJEC"), # use the second one if updateJetCollection is used
    #jetsPuppi = cms.InputTag("ak4PFJets"),
    #jetsAK8 = cms.InputTag("ak8PFJetsCHS"),
    #jetsAK8 = cms.InputTag("selectedPatJetsAK8PFCHS"),

    #jetsAK8 = cms.InputTag("ak8PFJetsPuppi"),

    #mets = cms.InputTag("slimmedMETs"),
    met = cms.InputTag("patMETs"), # this is used in the .cc file as "met"
    metPuppi = cms.InputTag("pfMetPuppi"), # this is used in the .cc file as "metPuppi"
    #metsNoHF = cms.InputTag("pfMet30"),
    # metsPuppi = cms.InputTag("pfMet"),
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

# ----- Jet Energy Resolution (only MC), and JEC uncertainties (data and MC) ----- # GK

if not options.isData: 
    process.DisplacedHcalJets.jer_PtResolution = cms.FileInPath("cms_lpc_llp/Run3-HCAL-LLP-NTupler/data/JEC_JER/JRDatabase/textFiles/"+JER_tag_name+"/"+JER_tag_name+"_PtResolution_AK4PFPuppi.txt")
    process.DisplacedHcalJets.jer_ScaleFactor = cms.FileInPath("cms_lpc_llp/Run3-HCAL-LLP-NTupler/data/JEC_JER/JRDatabase/textFiles/"+JER_tag_name+"/"+JER_tag_name+"_SF_AK4PFPuppi.txt")
    process.DisplacedHcalJets.btagSysSF = cms.FileInPath("cms_lpc_llp/Run3-HCAL-LLP-NTupler/data/BTag/btv-scale-factors/"+BTag_SF_tag_name+"/json/btagging_v2.json")
    # note that _AK4PFchs.txt is a symlink back to _AK4PFPuppi.txt, and there are issues when the symlinked version is used. So the puppi version is listed. 
    # based on https://cms-jerc.web.cern.ch/Recommendations/#2023_1
# process.DisplacedHcalJets.jec_Uncertainty = cms.FileInPath("cms_lpc_llp/Run3-HCAL-LLP-NTupler/data/JEC_JER/JECDatabase/textFiles/"+tag_name+"/"+tag_name+"_Uncertainty_AK4PFPuppi.txt")


# ------ MET Filter Handling for RECO From RAW ------ #

if options.recoFromRAW:

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

    process.DisplacedHcalJets.metFilterBits = cms.InputTag("")


# ---------------------------------------------------------------------------------------
# Final Processing: Path, EndPath, Schedule, Associate
# ---------------------------------------------------------------------------------------

# ------ Standard NTupler Processing ------ #

# Output
process.outputPath = cms.EndPath(process.output)

# NTuplization Process
process.p = cms.Path( process.primaryVertexAssociationLocal * process.egmGsfElectronIDSequence * process.egmPhotonIDSequence * process.NjettinessAK8CHS * process.metFiltersRecommended * process.DisplacedHcalJets )

# Schedule
process.schedule = cms.Schedule( process.p )

if options.recoFromRAW:

    # Include HLT Filter to NTuplization Process
    process.p = cms.Path( process.hltFilter * process.primaryVertexAssociationLocal * process.egmGsfElectronIDSequence * process.egmPhotonIDSequence * process.NjettinessAK8CHS * process.metFiltersRecommended * process.DisplacedHcalJets )

    # Schedule definition
    process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,
        process.Flag_HBHENoiseFilter,process.Flag_HBHENoiseIsoFilter,process.Flag_CSCTightHaloFilter,process.Flag_CSCTightHaloTrkMuUnvetoFilter,process.Flag_CSCTightHalo2015Filter,process.Flag_globalTightHalo2016Filter,process.Flag_globalSuperTightHalo2016Filter,process.Flag_HcalStripHaloFilter,process.Flag_hcalLaserEventFilter,process.Flag_EcalDeadCellTriggerPrimitiveFilter,process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_ecalBadCalibFilter,process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,process.Flag_trkPOGFilters,process.Flag_chargedHadronTrackResolutionFilter,process.Flag_muonBadTrackFilter,process.Flag_BadChargedCandidateFilter,process.Flag_BadPFMuonFilter,process.Flag_BadPFMuonDzFilter,process.Flag_hfNoisyHitsFilter,process.Flag_BadChargedCandidateSummer16Filter,process.Flag_BadPFMuonSummer16Filter,process.Flag_trkPOG_manystripclus53X,process.Flag_trkPOG_toomanystripclus53X,process.Flag_trkPOG_logErrorTooManyClusters,
        process.p) #, #process.endjob_step,process.AODoutput_step)

# ------ Associate PAT Task ------ #

process.schedule.associate(process.patTask)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# ------ Customization ------ #

# Automatic addition of the customisation function from Configuration.DataProcessing.RecoTLR
from Configuration.DataProcessing.RecoTLR import customisePrompt 

#call to customisation function customisePrompt imported from Configuration.DataProcessing.RecoTLR
process = customisePrompt(process)

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)