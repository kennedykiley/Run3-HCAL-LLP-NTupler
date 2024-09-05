# -----------------------------------------------------------------------------------------------------------------------------
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial  #**Up-to-date**
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile #
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3CheatSheet        #
#
# Environment setup:
#    cmsenv
#    source /cvmfs/cms.cern.ch/crab3/crab.sh
# To submit:
#    crab submit -c MyCrabConfig_Splash.py
# To check status:
#    crab status -d <CRAB-project-directory> [--jobids <comma-separated-list-of-jobs-and/or-job-ranges>]
# To kill jobs:
#    crab kill -d <CRAB-project-directory> [--jobids <comma-separated-list-of-jobs-and/or-job-ranges>]
# To retrieve output:
#    crab getoutput -d <CRAB-project-directory> [--jobids <comma-separated-list-of-jobs-and/or-job-ranges>]
# -----------------------------------------------------------------------------------------------------------------------------
from CRABClient.UserUtilities import config
#from CRABClient.UserUtilities import getUsernameFromSiteDB

# Select dataset to crab over
number = 0 # starting at 0 -> refers to datasetnames # number wrapper

# List of possible datasets
datasetnames = ['/Muon0/Run2023C-ZMu-PromptReco-v4/RAW-RECO']# dataset wrapper
"""
'/DisplacedJet/Run2023B-EXOLLPJetHCAL-PromptReco-v1/AOD', # 0
'/DisplacedJet/Run2023C-EXOLLPJetHCAL-PromptReco-v1/AOD',
'/DisplacedJet/Run2023C-EXOLLPJetHCAL-PromptReco-v2/AOD',
'/DisplacedJet/Run2023C-EXOLLPJetHCAL-PromptReco-v3/AOD',
'/DisplacedJet/Run2023C-EXOLLPJetHCAL-PromptReco-v4/AOD',
'/DisplacedJet/Run2023D-EXOLLPJetHCAL-PromptReco-v1/AOD',
'/DisplacedJet/Run2023D-EXOLLPJetHCAL-PromptReco-v2/AOD',
'/JetMET1/Run2023A-EXOHighMET-PromptReco-v2/RAW-RECO',    # 7
'/JetMET1/Run2023B-EXOHighMET-PromptReco-v1/RAW-RECO',
'/JetMET1/Run2023C-EXOHighMET-PromptReco-v1/RAW-RECO',
'/JetMET1/Run2023C-EXOHighMET-PromptReco-v2/RAW-RECO',
'/JetMET1/Run2023C-EXOHighMET-PromptReco-v3/RAW-RECO',
'/JetMET1/Run2023C-EXOHighMET-PromptReco-v4/RAW-RECO',
'/JetMET1/Run2023D-EXOHighMET-PromptReco-v1/RAW-RECO',
'/JetMET1/Run2023D-EXOHighMET-PromptReco-v2/RAW-RECO', # 14
'/Muon0/Run2023A-ZMu-PromptReco-v2/RAW-RECO', # no events, do not submit! 
'/Muon0/Run2023B-ZMu-PromptReco-v1/RAW-RECO',
'/Muon0/Run2023C-ZMu-PromptReco-v1/RAW-RECO',
'/Muon0/Run2023C-ZMu-PromptReco-v2/RAW-RECO',
'/Muon0/Run2023C-ZMu-PromptReco-v3/RAW-RECO',
'/Muon0/Run2023C-ZMu-PromptReco-v4/RAW-RECO',
'/Muon0/Run2023D-ZMu-PromptReco-v1/RAW-RECO', # 21
'/Muon0/Run2023D-ZMu-PromptReco-v2/RAW-RECO',
'/Muon1/Run2023A-ZMu-PromptReco-v2/RAW-RECO',
'/Muon1/Run2023B-ZMu-PromptReco-v1/RAW-RECO',
'/Muon1/Run2023C-ZMu-PromptReco-v1/RAW-RECO',
'/Muon1/Run2023C-ZMu-PromptReco-v2/RAW-RECO',
'/Muon1/Run2023C-ZMu-PromptReco-v3/RAW-RECO',
'/Muon1/Run2023C-ZMu-PromptReco-v4/RAW-RECO', # 28
'/Muon1/Run2023D-ZMu-PromptReco-v1/RAW-RECO',
'/Muon1/Run2023D-ZMu-PromptReco-v2/RAW-RECO'
]
"""

datasetblock = [
#'/DisplacedJet/Run2023C-EXOLLPJetHCAL-PromptReco-v4/AOD#0882cc9a-f2ab-4626-992e-c787a9d5017c' # in 2023C v4
]

# runrange = '362085,362087' # Nov2022 Phase Scan
runrange = ''

# JSON files for lumiMask are available at: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/
lumimask = ''

# Storage path for output files - EOS specific
#storagepath = '/store/user/'+getUsernameFromSiteDB()+'/HCALnoise2016'

# cmsRun file
psetname = 'DisplacedHcalJetNTuplizer.py'

# Output filename
#OutputFilename = '/eos/home-k/kikenned/HCALtupleMaker/CrabOutput_Run'+runrange+'.root'

# Storage site of output files
storageSite = 'T2_US_Wisconsin' # no write access to: 'T2_CH_CERN'

# White list sites
whiteList = ['T2_CH_CERN','T2_US_Caltech','T2_US_Florida', 'T2_US_MIT', 'T2_US_Nebraska', 'T2_US_Purdue', 'T2_US_UCSD', 'T2_US_Vanderbilt', 'T2_US_Wisconsin', 'T1_US_FNAL','T2_US_MIT','T1_FR_CCIN2P3']
# ['T2_US_UCSD']

# Black list sites
blackList = ['']

# -----------------------------------------------------------------------------------------------------------------------------
# No modifications below this line are necessary

import datetime
timestamp = datetime.datetime.now().strftime("_%Y%m%d_%H%M%S")
date = datetime.datetime.now().strftime("_%Y%m%d")

dataset = filter(None, datasetnames[number].split('/'))
dataset = list(dataset)

config = config()

# General
config.General.workArea        = '/afs/cern.ch/work/g/gkopp/2022_LLP_analysis/CMSSW_14_0_0/src/cms_lpc_llp/Run3-HCAL-LLP-NTupler/python/../../../../../crab_ZMu_20240905' # workArea wrapper
config.General.instance        = 'prod'
config.General.requestName     = 'Muon0_Run2023C-ZMu-PromptReco-v4_RAW-RECO_20240905_204857_version3' # requestName wrapper
config.General.transferOutputs = True
config.General.transferLogs    = True

# JobType
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = psetname
#config.JobType.outputFiles = [OutputFilename]
#config.JobType.pyCfgParams = ['outputFile='+OutputFilename]

# Data
# four below lines for standard dataset input
config.Data.inputDataset     = datasetnames[number]
config.Data.inputBlocks      = datasetblock
config.Data.inputDBS         = 'global'
config.Data.splitting        = 'Automatic' #'LumiBased'
# for single file test (3 below lines)
#config.Data.userInputFiles    = ['/store/data/Run2023C/DisplacedJet/AOD/EXOLLPJetHCAL-PromptReco-v4/000/367/881/00000/36ade28b-f320-4680-9dab-57ce2b536531.root']
#config.Data.splitting         = 'FileBased'
#config.Data.unitsPerJob       = 1
#config.Data.totalUnits       = 1
config.Data.ignoreLocality   = True
config.Data.publication      = False
config.Data.outputDatasetTag = 'Run2023C-ZMu-PromptReco-v4_RAW-RECO_20240905_204857_version3' # outputDatasetTag wrapper

config.Data.runRange        =  runrange
if lumimask != '':
  config.Data.lumiMask        = lumimask

config.Site.storageSite = storageSite

config.Site.whitelist = whiteList

if not blackList:
    config.Site.blacklist = blackList
