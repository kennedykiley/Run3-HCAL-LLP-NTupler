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
datasetnames = ['MYVAR_DATASET_NAME'] # dataset wrapper

datasetblock = [
#'/DisplacedJet/Run2023C-EXOLLPJetHCAL-PromptReco-v4/AOD#0882cc9a-f2ab-4626-992e-c787a9d5017c' # in 2023C v4
]

# runrange = '362085,362087' # Nov2022 Phase Scan
runrange = ''

# JSON files for lumiMask are available at: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/
lumimask = 'MYVAR_GOLDEN_JSON'

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
blackList = ['T2_BE_UCL']

# -----------------------------------------------------------------------------------------------------------------------------
# No modifications below this line are necessary

import datetime
timestamp = datetime.datetime.now().strftime("_%Y%m%d_%H%M%S")
date = datetime.datetime.now().strftime("_%Y%m%d")

dataset = filter(None, datasetnames[number].split('/'))
dataset = list(dataset)

config = config()

# Arguments
config.JobType.pyCfgParams = [
    'isData=MYVAR_ISDATA',
    'isSignal=MYVAR_ISSIGNAL',
    'recoFromRAW=MYVAR_RECO_FROM_RAW',
    'tagJEC=MYVAR_DATASET_NAME'
]

# General
config.General.workArea        = 'MYVAR_CRAB_OUTPUT_NAME' #'/afs/cern.ch/work/k/kikenned/Run3-HCAL-LLP-NTupler/CRAB_Workarea/NTuples_v4/'
config.General.instance        = 'prod'
config.General.requestName     = "MYVAR_REQUEST_NAME"
config.General.transferOutputs = True
config.General.transferLogs    = True

# JobType
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = psetname
config.JobType.inputFiles  = ["../data/JEC_JER"]
#config.JobType.outputFiles = [OutputFilename]
#config.JobType.pyCfgParams = ['outputFile='+OutputFilename]

# Data
# four below lines for standard dataset input
config.Data.inputDataset     = datasetnames[number]
config.Data.inputBlocks      = datasetblock
config.Data.inputDBS         = 'MY_VAR_INPUTDBS'
config.Data.splitting        = 'Automatic' #'LumiBased'
# for single file test (3 below lines)
#config.Data.userInputFiles    = ['/store/data/Run2023C/DisplacedJet/AOD/EXOLLPJetHCAL-PromptReco-v4/000/367/881/00000/36ade28b-f320-4680-9dab-57ce2b536531.root']
#config.Data.splitting         = 'FileBased' #'EventAwareLumiBased' #'FileBased'
#config.Data.unitsPerJob       = 400 #MYVAR_EVENTS_PER_FILE
#config.Data.totalUnits       = 1
config.Data.ignoreLocality   = True
config.Data.publication      = False
config.Data.outputDatasetTag = "MYVAR_DATASET_TAG"

# MYVAR_EXTRACONFIG

config.Data.runRange        =  runrange
if lumimask != '':
  config.Data.lumiMask        = lumimask

config.Site.storageSite = storageSite

config.Site.whitelist = whiteList

if not blackList:
    config.Site.blacklist = blackList
