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
number = 0 # starting at 0 -> refers to datasetnames

# List of possible datasets
datasetnames = [
'/JetMET/Run2022G-EXOHighMET-PromptReco-v1/RAW-RECO'
]

datasetblock = [
'/JetMET/Run2022G-EXOHighMET-PromptReco-v1/RAW-RECO#8469adfd-7170-448a-a270-211dc028c8d2'
]

# runrange = '362085,362087' # Nov2022 Phase Scan
runrange = '362696'

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
whiteList = ['T2_CH_CERN','T2_US_Caltech','T2_US_Florida', 'T2_US_MIT', 'T2_US_Nebraska', 'T2_US_Purdue', 'T2_US_UCSD', 'T2_US_Vanderbilt', 'T2_US_Wisconsin', 'T1_US_FNAL','T2_US_MIT']
# ['T2_US_UCSD']

# Black list sites
blackList = ['']

# -----------------------------------------------------------------------------------------------------------------------------
# No modifications below this line are necessary

import datetime
timestamp = datetime.datetime.now().strftime("_%Y%m%d_%H%M%S")


dataset = filter(None, datasetnames[number].split('/'))
dataset = list(dataset)

config = config()

# General
config.General.workArea        = 'crab_test_2023-05-08'
config.General.instance        = 'prod'
config.General.requestName     = dataset[0]+'_'+dataset[1]+'_'+dataset[2]+timestamp
config.General.transferOutputs = True
config.General.transferLogs    = True

# JobType
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = psetname
#config.JobType.outputFiles = [OutputFilename]
#config.JobType.pyCfgParams = ['outputFile='+OutputFilename]

# Data
config.Data.inputDataset     = datasetnames[number]
config.Data.inputBlocks      = datasetblock
config.Data.inputDBS         = 'global'
config.Data.splitting        = 'Automatic' #'LumiBased'
#config.Data.unitsPerJob      = 1
#config.Data.totalUnits       = 1
config.Data.ignoreLocality   = True
config.Data.publication      = False
config.Data.outputDatasetTag = dataset[1]+'_'+dataset[2]+timestamp

config.Data.runRange        =  runrange
if lumimask != '':
  config.Data.lumiMask        = lumimask

config.Site.storageSite = storageSite

config.Site.whitelist = whiteList

if not blackList:
    config.Site.blacklist = blackList
