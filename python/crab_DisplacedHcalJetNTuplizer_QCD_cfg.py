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
#'/QCD_PT-15to7000_TuneCP5_13p6TeV_pythia8/Run3Winter23Reco-FlatPU0to120GTv3_GTv3_126X_mcRun3_2023_forPU65_v3-v2/GEN-SIM-RECO'
'/QCD_PT-15to7000_TuneCP5_13p6TeV_pythia8/Run3Winter23Reco-FlatPU0to120_126X_mcRun3_2023_forPU65_v1-v2/GEN-SIM-RECO'
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
whiteList = ['T1_UK_RAL', 'T2_CH_CERN','T2_US_Caltech','T2_US_Florida', 'T2_US_MIT', 'T2_US_Nebraska', 'T2_US_Purdue', 'T2_US_UCSD', 'T2_US_Vanderbilt', 'T2_US_Wisconsin', 'T1_US_FNAL','T2_US_MIT']
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
#config.General.workArea        = 'crab_QCD-MC_2023-10-09'
config.General.workArea        = 'crab_QCD_MC'+date
config.General.instance        = 'prod'
config.General.requestName     = 'QCD_MC_'+timestamp #dataset[0]+'_'+dataset[1]+'_'+dataset[2]+timestamp
config.General.transferOutputs = True
config.General.transferLogs    = True

# JobType
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = psetname
#config.JobType.outputFiles = [OutputFilename]
#config.JobType.pyCfgParams = ['outputFile='+OutputFilename]

# Data
config.Data.inputDataset     = datasetnames[number]
config.Data.inputDBS         = 'global' # because input files are from another crab production run #'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/' #'global'
config.Data.splitting        = 'FileBased' #'Automatic' #'LumiBased'
config.Data.unitsPerJob      = 10 # with file based splitting, this is how many files in 1 job (5 * 200 events = 1k events)
nJobs                        = 150
config.Data.totalUnits       = config.Data.unitsPerJob * nJobs # total number of units (not number of jobs!). 100 jobs * 5 units per job = 500 units, with 200 events per unit * 500 units = 100k events
config.Data.ignoreLocality   = True
config.Data.publication      = False
config.Data.outputDatasetTag = 'QCD_MC_'+timestamp #dataset[1]+'_'+dataset[2]+timestamp

config.Data.runRange        =  runrange
if lumimask != '':
  config.Data.lumiMask        = lumimask

config.Site.storageSite = storageSite

config.Site.whitelist = whiteList

if not blackList:
    config.Site.blacklist = blackList
