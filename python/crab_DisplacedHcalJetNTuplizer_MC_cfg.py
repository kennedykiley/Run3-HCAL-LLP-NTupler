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
# '/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV/lpclonglived-crab_PrivateProduction_Summer22_DR_step2_RECOSIM_ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_batch1_v1-59a22edf0600a784f6c900595d24e883/USER'
'/HToSSTo4B_MH350_MS80_CTau500/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH350_MS80_CTau500_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER'
# /eos/uscms/store/user/lpclonglived/apresyan/privateProduction/DR/step2_RECOSIM/Run3Summer22/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV
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
whiteList = ['T2_CH_CERN','T2_US_Caltech','T2_US_Florida', 'T2_US_MIT', 'T2_US_Nebraska', 'T2_US_Purdue', 'T2_US_UCSD', 'T2_US_Vanderbilt', 'T2_US_Wisconsin', 'T1_US_FNAL','T2_US_MIT']
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
#config.General.workArea        = 'crab_signalMC_2023-06-29'
config.General.workArea        = 'crab_signalMC'+date
config.General.instance        = 'prod'
config.General.requestName     = 'LLP_MC_350_'+timestamp #dataset[0]+'_'+dataset[1]+'_'+dataset[2]+timestamp
config.General.transferOutputs = True
config.General.transferLogs    = True

# JobType
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = psetname
#config.JobType.outputFiles = [OutputFilename]
#config.JobType.pyCfgParams = ['outputFile='+OutputFilename]

# Data
config.Data.inputDataset     = datasetnames[number]
config.Data.inputDBS         = 'phys03' # because input files are from another crab production run #'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/' #'global'
config.Data.splitting        = 'FileBased' #'Automatic' #'LumiBased'
config.Data.unitsPerJob      = 50 # with file based splitting, this is how many files in 1 job (50 * 100 events = 5k events, for 350 GeV MC. 200 events for 125 GeV MC)
nJobs                        = 100
config.Data.totalUnits       = config.Data.unitsPerJob * nJobs # total number of units (not number of jobs!). 100 jobs * 50 units per job = 5000 units, with 100 events per unit * 5000 units = 500k events
config.Data.ignoreLocality   = True
config.Data.publication      = False
config.Data.outputDatasetTag = 'LLP_MC_350_'+timestamp #dataset[1]+'_'+dataset[2]+timestamp

config.Data.runRange        =  runrange
if lumimask != '':
  config.Data.lumiMask        = lumimask

config.Site.storageSite = storageSite

config.Site.whitelist = whiteList

if not blackList:
    config.Site.blacklist = blackList
