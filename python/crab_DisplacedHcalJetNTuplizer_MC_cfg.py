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
datasetnames = ['/HToSSTo4B_MH350_MS80_CTau500/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH350_MS80_CTau500_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER']# dataset wrapper 
'''
#'/HToSSTo4B_MH125_MS50_CTau3000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH125_MS50_CTau3000_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER'
#'/HToSSTo4B_MH125_MS50_CTau3000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH125_MS50_CTau3000_batch2_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER'
#'/HToSSTo4B_MH250_MS120_CTau10000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH250_MS120_CTau10000_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER'
#'/HToSSTo4B_MH250_MS120_CTau10000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH250_MS120_CTau10000_batch2_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER'
#'/HToSSTo4B_MH350_MS160_CTau10000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH350_MS160_CTau10000_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER'
##'/HToSSTo4B_MH350_MS160_CTau10000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH350_MS160_CTau10000_batch2_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER'
# '/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV/lpclonglived-crab_PrivateProduction_Summer22_DR_step2_RECOSIM_ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_batch1_v1-59a22edf0600a784f6c900595d24e883/USER'
'/HToSSTo4B_MH350_MS80_CTau500/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH350_MS80_CTau500_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER'
# /eos/uscms/store/user/lpclonglived/apresyan/privateProduction/DR/step2_RECOSIM/Run3Summer22/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV
]
'''

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
config.General.workArea        = '/afs/cern.ch/work/g/gkopp/2022_LLP_analysis/CMSSW_14_0_0/src/cms_lpc_llp/Run3-HCAL-LLP-NTupler/python/../../../../../crab_signalMC_20240905' # workArea wrapper
config.General.instance        = 'prod'
config.General.requestName     = 'HToSSTo4B_MH350_MS80_CTau500_USER_submission7_20240905_204933_version3' # requestName wrapper
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
config.Data.unitsPerJob      = 50 # with file based splitting, this is how many files in 1 job (50 * 100 events = 5k events, for 350 GeV MC. 25 * 200 events = 5k events for 125 GeV MC)
nJobs                        = 200
config.Data.totalUnits       = config.Data.unitsPerJob * nJobs # total number of units (not number of jobs!). 100 jobs * 50 units per job = 5000 units, with 100 events per unit * 5000 units = 500k events
config.Data.ignoreLocality   = True
config.Data.publication      = False
config.Data.outputDatasetTag = 'HToSSTo4B_MH350_MS80_CTau500_USER_submission7_20240905_204933_version3' # outputDatasetTag wrapper

config.Data.runRange        =  runrange
if lumimask != '':
  config.Data.lumiMask        = lumimask

config.Site.storageSite = storageSite

config.Site.whitelist = whiteList

if not blackList:
    config.Site.blacklist = blackList
