
from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = "CMSSW_9_4_7_WplusH_HToSSTobbbb_ms55_pl10000_DRstep1-Ntuple_CaltechT2"
config.General.workArea = "crab_prod"

config.section_("JobType")
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "/afs/cern.ch/work/c/christiw/public/LLP/CMSSW_9_4_7/src/cms_lpc_llp/llp_ntupler/python/displacedJetMuon_step2_ntupler_cfg.py"
#config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
config.Data.inputDataset = "/WplusH_HToSSTobbbb_ms55_pl10000_ev150000/sixie-crab_CMSSW_9_4_12_PrivateProduction_Fall17_DR_step1_WplusH_HToSSTobbbb_ms55_pl10000_v2_DR_CaltechT2-fb4de91b3672ad6012188656f7233fe2/USER"
config.Data.inputDBS = 'phys03'
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 1
config.Data.outputDatasetTag = "Run2"
#config.Data.outputDatasetTag = "Run2_displacedJetMuonNtupler_V1p9_MC_Summer16_jingyu-ggHdddd_M55_30mm_CP2_AODSIM-c3d6de13a4792afb4dd0c4ab58e49a3d_v1_v1"
config.Data.ignoreLocality = True

config.section_("Site")
config.Site.whitelist = ["T2_*"]
config.Site.storageSite = "T2_US_Caltech"
#config.Data.inputDBS = "https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader/"
config.Data.outLFNDirBase = '/store/group/phys_exotica/privateProduction/ntuple/RunIIFall17/WplusH_HToSSTobbbb_ms55_pl10000/v3/'
