datasets = {}


datasets["Data_DisplacedJet_Run2022"] = [
    #"/DisplacedJet/Run2022A-v1/RAW", #  900 Gev
    #"/DisplacedJet/Run2022B-v1/RAW", # Commissioning
    "/DisplacedJet/Run2022C-v1/RAW",
    "/DisplacedJet/Run2022D-v1/RAW",  # There are D-v2 and -v3, but they're not listed in DAS 
    "/DisplacedJet/Run2022E-v1/RAW",
    "/DisplacedJet/Run2022F-v1/RAW",
    "/DisplacedJet/Run2022G-v1/RAW",
]

datasets["Data_EXOLLPJetHCAL_Run2023"] = [
    "/DisplacedJet/Run2023B-EXOLLPJetHCAL-PromptReco-v1/AOD",
    "/DisplacedJet/Run2023C-EXOLLPJetHCAL-PromptReco-v1/AOD",
    "/DisplacedJet/Run2023C-EXOLLPJetHCAL-PromptReco-v2/AOD",
    "/DisplacedJet/Run2023C-EXOLLPJetHCAL-PromptReco-v3/AOD",
    "/DisplacedJet/Run2023C-EXOLLPJetHCAL-PromptReco-v4/AOD",
    "/DisplacedJet/Run2023D-EXOLLPJetHCAL-PromptReco-v1/AOD",
    "/DisplacedJet/Run2023D-EXOLLPJetHCAL-PromptReco-v2/AOD",

]

datasets["Data_EXOLLPJetHCAL_Run2024"] = [
    #"/DisplacedJet/Run2024A-EXOLLPJetHCAL-PromptReco-v1/AOD",
    "/DisplacedJet/Run2024B-EXOLLPJetHCAL-PromptReco-v1/AOD",
    "/DisplacedJet/Run2024C-EXOLLPJetHCAL-PromptReco-v1/AOD",
    "/DisplacedJet/Run2024D-EXOLLPJetHCAL-PromptReco-v1/AOD",
    "/DisplacedJet/Run2024E-EXOLLPJetHCAL-PromptReco-v1/AOD",
    "/DisplacedJet/Run2024E-EXOLLPJetHCAL-PromptReco-v2/AOD",
    "/DisplacedJet/Run2024F-EXOLLPJetHCAL-PromptReco-v1/AOD",
    "/DisplacedJet/Run2024G-EXOLLPJetHCAL-PromptReco-v1/AOD",
    "/DisplacedJet/Run2024H-EXOLLPJetHCAL-PromptReco-v1/AOD",
    "/DisplacedJet/Run2024I-EXOLLPJetHCAL-PromptReco-v1/AOD",
    "/DisplacedJet/Run2024I-EXOLLPJetHCAL-PromptReco-v2/AOD",
]

datasets["Data_ZMu_Run2022"] = [
    #"/Muon/Run2022C-ZMu-10Dec2022-v1/RAW-RECO",
    #"/Muon/Run2022C-ZMu-27Jun2023-v1/RAW-RECO",
    "/Muon/Run2022C-ZMu-PromptReco-v1/RAW-RECO",
    #"/Muon/Run2022D-ZMu-10Dec2022-v1/RAW-RECO",
    #"/Muon/Run2022D-ZMu-27Jun2023-v2/RAW-RECO",
    "/Muon/Run2022D-ZMu-PromptReco-v1/RAW-RECO",
    "/Muon/Run2022D-ZMu-PromptReco-v2/RAW-RECO",
    "/Muon/Run2022D-ZMu-PromptReco-v3/RAW-RECO",
    #"/Muon/Run2022E-ZMu-10Dec2022-v2/RAW-RECO",
    #"/Muon/Run2022E-ZMu-27Jun2023-v1/RAW-RECO",
    "/Muon/Run2022E-ZMu-PromptReco-v1/RAW-RECO",
    "/Muon/Run2022F-ZMu-PromptReco-v1/RAW-RECO",
    "/Muon/Run2022G-ZMu-PromptReco-v1/RAW-RECO",
]

datasets["Data_ZMu_Run2023"] = [
    '/Muon0/Run2023B-ZMu-PromptReco-v1/RAW-RECO',
    '/Muon0/Run2023C-ZMu-PromptReco-v1/RAW-RECO',
    '/Muon0/Run2023C-ZMu-PromptReco-v2/RAW-RECO',
    '/Muon0/Run2023C-ZMu-PromptReco-v3/RAW-RECO',
    '/Muon0/Run2023C-ZMu-PromptReco-v4/RAW-RECO',
    '/Muon0/Run2023D-ZMu-PromptReco-v1/RAW-RECO',
    '/Muon0/Run2023D-ZMu-PromptReco-v2/RAW-RECO',
    #'/Muon1/Run2023A-ZMu-PromptReco-v2/RAW-RECO',
    '/Muon1/Run2023B-ZMu-PromptReco-v1/RAW-RECO',
    '/Muon1/Run2023C-ZMu-PromptReco-v1/RAW-RECO',
    '/Muon1/Run2023C-ZMu-PromptReco-v2/RAW-RECO',
    '/Muon1/Run2023C-ZMu-PromptReco-v3/RAW-RECO',
    '/Muon1/Run2023C-ZMu-PromptReco-v4/RAW-RECO',
    '/Muon1/Run2023D-ZMu-PromptReco-v1/RAW-RECO',
    '/Muon1/Run2023D-ZMu-PromptReco-v2/RAW-RECO'
]

#datasets["Signal_HToSSTo4B_MH125_MS15_CTau1000"] = [
#    "/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV/lpclonglived-crab_PrivateProduction_Summer22_DR_step2_RECOSIM_ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_batch1_v1-59a22edf0600a784f6c900595d24e883/USER"
#]

# private datasets from here: https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset%3D%2FHToSSTo4B_MH*_MS*_CTau*%2Flpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH*_MS*_CTau*_batch*%2FUSER+instance%3Dprod%2Fphys03
datasets["Signal_HToSSTo4B_MH125_MS50_CTau3000"] = [
    "/HToSSTo4B_MH125_MS50_CTau3000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH125_MS50_CTau3000_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER",
    "/HToSSTo4B_MH125_MS50_CTau3000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH125_MS50_CTau3000_batch2_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER"
]

datasets["Signal_HToSSTo4B_MH125_MS50_CTau10000"] = [
    "/HToSSTo4B_MH125_MS50_CTau10000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH125_MS50_CTau10000_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER"
]

datasets["Signal_HToSSTo4B_MH125_MS15_CTau3000"] = [
    "/HToSSTo4B_MH125_MS15_CTau3000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH125_MS15_CTau3000_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER"
]
datasets["Signal_HToSSTo4B_MH125_MS15_CTau10000"] = [
    "/HToSSTo4B_MH125_MS15_CTau10000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH125_MS15_CTau10000_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER"
]

datasets["Signal_HToSSTo4B_MH250_MS120_CTau10000"] = [
    "/HToSSTo4B_MH250_MS120_CTau10000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH250_MS120_CTau10000_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER",
    "/HToSSTo4B_MH250_MS120_CTau10000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH250_MS120_CTau10000_batch2_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER"
]

datasets["Signal_HToSSTo4B_MH350_MS80_CTau500"] = [
    "/HToSSTo4B_MH350_MS80_CTau500/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH350_MS80_CTau500_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER"
]

datasets["Signal_HToSSTo4B_MH350_MS160_CTau10000"] = [
    "/HToSSTo4B_MH350_MS160_CTau10000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH350_MS160_CTau10000_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER",
    "/HToSSTo4B_MH350_MS160_CTau10000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH350_MS160_CTau10000_batch2_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER"
]

# central signal production samples, DAS = https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset%3D%2FHTo2LongLivedTo4b_MH-*_MFF-*_CTau-*mm_TuneCP5_13p6TeV_pythia8%2FRun3*-124X_mcRun3_*%2FGEN-SIM-RAW
datasets["CentralSignal_RAW_HToSSTo4B_MH125_MS15_CTau10000_2022"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-15_CTau-10000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]
datasets["CentralSignal_RAW_HToSSTo4B_MH125_MS15_CTau10000_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-15_CTau-10000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]
datasets["CentralSignal_RAW_HToSSTo4B_MH125_MS15_CTau3000_2022"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-15_CTau-3000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]
datasets["CentralSignal_RAW_HToSSTo4B_MH125_MS15_CTau3000_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-15_CTau-3000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]
datasets["CentralSignal_RAW_HToSSTo4B_MH125_MS23_CTau10000_2022"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-23_CTau-10000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_20２２_realistic_v1２-v２/GEN-SIM-RAW"
]
datasets["CentralSignal_RAW_HToSSTo4B_MH125_MS23_CTau10000_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-23_CTau-10000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]
datasets["CentralSignal_RAW_HToSSTo4B_MH125_MS23_CTau3000_2022"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-23_CTau-3000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]
datasets["CentralSignal_RAW_HToSSTo4B_MH125_MS23_CTau3000_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-23_CTau-3000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]
datasets["CentralSignal_RAW_HToSSTo4B_MH125_MS30_CTau10000_2022"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-30_CTau-10000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]
datasets["CentralSignal_RAW_HToSSTo4B_MH125_MS30_CTau10000_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-30_CTau-10000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]
datasets["CentralSignal_RAW_HToSSTo4B_MH125_MS30_CTau3000_2022"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-30_CTau-3000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]
datasets["CentralSignal_RAW_HToSSTo4B_MH125_MS30_CTau3000_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-30_CTau-3000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]
datasets["CentralSignal_RAW_HToSSTo4B_MH125_MS40_CTau10000_2022"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-40_CTau-10000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]
datasets["CentralSignal_RAW_HToSSTo4B_MH125_MS40_CTau10000_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-40_CTau-10000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]
datasets["CentralSignal_RAW_HToSSTo4B_MH125_MS40_CTau3000_2022"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-40_CTau-3000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]
datasets["CentralSignal_RAW_HToSSTo4B_MH125_MS40_CTau3000_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-40_CTau-3000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]
datasets["CentralSignal_RAW_HToSSTo4B_MH125_MS55_CTau10000_2022"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-55_CTau-10000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]
datasets["CentralSignal_RAW_HToSSTo4B_MH125_MS55_CTau10000_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-55_CTau-10000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]
datasets["CentralSignal_RAW_HToSSTo4B_MH125_MS55_CTau3000_2022"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-55_CTau-3000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]
datasets["CentralSignal_RAW_HToSSTo4B_MH125_MS55_CTau3000_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-55_CTau-3000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]

# data used for background studies 
datasets["Background_WPlusJets"] = [
    "/WJetsToLNu_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Winter23Reco-TRKRealistic_AlcaRecoRealisticTRK_preEE_126X_mcRun3_2022_realistic_v4-v2/GEN-SIM-RECO",
    "/WJetsToLNu_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Winter23Reco-TRKRealistic_AlcaRecoRealisticTRK_126X_mcRun3_2022_realistic_postEE_v2-v2/GEN-SIM-RECO"
]

datasets["Background_ZPlusJets"] = [
    "/DYJetsToMuMu_M-50_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Winter23Reco-TRKRealistic_AlcaRecoRealisticTRK_preEE_126X_mcRun3_2022_realistic_v4-v2/GEN-SIM-RECO",
    "/DYJetsToMuMu_M-50_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Winter23Reco-TRKRealistic_AlcaRecoRealisticTRK_REAL_126X_mcRun3_2022_realistic_postEE_v2-v3/GEN-SIM-RECO"
]

datasets["Background_QCD"] = [
    "/QCD_PT-15to7000_TuneCP5_13p6TeV_pythia8/Run3Winter23Reco-FlatPU0to80_126X_mcRun3_2023_forPU65_v1-v2/GEN-SIM-RECO"
]

# central signal samples for HToSSTo4B, these are the low stats samples
datasets["Signal_RAW_HToSSTo4B_MH1000_MS450_CTau100000_PU60_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-100000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson60KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH1000_MS450_CTau100000_PU70_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-100000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v1/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH1000_MS450_CTau100000_Premix_2022preEE"] = [
    "/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-100000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH1000_MS450_CTau100000_Premix_2023postBPix"] = [
    "/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-100000mm_TuneCP5_13p6TeV-pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW",
    "/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-100000mm_TuneCP5_13p6TeV_pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH1000_MS450_CTau10000_PU60_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-10000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson60KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH1000_MS450_CTau10000_PU70_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-10000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v1/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH1000_MS450_CTau10000_Premix_2022preEE"] = [
    "/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-10000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH1000_MS450_CTau10000_Premix_2023postBPix"] = [
    "/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-10000mm_TuneCP5_13p6TeV-pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW",
    "/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-10000mm_TuneCP5_13p6TeV_pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH125_MS12_CTau9000_PU60_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-12_CTau-9000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson60KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH125_MS12_CTau9000_PU70_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-12_CTau-9000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v1/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH125_MS12_CTau9000_Premix_2022preEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-12_CTau-9000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH125_MS12_CTau9000_Premix_2023postBPix"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-12_CTau-9000mm_TuneCP5_13p6TeV-pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW",
    "/HTo2LongLivedTo4b_MH-125_MFF-12_CTau-9000mm_TuneCP5_13p6TeV_pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH125_MS12_CTau900_PU60_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-12_CTau-900mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson60KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH125_MS12_CTau900_PU70_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-12_CTau-900mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v1/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH125_MS12_CTau900_Premix_2022preEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-12_CTau-900mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH125_MS12_CTau900_Premix_2023postBPix"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-12_CTau-900mm_TuneCP5_13p6TeV-pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW",
    "/HTo2LongLivedTo4b_MH-125_MFF-12_CTau-900mm_TuneCP5_13p6TeV_pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH125_MS25_CTau15000_PU60_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-25_CTau-15000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson60KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH125_MS25_CTau15000_PU70_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-25_CTau-15000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v1/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH125_MS25_CTau15000_Premix_2022preEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-25_CTau-15000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH125_MS25_CTau15000_Premix_2023postBPix"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-25_CTau-15000mm_TuneCP5_13p6TeV-pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW",
    "/HTo2LongLivedTo4b_MH-125_MFF-25_CTau-15000mm_TuneCP5_13p6TeV_pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH125_MS25_CTau1500_PU60_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-25_CTau-1500mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson60KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH125_MS25_CTau1500_PU70_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-25_CTau-1500mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v1/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH125_MS25_CTau1500_Premix_2022preEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-25_CTau-1500mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH125_MS25_CTau1500_Premix_2023postBPix"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-25_CTau-1500mm_TuneCP5_13p6TeV-pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW",
    "/HTo2LongLivedTo4b_MH-125_MFF-25_CTau-1500mm_TuneCP5_13p6TeV_pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH125_MS50_CTau30000_PU60_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-50_CTau-30000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson60KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH125_MS50_CTau30000_PU70_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-50_CTau-30000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v1/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH125_MS50_CTau30000_Premix_2022preEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-50_CTau-30000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH125_MS50_CTau30000_Premix_2023postBPix"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-50_CTau-30000mm_TuneCP5_13p6TeV-pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW",
    "/HTo2LongLivedTo4b_MH-125_MFF-50_CTau-30000mm_TuneCP5_13p6TeV_pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH125_MS50_CTau3000_PU60_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-50_CTau-3000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson60KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH125_MS50_CTau3000_PU70_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-50_CTau-3000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v1/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH125_MS50_CTau3000_Premix_2022preEE"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-50_CTau-3000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH125_MS50_CTau3000_Premix_2023postBPix"] = [
    "/HTo2LongLivedTo4b_MH-125_MFF-50_CTau-3000mm_TuneCP5_13p6TeV-pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW",
    "/HTo2LongLivedTo4b_MH-125_MFF-50_CTau-3000mm_TuneCP5_13p6TeV_pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH250_MS120_CTau10000_PU60_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-10000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson60KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH250_MS120_CTau10000_PU70_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-10000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v1/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH250_MS120_CTau10000_Premix_2022preEE"] = [
    "/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-10000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH250_MS120_CTau10000_Premix_2023postBPix"] = [
    "/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-10000mm_TuneCP5_13p6TeV-pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW",
    "/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-10000mm_TuneCP5_13p6TeV_pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH250_MS120_CTau1000_PU60_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson60KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH250_MS120_CTau1000_PU70_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v1/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH250_MS120_CTau1000_Premix_2022preEE"] = [
    "/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH250_MS120_CTau1000_Premix_2023postBPix"] = [
    "/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_13p6TeV-pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW",
    "/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_13p6TeV_pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH250_MS120_CTau500_PU60_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-500mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson60KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH250_MS120_CTau500_PU70_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-500mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v1/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH250_MS120_CTau500_Premix_2022preEE"] = [
    "/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-500mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH250_MS120_CTau500_Premix_2023postBPix"] = [
    "/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-500mm_TuneCP5_13p6TeV-pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW",
    "/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-500mm_TuneCP5_13p6TeV_pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH250_MS60_CTau10000_PU60_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-250_MFF-60_CTau-10000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson60KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH250_MS60_CTau10000_PU70_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-250_MFF-60_CTau-10000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v1/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH250_MS60_CTau10000_Premix_2022preEE"] = [
    "/HTo2LongLivedTo4b_MH-250_MFF-60_CTau-10000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH250_MS60_CTau10000_Premix_2023postBPix"] = [
    "/HTo2LongLivedTo4b_MH-250_MFF-60_CTau-10000mm_TuneCP5_13p6TeV-pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW",
    "/HTo2LongLivedTo4b_MH-250_MFF-60_CTau-10000mm_TuneCP5_13p6TeV_pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH250_MS60_CTau1000_PU60_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-250_MFF-60_CTau-1000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson60KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH250_MS60_CTau1000_PU70_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-250_MFF-60_CTau-1000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v1/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH250_MS60_CTau1000_Premix_2022preEE"] = [
    "/HTo2LongLivedTo4b_MH-250_MFF-60_CTau-1000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH250_MS60_CTau1000_Premix_2023postBPix"] = [
    "/HTo2LongLivedTo4b_MH-250_MFF-60_CTau-1000mm_TuneCP5_13p6TeV-pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW",
    "/HTo2LongLivedTo4b_MH-250_MFF-60_CTau-1000mm_TuneCP5_13p6TeV_pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH250_MS60_CTau500_PU60_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-250_MFF-60_CTau-500mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson60KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH250_MS60_CTau500_PU70_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-250_MFF-60_CTau-500mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v1/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH250_MS60_CTau500_Premix_2022preEE"] = [
    "/HTo2LongLivedTo4b_MH-250_MFF-60_CTau-500mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH250_MS60_CTau500_Premix_2023postBPix"] = [
    "/HTo2LongLivedTo4b_MH-250_MFF-60_CTau-500mm_TuneCP5_13p6TeV-pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW",
    "/HTo2LongLivedTo4b_MH-250_MFF-60_CTau-500mm_TuneCP5_13p6TeV_pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH350_MS160_CTau10000_PU60_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-10000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson60KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH350_MS160_CTau10000_PU70_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-10000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v1/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH350_MS160_CTau10000_Premix_2022preEE"] = [
    "/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-10000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH350_MS160_CTau10000_Premix_2023postBPix"] = [
    "/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-10000mm_TuneCP5_13p6TeV-pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW",
    "/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-10000mm_TuneCP5_13p6TeV_pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH350_MS160_CTau1000_PU60_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-1000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson60KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH350_MS160_CTau1000_PU70_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-1000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v1/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH350_MS160_CTau1000_Premix_2022preEE"] = [
    "/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-1000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH350_MS160_CTau1000_Premix_2023postBPix"] = [
    "/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-1000mm_TuneCP5_13p6TeV-pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW",
    "/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-1000mm_TuneCP5_13p6TeV_pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH350_MS160_CTau500_PU60_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson60KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH350_MS160_CTau500_PU70_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v1/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH350_MS160_CTau500_Premix_2022preEE"] = [
    "/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH350_MS160_CTau500_Premix_2023postBPix"] = [
    "/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_13p6TeV-pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW",
    "/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_13p6TeV_pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH350_MS80_CTau10000_PU60_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-350_MFF-80_CTau-10000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson60KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH350_MS80_CTau10000_PU70_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-350_MFF-80_CTau-10000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v1/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH350_MS80_CTau10000_Premix_2022preEE"] = [
    "/HTo2LongLivedTo4b_MH-350_MFF-80_CTau-10000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH350_MS80_CTau10000_Premix_2023postBPix"] = [
    "/HTo2LongLivedTo4b_MH-350_MFF-80_CTau-10000mm_TuneCP5_13p6TeV-pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW",
    "/HTo2LongLivedTo4b_MH-350_MFF-80_CTau-10000mm_TuneCP5_13p6TeV_pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH350_MS80_CTau1000_PU60_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-350_MFF-80_CTau-1000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson60KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH350_MS80_CTau1000_PU70_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-350_MFF-80_CTau-1000mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v1/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH350_MS80_CTau1000_Premix_2022preEE"] = [
    "/HTo2LongLivedTo4b_MH-350_MFF-80_CTau-1000mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH350_MS80_CTau1000_Premix_2023postBPix"] = [
    "/HTo2LongLivedTo4b_MH-350_MFF-80_CTau-1000mm_TuneCP5_13p6TeV-pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW",
    "/HTo2LongLivedTo4b_MH-350_MFF-80_CTau-1000mm_TuneCP5_13p6TeV_pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH350_MS80_CTau500_PU60_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-350_MFF-80_CTau-500mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson60KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH350_MS80_CTau500_PU70_2022postEE"] = [
    "/HTo2LongLivedTo4b_MH-350_MFF-80_CTau-500mm_TuneCP5_13p6TeV-pythia8/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v1/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH350_MS80_CTau500_Premix_2022preEE"] = [
    "/HTo2LongLivedTo4b_MH-350_MFF-80_CTau-500mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/GEN-SIM-RAW"
]

datasets["Signal_RAW_HToSSTo4B_MH350_MS80_CTau500_Premix_2023postBPix"] = [
    "/HTo2LongLivedTo4b_MH-350_MFF-80_CTau-500mm_TuneCP5_13p6TeV-pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW",
    "/HTo2LongLivedTo4b_MH-350_MFF-80_CTau-500mm_TuneCP5_13p6TeV_pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/GEN-SIM-RAW"
]

