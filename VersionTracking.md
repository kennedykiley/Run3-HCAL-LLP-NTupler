# Ntuples by version:

## V1
Updates: 
* Thread friendly EDAnalyzer version
* L1 jets added
* Track information added
* HBHE rechits > 0.5 GeV
* HLT paths updated
* TDC decoding in ntupler

Location of ntuples:
```
/store/user/gkopp/DisplacedJet/Run2023C-EXOLLPJetHCAL-PromptReco-v4_AOD_20230712_143615/230712_123655/0000/*.root				# 2023 Era C data
/store/user/gkopp/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV/LLP_MC_test__20230713_102955/230713_083027/0000/*.root			# LLP MC, MH125, MS15, cTau1000
```

## V2
Updates:
* HB rechits over 0.1 GeV within dR < 0.5 of a PF jet saved
* `FillTriggerBranches` fixed -- HLT now has correct order for name association
* Unprescaled HLT added, `HLT_L1SingleLLPJet`
* Added PU information
* Python script for QCD MC files added

Location of ntuples:
```
/store/user/kikenned/QCD_PT-15to7000_TuneCP5_13p6TeV_pythia8/QCD_MC__20231005_210638/231005_190646/0000/*.root					# QCD MC
/store/user/gkopp/DisplacedJet/Run2023C-EXOLLPJetHCAL-PromptReco-v4_AOD_20231017_175654/231017_155703/0000/*.root				# 2023 Era C data
/store/user/gkopp/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV/LLP_MC_test__20231017_175630/231017_155638/0000/*.root			# LLP MC, MH125, MS15, cTau1000. 500k events processed of 2M available 
```
About 500k events in LLP MC and LLP skim, and about 250k events in QCD MC. 

## V3
Updates:
* Add muon and electron isolation variables (pileup, charged, photon, and neutral hadron isolation)
* Save eta of electron super cluster
* Save muon IP 3D significance (for muon ID requirements) -- note, this is saved, but in all test files is the default value (max of a double for the vars used to calculate this)
* Save fixed grid rho fast jet all (for electron isolation calculation)

Processed 2023 DisplacedJet skim data.

Location of ntuples:
```
/store/user/gkopp/DisplacedJet/Run2023B-EXOLLPJetHCAL-PromptReco-v1_AOD_20231107_180123/231107_170140/0000/*.root
/store/user/gkopp/DisplacedJet/Run2023C-EXOLLPJetHCAL-PromptReco-v1_AOD_20231108_105637/231108_095708/0000/*.root
/store/user/gkopp/DisplacedJet/Run2023C-EXOLLPJetHCAL-PromptReco-v2_AOD_20231108_110936/231108_101012/0000/*.root and .../0001/*.root
/store/user/gkopp/DisplacedJet/Run2023C-EXOLLPJetHCAL-PromptReco-v3_AOD_20231108_112347/231108_102434/0000/*.root and .../0001/*.root
/store/user/gkopp/DisplacedJet/Run2023C-EXOLLPJetHCAL-PromptReco-v4_AOD_20231114_143556/231114_133610/0000/*.root
/store/user/gkopp/DisplacedJet/Run2023D-EXOLLPJetHCAL-PromptReco-v1_AOD_20231114_143525/231114_133531/0000/*.root
/store/user/gkopp/DisplacedJet/Run2023D-EXOLLPJetHCAL-PromptReco-v2_AOD_20231114_143510/231114_133516/0000/*.root
```

* Removed `GEMRecoGeometryRcd` as it is not a collection present in the 350 GeV H->LLP samples and not needed for the analysis.

Location of 350 GeV ntuples:
```
/hdfs/store/user/gkopp/HToSSTo4B_MH350_MS80_CTau500/LLP_MC_350__20231129_104033/231129_094141/0000/*.root (500k events, all available MC processed)
```

Location of 125 GeV ntuples:
```
/hdfs/store/user/gkopp/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV/LLP_MC_125__20240119_174604/240119_164620/0000/*.root (500k events, 25\% of available MC)
```

## V4
Updates:
* MET filters added
* All displaced jet HLT paths added
* Prescale for HLT paths added