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
/store/user/gkopp/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV/LLP_MC_test__20231017_175630/231017_155638/0000/*.root			# LLP MC, MH125, MS15, cTau1000
```
About 500k events in LLP MC and LLP skim, and about 250k events in QCD MC. 