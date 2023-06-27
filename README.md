# LLP Ntupler
Long Lived Particle Ntupler based on AOD, adapted for use with HBHE rechits for Run 3 LLP analysis.

# Setup Ntupler
```
cmsrel CMSSW_<>
mkdir cms_lpc_llp
cd cms_lpc_llp
git clone git@github.com:kennedykiley/Run3-HCAL-LLP-NTupler.git
cd Run3-HCAL-LLP-NTupler
git checkout -b <your-branch>
```

# Run the ntuples
Setup grid proxy
```
voms-proxy-init --rfc --voms cms --valid 48:00
cp /tmp/x509up_u101898 /afs/cern.ch/user/g/gkopp
chmod 777 /afs/cern.ch/user/g/gkopp/x509up_u101898
source /cvmfs/cms.cern.ch/common/crab-setup.sh
# OR
proxy
crab_setup

scram b -j 8
```
Ntupler runs on data and MC, specifying which as an argument to `cmsRun`. Many earlier files have been moved to `python/Archive`. 

The HCAL LLP skim is EXOLLPJetHCAL, in the dataset DisplacedJet. This can be found on [DAS](https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset%3D%2FDisplacedJet%2FRun2023*EXOLLPJetHCAL*%2FAOD) with the query `dataset=/DisplacedJet/Run2023*EXOLLPJetHCAL*/AOD`. The input file list can be made by running 
```
dasgoclient --query="file dataset=/DisplacedJet/Run2023C-EXOLLPJetHCAL-PromptReco-v4/AOD" > InputData_Run2023C-EXOLLPJetHCAL-PromptReco-v4.txt
```
Can also add specifications to ensure files are on DISK, such as `site=T1_US_FNAL_Disk`. Check event content with
```
edmDumpEventContent root://cmsxrootd.fnal.gov/</store/path/to/file.root> > EDM_content.txt
```

```
cd run
# High MET skim
cmsRun ../python/DisplacedHcalJetNTuplizer.py isData=True isSignal=False processEvents=200 inputFiles=InputDataMETSkimTest.txt debug=False outputFile=ntuple_output_test_data_METSkimtest.root

# HCAL LLP skim
cmsRun ../python/DisplacedHcalJetNTuplizer.py isData=True isSignal=False processEvents=200 inputFiles=InputData_Run2023C-EXOLLPJetHCAL-PromptReco-v4.txt debug=False outputFile=ntuple_output_data_Run2023C-EXOLLPJetHCAL-PromptReco-v4.root

# MC signal 
cmsRun ../python/DisplacedHcalJetNTuplizer.py isData=False isSignal=True processEvents=200 inputFiles=InputSignalFilesTest.txt debug=False outputFile=ntuple_output_test_signal1.root
```

Running with CRAB:
```
cmsenv
cd python
crab submit -c crab_DisplacedHcalJetNTuplizer_MC_cfg.py 
# note that for MC, the variables "signal" and "data" must be set by hand now in python/DisplacedHcalJetNTuplizer.py
crab submit -c crab_DisplacedHcalJetNTuplizer_cfg.py

# Useful commands
crab submit -c <crab_cfg.py file> --dryrun
crab status -d <crab_directory>/<crab_project> --long

crab checkwrite --site T2_US_Wisconsin
```

Output is in `/hdfs/store/user/gkopp/ggH_HToSSTobbbb_MH*`.


The High MET skim we start with from 2022 data are here on [DAS](https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset%3D%2FJetMET%2FRun2022G-EXOHighMET-PromptReco-v1%2FRAW-RECO). 

The H->XX->4b MC for 2022 are here:
```
dasgoclient --limit=100 --query="file dataset=/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV/lpclonglived-crab_PrivateProduction_Summer22_DR_step2_RECOSIM_ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_batch1_v1-59a22edf0600a784f6c900595d24e883/USER instance=prod/phys03 | grep file.name"
```

### Variables to check before running nTupler:
* isData: False if running signals, True if running data
* isSignal: True if running signals, False if running data
  ** in particular, isSignal and isData need to be set in the file for a crab submission! (until determine how to pass via config)
* readGenVertexTime: False for LLP samples

## Ntuple use
After ntuples are made, they are used in the LLP_NuplerAnalyzer, from [here](https://github.com/gk199/Run3-HCAL-LLP-Analysis/tree/main)

## Location 

lxplus location (Gillian): `/afs/cern.ch/work/g/gkopp/2022_LLP_analysis/CMSSW_12_4_6/src/cms_lpc_llp/Run3-HCAL-LLP-NTupler'

## Archive

### Initial setup of CMSSW & clone the ntupler from CMS LPC area:
```bash
# Done once to setup environment
cmsrel CMSSW_12_4_6
cd CMSSW_12_4_6/src
git clone -b run3 git@github.com:cms-lpc-llp/llp_ntupler.git cms_lpc_llp/llp_ntupler
rm cms_lpc_llp/llp_ntupler/plugins/displacedJetMuon.*
rm cms_lpc_llp/llp_ntupler/plugins/displacedJetTiming_ntupler.*
rm cms_lpc_llp/llp_ntupler/plugins/displacedJetTiming_aux.cc
rm cms_lpc_llp/llp_ntupler/plugins/llp_ntupler*

# Fix python print statements in cms_lpc_llp/llp_ntupler/python/crab_scripts/multi_crab_ntuples.py and cms_lpc_llp/llp_ntupler/python/crab_scripts/multi_crab_resubmit.py to make compatible with python 3

# crab scripts are giving errors, move them so they are still avaliable for reference:
mv cms_lpc_llp/llp_ntupler/python/crab_scripts/multi_crab_resubmit.py cms_lpc_llp/llp_ntupler/python/crab_scripts/multi_crab_resubmit.py.old
mv cms_lpc_llp/llp_ntupler/python/crab_scripts/multi_crab_ntuples.py cms_lpc_llp/llp_ntupler/python/crab_scripts/multi_crab_ntuples.py.old

scram b -j 8
cmsenv

cd cms_lpc_llp/llp_ntupler
```