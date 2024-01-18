# LLP Ntupler
Long-Lived Particle Ntupler based on AOD, adapted for use with HBHE rechits for Run 3 LLP analysis.

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

#### OR use shortcuts ####
proxy
crab_setup
```
Compile: 
```
scram b -j 8
```
If there is an error of the form `edmWriteConfigs: error while loading shared libraries: libssl.so.10: cannot open shared object file: No such file or directory`, try moving to lxplus8 and recompiling with `scram b clean; scram b -j 8`.

### Location of data and MC
Ntupler runs on data and MC, specifying which as an argument to `cmsRun`. The HCAL LLP skim is EXOLLPJetHCAL, in the dataset DisplacedJet. This can be found on [DAS](https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset%3D%2FDisplacedJet%2FRun2023*EXOLLPJetHCAL*%2FAOD) with the query `dataset=/DisplacedJet/Run2023*EXOLLPJetHCAL*/AOD`. The input file list can be made by running 
```
dasgoclient --query="file dataset=/DisplacedJet/Run2023C-EXOLLPJetHCAL-PromptReco-v4/AOD" > InputData_Run2023C-EXOLLPJetHCAL-PromptReco-v4.txt
```
Can also add specifications to ensure files are on DISK, such as `site=T1_US_FNAL_Disk`. CRAB job submission will cause a TAPE recall if the entire dataset is on TAPE, otherwise need to create a Rucio rule.

The H->XX->4b MC for 2022 are on DAS for [125 GeV](https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset%3D%2FggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV%2Flpclonglived-crab_PrivateProduction_Summer22_DR_step2_RECOSIM_ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_batch1_v1-59a22edf0600a784f6c900595d24e883%2FUSER+instance%3Dprod%2Fphys03) and [350 GeV](https://cmsweb.cern.ch/das/request?input=dataset%3D%2FHToSSTo4B_MH350_MS80_CTau500%2Flpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH350_MS80_CTau500_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76%2FUSER&instance=prod/phys03):
```
dataset=/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV/lpclonglived-crab_PrivateProduction_Summer22_DR_step2_RECOSIM_ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_batch1_v1-59a22edf0600a784f6c900595d24e883/USER instance=prod/phys03

dataset=/HToSSTo4B_MH350_MS80_CTau500/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH350_MS80_CTau500_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER/HToSSTo4B_MH350_MS80_CTau500/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH350_MS80_CTau500_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER
```
On command line, can be found via:
```
dasgoclient --limit=100 --query="file dataset=/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV/lpclonglived-crab_PrivateProduction_Summer22_DR_step2_RECOSIM_ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_batch1_v1-59a22edf0600a784f6c900595d24e883/USER instance=prod/phys03 | grep file.name"
```
### Running the ntupler
Test by running small numbers of events for testing:
```
cd run
# HCAL LLP skim
cmsRun ../python/DisplacedHcalJetNTuplizer.py isData=True isSignal=False processEvents=200 inputFiles=InputData_Run2023C-EXOLLPJetHCAL-PromptReco-v4.txt debug=False outputFile=ntuple_output_data_Run2023C-EXOLLPJetHCAL-PromptReco-v4.root

# MC signal 
cmsRun ../python/DisplacedHcalJetNTuplizer.py isData=False isSignal=True processEvents=200 inputFiles=InputSignalFilesTest.txt debug=False outputFile=ntuple_output_test_signal1.root
```
If there are issues with file access (''failed to open fiele at URL''), check on DAS if the data has been moved to TAPE. A temp workaround before they are accessible is to only use ones on disk.

Use CRAB to submit jobs (remember to change the Data / Signal flags in `DisplacedHcalJetNTuplizer.py`!). Note that the `number` in the CRAB python script must be changed for each job to run over all the datasets input. 
```
cmsenv
cd python
# note that for crab jobs (MC and data), the variables "signal" and "data" must be set by hand now in python/DisplacedHcalJetNTuplizer.py
crab submit -c crab_DisplacedHcalJetNTuplizer_MC_cfg.py 
crab submit -c crab_DisplacedHcalJetNTuplizer_QCD_cfg.py 
crab submit -c crab_DisplacedHcalJetNTuplizer_cfg.py
```
Useful CRAB commands:
```
crab submit -c <crab_cfg.py file> --dryrun
crab status -d <crab_directory>/<crab_project> --long
crab checkwrite --site T2_US_Wisconsin
```
Output ntuples are in `/hdfs/store/user/gkopp/ggH_HToSSTobbbb_MH*`.

Check event content with
```
edmDumpEventContent root://cmsxrootd.fnal.gov/</store/path/to/file.root> > EDM_content.txt
```

### Variables to check before running ntupler:
* isData: False if running signals, True if running data
* isSignal: True if running signals, False if running data
  * in particular, isSignal and isData need to be set in the file for a crab submission!
  * TODO: determine how to pass via config
* readGenVertexTime: False for LLP samples

## Ntuple use
After ntuples are made, they are used in the LLP_NuplerAnalyzer, from [here](https://github.com/gk199/Run3-HCAL-LLP-Analysis/tree/main)

## Location 

lxplus location (Gillian):
```
/afs/cern.ch/work/g/gkopp/2022_LLP_analysis/CMSSW_13_1_0/src/cms_lpc_llp/Run3-HCAL-LLP-NTupler
```

## Archive

Many earlier files have been moved to `python/Archive`. 

The High MET skim we start with from 2022 data are here on [DAS](https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset%3D%2FJetMET%2FRun2022G-EXOHighMET-PromptReco-v1%2FRAW-RECO). 

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