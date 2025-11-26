# LLP Ntupler
Long-Lived Particle Ntupler based on AOD, adapted for use with HBHE rechits for Run 3 LLP analysis.

# Setup Ntupler Code
```
cmsrel <CMSSW version> # Use CMSSW_13_2_0 for NTuples v4
mkdir cms_lpc_llp
cd cms_lpc_llp
git clone git@github.com:kennedykiley/Run3-HCAL-LLP-NTupler.git
cd Run3-HCAL-LLP-NTupler
git checkout -b <your-branch>
```

## Get JEC and JER files
[JER twiki](https://cms-jerc.web.cern.ch/Recommendations/#jet-energy-resolution): Get textfiles from the JRDatabase from github and put them in `Run3-HCAL-LLP-NTupler/data/JEC_JER/JRDatabase/textFiles/`. 

[JEC twiki](https://cms-jerc.web.cern.ch/Recommendations/#jet-energy-scale): Get the database files from the JECDatabase github, and put them in `Run3-HCAL-LLP-NTupler/data/JEC_JER/JECDatabase/SQLiteFiles/`.

At the time of v5 ntuple preparation, v3 is the recommended JES version and v1 is the recommended JER version. 

Make sure the `mapping` in L324 of `DisplacedHcalJetsTuplizer.py` is correct for all the data and MC processed, otherwise the JEC and JER tags will not be found (this will cause a runtime error). 

The saved branches are:
- `jetRaw`: no JEC applied, this is the uncorrected jet
- `jet_*_noJER`: no smearing applied (MC only)
- `jet_*_JER_up/down`: smeared up / down variation with JER
- `jet_E`, `jet_Pt`: corrected and (MC only) smeared jet. This should be used for analysis

# Run the ntuples
Setup grid proxy
```
voms-proxy-init --rfc --voms cms --valid 48:00
cp /tmp/x509up_u101898 /afs/cern.ch/user/<initial>/<name>
chmod 777 /afs/cern.ch/user/<initial>/<name>/x509up_u101898
source /cvmfs/cms.cern.ch/common/crab-setup.sh
```
Or use the bash shortcuts (if setup):
```
proxy
crab_setup
```
Compile: 
```
cmsenv
scram b -j 8
```
If there is an error of the form `edmWriteConfigs: error while loading shared libraries: libssl.so.10: cannot open shared object file: No such file or directory`, try moving to lxplus8 and recompiling with `scram b clean; scram b -j 8`.

## Running the ntupler
Test by running small numbers of events for testing:
```
cd run
# HCAL LLP skim
cmsRun ../python/DisplacedHcalJetNTuplizer.py isData=True isSignal=False processEvents=200 inputFiles=InputData_Run2023C-EXOLLPJetHCAL-PromptReco-v4.txt debug=False outputFile=ntuple_output_data_Run2023C-EXOLLPJetHCAL-PromptReco-v4.root

# MC signal 
cmsRun ../python/DisplacedHcalJetNTuplizer.py isData=False isSignal=True processEvents=200 inputFiles=InputSignalFilesTest.txt debug=False outputFile=ntuple_output_test_signal1.root
```
If there are issues with file access (''failed to open file at URL''), check on DAS if the data has been moved to TAPE. A temp workaround before they are accessible is to only use ones on disk.

Use CRAB to submit jobs (remember to change the Data / Signal flags in `DisplacedHcalJetNTuplizer.py`!). Note that the `number` in the CRAB python script must be changed for each job to run over all the datasets input. 
```
cmsenv
cd python
# note that for crab jobs (MC and data), the variables "signal" and "data" must be set by hand now in python/DisplacedHcalJetNTuplizer.py
# note that as of February 2024, the crab submissions ONLY work on lxplus 7, with crab-pre submit
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

### CRAB Wrapper 1 for Automating Submissions
First check that the dataset is on disk:
```
python3 checkDatasetAvailability.py <txt file of datasets to check>
```
If so, proceed with the crab submission. AOD-tier will be automaticall recalled, but RAW-RECO will not. 
```
proxy
crab_setup
python3 CrabSubmitWrapper.py
```
This handles changing the variables (isData and isSignal) in DisplacedHcalJetNTuplizer.py, as well as doing one crab submission per dataset listed. 

### CRAB Wrapper 2 for Automating Submissions

From src/
```
cmsenv
voms-proxy-init -voms cms
source /cvmfs/cms.cern.ch/crab3/crab.sh
# scram b -j 10 # if needed
cd cms_lpc_llp/Run3-HCAL-LLP-NTupler/python
```
Edit `makeBulkCrabSubmission.py`:
* Update line 15 to point to your crab output directory (should be outside of your CMSSW area). Make this directory if it does not exist

Run: 
```
# Setup crab jobs (setup to signal by default, but need to edit line 147 of `makeBulkCrabSubmission.py` to run over other datasets):
python3 makeBulkCrabSubmission.py Signal

# Submit crab jobs
source Signal/submit.sh
```

### Variables to check before running ntupler:
* isData: False if running signals, True if running data
* isSignal: True if running signals, False if running data
  * in particular, isSignal and isData need to be set in the file for a crab submission!
  * TODO: determine how to pass via config
* readGenVertexTime: False for LLP samples

### Location of data and MC
Ntupler runs on data and MC, specifying which as an argument to `cmsRun`. The HCAL LLP skim is EXOLLPJetHCAL, in the dataset DisplacedJet. This can be found on [DAS](https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset%3D%2FDisplacedJet%2FRun2023*EXOLLPJetHCAL*%2FAOD) with the query `dataset=/DisplacedJet/Run2023*EXOLLPJetHCAL*/AOD`. The input file list can be made by running 
```
dasgoclient --query="file dataset=/DisplacedJet/Run2023C-EXOLLPJetHCAL-PromptReco-v4/AOD" > InputData_Run2023C-EXOLLPJetHCAL-PromptReco-v4.txt
```
Can also add specifications to ensure files are on DISK, such as `site=T1_US_FNAL_Disk`. CRAB job submission will cause a TAPE recall if the entire dataset is on TAPE, otherwise need to create a Rucio rule.

Data:
```
dataset=/DisplacedJet/Run2023*-EXOLLPJetHCAL-PromptReco*/AOD
dataset=/JetMET*/Run*EXOHighMET-PromptReco*/RAW-RECO
```
which is the displaced jet skim (implemented for 2023 onwards) and the EXO high MET skim, for the background estimation. 

The H->XX->4b MC for 2022 are on DAS for [125 GeV, mX = 15](https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset%3D%2FggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV%2Flpclonglived-crab_PrivateProduction_Summer22_DR_step2_RECOSIM_ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_batch1_v1-59a22edf0600a784f6c900595d24e883%2FUSER+instance%3Dprod%2Fphys03) (2M total events) and [350 GeV, mX = 80](https://cmsweb.cern.ch/das/request?input=dataset%3D%2FHToSSTo4B_MH350_MS80_CTau500%2Flpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH350_MS80_CTau500_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76%2FUSER&instance=prod/phys03) (0.5M total events)
```
dataset=/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV/lpclonglived-crab_PrivateProduction_Summer22_DR_step2_RECOSIM_ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_batch1_v1-59a22edf0600a784f6c900595d24e883/USER instance=prod/phys03

dataset=/HToSSTo4B_MH350_MS80_CTau500/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH350_MS80_CTau500_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER/
```
The samples produced in 2024 have higher LLP masses, with 125 GeV (mX = 50), 250 GeV (mX = 120), and 350 GeV (mX = 160). Batch 1 has 1M events and 10k files per sample (100 events per file). Batch 2 has higher stats, with 3M events total and 10k files per sample (300 events per sample). In the crab config, set: `units per job = 50`, `n jobs = 200`, which will work for both batch 1 and batch 2.
```
dataset=/HToSSTo4B_MH125_MS50_CTau3000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH125_MS50_CTau3000_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER
dataset=/HToSSTo4B_MH125_MS50_CTau3000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH125_MS50_CTau3000_batch2_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER

dataset=/HToSSTo4B_MH250_MS120_CTau10000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH250_MS120_CTau10000_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER
dataset=/HToSSTo4B_MH250_MS120_CTau10000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH250_MS120_CTau10000_batch2_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER

dataset=/HToSSTo4B_MH350_MS160_CTau10000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH350_MS160_CTau10000_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER
dataset=/HToSSTo4B_MH350_MS160_CTau10000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH350_MS160_CTau10000_batch2_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER
```
On command line, can be found via:
```
dasgoclient --limit=100 --query="file dataset=/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV/lpclonglived-crab_PrivateProduction_Summer22_DR_step2_RECOSIM_ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_batch1_v1-59a22edf0600a784f6c900595d24e883/USER instance=prod/phys03" >> output_file_name.txt
```

## Ntuple use
After ntuples are made, they are used in the LLP_NuplerAnalyzer, from [here](https://github.com/gk199/Run3-HCAL-LLP-Analysis/tree/main)

## Location 

lxplus location (Gillian):
```
/afs/cern.ch/work/g/gkopp/2022_LLP_analysis/CMSSW_13_2_0/src/cms_lpc_llp/Run3-HCAL-LLP-NTupler
```

# Archive

Many earlier files have been moved to `python/Archive`. 

The High MET skim we start with from 2022 data are here on [DAS](https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset%3D%2FJetMET%2FRun2022G-EXOHighMET-PromptReco-v1%2FRAW-RECO). 

## Initial setup of CMSSW & clone the ntupler from CMS LPC area:
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