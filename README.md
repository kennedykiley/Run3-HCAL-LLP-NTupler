# llp_ntupler
Long Lived Particle Ntupler based on AOD 

lxplus location: `/afs/cern.ch/work/g/gkopp/2022_LLP_analysis/CMSSW_12_4_6/src/cms_lpc_llp/llp_ntupler`
Moved to: `/afs/cern.ch/work/g/gkopp/2022_LLP_analysis/CMSSW_12_4_6/src/cms_lpc_llp/Run3-HCAL-LLP-NTupler'

### Setup CMSSW & clone the ntuples:
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

# Run the ntuples
cd cms_lpc_llp/llp_ntupler

voms-proxy-init --rfc --voms cms --valid 48:00
cp /tmp/x509up_u101898 /afs/cern.ch/user/g/gkopp
chmod 777 /afs/cern.ch/user/g/gkopp/x509up_u101898
source /cvmfs/cms.cern.ch/common/crab-setup.sh
# OR
proxy
crab_setup

scram b -j 8

# below are files directly adapted from muon ntupler
cmsRun python/displacedJetMuon_ntupler_Data_2022_MuonShowerSkim.py
cmsRun python/displacedJetMuon_ntupler_Data_2022_MuonShowerSkim_small.py # this is for Run 3 data
cmsRun python/prod.py # this is for Run 3 MC

# moving to HCAL jets more specific files. Data input is working now, MC gives error (still troubleshooting)
cmsRun python/DisplacedHcalJetNTuplizer.py isData=1 inputFiles=2022Data.txt processEvents=500
cmsRun python/DisplacedHcalJetNTuplizer.py inputFiles=2022MC.txt processEvents=500

# full updates with working ntupler for data and MC! Many of the above files have now been moved to python/Archive 
cd run
cmsRun ../python/DisplacedHcalJetNTuplizer.py isData=True isSignal=False processEvents=1000 inputFiles=InputDataTest.txt debug=False outputFile=ntuple_output_test_data1.root
cmsRun ../python/DisplacedHcalJetNTuplizer.py isData=True isSignal=False processEvents=20000 inputFiles=InputDataMETSkimTest.txt debug=False outputFile=ntuple_output_test_data_METSkimtest.root

cmsRun ../python/DisplacedHcalJetNTuplizer.py isData=False isSignal=True processEvents=1000 inputFiles=InputSignalFilesTest.txt debug=False outputFile=ntuple_output_test_signal1.root
```

The files from the HMT dataset in 2022 are on Caltech T2: `/DisplacedJet/Run2022E-EXOCSCCluster-PromptReco-v1/USER`. This file is used for testing, it has 5k events: `/store/data/Run2022E/DisplacedJet/USER/EXOCSCCluster-PromptReco-v1/000/360/017/00000/eae65e97-9f58-4119-9806-a3226ecba729.root`. 

### Variables to check before running nTupler:
* isQCD: False if running signals, True if running QCD
  * Matching is done differently for QCD
* isFourJet: False for glueball model; True for fourjet
  * Difference is LLP particle IDs are different in the two models
* readGenVertexTime: False for glueball model; True for fourjet
  * Difference is glueball model used 2016 condition, doesn't have genVertexTime info for now. (will have them when the new samples are ready)
* Output file name:
```bash
process.TFileService = cms.Service("TFileService",
        fileName = cms.string('input.root'),
    closeFileFast = cms.untracked.bool(True)
)
```
* Input file name:
```bash
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('output.root'),
)
```

### AODSIM Location

They are all on Caltech tier2: (should be accessible through xrootd)
```/mnt/hadoop/store/user/christiw/RunII2016/```
The different directories are different h production mode
* ggh: ppTohToSS1SS2_SS1Tobb_SS2Toveve_MC_prod
* vbfh: ppTohjjToSS1SS2_SS1Tobb_SS2Toveve_MC_prod
* wh: ppTohwToSS1SS2_SS1Tobb_SS2Toveve_MC_prod
* zh: ppTohwToSS1SS2_SS1Tobb_SS2Toveve_MC_prod
* ggh with ISR: ppTohToSS1SS2_SS1Tobb_SS2Toveve-ppTojhToSS1SS2_SS1Tobb_SS2Toveve_MC_prod
* There are a few mass points / ctau points for each production modes, some production modes only have one point
* Once in a particular mass/ctau directory, AODSIM is in the directory that ends with *DR-AODSIM_CaltechT2/
* eg: the four ROOT files in ```/mnt/hadoop/store/user/christiw/RunII2016/ppTohToSS1SS2_SS1Tobb_SS2Toveve_MC_prod/ppTohToSS1SS2_SS1Tobb_SS2Toveve_run_m10_pl10000_ev10000/crab_CMSSW_8_0_31_ppTohToSS1SS2_SS1Tobb_SS2Toveve_run_m10_pl10000_ev10000_DR-AODSIM_CaltechT2/19012\
0_072347/0000/```
are AODSIM
