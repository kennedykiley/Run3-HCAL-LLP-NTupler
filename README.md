# llp_ntupler
Long Lived Particle Ntupler based on AOD 


### Setup CMSSW & clone the ntuples:
```bash
# Done once to setup environment
cmsrel CMSSW_12_4_6
cd CMSSW_12_4_6/src
git clone -b run3 git@github.com:cms-lpc-llp/llp_ntupler.git cms_lpc_llp/llp_ntupler
rm cms_lpc_llp/llp_ntupler/plugins/displacedJetMuon_dump.*
rm cms_lpc_llp/llp_ntupler/plugins/displacedJetMuon_rechit_studies.*
rm cms_lpc_llp/llp_ntupler/plugins/displacedJetTiming_ntupler.*
rm cms_lpc_llp/llp_ntupler/plugins/displacedJetTiming_aux.cc 
rm cms_lpc_llp/llp_ntupler/plugins/llp_ntupler*  

# Fix python print statements in cms_lpc_llp/llp_ntupler/python/crab_scripts/multi_crab_ntuples.py and cms_lpc_llp/llp_ntupler/python/crab_scripts/multi_crab_resubmit.py to make compatible with python 3

# crab scripts are giving errors, move them so they are still avaliable for reference: 
mv cms_lpc_llp/llp_ntupler/python/crab_scripts/multi_crab_resubmit.py cms_lpc_llp/llp_ntupler/python/crab_scripts/multi_crab_resubmit.py.old 
mv cms_lpc_llp/llp_ntupler/python/crab_scripts/multi_crab_ntuples.py cms_lpc_llp/llp_ntupler/python/crab_scripts/multi_crab_ntuples.py.old

scram b
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

cmsRun python/displacedJetMuon_ntupler_Data_2022_MuonShowerSkim.py
cmsRun python/displacedJetMuon_ntupler_Data_2022_MuonShowerSkim_small.py
#cmsRun python/jetNtupler_MC_AOD.py
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