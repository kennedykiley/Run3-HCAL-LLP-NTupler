# llp_ntupler
Long Lived Particle Ntupler based on AOD 


### Setup CMSSW & clone the ntuples:
```bash
# Done once to setup environment
cmsrel CMSSW_9_4_4
cd src
git clone git@github.com:RazorCMS/SUSYBSMAnalysis-JetNtupler.git SUSYBSMAnalysis/JetNtupler
scram b
cmsenv
# Run the ntuples
cmsRun python/jetNtupler_MC_AOD.py
```

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