# Test Data (AOD)
cmsRun ../python/DisplacedHcalJetNTuplizer.py isData=True isSignal=False processEvents=100 inputFiles=InputDataTest.txt debug=False outputFile=ntuple_output_test_data.root

# Test Data (Reco from RAW)
cmsRun ../python/DisplacedHcalJetNTuplizer.py isData=True isSignal=False processEvents=100 inputFiles=InputData_DisplacedJet_Run2023D-v1_RAW_369927.txt debug=False outputFile=ntuple_output_test_data_fromraw.root recoFromRAW=True

# Test Signal
cmsRun ../python/DisplacedHcalJetNTuplizer.py isData=False isSignal=True processEvents=100 inputFiles=InputSignal_Run3_ggH_HToSSTobbbb_MH-125_MS-50_CTau3000_13p6TeV.txt debug=False outputFile=ntuple_output_test_signal.root
