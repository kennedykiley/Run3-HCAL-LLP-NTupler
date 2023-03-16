cmsRun DisplacedHcalJetNTuplizer.py isData=True isSignal=False processEvents=1000 inputFiles=InputDataTest.txt debug=False outputFile=ntuple_output_test_data1.root
cmsRun DisplacedHcalJetNTuplizer.py isData=False isSignal=True processEvents=10000 inputFiles=InputSignalFiles.txt debug=False outputFile=ntuple_output_test_signal1.root
