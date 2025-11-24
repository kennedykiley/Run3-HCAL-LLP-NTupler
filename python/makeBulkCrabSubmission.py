""""
Script to mass-produce crab files to avoid annoying editing

Usage:
1. Edit dataset names (see lists and "EDIT ME!" comment below)
2. Run: `python3 makeBuleCrabSubmission.py <output_directory>`
"""

import os, re, sys

pwd = os.getcwd()
crab_filepath = os.path.join(pwd, "../python/crab_DisplacedHcalJetNTuplizer_DO-NOT-EDIT_cfg.py")

# Edit me:
crab_output_dir = '/afs/cern.ch/work/g/gkopp/2022_LLP_analysis/CRAB_Workarea/NTuples_v4/'

datasets = {}


datasets["Data_DisplacedJet_Run2022"] = [
    #"/DisplacedJet/Run2022A-v1/RAW", #  900 Gev
    #"/DisplacedJet/Run2022B-v1/RAW", # Commissioning
    "/DisplacedJet/Run2022C-v1/RAW",
    "/DisplacedJet/Run2022D-v1/RAW",  # There are D-v2 and -v3, but they're not listed in DAS 
    "/DisplacedJet/Run2022E-v1/RAW",
    "/DisplacedJet/Run2022F-v1/RAW",
    "/DisplacedJet/Run2022G-v1/RAW",
]

datasets["Data_EXOLLPJetHCAL_Run2023"] = [
    "/DisplacedJet/Run2023B-EXOLLPJetHCAL-PromptReco-v1/AOD",
    "/DisplacedJet/Run2023C-EXOLLPJetHCAL-PromptReco-v1/AOD",
    "/DisplacedJet/Run2023C-EXOLLPJetHCAL-PromptReco-v2/AOD",
    "/DisplacedJet/Run2023C-EXOLLPJetHCAL-PromptReco-v3/AOD",
    "/DisplacedJet/Run2023C-EXOLLPJetHCAL-PromptReco-v4/AOD",
    "/DisplacedJet/Run2023D-EXOLLPJetHCAL-PromptReco-v1/AOD",
    "/DisplacedJet/Run2023D-EXOLLPJetHCAL-PromptReco-v2/AOD",

]

datasets["Data_EXOLLPJetHCAL_Run2024"] = [
    #"/DisplacedJet/Run2024A-EXOLLPJetHCAL-PromptReco-v1/AOD",
    "/DisplacedJet/Run2024B-EXOLLPJetHCAL-PromptReco-v1/AOD",
    "/DisplacedJet/Run2024C-EXOLLPJetHCAL-PromptReco-v1/AOD",
    "/DisplacedJet/Run2024D-EXOLLPJetHCAL-PromptReco-v1/AOD",
    "/DisplacedJet/Run2024E-EXOLLPJetHCAL-PromptReco-v1/AOD",
    "/DisplacedJet/Run2024E-EXOLLPJetHCAL-PromptReco-v2/AOD",
    "/DisplacedJet/Run2024F-EXOLLPJetHCAL-PromptReco-v1/AOD",
    "/DisplacedJet/Run2024G-EXOLLPJetHCAL-PromptReco-v1/AOD",
    "/DisplacedJet/Run2024H-EXOLLPJetHCAL-PromptReco-v1/AOD",
    "/DisplacedJet/Run2024I-EXOLLPJetHCAL-PromptReco-v1/AOD",
    "/DisplacedJet/Run2024I-EXOLLPJetHCAL-PromptReco-v2/AOD",
]

datasets["Data_ZMu_Run2022"] = [
    "/Muon/Run2022C-ZMu-10Dec2022-v1/RAW-RECO",
    "/Muon/Run2022C-ZMu-27Jun2023-v1/RAW-RECO",
    "/Muon/Run2022C-ZMu-PromptReco-v1/RAW-RECO",
    "/Muon/Run2022D-ZMu-10Dec2022-v1/RAW-RECO",
    "/Muon/Run2022D-ZMu-27Jun2023-v2/RAW-RECO",
    "/Muon/Run2022D-ZMu-PromptReco-v1/RAW-RECO",
    "/Muon/Run2022D-ZMu-PromptReco-v2/RAW-RECO",
    "/Muon/Run2022D-ZMu-PromptReco-v3/RAW-RECO",
    "/Muon/Run2022E-ZMu-10Dec2022-v2/RAW-RECO",
    "/Muon/Run2022E-ZMu-27Jun2023-v1/RAW-RECO",
    "/Muon/Run2022E-ZMu-PromptReco-v1/RAW-RECO",
    "/Muon/Run2022F-ZMu-PromptReco-v1/RAW-RECO",
    "/Muon/Run2022G-ZMu-PromptReco-v1/RAW-RECO",
]

datasets["Data_ZMu_Run2023"] = [
    '/Muon0/Run2023B-ZMu-PromptReco-v1/RAW-RECO',
    '/Muon0/Run2023C-ZMu-PromptReco-v1/RAW-RECO',
    '/Muon0/Run2023C-ZMu-PromptReco-v2/RAW-RECO',
    '/Muon0/Run2023C-ZMu-PromptReco-v3/RAW-RECO',
    '/Muon0/Run2023C-ZMu-PromptReco-v4/RAW-RECO',
    '/Muon0/Run2023D-ZMu-PromptReco-v1/RAW-RECO',
    '/Muon0/Run2023D-ZMu-PromptReco-v2/RAW-RECO',
    #'/Muon1/Run2023A-ZMu-PromptReco-v2/RAW-RECO',
    '/Muon1/Run2023B-ZMu-PromptReco-v1/RAW-RECO',
    '/Muon1/Run2023C-ZMu-PromptReco-v1/RAW-RECO',
    '/Muon1/Run2023C-ZMu-PromptReco-v2/RAW-RECO',
    '/Muon1/Run2023C-ZMu-PromptReco-v3/RAW-RECO',
    '/Muon1/Run2023C-ZMu-PromptReco-v4/RAW-RECO',
    '/Muon1/Run2023D-ZMu-PromptReco-v1/RAW-RECO',
    '/Muon1/Run2023D-ZMu-PromptReco-v2/RAW-RECO'
]

#datasets["Signal_HToSSTo4B_MH125_MS15_CTau1000"] = [
#    "/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV/lpclonglived-crab_PrivateProduction_Summer22_DR_step2_RECOSIM_ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_batch1_v1-59a22edf0600a784f6c900595d24e883/USER"
#]

datasets["Signal_HToSSTo4B_MH125_MS50_CTau3000"] = [
    "/HToSSTo4B_MH125_MS50_CTau3000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH125_MS50_CTau3000_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER",
    "/HToSSTo4B_MH125_MS50_CTau3000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH125_MS50_CTau3000_batch2_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER"
]

datasets["Signal_HToSSTo4B_MH250_MS120_CTau10000"] = [
    "/HToSSTo4B_MH250_MS120_CTau10000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH250_MS120_CTau10000_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER",
    "/HToSSTo4B_MH250_MS120_CTau10000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH250_MS120_CTau10000_batch2_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER"
]

datasets["Signal_HToSSTo4B_MH350_MS80_CTau500"] = [
    "/HToSSTo4B_MH350_MS80_CTau500/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH350_MS80_CTau500_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER"
]

datasets["Signal_HToSSTo4B_MH350_MS160_CTau10000"] = [
    "/HToSSTo4B_MH350_MS160_CTau10000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH350_MS160_CTau10000_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER",
    "/HToSSTo4B_MH350_MS160_CTau10000/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH350_MS160_CTau10000_batch2_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER"
]

datasets["Background_WPlusJets"] = [
    "/WJetsToLNu_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Winter23Reco-TRKRealistic_AlcaRecoRealisticTRK_preEE_126X_mcRun3_2022_realistic_v4-v2/GEN-SIM-RECO",
    "/WJetsToLNu_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Winter23Reco-TRKRealistic_AlcaRecoRealisticTRK_126X_mcRun3_2022_realistic_postEE_v2-v2/GEN-SIM-RECO"
]

datasets["Background_ZPlusJets"] = [
    "/DYJetsToMuMu_M-50_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Winter23Reco-TRKRealistic_AlcaRecoRealisticTRK_preEE_126X_mcRun3_2022_realistic_v4-v2/GEN-SIM-RECO",
    "/DYJetsToMuMu_M-50_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Winter23Reco-TRKRealistic_AlcaRecoRealisticTRK_REAL_126X_mcRun3_2022_realistic_postEE_v2-v3/GEN-SIM-RECO"
]

datasets["Background_QCD"] = [
    "/QCD_PT-15to7000_TuneCP5_13p6TeV_pythia8/Run3Winter23Reco-FlatPU0to80_126X_mcRun3_2023_forPU65_v1-v2/GEN-SIM-RECO"
]

dataset_name_to_request_name = {}
dataset_name_to_request_name["/WJetsToLNu_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Winter23Reco-TRKRealistic_AlcaRecoRealisticTRK_preEE_126X_mcRun3_2022_realistic_v4-v2/GEN-SIM-RECO"]  = "WJetsToLNu_Run3Winter23Reco_preEE_126X_mcRun3_2022"
dataset_name_to_request_name["/WJetsToLNu_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Winter23Reco-TRKRealistic_AlcaRecoRealisticTRK_126X_mcRun3_2022_realistic_postEE_v2-v2/GEN-SIM-RECO"] = "WJetsToLNu_Run3Winter23Reco_postEE_126X_mcRun3_2022"
dataset_name_to_request_name["/DYJetsToMuMu_M-50_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Winter23Reco-TRKRealistic_AlcaRecoRealisticTRK_preEE_126X_mcRun3_2022_realistic_v4-v2/GEN-SIM-RECO"] = "DYJetsToMuMu_Run3Winter23Reco_preEE_126X_mcRun3_2022"
dataset_name_to_request_name["/DYJetsToMuMu_M-50_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Winter23Reco-TRKRealistic_AlcaRecoRealisticTRK_REAL_126X_mcRun3_2022_realistic_postEE_v2-v3/GEN-SIM-RECO"] = "DYJetsToMuMu_Run3Winter23Reco_postEE_126X_mcRun3_2022"
dataset_name_to_request_name["/QCD_PT-15to7000_TuneCP5_13p6TeV_pythia8/Run3Winter23Reco-FlatPU0to80_126X_mcRun3_2023_forPU65_v1-v2/GEN-SIM-RECO"] = "QCD_PT-15to7000_Run3Winter23Reco_FlatPU0to80_126X_mcRun3_2023"

# -------------------------------------------------------------------------------------------------
def main():

    output_directory = "Bulk_CRAB_Scripts"
    if len(sys.argv) > 1: output_directory = sys.argv[1]

    os.mkdir(output_directory)
    os.chdir(output_directory)

    print("Generating CRAB Scripts in directory:", output_directory)

    crab_script_list = [] 

    signal_tags = [ "Signal_HToSSTo4B_MH125_MS50_CTau3000", "Signal_HToSSTo4B_MH250_MS120_CTau10000", "Signal_HToSSTo4B_MH350_MS80_CTau500", "Signal_HToSSTo4B_MH350_MS160_CTau10000" ]
    data_tags = [ "Data_EXOLLPJetHCAL_Run2023" ]
    for dataset_tag in data_tags: #["Data_DisplacedJet_Run2022"]: #"Data_DisplacedJet_Run2022"]: #signal_tags: #["Background_QCD"]: #"Background_ZPlusJets"]: #signal_tags: 

        i = 0
        for dataset_name in datasets[dataset_tag]:

            replacements = {
                "MYVAR_CRAB_OUTPUT_NAME": crab_output_dir, 
                "MYVAR_DATASET_NAME": dataset_name,
                "MYVAR_GOLDEN_JSON": "",
                "MYVAR_ISDATA": "False",
                "MYVAR_ISSIGNAL": "False",
                "MYVAR_RECO_FROM_RAW": "False", 
                "MYVAR_REQUEST_NAME": dataset_name.replace("/","_")[1:]+"_v4",
                "MYVAR_DATASET_TAG": dataset_name.replace("/","_")[1:]+"_v4",
                "MY_VAR_INPUTDBS": "global",
            }

            if "Data_"   in dataset_tag: 
                replacements["MYVAR_ISDATA"]   = "True"
                #replacements["MYVAR_EVENTS_PER_FILE"] = 50000
                if "Run2022" in dataset_tag: replacements["MYVAR_GOLDEN_JSON"] = "/eos/user/c/cmsdqm/www/CAF/certification/Collisions22/Cert_Collisions2022_355100_362760_Golden.json"
                if "Run2023" in dataset_tag: replacements["MYVAR_GOLDEN_JSON"] = "/eos/user/c/cmsdqm/www/CAF/certification/Collisions23/Cert_Collisions2023_366442_370790_Golden.json"
                if "Run2024" in dataset_tag: replacements["MYVAR_GOLDEN_JSON"] = "/eos/user/c/cmsdqm/www/CAF/certification/Collisions24/Cert_Collisions2024_378981_386951_Golden.json"
                if "Data_DisplacedJet_Run2022" in dataset_tag: 
                    replacements["MYVAR_RECO_FROM_RAW"] = "True"
                    replacements["MYVAR_EXTRACONFIG"]  = "\nconfig.Data.splitting         = 'EventAwareLumiBased'"
                    replacements["MYVAR_EXTRACONFIG"] += "\nconfig.Data.unitsPerJob       = 35000"
            if "Signal_" in dataset_tag: 
                replacements["MYVAR_ISSIGNAL"] = "True"
                replacements["MY_VAR_INPUTDBS"] = "phys03"  #private mc
                replacements["MYVAR_REQUEST_NAME"] = dataset_tag.replace("Signal_","") + "_batch" + str(i+1) + "_v4"
                replacements["MYVAR_DATASET_TAG"]  = dataset_tag.replace("Signal_","") + "_batch" + str(i+1) + "_v4"
                replacements["MYVAR_EXTRACONFIG"]  = "\nconfig.Data.splitting         = 'FileBased'"
                replacements["MYVAR_EXTRACONFIG"] += "\nconfig.Data.unitsPerJob       = 200"
            if "Background_" in dataset_tag:
                replacements["MYVAR_EXTRACONFIG"]  = "\nconfig.Data.partialDataset = True" # Just run over what is available
                replacements["MYVAR_REQUEST_NAME"] = dataset_name_to_request_name[dataset_name] + "_v4"
                replacements["MYVAR_DATASET_TAG"]  = dataset_name_to_request_name[dataset_name] + "_v4"
                #replacements["MYVAR_EVENTS_PER_FILE"] = "100000" # check
            #if "Data_ZMu_" in dataset_tag:
            #    replacements["MYVAR_EVENTS_PER_FILE"] = 10000 

            # Read the file
            with open(crab_filepath, "r") as f:
                content = f.read()

            # Replace each MYVAR_XXX with the user-defined value
            for key, value in replacements.items():
                content = re.sub(rf"\b{re.escape(key)}\b", value, content)

            output_crab_filepath = crab_filepath.split("/")[-1].replace("DO-NOT-EDIT", dataset_tag+"-"+str(i))

            print( output_crab_filepath )
            crab_script_list.append( output_crab_filepath )

            # Write the modified content back to the file (or to a new file)
            with open(output_crab_filepath, "w") as f:
                f.write(content)

            i += 1

    os.chdir( pwd )

    print("Done") #. To submit, run: ")
    #print("source <TODO> ")

    with open( os.path.join( output_directory, "submit.sh"), "w") as f:
        f.write("SUBMITLOG="+os.path.join( output_directory, "SubmitLog.txt" )+ "\n")
        f.write("touch $SUBMITLOG \n")
        f.write("echo \"Submitting CRAB jobs, output will be written to $SUBMITLOG\" \n")
        for file in crab_script_list: 
            f.write( "echo \"Submitting: "+ os.path.join( output_directory, file ) + "\" \n" )
            f.write( "crab submit -c " + os.path.join( output_directory, file ) + " >> $SUBMITLOG \n" )
        f.write("grep \"crab status\" $SUBMITLOG | sed \"s/Please use '//g\" | sed \"s/' to check how the submission process proceeds.//g\"  > "+ os.path.join( output_directory, "status.sh") + "\n" )
        f.write("echo \"You can monitor your jobs by running: `source "+os.path.join( output_directory, "status.sh")+" `\" \n" )

    with open( os.path.join( output_directory, "submit_dryrun.sh"), "w") as f:
        for file in crab_script_list: f.write( "crab submit -c " + os.path.join( output_directory, file ) + " --dryrun \n" )

    print("\nTo submit all, run:")
    print("source "+os.path.join( output_directory, "submit.sh") )

    print("\nTo submit all in --dryrun mode, run:")
    print("source "+os.path.join( output_directory, "submit_dryrun.sh") )

    

# -------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    main()

