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
#crab_output_dir = '/afs/cern.ch/work/g/gkopp/2022_LLP_analysis/CRAB_Workarea/NTuples_v6/'
#crab_output_dir = '/afs/cern.ch/work/f/fsimpson/2022_LLP_analysis/CRAB_Workarea/NTuples_v6/'
crab_output_dir = '/afs/cern.ch/work/k/kikenned/Run3-HCAL-LLP-NTupler/CRAB_Workarea/NTuples_v5/'

#datasets = {}

from crab_formatted_datasets import datasets

dataset_name_to_request_name = {}
dataset_name_to_request_name["/WJetsToLNu_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Winter23Reco-TRKRealistic_AlcaRecoRealisticTRK_preEE_126X_mcRun3_2022_realistic_v4-v2/GEN-SIM-RECO"]  = "WJetsToLNu_Run3Winter23Reco_preEE_126X_mcRun3_2022"
dataset_name_to_request_name["/WJetsToLNu_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Winter23Reco-TRKRealistic_AlcaRecoRealisticTRK_126X_mcRun3_2022_realistic_postEE_v2-v2/GEN-SIM-RECO"] = "WJetsToLNu_Run3Winter23Reco_postEE_126X_mcRun3_2022"
dataset_name_to_request_name["/DYJetsToMuMu_M-50_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Winter23Reco-TRKRealistic_AlcaRecoRealisticTRK_preEE_126X_mcRun3_2022_realistic_v4-v2/GEN-SIM-RECO"] = "DYJetsToMuMu_Run3Winter23Reco_preEE_126X_mcRun3_2022"
dataset_name_to_request_name["/DYJetsToMuMu_M-50_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Winter23Reco-TRKRealistic_AlcaRecoRealisticTRK_REAL_126X_mcRun3_2022_realistic_postEE_v2-v3/GEN-SIM-RECO"] = "DYJetsToMuMu_Run3Winter23Reco_postEE_126X_mcRun3_2022"
dataset_name_to_request_name["/QCD_PT-15to7000_TuneCP5_13p6TeV_pythia8/Run3Winter23Reco-FlatPU0to80_126X_mcRun3_2023_forPU65_v1-v2/GEN-SIM-RECO"] = "QCD_PT-15to7000_Run3Winter23Reco_FlatPU0to80_126X_mcRun3_2023"
dataset_name_to_request_name["/QCD_Bin-Pt-15to7000_TuneCP5_13p6TeV_pythia8/Run3Winter25Reco-FlatPU0to120_142X_mcRun3_2025_realistic_v7-v3/GEN-SIM-RECO"] = "QCD_Bin-Pt-15to7000_Run3Winter25Reco-FlatPU0to120_142X_mcRun3_2025_realistic_v7-v3"
dataset_name_to_request_name["/QCD_Bin-Pt-15to7000_TuneCP5_13p6TeV_pythia8/Run3Winter25Reco-FEVTOUTPUT_142X_mcRun3_2025_realistic_v7-v1/GEN-SIM-RECO"] = "QCD_Bin-Pt-15to7000_Run3Winter25Reco-FEVTOUTPUT_142X_mcRun3_2025_realistic_v7-v1"
datasets["Background_QCD"] = [
    #"/QCD_PT-15to7000_TuneCP5_13p6TeV_pythia8/Run3Winter23Reco-FlatPU0to80_126X_mcRun3_2023_forPU65_v1-v2/GEN-SIM-RECO",
    "/QCD_Bin-Pt-15to7000_TuneCP5_13p6TeV_pythia8/Run3Winter25Reco-FlatPU0to120_142X_mcRun3_2025_realistic_v7-v3/GEN-SIM-RECO",
    "/QCD_Bin-Pt-15to7000_TuneCP5_13p6TeV_pythia8/Run3Winter25Reco-FEVTOUTPUT_142X_mcRun3_2025_realistic_v7-v1/GEN-SIM-RECO"]


# -------------------------------------------------------------------------------------------------
def main():

    output_directory = "Bulk_CRAB_Scripts"
    if len(sys.argv) > 1: output_directory = sys.argv[1]

    os.mkdir(output_directory)
    os.chdir(output_directory)

    print("Generating CRAB Scripts in directory:", output_directory)

    crab_script_list = [] 

    signal_tags = [ "Signal_HToSSTo4B_MH125_MS50_CTau3000", "Signal_HToSSTo4B_MH250_MS120_CTau10000", "Signal_HToSSTo4B_MH350_MS80_CTau500", "Signal_HToSSTo4B_MH350_MS160_CTau10000" ]
    data_tags = [ "Data_DisplacedJet_Run2022" ] #"Data_EXOLLPJetHCAL_Run2023" ]
    zmu_tags = [ "Data_ZMu_Run2022", "Data_ZMu_Run2023" ]

    signal_central_tags = []
    for tag in datasets: 
        if "Signal_RAW_" in tag: signal_central_tags.append(tag)

    for dataset_tag in signal_tags: #data_tags: # signal_central_tags, signal_tags, zmu_tags
        i = 0
        for dataset_name in datasets[dataset_tag]:
            replacements = {
                "MYVAR_CRAB_OUTPUT_NAME": crab_output_dir, 
                "MYVAR_DATASET_NAME": dataset_name,
                "MYVAR_GOLDEN_JSON": "",
                "MYVAR_ISDATA": "False",
                "MYVAR_ISSIGNAL": "False",
                "MYVAR_RECO_FROM_RAW": "False", 
                "MYVAR_REQUEST_NAME": dataset_name.replace("/","_")[1:]+"_v6",
                "MYVAR_DATASET_TAG": dataset_name.replace("/","_")[1:]+"_v6",
                "MY_VAR_INPUTDBS": "global",
            }

            if "Data_" in dataset_tag: 
                replacements["MYVAR_ISDATA"]   = "True"
                #replacements["MYVAR_EVENTS_PER_FILE"] = 50000
                if "Run2022" in dataset_tag: replacements["MYVAR_GOLDEN_JSON"] = "/eos/user/c/cmsdqm/www/CAF/certification/Collisions22/Cert_Collisions2022_355100_362760_Golden.json"
                if "Run2023" in dataset_tag: replacements["MYVAR_GOLDEN_JSON"] = "/eos/user/c/cmsdqm/www/CAF/certification/Collisions23/Cert_Collisions2023_366442_370790_Golden.json"
                if "Run2024" in dataset_tag: replacements["MYVAR_GOLDEN_JSON"] = "/eos/user/c/cmsdqm/www/CAF/certification/Collisions24/Cert_Collisions2024_378981_386951_Golden.json"
                if "Data_DisplacedJet_Run2022" in dataset_tag: # Basically, if RAW # or "ZMu" in dataset_tag: 
                    replacements["MYVAR_RECO_FROM_RAW"] = "True"
                    replacements["MYVAR_EXTRACONFIG"]  = "\nconfig.Data.splitting         = 'EventAwareLumiBased'"
                    replacements["MYVAR_EXTRACONFIG"] += "\nconfig.Data.unitsPerJob       = 50000"
                    replacements["MYVAR_EXTRACONFIG"] += "\nconfig.JobType.numCores = 2"
                    replacements["MYVAR_EXTRACONFIG"] += "\nconfig.JobType.maxMemoryMB = 5000"
            elif "Signal_" in dataset_tag: 
                replacements["MYVAR_ISSIGNAL"] = "True"
                replacements["MYVAR_REQUEST_NAME"] = dataset_tag.replace("Signal_","") + "_batch" + str(i+1) + "_v6"
                replacements["MYVAR_DATASET_TAG"]  = dataset_tag.replace("Signal_","") + "_batch" + str(i+1) + "_v6"
                if "RAW" in dataset_tag:
                    replacements["MYVAR_RECO_FROM_RAW"] = "True"
                    replacements["MYVAR_EXTRACONFIG"]  = "\nconfig.Data.splitting         = 'EventAwareLumiBased'"
                    replacements["MYVAR_EXTRACONFIG"] += "\nconfig.Data.unitsPerJob       = 1000"
                    replacements["MYVAR_EXTRACONFIG"] += "\nconfig.JobType.numCores = 2"
                    replacements["MYVAR_EXTRACONFIG"] += "\nconfig.JobType.maxMemoryMB = 5000"
                else:
                    replacements["MY_VAR_INPUTDBS"] = "phys03"  #private mc
                    replacements["MYVAR_EXTRACONFIG"]  = "\nconfig.Data.splitting         = 'FileBased'"
                    replacements["MYVAR_EXTRACONFIG"] += "\nconfig.Data.unitsPerJob       = 200"
            elif "Background_" in dataset_tag:
                replacements["MYVAR_EXTRACONFIG"]  = "\nconfig.Data.partialDataset = True" # Just run over what is available
                replacements["MYVAR_REQUEST_NAME"] = dataset_name_to_request_name[dataset_name] + "_v6"
                replacements["MYVAR_DATASET_TAG"]  = dataset_name_to_request_name[dataset_name] + "_v6"
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

