"""
CRAB job submission wrapper
Author: Gillian Kopp <gkopp@princeton.edu>

Submits CRAB jobs on lxplus for displaced jet analysis, running on RECO datasets or skims
"""

import argparse, sys, os, time
from datetime import datetime

import re
import shutil
from tempfile import mkstemp

# -----------------------------------------------------------------------------------------
# Section to modify datasets! 
datasets = ['/Muon0/Run2023B-ZMu-PromptReco-v1/RAW-RECO',
            '/Muon0/Run2023C-ZMu-PromptReco-v1/RAW-RECO',
            '/Muon0/Run2023C-ZMu-PromptReco-v2/RAW-RECO',
            '/Muon0/Run2023C-ZMu-PromptReco-v3/RAW-RECO',
            '/Muon0/Run2023C-ZMu-PromptReco-v4/RAW-RECO',
            '/Muon0/Run2023D-ZMu-PromptReco-v1/RAW-RECO', 
            '/Muon0/Run2023D-ZMu-PromptReco-v2/RAW-RECO',
            '/HToSSTo4B_MH350_MS80_CTau500/lpclonglived-crab_PrivateProduction_Summer23BPix_DR_step2_RECOSIM_HToSSTo4B_MH350_MS80_CTau500_batch1_v1-6c03a81f0d97498cab5c296ab3fa9a76/USER',
            '/QCD_PT-15to7000_TuneCP5_13p6TeV_pythia8/Run3Winter23Reco-FlatPU0to120_126X_mcRun3_2023_forPU65_v1-v2/GEN-SIM-RECO']
# LLP MC will have CTau*USER
# QCD MC will have QCD*GEN-SIM-RECO
# look at string in datasets to determine what type this is, and set isData, isSignal accordingly 

number = 0

version = "version3"

# -----------------------------------------------------------------------------------------
def sed(pattern, replace, source, dest=None, count=0):
    """
    From https://stackoverflow.com/questions/12714415/python-equivalent-to-sed
    Reads a source file and writes the destination file.
    In each line, replaces pattern with replace.

    Args:
        pattern (str): pattern to match (can be re.pattern)
        replace (str): replacement str
        source  (str): input filename
        count (int): number of occurrences to replace
        dest (str):   destination filename, if not given, source will be over written.        
    """

    fin = open(source, 'r')
    num_replaced = count

    if dest:
        fout = open(dest, 'w')
    else:
        fd, name = mkstemp()
        fout = open(name, 'w')

    for line in fin:
        out = re.sub(pattern, replace, line)
        fout.write(out)

        if out != line:
            num_replaced += 1
        if count and num_replaced > count:
            break
    try:
        fout.writelines(fin.readlines())
    except Exception as E:
        raise E

    fin.close()
    fout.close()

    if not dest:
        shutil.move(name, source) 

# -----------------------------------------------------------------------------------------
def edit_cfg(data_num = 0):
    """
    Things that need to change in the displaced jet ntuplizer python script:
    number which dataset to look at
    datasetnames from inputs given here
    config.General.workArea change name, make non-local
    config.General.requestName to add version number
    """
    
    # Parse string in datasets to determine what data type this is, and set isData, isSignal accordingly. Use this to determine which cfg file must be edited
    var_isData = 'False'
    var_isSignal = 'False'
    fileToEdit = 'crab_DisplacedHcalJetNTuplizer_cfg.py'
    if ("CTau" in datasets[data_num] and "USER" in datasets[data_num]): 
        var_isSignal = 'True'
        fileToEdit = 'crab_DisplacedHcalJetNTuplizer_MC_cfg.py'
    elif ("Run2023" in datasets[data_num] and "PromptReco" in datasets[data_num]): 
        var_isData = 'True'
    elif ("QCD" in datasets[data_num] and "GEN-SIM-RECO" in datasets[data_num]): 
        fileToEdit = 'crab_DisplacedHcalJetNTuplizer_QCD_cfg.py'

    # edit DisplacedHcalJetNTuplizer.py defaults
    sed('.*# isData wrapper','    ' + var_isData + ', # default value # isData wrapper','DisplacedHcalJetNTuplizer.py',count=1)
    sed('.*# isSignal wrapper','    ' + var_isSignal + ', # default value # isSignal wrapper','DisplacedHcalJetNTuplizer.py',count=1)
    # using .*$ enables replacing to the end of the line

    # setting up variables, string parsing, and directories / dates
    dataset_replace = "datasetnames = [\'" + datasets[data_num] + "\']"
    dataset = filter(None, datasets[data_num].split('/'))
    dataset = list(dataset)
    skim_type = filter(None, dataset[1].split('-'))
    skim_type = list(skim_type)
    pwd = os.getcwd()
    date = datetime.now().strftime("_%Y%m%d")
    timestamp = datetime.now().strftime("_%Y%m%d_%H%M%S")

    # edit crab_DisplacedHcalJetNTuplizer*_cfg.py
    sed('.*# dataset wrapper',dataset_replace + '# dataset wrapper',fileToEdit,count=1)
    sed('.*# number wrapper','number = 0 # starting at 0 -> refers to datasetnames # number wrapper',fileToEdit,count=1)

    workArea = "'" + pwd + "/../../../../../crab_" + skim_type[1] + date + "'"                # in 2022_LLP_analysis directory
    requestName = "'" + dataset[0]+'_'+dataset[1]+'_'+dataset[2]+timestamp+'_'+version + "'"
    outputDatasetTag = "'" + dataset[1]+'_'+dataset[2]+timestamp+'_'+version + "'"

    if (fileToEdit == 'crab_DisplacedHcalJetNTuplizer_MC_cfg.py'): 
        workArea = "'" + pwd + "/../../../../../crab_signalMC" + date + "'"
        requestName = "'" + dataset[0]+'_'+dataset[2]+'_submission'+str(data_num)+timestamp+'_'+version + "'"
        outputDatasetTag = "'" + dataset[0]+'_'+dataset[2]+'_submission'+str(data_num)+timestamp+'_'+version + "'"
    elif (fileToEdit == 'crab_DisplacedHcalJetNTuplizer_QCD_cfg.py'): 
        workArea = "'" + pwd + "/../../../../../crab_QCD_MC" + date + "'"
        requestName = "'QCD_PT_13p6TeV_"+dataset[2]+'_submission'+str(data_num)+timestamp+'_'+version + "'" # TODO needs work, dataset[0] has '-' and configs don't work with that
        outputDatasetTag = "'QCD_PT_13p6TeV_"+dataset[2]+'_submission'+str(data_num)+timestamp+'_'+version + "'"

    sed('.*# workArea wrapper','config.General.workArea        = ' + workArea + ' # workArea wrapper',fileToEdit,count=1)
    sed('.*# requestName wrapper','config.General.requestName     = ' + requestName + ' # requestName wrapper',fileToEdit,count=1)
    sed('.*# outputDatasetTag wrapper','config.Data.outputDatasetTag = ' + outputDatasetTag + ' # outputDatasetTag wrapper',fileToEdit,count=1)

    if (data_num == 4): submit_cfg(fileToEdit, workArea)

# -----------------------------------------------------------------------------------------
def submit_cfg(cfg_file = "", area = ""):
    """
    Submit the displaced jet ntuplizer python script:
    cmsenv
    scram b -j 8
    crab submit -c crab_DisplacedHcalJetNtuplizer_cfg.py
    """
    os.system("echo 'Setting up environment' ")
    os.system("cmsenv")
    os.system("scram b -j 8")
    os.system("pwd")
    os.system("echo 'Submitting crab jobs' ")
    os.system("crab submit -c " + cfg_file) # + " --dir=" + area)
    os.system("echo ' ' ")

# -----------------------------------------------------------------------------------------
if __name__ == "__main__":
    for i in range(len(datasets)):
        if (i == 4 or i > 6): edit_cfg(i)