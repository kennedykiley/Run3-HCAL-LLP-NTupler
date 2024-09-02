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
            '/Muon0/Run2023D-ZMu-PromptReco-v2/RAW-RECO']

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
def edit_cfg():
    for i in range(len(datasets)):
        print(i)
        print(datasets[i])
        dataset_replace = "datasetnames = [\'" + datasets[i] + "\']"
        print(dataset_replace)
        sed('datasetnames = \[.*$',dataset_replace,'test.py',count=1)
"""
Things that need to change in the displaced jet ntuplizer python script:
number which dataset to look at
datasetnames from inputs given here
config.General.workArea change name, make non-local
config.General.requestName to add version number
"""

# -----------------------------------------------------------------------------------------
if __name__ == "__main__":
    sed('regex.*$','regex foo replaced full','test.py',count=1)
    # using .*$ enables replacing to the end of the line
    edit_cfg()