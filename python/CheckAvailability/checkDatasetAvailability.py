#!/usr/bin/python

import os
import datetime
import time
import subprocess
import glob
import sys
import json

if (len(sys.argv) -1 < 1):
    print ("Error. Not enough arguments provided.\n")
    print ("Usage: python checkDatasetAvailability.py [DatasetListFile]  \n")
    exit()


datasetListFilename = sys.argv[1]
os.system("rm DiskStatus.txt -v")
outfile = open("DiskStatus.txt", "w")
blockfile = open("Blocks.txt", "w")
tempfile = open(datasetListFilename,"r")
templines = tempfile.readlines()
for line in templines:
    datasetName = line.strip()

    # datasetSize = 0
    # command = "dasgoclient -query=\"dataset dataset=" + datasetName + "\" -json > tmpOutput.json" 
    # #print command
    # os.system(command)
    # jsonFile = open("tmpOutput.json","r")
    # data = json.load(jsonFile)
    # for p in data[2]["dataset"]:    
    #     #print p
    #     datasetSize = p["size"] / (1024*1024*1024*1024.)
    # print datasetName + " : " + "{0:.2f}".format(datasetSize) + " TB"

    print (datasetName)

    #First get all the blocks of the dataset
    # commandPrintBlocks = "/afs/cern.ch/user/v/valya/public/dasgoclient/dasgoclient -query=\"block dataset=" + datasetName + "\" -json > tmpBlocks.json"
    commandPrintBlocks = "dasgoclient -query=\"block dataset=" + datasetName + "\" -json > tmpBlocks.json"
    #command = "dasgoclient -query=\"site dataset=" + datasetName + "\" -json > tmpOutput.json"
    #print command
    os.system(commandPrintBlocks)

    jsonFile = open("tmpBlocks.json","r")
    data = json.load(jsonFile)

    nBlocksTotal = 0
    nBlocksOnDisk = 0
    nFilesTotal = 0
    nFilesOnDisk = 0
    blocksNotOnDisk = []


    for p in data:        
        blockName = p["block"][0]["name"]
        commandRUCIOCheckBlockReplica = 'echo "source /cvmfs/cms.cern.ch/rucio/setup-py3.sh; rucio list-dataset-replicas cms:' + blockName + ' --deep --csv > tmpReplica.csv" > tmp.cmd; sh tmp.cmd'
        #print commandRUCIOCheckBlockReplica
        os.system(commandRUCIOCheckBlockReplica)
        blockfile.write(blockName+"\n")

        nBlocksTotal += 1        
        thisBlockFullyOnDisk = False
        nFilesInThisBlock = 0

        tmpfile = open('tmpReplica.csv', 'r')
        tmplines = tmpfile.readlines()
        for line in tmplines:
            tmpdata = line.split(",")
            sitename = tmpdata[0]
            tmpNFoundFiles = int(tmpdata[1])
            tmpNTotalFiles = int(tmpdata[2])
            #print (sitename + " : " + str(tmpNFoundFiles) + " : " + str(tmpNTotalFiles) )
            
            if (nFilesInThisBlock == 0):
                nFilesInThisBlock = tmpNTotalFiles

            #print ( "Disk ? " + str(not ("Tape" in tmpdata[0])))
            #print ( "Complete ? " + str(tmpNFoundFiles == tmpNTotalFiles))

            if (not ("Tape" in tmpdata[0]) and tmpNFoundFiles == tmpNTotalFiles):
                thisBlockFullyOnDisk = True
                
            #print ("Complete : " + str(thisBlockFullyOnDisk))            
        tmpfile.close()

        nFilesTotal  += nFilesInThisBlock
        if (thisBlockFullyOnDisk):
            nBlocksOnDisk += 1
            nFilesOnDisk += nFilesInThisBlock
        else:
            blocksNotOnDisk.append(blockName)

    print ( "Files On Disk : " + str(nFilesOnDisk) + " / " + str(nFilesTotal) + " ; Blocks On Disk : " + str(nBlocksOnDisk) + " / " + str(nBlocksTotal) )

    outfile.write(datasetName + " : " + str( 100 * float(nFilesOnDisk)/float(nFilesTotal)) + " , Blocks On Disk : " +  str(nBlocksOnDisk) + " / " + str(nBlocksTotal) + "\n" )

    print ("Blocks not on Disk:")
    for b in blocksNotOnDisk:
        print (b)

print ("DONE")