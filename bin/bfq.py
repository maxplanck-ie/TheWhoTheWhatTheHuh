#!/usr/bin/env python3
import sys
import os
import datetime
import time
from bcl2fastq_pipeline import getConfig
from bcl2fastq_pipeline import findFlowCells
from bcl2fastq_pipeline import makeFastq
from bcl2fastq_pipeline import afterFastq
from bcl2fastq_pipeline import misc

def sleep(config) :
    time.sleep(float(config['Options']['sleepTime'])*60*60)

while True:
    #Read the config file
    config = getConfig.getConfig()
    if(config is None) :
        #There's no recovering from this!
        sys.exit("Error: couldn't read the config file!")

    #Get the next flow cell to process, or sleep
    config = findFlowCells.newFlowCell(config)
    if(config.get('Options','runID') == "") :
        sleep(config)
        continue

    #Ensure we have sufficient space
    if(misc.enoughFreeSpace(config) == False) :
        sys.stderr.write("Error: insufficient free space!\n")
        sys.stderr.flush()
        misc.errorEmail(config,"Error: insufficient free space!")
        sleep(config)
        continue

    startTime=datetime.datetime.now()

    #Make the fastq files
    try:
        makeFastq.bcl2fq(config)
    except :
        sys.stderr.write("Got an error in bcl2fq\n")
        sys.stderr.flush()
        misc.errorEmail(config,"Got an error in bcl2fq")
        sleep(config)
        continue

    #Run post-processing steps
    try :
        message = afterFastq.postMakeSteps(config)
    except :
        sys.stderr.write("Got an error during postMakeSteps\n")
        sys.stderr.flush()
        misc.errorEmail(config, "Got an error during postMakeSteps")
        sleep(config)
        continue

    #Copy over xml and FastQC stuff
    try:
        makeFastq.cpSeqFac(config)
    except :
        sys.stderr.write("Got an error in cpSeqFac\n")
        sys.stderr.flush()
        misc.errorEmail(config,"Got an error in cpSeqFac")
        sleep(config)
        continue

    message = misc.parseConversionStats(config)
    #Get more statistics and create PDFs
    try :
        message += "\n\n"+misc.parseConversionStats(config)
    except :
        sys.stderr.write("Got an error during parseConversionStats\n")
        sys.stderr.flush()
        misc.errorEmail(config, "Got an error during parseConversionStats")
        sleep(config)
        continue

    runTime = datetime.datetime.now()-startTime

    #Email finished message
    try :
        misc.finishedEmail(config, message, runTime)
    except :
        #Unrecoverable error
        sys.exit("Couldn't send the finished email! Quiting")

    #Mark the flow cell as having been processed
    findFlowCells.markFinished(config)

    #DEBUGGING!
    break
