#!/usr/bin/env python3
import sys
import os
import datetime
from bcl2fastq_pipeline import getConfig
from bcl2fastq_pipeline import findFlowCells
from bcl2fastq_pipeline import makeFastq
from bcl2fastq_pipeline import afterFastq
from bcl2fastq_pipeline import misc

while True:
    #Read the config file
    config = getConfig.getConfig()
    if(config is None) :
        #There's no recovering from this!
        sys.exit("Error: couldn't read the config file!")

    #Get the next flow cell to process, or sleep
    config = findFlowCells.newFlowCell(config)
    if(config.get('Options','runID') == None) :
        sleep(config['Options']['sleepTime']*60*60)
        continue

    #Ensure we have sufficient space
    if(misc.enoughFreeSpace(config) == False) :
        print("Error: insufficient free space!")
        errorEmail(config,"Error: insufficient free space!")
        sleep(config['Options']['sleepTime']*60*60)
        continue

    startTime=datetime.datetime.now()

    #Make the fastq files
    try:
        makeFastq.bcl2fq(config)
    except :
        print("Got an error in bcl2fq")
        errorEmail(config,"Got an error in bcl2fq")
        sleep(config['Options']['sleepTime']*60*60)
        continue

    #Copy over xml and InterOp/
    try:
        makeFastq.cpXmlInterOp(config)
    except :
        print("Got an error in cpXmlInteOp")
        errorEmail(config,"Got an error in cpXmlInteOp")
        sleep(config['Options']['sleepTime']*60*60)
        continue

    #Run post-processing steps
    try :
        message = afterFastq.postMakeSteps(config)
    except :
        print("Got an error during postMakeSteps")
        errorEmail(config, "Got an error during postMakeSteps")
        sleep(config['Options']['sleepTime']*60*60)
        continue

    runTime = datetime.datetime.now()-startTime

    #Email finished message
    #try :
    misc.finishedEmail(config, message, runTime)
    #except :
    #    #Unrecoverable error
    #    sys.exit("Couldn't send the finished email! Quiting")

    #Mark the flow cell as having been processed
    open("%s/%s/fastq.made" % (config.get("Paths","outputDir"),config.get("Options","runID")), "w").close()

    break
