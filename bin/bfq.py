#!/usr/bin/env python3
import sys
import os
import datetime
import time
import syslog
import bcl2fastq_pipeline.getConfig
import bcl2fastq_pipeline.findFlowCells
import bcl2fastq_pipeline.makeFastq
import bcl2fastq_pipeline.afterFastq
import bcl2fastq_pipeline.misc
import bcl2fastq_pipeline.galaxy
import importlib

def sleep(config) :
    time.sleep(float(config['Options']['sleepTime'])*60*60)

while True:
    #Reimport to allow reloading a new version
    importlib.reload(bcl2fastq_pipeline.getConfig)
    importlib.reload(bcl2fastq_pipeline.findFlowCells)
    importlib.reload(bcl2fastq_pipeline.makeFastq)
    importlib.reload(bcl2fastq_pipeline.afterFastq)
    importlib.reload(bcl2fastq_pipeline.misc)
    importlib.reload(bcl2fastq_pipeline.galaxy)

    #Read the config file
    config = bcl2fastq_pipeline.getConfig.getConfig()
    if(config is None) :
        #There's no recovering from this!
        sys.exit("Error: couldn't read the config file!")

    #Get the next flow cell to process, or sleep
    config = bcl2fastq_pipeline.findFlowCells.newFlowCell(config)
    if(config.get('Options','runID') == "") :
        sleep(config)
        continue

    #Ensure we have sufficient space
    if(bcl2fastq_pipeline.misc.enoughFreeSpace(config) == False) :
        syslog.syslog("Error: insufficient free space!\n")
        bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), "Error: insufficient free space!")
        sleep(config)
        continue

    startTime=datetime.datetime.now()

    #Make the fastq files, if not already done
    lanes = config["Options"]["lanes"]
    if lanes != "":
        lanes = "_lanes{}".format(lanes)
    if not os.path.exists("{}/{}{}/bcl.done".format(config["Paths"]["outputDir"], config["Options"]["runID"], lanes)):
        try:
            bcl2fastq_pipeline.makeFastq.bcl2fq(config)
            open("{}/{}{}/bcl.done".format(config["Paths"]["outputDir"], config["Options"]["runID"], lanes), "w").close()
        except :
            syslog.syslog("Got an error in bcl2fq\n")
            bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), "Got an error in bcl2fq")
            sleep(config)
            continue

    #Run post-processing steps
    try :
        message = bcl2fastq_pipeline.afterFastq.postMakeSteps(config)
    except :
        syslog.syslog("Got an error during postMakeSteps\n")
        bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), "Got an error during postMakeSteps")
        sleep(config)
        continue

    #Get more statistics and create PDFs
    try :
        message += "\n\n"+bcl2fastq_pipeline.misc.parseConversionStats(config)
    except :
        syslog.syslog("Got an error during parseConversionStats\n")
        bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), "Got an error during parseConversionStats")
        sleep(config)
        continue

    #Copy over xml, FastQC, and PDF stuff
    try:
        message += bcl2fastq_pipeline.makeFastq.cpSeqFac(config)
    except :
        syslog.syslog("Got an error in cpSeqFac\n")
        bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), "Got an error in cpSeqFac")
        sleep(config)
        continue

    runTime = datetime.datetime.now()-startTime
    startTime = datetime.datetime.now()

    #Transfer data to groups
    try : 
        message += bcl2fastq_pipeline.misc.transferData(config)
    except :
        syslog.syslog("Got an error during transferData\n")
        bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), "Got an error during transferData")
        sleep(config)
        continue

    #Upload to Galaxy
    try :
        message += bcl2fastq_pipeline.galaxy.linkIntoGalaxy(config)
    except:
        syslog.syslog("Got an error while uploading to Galaxy!\n")
        bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), "Got an error while uploading to Galaxy!")
        sleep(config)

    transferTime = datetime.datetime.now()-startTime

    #Email finished message
    try :
        bcl2fastq_pipeline.misc.finishedEmail(config, message, runTime, transferTime)
    except :
        #Unrecoverable error
        syslog.syslog("Couldn't send the finished email! Quiting")
        sys.exit()

    #Mark the flow cell as having been processed
    bcl2fastq_pipeline.findFlowCells.markFinished(config)
