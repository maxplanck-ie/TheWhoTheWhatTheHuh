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
import signal
from threading import Event
import urllib3

# Disable excess warning messages if we disable SSL checks
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
gotHUP = Event()

def breakSleep(signo, _frame):
    gotHUP.set()


def sleep(config) :
    gotHUP.wait(timeout=float(config['Options']['sleepTime'])*60*60)
    gotHUP.clear()

signal.signal(signal.SIGHUP, breakSleep)


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
    config, parkour_dict = bcl2fastq_pipeline.findFlowCells.newFlowCell(config)
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
    message = ""
    if not os.path.exists("{}/{}{}/bcl.done".format(config["Paths"]["outputDir"], config["Options"]["runID"], lanes)):
        try:
            print("pre_bcl2fastq")
            message += bcl2fastq_pipeline.makeFastq.bcl2fq(config, parkour_dict)
            print(message)
            open("{}/{}{}/bcl.done".format(config["Paths"]["outputDir"], config["Options"]["runID"], lanes), "w").close()
        except :
            syslog.syslog("Got an error in bcl2fq\n")
            #bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), "Got an error in bcl2fq")
            sleep(config)
            continue

    #Fix the file names (prepend "Project_" and "Sample_" and such as appropriate)
    if not os.path.exists("{}/{}{}/files.renamed".format(config["Paths"]["outputDir"], config["Options"]["runID"], lanes)):
        try:
            bcl2fastq_pipeline.makeFastq.fixNames(config)
            open("{}/{}{}/files.renamed".format(config["Paths"]["outputDir"], config["Options"]["runID"], lanes), "w").close()
        except :
            syslog.syslog("Got an error in fixNames\n")
            #bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), "Got an error in fixNames")
            sleep(config)
            continue

    #Run post-processing steps
    try :
         message += bcl2fastq_pipeline.afterFastq.postMakeSteps(config)
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
        print("skip transfer")
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

    #Update parkour, errors are non-fatal here
    try:
        message += bcl2fastq_pipeline.misc.jsonParkour(config, message)
    except:
        syslog.syslog("Received an error while updating Parkour!\n")
        bcl2fastq_pipeline.misc.errorEmail(config, sys.exc_info(), "Got an error while updating Parkour!")

    #Email finished message
    try :
         bcl2fastq_pipeline.misc.finishedEmail(config, message, runTime, transferTime)
    except :
         #Unrecoverable error
         syslog.syslog("Couldn't send the finished email! Quiting")
         sys.exit()

    #Mark the flow cell as having been processed
    bcl2fastq_pipeline.findFlowCells.markFinished(config)
