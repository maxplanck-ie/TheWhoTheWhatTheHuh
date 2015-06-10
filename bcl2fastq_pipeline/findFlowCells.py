'''
This file includes anything involved in finding new flow cells to process.

Note that this includes anything done after a flow cell has been processed,
such as marking it as having been processed and sending emails.
'''

import os
import sys
import smtplib
import glob
from email.mime.text import MIMEText

#Returns 1 on processed, 0 on unprocessed
def flowCellProcessed(config) :
    if(os.access("%s/%s/casava.finished" % (config.get("Paths","outputDir"), config.get("Options","runID")), os.F_OK)) :
        return True
    if(os.access("%s/%s/fastq.made" % (config.get("Paths","outputDir"), config.get("Options","runID")), os.F_OK)) :
        return True
    return False

'''
Iterate over all folders in config.baseDir from machine SN7001180. For each,
see if it contains an RTAComplete.txt file, as this signifies the completion
of a sequencing run.

If a run has finished, we then check to see if it's already been processed.
Runs are marked as having been processed if they appear in config.finalDir
and have a file called casava.finished or fastq.made. casava.finished is
produced by the old pipeline, while this one creates fastq.made.

This function always returns its configuration. If there's a new flow cell to
process, then the runID is filled in. Otherwise, that's set to None.
'''
def newFlowCell(config) :
    dirs = glob.glob("%s/*SN7001180*/RTAComplete.txt" % config.get("Paths","baseDir"))
    for d in dirs :
        #Get the flow cell ID (e.g., 150416_SN7001180_0196_BC605HACXX)
        config.set('Options','runID',d.split("/")[-2])

        if(flowCellProcessed(config) is False) :
            sys.stderr.write("Found a new flow cell: %s\n" % config.get("Options","runID"))
            sys.stderr.flush()
            return config
        else :
            config.set("Options","runID","")
    return config

def markFinished(config) :
    open("%s/%s/fastq.made" % (config["Paths"]["outputDir"], config["Options"]["runID"]), "w").close()

'''
This function needs to be run after newFlowCell() returns with config.runID
filled in. It creates the output directories.
'''
def MakeTargetDirs(config) :
    assert(config["Paths"]["runID"] != None)
    os.mkdirs("%s/%s" % (config["Paths"]["outputDir"], config["Options"]["runID"]))
