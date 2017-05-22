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
import syslog

#Returns True on processed, False on unprocessed
def flowCellProcessed(config) :
    lanes = config.get("Options", "lanes")
    if lanes != "": 
        lanes = "_lanes{}".format(lanes)

    path = "%s/%s%s/fastq.made" % (config.get("Paths","outputDir"), config.get("Options","runID"), lanes)
    if os.access(path, os.F_OK):
        return True
    return False


def getSampleSheets(d):
    """
    Provide a list of output directories and sample sheets
    """
    ss = glob.glob("%s/SampleSheet*.csv" % d)

    if len(ss) == 0:
        return ([None], [None], [''])

    laneOut = []
    bcLens = []
    for sheet in ss:
        # get the lanes
        lanes = []
        f = open(sheet)
        inData = False
        lastLine = None
        colNum = None
        indexCols = [None, None]
        bcLensAppended = False
        for line in f:
            bcLen = '0,0'
            if inData is False:
                if line.startswith("[Data]"):
                    inData = True
                    continue
            else:
                cols = line.strip().split(",")
                if lastLine is True:
                    if indexCols[0] is not None:
                        bcLen = "{}".format(len(cols[indexCols[0]]))
                    if indexCols[1] is not None:
                        bcLen = "{},{}".format(bcLen, len(cols[indexCols[1]]))
                    else:
                        bcLen += ",0"
                    bcLens.append(bcLen)
                    bcLensAppended = True
                    lastLine = False

                # Handle barcodes (once)
                if "index" in cols and lastLine is None:
                    indexCols[0] = cols.index("index")
                    if "index2" in cols:
                        indexCols[1] = cols.index("index2")
                        lastLine = True
                    lastLine = True

                if "Lane" in cols:
                    colNum = cols.index("Lane")
                    lastLine = True
                    continue

                if colNum is not None:
                    if cols[colNum] != '' and int(cols[colNum]) not in lanes:
                        lanes.append(int(cols[colNum]))
        if len(lanes) > 0:
            lanes = sorted(lanes)
            laneOut.append("_".join(["{}".format(x) for x in lanes]))
        if bcLensAppended == False:
            # This only happens with runs having no Lane or index in the sample sheet
            bcLens.append(None)
    if len(ss) == 1:
        laneOut = [None]
    return ss, laneOut, bcLens


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
    dirs = glob.glob("%s/*_SN*_*/RTAComplete.txt" % config.get("Paths","baseDir"))
    dirs.extend(glob.glob("%s/*_NB*_*/RTAComplete.txt" % config.get("Paths","baseDir")))
    dirs.extend(glob.glob("%s/*_M*_*/RTAComplete.txt" % config.get("Paths","baseDir")))
    dirs.extend(glob.glob("%s/*_J*_*/RTAComplete.txt" % config.get("Paths","baseDir")))
    for d in dirs :
        #Get the flow cell ID (e.g., 150416_SN7001180_0196_BC605HACXX)
        config.set('Options','runID',d.split("/")[-2])

        # Before 1703 only a single sample sheet was supported
        if config.get("Options","runID")[:4] < "1703":
            continue

        sampleSheet, lanes, bcLens = getSampleSheets(os.path.dirname(d))

        for ss, lane, bcLen in zip(sampleSheet, lanes, bcLens):
            config.set('Options','runID',d.split("/")[-2])
            if lane is not None:
                config.set("Options","lanes",lane)
            else:
                config.set("Options","lanes","")
            if ss is None:
                ss = ''
            if bcLen is not None and bcLen is not '':
                config.set("Options","bcLen",bcLen)
            else:
                config.set("Options","bcLen","0,0")
    
            if flowCellProcessed(config) is False:
                syslog.syslog("Found a new flow cell: %s\n" % config.get("Options","runID"))
                config.set("Options","sampleSheet",ss)
                return config
            else :
                config.set("Options","runID","")
    config.set("Options","runID","")
    return config


def markFinished(config) :
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes = "_lanes{}".format(lanes)

    open("%s/%s%s/fastq.made" % (config["Paths"]["outputDir"], config["Options"]["runID"], lanes), "w").close()

'''
This function needs to be run after newFlowCell() returns with config.runID
filled in. It creates the output directories.
'''
def MakeTargetDirs(config) :
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes = "_lanes{}".format(lanes)

    assert(config["Paths"]["runID"] != None)
    os.mkdirs("%s/%s%s" % (config["Paths"]["outputDir"], config["Options"]["runID"], lanes))
