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
import xml.etree.ElementTree as ET
from pyBarcodes import getStats


#Returns True on processed, False on unprocessed
def flowCellProcessed(config) :
    lanes = config.get("Options", "lanes")
    if lanes != "": 
        lanes = "_lanes{}".format(lanes)

    path = "%s/%s%s/fastq.made" % (config.get("Paths","outputDir"), config.get("Options","runID"), lanes)
    if os.access(path, os.F_OK):
        return True
    return False


# Get the number of lanes in the run. This might not match the number of lanes in the sampleSheet
def getNumLanes(d):
    try:
        tree = ET.parse("{}/RunInfo.xml".format(d))
        root = tree.getroot()[0]
        numLanes = root.findall("FlowcellLayout")[0]
        return int(numLanes.get("LaneCount"))
    except:
        return 1


def revComp(s):
    """Reverse complement a primer sequence"""
    d = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    s = s[::-1]
    x = [d[c] for c in s]
    return "".join(x)


def formatHeaderLine(cols, colLabs, indexCols, storeLanes):
    """Format a single header line in a sample sheet"""
    l = []
    if storeLanes is True:
        l.append("Lane")
    if colLabs[1] is not None:
        l.append("Sample_ID")
    if colLabs[2] is not None:
        l.append("Sample_Name")
    if indexCols[0] is not None and len(cols[indexCols[0]]) > 0:
        l.append("index")
    if indexCols[1] is not None and len(cols[indexCols[1]]) > 0:
        l.append("index2")
    if colLabs[3] is not None:
        l.append("Sample_Project")
    return ",".join(l)


def formatLine(cols, colLabs, indexCols, storeLanes):
    """Format a single line in a sample sheet"""
    l = []
    if storeLanes is True and colLabs[0] is not None:
        l.append(cols[colLabs[0]])
    if colLabs[1] is not None:
        l.append(cols[colLabs[1]])
    if colLabs[2] is not None:
        l.append(cols[colLabs[2]])
    if indexCols[0] is not None and len(cols[indexCols[0]]) > 0:
        l.append(cols[indexCols[0]])
    if indexCols[1] is not None and len(cols[indexCols[1]]) > 0:
        l.append(cols[indexCols[1]])
    if colLabs[3] is not None:
        l.append(cols[colLabs[3]])
    return ",".join(l)


def reformatSS(rv):
    """This is used in parseSampleSheet to reformat the output for upstream use"""
    ss = []
    laneOut = []
    bcLens = []
    nLanes = 0

    for k, v in rv.items():
        ss.append("\n".join(v[0]))
        lanes = ""
        if len(v[1]) > 0:
            nLanes += len(v[1])
            lanes = "_".join(["{}".format(x) for x in sorted(list(v[1]))])
        laneOut.append(lanes)
        bcLens.append(v[2])

    if len(ss) < 2 and nLanes == 8:
        laneOut = None
    return ss, laneOut, bcLens


def getReadLengths(basePath):
    """
    Parse RunInfo.xml to get the read lengths
    Return them as a list
    """
    tree = ET.parse("{}/RunInfo.xml".format(basePath))
    root = tree.getroot()[0]
    offsets = [1]
    for node in root.iter("Read"):
        offsets.append(int(node.get("NumCycles")) + offsets[-1])
    return offsets[1:]


def handleRevComp(d, basePath):
    """
    Input is a dictionary with masks as keys and values as lists with 3 items: output sample sheet(s) (list of lines), lane(s) (set), barcode lengths (string)

    If there's no barcode 2, simply return the lists as is. Otherwise, see if the barcodes match better what the sequencer saw or the rev. comp.
    In the latter case, rev. comp. and then return.
    """
    # Empty sample sheet
    if not d or not len(d):
        return d

    # No second barcode
    tot = 0
    for v in d.keys():  # A list of barcode lengths, e.g., ['6,0', '8,8']
        tot += int(v.split(",")[1])
    if tot == 0:
        return d

    # Get the flow cell type
    machine = os.path.split(basePath.rstrip("/"))[1].split("_")[1]
    if machine[0] == 'N':
        runType = "NextSeq"
    elif machine[0] == 'M':
        runType = "MiSeq"
    elif machine[0] == 'S':
        runType = "HiSeq2500"
    else:
        runType = "HiSeq3000"

    # At least 1 lane has a barcode 2
    readOffsets = getReadLengths(basePath)
    sampleSheets = []
    lanes = []
    masks = []
    vals = list(d.values())
    for i in range(len(vals)):  # Iterate over each mask possibility
        ss = vals[i][0][2:]  # The header is stripped
        localLanes = vals[i][1]  # This is actually a set
        mask = vals[i][2]
        if int(mask.split(",")[1]) == 0:  #Skip, no second barcode
            sampleSheets.append(vals[i][0])
            lanes.append(localLanes)
            masks.append(mask)
            continue

        hasLane = True
        if localLanes == '' or len(localLanes) == 0:
            hasLane = False
            lane = 1
            localLanes = [1]
        else:
            localLanes = sorted(list(localLanes))

        # Make a dictionary with key the lane and values the entries (as lists of lists)
        d2 = dict()
        for line in ss:
            if hasLane:
                lane = int(line[0])
            if lane not in d2:
                d2[lane] = list()
            d2[lane].append(line)

        # for each lane, get the observed barcode frequency
        cycles = list(range(readOffsets[0], readOffsets[0] + int(mask.split(",")[0])))
        cycles.extend(list(range(readOffsets[1], readOffsets[1] + int(mask.split(",")[1]))))
        finalSS = []
        outputLanes = set()
        for lane in localLanes:
            barcodes = getStats(basePath, runType, cycles, lane)

            totF = 0.0
            totR = 0.0
            # See what the total is if we used the barcodes as given
            for line in d2[lane]:
                line = line.split(",")
                # The .strip() stuff is due to one flow cell having extra spaces in it.
                if hasLane:
                    bcF = "{}{}".format(line[3].strip(), line[4].strip())
                    bcR = "{}{}".format(line[3].strip(), revComp(line[4].strip()))
                else:
                    bcF = "{}{}".format(line[2].strip(), line[3].strip())
                    bcR = "{}{}".format(line[2].strip(), revComp(line[3].strip()))
                if bcF in barcodes:
                    totF += barcodes[bcF]
                if bcR in barcodes:
                    totR += barcodes[bcR]
            if totR > totF:
                for idx in range(len(d2[lane])):
                    cols = d2[lane][idx].split(",")
                    if hasLane:
                        cols[4] = revComp(cols[4])
                    else:
                        cols[3] = revComp(cols[3])
                    d2[lane][idx] = ",".join(cols)
            if hasLane:
                outputLanes.add(lane)
            finalSS.extend(d2[lane])

        # Add the header lines to finalSS and then make it a big string
        if hasLane:
            lanes.append(localLanes)
        masks.append(mask)
        finalSS.insert(0, "[Data]")
        if hasLane:
            finalSS.insert(1, "Lane,Sample_ID,Sample_Name,index,index2,Sample_Project")
        else:
            finalSS.insert(1, "Sample_ID,Sample_Name,index,index2,Sample_Project")
        sampleSheets.append(finalSS)

    rv = dict()
    for i in range(len(masks)):
        if hasLane:
            rv[masks[i]] = [sampleSheets[i], lanes[i], masks[i]]
        else:
            rv[masks[i]] = [sampleSheets[i], set(), masks[i]]
    return rv


def parseSampleSheet(ss, fullSheets=False):
    """
    Return a dictionary with keys: (Barcode length 1, Barcode length 2)

    return ss, laneOut, bcLens
    """
    rv = dict()

    # If this is a NextSeq or a HiSeq 2500 rapid run, then don't store the incorrect Lane column
    storeLanes = True
    if getNumLanes(os.path.dirname(ss)) < 8:
        storeLanes = False

    f = open(ss)
    inData = False
    lastLine = None
    colLabs = [None, None, None, None] # Lane, Sample_ID, Sample_Name, Sample_Project
    indexCols = [None, None] # index, index2
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

                # Append to rv, with header
                if bcLen not in rv:
                    rv[bcLen] = [["[Data]", formatHeaderLine(cols, colLabs, indexCols, storeLanes)], set(), bcLen]

                # Add the lane to the set, if relevant
                if colLabs[0] is not None and storeLanes is True:
                    rv[bcLen][1].add(int(cols[colLabs[0]]))

                rv[bcLen][0].append(formatLine(cols, colLabs, indexCols, storeLanes))

            # Set columns for barcodes, etc.
            if lastLine is None:
                lastLine = True
                if "index" in cols:
                    indexCols[0] = cols.index("index")
                    if "index2" in cols:
                        indexCols[1] = cols.index("index2")

                if "Lane" in cols:
                    colLabs[0] = cols.index("Lane")

                if "Sample_ID" in cols:
                    colLabs[1] = cols.index("Sample_ID")

                if "Sample_Name" in cols:
                    colLabs[2] = cols.index("Sample_Name")

                if "Sample_Project" in cols:
                    colLabs[3] = cols.index("Sample_Project")
                continue

    if fullSheets:
        return reformatSS(handleRevComp(rv, os.path.dirname(ss)))
    else:
        return reformatSS(rv)


def getSampleSheets(d, fullSheets=False):
    """
    Provide a list of output directories and sample sheets
    """
    ss = glob.glob("%s/SampleSheet*.csv" % d)

    if len(ss) == 0:
        return ([None], [None], [''])

    laneOut = []
    bcLens = []
    ssUse = []
    for sheet in ss:
        ss_, laneOut_, bcLens_ = parseSampleSheet(sheet, fullSheets=False)
        nSS = 0
        if ss_ is not None and len(ss_) > 0:
            ssUse.extend(ss_)
            nSS = len(ss_)
        if nSS > 0 and laneOut_ is not None and len(laneOut_) > 0:
            laneOut.extend(laneOut_)
        elif nSS > 0:
            laneOut.extend([None] * nSS)
        if nSS > 0 and bcLens_ is not None and len(bcLens_) > 0:
            bcLens.extend(bcLens_)
        elif nSS > 0:
            bcLens.extend([None] * nSS)

    return ssUse, laneOut, bcLens


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
        # Before 1706, parkour wasn't being used, so barcode revComp might be wrong
        if config.get("Options","runID")[:4] < "1706":
            continue

        gotHits = False
        sampleSheet, lanes, bcLens = getSampleSheets(os.path.dirname(d))

        for ss, lane, bcLen in zip(sampleSheet, lanes, bcLens):
            config.set('Options','runID',d.split("/")[-2])
            lanesUse = ""
            if lane is not None and lane != "":
                config.set("Options","lanes",lane)
                lanesUse = "_lanes{}".format(lane)
            else:
                config.set("Options","lanes","")
            if ss is None:
                ss = ''
            if bcLen is not None and bcLen is not '':
                config.set("Options","bcLen",bcLen)
            else:
                config.set("Options","bcLen","0,0")
    
            if flowCellProcessed(config) is False:
                gotHits = True
            else :
                config.set("Options","runID","")

        # This may seem like code duplication, but for things like a MiSeq it takes a long time to parse the BCL files. This skips that unless needed
        if gotHits:
            sampleSheet, lanes, bcLens = getSampleSheets(os.path.dirname(d), fullSheet=True)
            for ss, lane, bcLen in zip(sampleSheet, lanes, bcLens):
                config.set('Options','runID',d.split("/")[-2])
                lanesUse = ""
                if lane is not None and lane != "":
                    config.set("Options","lanes",lane)
                    lanesUse = "_lanes{}".format(lane)
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
                    odir = "{}/{}{}".format(config.get("Paths", "outputDir"), config.get("Options", "runID"), lanesUse)
                    if not os.path.exists(odir):
                        os.makedirs(odir)
                    if ss is not None and not os.path.exists("{}/SampleSheet.csv".format(odir)):
                        o = open("{}/SampleSheet.csv".format(odir), "w")
                        o.write(ss)
                        o.close()
                        ss = "{}/SampleSheet.csv".format(odir)
                    config.set("Options","sampleSheet",ss)
                    return config

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
