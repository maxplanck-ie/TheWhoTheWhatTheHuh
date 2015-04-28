'''
This file includes code that actually runs FastQC and any other tools after the fastq files have actually been made. This uses a pool of workers to process each request.
'''
import multiprocessing as mp
import glob
import sys
import subprocess
import os
import shutil
import xml.etree.ElementTree as ET

'''
Do we really need the md5sum?
'''

localConfig = None

def FastQC_worker(fname) :
    global localConfig
    config = localConfig
    projectName = fname.split("/")[-3] #It's the penultimate directory
    #--noextract -f fastq -t 1 -q
    cmd = "%s %s -o %s/%s/FASTQC_%s %s" % (
          config.get("FastQC","fastqc_command"),
          config.get("FastQC","fastqc_options"),
          config.get("Paths","outputDir"),
          config.get("Options","runID"),
          projectName,
          fname)
    os.makedirs("%s/%s/FASTQC_%s" % (config.get("Paths","outputDir"),
          config.get("Options","runID"),
          projectName), exist_ok=True)
    print("[FastQC_worker] Running %s" % cmd)
    subprocess.check_call(cmd, shell=True)

def toDirs(files) :
    s = set()
    for f in files :
        d = os.path.dirname(f)
        s.add(d[:d.rfind('/')]) #We just want projects, not individual libraries
    return s

def md5sum_worker(d) :
    global localConfig
    config = localConfig
    oldWd = os.getcwd()
    os.chdir(d)
    cmd = "md5sum */*.fastq.gz > md5sums.txt"
    print("[md5sum_worker] Processing %s" % d)
    subprocess.check_call(cmd, shell=True)
    os.chdir(oldWd)

def parserDemultiplexStats(config) :
    '''
    Parse DemultiplexingStats.xml under outputDir/Stats/ to get the
    number/percent of undetermined indices.

    In particular, we extract the BarcodeCount values from Project "default"
    Sample "all" and Project "all" Sample "all", as the former gives the total
    undetermined and the later simply the total clusters
    '''
    totals = [0,0,0,0,0,0,0,0]
    undetermined = [0,0,0,0,0,0,0,0]
    tree = ET.parse("%s/%s/Stats/DemultiplexingStats.xml" % (config.get("Paths","outputDir"),config.get("Options","runID")))
    root = tree.getroot()
    for child in root[0].findall("Project") :
        if(child.get("name") == "default") :
            break
    for sample in child.findall("Sample") :
        if(sample.get("name") == "all") :
            break
    child = sample[0] #Get inside Barcode
    for lane in child.findall("Lane") :
        lnum = int(lane.get("number"))
        undetermined[lnum-1] += int(lane[0].text)

    for child in root[0].findall("Project") :
        if(child.get("name") == "all") :
            break
    for sample in child.findall("Sample") :
        if(sample.get("name") == "all") :
            break
    child = sample[0] #Get Inside Barcode
    for lane in child.findall("Lane") :
        lnum = int(lane.get("number"))
        totals[lnum-1] += int(lane[0].text)

    out = ""
    for i in range(8) :
        if(totals[i] == 0) :
            continue
        out += "\nLane %i: %i of %i reads/pairs had undetermined indices (%5.2f%%)" % (
            i+1,undetermined[i],totals[i],100*undetermined[i]/totals[i])
    return out

#All steps that should be run after `make` go here
def postMakeSteps(config) :
    '''
    Current steps are:
      1) Run FastQC on each fastq.gz file
      2) Run md5sum on the files in each project directory
    Other steps could easily be added to follow those. Note that this function
    will try to use a pool of threads. The size of the pool is set by config.postMakeThreads
    '''

    projectDirs = glob.glob("%s/%s/*/*/*.fastq.gz" % (config.get("Paths","outputDir"), config.get("Options","runID")))
    projectDirs = toDirs(projectDirs)
    sampleFiles = glob.glob("%s/%s/*/*/*.fastq.gz" % (config.get("Paths","outputDir"),config.get("Options","runID")))
    global localConfig
    localConfig = config

    #FastQC
    p = mp.Pool(int(config.get("Options","postMakeThreads")))
    p.map(FastQC_worker, sampleFiles)

    #md5sum
    p = mp.Pool(int(config.get("Options","postMakeThreads")))
    p.map(md5sum_worker, projectDirs)

    #disk usage
    (tot,used,free) = shutil.disk_usage(config.get("Paths","outputDir"))
    tot /= 1024*1024*1024 #Convert to gigs
    used /= 1024*1024*1024
    free /= 1024*1024*1024

    #Undetermined indices
    undeter = parserDemultiplexStats(config)

    message = "\nCurrent free space: %i of %i gigs (%5.2f%%)" % (
        free,tot,100*free/tot)
    message += undeter
    return(message)
