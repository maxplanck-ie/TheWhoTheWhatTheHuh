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
    libName = fname.split("/")[-2] #The last directory
    cmd = "%s %s -o %s/%s/FASTQC_%s/%s %s" % (
          config.get("FastQC","fastqc_command"),
          config.get("FastQC","fastqc_options"),
          config.get("Paths","outputDir"),
          config.get("Options","runID"),
          projectName,
          libName,
          fname)
    os.makedirs("%s/%s/FASTQC_%s/%s" % (config.get("Paths","outputDir"),
          config.get("Options","runID"),
          projectName,
          libName), exist_ok=True)
    sys.stderr.write("[FastQC_worker] Running %s\n" % cmd)
    sys.stderr.flush()
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
    sys.stderr.write("[md5sum_worker] Processing %s\n" % d)
    sys.stderr.flush()
    subprocess.check_call(cmd, shell=True)
    os.chdir(oldWd)

def parseFlowcell_demux_summary(config) :
    '''
    Parse Basecall_Stats_XXX/Flowcell_demux_summary.xml to get the
    number/percent of undetermined indices.

    In particular, we extract the BarcodeCount values from Project "default"
    Sample "all" and Project "all" Sample "all", as the former gives the total
    undetermined and the later simply the total clusters
    '''
    determined = [0,0,0,0,0,0,0,0]
    undetermined = [0,0,0,0,0,0,0,0]
    FCID = config.get("Options","runID").split("_")[-1][1:]
    tree = ET.parse("%s/%s/Basecall_Stats_%s/Flowcell_demux_summary.xml" % (config.get("Paths","outputDir"),config.get("Options","runID"), FCID))
    root = tree.getroot()
    for lane in root.findall("Lane") :
        lnum = int(lane.get("index"))-1
        for sample in lane.findall("Sample") :
            for barcode in sample.findall("Barcode") :
                if(barcode.get("index") == "Undetermined") :
                    for tile in barcode.findall("Tile") :
                        undetermined[lnum] += int(tile[0][1][2].text) #Read 1->Pf->ClusterCount
                else :
                    for tile in barcode.findall("Tile") :
                        determined[lnum] += int(tile[0][1][2].text) #Read 1->Pf->ClusterCount

    out = ""
    for i in range(8) :
        if(determined[i]+undetermined[i] == 0) :
            continue
        out += "\nLane %i: %i of %i reads/pairs had undetermined indices (%5.2f%%)" % (
            i+1,undetermined[i],determined[i]+undetermined[i],
            100*undetermined[i]/(determined[i]+undetermined[i]))
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

    projectDirs = glob.glob("%s/%s/Project_*/*/*.fastq.gz" % (config.get("Paths","outputDir"), config.get("Options","runID")))
    projectDirs = toDirs(projectDirs)
    sampleFiles = glob.glob("%s/%s/Project_*/*/*.fastq.gz" % (config.get("Paths","outputDir"),config.get("Options","runID")))
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
    undeter = parseFlowcell_demux_summary(config)

    message = "Current free space: %i of %i gigs (%5.2f%%)\n" % (
        free,tot,100*free/tot)
    message += undeter
    return(message)
