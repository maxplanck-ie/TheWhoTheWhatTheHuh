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
import syslog
import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as img

'''
Do we really need the md5sum?
'''

localConfig = None

def bgzip_worker(fname) :
    global localConfig
    config = localConfig
    cmd = "%s -r %s" % (
        config.get("bgzip","bgzip_command"),
        fname)
    syslog.syslog("[bgzip_worker] Running %s\n" % cmd)
    subprocess.check_call(cmd, shell=True)

def plotFastqScreen(fname) :
    species=[]
    ohol=[]
    mhol=[]
    ohml=[]
    mhml=[]
    for line in csv.reader(open(fname, "r"), dialect="excel-tab") :
        if(len(line) == 0) :
            break
        if(line[0].startswith("#")) :
            continue
        if(line[0].startswith("Library")) :
            continue
        species.append(line[0])
        ohol.append(float(line[5]))
        mhol.append(float(line[7]))
        ohml.append(float(line[9]))
        mhml.append(float(line[11]))

    ohol = np.array(ohol)
    mhol = np.array(mhol)
    ohml = np.array(ohml)
    mhml = np.array(mhml)

    ind = np.arange(len(species))
    p1 = plt.bar(ind, tuple(ohol), color="#0000FF")
    p2 = plt.bar(ind, tuple(mhol), color="#6699FF", bottom=tuple(ohol))
    p3 = plt.bar(ind, tuple(ohml), color="#FF0000", bottom=tuple(ohol+mhol))
    p4 = plt.bar(ind, tuple(mhml), color="#FF6699", bottom=tuple(ohol+mhol+ohml))

    plt.title("%s" % fname.replace("_R1_screen.txt","").split("/")[-1])
    plt.ylabel("%")
    plt.ylim((0,105))
    plt.xticks(ind+0.4, species, rotation="vertical")
    plt.yticks(np.arange(0,110,10))
    plt.legend((p4[0], p3[0], p2[0], p1[0]), ("repeat", "conserved", "multimap", "unique"))
    plt.tight_layout()
    plt.savefig("%s.png" % fname.replace("_screen.txt","_screen"))
    plt.close()

def fastq_screen_worker(fname) :
    global localConfig
    config = localConfig

    #Skip read #2
    bname = fname.split("/")[-1]
    if(bname[-12:] == "_R2.fastq.gz") :
        return

    #Subsample
    ofile=fname.replace("_R1.fastq.gz","subsampled.fastq")
    cmd = "%s sample %s %s %s" % (
        config.get("fastq_screen","seqtk_command"),
        config.get("fastq_screen","seqtk_options"),
        fname,
        config.get("fastq_screen","seqtk_size"))
    syslog.syslog("[fastq_screen_worker] Running %s\n" % cmd)
    o = open(ofile, "w")
    subprocess.check_call(cmd, shell=True, stdout=o)
    o.close()

    #fastq_screen
    cmd = "%s %s %s" % (
        config.get("fastq_screen", "fastq_screen_command"),
        config.get("fastq_screen", "fastq_screen_options"),
        ofile)
    syslog.syslog("[fastq_screen_worker] Running %s\n" % cmd)
    subprocess.check_call(cmd, shell=True)

    #Unlink/rename
    os.unlink(ofile)
    os.rename(ofile.replace(".fastq","_screen.txt"), fname.replace("_R1.fastq.gz", "_R1_screen.txt"))

    #Create the images
    plotFastqScreen(fname.replace("_R1.fastq.gz", "_R1_screen.txt"))

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
    syslog.syslog("[FastQC_worker] Running %s\n" % cmd)
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
    syslog.syslog("[md5sum_worker] Processing %s\n" % d)
    subprocess.check_call(cmd, shell=True)
    os.chdir(oldWd)

def multiqc_worker(d) :
    global localConfig
    config = localConfig
    oldWd = os.getcwd()
    os.chdir(d)
    dname = d.split("/")
    dname[-1] = "FASTQC_{}".format(dname[-1])
    dname = "/".join(dname)
    cmd = "{} {} {}/*/*.zip".format(config.get("MultiQC", "multiqc_command"), config.get("MultiQC", "multiqc_options"), dname)
    syslog.syslog("[multiqc_worker] Processing %s\n" % d)
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

    projectDirs = glob.glob("%s/%s/Project_[ABC][0-9]*/*/*.fastq.gz" % (config.get("Paths","outputDir"), config.get("Options","runID")))
    projectDirs = toDirs(projectDirs)
    sampleFiles = glob.glob("%s/%s/Project_[ABC][0-9]*/*/*.fastq.gz" % (config.get("Paths","outputDir"),config.get("Options","runID")))
    global localConfig
    localConfig = config

    #FastQC
    p = mp.Pool(int(config.get("Options","postMakeThreads")))
    p.map(FastQC_worker, sampleFiles)
    p.close()
    p.join()

    #md5sum
    p = mp.Pool(int(config.get("Options","postMakeThreads")))
    p.map(md5sum_worker, projectDirs)
    p.close()
    p.join()

    #fastq_screen
    p = mp.Pool(int(config.get("Options", "postMakeThreads")))
    p.map(fastq_screen_worker, sampleFiles)
    p.close()
    p.join()

    # multiqc
    p = mp.Pool(int(config.get("Options","postMakeThreads")))
    p.map(multiqc_worker, projectDirs)
    p.close()
    p.join()

    #disk usage
    (tot,used,free) = shutil.disk_usage(config.get("Paths","outputDir"))
    tot /= 1024*1024*1024 #Convert to gigs
    used /= 1024*1024*1024
    free /= 1024*1024*1024

    #Undetermined indices
    undeter = parserDemultiplexStats(config)

    message = "Current free space: %i of %i gigs (%5.2f%%)\n" % (
        free,tot,100*free/tot)
    message += undeter
    return(message)
