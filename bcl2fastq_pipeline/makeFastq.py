'''
This file contains functions required to actually convert the bcl files to fastq
'''
import subprocess
import os
import sys
import shutil
import glob
import csv

def reformatSampleSheet(config) :
    '''
    I have no idea why the sample sheet that's input into the machine doesn't match
    what's being expected by Illumina's software. It works with the newer bcl2fastq
    software, but not the older stuff.
    '''
    newSS = open("%s/%s/SampleSheet.csv" % (config.get("Paths","outputDir"),config.get("Options","runID")), "w")
    newSS.write("FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject\n")
    inLane=False
    FCID = config.get("Options","runID").split("_")[-1][1:]
    for line in csv.reader(open("%s/%s/SampleSheet.csv" % (config.get("Paths","baseDir"), config.get("Options","runID")),"r")) :
        if(inLane) :
            newSS.write("%s,%s,%s,,%s,%s,N,,,%s\n" % (
              FCID,
              line[0],
              line[1],
              line[6],
              line[8],
              line[7]
            ))
        else :
            if(line[0] == "Lane") :
                inLane=True
    newSS.close()

def bcl2fq(config) :
    '''
    takes things from /dont_touch_this/solexa_runs/XXX/Data/Intensities/BaseCalls
    and writes most output into outputDir/XXX, where XXX is the run ID.
    
    config.baseDir = /dont_touch_this/solexa_runs
    config.outputDir = /dont_touch_this/solexa_data/final
    config.runID = XXX
    '''
    #Make the output directories
    os.makedirs("%s/%s" % (config.get("Paths","outputDir"),config.get("Options","runID")), exist_ok=True)
    os.makedirs("%s/%s/InterOp" % (config.get("Paths","seqFacDir"),config.get("Options","runID")), exist_ok=True)
    reformatSampleSheet(config)
    cmd = "%s %s -o %s/%s --input-dir %s/%s/Data/Intensities/BaseCalls --sample-sheet %s/%s/SampleSheet.csv" % (
        config.get("bcl2fastq","configure"),
        config.get("bcl2fastq","configure_options"),
        config.get("Paths","outputDir"),
        config.get("Options","runID"),
        config.get("Paths","baseDir"),
        config.get("Options","runID"),
        config.get("Paths","outputDir"),
        config.get("Options","runID")
    )
    sys.stderr.write("[bcl2fq] Running: %s\n" % cmd)
    sys.stderr.flush()
    subprocess.check_call(cmd, shell=True)
    cmd = "cd %s/%s && make -j %s" % (
        config.get("Paths","outputDir"),
        config.get("Options","runID"),
        config.get("bcl2fastq","make_threads")
    )
    sys.stderr.write("[bcl2fq] Running: %s\n" % cmd)
    sys.stderr.flush()
    subprocess.check_call(cmd, shell=True)

def cpSeqFac(config) :
    '''
    Copy over Xml and FastQC files
    '''
    shutil.rmtree("%s/%s" % (config.get("Paths","seqFacDir"),config.get("Options","runID")), ignore_errors=True)
    os.makedirs("%s/%s/InterOp" % (config.get("Paths","seqFacDir"),config.get("Options","runID")), exist_ok=True)
    #Xml
    shutil.copy2("%s/%s/RunInfo.xml" % (config.get("Paths","baseDir"), config.get("Options","runID")), "%s/%s/" % (config.get("Paths","seqFacDir"), config.get("Options","runID")))
    shutil.copy2("%s/%s/runParameters.xml" % (config.get("Paths","baseDir"), config.get("Options","runID")), "%s/%s/" % (config.get("Paths","seqFacDir"), config.get("Options","runID")))

    #FastQC
    dirs = glob.glob("%s/%s/FASTQC_*" % (config.get("Paths","outputDir"), config.get("Options","runID")))
    for d in dirs :
        dname = d.split("/")[-1]
        shutil.rmtree("%s/%s/%s" % (config.get("Paths","seqFacDir"), config.get("Options","runID"), dname), ignore_errors=True)
        shutil.copytree(d, "%s/%s/%s" % (config.get("Paths","seqFacDir"), config.get("Options","runID"), dname))
