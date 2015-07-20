'''
This file contains functions required to actually convert the bcl files to fastq
'''
import subprocess
import os
import sys
import shutil
import glob
import csv
import glob
import syslog
import codecs

def mergeLanesRename(config) :
    '''
    Merge samples that were sequenced on multiple lanes.

    Also, rename files such that they use sample names rather than IDs
    '''
    syslog.syslog("[mergeLanesRename] Merging and renaming samples\n")
    inLane=False
    for line in csv.reader(codecs.open("%s/%s/SampleSheet.csv" % (config.get("Paths","baseDir"), config.get("Options","runID")),"r","iso-8859-1")) :
        if(inLane) :
            #Read1
            if(len(line) == 7) :
                #Single sample, no demultiplexing needed
                files = glob.glob("%s/%s/Project_%s/Sample_%s/%s_NoIndex_L???_R1_???.fastq.gz" % (
                     config.get("Paths","outputDir"),
                     config.get("Options","runID"),
                     line[5],
                     line[1],
                     line[1]))
            else :
                files = glob.glob("%s/%s/Project_%s/Sample_%s/%s_%s_L???_R1_???.fastq.gz" % (
                     config.get("Paths","outputDir"),
                     config.get("Options","runID"),
                     line[7],
                     line[1],
                     line[1],
                     line[6]))
            if(len(files) == 0) :
                continue
            cmd = "cat "
            cmd2 = "rm "
            for f in files :
                cmd += "%s " % f
                cmd2 += "%s " % f
            if(len(line) == 7) :
                cmd += "> %s/%s/Project_%s/Sample_%s/%s_R1.fastq.gz" % (
                     config.get("Paths","outputDir"),
                     config.get("Options","runID"),
                     line[5],
                     line[1],
                     line[1])
            else :
                cmd += "> %s/%s/Project_%s/Sample_%s/%s_R1.fastq.gz" % (
                     config.get("Paths","outputDir"),
                     config.get("Options","runID"),
                     line[7],
                     line[1],
                     line[2])
            subprocess.check_call(cmd, shell=True)
            subprocess.check_call(cmd2, shell=True)
            #Read2
            if(len(line) == 7) :
                files = glob.glob("%s/%s/Project_%s/Sample_%s/%s_NoIndex_L???_R2_???.fastq.gz" % (
                     config.get("Paths","outputDir"),
                     config.get("Options","runID"),
                     line[5],
                     line[1],
                     line[1]))
            else :
                files = glob.glob("%s/%s/Project_%s/Sample_%s/%s_%s_L???_R2_???.fastq.gz" % (
                     config.get("Paths","outputDir"),
                     config.get("Options","runID"),
                     line[7],
                     line[1],
                     line[1],
                     line[6]))
            if(len(files) == 0) :
                continue
            cmd = "cat "
            cmd2 = "rm "
            for f in files :
                cmd += "%s " % f
                cmd2 += "%s " % f
            if(len(line) == 7) :
                files = glob.glob("%s/%s/Project_%s/Sample_%s/%s_NoIndex_L???_R2_???.fastq.gz" % (
                     config.get("Paths","outputDir"),
                     config.get("Options","runID"),
                     line[5],
                     line[1],
                     line[1]))
            else :
                files = glob.glob("%s/%s/Project_%s/Sample_%s/%s_%s_L???_R2_???.fastq.gz" % (
                     config.get("Paths","outputDir"),
                     config.get("Options","runID"),
                     line[7],
                     line[1],
                     line[1],
                     line[6]))
            if(len(line) == 7) :
                cmd += "> %s/%s/Project_%s/Sample_%s/%s_R1.fastq.gz" % (
                     config.get("Paths","outputDir"),
                     config.get("Options","runID"),
                     line[5],
                     line[1],
                     line[1])
            else :
                cmd += "> %s/%s/Project_%s/Sample_%s/%s_R1.fastq.gz" % (
                     config.get("Paths","outputDir"),
                     config.get("Options","runID"),
                     line[7],
                     line[1],
                     line[2])
            subprocess.check_call(cmd, shell=True)
            subprocess.check_call(cmd2, shell=True)
        else :
            if(len(line) == 0) :
                continue
            if(line[0] == "Lane") :
                inLane=True

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
    for line in csv.reader(codecs.open("%s/%s/SampleSheet.csv" % (config.get("Paths","baseDir"), config.get("Options","runID")),"r","iso-8859-1")) :
        if(inLane) :
            if(len(line) > 8) :
                newSS.write("%s,%s,%s,,%s,%s,N,,,%s\n" % (
                  FCID,
                  line[0],
                  line[1],
                  line[6],
                  line[8],
                  line[7]
                ))
            elif(len(line) == 7) :
                #Single sample format
                newSS.write("%s,%s,%s,,,%s,N,,,%s\n" % (
                  FCID,
                  line[0],
                  line[1],
                  line[6], 
                  line[5]
                ))
        else :
            if(len(line) == 0) :
                continue
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
    reformatSampleSheet(config)
    logOut = open("%s/%s.stdout" % (config.get("Paths","logDir"), config.get("Options","runID")), "w")
    logErr = open("%s/%s.stderr" % (config.get("Paths","logDir"), config.get("Options","runID")), "w")
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
    syslog.syslog("[bcl2fq] Running: %s\n" % cmd)
    subprocess.check_call(cmd, stdout=logOut, stderr=logErr, shell=True)
    cmd = "cd %s/%s && make -j %s" % (
        config.get("Paths","outputDir"),
        config.get("Options","runID"),
        config.get("bcl2fastq","make_threads")
    )
    syslog.syslog("[bcl2fq] Running: %s\n" % cmd)
    subprocess.check_call(cmd, stdout=logOut, stderr=logErr, shell=True)
    logOut.close()
    logErr.close()

    #This version of Illumina's software always splits by lane and names according to sampleID
    mergeLanesRename(config)

def cpSeqFac(config) :
    '''
    Copy over Xml and FastQC files
    '''
    #Ensure the directory doesn't exist
    shutil.rmtree("%s/%s" % (config.get("Paths","seqFacDir"),config.get("Options","runID")), ignore_errors=True)
    #InterOp
    shutil.copytree("%s/%s/InterOp" % (config.get("Paths","baseDir"), config.get("Options","runID")), "%s/%s/InterOp" % (config.get("Paths","seqFacDir"), config.get("Options","runID")))
    #Xml
    shutil.copy2("%s/%s/RunInfo.xml" % (config.get("Paths","baseDir"), config.get("Options","runID")), "%s/%s/" % (config.get("Paths","seqFacDir"), config.get("Options","runID")))
    shutil.copy2("%s/%s/runParameters.xml" % (config.get("Paths","baseDir"), config.get("Options","runID")), "%s/%s/" % (config.get("Paths","seqFacDir"), config.get("Options","runID")))

    #FastQC
    dirs = glob.glob("%s/%s/FASTQC_*" % (config.get("Paths","outputDir"), config.get("Options","runID")))
    for d in dirs :
        dname = d.split("/")[-1]
        shutil.rmtree("%s/%s/%s" % (config.get("Paths","seqFacDir"), config.get("Options","runID"), dname), ignore_errors=True)
        shutil.copytree(d, "%s/%s/%s" % (config.get("Paths","seqFacDir"), config.get("Options","runID"), dname))
