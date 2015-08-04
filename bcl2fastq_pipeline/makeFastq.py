'''
This file contains functions required to actually convert the bcl files to fastq
'''
import subprocess
import os
import sys
import shutil
import glob
import syslog
import csv
import codecs
import tempfile
import xml.etree.ElementTree as ET

def rewriteSampleSheet(config) :
    '''
    If it exists, make a modified copy of the sample sheet to ensure that:
     A) It contains no special characters
     B) It contains no adapter sequences (so adapter trimming is disabled
    An appropriate "--sample-sheet blah --use-bases-mask foo" is returned
    '''

    if(os.path.isfile("%s/%s/SampleSheet.csv" % (
        config.get("Paths", "baseDir"),
        config.get("Options", "runID")
    ))) :
        od, oname = tempfile.mkstemp()
        config.set("Options", "sampleSheet", oname)
        of = open(oname, "w")
        inData = False
        for line in codecs.open("%s/%s/SampleSheet.csv" % (config.get("Paths","baseDir"),config.get("Options","runID")), "r", "iso-8859-1") :
            if(line.startswith("Lane")) :
                inData = True
            elif(inData) :
                #+ to _plus_
                line = line.replace("+", "_plus_")
                #ö to oe
                line = line.replace("ö", "oe")
                line = line.replace("Ö", "Oe")
                #ä to ae
                line = line.replace("ä", "ae")
                line = line.replace("Ä", "Ae")
                #ü to ue
                line = line.replace("ü", "ue")
                line = line.replace("Ü", "Ue")
            else :
                if(line.startswith("Adapter")) :
                    continue
            of.write(line)
        of.close()
        os.close(od)
        return "--sample-sheet %s" % oname
    else :
        config.set("Options", "sampleSheet", "")
        return None

def fixNames(config) :
    fnames = glob.glob("%s/%s/[ABC][0-9]*/*/*.fastq.gz" % (config.get("Paths","outputDir"), config.get("Options","runID")))
    for fname in fnames :
        idx = fname.rindex("_")
        fnew = fname[0:idx]+".fastq.gz"
        syslog.syslog("Moving %s to %s\n" % (fname, fnew))
        shutil.move(fname, fnew)

def bcl2fq(config) :
    '''
    takes things from /dont_touch_this/solexa_runs/XXX/Data/Intensities/BaseCalls
    and writes most output into config.outputDir/XXX, where XXX is the run ID.
    '''
    #Make the output directories
    os.makedirs("%s/%s" % (config.get("Paths","outputDir"),config.get("Options","runID")), exist_ok=True)
    os.makedirs("%s/%s/InterOp" % (config.get("Paths","seqFacDir"),config.get("Options","runID")), exist_ok=True)

    #If there's no sample sheet then we need to not mask the last index base!
    rv = rewriteSampleSheet(config)
    mask = ""
    if(rv is not None) :
        mask = rv
    cmd = "%s %s %s -o %s/%s -R %s/%s --interop-dir %s/%s/InterOp" % (
        config.get("bcl2fastq","bcl2fastq"),
        config.get("bcl2fastq","bcl2fastq_options"),
        mask,
        config.get("Paths","outputDir"),
        config.get("Options","runID"),
        config.get("Paths","baseDir"),
        config.get("Options","runID"),
        config.get("Paths","seqFacDir"),
        config.get("Options","runID")
    )
    syslog.syslog("[bcl2fq] Running: %s\n" % cmd)
    logOut = open("%s/%s.stdout" % (config.get("Paths","logDir"), config.get("Options","runID")), "w")
    logErr = open("%s/%s.stderr" % (config.get("Paths","logDir"), config.get("Options","runID")), "w")
    subprocess.check_call(cmd, stdout=logOut, stderr=logErr, shell=True)
    logOut.close()
    logErr.close()
    fixNames(config)

def cpSeqFac(config) :
    '''
    Copy over Xml and FastQC files
    '''
    shutil.rmtree("%s/%s" % (config.get("Paths","seqFacDir"),config.get("Options","runID")), ignore_errors=True)
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
