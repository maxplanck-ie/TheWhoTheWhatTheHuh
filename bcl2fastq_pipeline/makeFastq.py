'''
This file contains functions required to actually convert the bcl files to fastq
'''
import subprocess
import os
import sys
import shutil
import glob

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
    cmd = "%s %s -o %s/%s -R %s/%s --interop-dir %s/%s/InterOp" % (
        config.get("bcl2fastq","bcl2fastq"),
        config.get("bcl2fastq","bcl2fastq_options"),
        config.get("Paths","outputDir"),
        config.get("Options","runID"),
        config.get("Paths","baseDir"),
        config.get("Options","runID"),
        config.get("Paths","seqFacDir"),
        config.get("Options","runID")
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
