'''
This file contains functions required to actually convert the bcl files to fastq
'''
import subprocess
import os
import shutil

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
    os.makedirs("%s/%s/InterOp" % (config.get("Paths","outputDir"),config.get("Options","runID")), exist_ok=True)
    cmd = "%s %s -o %s/%s -R %s/%s --interop-dir %s/%s/InterOp" % (
        config.get("bcl2fastq","bcl2fastq"),
        config.get("bcl2fastq","bcl2fastq_options"),
        config.get("Paths","outputDir"),
        config.get("Options","runID"),
        config.get("Paths","baseDir"),
        config.get("Options","runID"),
        config.get("Paths","outputDir"),
        config.get("Options","runID")
    )
    print("[bcl2fq] Running: %s" % cmd)
    subprocess.check_call(cmd, shell=True)

def cpXmlInterOp(config) :
    '''
    Copy over the InterOp directory and the base .xml files
    '''
    #Remove the tree if it exists
    shutil.rmtree("%s/%s/InterOp" % (config.get("Paths","outputDir"),config.get("Options","runID")), ignore_errors=True)
    shutil.copytree("%s/%s/InterOp" % (config.get("Paths","baseDir"),config.get("Options","runID")), "%s/%s/InterOp" % (config.get("Paths","outputDir"),config.get("Options","runID")))
    shutil.copy2("%s/%s/RunInfo.xml" % (config.get("Paths","baseDir"), config.get("Options","runID")), "%s/%s/" % (config.get("Paths","outputDir"), config.get("Options","runID")))
    shutil.copy2("%s/%s/runParameters.xml" % (config.get("Paths","baseDir"), config.get("Options","runID")), "%s/%s/" % (config.get("Paths","outputDir"), config.get("Options","runID")))
