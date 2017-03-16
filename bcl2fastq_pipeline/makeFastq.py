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
import re
from reportlab.lib import colors, utils
from reportlab.platypus import BaseDocTemplate, Table, Preformatted, Paragraph, Spacer, Image, Frame, NextPageTemplate, PageTemplate, TableStyle
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas

def determineMask(config):
    '''
    If there's already a mask set in the config file then return it.

    Otherwise:
     1. Check for RunInfo.xml
     2. Parse each <Read> child, adding it to a list.
     3. Join the list by commas
     4. If there's no mask then return nothing
    '''
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes2 = []
        for l in lanes.split("_"):
            lanes2.append("s_{}".format(l))
        lanes = "--tiles {}".format(",".join(lanes2))

    mask = config.get("Options", "index_mask")
    bcLens = [int(x) for x in config.get("Options","bcLen").split(",")]
    bcNum = 0
    if mask != "":
        return "--use-bases-mask {} {}".format(mask, lanes)
    elif os.path.isfile("{}/{}/RunInfo.xml".format(config.get("Paths","baseDir"),config.get("Options","runID"))):
        xml = ET.parse("{}/{}/RunInfo.xml".format(config.get("Paths","baseDir"),config.get("Options","runID")))
        root = xml.getroot()[0][3]
        l = []
        for read in root.findall("Read"):
            if read.get("IsIndexedRead") == "N":
                l.append("Y*")
            else:
                nc = int(read.get("NumCycles"))
                if nc > bcLens[bcNum]:
                    if bcLens[bcNum] > 0:
                        l.append("I{}{}".format(bcLens[bcNum], "n" * (nc - bcLens[bcNum])))
                    else:
                        l.append("{}".format("n" * nc))
                else:
                    l.append("I{}".format(bcLens[bcNum]))
                bcNum += 1
        if len(l) > 0:
            return "--use-bases-mask {} {}".format(",".join(l), lanes)
    return lanes

def rewriteSampleSheet(config) :
    '''
    If it exists, make a modified copy of the sample sheet to ensure that:
     A) It contains no special characters
     B) It contains no adapter sequences (so adapter trimming is disabled
    An appropriate "--sample-sheet blah --use-bases-mask foo" is returned
    '''

    ssheet = config.get("Options", "sampleSheet")
    if ssheet is None or ssheet == "":
        ssheet = "%s/%s/SampleSheet.csv" % (config.get("Paths", "baseDir"), config.get("Options", "runID"))

    if os.path.isfile(ssheet):
        od, oname = tempfile.mkstemp()
        config.set("Options", "sampleSheet", oname)
        of = open(oname, "w")
        inData = False
        inReads = 0
        for line in codecs.open(ssheet, "r", "iso-8859-1"):
            if((line.startswith("Lane") or line.startswith("Sample_ID")) and (inData is False)) :
                inData = True
            elif(line.startswith("[Reads]")) :
                inReads = 1
            elif(inReads == 1) :
                inReads = 2
            elif(inReads==2) :
                inReads = 0
            elif(inData) :
                #. to _dot_
                line = line.replace(".", "_dot_")
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
                #ß to sz
                line = line.replace("ß", "sz")
                #& to _and_
                line = line.replace("&", "_and_")
            else :
                if(line.startswith("Adapter")) :
                    continue
            of.write(line)
        of.close()
        os.close(od)
        return "--sample-sheet {} {}".format(oname, determineMask(config))
    else :
        config.set("Options", "sampleSheet", "")
        return None

def fixNames(config) :
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes = "_lanes{}".format(lanes)

    fnames = glob.glob("%s/%s%s/[ABC][0-9]*/*/*.fastq.gz" % (config.get("Paths","outputDir"), config.get("Options","runID"), lanes))
    for fname in fnames :
        idx = fname.rindex("_")
        fnew = fname[0:idx]
        fnew = re.sub(r"_S[0-9]+_R([12])$",r'_R\1', fnew) + ".fastq.gz"
        syslog.syslog("Moving %s to %s\n" % (fname, fnew))
        shutil.move(fname, fnew)

    snames = glob.glob("%s/%s%s/[ABC][0-9]*/*" % (config.get("Paths","outputDir"), config.get("Options","runID"), lanes))
    for sname in snames :
        idx = sname.rindex("/")
        snew = "%s/Sample_%s" % (sname[:idx], sname[idx+1:])
        syslog.syslog("Moving %s to %s\n" % (sname, snew))
        shutil.move(sname, snew)

    pnames = glob.glob("%s/%s%s/[ABC][0-9]*" % (config.get("Paths","outputDir"), config.get("Options","runID"), lanes))
    for pname in pnames :
        idx = pname.rindex("/")
        pnew = "%s/Project_%s" % (pname[:idx], pname[idx+1:])
        syslog.syslog("Moving %s to %s\n" % (pname, pnew))
        shutil.move(pname, pnew)

def bcl2fq(config) :
    '''
    takes things from /dont_touch_this/solexa_runs/XXX/Data/Intensities/BaseCalls
    and writes most output into config.outputDir/XXX, where XXX is the run ID.
    '''
    lanes = config.get("Options", "lanes")
    if lanes != '':
        lanes = '_lanes{}'.format(lanes)
        
    #Make the output directories
    os.makedirs("%s/%s%s" % (config.get("Paths","outputDir"), config.get("Options","runID"), lanes), exist_ok=True)
    os.makedirs("%s/%s%s/InterOp" % (config.get("Paths","seqFacDir"), config.get("Options","runID"), lanes), exist_ok=True)

    #If there's no sample sheet then we need to not mask the last index base!
    rv = rewriteSampleSheet(config)
    mask = ""
    if(rv is not None) :
        mask = rv
    cmd = "%s %s %s -o %s/%s%s -R %s/%s --interop-dir %s/%s%s/InterOp" % (
        config.get("bcl2fastq","bcl2fastq"),
        config.get("bcl2fastq","bcl2fastq_options"),
        mask,
        config.get("Paths","outputDir"),
        config.get("Options","runID"),
        lanes,
        config.get("Paths","baseDir"),
        config.get("Options","runID"),
        config.get("Paths","seqFacDir"),
        config.get("Options","runID"),
        lanes
    )
    syslog.syslog("[bcl2fq] Running: %s\n" % cmd)
    logOut = open("%s/%s%s.stdout" % (config.get("Paths","logDir"), config.get("Options","runID"), lanes), "w")
    logErr = open("%s/%s%s.stderr" % (config.get("Paths","logDir"), config.get("Options","runID"), lanes), "w")
    subprocess.check_call(cmd, stdout=logOut, stderr=logErr, shell=True)
    logOut.close()
    logErr.close()

def getOffSpecies(fname) :
    total = 0
    species=[]
    ohol=[]
    mhol=[]
    i = 0
    maxi = 0
    for line in csv.reader(open(fname, "r"), dialect="excel-tab") :
        if(len(line) == 0) :
            break
        if(line[0].startswith("#")) :
            continue
        if(line[0].startswith("Library")) :
            continue
        if(line[0].startswith("PhiX") or line[0].startswith("Adapters") or line[0].startswith("Vectors") or line[0].startswith("rRNA")):
            continue
        species.append(line[0])
        ohol.append(float(line[5]))

        if(ohol[maxi] < ohol[i]) :
            maxi = i
        i += 1

    off = 0
    for i in range(len(ohol)) :
        if(i != maxi) :
            off += ohol[i]
    return off

def MakeTotalPDF(config) :
    '''
    Make a PDF containing the fastq_screen images from each sample

    Also, parse the fastq_screen .txt files and calculate the per-sample off-species rate
    '''
    lanes = config.get("Options", "lanes")
    if lanes != '':
        lanes = '_lanes{}'.format(lanes)

    stylesheet=getSampleStyleSheet()

    pdf = BaseDocTemplate("%s/%s%s/ContaminationReport.pdf" % (
        config.get("Paths","outputDir"),
        config.get("Options","runID"),
        lanes), pagesize=A4)
    fM = Frame(pdf.leftMargin, pdf.bottomMargin, pdf.width, pdf.height, id="main")

    tab = [["Project", "Sample", "confident off-species reads/sample", "% Optical Duplication"]]
    txt = "\nProject\tSample\tconfident off-species reads/sample\t% Optical Duplicates\n"
    elements = []

    projs = glob.glob("%s/%s%s/Project_*" % (
        config.get("Paths","outputDir"),
        config.get("Options","runID"),
        lanes))
    projs.sort()

    #Make the table
    for proj in projs :
        pname=proj.split("/")[-1]
        txts = glob.glob("%s/Sample_*/*_screen.txt" % proj)
        txts.sort()
        for i in range(len(txts)) :
            # Get the percent optical duplication
            _ = "{}.duplicate.txt".format(txts[i][:-14])
            if os.path.exists(_):
                _ = open(_).read().strip().split("\t")
                duplicationRate = "{:5.2f}%".format((100.0 * int(_[0])) / int(_[1]))
            else:
                duplicationRate = "NA"
            tab.append([pname, txts[i].split("/")[-2], "%5.2f" % getOffSpecies(txts[i]), duplicationRate])
            txt += "%s\t%s\t%5.2f\t%s\n" % (pname, txts[i].split("/")[-2], getOffSpecies(txts[i]), duplicationRate)
    txt += "\n"

    #Add the table
    t = Table(tab, style=[
        ('ROWBACKGROUNDS', (0, 0), (-1, -1), (0xD3D3D3, None)) #Light grey
        ], repeatRows=1)
    elements.append(t)
    elements.append(Spacer(0,30))

    #Add the images
    for proj in projs :
        pname=proj.split("/")[-1]
        elements.append(Paragraph(pname, stylesheet['title']))
        imgs = glob.glob("%s/Sample_*/*.png" % proj)
        imgs.sort()
        for i in range(len(imgs)) :
            TmpImg = utils.ImageReader(imgs[i])
            iw, ih = TmpImg.getSize()
            iw = 0.5*iw
            ih = 0.5*ih
            elements.append(Image(imgs[i], width=iw, height=ih, hAlign="LEFT"))

    pdf.addPageTemplates([PageTemplate(id="foo", frames=[fM])])
    pdf.build(elements)

    return txt

def cpSeqFac(config) :
    '''
    Copy over Xml and FastQC files
    '''
    lanes = config.get("Options", "lanes")
    if lanes != '':
        lanes = '_lanes{}'.format(lanes)

    shutil.rmtree("%s/%s%s" % (config.get("Paths","seqFacDir"),config.get("Options","runID"), lanes), ignore_errors=True)
    shutil.copytree("%s/%s/InterOp" % (config.get("Paths","baseDir"), config.get("Options","runID")), "%s/%s%s/InterOp" % (config.get("Paths","seqFacDir"), config.get("Options","runID"), lanes))
    #Xml
    shutil.copy2("%s/%s/RunInfo.xml" % (config.get("Paths","baseDir"), config.get("Options","runID")), "%s/%s%s/" % (config.get("Paths","seqFacDir"), config.get("Options","runID"), lanes))
    try:
        shutil.copy2("%s/%s/runParameters.xml" % (config.get("Paths","baseDir"), config.get("Options","runID")), "%s/%s%s/" % (config.get("Paths","seqFacDir"), config.get("Options","runID"), lanes))
    except:
        #renamed on a Nextseq
        shutil.copy2("%s/%s/RunParameters.xml" % (config.get("Paths","baseDir"), config.get("Options","runID")), "%s/%s%s/" % (config.get("Paths","seqFacDir"), config.get("Options","runID"), lanes))

    #Make the PDF
    txt = MakeTotalPDF(config)
    shutil.copy2("%s/%s%s/ContaminationReport.pdf" % (config.get("Paths","outputDir"), config.get("Options","runID"), lanes), "%s/%s%s/" % (config.get("Paths","seqFacDir"), config.get("Options","runID"), lanes))

    #FastQC
    dirs = glob.glob("%s/%s%s/FASTQC_*" % (config.get("Paths","outputDir"), config.get("Options","runID"), lanes))
    for d in dirs :
        dname = d.split("/")[-1]
        shutil.rmtree("%s/%s%s/%s" % (config.get("Paths","seqFacDir"), config.get("Options","runID"), lanes, dname), ignore_errors=True)
        shutil.copytree(d, "%s/%s%s/%s" % (config.get("Paths","seqFacDir"), config.get("Options","runID"), lanes, dname))

    return txt
