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
from reportlab.lib import colors, utils
from reportlab.platypus import BaseDocTemplate, Table, Preformatted, Paragraph, Spacer, Image, Frame, NextPageTemplate, PageTemplate, TableStyle
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas

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
        PE = False
        BC = True
        inReads = 0
        for line in codecs.open("%s/%s/SampleSheet.csv" % (config.get("Paths","baseDir"),config.get("Options","runID")), "r", "iso-8859-1") :
            if(line.startswith("Lane") and (inData is False)) :
                inData = True
                #Do we have a barcode?
                if(line.split(",")[6] != "index") :
                    BC = False
            elif(line.startswith("[Reads]")) :
                inReads = 1
            elif(inReads == 1) :
                inReads = 2
            elif(inReads==2) :
                if(line.startswith(",") is False) :
                    PE = True
                inReads = 0
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
                #ß to sz
                line = line.replace("ß", "sz")
            else :
                if(line.startswith("Adapter")) :
                    continue
            of.write(line)
        of.close()
        os.close(od)
        if(BC is False) :
            return "--sample-sheet %s" % oname
        if(PE is True) :
            return "--sample-sheet %s --use-bases-mask Y*,I6n,Y*" % oname
        else :
            return "--sample-sheet %s --use-bases-mask Y*,I6n" % oname
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

    snames = glob.glob("%s/%s/[ABC][0-9]*/*" % (config.get("Paths","outputDir"), config.get("Options","runID")))
    for sname in snames :
        idx = sname.rindex("/")
        snew = "%s/Sample_%s" % (sname[:idx], sname[idx+1:])
        syslog.syslog("Moving %s to %s\n" % (sname, snew))
        shutil.move(sname, snew)

    pnames = glob.glob("%s/%s/[ABC][0-9]*" % (config.get("Paths","outputDir"), config.get("Options","runID")))
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
        species.append(line[0])
        ohol.append(float(line[5]))
        mhol.append(float(line[7]))
        if(ohol[maxi]+mhol[maxi] < ohol[i]+mhol[i]) :
            maxi = i
        i += 1

    off = 0
    for i in range(len(ohol)) :
        if(i != maxi) :
            off += ohol[i] + mhol[i]
    return off

def MakeTotalPDF(config) :
    '''
    Make a PDF containing the fastq_screen images from each sample

    Also, parse the fastq_screen .txt files and calculate the per-sample off-species rate
    '''

    stylesheet=getSampleStyleSheet()

    pdf = BaseDocTemplate("%s/%s/ContaminationReport.pdf" % (
        config.get("Paths","outputDir"),
        config.get("Options","runID")), pagesize=A4)
    fM = Frame(pdf.leftMargin, pdf.bottomMargin, pdf.width, pdf.height, id="main")

    tab = [["Project", "Sample", "off-species reads/sample"]]
    txt = "\nProject\tSample\toff-species reads/sample\n"
    elements = []

    projs = glob.glob("%s/%s/Project_*" % (
        config.get("Paths","outputDir"),
        config.get("Options","runID")))
    projs.sort()

    #Make the table
    for proj in projs :
        pname=proj.split("/")[-1]
        txts = glob.glob("%s/Sample_*/*.txt" % proj)
        txts.sort()
        for i in range(len(txts)) :
            tab.append([pname, txts[i].split("/")[-2], "%5.2f" % getOffSpecies(txts[i])])
            txt += "%s\t%s\t%5.2f\n" % (pname, txts[i].split("/")[-2], getOffSpecies(txts[i]))
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
    shutil.rmtree("%s/%s" % (config.get("Paths","seqFacDir"),config.get("Options","runID")), ignore_errors=True)
    shutil.copytree("%s/%s/InterOp" % (config.get("Paths","baseDir"), config.get("Options","runID")), "%s/%s/InterOp" % (config.get("Paths","seqFacDir"), config.get("Options","runID")))
    #Xml
    shutil.copy2("%s/%s/RunInfo.xml" % (config.get("Paths","baseDir"), config.get("Options","runID")), "%s/%s/" % (config.get("Paths","seqFacDir"), config.get("Options","runID")))
    shutil.copy2("%s/%s/runParameters.xml" % (config.get("Paths","baseDir"), config.get("Options","runID")), "%s/%s/" % (config.get("Paths","seqFacDir"), config.get("Options","runID")))

    #Make the PDF
    txt = MakeTotalPDF(config)
    shutil.copy2("%s/%s/ContaminationReport.pdf" % (config.get("Paths","outputDir"), config.get("Options","runID")), "%s/%s/" % (config.get("Paths","seqFacDir"), config.get("Options","runID")))

    #FastQC
    dirs = glob.glob("%s/%s/FASTQC_*" % (config.get("Paths","outputDir"), config.get("Options","runID")))
    for d in dirs :
        dname = d.split("/")[-1]
        shutil.rmtree("%s/%s/%s" % (config.get("Paths","seqFacDir"), config.get("Options","runID"), dname), ignore_errors=True)
        shutil.copytree(d, "%s/%s/%s" % (config.get("Paths","seqFacDir"), config.get("Options","runID"), dname))

    return txt
