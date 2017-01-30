"""
Misc. functions
"""

import configparser
import shutil
import smtplib
from email.mime.text import MIMEText
import xml.etree.ElementTree as ET
from reportlab.lib import colors, utils
from reportlab.platypus import BaseDocTemplate, Table, Preformatted, Paragraph, Spacer, Image, Frame, NextPageTemplate, PageTemplate, TableStyle, PageBreak, ListFlowable, ListItem
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.pagesizes import A4, landscape
from time import strftime
from reportlab.pdfgen import canvas
import csv
import sys
import glob
import pathlib
import os
import os.path
import syslog
import codecs

def transferData(config) :
    """
    Distribute fastq and fastQC files to users.
    """
    message = ""
    projects = glob.glob("%s/%s/Project_*" % (config.get("Paths","outputDir"),config.get("Options","runID")))
    for project in projects :
        pname = project.split("/")[-1][8:]
        if(pname[0] == "A") :
            #Copy
            group = pname.split("_")[1].lower()
            syslog.syslog("[transferData] Transferring %s\n" % pname)
            try :
                p = pathlib.Path("%s/%s/sequencing_data/%s" % (
                    config.get("Paths","groupDir"),
                    group,
                    config.get("Options","runID")))
                if(p.exists() == False) :
                    p.mkdir(mode=0o750, parents=True)

                shutil.copytree(project, "%s/%s/sequencing_data/%s/%s" % (
                    config.get("Paths","groupDir"),
                    group,
                    config.get("Options","runID"),
                    project.split("/")[-1]))

                shutil.copytree("%s/%s/FASTQC_%s" % (
                    config.get("Paths","outputDir"),
                    config.get("Options","runID"),
                    project.split("/")[-1])
                    , "%s/%s/sequencing_data/%s/FASTQC_%s" % (
                    config.get("Paths","groupDir"),
                    group,
                    config.get("Options","runID"),
                    project.split("/")[-1]))

                for r, dirs, files in os.walk("%s/%s/sequencing_data/%s" % (
                    config.get("Paths","groupDir"),
                    group,
                    config.get("Options","runID"))):
                    for d in dirs:
                        os.chmod(os.path.join(r, d), stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP)
                    for f in files:
                        os.chmod(os.path.join(r, f), stat.S_IRWXU | stat.S_IRGRP)

                message += "\n%s\ttransferred" % pname
            except :
                e = sys.exc_info()
                message += "\n%s\tError during transfer (%s: %s)!" % (pname, e[0], e[1])
        elif(pname[0] == "B") :
            syslog.syslog("[transferData] Transferring %s\n" % pname)
            try :
                # The Schuele group has its own person that should get these
                if project.split("/")[-1].startswith("B01Schuele_"):
                    recipient = config.get("Uni","Schuele")
                else:
                    recipient = config.get("Uni","default")
                cmd = "tar cf - %s/%s/FASTQC_%s %s/%s/%s | fexsend -s %s_%s.tar %s" % (
                    config.get("Paths","outputDir"),
                    config.get("Options","runID"),
                    project.split("/")[-1],
                    config.get("Paths","outputDir"),
                    config.get("Options","runID"),
                    project.split("/")[-1],
                    config.get("Options", "runID"),
                    project.split("/")[-1],
                    recipient)
                rv = os.system(cmd)
                #if rv != 0:
                #    assert(1==0)
                message += "\n%s\ttransferred (return code %s from command '%s')" % (pname, rv, cmd)
            except :
                # fexsend doesn't return 0 on success
                message += "\n%s\ttransferred (return code %s from command '%s')" % (pname, rv, cmd)
        elif(pname[0] == "C") :
            syslog.syslog("[transferData] Transferring %s\n" % pname)
            try :
                p = pathlib.Path("%s/sequencing_data/%s" % (
                    config.get("Paths","DEEPDir"),
                    config.get("Options","runID")))
                if(p.exists() == False) :
                    p.mkdir(mode=0o770, parents=True)

                shutil.copytree("%s/%s/FASTQC_%s" % (
                    config.get("Paths","outputDir"),
                    config.get("Options","runID"),
                    project.split("/")[-1]),
                    "%s/sequencing_data/%s/FASTQC_%s" % (
                    config.get("Paths","DEEPDir"),
                    config.get("Options","runID"),
                    project.split("/")[-1]))

                shutil.copytree("%s/%s/%s" % (
                    config.get("Paths","outputDir"),
                    config.get("Options","runID"),
                    project.split("/")[-1])
                    , "%s/sequencing_data/%s/%s" % (
                    config.get("Paths","DEEPDir"),
                    config.get("Options","runID"),
                    project.split("/")[-1]))
                message += "\n%s\ttransferred" % pname
            except :
                message += "\n%s\tError during transfer!" % pname
        else :
            message += "\n%s\tskipped" % pname
    return message

def getSampleIDNameProjectLaneTuple(config) :
    """
    Parse a sample sheet to get a tuple of:
      * Sample ID
      * Sample Name
      * Lane
      * Project Name
    This can then be used in the project PDFs
    """
    samples = []
    inBottom = False
    if(config.get("Options", "sampleSheet") == "" or os.path.isfile(config.get("Options", "sampleSheet")) == False) :
        syslog.syslog("[getSampleIDNameProjectLaneTuple] No sample sheet! This *must* be an unindexed project.\n")
        return None

    for line in csv.reader(codecs.open("%s" % config.get("Options", "sampleSheet"), "r", "iso-8859-1")) :
        if len(line) == 0:
            continue
        if(inBottom) :
            samples.append([line[1],line[2],line[0],line[7]])
        else :
            if(line[0] == "Lane") :
                inBottom = True
    if(inBottom is False) :
        syslog.syslog("[getSampleIDNameProjectLaneTuple] Apparently the sample sheet couldn't properly be parsed.\n")
        return None
    return samples

def getSampleID(sampleTuple, project, lane, sampleName) :
    if(sampleTuple is None) :
        return " "
    for item in sampleTuple :
        if(sampleName == item[1] and
            lane == item[2] and
            project == item[3]) :
            return item[0]
    return " "

def makeProjectPDF(node, project, config) :
    """
    Uses reportlab to generate a per-project PDF file.

    Contents are currently a table containing:
      * Sample
      * Barcode
      * Lane
      * # Reads (passing filter)
      * % Bases > Q30
      * Average base quality
    For paired-end datasets, the last two columns are repeated and named differently.
    """
    st = getSampleIDNameProjectLaneTuple(config)

    stylesheet=getSampleStyleSheet()

    pdf = BaseDocTemplate("%s/%s/Project_%s/SequencingReport.pdf" % (
        config.get("Paths","outputDir"),config.get("Options","runID"),
        project), pagesize=landscape(A4))
    topHeight=100 #The image is 86 pixels tall
    fTL = Frame(pdf.leftMargin, pdf.height, width=pdf.width/2, height=topHeight, id="col1") #Fixed height
    fTR = Frame(pdf.leftMargin+pdf.width/2, pdf.height, width=pdf.width/2, height=topHeight, id="col2")
    fB = Frame(pdf.leftMargin, pdf.bottomMargin, pdf.width, pdf.height-topHeight, id="bottom")
    fM = Frame(pdf.leftMargin, pdf.bottomMargin, pdf.width, pdf.height, id="main")
    
    elements = []
    PE = False
    if(len(node[0][0][0][0][1]) == 3) :
        PE = True
        data = [["Sample ID","Sample Name", "Barcode","Lane","# Reads","% Bases\n>= Q30\nRead #1","Ave. Qual.\nRead #1","% Bases\n>= Q30\nRead #2","Ave. Qual.\nRead #2"]]
    else :
        data = [["Sample ID","Sample Name", "Barcode","Lane","# Reads","% Bases\n>= Q30","Ave. Qual."]]

    #A text blurb
    string = "Project: %s" % project
    p = Paragraph(string, style=stylesheet['Normal'])
    elements.append(p)
    string = "Report generated: %s" % (strftime("%d-%m-%Y %H:%M:%S"))
    p = Paragraph(string, style=stylesheet['Normal'])
    elements.append(p)
    string = "Flow cell ID: %s" % (config.get("Options","runID"))
    p = Paragraph(string, style=stylesheet['Normal'])
    elements.append(p)
    string = "BCL2Fastq pipeline version: %s" % (config.get("Version","pipeline"))
    p = Paragraph(string, style=stylesheet['Normal'])
    elements.append(p)
    string = "bcl2fastq version: %s" % (config.get("Version","bcl2fastq"))
    p = Paragraph(string, style=stylesheet['Normal'])
    elements.append(p)
    try:
        readLength = int(int(node[0][0][0][0][0][1][0].text)/int(node[0][0][0][0][0][0].text))
    except:
        readLength = 0
    string = "FastQC version: %s" % (config.get("Version","fastQC"))
    p = Paragraph(string, style=stylesheet['Normal'])
    elements.append(p)
    if(PE) :
        string = "%i base paired-end reads" % readLength
    else :
        string = "%i base single-end reads" % readLength
    p = Paragraph(string, style=stylesheet['Normal'])
    elements.append(p)

    #Image
    #Scale things
    img = utils.ImageReader(config.get("Options","imagePath"))
    iw, ih = img.getSize()
    iw = 0.88*iw*topHeight/ih
    elements.append(Image(config.get("Options","imagePath"), width=iw, height=0.88*topHeight, hAlign="RIGHT"))
    elements.append(NextPageTemplate("RemainingPages"))

    #Data table
    for sample in node.findall("Sample") :
        e = [None,None,-1,0,0,0,0,0,0]
        e[0] = sample.get("name")
        if(e[0] == "all") :
            continue
        for barcode in sample.findall("Barcode") :
            e[1] = barcode.get("name")
            if(e[1] == "all") :
                continue
            for lane in barcode.findall("Lane") :
                e[2] = lane.get("number")
                e[3] = 0
                e[4] = 0
                e[5] = 0
                e[6] = 0
                e[7] = 0
                e[8] = 0
                for tile in lane.findall("Tile") :
                    e[3] += int(tile[1][0].text) #Pf->ClusterCount
                    e[4] += int(tile[1][1][0].text) #Pf->Read1->Yield
                    e[5] += int(tile[1][1][1].text) #Pf->Read1->YieldQ30
                    e[6] += int(tile[1][1][2].text) #Pf->Read1->QualSum
                    if(PE) :
                        e[7] += int(tile[1][2][1].text) #Pf->Read2->YieldQ30
                        e[8] += int(tile[1][2][2].text) #Pf->Read2->QualSum
                if(PE) :
                    try:
                        data.append([getSampleID(st, project, e[2], e[0]),
                                     e[0],
                                     e[1],
                                     e[2],
                                     e[3],
                                     "%5.2f" % (100*(e[5]/e[4])),
                                     "%5.2f" % (e[6]/e[4]),
                                     "%5.2f" % (100*(e[7]/e[4])),
                                     "%5.2f" % (e[8]/e[4])
                            ])
                    except:
                        data.append([getSampleID(st, project, e[2], e[0]),
                                     e[0],
                                     e[1],
                                     e[2],
                                     e[3],
                                     "NA",
                                     "NA",
                                     "NA",
                                     "NA"
                            ])
                else :
                    try:
                        data.append([getSampleID(st, project, e[2], e[0]),
                                     e[0],
                                     e[1],
                                     e[2],
                                     e[3],
                                     "%5.2f" % (100*(e[5]/e[4])),
                                     "%5.2f" % (e[6]/e[4])
                            ])
                    except:
                        data.append([getSampleID(st, project, e[2], e[0]),
                                     e[0],
                                     e[1],
                                     e[2],
                                     e[3],
                                     "NA",
                                     "NA"
                            ])

    t = Table(data, style=[
        ('ROWBACKGROUNDS', (0, 0), (-1, -1), (0xD3D3D3, None)) #Light grey
        ], repeatRows=1)
    elements.append(t)

    #Add the key
    elements.append(Spacer(0,30))
    key = []
    key.append([Paragraph("Sample ID",
            stylesheet['BodyText']),
        Paragraph("The sample ID as provided to the sequencer in the sample sheet. This may not match the final file name, but will match the directory in which it's held.", 
            stylesheet['BodyText'])])
    key.append([Paragraph("Sample Name",
            stylesheet['BodyText']),
        Paragraph("The sample name as provided to the sequencer in the sample sheet. This should match the final file name.",
            stylesheet['BodyText'])])
    key.append([Paragraph("Barcode",
            stylesheet['BodyText']),
        Paragraph("The sample barcode added by the sequencing facility (or you, if you created the libraries yourself). This will generally be 6 nucleotides long.",
            stylesheet['BodyText'])])
    key.append([Paragraph("Lane", 
            stylesheet['BodyText']),
        Paragraph("The lane number on the flow cell (there are 8 of them).",
            stylesheet['BodyText'])])
    key.append([Paragraph("# Reads", 
            stylesheet['BodyText']),
        Paragraph("The number of reads in a given file. For paired-end datasets, this is equivalent to the number of fragments sequenced, rather than summing the counts for read #1 and read #2. Note that this includes only reads passing the quality filter.",
            stylesheet['BodyText'])])
    key.append([Paragraph("% Bases >= Q30 Read #1", 
            stylesheet['BodyText']),
        Paragraph("The percentage of bases in read #1 of a pair having a Phred-scaled score of at least 30, meaning that the 0.1% or less chance that they're incorrect.",
            stylesheet['BodyText'])])
    key.append([Paragraph("Ave. Qual. Read #1", 
            stylesheet['BodyText']),
        Paragraph("The average Phred-scaled base quality of bases in read #1 of a pair. This number of -10*log10(Probability that the call is incorrect). In other words, if a call is 100% likely to be wrong, the score is 0 (or 10 for 10% likelihood, 20 for 1% likelihood, etc.).",
            stylesheet['BodyText'])])
    key.append([Paragraph("% Bases >= Q30 Read #2", 
            stylesheet['BodyText']),
        Paragraph("Identical to '% Bases >= Q30 Read #1', but for read #2 of a pair.",
            stylesheet['BodyText'])])
    key.append([Paragraph("Ave. Qual. Read #2", 
            stylesheet['BodyText']),
        Paragraph("Identical to 'Ave. Qual. Read #1', but for read #1 of a pair.",
            stylesheet['BodyText'])])
    key.append([Paragraph("# Reads", 
            stylesheet['BodyText']),
        Paragraph("Identical to '% Bases >= Q30 Read #1', but for single-end datasets.",
            stylesheet['BodyText'])])
    key.append([Paragraph("Ave. Qual.", 
            stylesheet['BodyText']),
        Paragraph("Identical to 'Ave. Qual. Read #1', but for single-end datasets.",
            stylesheet['BodyText'])])
    t2 = Table(key, colWidths=(80, None))
    t2.setStyle(TableStyle([('VALIGN',(0,0),(-1,-1),'TOP')]))
    elements.append(t2)

    #fastq_screen images
    elements.append(PageBreak())
    elements.append(Paragraph("Contaminant screen", stylesheet['title']))
    elements.append(Paragraph("Below are images generated on the output of fastq_screen. In short, 1 million reads are randomly taken from each indicated sample. These reads are then aligned against a variety of organisms (mouse, human, etc.). The resulting alignments are then categorized as follows:", stylesheet['Normal']))
    elements.append(ListFlowable([
        Paragraph("unique: aligns only a single time within a single species.", stylesheet['Normal']),
        Paragraph("multimap: aligns multiple times, but only within a single species.", stylesheet['Normal']),
        Paragraph("conserved: aligns a single time to each of two or more species.", stylesheet['Normal']),
        Paragraph("repeat: aligns multiple times to each of two or more species.", stylesheet['Normal'])],
        bulletType='bullet',
        start='circle'))
    elements.append(Spacer(0,30))
    elements.append(Paragraph("Ideally, the 'unique' and 'multimap' values will only be appreciably present in the species from which your sample should have arisen.", stylesheet['Normal']))
    elements.append(Spacer(0,30))
    elements.append(Paragraph("Note that as the mouse and human reference genomes are the best quality, many low complexity reads will align to them.", stylesheet['Normal']))
    elements.append(Spacer(0,30))
    fqs = glob.glob("%s/%s/Project_%s/*/*.png" % (
        config.get("Paths","outputDir"),
        config.get("Options","runID"),
        project))
    fqs.sort()
    for fq in fqs:
        img = utils.ImageReader(fq)
        iw, ih = img.getSize()
        iw = 0.7*iw
        ih = 0.7*ih
        elements.append(Image(fq, width=iw, height=ih, hAlign="LEFT"))

    pdf.addPageTemplates([PageTemplate(id="FirstPage", frames=[fTL, fTR, fB]),
        PageTemplate(id="RemainingPages", frames=[fM])]),
    pdf.build(elements)

def getFCmetrics(root) :
    barcode = root[0][0] #Sample "all", barcode "all"
    message = "Lane\t# Clusters (% pass)\t% Bases >=Q30\tAve. base qual.\n"
    for lane in barcode.findall("Lane") :
        message += "Lane %s" % lane.get("number")
        clusterCount = 0
        clusterCountPass = 0
        baseYield = [0,0]
        baseYieldQ30 = [0,0]
        QualSum = [0,0]
        rlens=[0,0]
        for tile in lane :
            clusterCount += int(tile[0][0].text)
            clusterCountPass += int(tile[1][0].text)
            #Yield
            baseYield[0] += int(tile[1][1][0].text)
            if(len(tile[1]) == 3) :
                baseYield[1] += int(tile[1][2][0].text)
            #YieldQ30
            baseYieldQ30[0] += int(tile[1][1][1].text)
            if(len(tile[1]) == 3) :
                baseYieldQ30[1] += int(tile[1][2][1].text)
            #QualSum
            QualSum[0] += int(tile[1][1][2].text)
            if(len(tile[1]) == 3) :
                QualSum[1] += int(tile[1][2][2].text)
        #Number of clusters (%passing filter)
        try:
            message += "\t%i (%5.2f%%)" % (clusterCount,100*clusterCountPass/clusterCount)
        except:
            message += "\t%i (NA)" % (clusterCount)
        #%bases above Q30
        if(baseYield[1] > 0) :
            try:
                message += "\t%5.2f%%/%5.2f%%" % (100*(baseYieldQ30[0]/baseYield[0]),
                    100*(baseYieldQ30[1]/baseYield[1]))
            except:
                message += "\tNA/NA"
        else :
            try:
                message += "\t%5.2f%%" % (100*(baseYieldQ30[0]/baseYield[0]))
            except:
                message += "\tNA"
        #Average base quality
        if(baseYield[1] > 0) :
            try:
                message += "\t%4.1f/%4.1f\n" % (QualSum[0]/float(baseYield[0]),
                    QualSum[1]/float(baseYield[1]))
            except:
                message += "\tNA/NA\n"
        else :
            try:
                message += "\t%4.1f\n" % (QualSum[0]/float(baseYield[0]))
            except:
                message += "\tNA\n"

    return message

def parseConversionStats(config) :
    """
    Parse ConversionStats.xml, producing:
     1) A PDF file for each project
     2) A message that will be included in the email message
    """
    try :
        tree = ET.parse("%s/%s/Stats/ConversionStats.xml" % (config.get("Paths","outputDir"),config.get("Options","runID")))
        root = tree.getroot()[0] #We only ever have a single flow cell
    except :
        return None
    metrics = None
    #Per-project PDF files
    for project in root.findall("Project") :
        if(project.get("name") == "default") :
            continue
        if(project.get("name") == "all") :
            metrics = getFCmetrics(project)
        else :
            makeProjectPDF(project, project.get("name"), config)
    return metrics

def enoughFreeSpace(config) :
    """
    Ensure that outputDir has at least minSpace gigs
    """
    (tot,used,free) = shutil.disk_usage(config.get("Paths","outputDir"))
    free /= 1024*1024*1024
    if(free >= float(config.get("Options","minSpace"))) :
        return True
    return False

def errorEmail(config, errTuple, msg) :
    msg = MIMEText(msg + "\nError type: %s\nError value: %s\n%s\n" % (errTuple[0], errTuple[1], errTuple[2]))
    msg['Subject'] = "[bcl2fastq_pipeline] Error"
    msg['From'] = config.get("Email","fromAddress")
    msg['To'] = config.get("Email","errorTo")

    s = smtplib.SMTP(config.get("Email","host"))
    s.send_message(msg)
    s.quit()

def finishedEmail(config, msg, runTime, transferTime) :
    message = "Flow cell: %s\n" % config.get("Options","runID")
    message += "Run time: %s\n" % runTime
    message += "Data transfer: %s\n" % transferTime
    message += msg

    msg = MIMEText(message)
    msg['Subject'] = "[bcl2fastq_pipeline] %s processed" % config.get("Options","runID")
    msg['From'] = config.get("Email","fromAddress")
    msg['To'] = config.get("Email","finishedTo")

    s = smtplib.SMTP(config.get("Email","host"))
    s.send_message(msg)
    s.quit()
