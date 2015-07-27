"""
Misc. functions
"""

import configparser
import shutil
import smtplib
from email.mime.text import MIMEText
import xml.etree.ElementTree as ET
from reportlab.lib import colors, utils
from reportlab.platypus import BaseDocTemplate, Table, Preformatted, Paragraph, Spacer, Image, Frame, NextPageTemplate, PageTemplate, TableStyle
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.pagesizes import A4, landscape
from time import strftime
from reportlab.pdfgen import canvas
import csv
import sys
import glob
import pathlib
import subprocess
import syslog
import codecs

def transferData(config) :
    """
    Distribute fastq and fastQC files to users.
    """
    message = ""
    projects = glob.glob("%s/%s/FASTQC_*" % (config.get("Paths","outputDir"),config.get("Options","runID")))
    for project in projects :
        pname = project.split("/")[-1][7:]
        if(pname[8] == "A") :
            #Copy
            group = pname.split("_")[2].lower()
            syslog.syslog("[transferData] Transferring %s\n" % pname)
            try :
                p = pathlib.Path("%s/%s/sequencing_data/%s" % (
                    config.get("Paths","groupDir"),
                    group,
                    config.get("Options","runID")))
                if(p.exists() == False) :
                    p.mkdir(mode=0o770, parents=True)
                shutil.copytree(project, "%s/%s/sequencing_data/%s/FASTQC_%s" % (
                    config.get("Paths","groupDir"),
                    group,
                    config.get("Options","runID"),
                    pname))
                shutil.copytree("%s/%s/%s" % (
                    config.get("Paths","outputDir"),
                    config.get("Options","runID"),
                    pname)
                    , "%s/%s/sequencing_data/%s/%s" % (
                    config.get("Paths","groupDir"),
                    group,
                    config.get("Options","runID"),
                    pname))
                subprocess.call(['chmod','-R','g-w', "%s/%s/sequencing_data/%s" % (
                    config.get("Paths","groupDir"),
                    group,
                    config.get("Options","runID"))])
                message += "\n%s\ttransferred" % pname
            except :
                message += "\n%s\tError during transfer!" % pname
        elif(pname[8] == "B") :
            syslog.syslog("[transferData] Transferring %s\n" % pname)
            try :
                p = pathlib.Path("%s/sequencing_data/%s" % (
                    config.get("Paths","UniDir"),
                    config.get("Options","runID")))
                if(p.exists() == False) :
                    p.mkdir(mode=0o770, parents=True)
                shutil.copytree(project, "%s/sequencing_data/%s/FASTQC_%s" % (
                    config.get("Paths","UniDir"),
                    config.get("Options","runID"),
                    pname))
                shutil.copytree("%s/%s/%s" % (
                    config.get("Paths","outputDir"),
                    config.get("Options","runID"),
                    pname)
                    , "%s/sequencing_data/%s/%s" % (
                    config.get("Paths","UniDir"),
                    config.get("Options","runID"),
                    pname))
                message += "\n%s\ttransferred" % pname
            except :
                message += "\n%s\tError during transfer!" % pname
        elif(pname[8] == "C") :
            syslog.syslog("[transferData] Transferring %s\n" % pname)
            try :
                p = pathlib.Path("%s/sequencing_data/%s" % (
                    config.get("Paths","DEEPDir"),
                    config.get("Options","runID")))
                if(p.exists() == False) :
                    p.mkdir(mode=0o770, parents=True)
                shutil.copytree(project, "%s/sequencing_data/%s/FASTQC_%s" % (
                    config.get("Paths","DEEPDir"),
                    config.get("Options","runID"),
                    pname))
                shutil.copytree("%s/%s/%s" % (
                    config.get("Paths","outputDir"),
                    config.get("Options","runID"),
                    pname)
                    , "%s/sequencing_data/%s/%s" % (
                    config.get("Paths","DEEPDir"),
                    config.get("Options","runID"),
                    pname))
                subprocess.call(['chmod','-R','g-w', "%s/sequencing_data/%s/%s" % (
                    config.get("Paths","DEEPDir"),
                    config.get("Options","runID"),
                    pname)])
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
    for line in csv.reader(codecs.open("%s/%s/SampleSheet.csv" % (config.get("Paths","baseDir"),config.get("Options","runID")), "r","iso-8859-1")) :
        if(inBottom) :
            if(len(line) == 7) :
                samples.append([line[1], line[6], line[0], line[5]])
            else :
                samples.append([line[1], line[2], line[0], line[7]])
        else :
            if(len(line) == 0) :
                continue
            if(line[0] == "Lane") :
                inBottom = True
    if(inBottom is False) :
        syslog.syslog("[getSampleIDNameProjectLaneTuple] Apparently the sample sheet couldn't properly be parsed.\n")
        return None
    return samples

def getSampleName(sampleTuple, project, lane, sampleID) :
    if(sampleTuple is None) :
        return " "
    for item in sampleTuple :
        if(sampleID == item[0] and
            lane == item[2] and
            project == item[3]) :
            return item[1]
    return " "

def makeProjectPDF(project, matrix, config) :
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
    if(matrix[0][9]>0) :
        PE = True
    if(PE) :
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
    readLength = int(matrix[0][12])
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
    for line in matrix :
        if(project != line[4]) :
            continue
        if(PE) :
            data.append([line[1],
                line[2],
                line[3],
                line[0],
                line[6],
                "%5.2f" % (100*(line[8]/line[7])),
                "%5.2f" % (line[10]/line[7]),
                "%5.2f" % (100*(line[9]/line[7])),
                "%5.2f" % (line[11]/line[7])])
        else : 
            data.append([line[1],
                line[2],
                line[3],
                line[0],
                line[6],
                "%5.2f" % (100*(line[8]/line[7])),
                "%5.2f" % (line[10]/line[7])])

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

    pdf.addPageTemplates([PageTemplate(id="FirstPage", frames=[fTL, fTR, fB]),
        PageTemplate(id="RemainingPages", frames=[fM])]),
    pdf.build(elements)

def getID(matrix, lane, sampleID, barcode) :
    idx = 0
    for row in matrix :
        if(row[3] == "" or barcode == "NoIndex") :
            if(row[0] == lane and row[1] == sampleID) :
                return idx
        else :
            if(row[0] == lane and row[1] == sampleID and row[3] == barcode) :
                return idx
        idx+=1
    return None

def getFCmetrics(root, matrix) :
    for lane in root.findall("Lane") :
        lnum = lane.get("index")
        for sample in lane.findall("Sample") :
            sampleID = sample.get("index")
            for barcode in sample.findall("Barcode") :
                BC = barcode.get("index")
                idx = getID(matrix, lnum, sampleID, BC)
                if idx is None:
                    continue
                for tile in barcode.findall("Tile") :
                    matrix[idx][5] += int(tile[0][0][2].text) #clusterCount
                    matrix[idx][6] += int(tile[0][1][2].text) #clusterCountPass
                    matrix[idx][7] += int(tile[0][1][0].text) #baseYield #1
                    matrix[idx][8] += int(tile[0][1][1].text) #baseYieldQ30 #1
                    matrix[idx][10] += int(tile[0][1][5].text) #QualSum #1
                    matrix[idx][12] = int(tile[0][1][0].text)/int(tile[0][1][2].text) #rlen #1
                    if(len(tile) > 1) :
                        matrix[idx][9] += int(tile[1][1][1].text) #baseYieldQ30 #2
                        matrix[idx][11] += int(tile[1][1][5].text) #QualSum #2
    return matrix

def getFCLaneMetrics(root) :
    message = "Lane\t# Clusters (% pass)\t% Bases >=Q30\tAve. base qual.\n"
    for lane in root.findall("Lane") :
        lnum = lane.get("index")
        clusterCount = 0
        clusterCountPass = 0
        baseYield = [0,0]
        baseYieldQ30 = [0,0]
        QualSum = [0,0]
        rlens = [0,0]
        for sample in lane.findall("Sample") :
            sampleID = sample.get("index")
            for barcode in sample.findall("Barcode") :
                BC = barcode.get("index")
                if(BC == "Undetermined") :
                    continue
                for tile in barcode.findall("Tile") :
                    clusterCount += int(tile[0][0][2].text) #clusterCount
                    clusterCountPass += int(tile[0][1][2].text) #clusterCountPass
                    baseYield[0] += int(tile[0][1][0].text) #baseYield #1
                    baseYieldQ30[0] += int(tile[0][1][1].text) #baseYieldQ30 #1
                    QualSum[0] += int(tile[0][1][5].text) #QualSum #1
                    rlens[0] = int(tile[0][1][0].text)/int(tile[0][1][2].text) #rlen #1
                    if(len(tile) > 1) :
                        baseYield[1] += int(tile[1][1][0].text) #baseYield #2
                        baseYieldQ30[1] += int(tile[1][1][1].text) #baseYieldQ30 #2
                        QualSum[1] += int(tile[1][1][5].text) #QualSum #2
                        rlens[1] = baseYield[1]/clusterCountPass 
        message += "Lane %s" % lnum
        message += "\t%i (%5.2f%%)" % (clusterCount,100*clusterCountPass/clusterCount)
        if(baseYieldQ30[1]>0) :
            message += "\t%5.2f%%/%5.2f%%" % (100*(baseYieldQ30[0]/baseYield[0]),
                100*(baseYieldQ30[1]/baseYield[1]))
            message += "\t%4.1f/%4.1f\n" % (QualSum[0]/baseYield[0],
                QualSum[1]/baseYield[1])
        else :
            message += "\t%5.2f%%" % (100*(baseYieldQ30[0]/baseYield[0]))
            message += "\t%4.1f\n" % (QualSum[0]/baseYield[0])
    return message

def parseConversionStats(config) :
    """
    The file name hasn't really beem changed from the bcl2fastq2 version.
    This actually creates a per-project PDF file by parsing Flowcell_demux_summary.xml
     1) A PDF file for each project
     2) A message that will be included in the email message
    """
    #Make an array of tuples to hold the metrics
    matrix = []
    inLane=False
    projects = set()
    for line in codecs.open("%s/%s/SampleSheet.csv" % (config.get("Paths","baseDir"), config.get("Options","runID")),"r","iso-8859-1") :
        line = line.strip().split(",")
        if(inLane) :
            #Lane, SampleID, SampleName, Barcode, Project,
            #clusterCount, clusterCountPass, baseYield*2,
            #baseYieldQ30*2, QualSum*2, rlens*2
            if(len(line) == 7) :
                matrix.append([line[0], line[1], line[6], line[2], line[5], 0,0,0,0,0,0,0,0,0,0])
                if line[5] not in projects :
                    projects.add(line[5])
            else :
                matrix.append([line[0],line[1],line[2],line[6],line[7], 0,0,0,0,0,0,0,0,0,0])
                if line[7] not in projects :
                    projects.add(line[7])
        else :
            if(line[0] == "Lane") :
                inLane = True
        
    FCID = config.get("Options","runID").split("_")[-1][1:]
    tree = ET.parse("%s/%s/Basecall_Stats_%s/Flowcell_demux_summary.xml" % (config.get("Paths","outputDir"),config.get("Options","runID"), FCID))
    root = tree.getroot()
    matrix = getFCmetrics(root, matrix)
    
    #Per-project PDF files
    for project in projects :
        makeProjectPDF(project, matrix, config)

    return getFCLaneMetrics(root)

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
