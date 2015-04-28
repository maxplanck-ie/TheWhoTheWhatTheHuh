"""
Misc. functions
"""

import configparser
import shutil
import smtplib
from email.mime.text import MIMEText
import xml.etree.ElementTree as ET
from reportlab.lib import colors
from reportlab.platypus import SimpleDocTemplate, Table, Paragraph, Spacer, Image, Frame
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.pagesizes import A4, landscape
from time import strftime

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

    To Do:
    Add an image?
    Footer?
    """
    pdf = SimpleDocTemplate("%s/%s/%s/SequencingReport.pdf" % (
        config.get("Paths","outputDir"),config.get("Options","runID"),
        project), pagesize=A4)
    f1 = Frame(0, 0, width=5, height=4)
    f2 = Frame(5, 4, width=5, height=4)
    elements = []
    PE = False
    if(len(node[0][0][0][0][1]) == 3) :
        PE = True
        data = [["Sample ID","Barcode","Lane","# Reads","% Bases\n>= Q30\nRead #1","Ave. Qual.\nRead #1","% Bases\n>= Q30\nRead #2","Ave. Qual.\nRead #2"]]
    else :
        data = [["Sample ID","Barcode","Lane","# Reads","% Bases\n>= Q30","Ave. Qual."]]

    #A text blurb
    stylesheet=getSampleStyleSheet()
    string = "Project: %s" % project
    p = Paragraph(string, style=stylesheet['Normal'])
    elements.append(p)
    string = "Report generated: %s" % (strftime("%d-%m-%Y %H:%M:%S"))
    p = Paragraph(string, style=stylesheet['Normal'])
    elements.append(p)
    string = "BCL2Fastq pipeline version: %s" % (config.get("Version","pipeline"))
    p = Paragraph(string, style=stylesheet['Normal'])
    elements.append(p)
    string = "bcl2fastq version: %s" % (config.get("Version","bcl2fastq"))
    p = Paragraph(string, style=stylesheet['Normal'])
    elements.append(p)
    string = "FastQC version: %s" % (config.get("Version","fastQC"))
    p = Paragraph(string, style=stylesheet['Normal'])
    elements.append(p)
    elements.append(Spacer(1,72))

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
                for tile in lane.findall("Tile") :
                    e[3] += int(tile[1][0].text) #Pf->ClusterCount
                    e[4] += int(tile[1][1][0].text) #Pf->Read1->Yield
                    e[5] += int(tile[1][1][1].text) #Pf->Read1->YieldQ30
                    e[6] += int(tile[1][1][2].text) #Pf->Read1->QualSum
                    if(PE) :
                        e[7] += int(tile[1][2][1].text) #Pf->Read2->YieldQ30
                        e[8] += int(tile[1][2][2].text) #Pf->Read2->QualSum
            if(PE) :
                data.append([e[0],
                             e[1],
                             e[2],
                             e[3],
                             "%5.2f" % (100*(e[5]/e[4])),
                             "%5.2f" % (e[6]/e[4]),
                             "%5.2f" % (100*(e[7]/e[4])),
                             "%5.2f" % (e[8]/e[4])
                    ])
            else :
                data.append([e[0],
                             e[1],
                             e[2],
                             e[3],
                             "%5.2f" % (100*(e[5]/e[4])),
                             "%5.2f" % (e[6]/e[4])
                    ])
                data.append([e[0],e[1],e[2],e[3],100*(e[5]/e[4]),e[6]/e[4]])
    t = Table(data, style=[
        ('ROWBACKGROUNDS', (0, 0), (-1, -1), (0xD3D3D3, None)) #Light grey
        ], repeatRows=1)
    elements.append(t)
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
        message += "\t%i (%5.2f%%)" % (clusterCount,100*clusterCountPass/clusterCount)
        #%bases above Q30
        if(baseYield[1] > 0) :
            message += "\t%5.2f%%/%5.2f%%" % (100*(baseYieldQ30[0]/baseYield[0]),
                100*(baseYieldQ30[1]/baseYield[1]))
        else :
            message += "\t%5.2f%%" % (100*(baseYieldQ30[0]/baseYield[0]))
        #Average base quality
        if(baseYield[1] > 0) :
            message += "\t%4.1f/%4.1f\n" % (QualSum[0]/clusterCountPass/(baseYield[0]/clusterCount),
                QualSum[1]/clusterCountPass/(baseYield[1]/clusterCount))
        else :
            message += "\t%4.1f\n" % (QualSum[0]/clusterCountPass/(baseYield[0]/clusterCount))

    return message

def parseConversionStats(config) :
    """
    Parse ConversionStats.xml, producing:
     1) A PDF file for each project
     2) A message that will be included in the email message
    """
    tree = ET.parse("%s/%s/Stats/ConversionStats.xml" % (config.get("Paths","outputDir"),config.get("Options","runID")))
    root = tree.getroot()[0] #We only ever have a single flow cell
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

def errorEmail(config, msg) :
    msg = MIMEText(msg)
    msg['Subject'] = "[bcl2fastq_pipeline] Error"
    msg['From'] = config.get("Email","fromAddress")
    msg['To'] = config.get("Email","errorTo")

    s = smtplib.SMTP(config.get("Email","host"))
    s.send_message(msg)
    s.quit()

def finishedEmail(config, msg, runTime) :
    message = "Flow cell: %s\n" % config.get("Options","runID")
    message += "Run time: %s\n" % runTime
    message += msg

    print(message)
    msg = MIMEText(message)
    msg['Subject'] = "[bcl2fastq_pipeline] %s processed" % config.get("Options","runID")
    msg['From'] = config.get("Email","fromAddress")
    msg['To'] = config.get("Email","finishedTo")

    s = smtplib.SMTP(config.get("Email","host"))
    s.send_message(msg)
    s.quit()
