"""
Misc. functions
"""

import configparser
import shutil
import smtplib
from email.mime.text import MIMEText
import xml.etree.ElementTree as ET

def makeProjectPDF(node, project, config) :
    return

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
