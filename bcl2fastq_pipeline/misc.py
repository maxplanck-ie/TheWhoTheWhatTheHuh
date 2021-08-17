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
import stat
import codecs
import requests
import json

def fetchGalaxyUsers(userFile):
    l = list()
    with open(userFile) as f:
        for line in f:
            l.append(line.strip().split()[1]) #Grab last name
    return l

def getLatestSeqdir(groupData, PI):
    seqDirNum = 0
    for dirs in os.listdir(os.path.join(groupData, PI)):
        if 'sequencing_data' in dirs:
            seqDirStrip = dirs.replace('sequencing_data','')
            if seqDirStrip is not '':
                if int(seqDirStrip) > seqDirNum:
                    seqDirNum = int(seqDirStrip)
    if seqDirNum == 0:
        return 'sequencing_data'
    else:
        return 'sequencing_data' + str(seqDirNum)

def transferData(config) :
    """
    Distribute fastq and fastQC files to users.
    """
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes = "_lanes{}".format(lanes)

    message = ""
    projects = glob.glob("%s/%s%s/Project_*" % (config.get("Paths","outputDir"),config.get("Options","runID"), lanes))
    for project in projects :
        pname = project.split("/")[-1][8:]
        group = pname.split("_")[-1].lower()
        if "-" in group:
            # Handle things like cabezas-wallschied -> cabezas
            group =  group.split("-")[0]
        syslog.syslog("[transferData] Transferring %s\n" % pname)

        # Get the latest sequencing data folder.
        if os.path.exists("{}/{}/sequencing_data".format(config.get("Paths","groupDir"), group)):
            # Get latest sequencing data folder (after the check for regular sequencing data to discriminate university runs.
            latestSeqDir = getLatestSeqdir(config.get("Paths","groupDir"), group)
            #Local group
            syslog.syslog("[transferData] Transferring %s\n" % pname)
            try :
                p = pathlib.Path("%s/%s/%s/%s%s" % (
                    config.get("Paths","groupDir"),
                    group,
                    latestSeqDir,
                    config.get("Options","runID"),
                    lanes))
                if(p.exists() == False) :
                    p.mkdir(mode=0o750, parents=True)

                shutil.copytree(project, "%s/%s/%s/%s%s/%s" % (
                    config.get("Paths","groupDir"),
                    group,
                    latestSeqDir,
                    config.get("Options","runID"),
                    lanes,
                    project.split("/")[-1]))

                shutil.copytree("%s/%s%s/FASTQC_%s" % (
                    config.get("Paths","outputDir"),
                    config.get("Options","runID"),
                    lanes,
                    project.split("/")[-1])
                    , "%s/%s/%s/%s%s/FASTQC_%s" % (
                    config.get("Paths","groupDir"),
                    group,
                    latestSeqDir,
                    config.get("Options","runID"),
                    lanes,
                    project.split("/")[-1]))

                for r, dirs, files in os.walk("%s/%s/%s/%s%s" % (
                    config.get("Paths","groupDir"),
                    group,
                    latestSeqDir,
                    config.get("Options","runID"),
                    lanes)):
                    for d in dirs:
                        os.chmod(os.path.join(r, d), 0o700)
                    for f in files:
                        os.chmod(os.path.join(r, f), 0o700)

                message += "\n%s\ttransferred" % pname
            except :
                e = sys.exc_info()
                message += "\n%s\tError during transfer (%s: %s)!" % (pname, e[0], e[1])
        else:
            #Upload with FEX
            syslog.syslog("[transferData] Transferring %s\n" % pname)
            try :
                # The Schuele group has its own person that should get these
                if project.split("/")[-1].startswith("B01Schuele_"):
                    recipient = config.get("Uni","Schuele")
                else:
                    recipient = config.get("Uni","default")
                cmd = "tar cf - %s/%s%s/FASTQC_%s %s/%s%s/%s | fexsend -s %s%s_%s.tar %s" % (
                    config.get("Paths","outputDir"),
                    config.get("Options","runID"),
                    lanes,
                    project.split("/")[-1],
                    config.get("Paths","outputDir"),
                    config.get("Options","runID"),
                    lanes,
                    project.split("/")[-1],
                    config.get("Options", "runID"),
                    lanes,
                    project.split("/")[-1],
                    recipient)
                rv = os.system(cmd)
                message += "\n%s\ttransferred (return code %s from command '%s')" % (pname, rv, cmd)
            except :
                # fexsend doesn't return 0 on success
                message += "\n%s\ttransferred (return code %s from command '%s')" % (pname, rv, cmd)
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
    sampleIDcol = None
    sampleNameCol = None
    laneCol = None
    projectNameCol = None
    if(config.get("Options", "sampleSheet") == "" or os.path.isfile(config.get("Options", "sampleSheet")) == False) :
        syslog.syslog("[getSampleIDNameProjectLaneTuple] No sample sheet! This *must* be an unindexed project.\n")
        return None

    for line in csv.reader(codecs.open("%s" % config.get("Options", "sampleSheet"), "r", "iso-8859-1")) :
        if len(line) == 0:
            continue
        if(inBottom) :
            samples.append([line[sampleIDcol],line[sampleNameCol],line[laneCol],line[projectNameCol]])
        else :
            if "Lane" in line:
                try:
                    laneCol = line.index("Lane")
                    sampleIDcol = line.index("Sample_ID")
                    sampleNameCol = line.index("Sample_Name")
                    projectNameCol = line.index("Sample_Project")
                    inBottom = True
                except:
                    pass
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
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes = "_lanes{}".format(lanes)

    st = getSampleIDNameProjectLaneTuple(config)

    stylesheet=getSampleStyleSheet()

    pdf = BaseDocTemplate("%s/%s%s/Project_%s/SequencingReport.pdf" % (
        config.get("Paths","outputDir"),config.get("Options","runID"), lanes,
        project), pagesize=landscape(A4))
    topHeight=120 #The image is 86 pixels tall
    fTL = Frame(pdf.leftMargin, pdf.height, width=pdf.width/2, height=topHeight, id="col1") #Fixed height
    fTR = Frame(pdf.leftMargin+pdf.width/2, pdf.height, width=pdf.width/2, height=topHeight, id="col2")
    fB = Frame(pdf.leftMargin, pdf.bottomMargin, pdf.width, pdf.height-topHeight, id="bottom")
    fM = Frame(pdf.leftMargin, pdf.bottomMargin, pdf.width, pdf.height, id="main")
    
    elements = []
    PE = False
    if(len(node[0][0][0][0][1]) == 3) :
        PE = True
#        data = [["Sample ID","Sample Name", "Barcode","Lane","# Reads","% Bases\n>= Q30\nRead #1","Ave. Qual.\nRead #1","% Bases\n>= Q30\nRead #2","Ave. Qual.\nRead #2"]]
        data = [["Sample ID","Sample Name", "Barcode(s)","Lane(s)","# Reads (million)","% Bases\n>= Q30\nRead #1","% Bases\n>= Q30\nRead #2"]]
    else :
#        data = [["Sample ID","Sample Name", "Barcode","Lane","# Reads","% Bases\n>= Q30","Ave. Qual."]]
        data = [["Sample ID","Sample Name", "Barcode(s)","Lane(s)","# Reads (million)","% Bases\n>= Q30"]]
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
    FC = config.get("Options", "runID")
    if FC[7] == 'M':
        string = "Sequencer type: MiSeq"
    elif FC[7:9] == 'NB':
        string = "Sequencer type: NextSeq 500"
    elif FC[7:9] == 'SN':
        string = "Sequencer type: HiSeq 2500"
    elif FC[7] == 'J':
        string = "Sequencer type: HiSeq 3000"
    elif FC[7] == 'A':
        string = "Sequencer type: NovaSeq 6000"
    else:
        string = "Sequencer type: Unknown"
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
    #string = "FastQC version: %s" % (config.get("Version","fastQC"))
    #p = Paragraph(string, style=stylesheet['Normal'])
    #elements.append(p)
    if(PE) :
        string = "%i base paired-end reads" % readLength
    else :
        string = "%i base single-end reads" % readLength
    p = Paragraph(string, style=stylesheet['Normal'])
    elements.append(p)
    # Try to fetch the libtype and protocol.
    FCID = FC.split("_")[3][1:]
    if '-' in FCID:
        FCID = FCID.split('-')[-1]
    d = {'flowcell_id': FCID}
    res = requests.get( config.get("parkour", "QueryURL"), auth=(config.get("parkour", "user"), config.get("parkour", "password")), params=d )
    if res.status_code == 200:
        resDic = res.json()
        projectIndex = project.split('_')[0]
        for p in resDic:
            if projectIndex in p:
                projMatch = p
        if not projMatch:
            p = Paragraph("Lib. Type: Unknown", style=stylesheet['Normal'])
            elements.append(p)
            p = Paragraph("Protocol: Unknown", style=stylesheet['Normal'])
            elements.append(p)
        else:
            libType = set( [ resDic[projMatch][sam][1] for sam in resDic[projMatch] ] )
            p = Paragraph("Lib. Type: {}".format(','.join(libType)), style=stylesheet['Normal'])
            elements.append(p)
            protType = set( [ resDic[projMatch][sam][2] for sam in resDic[projMatch] ] )
            p = Paragraph("Protocol: {}".format(','.join(protType)), style=stylesheet['Normal'])
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
                                     #"%5.2f" % (e[6]/e[4]),
                                     "%5.2f" % (100*(e[7]/e[4]))
                                     #"%5.2f" % (e[8]/e[4])
                            ])
                    except:
                        data.append([getSampleID(st, project, e[2], e[0]),
                                     e[0],
                                     e[1],
                                     e[2],
                                     e[3],
                                     "NA",
                                     #"NA",
                                     "NA"
                                     #"NA"
                            ])
                else :
                    try:
                        data.append([getSampleID(st, project, e[2], e[0]),
                                     e[0],
                                     e[1],
                                     e[2],
                                     e[3],
                                     "%5.2f" % (100*(e[5]/e[4]))
                                     #"%5.2f" % (e[6]/e[4])
                            ])
                    except:
                        data.append([getSampleID(st, project, e[2], e[0]),
                                     e[0],
                                     e[1],
                                     e[2],
                                     e[3],
                                     "NA"
                                     #"NA"
                            ])
    # Iterate over data to collapse samples on different lanes.
    print("DataEntries:")
    for i in data:
        print(i)
    dataDic = {}
    for row in data:
        if row[0] == 'Sample ID':
            dataDic['Header'] = row
        else:
            if row[1] in dataDic:
                dataDic[ row[1] ].append(row)
            else:
                dataDic[ row[1] ] = [row]

    # Fetch the barcode ID from the original sampleSheet.
    fulSS =  os.path.join(config.get("Paths", "baseDir"), config.get("Options","runID"), 'SampleSheet.csv' )
    sampleBarCodeIxDic = {}
    if os.path.exists(fulSS):
        with open(fulSS) as f:
            for line in f:
                if project in line:
                    sampleBarCodeIxDic[ line.strip().split(',')[2] ] = ','.join([ line.strip().split(',')[5], line.strip().split(',')[7]] )

    dataSquashed = []
    for entry in dataDic:
        if entry == 'Header':
            headerUpd = dataDic[entry][0:3]
            headerUpd.append("Barcode ID(s)")
            headerUpd = headerUpd + dataDic[entry][3::]
            dataSquashed.append( headerUpd )
            headerLen = len( dataDic[entry] )
        else:
            collapsedSample = []
            collapsedSample = [ dataDic[entry][0][0], dataDic[entry][0][1], dataDic[entry][0][2] ] #Add sample ID, sample Name and Barcode
            collapsedSample.append( sampleBarCodeIxDic[ dataDic[entry][0][1] ]  )
            collapsedLanes = ', '.join( [laneNum[3] for laneNum in dataDic[entry] ] )
            collapsedSample.append(collapsedLanes)
            collapsedReads = round( sum( [int( readNum[4] ) for readNum in dataDic[entry]] ) / 1000000 , 2)
            collapsedSample.append(collapsedReads)
            try:
                avQualR1 = [ float(qual[5]) for qual in dataDic[entry]  ]
                avQualR1 = round( ( sum(avQualR1) / len(avQualR1) ), 2)
            except:
                avQualR1 = 0
            collapsedSample.append(avQualR1)
            if headerLen == 7:
                try:
                    avQualR2 = [ float(qual[6]) for qual in dataDic[entry]  ]
                    avQualR2 = round( ( sum(avQualR2) / len(avQualR2) ), 2)
                except:
                    avQualR2 = 0
                collapsedSample.append(avQualR2)
            dataSquashed.append(collapsedSample)
    t = Table(dataSquashed, style=[
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
        Paragraph("The sample barcode added by the sequencing facility (or you, if you created the libraries yourself).",
            stylesheet['BodyText'])])
    key.append([Paragraph("Barcode ID(s)",
            stylesheet['BodyText']),
        Paragraph("The ID of the barcode added by the sequencing facility (or you, if you created the libraries yourself).",
            stylesheet['BodyText'])])
    key.append([Paragraph("Lane", 
            stylesheet['BodyText']),
        Paragraph("The lane number on the flow cell",
            stylesheet['BodyText'])])
    key.append([Paragraph("# Reads", 
            stylesheet['BodyText']),
        Paragraph("The number of reads in a given file. For paired-end datasets, this is equivalent to the number of fragments sequenced, rather than summing the counts for read #1 and read #2. Note that this includes only reads passing the quality filter.",
            stylesheet['BodyText'])])
    key.append([Paragraph("% Bases >= Q30 Read #1", 
            stylesheet['BodyText']),
        Paragraph("The percentage of bases in read #1 of a pair having a Phred-scaled score of at least 30, meaning that the 0.1% or less chance that they're incorrect.",
            stylesheet['BodyText'])])
#    key.append([Paragraph("Ave. Qual. Read #1", 
#            stylesheet['BodyText']),
#        Paragraph("The average Phred-scaled base quality of bases in read #1 of a pair. This number of -10*log10(Probability that the call is incorrect). In other words, if a call is 100% likely to be wrong, the score is 0 (or 10 for 10% likelihood, 20 for 1% likelihood, etc.).",
#            stylesheet['BodyText'])])
    key.append([Paragraph("% Bases >= Q30 Read #2", 
            stylesheet['BodyText']),
        Paragraph("Identical to '% Bases >= Q30 Read #1', but for read #2 of a pair.",
            stylesheet['BodyText'])])
#    key.append([Paragraph("Ave. Qual. Read #2", 
#            stylesheet['BodyText']),
#        Paragraph("Identical to 'Ave. Qual. Read #1', but for read #1 of a pair.",
#            stylesheet['BodyText'])])
#    key.append([Paragraph("# Reads", 
#            stylesheet['BodyText']),
#        Paragraph("Identical to '% Bases >= Q30 Read #1', but for single-end datasets.",
#            stylesheet['BodyText'])])
#    key.append([Paragraph("Ave. Qual.", 
#            stylesheet['BodyText']),
#        Paragraph("Identical to 'Ave. Qual. Read #1', but for single-end datasets.",
#            stylesheet['BodyText'])])
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
    fqs = glob.glob("%s/%s%s/Project_%s/*/*.png" % (
        config.get("Paths","outputDir"),
        config.get("Options","runID"),
        lanes,
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
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes = "_lanes{}".format(lanes)

    try :
        tree = ET.parse("%s/%s%s/Stats/ConversionStats.xml" % (config.get("Paths","outputDir"),config.get("Options","runID"), lanes))
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
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes = "_lanes{}".format(lanes)

    message = "Flow cell: %s%s\n" % (config.get("Options","runID"), lanes)
    message += "Run time: %s\n" % runTime
    message += "Data transfer: %s\n" % transferTime
    message += msg

    msg = MIMEText(message)
    msg['Subject'] = "[bcl2fastq_pipeline] %s%s processed" % (config.get("Options","runID"), lanes)
    msg['From'] = config.get("Email","fromAddress")
    msg['To'] = config.get("Email","finishedTo")

    s = smtplib.SMTP(config.get("Email","host"))
    s.send_message(msg)
    s.quit()

def jsonParkour(config, msg):
    d = dict()
    d['flowcell_id'] = config.get("Options", "runID").split("_")[3][1:]
    if "-" in d['flowcell_id']:
        # MiSeq runs need to have 000000- stripped from the front
        d['flowcell_id'] = d['flowcell_id'].split('-')[-1]

    # For each lane:
    # - % clusters passing filter (OK)
    # - number reads passing filter (OK)
    # - undetermined indices (%) (OK)
    # - % PhiX (unknowable?)
    # - read 1 %bases >=Q30 (OK)
    # - read 2 %bases >=Q30 (OK)

    laneDict = dict()
    for line in msg.split("\n"):
        if line.startswith("Lane\t# Clusters (% pass)"):
            continue
        if not line.startswith("Lane"):
            continue
        if "had undetermined" in line:
            cols = line.split(" ")
            lane = "Lane {}".format(cols[1][:-1])
            if lane not in laneDict:
                laneDict[lane] = dict()
            percent = cols[-1].strip("(").strip(")")
            readsPF = cols[4]
            laneDict[lane]["undetermined_indices"] = percent
            laneDict[lane]["reads_pf"] = readsPF
        else:
            cols = line.split("\t")
            lane = cols[0]
            clusterPF = cols[1].split()[-1].strip("(").strip(")").strip("%")
            cols2 = cols[2].split("/")
            read1q30 = cols2[0][:-1]
            read2q30 = None
            if len(cols2) == 2:
                read2q30 = cols2[1][:-1]
            # if there's only one sample, there will be no undetermined indices
            if lane not in laneDict:
                readsPF = cols[1].split(0)
                laneDict[lane] = {'undetermined_indices': 0.0,
                                  'reads_pf': readsPF}
            laneDict[lane]["read_1"] = read1q30
            laneDict[lane]["read_2"] = read2q30
            laneDict[lane]["cluster_pf"] = clusterPF
            laneDict[lane]["name"] = lane

    d['matrix'] = json.dumps(list(laneDict.values()))
    res = requests.post(config.get("parkour", "URL"), auth=(config.get("parkour", "user"), config.get("parkour", "password")), data=d)
    if res.status_code == 200:
        return "\nParkour: updated\n"
    else:
        return "\nParkour: parkour returned {}, status {}\nSent: {}".format(res.text, res.status_code, d)
