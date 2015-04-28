"""
Misc. functions
"""

import configparser
import shutil
import smtplib
from email.mime.text import MIMEText

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
