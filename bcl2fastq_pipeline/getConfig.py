import configparser
import os

def getConfig() :
    config = configparser.ConfigParser()
    config.read_file(open("%s/bcl2fastq.ini" % os.path.expanduser("~"),"r"))
    if("Paths" in config.sections()) :
        return config
    return None
