#!/usr/bin/env python
import argparse
import glob
import shutil
import re
import os

parser = argparse.ArgumentParser(description="Rename files and folders in a given sequencing data directory appropriately. This is normally done by bfq.py.")
parser.add_argument("projects", help="Project/directory name to rename. You can specify more than one if you'd like", nargs="+")
args = parser.parse_args()

for project in args.projects:
    # Fix sample names
    fnames = glob.glob("{}/*/*.fastq.gz".format(project))
    for fname in fnames :
        idx = fname.rindex("_")
        fnew = fname[0:idx]
        fnew = re.sub(r"_S[0-9]+_R([12])$",r'_R\1', fnew) + ".fastq.gz"
        shutil.move(fname, fnew)

    # Fix library names
    snames = glob.glob("{}/*".format(project))
    for sname in snames :
        idx = sname.rindex("/")
        snew = "{}/Sample_{}".format(sname[:idx], sname[idx+1:])
        shutil.move(sname, snew)

    # Fix project name
    if not os.path.isdir(project):
        continue
    pnew = "Project_{}".format(project)
    shutil.move(project, pnew)
