import configparser
from bioblend.galaxy import GalaxyInstance
from bioblend.galaxy.objects import GalaxyInstance as GI
import sys
import os
import glob
import syslog


def getLibID(gi, libName):
    """
    Given a library name, like "foo bar", return the ID for the first library matching that
    """
    lib = gi.libraries.get_libraries(name=libName)
    if not lib or len(lib) == 0:
        libs = gi.libraries.get_libraries()
        # Handle differences in capitalization
        for l in libs:
            if l["name"].lower() == libName.lower():
                return l["id"]
        raise RuntimeError("No library named {}".format(libName))
    return lib[0]["id"]


def getFolderID(gi, libID, path):
    """
    Given a library ID (lib) and a path (e.g., "/foo/bar/sniggly3"), return the folder ID for sniggly3.
    
    If the path doesn't exist (in part or in total), then create it.
    """
    # Does the path already exist?
    folders = gi.libraries.get_folders(libID, name=path)
    if folders is not None and len(folders) > 0:
        return folders[0]["id"]
    
    # Get the closest base folder
    longest = gi.libraries.get_folders(libID, name="/")[0]
    folders = gi.libraries.get_folders(libID)
    for folder in folders:
        if path.startswith(folder["name"]):
            # Look for the longest pathname overlap. The next character in path MUST be "/",
            # since otherwise adding "/a/b/c2" when "/a/b/c" exists would result in "/a/b/c/2"!
            if len(folder["name"]) > len(longest["name"]) and path[len(folder["name"])] == "/":
                longest = folder

    # shorten the path name if relevant
    idx = len(longest["name"])
    pathLeft = path[idx:]
    if pathLeft.startswith("/"):
        pathLeft = pathLeft[1:]

    for fName in pathLeft.split("/"):
        gi.libraries.create_folder(libID, fName, base_folder_id=longest["id"]) # returns None
        if longest["name"] != "/":
            newFName = "{}/{}".format(longest["name"], fName)
        else:
            newFName = "/{}".format(fName)
        longest = gi.libraries.get_folders(libID, name=newFName)[0]

    return longest["id"]


def addFileToLibraryFolder(gi, libID, folderID, fileName, file_type='auto', dbkey='?', link=True, roles=''):
    """
    Link/copy "fname" into the library and folder specified by libID and folderID. These MUST exist.
    
    file_type, dbkey, and roles are pass through to upload_from_galaxy_filesystem().
    
    link must be True (default) or False. If it's False then files are copied in.
    
    This returns a dictionary with keys: name, url, and id (or presumably None on error).
    """
    if link == True:
        link_data_only = 'link_to_files'
    else:
        link_data_only = 'copy_files'

    rv = gi.libraries.upload_from_galaxy_filesystem(libID,
                                                    fileName,
                                                    folder_id=folderID,
                                                    file_type=file_type,
                                                    dbkey=dbkey,
                                                    link_data_only=link_data_only,
                                                    roles=roles)
    return rv


def addFileToLibrary(gi, libraryName, path, fileName, file_type='auto', dbkey='?', link=True, roles=''):
    """
    Add fileName to path in libraryName. gi is a GalaxyInstance
    
    The other parameters are passed to upload_from_galaxy_filesystem()
    """
    libID = getLibID(gi, libraryName)
    folderID = getFolderID(gi, libID, path)
    dataset = addFileToLibraryFolder(gi, libID, folderID, fileName, file_type=file_type, dbkey=dbkey, link=link, roles='')
    if not dataset:
        raise RuntimeError("Error adding '{}' to '{}'".format(fileName, path))
    return dataset


def getFileType(fName):
    """
    If the file name ends with .fastq.gz then return 'fastqsanger'. Otherwise, return 'auto'.
    """
    if fName.endswith(".fastq.gz") or fName.endswith(".fq.gz"):
        return "fastqsanger"
    return "auto"


def checkExists(datasets, folderID, fName):
    """
    Return true if there's a file named fName in a folder with folderID in side the library with ID libID. Otherwise, return False

    For fastq.gz files, the file name in Galaxy might be lacking the .gz extension when we check...
    """
    for dataset in datasets:
        if dataset.wrapped["folder_id"] == folderID:
            if dataset.name == os.path.basename(fName):
                return True
            if fName.endswith(".fastq.gz") and os.path.basename(fName)[:-3] == dataset.name:
                return True
    return False


def getFiles(d):
    """
    Given a directory, return a list of all files in all subdirectories
    """
    files = []
    for (dpath, dnames, filenames) in os.walk(d):
        for fname in filenames:
            files.append("{}/{}".format(dpath, fname))
    return files


def linkIntoGalaxy(config):
    """
    Given a project, connect to a running Galaxy instance and do the following:
      1) If the associated group exists, create a "sequencing data"/"flow cell"/"project"
         shared library folder structure
      2) Link all samples, preserving directory structures
      3) TODO: Handle Galaxy stripping off the .gz file extension
    """
    lanes = config.get("Options", "lanes")
    if lanes != "":
        lanes = "_lanes{}".format(lanes)

    url = config.get("Galaxy", "URL")
    userKey = config.get("Galaxy", "API key")
    gi = GalaxyInstance(url=url, key=userKey)
    gi2 = GI(url=url, api_key=userKey)

    message = "\n"
    projects = glob.glob("%s/%s%s/Project_*" % (config.get("Paths","outputDir"),config.get("Options","runID"), lanes))
    for project in projects :
        pname = project.split("/")[-1][8:]
        if(pname[0] == "A") :
            group = pname.split("_")[1].lower()
            basePath = "{}/{}/sequencing_data".format(config.get("Paths","groupDir"), group)

            try:
                libID = getLibID(gi, "{} sequencing runs".format(group))
            except:
                message += "\n{}\tNo sequencing data folder!".format(group)
                continue

            try:
                # memoize the datasets already in the library
                _ = gi.libraries.get_libraries(library_id=libID)[0]["name"]
                l = gi2.libraries.list(name=_)[0]
                currentDatasets = l.get_datasets()

                fileList = getFiles("{}/{}{}".format(basePath, config.get("Options", "runID"), lanes))
                for fName in fileList:
                    basePath2 = os.path.dirname(fName)[len(basePath):]
                    folderID = getFolderID(gi, libID, basePath2)
                    if checkExists(currentDatasets, folderID, fName):
                        syslog.syslog("[linkIntoGalaxy] Skipping {}, already added\n".format(os.path.basename(fName)))
                        continue
                    f = addFileToLibraryFolder(gi, libID, folderID, fName, file_type=getFileType(fName))
                message += "\n{}\tSuccessfully uploaded to Galaxy".format(pname)
            except:
                message += "\n{}\tError during Galaxy upload".format(pname)
    return message
