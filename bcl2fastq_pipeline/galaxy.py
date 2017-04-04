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


def getFolderID(gi, folderDict, libID, path):
    """
    Given a library ID (lib) and a path (e.g., "/foo/bar/sniggly3"), return the folder ID for sniggly3.
    
    If the path doesn't exist (in part or in total), then create it.
    """
    # Does the path already exist?
    if path in folderDict:
        return folderDict[path][0]
    
    # Get the closest base folder
    longest = None
    longestID = None
    dirs = [x for x in path.split("/") if x != ""]  # Strip out "" in the split
    l = len(dirs)
    while l >= 0:
        foo = "/{}".format("/".join(dirs[:l]))
        if foo in folderDict:
            longest = foo
            longestID = folderDict[foo][0]
            break
        l = l - 1
    
    
    pathLeft = path[len(longest):]
    for fName in pathLeft.strip("/").split("/"):
        rv = gi.libraries.create_folder(libID, fName, base_folder_id=longestID) # returns a dict!
        if longest != "/":
            longest = "{}/{}".format(longest, fName)
        else:
            longest = "/{}".format(fName)
        longestID = rv[0]['id']
        folderDict[longest] = [longestID]

    return longestID


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


def getFileType(fName):
    """
    If the file name ends with .fastq.gz then return 'fastqsanger'. Otherwise, return 'auto'.
    """
    if fName.endswith(".fastq.gz") or fName.endswith(".fq.gz"):
        return "fastqsanger"
    return "auto"


def checkExists(datasetDict, folderID, fName):
    """
    Return true if there's a file named fName in a folder with folderID in side the library with ID libID. Otherwise, return False
    """
    if fName in datasetDict:
        if folderID in datasetDict[fName]:
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

                currentDatasets = l.get_datasets() # A list of LibraryDataset() objects, which takes a while to return
                                                   # Has a folder_id and file_name
                                                   # Can the folders themselves be memoized? Then we can also hash folder IDs and add file names as a list.
                # Make a dict() out of the datasets
                dataDict = {}
                for dataset in currentDatasets:
                    if dataset.wrapped['file_name'] in dataDict:
                        dataDict[dataset.wrapped['file_name']].append(dataset.wrapped['folder_id'])
                    else:
                        dataDict[dataset.wrapped['file_name']] = [dataset.wrapped['folder_id']]
                folders = gi.libraries.get_folders(libID)

                # Make a dict() out of the folders
                folderDict = {}
                for folder in folders:
                    if folder['name'] == '/':
		        # Otherwise / becomes ""
                        folderDict[folder['name']] = [folder['id']]
                        continue
                    if folder['name'].rstrip("/") in folderDict:
                        folderDict[folder['name'].rstrip("/")].append(folder['id'])
                    else:
                        folderDict[folder['name'].rstrip("/")] = [folder['id']]

                fileList = getFiles("{}/{}{}".format(basePath, config.get("Options", "runID"), lanes))
                for fName in fileList:
                    basePath2 = os.path.dirname(fName)[len(basePath):]
                    folderID = getFolderID(gi, folderDict, libID, basePath2)
                    if checkExists(dataDict, folderID, fName):
                        syslog.syslog("[linkIntoGalaxy] Skipping {}, already added\n".format(os.path.basename(fName)))
                        continue
                    f = addFileToLibraryFolder(gi, libID, folderID, fName, file_type=getFileType(fName))
                    if fName in dataDict:
                        dataDict[fName].append(folderID)
                    else:
                        dataDict[fName] = [folderID]
                message += "\n{}\tSuccessfully uploaded to Galaxy".format(pname)
            except:
                message += "\n{}\tError during Galaxy upload".format(pname)
    return message
