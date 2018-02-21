#!/usr/bin/env python3
from setuptools import setup, Extension
from distutils import sysconfig

srcs = ["pyBarcodes.c"]
libs=["z"]
if sysconfig.get_config_vars('BLDLIBRARY') is not None:
    #Note the "-l" prefix!
    for e in sysconfig.get_config_vars('BLDLIBRARY')[0].split():
        if e[0:2] == "-l":
            libs.append(e[2:])
elif(sys.version_info[0] >= 3 and sys.version_info[1] >= 3) :
    libs.append("python%i.%im" % (sys.version_info[0], sys.version_info[1]))
else :
    libs.append("python%i.%i" % (sys.version_info[0], sys.version_info[1]))

additional_libs = [sysconfig.get_config_var("LIBDIR"), sysconfig.get_config_var("LIBPL")]

module1 = Extension('pyBarcodes',
                    sources = srcs,
                    libraries = libs,
                    library_dirs = additional_libs)

setup(name = 'bcl2fastq_pipeline',
       version = '0.3.1',
       description = 'bcl2fastq_pipeline',
       author = "Devon P. Ryan",
       author_email = "ryan@ie-freiburg.mpg.de",
       scripts = ['bin/bfq.py'],
       packages = ['bcl2fastq_pipeline'],
       include_package_data = False,
       install_requires = ['configparser',
                           'reportlab',
                           'numpy',
                           'matplotlib',
                           'bioblend'],
       ext_modules = [module1])

