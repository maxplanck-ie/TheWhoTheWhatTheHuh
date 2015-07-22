#!/usr/bin/env python3
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'bcl2fastq_pipeline',
    'author': 'Devon Ryan',
    'url': 'github.com/',
    'author_email': 'ryan@ie-freiburg.mpg.de.',
    'version': '0.2',
    'packages': ['bcl2fastq_pipeline'],
    'scripts': ['bin/bfq.py'],
    'name': 'bcl2fastq_pipeline'
}

setup(**config)
