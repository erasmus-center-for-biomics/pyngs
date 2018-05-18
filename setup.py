#!/bin/env python

#
from setuptools import setup

# set the variables for the installation
DEPENDENCIES = ["setuptools"]
VERSION = "Undefined"
DOC = ""
ADD_TO_DOC = True
with open("pyngs/__init__.py") as infile:
    for line in infile:
        line = line.rstrip()
        if line.startswith("VERSION:"):
            VERSION = line.split(":", 1)[1]
        if "\"\"\"" in line:
            ADD_TO_DOC = True
        elif "\"\"\"" in line and ADD_TO_DOC:
            ADD_TO_DOC = False

        if ADD_TO_DOC:
            DOC += "{line}\n".format(line=line.replace("\"\"\"", ""))

# prepare the setup
setup(
    name="pyngs",
    packages=["pyngs", "pyngs.tests"],
    scripts=["combine_flagstat"],
    author="r.w.w.brouwer",
    author_email="@erasmusmc.nl",
    description="A NGS library for Python",
    long_description=DOC,
    test_suite='pyngs.tests',
    install_requires=DEPENDENCIES,
    entry_points={
    },
    url="",
    version=VERSION,
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    keywords='bioinformatics',
    use_2to3=False,
    include_package_data=True,
    package_data={
        
    }
)
