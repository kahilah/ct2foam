import setuptools

__title__ = 'ct2foam'
__version__ = '1.0.0'
__author__ = 'Heikki Kahila'
__license__ = 'XXX'
__copyright__ = 'Copyright 2021 by Heikki Kahila'


with open("README.md", "r") as fh:
    long_description = fh.read()
with open("requirements.txt", "r") as fh:
    requirements = [line.strip() for line in fh]

setuptools.setup(
    name="ct2foam",
    version="0.0.1",
    author="Heikki Kahila",
    author_email="heikki.kahila@wartsila.com",
    description="""A package including utilities for generating OpenFoam dictionaries for thermodynamic and transport entries. 
    In addition, utilities to use pyJac libraries within the DLBFoam is provided.""",
    long_description="""A library providing utilities for converting chemical mechanism data into  
    OpenFoam dictionaries for thermodynamic and transport entries. In addition, provides functionality to convert 
    pyJac [Niemeyer et al. 2017] generated C source files into a library usable together 
    with DLBFoam https://github.com/blttkgl/DLBFoam_dev.""",
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: XXX License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=requirements,
    # files treated as global executables
    #scripts=['ct2Foam/xdmfSurfaceWriter.py']
)