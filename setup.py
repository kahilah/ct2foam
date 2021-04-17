import setuptools

__title__ = 'ct2foam'
__version__ = '1.0.0'
__author__ = 'Heikki Kahila'
__license__ = 'XXX'
__copyright__ = 'Copyright 2021 by Heikki Kahila'


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ct2foam",
    version="0.0.1",
    author="Heikki Kahila",
    author_email="heikki.kahila@wartsila.com",
    description="""A package including utilities for generating OpenFoam dictionaries for thermodynamic and transport entries. 
    In addition, utilities to use pyJac libraries within the DLBFoam is provided.""",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=["ct2foam", "ct2foam.thermo_transport", "ct2foam.scripts"],
    include_package_data=True,    
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: XXX License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        'cvxopt>=1.2.5',
        'numpy>=1.19.2',
        'scipy>=1.6.2',
        'matplotlib>=3.3.4'
    ],
    # files treated as global executables
    scripts=["pyjac2foam/pyjac2foam.sh"],
    entry_points = {
        'console_scripts': 
            ['ct2foam=ct2foam.scripts.mech2foam:main',
             'ctmix2foam=ct2foam.scripts.mix2foam:main'],
    }
)