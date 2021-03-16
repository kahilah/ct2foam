import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()
with open("requirements.txt", "r") as fh:
    requirements = [line.strip() for line in fh]

setuptools.setup(
    name="pyjac2foam",
    version="0.0.1",
    author="Heikki Kahila",
    author_email="heikki.kahila@wartsila.com",
    description="""A package including utilities for generating pyJac libraries for the use in DLBFoam.""",
    long_description="""A Python library providing utilities for exporting pyJac [Niemeyer et al. 2017] 
    generated C source files into OpenFoam format. In particular, this package is designed together 
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
    #scripts=['pywice/xdmfSurfaceWriter.py']
)