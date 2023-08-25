# ct2foam -- Convert Cantera based thermophysical data To OpenFoam format 
ct2foam enables user to generate OpenFoam dictionary entries for NASA-polynomials, Sutherland and other transport models required in <code>thermophysicalProperties</code> file. This python package utilises [Cantera](https://cantera.org/) to generate such entries for all species in a given chemical mechanism or for gaseous mixtures defined by the user. With the general functions provided in this package, it is also possible to create NASA-polynomial, Sutherland and other polynomial type fits for thermophysical data based on experiments.

Furthermore, this package supports the users of [DLBFoam](https://github.com/Aalto-CFD/DLBFoam) and [pyJac](https://github.com/SLACKHA/pyJac) by introducing an automated pyjac2foam script which builds a compilation environment for pyjac routines as well as utilises ct2foam to generate consistent thermodynamics dictionaries with pyJac. See pyjac2foam module instructions below for further information.

## Installation
Package relies on Cantera installation, for which the recommended installation principle is via [conda](https://cantera.org/install/index.html). Hence, the recommended and easiest installation path is via conda as well. In particular, we recommend using the Miniconda package manager for this. 

For users interested to use the pyjac2foam module, python version must be set <code>python<=3.6</code>. Otherwise, the latest releases can be utilised (tested up to python 3.8.8).

- To install in conda environment:
```
conda create --name ct2foam_env --channel cantera cantera python=3.6 numpy scipy matplotlib
cd my/installation/path
git clone git@github.com:kahilah/ct2foam.git
cd ct2foam
conda activate ct2foam_env
pip install .
```

- To install on ubuntu without conda:
    - follow the [instructions for cantera](https://cantera.org/install/ubuntu-install.html)
    - use <code>pip install .</code>, as above.

- To install in other environments:
    - See https://cantera.org/install/index.html and ensure that dependencies mentioned in <code>setup.py</code> are appreciated.

- Dependencies are installed automatically:
    - cantera, numpy, scipy, cvxopt, matplotlib
    - Distribution is tested with Cantera 2.5.1
    - <code>setup.py</code> lists the dependencies.

## Run tests:
```
cd ct2foam
python -m unittest discover
```
- Note that <code>ct2foam/test_data/OF_reference</code> includes C++ code requiring OpenFoam based compilation. However, this is mainly included here for development / testing purposes and standard user is not required to compile anything here.

## Running ct2foam:
- Installation generates two global executables: <code>ct2foam</code> and <code>ctmix2foam</code> which can be run as follows:
- <code>ct2foam --input h2o2.cti --output test_output --Tmid 1000.0 --plot </code>
    - outputs OpenFoam compatible thermodynamical and transport dictionary entries under output directory such that NASA-polynomials have a middle temperature of 1000 K. Thermodynamical fits with error larger than user given (or here default) tolerances, will be plotted under <code>output/Figures</code> directory. Supports
    input arguments <code>--Tlow</code> and <code>--Thigh</code> for low and high temperature limits of the
    NASA-polynomials, respectively.  
-  <code> ctmix2foam --input h2o2.cti --output test_output --name air --mixture "O2:1,AR:4" --Tmid 1000.0 --plot  </code>
    - outputs data for a mixture, defines similar to standard input format in Cantera.
- See <code>ct2foam -h</code> for help.
- Can be executed as a python modules:
    - python -m ct2foam.scripts.mech2foam
    - python -m ct2foam.scripts.mix2foam

- Outputs <code>log.txt</code> including all the relevant information related to the particular execution.
- Note that if you want plots for all species, you need to provide strict tolerance values for <code>transport_fit_quality</code> and <code>nasa7_fit_quality</code>

## pyjac2foam
- <code>pyjac2foam</code> module provides a consistent input of thermochemical data for users utilising [DLBFoam](https://github.com/Aalto-CFD/DLBFoam) and [pyJac](https://github.com/SLACKHA/pyJac) in OpenFoam environment. In particular, <code>pyjac2foam.sh</code> executable generates:
    - shared libray *.so, compiled with cmake.
    - consistent species ordering between OpenFOAM and pyjac inputs.
    - consistent <code>thermophysicalProperties</code>, <code>chemistsryProperties</code> and <code>controlDict</code> entries for an OpenFOAM case setup.


- Note that [pyJac](https://github.com/SLACKHA/pyJac) is not set as a requirement for this python package. To install, ensure that you have <code>python=3.6</code> installed. Then,
    - <code>conda install -c slackha pyjac</code> (or)
    - <code>pip install pyjac</code>

- Furthermore, you need cmake to use the automatic library compilation. 
- For a simple test:
```
cd pyjac2foam/test_data
./run_test.sh
```

- How-to:
    - Generate pyjac C-files (requires pyjac):
        - <code>python -m pyjac --lang c --last_species N2 --input my_mechanism.cti</code>
    - Compile + generate *.foam files under the same directory where pyjac output files are located.
        - <code>pyjac2foam -m my_mechanism.cti -i pyjac/path -c</code>
        - Example dictionary entries with correct absolute paths are printed after execution + written into file <code>foam/include_example.txt</code>

- Notes:    
    - note that compilation to a shared library is now with rather standard flags --> depending on the architecture, you may want to optimise CMakeLists
    - you can use  <code>pyjac2foam</code> for library compilation with <code>-c / --compile</code> flag, or copy <code>pyjac2foam/cmake_directives</code> directory, modify them and run <code>./runCmake.sh</code>


## Notes
- See TODO.txt for development to be considered in the future.

## Acknowledgements
- The constrained least-squares fitting procedure in <code>ct2foam/thermo_transport/lsqlin.py</code> is based on the work by Valeriy Vishnevskiy and Michael Hirsch, originally published in http://maggotroot.blogspot.ch/2013/11/constrained-linear-least-squares-in.html
