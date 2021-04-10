import argparse
import cantera as ct
import os
from thermo_transport import  ct_properties
from thermo_transport import  transport_fitter as tr_fitter
from thermo_transport import  thermo_fitter as th_fitter
from thermo_transport import  ct2foam_utils as utils


def init_dirs(output_dir=None):
    if(output_dir==None):
        output_dir = os.getcwd()
        print("Using current directory as an output directory.")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir,"Figures"), exist_ok=True)

def init_dirs2(self):

    if(self.outputDir == None):
        self.outputDir = os.getcwd()
    thermoDir = os.path.join(self.outputDir, 'foamThermoTransport')
    try:
        os.makedirs(thermoDir, exist_ok=True)
    except OSError:  
        None
    thermoFile= os.path.join(thermoDir, 'thermo.foam')
    try:
        os.remove(thermoFile) 
    except OSError:
        None
    # Add the chem.foam include to the file first
    #with open(thermoFN,'a') as output:
    #    output.write('#include "chem.foam"')
    #    output.write('\n')
    #This would write the species names list in OF style. However, we want to include that into the reaction file in pyJac style. 
    #See makeMechFileForOF.sh at pyjac folder
    #wr.write_sp_list(gas.species_names,thermoFN)
    # --------------------------------------------- #


if(__name__ == '__main__'):
    parser = argparse.ArgumentParser(description='Convert/refit cantera .cti -based transport and thermodynamic data into OpenFOAM format.')
    parser.add_argument('-i','--input', type=str, help='Mechanism path.', default=None, required=True)
    parser.add_argument('-o','--output', type=str, help='Output directory path.', default=None, required=False)
    parser.add_argument('-T','--Tmid', type=float, help='Common temperature for NASA-7 thermodynamical fits.', default=1000.0, required=False)
    parser.add_argument('-p','--plot', type=bool, help='Generate plots when available.', default=False, required=False)
    args = parser.parse_args()


    init_dirs(args.output)

    thermo = ct_properties.ctThermoTransport(args.input, Tcommon=args.Tmid)
    thermo.evaluate_properties()
   
    print("\nFitting transport properties:")
    As, Ts, std_err, poly_mu, poly_kappa = utils.fit_ct_transport(thermo)
    success = utils.transport_fit_quality(thermo.gas.species_names, thermo.T, thermo.mu, thermo.kappa, poly_mu, poly_kappa, As, Ts, std_err, thermo.cv_mole, thermo.W, plot=args.plot)
    if(not success):
        print("\tSome transport fits have failed.\n")
    else:
        print("\tSuccess.")


    print("\nFitting thermodynamic properties:")
    nasa_coeffs = utils.refit_ct_thermo(thermo, thermo.Tcommon)
    success = utils.nasa7_fit_quality(thermo.gas.species_names, thermo.T, thermo.cp, thermo.h, thermo.s, thermo.Tcommon, nasa_coeffs[0], nasa_coeffs[1], plot=args.plot)
    if(not success):
        print("\tSome NASA7 polynomial fits have failed.\n")
    else:
        print("\tSuccess.")
    