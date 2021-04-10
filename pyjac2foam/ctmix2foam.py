import argparse
import cantera as ct
import numpy as np
from thermo_transport import  ct_properties
from thermo_transport import  transport_fitter as tr_fitter
from thermo_transport import  thermo_fitter as th_fitter
from thermo_transport import  ct2foam_utils as utils


if(__name__ == '__main__'):
    parser = argparse.ArgumentParser(description='Convert/refit cantera .cti -based transport and thermodynamic data into OpenFOAM format.')
    parser.add_argument('-i','--input', type=str, help='Mechanism path.', default=None, required=True)
    parser.add_argument('-o','--output', type=str, help='Output directory path.', default=None, required=False)
    parser.add_argument('-n','--name', type=str, help='Mixture name, e.g. "air".', default=None, required=True)
    parser.add_argument('-m','--mixture', type=str, help='Molecular mixture ratio in cantera style: "O2:1, N2:3.76" ', default=None, required=True)
    parser.add_argument('-T','--Tmid', type=float, help='Common temperature for NASA-7 thermodynamical fits.', default=1000.0, required=False)
    parser.add_argument('-p','--plot', type=bool, help='Generate plots when available.', default=False, required=False)
    args = parser.parse_args()

    mech = args.input
    mix_name = args.name
    ct_mixture = args.mixture
    nasa7_Tmid = args.Tmid

    thermo = ct_properties.ctThermoTransport(mech)
    mix_data = thermo.evaluate_mixture_properties(mix_name, ct_mixture)

    print("\nFitting transport properties:")
    As, Ts, std_err = tr_fitter.fit_sutherland(thermo.T, mix_data[0]["mu"])
    poly_mu, poly_kappa = tr_fitter.fit_polynomial(thermo.T, mix_data[0]["mu"], mix_data[0]["kappa"])

    success = utils.transport_fit_quality(
        np.array([mix_name]), thermo.T,  mix_data[0]["mu"],  mix_data[0]["kappa"],
        poly_mu, poly_kappa, As, Ts, std_err, mix_data[1]["cv_mole"], mix_data[0]["W"],
        plot=args.plot, rel_tol_sutherland=2e-2, rel_tol_Euken=5e-2, rel_tol_poly=5e-3
        )
    if(not success):
        print("\nMixture transport fit has failed.\n")
    else:
        print("\tSuccess.")

    print("\nFitting thermodynamic properties:")
    nasa_coeffs_lo, nasa_coeffs_hi = utils.fit_mixture_thermo(
        thermo.T, nasa7_Tmid, mix_data[1]["cp"], mix_data[1]["h"], mix_data[1]["s"], 
        mix_data[1]["cp0_over_R"], mix_data[1]["dhf_over_R"], mix_data[1]["s0_over_R"]
        )
    success = utils.nasa7_fit_quality(
        np.array([mix_name]), thermo.T, mix_data[1]["cp"], mix_data[1]["h"], mix_data[1]["s"],
        nasa7_Tmid, nasa_coeffs_lo, nasa_coeffs_hi, plot=args.plot)
    if(not success):
        print("\tThe NASA7 polynomial fit has failed.\n")
    else:
        print("\tSuccess.")
    
    
