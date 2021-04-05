import argparse
import numpy as np
import os

from thermo_transport import  ct_properties
from thermo_transport import  transport_fitter as tr_fitter

        

if(__name__ == '__main__'):
    parser = argparse.ArgumentParser(description='Solve adiabatic free flame over multiple initial conditions.')
    parser.add_argument('-n','--ncores',  type=int, help='Number available cores for parallel execution.', default=1, required=False)
    parser.add_argument('-id','--rankID', type=int, help='Individual rank/task id.', default=0, required=False)
    args = parser.parse_args()


    
    thermo = ct_properties.ctThermoTransport("h2o2.cti")
    thermo.evaluate_properties()
    
    As, Ts, std_err = tr_fitter.fit_sutherland(thermo.T, thermo.mu)
    err_mu_sutherland, err_kappa_sutherland = tr_fitter.error_sutherland(thermo.mu, thermo.kappa, thermo.T, As, Ts, thermo.cv_mass, thermo.cv_mole, thermo.R)

    poly_mu, poly_kappa = tr_fitter.fit_polynomial(thermo.T, thermo.mu, thermo.kappa)
    err_mu_poly, err_kappa_poly = tr_fitter.error_polynomial(thermo.mu, thermo.kappa, poly_mu, poly_kappa, thermo.T)

    # mixture example
    #mix_data = thermo.evaluate_mixture_properties("air", "O2:1,AR:3.67")
    #As, Ts, std_err = tr_fitter.fit_sutherland(thermo.T, mix_data[0]["mu"])
    #err_mu_sutherland, err_kappa_sutherland = tr_fitter.error_sutherland(mix_data[0]["mu"], mix_data[0]["kappa"], thermo.T, As, Ts, mix_data[1]["cv_mass"], mix_data[1]["cv_mole"], thermo.R)


