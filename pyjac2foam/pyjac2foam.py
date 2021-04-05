import argparse
import numpy as np
import os

from thermo_transport import  ct_properties
from thermo_transport import  transport_fitter as tr_fitter
from thermo_transport import  thermo_fitter as th_fitter

        

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


    nasa_polys = th_fitter.nasaPolynomials()
    test=nasa_polys.refit_nasapolys(thermo.nasa_coeffs, 1000.0)
    print(test)

    Tc_i = np.where(thermo.T==1000.0)[0][0]
    # - enthalpy of formation and entropy at standard conditions are required by definition
    cp0_over_R =  thermo.gas.species(thermo.gas.species_index("H2")).thermo.cp(298.15)/thermo.R # this calls molar cp <==> consistent
    dhf_by_R =  thermo.gas.species(thermo.gas.species_index("H2")).thermo.h(298.15)/thermo.R
    s0_by_R =  thermo.gas.species(thermo.gas.species_index("H2")).thermo.s(298.15)/thermo.R
    cp_over_R = thermo.cp[0,:]/thermo.R
    h_over_RT = thermo.h[0,:]/(thermo.T*thermo.R)
    s_over_R = thermo.s[0,:]/thermo.R
    coeffs = nasa_polys.fit_nasapolys_0(thermo.T,Tc_i,cp_over_R,h_over_RT,s_over_R,cp0_over_R,dhf_by_R,s0_by_R) 

