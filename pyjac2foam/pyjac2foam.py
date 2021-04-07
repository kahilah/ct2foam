import argparse
import sys
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
from thermo_transport import  ct_properties
from thermo_transport import  transport_fitter as tr_fitter
from thermo_transport import  thermo_fitter as th_fitter


def fit_transport(thermo):
    print("\nFitting transport properties:")
    As, Ts, std_err = tr_fitter.fit_sutherland(thermo.T, thermo.mu)
    err_mu_sutherland, err_kappa_sutherland = tr_fitter.error_sutherland(thermo.mu, thermo.kappa, thermo.T, As, Ts, thermo.cv_mass, thermo.cv_mole, thermo.R)

    poly_mu, poly_kappa = tr_fitter.fit_polynomial(thermo.T, thermo.mu, thermo.kappa)
    err_mu_poly, err_kappa_poly = tr_fitter.error_polynomial(thermo.mu, thermo.kappa, poly_mu, poly_kappa, thermo.T)
    print("\tSuccess.")


def ct_thermo_consistency(thermo, sp_i):
    nasa7_coeffs = thermo.get_nasa7_coeffs(sp_i)
    cp_over_R, h_over_RT, s_over_R = thermo.nasa_normalisation(thermo.T, thermo.cp, thermo.h, thermo.s)
    i = thermo.gas.species_index(sp_i)
    consistent_data = th_fitter.consistency(thermo.T, nasa7_coeffs[0], nasa7_coeffs[8:15], nasa7_coeffs[1:8], cp_over_R[i,:], h_over_RT[i,:], s_over_R[i,:])
    if(not consistent_data):
        print("\tError, " + sp_i + " : there is a mismatch between Cantera thermo data and NASA polynomials.")
        print("\tExiting...")
        sys.exit(1)
    else:
        return True

def ct_thermo_continuity(thermo, sp_i):
    nasa7_coeffs = thermo.get_nasa7_coeffs(sp_i)
    continuous_fit = th_fitter.continuity(nasa7_coeffs[0], nasa7_coeffs[8:15], nasa7_coeffs[1:8])
    if(not continuous_fit):
        print("\t" + sp_i + " : original NASA polynomial is not continuous, a re-fit is required")
        return False
    else:
        return True


def thermo_fit_error(T, cp, h, s, Tmid, coeffs_lo, coeffs_hi, rel_tol=5e-3):

    cp_fit = th_fitter.cp_nasa7(T, Tmid, coeffs_lo, coeffs_hi)*ct.gas_constant
    dcp = np.abs(cp - cp_fit)
    err_cp = np.linalg.norm(dcp)/np.linalg.norm(cp)

    h_fit = th_fitter.h_nasa7(T, Tmid, coeffs_lo, coeffs_hi)*ct.gas_constant*T
    dh = np.abs(h - h_fit)
    err_h = np.linalg.norm(dh)/np.linalg.norm(h)

    s_fit = th_fitter.s_nasa7(T, Tmid, coeffs_lo, coeffs_hi)*ct.gas_constant
    ds = np.abs(s - s_fit)
    err_s = np.linalg.norm(ds)/np.linalg.norm(s)

    failed = err_cp > rel_tol or err_h > rel_tol or err_s > rel_tol
    if(failed):
        print("\tWarning : fit error is larger than tolerance " + repr(rel_tol))

    continuous_fit = th_fitter.continuity(Tmid, coeffs_lo, coeffs_hi)
    if(not continuous_fit):
            print("\tWarning : the new NASA polynomial fit is not continuous.")

    if( (not failed) and continuous_fit):
        return True
    else:
        return False


def plot_thermo_fits(sp_i, T, cp_over_R, h_over_RT, s_over_R, Tmid, coeffs_lo, coeffs_hi):
    
    cp_fit = th_fitter.cp_nasa7(T, Tmid, coeffs_lo, coeffs_hi)
    h_fit = th_fitter.h_nasa7(T, Tmid, coeffs_lo, coeffs_hi)
    s_fit = th_fitter.s_nasa7(T, Tmid, coeffs_lo, coeffs_hi)

    fig=plt.figure(num=2,figsize=(7.5,10))
    ax1 = plt.subplot(311)
    plt.plot(T, cp_over_R, '-', color='r', label='Orig.')
    plt.plot(T, cp_fit, '--', color='b', label='Fit')
    ax1.set_ylabel('cp/R')
    plt.legend(loc=4)

    ax2 = plt.subplot(312)
    plt.plot(T, h_over_RT, '-', color='r')
    plt.plot(T, h_fit, '--', color='b')
    ax2.set_ylabel('h/RT')

    ax3 = plt.subplot(313)
    plt.plot(T, s_over_R, '-', color='r')
    plt.plot(T, s_fit, '--', color='b')
    ax3.set_ylabel('s/R')
    
    fig.savefig('Figures/'+sp_i+'_thermo.png', bbox_inches='tight')
    fig.clf()


def refit_ct_thermo(thermo, Tmid):
    
    print("\nFitting thermodynamic properties:\n")
    
    Tci = np.where(thermo.T==1000.0)[0][0]

    for sp_i in thermo.gas.species_names:
        
        i = thermo.gas.species_index(sp_i)
        cp_over_R, h_over_RT, s_over_R = thermo.nasa_normalisation(thermo.T, thermo.cp, thermo.h, thermo.s)

        if(thermo.is_nasa7(sp_i)):

            consistent = ct_thermo_consistency(thermo, sp_i)
            continuos = ct_thermo_continuity(thermo, sp_i)
            equal_Tmid = (abs(thermo.get_nasa7_coeffs(sp_i)[0] - Tmid)/Tmid < 1e-12)

            if(consistent and continuos and equal_Tmid):
                coeffs = thermo.get_nasa7_coeffs(sp_i)
                coeffs_lo, coeffs_hi = coeffs[8:15], coeffs[1:8]
                continue
            elif(consistent and continuos):
                coeffs_lo, coeffs_hi = th_fitter.fit_nasapolys_cp(thermo.T, Tci, cp_over_R[i,:], thermo.cp0_over_R[i], thermo.dhf_over_R[i], thermo.s0_over_R[i])
            else:
                coeffs_lo, coeffs_hi = th_fitter.fit_nasapolys_full(thermo.T, Tci, cp_over_R[i,:], h_over_RT[i,:], s_over_R[i,:], thermo.cp0_over_R[i], thermo.dhf_over_R[i], thermo.s0_over_R[i])

        else:
            type_name = thermo.get_thermo_fit_type(sp_i)
            print("\tWarning, " + sp_i + " : "+ type_name + " thermo type will be fitted to NASA 7-coefficient format.")
            print("\tUser is recommended to visually check the fit quality.")           
            coeffs_lo, coeffs_hi = th_fitter.fit_nasapolys_full(thermo.T, Tci, cp_over_R[i,:], h_over_RT[i,:], s_over_R[i,:], thermo.cp0_over_R[i], thermo.dhf_over_R[i], thermo.s0_over_R[i])

        success = thermo_fit_error(thermo.T, thermo.cp[i,:], thermo.h[i,:], thermo.s[i,:], Tmid, coeffs_lo, coeffs_hi)
        #plot_thermo_fits(sp_i, thermo.T, cp_over_R[i,:], h_over_RT[i,:], s_over_R[i,:], Tmid, coeffs_lo, coeffs_hi)

        if(not success):
            print("\tFailure.\n")



if(__name__ == '__main__'):
    parser = argparse.ArgumentParser(description='Solve adiabatic free flame over multiple initial conditions.')
    parser.add_argument('-n','--ncores',  type=int, help='Number available cores for parallel execution.', default=1, required=False)
    parser.add_argument('-id','--rankID', type=int, help='Individual rank/task id.', default=0, required=False)
    args = parser.parse_args()


    
    thermo = ct_properties.ctThermoTransport("h2o2.yaml")
    thermo = ct_properties.ctThermoTransport("gri30.yaml")

    # test NASA9 polynomial based thermo
    #thermo = ct_properties.ctThermoTransport("h2o2_mod.yaml")
    thermo.evaluate_properties()
    
    fit_transport(thermo)

    Tmid = 1000.0
    refit_ct_thermo(thermo, Tmid)

    # mixture example
    #mix_data = thermo.evaluate_mixture_properties("air", "O2:1,AR:3.67")
    #As, Ts, std_err = tr_fitter.fit_sutherland(thermo.T, mix_data[0]["mu"])
    #err_mu_sutherland, err_kappa_sutherland = tr_fitter.error_sutherland(mix_data[0]["mu"], mix_data[0]["kappa"], thermo.T, As, Ts, mix_data[1]["cv_mass"], mix_data[1]["cv_mole"], thermo.R)

    '''
    Tc_i = np.where(thermo.T==1000.0)[0][0]
    # - enthalpy of formation and entropy at standard conditions are required by definition
    cp0_over_R =  thermo.gas.species(thermo.gas.species_index("H2")).thermo.cp(298.15)/thermo.R # this calls molar cp <==> consistent
    dhf_by_R =  thermo.gas.species(thermo.gas.species_index("H2")).thermo.h(298.15)/thermo.R
    s0_by_R =  thermo.gas.species(thermo.gas.species_index("H2")).thermo.s(298.15)/thermo.R
    cp_over_R = thermo.cp[0,:]/thermo.R
    h_over_RT = thermo.h[0,:]/(thermo.T*thermo.R)
    s_over_R = thermo.s[0,:]/thermo.R
    coeffs = nasa_polys.fit_nasapolys_cp(thermo.T,Tc_i,cp_over_R,cp0_over_R,dhf_by_R,s0_by_R) 
    '''
