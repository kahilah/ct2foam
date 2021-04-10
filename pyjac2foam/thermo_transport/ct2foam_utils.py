import sys
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
from thermo_transport import  transport_fitter as tr_fitter
from thermo_transport import  thermo_fitter as th_fitter



def nasa_normalisation(T, cp_mole, h_mole, s_mole):
    """
    Divide by gas constant (and T) according to NASA JANAF definitions.
    Assuming data structure [M,] or [M,N] shaped data, where M is the number of species and N is data size      
    """
    cp_over_R = cp_mole / ct.gas_constant
    h_over_RT = h_mole[None,:] / (T * ct.gas_constant)
    s_over_R  = s_mole / ct.gas_constant

    return cp_over_R, h_over_RT[0], s_over_R


def fit_ct_transport(thermo):
    As, Ts, std_err = tr_fitter.fit_sutherland(thermo.T, thermo.mu)
    poly_mu, poly_kappa = tr_fitter.fit_polynomial(thermo.T, thermo.mu, thermo.kappa)
    return As, Ts, std_err, poly_mu, poly_kappa


def transport_fit_quality(names, T, mu, kappa, poly_mu, poly_kappa, As, Ts, suth_std_err, cv_mole, W, plot=True, rel_tol_sutherland=2e-2, rel_tol_Euken=2e-2, rel_tol_poly=5e-3):

    N = mu.shape[0]
    success = True
    err_mu_sutherland, err_kappa_sutherland = tr_fitter.error_sutherland(mu, kappa, T, As, Ts, cv_mole, W, ct.gas_constant)
    err_mu_poly, err_kappa_poly = tr_fitter.error_polynomial(mu, kappa, poly_mu, poly_kappa, T)

    for i in range(N):
        
        fit_ok = (err_mu_sutherland[i] < rel_tol_sutherland) and (err_kappa_sutherland[i] < rel_tol_Euken) and (err_mu_poly[i] < rel_tol_poly) and (err_kappa_poly[i] < rel_tol_poly)
        if(not fit_ok):
            success = False
            print("\n\tWarning : " + repr(names[i]) + " fit error is larger than tolerance.")
            print("\tSutherland mu error: " + repr(err_mu_sutherland[i]))
            print("\tEuken kappa error: " + repr(err_kappa_sutherland[i]))
            print("\tSutherland STD: " + repr(suth_std_err[i]))
            print("\tPolynomial mu error: " + repr(err_mu_poly[i]))           
            print("\tPolynomial kappa error: " + repr(err_kappa_poly[i]))           

            if(plot):
                plot_transport_comparison(names[i], T, mu[i,:], kappa[i,:], poly_mu[:,i], poly_kappa[:,i], As[i], Ts[i], cv_mole[i], W[i])

    return success


def plot_transport_comparison(name, T, mu, kappa, poly_coeffs_mu, poly_coeffs_kappa, As, Ts, cv_mole, W):
    
    mu_sutherland = tr_fitter.sutherland(T, As, Ts)
    kappa_euken = tr_fitter.euken(mu_sutherland, cv_mole, W, ct.gas_constant)
    mu_poly, kappa_poly = tr_fitter.eval_polynomial(poly_coeffs_mu, poly_coeffs_kappa, T, 3, True)

    fig=plt.figure(num=2,figsize=(7.5,10))
    ax1 = plt.subplot(211)
    plt.plot(T, mu, '-', color='r', label='Orig.')
    plt.plot(T, mu_sutherland, '--', color='b', label='Sutherland')
    plt.plot(T, mu_poly, '--', color='k', label='Polynomial')
    ax1.set_ylabel('$\mu$')
    ax1.set_xlabel('$T$[K]')
    plt.legend(loc=4)

    ax2 = plt.subplot(212)
    plt.plot(T, kappa, '-', color='r', label='Orig.')
    plt.plot(T, kappa_euken, '--', color='b', label='Euken')
    plt.plot(T, kappa_poly, '--', color='k', label='Polynomial')
    ax2.set_ylabel('$\kappa$')
    ax2.set_xlabel('$T$[K]')
    plt.legend(loc=4)

    fig.savefig('Figures/'+name+'_transport.png', bbox_inches='tight')
    fig.clf()


def ct_thermo_consistency(thermo, sp_i):
    nasa7_coeffs = thermo.get_nasa7_coeffs(sp_i)
    cp_over_R, h_over_RT, s_over_R = nasa_normalisation(thermo.T, thermo.cp, thermo.h, thermo.s)
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
        return False
    else:
        return True


def plot_nasa7_comparison(name, T, cp_over_R, h_over_RT, s_over_R, Tmid, coeffs_lo, coeffs_hi):
    
    cp_fit = th_fitter.cp_nasa7(T, Tmid, coeffs_lo, coeffs_hi)
    h_fit = th_fitter.h_nasa7(T, Tmid, coeffs_lo, coeffs_hi)
    s_fit = th_fitter.s_nasa7(T, Tmid, coeffs_lo, coeffs_hi)

    fig=plt.figure(num=2,figsize=(7.5,10))
    ax1 = plt.subplot(311)
    plt.plot(T, cp_over_R, '-', color='r', label='Orig.')
    plt.plot(T, cp_fit, '--', color='b', label='Fit')
    ax1.set_ylabel('cp/R')
    plt.legend(loc=4)
    ax1.set_xlabel('$T$[K]')

    ax2 = plt.subplot(312)
    plt.plot(T, h_over_RT, '-', color='r')
    plt.plot(T, h_fit, '--', color='b')
    ax2.set_ylabel('h/RT')
    ax2.set_xlabel('$T$[K]')

    ax3 = plt.subplot(313)
    plt.plot(T, s_over_R, '-', color='r')
    plt.plot(T, s_fit, '--', color='b')
    ax3.set_ylabel('s/R')
    ax3.set_xlabel('$T$[K]')

    fig.savefig('Figures/'+name+'_thermo.png', bbox_inches='tight')
    fig.clf()


def refit_ct_thermo(thermo, Tmid):
       
    nasa_coeffs_lo = np.zeros((thermo.gas.n_species, 7))
    nasa_coeffs_hi = np.zeros((thermo.gas.n_species, 7))
    Tci = np.where(thermo.T==Tmid)[0][0]

    for sp_i in thermo.gas.species_names:
        
        i = thermo.gas.species_index(sp_i)
        cp_over_R, h_over_RT, s_over_R = nasa_normalisation(thermo.T, thermo.cp, thermo.h, thermo.s)
        if(thermo.is_nasa7(sp_i)):

            ct_thermo_consistency(thermo, sp_i)
            continuos = ct_thermo_continuity(thermo, sp_i)
            equal_Tmid = (abs(thermo.get_nasa7_coeffs(sp_i)[0] - Tmid)/Tmid < 1e-12)

            if(continuos and equal_Tmid):
                coeffs = thermo.get_nasa7_coeffs(sp_i)
                coeffs_lo, coeffs_hi = coeffs[8:15], coeffs[1:8]
            elif(continuos):
                coeffs_lo, coeffs_hi = th_fitter.fit_nasapolys_cp(thermo.T, Tci, cp_over_R[i,:], thermo.cp0_over_R[i], thermo.dhf_over_R[i], thermo.s0_over_R[i])
            else:
                coeffs_lo, coeffs_hi = th_fitter.fit_nasapolys_full(thermo.T, Tci, cp_over_R[i,:], h_over_RT[i,:], s_over_R[i,:], thermo.cp0_over_R[i], thermo.dhf_over_R[i], thermo.s0_over_R[i])
        else:
            type_name = thermo.get_thermo_fit_type(sp_i)
            print("\tWarning, " + sp_i + " : "+ type_name + " thermo type will be fitted to NASA 7-coefficient format.")
            print("\tUser is recommended to visually check the fit quality.")           
            coeffs_lo, coeffs_hi = th_fitter.fit_nasapolys_full(thermo.T, Tci, cp_over_R[i,:], h_over_RT[i,:], s_over_R[i,:], thermo.cp0_over_R[i], thermo.dhf_over_R[i], thermo.s0_over_R[i])

        nasa_coeffs_lo[i] = coeffs_lo
        nasa_coeffs_hi[i] = coeffs_hi

    return nasa_coeffs_lo, nasa_coeffs_hi


def fit_mixture_thermo(T, Tmid, cp, h, s, cp0_over_R, dhf_over_R, s0_over_R):
       
    Tci = np.where(T==Tmid)[0][0]
    cp_over_R, h_over_RT, s_over_R = nasa_normalisation(T, cp, h, s)
    coeffs_lo, coeffs_hi = th_fitter.fit_nasapolys_full(T, Tci, cp_over_R[0], h_over_RT[0], s_over_R[0], cp0_over_R, dhf_over_R, s0_over_R)
    return coeffs_lo, coeffs_hi


def nasa7_fit_quality(name, T, cp_ref, h_ref, s_ref, Tmid, nasa_coeffs_lo, nasa_coeffs_hi, plot=True):

    cp_ref = np.atleast_2d(cp_ref)
    h_ref = np.atleast_2d(h_ref)
    s_ref = np.atleast_2d(s_ref)
    nasa_coeffs_lo = np.atleast_2d(nasa_coeffs_lo)
    nasa_coeffs_hi = np.atleast_2d(nasa_coeffs_hi)    
    rel_tol = 5e-3
    N = cp_ref.shape[0]
    success = True
    for i in range(N):
        err_cp, err_h, err_s = th_fitter.error_nasa7(T, cp_ref[i,:], h_ref[i,:], s_ref[i,:], Tmid, nasa_coeffs_lo[i,:], nasa_coeffs_hi[i,:], ct.gas_constant)
        
        fit_ok = (err_cp < rel_tol) and (err_h < rel_tol) and (err_s < rel_tol)
        if(not fit_ok):
            success = False
            print("\tWarning : " + repr(name[i]) + " fit error is larger than tolerance " + repr(rel_tol))
            print("\tcp error: " + repr(err_cp))
            print("\th error: " + repr(err_h))
            print("\ts error: " + repr(err_s))           

        continuous_fit = th_fitter.continuity(Tmid, nasa_coeffs_lo[i,:], nasa_coeffs_hi[i,:])
        if(not continuous_fit):
            success = False
            print("\tWarning : " + repr(name[i]) + "the new NASA polynomial fit is not continuous.")

        if(plot):
            cp_over_R, h_over_RT, s_over_R = nasa_normalisation(T, cp_ref, h_ref, s_ref)
            plot_nasa7_comparison(name[i], T, cp_over_R[i,:], h_over_RT[i,:], s_over_R[i,:], Tmid, nasa_coeffs_lo[i,:], nasa_coeffs_hi[i,:])


    return success
