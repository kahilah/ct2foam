import sys
import os
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt

from ct2foam.thermo_transport import warning_msg
from ct2foam.thermo_transport import transport_fitter as tr_fitter
from ct2foam.thermo_transport import thermo_fitter as th_fitter
from ct2foam.thermo_transport import foam_writer as writer



def nasa_normalisation(T, cp_mole, h_mole, s_mole, R=ct.gas_constant):
    """
    Divide by gas constant (and T) according to NASA JANAF definitions.
    Assuming data structure [M,] or [M,N] shaped data, where M is the number of species and N is data size      
    """
    cp_over_R = cp_mole / R
    h_over_RT = h_mole[None,:] / (T * R)
    s_over_R  = s_mole / R

    return cp_over_R, h_over_RT[0], s_over_R


def fit_ct_transport(data, poly_order=3):
    """
    data: a class object including T, mu, kappa entries.
    By default, a chemkin type order 3 fit is utilised.
    return: Sutherland, polynomial and log-polynomial coefficients.
    """
    N = data.mu.shape[0]
    As = np.zeros(N)
    Ts = np.zeros(N)
    std = np.zeros((N, 2))
    poly_mu = np.zeros((N, poly_order+1))
    poly_kappa = np.zeros((N, poly_order+1))
    logpoly_mu = np.zeros((N, poly_order+1))
    logpoly_kappa = np.zeros((N, poly_order+1))

    for i in range(N):
        As[i], Ts[i], std[i] = tr_fitter.fit_sutherland(data.T, data.mu[i,:])
        poly_mu[i], poly_kappa[i] = tr_fitter.fit_polynomial(data.T, data.mu[i,:], data.kappa[i,:], poly_order)
        logpoly_mu[i], logpoly_kappa[i] = tr_fitter.fit_log_polynomial(data.T, data.mu[i,:], data.kappa[i,:], poly_order)

    return As, Ts, std, poly_mu, poly_kappa, logpoly_mu, logpoly_kappa


def transport_fit_quality(data, transport_fits, output_dir, plot=True, rel_tol_sutherland=2e-2, rel_tol_Euken=2e-2, rel_tol_poly=5e-3):
    """
    Inspect the quality of transport fits againts the original data.
    data: a class object including species names, mu, kappa, cv_mole and W.
    transport_fits: input corresponds to the return tuple from fit_ct_tranport()
    return: boolean of success.
    """
    names, T, mu, kappa, cv_mole, W = data.names, data.T, data.mu, data.kappa, data.cv_mole, data.W
    As, Ts, std, poly_mu, poly_kappa, logpoly_mu, logpoly_kappa = transport_fits
    
    N = mu.shape[0]
    success = True

    for i in range(N):

        err_mu_sutherland, err_kappa_sutherland = tr_fitter.error_sutherland(mu[i,:], kappa[i,:], T, As[i], Ts[i], cv_mole[i], W[i], ct.gas_constant)
        err_mu_poly, err_kappa_poly = tr_fitter.error_polynomial(mu[i,:], kappa[i,:], poly_mu[i,:], poly_kappa[i,:], T)
        err_mu_logpoly, err_kappa_logpoly = tr_fitter.error_log_polynomial(mu[i,:], kappa[i,:], logpoly_mu[i,:], logpoly_kappa[i,:], T)

        sutherland_ok = (err_mu_sutherland < rel_tol_sutherland) and (err_kappa_sutherland < rel_tol_Euken) 
        polynomial_ok =  (err_mu_poly < rel_tol_poly) and (err_kappa_poly < rel_tol_poly)
        logpolynomial_ok =  (err_mu_logpoly < rel_tol_poly) and (err_kappa_logpoly < rel_tol_poly)

        if(not sutherland_ok):
            success = False
            warning_msg.warning_sutherland(output_dir, names[i], err_mu_sutherland, std[i], err_kappa_sutherland)

        if(not polynomial_ok):
            success = False
            warning_msg.warning_polynomial(output_dir, names[i], err_mu_sutherland, err_kappa_sutherland, False)

        if(not logpolynomial_ok):
            success = False
            warning_msg.warning_polynomial(output_dir, names[i], err_mu_sutherland, err_kappa_sutherland, True)

        if( (not (sutherland_ok and polynomial_ok and logpolynomial_ok)) and plot):
            fig_dir = os.path.join(output_dir, "Figures")
            os.makedirs(fig_dir, exist_ok=True)
            plot_transport_comparison(names[i], T, mu[i,:], kappa[i,:], poly_mu[i,:], poly_kappa[i,:], As[i], Ts[i], cv_mole[i], W[i], fig_dir)

    return success


def plot_transport_comparison(name, T, mu, kappa, poly_coeffs_mu, poly_coeffs_kappa, As, Ts, cv_mole, W, output_dir):
    
    mu_sutherland = tr_fitter.sutherland(T, As, Ts)
    kappa_euken = tr_fitter.euken(mu_sutherland, cv_mole, W, ct.gas_constant)
    mu_poly, kappa_poly = tr_fitter.eval_polynomial(poly_coeffs_mu, poly_coeffs_kappa, T)

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

    file = os.path.join(output_dir, (name+'_transport.png'))
    fig.savefig(file, bbox_inches='tight')
    fig.clf()


def ct_thermo_consistency(thermo, sp_i, output_dir):
    """
    Check if thermodynamical data is consistent with the original data (==original fit).
    sp_i: species name
    output_dir: output directory into which a potential warning log is saved.
    return: Boolean of success or sys.exit()
    """
    nasa7_coeffs = thermo.get_nasa7_coeffs(sp_i)
    cp_over_R, h_over_RT, s_over_R = nasa_normalisation(thermo.T, thermo.cp, thermo.h, thermo.s)
    i = thermo.gas.species_index(sp_i)
    consistent_data = th_fitter.consistency(thermo.T, nasa7_coeffs[0], nasa7_coeffs[8:15], nasa7_coeffs[1:8], cp_over_R[i,:], h_over_RT[i,:], s_over_R[i,:])
    if(not consistent_data):
        warning_msg.error_consistency(sp_i, output_dir)
        sys.exit(1)
    else:
        return True


def ct_thermo_continuity(thermo, sp_i):
    """
    Check whether the original Cantera-based thermodynamical coefficients are C0/C1 continuous.
    thermo: a class object compatible with the get_nasa7_coeffs function
    return: Boolean
    """
    nasa7_coeffs = thermo.get_nasa7_coeffs(sp_i)
    continuous_fit = th_fitter.continuity(nasa7_coeffs[0], nasa7_coeffs[8:15], nasa7_coeffs[1:8])
    if(not continuous_fit):
        return False
    else:
        return True


def plot_nasa7_comparison(name, T, cp_over_R, h_over_RT, s_over_R, Tmid, coeffs_lo, coeffs_hi, output_dir):
    
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

    file = os.path.join(output_dir, (name+'_transport.png'))
    fig.savefig(file, bbox_inches='tight')
    fig.clf()


def refit_ct_thermo(thermo, Tmid, output_dir):
    """
    Refit Cantera based thermodynamics (=NASA7-polynomials). If the original NASA polynomials
    are consistent, continuous and share the desired Tmid value, no re-fit is required.
    If data is continuous but does not have a desired Tmid, a new fit is made, but only for Cp. (See thermo_fitter.py)
    If data is not continuous and not have a desired Tmid, a new fit is made for a full cp, h, s matrix system.
    thermo: a class object including access to the Cantera data.
    Tmid: the desired common middle temperature value for NASA7-polynomials
    output_dir: directory into which warning logs are saved.
    return: NASA coefficients (low and high T range).
    """
    nasa_coeffs_lo = np.zeros((thermo.gas.n_species, 7))
    nasa_coeffs_hi = np.zeros((thermo.gas.n_species, 7))
    Tci = np.where(thermo.T==Tmid)[0][0]

    for sp_i in thermo.gas.species_names:
        
        i = thermo.gas.species_index(sp_i)
        cp_over_R, h_over_RT, s_over_R = nasa_normalisation(thermo.T, thermo.cp, thermo.h, thermo.s)
        if(thermo.is_nasa7(sp_i)):

            ct_thermo_consistency(thermo, sp_i, output_dir)
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
            warning_msg.warning_not_nasa7(thermo, output_dir)
            coeffs_lo, coeffs_hi = th_fitter.fit_nasapolys_full(thermo.T, Tci, cp_over_R[i,:], h_over_RT[i,:], s_over_R[i,:], thermo.cp0_over_R[i], thermo.dhf_over_R[i], thermo.s0_over_R[i])

        nasa_coeffs_lo[i] = coeffs_lo
        nasa_coeffs_hi[i] = coeffs_hi

    return nasa_coeffs_lo, nasa_coeffs_hi


def fit_mixture_thermo(data):
    """
    Refit Cantera based thermodynamics (=NASA7-polynomials).
    """
    Tci = np.where(data.T==data.Tmid)[0][0]
    cp_over_R, h_over_RT, s_over_R = nasa_normalisation(data.T, data.cp, data.h, data.s)
    coeffs_lo, coeffs_hi = th_fitter.fit_nasapolys_full(data.T, Tci, cp_over_R[0], h_over_RT[0], s_over_R[0], data.cp0_over_R, data.dhf_over_R, data.s0_over_R)
    return coeffs_lo, coeffs_hi


def nasa7_fit_quality(data, thermo_fits, output_dir, plot=True):
    """
    Ensure that the newly generated NASA7-polynomial coefficients fullfil their 
    continuity-property and that the error to original data is not large.
    data: a class object with an access to original Cantera data.
    thermo_fits: output of the refit_ct_thermo() function
    return: Boolean success
    """
    nasa_coeffs_lo = np.atleast_2d(thermo_fits[0])
    nasa_coeffs_hi = np.atleast_2d(thermo_fits[1])    
    rel_tol = 5e-3
    N = data.cp.shape[0]
    success = True
    for i in range(N):
        err_cp, err_h, err_s = th_fitter.error_nasa7(data.T, data.cp[i,:], data.h[i,:], data.s[i,:], data.Tmid, nasa_coeffs_lo[i,:], nasa_coeffs_hi[i,:], ct.gas_constant)
        
        fit_ok = (err_cp < rel_tol) and (err_h < rel_tol) and (err_s < rel_tol)
        if(not fit_ok):
            success = False
            warning_msg.warning_nasa7_err(data.names[i], rel_tol, err_cp, err_h, err_s, output_dir)

        continuous_fit = th_fitter.continuity(data.Tmid, nasa_coeffs_lo[i,:], nasa_coeffs_hi[i,:])
        if(not continuous_fit):
            success = False
            warning_msg.warning_nasa7_continuity(data.names[i], output_dir)

        if( (not (fit_ok and continuous_fit)) and plot):
            fig_dir = os.path.join(output_dir, "Figures")
            os.makedirs(fig_dir, exist_ok=True)
            cp_over_R, h_over_RT, s_over_R = nasa_normalisation(data.T, data.cp, data.h, data.s)
            plot_nasa7_comparison(data.names[i], data.T, cp_over_R[i,:], h_over_RT[i,:], s_over_R[i,:], data.Tmid, nasa_coeffs_lo[i,:], nasa_coeffs_hi[i,:], fig_dir)

    return success

def get_elements(thermo, sp_i):
    """
    The elemental composition entry for species
    return: dictionary elements = {"C":2, "H":2}
    """
    elements = {}
    for elem_i in thermo.gas.element_names:
        na = thermo.gas.n_atoms(sp_i, elem_i)
        if(na > 0):
            elements[elem_i] = na
    return elements


def ct2foam_thermo_writer(species_file, thermo_file, reactions_file, data, transport_fits, thermo_fits):
    """
    Wrapper to call the OpenFOAM writer to generate entries for thermophysicalProperties file.
    species_file: path to species.foam including a list of species in mechanism based order
    thermo_file: thermo.foam including thermodynamics and transport data entries.
    reactions_file: reactions.foam including an empty dictionary entry for reactions.
    transport_fits: output of fit_ct_transport()
    thermo_fits: output of refit_ct_thermo()
    """
    if os.path.exists(thermo_file):
        os.remove(thermo_file)
    
    if os.path.exists(reactions_file):
        os.remove(reactions_file)
    
    if os.path.exists(species_file):
        os.remove(species_file)

    writer.write_reactions(reactions_file)

    T, Tmid, W = data.T, data.Tmid, data.W
    As, Ts, _, poly_mu, poly_kappa, logpoly_mu, logpoly_kappa  = transport_fits
    nasa7_lo, nasa7_hi = thermo_fits

    writer.write_species_list(species_file, data.gas.species_names)

    for sp_i in data.gas.species_names:
        
        i = data.gas.species_index(sp_i)
        elements = get_elements(data, sp_i)
        writer.write_thermo_transport(thermo_file, sp_i, W[i], As[i], Ts[i], poly_mu[i,:], poly_kappa[i,:], logpoly_mu[i,:], logpoly_kappa[i,:], Tmid, T[0], T[-1], nasa7_lo[i,:], nasa7_hi[i,:], elements=elements)


def ctmix2foam_thermo_writer(species_file, thermo_file, reactions_file, data, transport_fits, thermo_fits):
    """
    Wrapper to call the OpenFOAM writer to generate entries for thermophysicalProperties file.
    species_file: path to species.foam including a list of species in mechanism based order
    thermo_file: thermo.foam including thermodynamics and transport data entries.
    reactions_file: reactions.foam including an empty dictionary entry for reactions.
    transport_fits: output of fit_ct_transport()
    thermo_fits: output of refit_ct_thermo()
    """

    if os.path.exists(thermo_file):
        os.remove(thermo_file)
    
    if os.path.exists(reactions_file):
        os.remove(reactions_file)
    
    if os.path.exists(species_file):
        os.remove(species_file)

    writer.write_reactions(reactions_file)

    T, Tmid, W = data.T, data.Tmid, data.W
    As, Ts, _, poly_mu, poly_kappa, logpoly_mu, logpoly_kappa  = transport_fits
    nasa7_lo, nasa7_hi = thermo_fits
      
    elements = None
    writer.write_thermo_transport(thermo_file, data.names[0], W[0], As[0], Ts[0], poly_mu[0], poly_kappa[0], logpoly_mu[0], logpoly_kappa[0], Tmid, T[0], T[-1], nasa7_lo, nasa7_hi, elements=elements)
