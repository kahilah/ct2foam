import numpy as np
from scipy.optimize import curve_fit
import cantera as ct


def sutherland(T, As, Ts):
    '''
    - Original formulation: mu = mu0*(T0+C)/(T+C)*(T/T0)^(3/2)
    - Here a simplified version is used (as in OpenFOAM): mu = As*sqrt(T)/(1.0 + Ts/T);
    '''
    return  As*np.sqrt(T)/(1.0 + Ts/T)


def euken(mu, cv_mole, W, R):
    """
    The Euken formulation for conductivity based on Sutherland coefficients used in Openfoam
    Foam::scalar Foam::sutherlandTransport<Thermo>::kappa() = mu(p, T)*Cv_*(1.32 + 1.77*this->R()/Cv_);
    """
    Cv = cv_mole/W
    Rspecific = R/W
    return mu*Cv*(1.32 + 1.77*Rspecific/Cv)


def eval_polynomial(poly_coeffs_mu, poly_coeffs_kappa, T):
    '''
    Evaluate polynomial transport fit.
    return: viscosity and conductivity values at given temperature.
    '''
    f_mu = np.poly1d(poly_coeffs_mu)
    f_kappa = np.poly1d(poly_coeffs_kappa)

    return f_mu(T), f_kappa(T)

def eval_log_polynomial(poly_coeffs_mu, poly_coeffs_kappa, T):
    '''
    Evaluate polynomial transport fit.
    return: viscosity and conductivity values at given temperature.
    '''
    f_mu = np.poly1d(poly_coeffs_mu)
    f_kappa = np.poly1d(poly_coeffs_kappa)

    T_ = np.log(T)
    return np.exp(f_mu(T_)), np.exp(f_kappa(T_))


def fit_polynomial(T, mu, kappa, poly_order=3):
    '''
    Make a least-squares polynomial transport fit.
    return: fit coefficients according to numpy standards.
    mu is [N,M]
    kappa is [N,M]
    '''
    # transpose to achieve vectorised fittinh
    mu = np.transpose(mu)
    kappa = np.transpose(kappa)

    poly_coeffs_mu = np.polyfit(T, mu, poly_order)
    poly_coeffs_kappa = np.polyfit(T, kappa, poly_order)

    return poly_coeffs_mu, poly_coeffs_kappa


def fit_log_polynomial(T, mu, kappa, poly_order=3):
    '''
    Make a least-squares polynomial transport fit.
    return: fit coefficients according to numpy standards.
    mu is [N,M]
    kappa is [N,M]
    '''
    # transpose to achieve vectorised fittinh
    mu = np.transpose(mu)
    kappa = np.transpose(kappa)

    mu_ = np.log(mu)
    kappa_ = np.log(kappa)
    T_ = np.log(T)

    poly_coeffs_mu = np.polyfit(T_, mu_, poly_order)
    poly_coeffs_kappa = np.polyfit(T_, kappa_, poly_order)

    return poly_coeffs_mu, poly_coeffs_kappa



def fit_sutherland(T, mu, p0=np.array([1.0, 1.0])):
    """
    Make a least-squares polynomial transport fit.
    return: fit coefficients according to curve_fit standards.
    mu is [N,M]
    kappa is [N,M]
    """
    popt, pcov = curve_fit(sutherland, T, mu, p0=p0)
    As = popt[0]
    Ts = popt[1]
    # one standard deviation errors
    std_err = np.sqrt(np.diag(pcov))

    return As, Ts, std_err
    

def error_sutherland(mu, kappa, T, As, Ts, cv_mole, W, R=ct.gas_constant):
    '''
    Assuming [M,] shaped arrays (column vector data)
    calculate species-wise L2 error for the sutherland fit
    T is array of sample rate
    As and Ts are sutherland coefficients with size N corresponding number of species
    '''

    mu_sutherland = sutherland(T, As, Ts)
    kappa_euken = euken(mu_sutherland, cv_mole, W, R)

    diff_mu = np.abs(mu - mu_sutherland)
    diff_kappa = np.abs(kappa - kappa_euken)
    err_mu = np.linalg.norm(diff_mu)/np.linalg.norm(mu)
    err_kappa = np.linalg.norm(diff_kappa)/np.linalg.norm(kappa)

    return err_mu, err_kappa


def error_polynomial(mu, kappa, poly_coeffs_mu, poly_coeffs_kappa, T):
    '''
    calculate species-wise L2 error for the sutherland fit
    T is array of sample rate
    As and Ts are sutherland coefficients with size N corresponding number of species
    '''
    mu_poly, kappa_poly = eval_polynomial(poly_coeffs_mu, poly_coeffs_kappa, T)
    
    diff_mu = np.abs(mu - mu_poly)
    diff_kappa = np.abs(kappa - kappa_poly)
    err_mu = np.linalg.norm(diff_mu)/np.linalg.norm(mu)
    err_kappa = np.linalg.norm(diff_kappa)/np.linalg.norm(kappa)

    return err_mu, err_kappa


def error_log_polynomial(mu, kappa, poly_coeffs_mu, poly_coeffs_kappa, T):
    '''
    calculate species-wise L2 error for the sutherland fit
    T is array of sample rate
    poly* are coefficients with size N corresponding number of species
    '''
    mu_poly, kappa_poly = eval_log_polynomial(poly_coeffs_mu, poly_coeffs_kappa, T)
    
    diff_mu = np.abs(mu - mu_poly)
    diff_kappa = np.abs(kappa - kappa_poly)
    err_mu = np.linalg.norm(diff_mu)/np.linalg.norm(mu)
    err_kappa = np.linalg.norm(diff_kappa)/np.linalg.norm(kappa)

    return err_mu, err_kappa    