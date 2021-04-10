import numpy as np
from scipy.optimize import curve_fit


def sutherland(T, As, Ts):
    '''
    - Original formulation: mu = mu0*(T0+C)/(T+C)*(T/T0)^(3/2)
    - Here a simplified version is used (as in OpenFOAM): mu = As*sqrt(T)/(1.0 + Ts/T);
    '''
    return  As*np.sqrt(T)/(1.0 + Ts/T)


def euken(mu, cv_mole, W, R):
    # The Euken formulation for conductivity based on Sutherland coefficients used in Openfoam
    # Foam::scalar Foam::sutherlandTransport<Thermo>::kappa() = mu(p, T)*Cv_*(1.32 + 1.77*this->R()/Cv_);
    Cv = cv_mole/W
    Rspecific = R/W
    return mu*Cv*(1.32 + 1.77*Rspecific/Cv)


def eval_polynomial(poly_coeffs_mu, poly_coeffs_kappa, T, poly_order=3, logBasis=True):
    '''
    poly_order=3 yields Chemkin type fit for polynomial for viscosity and conductivity
    '''
    f_mu = np.poly1d(poly_coeffs_mu)
    f_kappa = np.poly1d(poly_coeffs_kappa)

    if(logBasis):
        T_ = np.log(T)
        return np.exp(f_mu(T_)), np.exp(f_kappa(T_))
    else:
        T_ = T
        return f_mu(T_), f_kappa(T_)

    


def fit_polynomial(T, mu, kappa, poly_order=3, logBasis=True):
    '''
    poly_order=3 yields Chemkin type fit for polynomial fit for viscosity and conductivity
    mu is [N,M]
    kappa is [N,M]
    '''
    # transpose to achieve vectorised fittinh
    mu = np.transpose(mu)
    kappa = np.transpose(kappa)

    if(logBasis):
        mu_ = np.log(mu)
        kappa_ = np.log(kappa)
        T_ = np.log(T)
    else:
        mu_ = mu
        kappa_ = kappa
        T_ = T

    poly_coeffs_mu = np.polyfit(T_, mu_, poly_order)
    poly_coeffs_kappa = np.polyfit(T_, kappa_, poly_order)

    return poly_coeffs_mu, poly_coeffs_kappa



def fit_sutherland(T, mu, p0=np.array([1.0, 1.0])):
    """
        mu is M or [M,N] array where N is the number of mixture states (species) and M is the temperature sampling size
        p0 #initial guess
    """
    mu = np.atleast_2d(mu) # consistent behavior with 1d and 2d arrays
    N = mu.shape[0]
    As = np.zeros(N)
    Ts = np.zeros(N)
    std_err = np.zeros((N, 2))

    for i in range(N):
        mu_i = mu[i,:]
        popt, pcov = curve_fit(sutherland, T, mu_i, p0=p0)
        As[i] = popt[0]
        Ts[i] = popt[1]

        # one standard deviation errors
        std_err[i,:] = np.sqrt(np.diag(pcov))

    return As, Ts, std_err
    

def error_sutherland(mu, kappa, T, As, Ts, cv_mole, W, R):
    '''
    Assuming [M,] or [M,K] shaped arrays (column vector data)
    calculate species-wise L2 error for the sutherland fit
    T is array of sample rate
    As and Ts are sutherland coefficients with size N corresponding number of species
    '''
    # consistent behavior with 1d and 2d arrays
    mu = np.atleast_2d(mu) 
    kappa = np.atleast_2d(kappa) 
    cv_mole = np.atleast_2d(cv_mole) 
    As = np.atleast_1d(As) 
    Ts = np.atleast_1d(Ts) 
    W = np.atleast_1d(W) 

    N = len(As)
    err_mu = np.zeros(N)
    err_kappa = np.zeros(N)
    for i in range(N):
        mu_sutherland = sutherland(T, As[i], Ts[i])
        kappa_euken = euken(mu_sutherland, cv_mole[i,:], W[i], R)
        mu_err = np.abs(mu[i,:] - mu_sutherland)
        kappa_err = np.abs(kappa[i,:] - kappa_euken)

        err_mu[i] = np.linalg.norm(mu_err)/np.linalg.norm(mu[i,:])
        err_kappa[i] = np.linalg.norm(kappa_err)/np.linalg.norm(kappa[i,:])

    if((np.max(err_mu) > 0.05) or (np.max(err_kappa) > 0.1)):
        print("\nWARNING: Sutherland transport fit quality may be poor:")
        print("\tmax(err(mu))="+repr(np.max(err_mu)))
        print("\tmax(err(kappa))="+repr(np.max(err_kappa))+"\n")

    return err_mu, err_kappa


def error_polynomial(mu, kappa, poly_coeffs_mu, poly_coeffs_kappa, T, poly_order=3, logBasis=True):
    '''
    calculate species-wise L2 error for the sutherland fit
    T is array of sample rate
    As and Ts are sutherland coefficients with size N corresponding number of species
    '''
    # consistent behavior with 1d and 2d arrays
    mu = np.atleast_2d(mu) 
    kappa = np.atleast_2d(kappa) 

    N = mu.shape[0]
    err_mu = np.zeros(N)
    err_kappa = np.zeros(N)

    for i in range(N):
        mu_poly, kappa_poly = eval_polynomial(poly_coeffs_mu[:,i], poly_coeffs_kappa[:,i], T, poly_order, logBasis)
        mu_err = np.abs(mu[i,:] - mu_poly)
        kappa_err = np.abs(kappa[i,:] - kappa_poly)

        err_mu[i] = np.linalg.norm(mu_err)/np.linalg.norm(mu[i,:])
        err_kappa[i] = np.linalg.norm(kappa_err)/np.linalg.norm(kappa[i,:])

    if((np.max(err_mu) > 0.01) or (np.max(err_kappa) > 0.01)):
        print("\nWARNING: Polynomial transport fit quality may be poor:")
        print("\tmax(err(mu))="+repr(np.max(err_mu)))
        print("\tmax(err(kappa))="+repr(np.max(err_kappa))+"\n")

    return err_mu, err_kappa