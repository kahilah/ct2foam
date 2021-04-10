import numpy as np
from . import lsqlin #3rd party matlab like lsqlin, providing better performance than scipy (in 2016)

"""
See the following links for further information on the NASA polynomial fitting procedure:
- https://pdfs.semanticscholar.org/4920/6eb96b41fcbc8b19526b2cf2a3b10e02e0b1.pdf
- http://shepherd.caltech.edu/EDL/publications/reprints/RefittingThermoDataNew.pdf
- https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19930003779.pdf
- https://ntrs.nasa.gov/api/citations/19940013151/downloads/19940013151.pdf

This function fits polynomials to specific heat, using NASA 7 term format from SP-272. Gordon and McBride 1971.
- Eq. 90: cp0/R=a1+a2T+a3*T^2+a4*T^3+a5*T^4 
- Eq. 91: H0_T/(RT)=a1+a2/2*T+a3/3*T^2+a4/4*T^3+a5/5*T^4+a6/T 
- Eq. 92: S0_T/R=a1*ln(T)+a2*T+a3/2*T^2+a4/3*T^3+a5/4*T^4+a7

Depending on the initial data, the fitting is made by either of the following methods:

1) Original data is from high-quality fit or from partition functions.
- Fit the specific heat only and use the analytical integration to obtain the additional coefficients for entropy and enthalpy.  
The specific heat fit constrains the derivative to be continuous at Tcommon and the value of 
specific heat at 298.15 K is constrained to be the specified value of the standard state (enthalpy of formation). 
The entropy and enthalpy values at the midpoints are forced to be continuous by analytical constaining.

- This procedure produces high quality fits as long as the input entropy and enthalpy are thermodynamically 
consistent with the specific heat data.  This will always be the case if the thermodynamic data have been 
derived from partition functions and presumably when refitting existing fit functions. Be cautious in refitting!

2) Original data from low-quality fit or from experiments
- Create a simultaneous least squares fit considering the whole system with prescribed constraints.


Fitting procedure: 

Solve linear constrained l2-regularized least squares by a Matlab's lsqlin equivalent algorithm. It is actually wrapper around CVXOPT QP solver.
    min_x ||C*x  - d||^2_2 + reg * ||x||^2_2
    s.t.  A * x <= b
            Aeq * x = beq
            lb <= x <= ub

NOTES:
- Note that multiplying with a fraction prior to the exponent ((1./2.)**(1./3.)*T)**3 ensures higher numerical arithmetic accuracy
"""

# global variables in this module
__T_std = 298.15 # private

def _cp_nasa7(c, T):
    """
    Private function for direct call of nasa 7-coefficient polynomial.
    No Tcommon checks.
    c is a coefficient array of size 4.
    return cp(T)/R
    """
    return c[0] + c[1]*T + c[2]*pow(T,2) + c[3]*pow(T,3) + c[4]*pow(T,4)

def cp_nasa7(T, Tmid, c_lo, c_hi):
    """
    Full treatment of all temperature ranges in a numpy efficient manner.
    coeffs[0] is Tcommon
    T is [M,] shaped array
    return cp(T)/R
    """
    T = np.atleast_1d(T)
    cp_over_R = np.zeros(np.shape(T))
    i_lo = (T <= Tmid)
    i_hi = (T > Tmid)
    cp_over_R[i_lo] = _cp_nasa7(c_lo, T[i_lo])
    cp_over_R[i_hi] = _cp_nasa7(c_hi, T[i_hi])   
    return cp_over_R


def _h_nasa7(c, T):
    """
    Private function for direct call of nasa 7-coefficient polynomial.
    No Tcommon checks.
    c is a coefficient array of size 4.
    return h(T)/TR
    """
    return c[0] + c[1]*T/2 + c[2]*pow(T,2)/3 + c[3]*pow(T,3)/4 + c[4]*pow(T,4)/5 + c[5]/T    

def h_nasa7(T, Tmid, c_lo, c_hi):
    """
    h/RT
    return h(T)/TR
    """
    T = np.atleast_1d(T)
    h_over_RT = np.zeros(len(T))
    i_lo = (T <= Tmid)
    i_hi = (T > Tmid)
    h_over_RT[i_lo] = _h_nasa7(c_lo, T[i_lo])
    h_over_RT[i_hi] = _h_nasa7(c_hi, T[i_hi])   
    return h_over_RT


def _s_nasa7(c, T):
    """
    Private function for direct call of nasa 7-coefficient polynomial.
    No Tcommon checks.
    c is a coefficient array of size 4.
    return s(T)/R
    """
    return c[0]*np.log(T) + c[1]*T + c[2]*pow(T,2)/2 + c[3]*pow(T,3)/3 + c[4]*pow(T,4)/4 + c[6]

def s_nasa7(T, Tmid, c_lo, c_hi):
    """
    return s(T)/R
    """
    T = np.atleast_1d(T)
    s_over_R = np.zeros(len(T))
    i_lo = (T <= Tmid)
    i_hi = (T > Tmid)
    s_over_R[i_lo] = _s_nasa7(c_lo, T[i_lo])
    s_over_R[i_hi] = _s_nasa7(c_hi, T[i_hi])   
    return s_over_R


def _dcpdT_nasa7(c, T):
    """
    Private function for direct call of nasa 7-coefficient polynomial.
    No Tcommon checks.
    c is a coefficient array of size 4.
    return d cp(T)/dT * 1/R
    """
    return c[1] + 2.*c[2]*T + 3.*c[3]*pow(T,2) + 4.*c[4]*pow(T,3)

def dcpdT_nasa7(T, Tmid, c_lo, c_hi):
    """
    Derivative of cp/R w.r.t. temperature == 1/R * dcp/dT
    """
    T = np.atleast_1d(T)
    dcpdT_over_R = np.zeros(len(T))
    i_lo = (T <= Tmid)
    i_hi = (T > Tmid)
    dcpdT_over_R[i_lo] = _dcpdT_nasa7(c_lo, T[i_lo])
    dcpdT_over_R[i_hi] = _dcpdT_nasa7(c_hi, T[i_hi])   
    return dcpdT_over_R


def consistency(T, Tmid, c_lo, c_hi, cp_over_R, h_over_RT, s_over_R, abs_tol=1e-6):
    """
    cp: [M,] sized array where M is the temperature sampling ratio.
    coeffs: coeffs[0] is Tcommon, rest are NASA 7 coefficients (2*7) in chemkin order (first high then low)
    """
    dcp = np.abs(cp_over_R - cp_nasa7(T, Tmid, c_lo, c_hi))
    err_cp = np.linalg.norm(dcp)
    dh = np.abs(h_over_RT - h_nasa7(T, Tmid, c_lo, c_hi))
    err_h = np.linalg.norm(dh)
    ds = np.abs(s_over_R - s_nasa7(T, Tmid, c_lo, c_hi))
    err_s = np.linalg.norm(ds)

    if((err_cp < abs_tol) and (err_h < abs_tol) and (err_s < abs_tol)):
        return True
    else:
        return False



def continuity(Tmid, c_lo, c_hi, cp_tol=1e-6, cpdT_tol=0.01, h_tol=1e-6, s_tol=1e-6):
    """
    check_quality checks whether the nasa polynomial coefficients fulfill
    the original requirements on C0 and C1 continuity.
    - C0 continuity: 1e-6 error is acceptable
    - C1 continuity (not mandatory but recommended): For some mechanism fits C1 continuity is not guaranteed --> a loose tolerance is chosen
    :coeffs: 
    :return: boolean whether a good fit or not.
    """    
    # C0 continuity of cp
    cp_mid_lo = _cp_nasa7(c_lo, Tmid)
    cp_mid_hi = _cp_nasa7(c_hi, Tmid)
    abs_err = np.abs(cp_mid_lo - cp_mid_hi)
    C0_cp = (abs_err < cp_tol)
    # C1 continuity of cp
    dcpdT_mid_lo = _dcpdT_nasa7(c_lo, Tmid)
    dcpdT_mid_hi = _dcpdT_nasa7(c_hi, Tmid)
    abs_err = np.abs(dcpdT_mid_lo - dcpdT_mid_hi)
    C1_cp = (abs_err < cpdT_tol)
    # C0 continuity of h
    h_mid_lo = _h_nasa7(c_lo, Tmid)
    h_mid_hi = _h_nasa7(c_hi, Tmid)
    abs_err = np.abs(h_mid_lo - h_mid_hi)
    C0_h = (abs_err < h_tol)
    # C0 continuity of s
    s_mid_lo = _s_nasa7(c_lo, Tmid)
    s_mid_hi = _s_nasa7(c_hi, Tmid)
    abs_err = np.abs(s_mid_lo - s_mid_hi)
    C0_s =  (abs_err < s_tol)

    continuous = (C0_cp and C0_h and C0_s and C1_cp)

    return continuous 


def correct_coeffs(coeffs, Tcommon, dhf_over_R, s0_over_R):
    """
    Solving additional constant over means of conservation of enthalpy and entropy + ensuring C0 continuity at T=Tcommon 
    coeffs = 14 size nasa coeffs without Tcommon in 0
    """
    # Define coeff[5] in terms of enthalpy of formation at standard conditions
    coeffs[5] = dhf_over_R - ( coeffs[0]*__T_std + coeffs[1]*((1./2.)**(1./2.)*__T_std)**2 \
            + coeffs[2]*((1./3.)**(1./3.)*__T_std)**3 + coeffs[3]*((1./4.)**(1./4.)*__T_std)**4 + coeffs[4]*((1./5.)**(1./5.)*__T_std)**5 )
    # Define coeffs[5+7] by ensuring continuity at Tcommon
    h_sum_term_l = lambda T:  coeffs[0] + coeffs[1]*(1./2.)*T + coeffs[2]*((1./3.)**(1./2.)*T)**2 + coeffs[3]*((1./4.)**(1./3.)*T)**3 + coeffs[4]*((1./5.)**(1./4.)*T)**4 
    h_sum_term_h = lambda T:  coeffs[0+7] + coeffs[1+7]*(1./2.)*T + coeffs[2+7]*((1./3.)**(1./2.)*T)**2 + coeffs[3+7]*((1./4.)**(1./3.)*T)**3 + coeffs[4+7]*((1./5.)**(1./4.)*T)**4 
    coeffs[5+7] = Tcommon*( h_sum_term_l(Tcommon) + coeffs[5]/Tcommon - h_sum_term_h(Tcommon) )
    
    # Define coeff[6] in terms of entropy at standard conditions
    coeffs[6] = s0_over_R - ( coeffs[0]*np.log(__T_std) + coeffs[1]*__T_std + coeffs[2]*((1./2.)**(1./2.)*__T_std)**2 \
            + coeffs[3]*((1./3.)**(1./3.)*__T_std)**3 + coeffs[4]*((1./4.)**(1./4.)*__T_std)**4  )
    # Define coeffs[6+7] by ensuring continuity at Tcommon
    s_sum_term_l = lambda T:  coeffs[0]*np.log(T) + coeffs[1]*T + coeffs[2]*((1./2.)**(1./2.)*T)**2 + coeffs[3]*((1./3.)**(1./3.)*T)**3 + coeffs[4]*((1./4.)**(1./4.)*T)**4 
    s_sum_term_h = lambda T:  coeffs[0+7]*np.log(T) + coeffs[1+7]*T + coeffs[2+7]*((1./2.)**(1./2.)*T)**2 + coeffs[3+7]*((1./3.)**(1./3.)*T)**3 + coeffs[4+7]*((1./4.)**(1./4.)*T)**4 
    coeffs[6+7] = s_sum_term_l(Tcommon) + coeffs[6] - s_sum_term_h(Tcommon) 

    return coeffs


def fit_nasapolys_cp(T0, Tc_i, cp_over_R, cp0_over_R, dhf_over_R, s0_over_R, verbose=False):
    """
    take notes from main notes
    Supports 1d-arrays only [M,] sized, where M refers to number of T0 array (sample) size.
    Tc_i is the common temperature index in the T0 array
    """
    if(verbose):
        print('\tFitting cp/R only and derive the remaining coeffcients analytically...')

    Tcommon = T0[Tc_i]
    T_low = T0[0:Tc_i+1]  
    T_high = T0[Tc_i:]
    T = np.concatenate((T_low,T_high))

    Nl=len(T_low)
    Nh=len(T_high)
    N=len(T)
    
    M = 5
    C = np.zeros((N,2*M))
    d = np.zeros(N)

    # cp_over_R formulation
    for i in range(0,M):
        C[:Nl,i] = pow(T_low,i)
        C[Nl:Nl+Nh,i+M] = pow(T_high,i)

    #The right hand side
    d[:Nl]      = cp_over_R[0:Tc_i+1] # low part
    d[Nl:Nl+Nh] = cp_over_R[Tc_i:]    # high part

    # Constraint matrix and rhs
    Aeq = np.zeros((5,2*M))
    beq = np.zeros(5)
    
    # cp_over_R : C0 continuity at Tcommon and equal first and last element values
    for i in range(0,M):
        Aeq[0,i]    = Tcommon**i
        Aeq[0,i+M]  = -Tcommon**i
        Aeq[1,i]    = Tcommon**i
        Aeq[2,i]    = __T_std**i
        Aeq[3,i+M]  = T[-1]**i

    beq[0] = 0.0
    beq[1] = cp_over_R[Tc_i]
    beq[2] = cp0_over_R
    beq[3] = cp_over_R[-1]
    
    # cp_over_R : C1 continuity at Tcommon
    Aeq[4,1]    = 1.
    Aeq[4,2]    = 2.*Tcommon
    Aeq[4,3]    = ((3.**(1./2.))*Tcommon)**2
    Aeq[4,4]    = ((4.**(1./3.))*Tcommon)**3
    Aeq[4,1+M]  = -1.
    Aeq[4,2+M]  = -2.*Tcommon
    Aeq[4,3+M]  = -((3.**(1./2.))*Tcommon)**2
    Aeq[4,4+M]  = -((4.**(1./3.))*Tcommon)**3
    
    # Find the least-squares solution
    sol = lsqlin.lsqlin(C, d, 0, None, None, Aeq, beq,-1e9,1e9,None,{'show_progress': False, 'abstol': 1e-12, 'reltol': 1e-8})
    coeffs_tmp = sol['x']

    # fill the coefficient matrix
    coeffs = np.zeros(14) 
    for i in range(0,5):
        coeffs[i] = coeffs_tmp[i]
        coeffs[i+7] = coeffs_tmp[i+M]


    # solving additional constant by means of conservation of enthalpy and entropy + ensuring C0 continuity at T=Tcommon 
    coeffs_corrected = correct_coeffs(coeffs, Tcommon, dhf_over_R, s0_over_R)
    coeffs_lo = coeffs_corrected[:7]
    coeffs_hi = coeffs_corrected[7:]

    return coeffs_lo, coeffs_hi


def fit_nasapolys_full(T0, Tc_i, cp_over_R, h_over_RT, s_over_R, cp0_over_R, dhf_over_R, s0_over_R, verbose=False):
    """
    take notes from main notes
    Supports 1d-arrays only [M,] sized, where M refers to number of T0 array (sample) size.
    Tc_i is the common temperature index in the T0 array
    """
    if(verbose):
        print( '\tFitting cp/R, h/RT and s/R simultaneously...' )

    Tcommon = T0[Tc_i]
    T_low = T0[0:Tc_i+1]  
    T_high = T0[Tc_i:]
    T = np.concatenate((T_low,T_high))

    cp_over_R_L = cp_over_R[0:Tc_i+1]
    cp_over_R_H = cp_over_R[Tc_i:]

    h_over_RT_L = h_over_RT[0:Tc_i+1]
    h_over_RT_H = h_over_RT[Tc_i:]

    s_over_R_L = s_over_R[0:Tc_i+1]
    s_over_R_H = s_over_R[Tc_i:]

    Nl=len(T_low)
    Nh=len(T_high)
    N=len(T)

    M = 5
    C = np.zeros((3*N,2*M))
    d = np.zeros(3*N)

    # cp formulation
    for i in range(0,5):
        C[:Nl,i] = pow(T_low,i)
        C[Nl:Nl+Nh,i+M] = pow(T_high,i)           
    
    # h equation
    c1 = 1./2.
    c2 = (1./3.)**(1./2.)
    c3 = (1./4.)**(1./3.)
    c4 = (1./5.)**(1./4.)
    
    i1=Nl+Nh
    i2=i1+Nl
    C[i1:i2,0] = np.ones((1,Nl)) - __T_std/T_low
    C[i1:i2,1] = (c1*T_low)    - (__T_std/T_low)*(c1*__T_std)   
    C[i1:i2,2] = (c2*T_low)**2 - (__T_std/T_low)*(c2*__T_std)**2
    C[i1:i2,3] = (c3*T_low)**3 - (__T_std/T_low)*(c3*__T_std)**3
    C[i1:i2,4] = (c4*T_low)**4 - (__T_std/T_low)*(c4*__T_std)**4

    i3=i2
    i4=i2+Nh
    C[i3:i4,0+M] = np.ones((1,Nh)) - __T_std/T_high
    C[i3:i4,1+M] = (c1*T_high)    - (__T_std/T_high)*(c1*__T_std)
    C[i3:i4,2+M] = (c2*T_high)**2 - (__T_std/T_high)*(c2*__T_std)**2
    C[i3:i4,3+M] = (c3*T_high)**3 - (__T_std/T_high)*(c3*__T_std)**3
    C[i3:i4,4+M] = (c4*T_high)**4 - (__T_std/T_high)*(c4*__T_std)**4

    
    # s equation
    c2 = (1./2.)**(1./2.)
    c3 = (1./3.)**(1./3.)
    c4 = (1./4.)**(1./4.)

    i5=i4
    i6=i4+Nl
    C[i5:i6,0] = np.log(T_low/__T_std)
    C[i5:i6,1] = T_low - __T_std
    C[i5:i6,2] = (c2*T_low)**2 - (c2*__T_std)**2
    C[i5:i6,3] = (c3*T_low)**3 - (c3*__T_std)**3
    C[i5:i6,4] = (c4*T_low)**4 - (c4*__T_std)**4
    
    i7=i6
    i8=i6+Nh
    C[i7:i8,0+M] = np.log(T_high/__T_std)
    C[i7:i8,1+M] = T_high - __T_std
    C[i7:i8,2+M] = (c2*T_high)**2 - (c2*__T_std)**2
    C[i7:i8,3+M] = (c3*T_high)**3 - (c3*__T_std)**3
    C[i7:i8,4+M] = (c4*T_high)**4 - (c4*__T_std)**4
    
    #The right hand side
    d[:Nl]=cp_over_R_L
    d[Nl:Nl+Nh]=cp_over_R_H
    
    d[i1:i2] = h_over_RT_L - dhf_over_R/T_low
    d[i3:i4] = h_over_RT_H - dhf_over_R/T_high
    
    d[i5:i6] = s_over_R_L - s0_over_R
    d[i7:i8] = s_over_R_H - s0_over_R

    # Constraints
    # Note that you can add here constraints if e.g. maximum values are known etc.
    Aeq=np.zeros((3,2*M))
    beq=np.zeros(3)
    
    # cp_over_R : C0 continuity at Tcommon (presumably no knowledge on min/max/std values)
    for i in range(0,5):
        Aeq[0,i] = Tcommon**i
        Aeq[0,i+M] = -Tcommon**i
        Aeq[2,i]    = __T_std**i    # not mandatory

    # cp_over_R : C1 continuity at Tcommon
    Aeq[1,1] = 1.
    Aeq[1,2] = 2.*Tcommon
    Aeq[1,3] = ((3.**(1./2.))*Tcommon)**2
    Aeq[1,4] = ((4.**(1./3.))*Tcommon)**3
    Aeq[1,1+M] = -1.
    Aeq[1,2+M] = -2.*Tcommon
    Aeq[1,3+M] = -((3.**(1./2.))*Tcommon)**2
    Aeq[1,4+M] = -((4.**(1./3.))*Tcommon)**3
    
    # constraint rhs
    beq[2] = cp0_over_R     # not mandatory

    '''       
    # - h and C0 continuity is guaranteed after solution by analytical consideration
    # - if forcing the constraints to the linear solution (below), result are typically less good
    #   and the other coefficients do the job anyways as suggested by pen and paper
    # h : C0 continuity 
    Aeq[2,0] = 1. - __T_std/Tcommon
    Aeq[2,M] = -(1. - __T_std/Tcommon)
    for i in range(1,5):
        cf = (1./(i+1))**(1./i)
        Aeq[2,i] = (cf*Tcommon)**i - (__T_std/Tcommon)*(cf*__T_std)**i  
        Aeq[2,i+M] = -((cf*Tcommon)**i - (__T_std/Tcommon)*(cf*__T_std)**i)  
    # s : C0 continuity 
    Aeq[3,0] = np.log(Tcommon/__T_std)
    Aeq[3,M] = -np.log(Tcommon/__T_std)
    for i in range(1,5):
        cf = (1./i)**(1./i)
        Aeq[3,i] = (cf*Tcommon)**i - (cf*__T_std)**i
        Aeq[3,i+M] = -(cf*Tcommon)**i + (cf*__T_std)**i
    '''

    # Find the least-squares solution
    sol=lsqlin.lsqlin(C, d, 0, None, None, Aeq, beq,-1e9,1e9,None,{'show_progress': False, 'abstol': 1e-12, 'reltol': 1e-8})
    coeffs_tmp=sol['x']

    # fill the coefficient matrix
    coeffs=np.zeros(14) 
    for i in range(0,5):
        coeffs[i] = coeffs_tmp[i]
        coeffs[i+7] = coeffs_tmp[i+M]
    
    # solving additional constant by means of conservation of enthalpy and entropy + ensuring C0 continuity at T=Tcommon 
    coeffs_corrected = correct_coeffs(coeffs, Tcommon, dhf_over_R, s0_over_R)
    coeffs_lo = coeffs_corrected[:7]
    coeffs_hi = coeffs_corrected[7:]

    return coeffs_lo, coeffs_hi


def error_nasa7(T, cp, h, s, Tmid, coeffs_lo, coeffs_hi, R=8314.46261815324):

    cp_fit = cp_nasa7(T, Tmid, coeffs_lo, coeffs_hi)*R
    dcp = np.abs(cp - cp_fit)
    err_cp = np.linalg.norm(dcp)/np.linalg.norm(cp)

    h_fit = h_nasa7(T, Tmid, coeffs_lo, coeffs_hi)*R*T
    dh = np.abs(h - h_fit)
    err_h = np.linalg.norm(dh)/np.linalg.norm(h)

    s_fit = s_nasa7(T, Tmid, coeffs_lo, coeffs_hi)*R
    ds = np.abs(s - s_fit)
    err_s = np.linalg.norm(ds)/np.linalg.norm(s)

    return err_cp, err_h, err_s
