import numpy as np
from . import lsqlin #3rd party matlab like lsqlin, providing better performance than scipy (in 2016)

class nasaPolynomials:
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


    def __init__(self):
        """
        Tcommon:  
        """
        self.T_std = 298.15 


    def cp_over_R(self, coeffs, T):
        """
        coeffs[0] is Tcommon
        """
        if(T <= coeffs[0]):
            c = coeffs[1:8]
        else: 
            c = coeffs[8:15]

        cp_over_R = c[0] + c[1]*T + c[2]*pow(T,2) + c[3]*pow(T,3) + c[4]*pow(T,4)
        
        return cp_over_R


    def cpdT_over_R(self, coeffs, T):
        """
        Derivative of cp/R w.r.t. temperature == 1/R * dcp/dT
        """
        if(T <= coeffs[0]):
            c = coeffs[1:8]
        else: 
            c = coeffs[8:15]

        cpdT_over_R = c[1] + 2.*c[2]*T + 3.*c[3]*pow(T,2) + 4.*c[4]*pow(T,3)
        
        return cpdT_over_R


    def h_over_RT(self, coeffs, T):
        """
        h/RT
        """
        if(T <= coeffs[0]):
            c = coeffs[1:8]
        else: 
            c = coeffs[8:15]

        h_over_RT = c[0] + c[1]*T/2 + c[2]*pow(T,2)/3 + c[3]*pow(T,3)/4 + c[4]*pow(T,4)/5 + c[5]/T
        
        return h_over_RT


    def s_over_R(self, coeffs, T):
        """
        s/R
        """
        if(T <= coeffs[0]):
            c = coeffs[1:8]
        else: 
            c = coeffs[8:15]

        s_over_R = c[0]*np.log(T) + c[1]*T + c[2]*pow(T,2)/2 + c[3]*pow(T,3)/3 + c[4]*pow(T,4)/4 + c[6]
        
        return s_over_R


    def check_quality(self, coeffs, cp_tol=1e-6, cpdT_tol=0.01, h_tol=1e-6, s_tol=1e-6):
        """
        check_quality checks whether the nasa polynomial coefficients fulfill
        the original requirements on C0 and C1 continuity.
        - C0 continuity: 1e-6 error is acceptable
        - C1 continuity (not mandatory but recommended): For some mechanism fits C1 continuity is not guaranteed --> a loose tolerance is chosen
        :coeffs: array with common temperature and low and high temperature coefficients
        :return: boolean whether a good fit or not.
        """ 
        
        # For the following continuity test, fill low/high coefficients to ensure evaluation with both coefficients independent of temperature
        coeffs_low = np.copy(coeffs)
        coeffs_high = np.copy(coeffs)
        coeffs_low[8:15] = coeffs[1:8]
        coeffs_high[1:8] = coeffs[8:15]

        Tcommon = coeffs[0]

        # C0 continuity of cp
        cp_at_Tc_low = self.cp_over_R(coeffs_low, Tcommon)
        cp_at_Tc_high = self.cp_over_R(coeffs_high, Tcommon)
        C0_cp = np.abs(cp_at_Tc_low - cp_at_Tc_high)/(np.abs(cp_at_Tc_low)+1e-15) < cp_tol
        # C1 continuity of cp
        cpdT_at_Tc_low = self.cpdT_over_R(coeffs_low, Tcommon)
        cpdT_at_Tc_high = self.cpdT_over_R(coeffs_high, Tcommon)
        C1_cp = np.abs(cpdT_at_Tc_low - cpdT_at_Tc_high)/(np.abs(cpdT_at_Tc_low)+1e-15) < cpdT_tol
        # C0 continuity of h
        h_at_Tc_low = self.h_over_RT(coeffs_low, Tcommon)
        h_at_Tc_high = self.h_over_RT(coeffs_high, Tcommon)
        C0_h = np.abs(h_at_Tc_low - h_at_Tc_high)/(np.abs(h_at_Tc_low)+1e-15) < h_tol
        # C0 continuity of s
        s_at_Tc_low = self.s_over_R(coeffs_low, Tcommon)
        s_at_Tc_high = self.s_over_R(coeffs_high, Tcommon)
        C0_s = np.abs(s_at_Tc_low - s_at_Tc_high)/(np.abs(s_at_Tc_low)+1e-15) < s_tol

        good_fit = (C0_cp and C0_h and C0_s and C1_cp)

        return good_fit 



    def fit_nasapolys_0(self, T0, Tc_i, cp_over_R, h_over_RT, s_over_R, cp0_by_R, dhf_by_R, s0_by_R):
        """
        take notes from main notes
        Supports 1d-arrays only [M,] sized, where M refers to number of T0 array (sample) size.
        Tc_i is the common temperature index in the T0 array
        """
        Tcommon = T0[Tc_i]
        T_low = T0[0:Tc_i+1]  
        T_high = T0[Tc_i:]
        T = np.concatenate((T_low,T_high))

        Nl=len(T_low)
        Nh=len(T_high)
        N=len(T)
        

        print( '\nFitting cp/R only and derive the remaining coeffcients analytically:\n' )

        M = 5
        C = np.zeros((N,2*M))
        d = np.zeros(N)

        # cp equation
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
            Aeq[2,i]    = self.T_std**i
            Aeq[3,i+M]  = T[-1]**i

        beq[0] = 0.0
        beq[1] = cp_over_R[Tc_i]
        beq[2] = cp0_by_R
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

        # Define coeff[5] in terms of enthalpy of formation at standard conditions
        coeffs[5] = dhf_by_R - ( coeffs[0]*self.T_std + coeffs[1]*((1./2.)**(1./2.)*self.T_std)**2 \
                + coeffs[2]*((1./3.)**(1./3.)*self.T_std)**3 + coeffs[3]*((1./4.)**(1./4.)*self.T_std)**4 + coeffs[4]*((1./5.)**(1./5.)*self.T_std)**5 )
        # Define coeffs[5+7] by ensuring continuity at Tcommon
        h_sum_term_l = lambda T:  coeffs[0] + coeffs[1]*(1./2.)*T + coeffs[2]*((1./3.)**(1./2.)*T)**2 + coeffs[3]*((1./4.)**(1./3.)*T)**3 + coeffs[4]*((1./5.)**(1./4.)*T)**4 
        h_sum_term_h = lambda T:  coeffs[0+7] + coeffs[1+7]*(1./2.)*T + coeffs[2+7]*((1./3.)**(1./2.)*T)**2 + coeffs[3+7]*((1./4.)**(1./3.)*T)**3 + coeffs[4+7]*((1./5.)**(1./4.)*T)**4 
        coeffs[5+7] = Tcommon*( h_sum_term_l(Tcommon) + coeffs[5]/Tcommon - h_sum_term_h(Tcommon) )
        
        # Define coeff[6] in terms of entropy at standard conditions
        coeffs[6] = s0_by_R - ( coeffs[0]*np.log(self.T_std) + coeffs[1]*self.T_std + coeffs[2]*((1./2.)**(1./2.)*self.T_std)**2 \
                + coeffs[3]*((1./3.)**(1./3.)*self.T_std)**3 + coeffs[4]*((1./4.)**(1./4.)*self.T_std)**4  )
        # Define coeffs[6+7] by ensuring continuity at Tcommon
        s_sum_term_l = lambda T:  coeffs[0]*np.log(T) + coeffs[1]*T + coeffs[2]*((1./2.)**(1./2.)*T)**2 + coeffs[3]*((1./3.)**(1./3.)*T)**3 + coeffs[4]*((1./4.)**(1./4.)*T)**4 
        s_sum_term_h = lambda T:  coeffs[0+7]*np.log(T) + coeffs[1+7]*T + coeffs[2+7]*((1./2.)**(1./2.)*T)**2 + coeffs[3+7]*((1./3.)**(1./3.)*T)**3 + coeffs[4+7]*((1./4.)**(1./4.)*T)**4 
        coeffs[6+7] = s_sum_term_l(Tcommon) + coeffs[6] - s_sum_term_h(Tcommon) 
        
        return coeffs


    def fit_nasapolys_simultaneous(self, T0, Tc_i, cp_over_R, h_over_RT, s_over_R, cp0_by_R, dhf_by_R, s0_by_R, simultaneous_fit):
        """
        take notes from main notes
        Supports 1d-arrays only [M,] sized, where M refers to number of T0 array (sample) size.
        Tc_i is the common temperature index in the T0 array
        """
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



        print( '\nFitting cp/R, h/RT and s/R simultaneously:\n' )

        M = 5
        C = np.zeros((3*N,2*M))
        d = np.zeros(3*N)

        # cp equation
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
        C[i1:i2,0] = np.ones((1,Nl)) - self.T_std/T_low
        C[i1:i2,1] = (c1*T_low)    - (self.T_std/T_low)*(c1*self.T_std)   
        C[i1:i2,2] = (c2*T_low)**2 - (self.T_std/T_low)*(c2*self.T_std)**2
        C[i1:i2,3] = (c3*T_low)**3 - (self.T_std/T_low)*(c3*self.T_std)**3
        C[i1:i2,4] = (c4*T_low)**4 - (self.T_std/T_low)*(c4*self.T_std)**4

        i3=i2
        i4=i2+Nh
        C[i3:i4,0+M] = np.ones((1,Nh)) - self.T_std/T_high
        C[i3:i4,1+M] = (c1*T_high)    - (self.T_std/T_high)*(c1*self.T_std)
        C[i3:i4,2+M] = (c2*T_high)**2 - (self.T_std/T_high)*(c2*self.T_std)**2
        C[i3:i4,3+M] = (c3*T_high)**3 - (self.T_std/T_high)*(c3*self.T_std)**3
        C[i3:i4,4+M] = (c4*T_high)**4 - (self.T_std/T_high)*(c4*self.T_std)**4

        
        # s equation
        c2 = (1./2.)**(1./2.)
        c3 = (1./3.)**(1./3.)
        c4 = (1./4.)**(1./4.)

        i5=i4
        i6=i4+Nl
        C[i5:i6,0] = np.log(T_low/self.T_std)
        C[i5:i6,1] = T_low - self.T_std
        C[i5:i6,2] = (c2*T_low)**2 - (c2*self.T_std)**2
        C[i5:i6,3] = (c3*T_low)**3 - (c3*self.T_std)**3
        C[i5:i6,4] = (c4*T_low)**4 - (c4*self.T_std)**4
        
        i7=i6
        i8=i6+Nh
        C[i7:i8,0+M] = np.log(T_high/self.T_std)
        C[i7:i8,1+M] = T_high - self.T_std
        C[i7:i8,2+M] = (c2*T_high)**2 - (c2*self.T_std)**2
        C[i7:i8,3+M] = (c3*T_high)**3 - (c3*self.T_std)**3
        C[i7:i8,4+M] = (c4*T_high)**4 - (c4*self.T_std)**4
        
        #The right hand side
        d[:Nl]=cp_over_R_L
        d[Nl:Nl+Nh]=cp_over_R_H
        
        d[i1:i2] = h_over_RT_L - dhf_by_R/T_low
        d[i3:i4] = h_over_RT_H - dhf_by_R/T_high
        
        d[i5:i6] = s_over_R_L - s0_by_R
        d[i7:i8] = s_over_R_H - s0_by_R




        # Constraints
        Aeq=np.zeros((2,2*M))
        beq=np.zeros(2)
        
        # Note that you can add here constraints if e.g. maximum values are known etc.

        # cp_over_R : C0 continuity at Tcommon (presumably no knowledge on min/max/std values)
        for i in range(0,5):
            Aeq[0,i] = Tcommon**i
            Aeq[0,i+M] = -Tcommon**i

        # cp_over_R : C1 continuity at Tcommon
        Aeq[1,1] = 1.
        Aeq[1,2] = 2.*Tcommon
        Aeq[1,3] = ((3.**(1./2.))*Tcommon)**2
        Aeq[1,4] = ((4.**(1./3.))*Tcommon)**3
        Aeq[1,1+M] = -1.
        Aeq[1,2+M] = -2.*Tcommon
        Aeq[1,3+M] = -((3.**(1./2.))*Tcommon)**2
        Aeq[1,4+M] = -((4.**(1./3.))*Tcommon)**3
        
        # - h and C0 continuity is guaranteed after solution by analytical consideration
        # - if forcing the constraints to the linear solution (below), result are typically less good
        #   and the other coefficients do the job anyways as suggested by pen and paper

        '''       
        # h : C0 continuity 
        Aeq[2,0] = 1. - self.T_std/Tcommon
        Aeq[2,M] = -(1. - self.T_std/Tcommon)
        for i in range(1,5):
            cf = (1./(i+1))**(1./i)
            Aeq[2,i] = (cf*Tcommon)**i - (self.T_std/Tcommon)*(cf*self.T_std)**i  
            Aeq[2,i+M] = -((cf*Tcommon)**i - (self.T_std/Tcommon)*(cf*self.T_std)**i)  
        
        # s : C0 continuity 
        
        Aeq[3,0] = np.log(Tcommon/self.T_std)
        Aeq[3,M] = -np.log(Tcommon/self.T_std)
        for i in range(1,5):
            cf = (1./i)**(1./i)
            Aeq[3,i] = (cf*Tcommon)**i - (cf*self.T_std)**i
            Aeq[3,i+M] = -(cf*Tcommon)**i + (cf*self.T_std)**i
        '''


        # Find the least-squares solution
        sol=fitter.lsqlin(C, d, 0, None, None, Aeq, beq,-1e9,1e9,None,{'show_progress': False, 'abstol': 1e-12, 'reltol': 1e-8})
        coeffs_tmp=sol['x']

        # fill the coefficient matrix
        coeffs=np.zeros(14) 
        for i in range(0,5):
            coeffs[i] = coeffs_tmp[i]
            coeffs[i+7] = coeffs_tmp[i+M]
        
        # solving additional constant by means of conservation of enthalpy and entropy + ensuring C0 continuity at T=Tcommon 
        
        # Define coeff[5] in terms of the enthalpy of formation at standard conditions
        coeffs[5] = dhf_by_R - ( coeffs[0]*self.T_std + coeffs[1]*((1./2.)**(1./2.)*self.T_std)**2 \
                + coeffs[2]*((1./3.)**(1./3.)*self.T_std)**3 + coeffs[3]*((1./4.)**(1./4.)*self.T_std)**4 + coeffs[4]*((1./5.)**(1./5.)*self.T_std)**5 )
        
        # Define coeffs[5+7] by ensuring continuity at Tcommon
        h_sum_term_l = lambda T:  coeffs[0] + coeffs[1]*(1./2.)*T + coeffs[2]*((1./3.)**(1./2.)*T)**2 + coeffs[3]*((1./4.)**(1./3.)*T)**3 + coeffs[4]*((1./5.)**(1./4.)*T)**4 
        h_sum_term_h = lambda T:  coeffs[0+7] + coeffs[1+7]*(1./2.)*T + coeffs[2+7]*((1./3.)**(1./2.)*T)**2 + coeffs[3+7]*((1./4.)**(1./3.)*T)**3 + coeffs[4+7]*((1./5.)**(1./4.)*T)**4 
        coeffs[5+7] = Tcommon*( h_sum_term_l(Tcommon) + coeffs[5]/Tcommon - h_sum_term_h(Tcommon) )

        # Define coeff[6] in terms of entropy at standardconditions
        coeffs[6] = s0_by_R - ( coeffs[0]*np.log(self.T_std) + coeffs[1]*self.T_std + coeffs[2]*((1./2.)**(1./2.)*self.T_std)**2 \
                + coeffs[3]*((1./3.)**(1./3.)*self.T_std)**3 + coeffs[4]*((1./4.)**(1./4.)*self.T_std)**4  )
        # Define coeffs[6+7] by ensuring continuity at Tcommon
        s_sum_term_l = lambda T:  coeffs[0]*np.log(T) + coeffs[1]*T + coeffs[2]*((1./2.)**(1./2.)*T)**2 + coeffs[3]*((1./3.)**(1./3.)*T)**3 + coeffs[4]*((1./4.)**(1./4.)*T)**4 
        s_sum_term_h = lambda T:  coeffs[0+7]*np.log(T) + coeffs[1+7]*T + coeffs[2+7]*((1./2.)**(1./2.)*T)**2 + coeffs[3+7]*((1./3.)**(1./3.)*T)**3 + coeffs[4+7]*((1./4.)**(1./4.)*T)**4 

        coeffs[6+7] = s_sum_term_l(Tcommon) + coeffs[6] - s_sum_term_h(Tcommon) 
            
        return coeffs


    def refit_nasapolys(self, nasa_coeffs, Tcommon):
        """
        Refit NASA polynomials based on existing coefficients. This function should be
        utilised when a different Tcommon value is desired.
        nasa_coeffs: [M,] or [M,N] sized array following the cantera notation, where M is 15 (Tcommon + 2*7 coeffs)
        and N is the number of species
        Tcommon
        """        

        if(type(Tcommon) is not np.ndarray):
            Tcommon = Tcommon*np.ones(nasa_coeffs.shape[0])

        nasa_coeffs = np.atleast_2d(nasa_coeffs) # consistent behavior with 1d and 2d arrays
        N = nasa_coeffs.shape[0]

        for i in range(N):

            good_fit = self.check_quality(nasa_coeffs[i,:])
            Tc_shared = nasa_coeffs[i,0] == Tcommon[i]
            print(nasa_coeffs[i,0])
            print(good_fit)
            if( (not good_fit) or (not Tc_shared) ):
                print("gg")

        return "uli"