import numpy as np

class nasaPolynomials:
    def __init__(self,   Tcommon=1000.0):
        """
        """
        self.mechanismFile =  mechanismFile
        self.outputDir =  outputDir

        self.T     = np.sort(T)


    def check_quality(self, coeffs):
        """
        check_quality checks whether the nasa polynomial coefficients fulfill
        the original requirements on C0 and C1 continuity.

        :coeffs: array with common temperature and low and high temperature coefficients
        :return: boolean whether a good fit or not.
        """ 


        # Check the quality of the existing fit (w.r.t Tcommon)
        low_T_coeffs = gas.species(gas.species_index(sp_i)).thermo.coeffs[8:15]
        high_T_coeffs = gas.species(gas.species_index(sp_i)).thermo.coeffs[1:8]
        Tcommon_orig = gas.species(gas.species_index(sp_i)).thermo.coeffs[0]
        coeffs[:7] = low_T_coeffs
        coeffs[7:14] = high_T_coeffs

        desired_Tcommon = coeffs[0] == Tcommon

        # relative error limits for consistency check are chosen by experience: 
        # - C0 continuity: 1e-6 error is acceptable 
        # - C1 continuityi (not mandatory but recommended): For some mechanism fits C1 continuity is not guaranteed --> a loose tolerance is chosen
        C0_at_Tc_cp = np.abs(cpFunc_low(Tcommon_orig) - cpFunc_high(Tcommon_orig))/np.abs(cpFunc_low(Tcommon_orig)) < 1e-6
        C1_at_Tc_cp = np.abs(cpdTFunc_low(Tcommon_orig) - cpdTFunc_high(Tcommon_orig))/np.abs(cpdTFunc_low(Tcommon_orig)) < 0.01 
        C0_at_Tc_h = np.abs(hFunc_low(Tcommon_orig) - hFunc_high(Tcommon_orig))/np.abs(hFunc_low(Tcommon_orig)) < 1e-6
        C0_at_Tc_s = np.abs(sFunc_low(Tcommon_orig) - sFunc_high(Tcommon_orig))/np.abs(sFunc_low(Tcommon_orig)) < 1e-6
        
        good_fit = C0_at_Tc_cp and C0_at_Tc_h and C0_at_Tc_s and C1_at_Tc_cp


    def refit_nasapoly(self, nasa_coeffs, Tcommon):
        
        check_quality
        if coeffs[0] = Tcommon
        fit 
        compare


        cp = np.atleast_2d(cp) # consistent behavior with 1d and 2d arrays
        N = cp.shape[0]

        As = np.zeros(N)
        Ts = np.zeros(N)
        std_err = np.zeros((N, 2))

        for i in range(N):
            cp_i = cp[i,:]
            popt, pcov = curve_fit(sutherland, T, mu_i, p0=p0)
