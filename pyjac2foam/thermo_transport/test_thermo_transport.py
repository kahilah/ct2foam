import unittest
import numpy as np
from . import  ct_properties
from . import  transport_fitter as tr_fitter
from . import  thermo_fitter as th_fitter
from . import  ct2foam_utils 

"""
OF_referece/Test-thermoMixture.C is used as a source of the reference data tested here.
"""
eps = 1e-12
R = 8314.47006650545 # taken from openFoam -- differs from standards

class TestThermoTransport(unittest.TestCase):

    def test_sutherland0(self):
        T, As, Ts = 1, 1, 1
        mu=tr_fitter.sutherland(T, As, Ts)
        mu_foam = 0.5
        self.assertTrue(np.abs(mu-mu_foam) < eps)

    def test_euken0(self):
        T, As, Ts, W = 1, 1, 1, 2
        cv_mole = -R + 4
        kappa_foam = 936.697882481863 
        mu=tr_fitter.sutherland(T, As, Ts)
        kappa=tr_fitter.euken(mu, cv_mole, W, R)
        self.assertTrue(np.abs(kappa-kappa_foam) < eps)

    def test_sutherland1(self):
        T, As, Ts = 400, 1.67212e-06, 170.672 #H2O
        mu=tr_fitter.sutherland(T, As, Ts)
        mu_foam = 2.34407155073317e-05
        self.assertTrue(np.abs(mu-mu_foam) < eps)

    def test_euken1(self):
        T, As, Ts, W = 400, 1.67212e-06, 170.672, 18.0153 #H2O
        cv_mole = 26123.236960773
        kappa_foam = 0.0640159441308283
        mu=tr_fitter.sutherland(T, As, Ts)
        kappa=tr_fitter.euken(mu, cv_mole, W, R)
        self.assertTrue(np.abs(kappa-kappa_foam) < eps)

    def test_poly0(self):
        T = 400
        poly_coeffs_mu = np.flip(np.array([1000, -0.05, 0.003, 0]))
        poly_coeffs_kappa = np.flip(np.array([2000, -0.15, 0.023, 0]))
        mu_foam =  1460.0
        kappa_foam = 5620.0
        mu, kappa = tr_fitter.eval_polynomial(poly_coeffs_mu, poly_coeffs_kappa, T, poly_order=3, logBasis=False)
        self.assertTrue(np.abs(mu-mu_foam)/np.abs(mu_foam) < eps)
        self.assertTrue(np.abs(kappa-kappa_foam)/np.abs(kappa_foam) < eps)

    def test_transport_fit(self):
        T = 400
        thermo = ct_properties.ctThermoTransport("h2o2.cti", verbose=False)
        thermo.evaluate_properties()

        As, Ts, std_err, poly_mu, poly_kappa = ct2foam_utils.fit_ct_transport(thermo)        
        success = ct2foam_utils.transport_fit_quality(thermo.gas.species_names, thermo.T, thermo.mu, thermo.kappa, poly_mu, poly_kappa, As, Ts, std_err, thermo.cv_mole, thermo.W, False, 0.02, 0.08, 5e-3)
        self.assertTrue(success)

    def test_transport_sanity(self):
        """
        Reference values taken from NIST data base
        """
        T = 400
        cv_mole, W = 21005.045895231186, 28.014
        species_name = "N2"

        thermo = ct_properties.ctThermoTransport("gri30.cti", verbose=False)
        thermo.evaluate_properties()
        i = thermo.gas.species_index(species_name)

        As, Ts, _, poly_mu, poly_kappa = ct2foam_utils.fit_ct_transport(thermo)        
        mu_s = tr_fitter.sutherland(T, As[i], Ts[i])
        kappa_s=tr_fitter.euken(mu_s, cv_mole, W, R)
        mu_p, kappa_p = tr_fitter.eval_polynomial(poly_mu[:,i], poly_kappa[:,i], T, poly_order=3, logBasis=True)


        # rough test whether they are in the same scale...
        mu_ref = 2.2217e-5
        kappa_ref = 0.032205

        self.assertTrue(np.abs(mu_s-mu_ref)/np.abs(mu_ref) < 0.07)
        self.assertTrue(np.abs(mu_p-mu_ref)/np.abs(mu_ref) < 0.01)
        self.assertTrue(np.abs(kappa_s-kappa_ref)/np.abs(kappa_ref) < 0.05)
        self.assertTrue(np.abs(kappa_p-kappa_ref)/np.abs(kappa_ref) < 0.05)

    def test_nasa7(self):
        T = 400.0
        # H2O taken from OF8 tutorials
        Tmid = 1000.0
        c_lo = np.array((3.38684, 0.00347498, -6.3547e-06, 6.96858e-09, -2.50659e-12, -30208.1, 2.59023))
        c_hi = np.array((2.67215, 0.00305629, -8.73026e-07, 1.201e-10, -6.39162e-15, -29899.2, 6.86282 ))
        cp = R * th_fitter.cp_nasa7(T, Tmid, c_lo, c_hi)
        h = R * T * th_fitter.h_nasa7(T, Tmid, c_lo, c_hi)
        s = R * th_fitter.s_nasa7(T, Tmid, c_lo, c_hi)
        cp_foam = 34437.7070272785
        h_foam = -238388055.112524
        s_foam = 198687.554172766
        self.assertTrue(np.abs(cp-cp_foam)/np.abs(cp_foam) < eps)
        self.assertTrue(np.abs(h-h_foam)/np.abs(h_foam) < eps)
        self.assertTrue(np.abs(s-s_foam)/np.abs(s_foam) < eps)

    def test_nasa9(self):
        import os
        cwd=os.getcwd()
        thermo = ct_properties.ctThermoTransport(os.path.join(cwd, "test_data/h2o2_mod.yaml"), verbose=False)
        nasa7 = thermo.is_nasa7(0)
        self.assertFalse(nasa7)

    def test_mixture_transport(self):
        mech = "gri30.cti"
        mix_name = "test"
        ct_mixture = "O2: 1, N2: 3.76"

        thermo = ct_properties.ctThermoTransport(mech, verbose=False)
        mix_data = thermo.evaluate_mixture_properties(mix_name, ct_mixture)

        As, Ts, std_err = tr_fitter.fit_sutherland(thermo.T, mix_data[0]["mu"])
        poly_mu, poly_kappa = tr_fitter.fit_polynomial(thermo.T, mix_data[0]["mu"], mix_data[0]["kappa"])
        success = ct2foam_utils.transport_fit_quality(
            np.array([mix_name]), thermo.T,  mix_data[0]["mu"],  mix_data[0]["kappa"],
            poly_mu, poly_kappa, As, Ts, std_err, mix_data[1]["cv_mole"], mix_data[0]["W"],
            plot=False, rel_tol_sutherland=2e-2, rel_tol_Euken=5e-2, rel_tol_poly=5e-3
            )
        # rough test whether they are in the right scale...
        Tref = 980
        mu_ref = 4.25e-5
        kappa_ref = 0.0685
        Ti = np.argmin(np.abs(thermo.T-Tref))
        mu_p, kappa_p = tr_fitter.eval_polynomial(poly_mu[:,0], poly_kappa[:,0], thermo.T, poly_order=3, logBasis=True)

        self.assertTrue(success)
        self.assertTrue(np.abs(mu_p[Ti] - mu_ref)/np.abs(mu_ref) < 0.01)
        self.assertTrue(np.abs(kappa_p[Ti] - kappa_ref)/np.abs(kappa_ref) < 0.01)

    def test_mixture_thermo(self):
        mech = "gri30.cti"
        mix_name = "test"
        ct_mixture = "O2: 1, N2: 3.76"
        nasa7_Tmid = 1000.0
        thermo = ct_properties.ctThermoTransport(mech, verbose=False)
        mix_data = thermo.evaluate_mixture_properties(mix_name, ct_mixture)
        nasa_coeffs_lo, nasa_coeffs_hi = ct2foam_utils.fit_mixture_thermo(
            thermo.T, nasa7_Tmid, mix_data[1]["cp"], mix_data[1]["h"], mix_data[1]["s"], 
            mix_data[1]["cp0_over_R"], mix_data[1]["dhf_over_R"], mix_data[1]["s0_over_R"]
            )
        success = ct2foam_utils.nasa7_fit_quality(
            np.array([mix_name]), thermo.T, mix_data[1]["cp"], mix_data[1]["h"], mix_data[1]["s"],
            nasa7_Tmid, nasa_coeffs_lo, nasa_coeffs_hi, plot=False)

        cp_ref = 1.15e3
        cp = R * th_fitter.cp_nasa7(nasa7_Tmid, nasa7_Tmid, nasa_coeffs_lo, nasa_coeffs_hi) / mix_data[0]["W"]

        self.assertTrue(success)
        self.assertTrue(np.abs(cp - cp_ref)/np.abs(cp_ref) < 0.01)



if __name__ == '__main__':
    unittest.main()