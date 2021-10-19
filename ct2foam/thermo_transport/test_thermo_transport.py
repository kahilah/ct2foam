import unittest
import numpy as np
import os
from pathlib import Path

from ct2foam.thermo_transport import ct_properties
from ct2foam.thermo_transport import  transport_fitter as tr_fitter
from ct2foam.thermo_transport import  thermo_fitter as th_fitter
from ct2foam.thermo_transport import ct2foam_utils 
from ct2foam.thermo_transport import  foam_writer 

"""
OF_referece/Test-thermoMixture.C is used as a source of the reference data tested here.
"""

eps = 1e-12
R = 8314.47006650545  # taken from openFoam -- differs slightly from standards
test_data_dir = Path(Path(__file__).parent.parent, "test_data")


class TestThermoTransport(unittest.TestCase):

    def test_sutherland0(self):
        T, As, Ts = 1, 1, 1
        mu = tr_fitter.sutherland(T, As, Ts)
        mu_foam = 0.5
        self.assertTrue(np.abs(mu-mu_foam) < eps)

    def test_euken0(self):
        T, As, Ts, W = 1, 1, 1, 2
        cv_mole = -R + 4
        kappa_foam = 936.697882481863 
        mu = tr_fitter.sutherland(T, As, Ts)
        kappa = tr_fitter.euken(mu, cv_mole, W, R)
        self.assertTrue(np.abs(kappa-kappa_foam) < eps)

    def test_sutherland1(self):
        T, As, Ts = 400, 1.67212e-06, 170.672  # H2O
        mu = tr_fitter.sutherland(T, As, Ts)
        mu_foam = 2.34407155073317e-05
        self.assertTrue(np.abs(mu-mu_foam) < eps)

    def test_euken1(self):
        T, As, Ts, W = 400, 1.67212e-06, 170.672, 18.0153  # H2O
        cv_mole = 26123.236960773
        kappa_foam = 0.0640159441308283
        mu = tr_fitter.sutherland(T, As, Ts)
        kappa = tr_fitter.euken(mu, cv_mole, W, R)
        self.assertTrue(np.abs(kappa-kappa_foam) < eps)

    def test_poly0(self):
        T = 400
        poly_coeffs_mu = np.flip(np.array([1000, -0.05, 0.003, 0]))
        poly_coeffs_kappa = np.flip(np.array([2000, -0.15, 0.023, 0]))
        mu_foam =  1460.0
        kappa_foam = 5620.0
        mu, kappa = tr_fitter.eval_polynomial(poly_coeffs_mu, poly_coeffs_kappa, T)
        self.assertTrue(np.abs(mu - mu_foam) / np.abs(mu_foam) < eps)
        self.assertTrue(np.abs(kappa - kappa_foam) / np.abs(kappa_foam) < eps)

    def test_logpoly0(self):
        T = 400
        poly_coeffs_mu = np.flip(np.array([0.1, 0.1, 0.1, 0]))
        poly_coeffs_kappa = np.flip(np.array([0.1, 0.1, 0.1, 0]))
        mu_foam =  72.8870655981874
        kappa_foam = 72.8870655981874
        mu, kappa = tr_fitter.eval_log_polynomial(poly_coeffs_mu, poly_coeffs_kappa, T)
        self.assertTrue(np.abs(mu - mu_foam) / np.abs(mu_foam) < eps)
        self.assertTrue(np.abs(kappa-kappa_foam)/np.abs(kappa_foam) < eps)

    def test_transport_fit(self):

        data = ct_properties.ctThermoTransport("h2o2.cti", verbose=False)
        data.evaluate_properties()

        transport_fits = ct2foam_utils.fit_ct_transport(data)        
        success = ct2foam_utils.transport_fit_quality(data, transport_fits, test_data_dir, False, 0.02, 0.08, 5e-3)
        self.assertTrue(success)
    
    def test_transport_sanity(self):
        """
        Reference values taken from NIST data base
        """
        T = 400
        cv_mole, W = 21005.045895231186, 28.014
        species_name = "N2"

        data = ct_properties.ctThermoTransport("gri30.cti", verbose=False)
        data.evaluate_properties()
        i = data.gas.species_index(species_name)

        As, Ts, _, poly_mu, poly_kappa, log_poly_mu, log_poly_kappa= ct2foam_utils.fit_ct_transport(data)

        mu_s = tr_fitter.sutherland(T, As[i], Ts[i])
        kappa_s = tr_fitter.euken(mu_s, cv_mole, W, R)
        mu_logp, kappa_logp = tr_fitter.eval_log_polynomial(log_poly_mu[i,:], log_poly_kappa[i,:], T)
        mu_p, kappa_p = tr_fitter.eval_polynomial(poly_mu[i,:], poly_kappa[i,:], T)


        # rough test whether they are in the same scale...
        mu_ref = 2.2217e-5
        kappa_ref = 0.032205

        self.assertTrue(np.abs(mu_s-mu_ref)/np.abs(mu_ref) < 0.07)
        self.assertTrue(np.abs(mu_p-mu_ref)/np.abs(mu_ref) < 0.01)
        self.assertTrue(np.abs(mu_logp-mu_ref)/np.abs(mu_ref) < 0.01)
        self.assertTrue(np.abs(kappa_s-kappa_ref)/np.abs(kappa_ref) < 0.05)
        self.assertTrue(np.abs(kappa_p-kappa_ref)/np.abs(kappa_ref) < 0.05)
        self.assertTrue(np.abs(kappa_logp-kappa_ref)/np.abs(kappa_ref) < 0.05)

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
        thermo = ct_properties.ctThermoTransport(Path(test_data_dir, "h2o2_mod.yaml"), verbose=False)
        nasa7 = thermo.is_nasa7(0)
        self.assertFalse(nasa7)
    
    def test_mixture_transport(self):
        mech = "gri30.cti"
        mix_name = "test"
        ct_mixture = "O2: 1, N2: 3.76"

        data = ct_properties.ctThermoTransport(mech, verbose=False)
        data.evaluate_mixture_properties(mix_name, ct_mixture)
        transport_fits = ct2foam_utils.fit_ct_transport(data)
        success = ct2foam_utils.transport_fit_quality(data, transport_fits, test_data_dir, plot=False, rel_tol_sutherland=2e-2, rel_tol_Euken=5e-2, rel_tol_poly=5e-3)

        # rough test whether they are in the right scale...
        Tref = 980
        mu_ref = 4.25e-5
        kappa_ref = 0.0685
        Ti = np.argmin(np.abs(data.T-Tref))
        poly_mu = transport_fits[3]
        poly_kappa = transport_fits[4]
        mu_p, kappa_p = tr_fitter.eval_polynomial(poly_mu[0,:], poly_kappa[0,:], data.T)

        self.assertTrue(success)
        self.assertTrue(np.abs(mu_p[Ti] - mu_ref)/np.abs(mu_ref) < 0.01)
        self.assertTrue(np.abs(kappa_p[Ti] - kappa_ref)/np.abs(kappa_ref) < 0.01)

    def test_mixture_thermo(self):
        mech = "gri30.cti"
        mix_name = "test"
        ct_mixture = "O2: 1, N2: 3.76"
        nasa7_Tmid = 1000.0
        data = ct_properties.ctThermoTransport(mech, verbose=False)
        data.evaluate_mixture_properties(mix_name, ct_mixture)
        thermo_fits = ct2foam_utils.fit_mixture_thermo(data)
        success = ct2foam_utils.nasa7_fit_quality(data, thermo_fits, test_data_dir, plot=False)

        cp_ref = 1.15e3
        cp = R * th_fitter.cp_nasa7(nasa7_Tmid, nasa7_Tmid, thermo_fits[0], thermo_fits[1]) / data.W

        self.assertTrue(success)
        self.assertTrue(np.abs(cp - cp_ref)/np.abs(cp_ref) < 0.01)

    def test_thermo_foam_writer(self):
        
        test_file = Path(test_data_dir, "testDict")
        foam_file_ref = Path(test_data_dir, "OF_reference", "refDict")

        name = "C2H2"
        W = 1.123
        As = 1.123
        Ts = 2.234
        poly_mu = np.flip([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
        poly_kappa = np.flip([0.9, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7])
        logpoly_mu = np.flip([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
        logpoly_kappa = np.flip([0.9, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7])
        nasa7_Tmid = 1000.0
        nasa7_Tlo, nasa7_Thi = 200, 5000
        nasa7_lo = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
        nasa7_hi = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
        elements = {"C": 2, "H": 2}
        foam_writer.write_thermo_transport(test_file, name, W, As, Ts, poly_mu, poly_kappa, logpoly_mu, logpoly_kappa, nasa7_Tmid, nasa7_Tlo, nasa7_Thi, nasa7_lo, nasa7_hi, elements=elements)
        value = os.system("diff -q " + str(test_file) + " " + str(foam_file_ref))
        self.assertTrue(value == 0)
        # clean-up
        test_file.unlink()

        # The following more elaborate test fails if ran on different architechtures (floating point differences)
        # Writes a dictionary which is used by test_data/OF_reference/Test_thermoMixture.C
        """
        from shutil import copyfile
        data = ct_properties.ctThermoTransport("h2o2.cti", verbose=False)
        data.evaluate_properties()
    
        transport_fits = ct2foam_utils.fit_ct_transport(data, poly_order=3)
        cwd = os.getcwd()
        thermo_fits = ct2foam_utils.refit_ct_thermo(data, data.Tmid, test_data_dir)

        foam_file = os.path.join( test_data_dir, "thermoDict_H2")
        foam_file_tmp = os.path.join(test_data_dir, "thermoDict_H2_orig")
        r_file = os.path.join(test_data_dir, "reactions.foam")
        sp_file = os.path.join(test_data_dir, "species.foam")
        copyfile(foam_file, foam_file_tmp)
        ct2foam_utils.ct2foam_thermo_writer(sp_file, foam_file, r_file, data, transport_fits, thermo_fits)
        value = os.system("diff -q " + repr(foam_file) + " " + repr(foam_file_tmp))
        self.assertTrue(value==0)
        # clean-up
        copyfile(foam_file_tmp, foam_file)
        os.remove(foam_file_tmp)
        os.remove(r_file)
        os.remove(sp_file)
        """



if __name__ == '__main__':
    unittest.main()