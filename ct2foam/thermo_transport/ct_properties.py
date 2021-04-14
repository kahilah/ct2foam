import cantera as ct
import numpy as np
import os

class ctThermoTransport:
    def __init__(self, mechanismFile, outputDir=None, T=np.linspace(280,3000,128), Tmid=1000.0, verbose=True):
        """
        A class providing functions to access cantera mechanism and thermo and transport data required in polynomial fitting.
        - mechanismFile: path to file
        - outputDir: Where you want the mechanism file output to be placed.
        - T: sampling ratio of temperature for evaluating thermophysical properties for fitting.
        - Tcommon: utilised by thermo_fitter w.r.t NASA polynomial fitting process.
        """
        self.mechanismFile =  mechanismFile
        self.outputDir =  outputDir
        self.verbose = verbose

        self.T_std = 298.15
        self.T     = np.sort(T)
        self.Tmid  =   Tmid
        # insert Tcommon into the T array sustaining sorting
        idx = self.T.searchsorted(self.Tmid)
        self.T = np.concatenate((self.T[:idx], [self.Tmid], self.T[idx:]))

        self.gas        = self.check_mechanism()
        self.R          = ct.gas_constant
        self.W          = self.gas.molecular_weights
        self.names      = self.gas.species_names

        # sampling matrices of thermophysical variables (data-set per row)
        self.mu = np.zeros((self.gas.n_species, len(self.T)))
        self.kappa = np.zeros((self.gas.n_species, len(self.T)))
        self.cv_mole = np.zeros((self.gas.n_species, len(self.T)))
        self.cp = np.zeros((self.gas.n_species, len(self.T)))       # molar
        self.h = np.zeros((self.gas.n_species, len(self.T)))        # molar
        self.s = np.zeros((self.gas.n_species, len(self.T)))        # molar

        # properties at standard temperature, required by nasa polynomial fitting
        self.cp0_over_R = np.zeros(self.gas.n_species)
        self.dhf_over_R = np.zeros(self.gas.n_species)
        self.s0_over_R = np.zeros(self.gas.n_species)


    def check_mechanism(self):
        "Should catch any error related to mechanism."
        gas = ct.Solution(self.mechanismFile)
        return gas


    def get_thermo_fit_type(self, species_name):
        """
        sp_i: species name
        return: - type_name = [NasaPoly2 (NASA7), Nasa9PolyMultiTempRegion (NASA9), Shomate], 
        """
        i = self.gas.species_index(species_name)
        type_name = type(self.gas.species(i).thermo).__name__
        return type_name


    def get_nasa7_coeffs(self, species_name):
        """
        gas: cantera gas object
        sp_i: species name
        return: 
                - coeffs is the corresponding coeffcients for the analytical polynomial.
        """
        i = self.gas.species_index(species_name)
        return self.gas.species(i).thermo.coeffs



    def is_nasa7(self, species_name):
        fit_type = self.get_thermo_fit_type(species_name)
        if(fit_type=="NasaPoly2" ):
            return True
        else:
            return False


    def evaluate_properties(self):
        '''
        Evaluate according to the temperature sampling, defined in class declaration.
        '''
        # Initialise a free flame object
        gas = ct.Solution(self.mechanismFile)
        p0 = ct.one_atm     # - not utilised
        gas.TP = 300.0, p0  # - not utilised
        gas.transport_model = 'Multi'

        if(self.verbose):
            print("Evaluating thermophysical properties over species.")

        for sp_i in gas.species_names:
            
            self.check_temperature_limits(sp_i)

            reactants = sp_i + ':1.0'
            i = gas.species_index(sp_i)

            # standard properties
            self.cp0_over_R[i] =  gas.species(i).thermo.cp(self.T_std)/ct.gas_constant # this calls molar cp <==> consistent
            self.dhf_over_R[i] =  gas.species(i).thermo.h(self.T_std)/ct.gas_constant
            self.s0_over_R[i]  =  gas.species(i).thermo.s(self.T_std)/ct.gas_constant

            for j in range(len(self.T)):

                gas.TPX = self.T[j], p0, reactants

                self.mu[i][j] = gas.viscosity
                self.kappa[i][j] = gas.thermal_conductivity

                self.cp[i][j] = gas.cp_mole
                self.h[i][j]  = gas.enthalpy_mole
                self.s[i][j]  = gas.entropy_mole

                # both Cv definitions are required in Sutherland formulation
                self.cv_mole[i][j] = gas.cv_mole

        return 0

    def evaluate_mixture_properties(self, mixture_name, X):
        '''
        Evaluate according to the temperature sampling, defined in class declaration.
        X is molar concentration in cantera style formatting
        '''
        # Initialise a free flame object
        gas = ct.Solution(self.mechanismFile)
        p0 = ct.one_atm     # - not utilised
        gas.TPX = 300.0, p0, X 
        gas.transport_model = 'Multi'

        # return 2d arrays to be consistent with the full set of species
        mu = np.atleast_2d(np.zeros(len(self.T)))
        kappa = np.atleast_2d(np.zeros(len(self.T)))
        cp = np.atleast_2d(np.zeros(len(self.T)))
        cv_mole = np.atleast_2d(np.zeros(len(self.T)))
        h = np.atleast_2d(np.zeros(len(self.T)))
        s = np.atleast_2d(np.zeros(len(self.T)))

        if(self.verbose):
            print("Evaluating thermophysical properties for mixture " + mixture_name + ": " + repr(X))

        # standard property, potentially required in cp re-fitting
        gas0 = ct.Solution(self.mechanismFile)
        gas0.TPX = self.T_std, p0, X 
        cp0_over_R = np.atleast_1d(gas0.cp_mole/ct.gas_constant)
        dhf_over_R = np.atleast_1d(gas0.enthalpy_mole/ct.gas_constant)
        s0_over_R = np.atleast_1d(gas0.entropy_mole/ct.gas_constant)
        W = np.atleast_1d(gas0.mean_molecular_weight)

        for i in range(len(self.T)):

            gas.TPX = self.T[i], p0, X

            mu[0][i] = gas.viscosity
            kappa[0][i] = gas.thermal_conductivity

            # divide by gas constant according to NASA JANAF definitions
            cp[0][i] = gas.cp_mole
            h[0][i]  = gas.enthalpy_mole
            s[0][i]  = gas.entropy_mole

            # both Cv definitions are required in Sutherland formulation
            cv_mole[0][i] = gas.cv_mole

        #override the class initialisation
        self.W          = W
        self.names      = [mixture_name]
        self.mu = mu
        self.kappa = kappa
        self.cv_mole = cv_mole
        self.cp = cp
        self.h = h
        self.s = s
        self.cp0_over_R = cp0_over_R
        self.dhf_over_R = dhf_over_R
        self.s0_over_R = s0_over_R
        
        return 0


    def check_temperature_limits(self, sp_i):

        minT_ref = self.gas.species(self.gas.species_index(sp_i)).thermo.min_temp
        minT_user = np.min(self.T)
        maxT_ref = self.gas.species(self.gas.species_index(sp_i)).thermo.max_temp
        maxT_user = np.max(self.T)

        if(self.verbose):
            if( minT_ref > minT_user ):
                print("\t" + sp_i + ": Warning, given min(T)=" + repr(minT_user) +  "K is out of original bounds. (" + repr(minT_ref) + 'K)' )
            if( maxT_ref < maxT_user ):
                print("\t" + sp_i + ": Warning, given max(T)=" + repr(maxT_user) +  "K is out of original bounds. (" + repr(maxT_ref) + 'K)' )

