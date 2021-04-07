import cantera as ct
import numpy as np
import os

class ctThermoTransport:
    def __init__(self, mechanismFile, outputDir=None, T=np.linspace(280,3000,128), Tcommon=1000.0):
        """
        A class providing functions to access cantera mechanism and thermo and transport data required in polynomial fitting.
        - mechanismFile: path to file
        - outputDir: Where you want the mechanism file output to be placed.
        - T: sampling ratio of temperature for evaluating thermophysical properties for fitting.
        - Tcommon: utilised by thermo_fitter w.r.t NASA polynomial fitting process.
        """
        self.mechanismFile =  mechanismFile
        self.outputDir =  outputDir

        self.T_std = 298.15
        self.T     = np.sort(T)
        self.Tcommon  =   Tcommon
        # insert Tcommon into the T array sustaining sorting
        idx = self.T.searchsorted(self.Tcommon)
        self.T = np.concatenate((self.T[:idx], [self.Tcommon], self.T[idx:]))

        self.gas        = self.check_mechanism()

        self.R          = ct.gas_constant

        # sampling matrices of thermophysical variables (data-set per row)
        self.mu = np.zeros((self.gas.n_species, len(self.T)))
        self.kappa = np.zeros((self.gas.n_species, len(self.T)))
        self.cv_mass = np.zeros((self.gas.n_species, len(self.T)))
        self.cv_mole = np.zeros((self.gas.n_species, len(self.T)))
        self.cp = np.zeros((self.gas.n_species, len(self.T)))       # molar
        self.h = np.zeros((self.gas.n_species, len(self.T)))        # molar
        self.s = np.zeros((self.gas.n_species, len(self.T)))        # molar

        # properties at standard temperature, required by nasa polynomial fitting
        self.cp0_over_R = np.zeros(self.gas.n_species)
        self.dhf_over_R = np.zeros(self.gas.n_species)
        self.s0_over_R = np.zeros(self.gas.n_species)

        #self.test = dict.fromkeys(self.gas.species_names, np.zeros((len(T))))

        self.init_dirs()


    def check_mechanism(self):
        gas = ct.Solution(self.mechanismFile)
        return gas


    def init_dirs(self):

        if(self.outputDir == None):
            self.outputDir = os.getcwd()
        thermoDir = os.path.join(self.outputDir, 'foamThermoTransport')
        try:
            os.makedirs(thermoDir, exist_ok=True)
        except OSError:  
            None
        thermoFile= os.path.join(thermoDir, 'thermo.foam')
        try:
            os.remove(thermoFile) 
        except OSError:
            None
        # Add the chem.foam include to the file first
        #with open(thermoFN,'a') as output:
        #    output.write('#include "chem.foam"')
        #    output.write('\n')
        #This would write the species names list in OF style. However, we want to include that into the reaction file in pyJac style. 
        #See makeMechFileForOF.sh at pyjac folder
        #wr.write_sp_list(gas.species_names,thermoFN)
        # --------------------------------------------- #


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


    def nasa_normalisation(self, T, cp_mole, h_mole, s_mole):
        """
        Divide by gas constant (and T) according to NASA JANAF definitions.
        Assuming data structure [M,] or [M,N] shaped data, where M is the number of species and N is data size      
        """
        cp_over_R = cp_mole / ct.gas_constant
        h_over_RT = h_mole[None,:] / (T * ct.gas_constant)
        s_over_R  = s_mole / ct.gas_constant

        return cp_over_R, h_over_RT[0], s_over_R


    def evaluate_properties(self):
        '''
        Evaluate according to the temperature sampling, defined in class declaration.
        '''
        # Initialise a free flame object
        gas = ct.Solution(self.mechanismFile)
        p0 = ct.one_atm     # - not utilised
        gas.TP = 300.0, p0  # - not utilised

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
                self.cv_mass[i][j] = gas.cv_mass
                self.cv_mole[i][j] = gas.cv_mole


    def evaluate_mixture_properties(self, mixture_name, X):
        '''
        Evaluate according to the temperature sampling, defined in class declaration.
        X is molar concentration in cantera style formatting
        '''
        # Initialise a free flame object
        gas = ct.Solution(self.mechanismFile)
        p0 = ct.one_atm     # - not utilised
        gas.TPX = 300.0, p0, X 

        mu = np.zeros(len(self.T))
        kappa = np.zeros(len(self.T))
        cp = np.zeros(len(self.T))
        cv_mass = np.zeros(len(self.T))
        cv_mole = np.zeros(len(self.T))
        h = np.zeros(len(self.T))
        s = np.zeros(len(self.T))

        print("Evaluating thermophysical properties for mixture " + mixture_name + ": " + repr(X))

        # standard property, potentially required in cp re-fitting
        gas0 = ct.Solution(self.mechanismFile)
        gas0.TPX = 298.15, p0, X 
        cp0_over_R =  gas0.cp_mole/ct.gas_constant
        
        for i in range(len(self.T)):

            gas.TPX = self.T[i], p0, X

            mu[i] = gas.viscosity
            kappa[i] = gas.thermal_conductivity

            # divide by gas constant according to NASA JANAF definitions
            cp[i] = gas.cp_mole
            h[i]  = gas.enthalpy_mole
            s[i]  = gas.entropy_mole

            # both Cv definitions are required in Sutherland formulation
            cv_mass[i] = gas.cv_mass
            cv_mole[i] = gas.cv_mole


        mix_tran_data = {"mu": mu, "kappa": kappa}
        mix_thermo_data = {"cp": cp, "h": h, "s": s, "cv_mass": cv_mass, "cv_mole": cv_mole, "cp0_over_R": cp0_over_R}
        
        return mix_tran_data, mix_thermo_data


    def check_temperature_limits(self, sp_i):

        minT_ref = self.gas.species(self.gas.species_index(sp_i)).thermo.min_temp
        minT_user = np.min(self.T)
        maxT_ref = self.gas.species(self.gas.species_index(sp_i)).thermo.max_temp
        maxT_user = np.max(self.T)

        if( minT_ref > minT_user ):
            print("\t" + sp_i + ": Warning, given min(T)=" + repr(minT_user) +  "K is out of original bounds. (" + repr(minT_ref) + 'K)' )
        if( maxT_ref < maxT_user ):
            print("\t" + sp_i + ": Warning, given max(T)=" + repr(maxT_user) +  "K is out of original bounds. (" + repr(maxT_ref) + 'K)' )

