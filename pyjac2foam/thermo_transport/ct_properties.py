import cantera as ct
import numpy as np
import os
import bisect

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
        self.cp = np.zeros((self.gas.n_species, len(self.T)))
        self.cv_mass = np.zeros((self.gas.n_species, len(self.T)))
        self.cv_mole = np.zeros((self.gas.n_species, len(self.T)))
        self.h = np.zeros((self.gas.n_species, len(self.T)))
        self.s = np.zeros((self.gas.n_species, len(self.T)))

        self.nasa_coeffs = np.zeros((self.gas.n_species, 15))

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

    

    def evaluate_properties(self):
        '''
        Evaluate according to the temperature sampling, defined in class declaration.
        '''
        # Initialise a free flame object
        gas = ct.Solution(self.mechanismFile)
        p0 = ct.one_atm     # - not utilised
        width = 0.03        # - not utilised
        gas.TP = 300.0, p0  # - not utilised
        f = ct.FreeFlame(gas, width=width)

        print("Evaluating thermophysical properties over species.")

        for sp_i in gas.species_names:
            
            self.check_temperature_limits(sp_i)

            reactants = sp_i + ':1.0'
            i = gas.species_index(sp_i)

            self.nasa_coeffs[i,:] = self.gas.species(gas.species_index(sp_i)).thermo.coeffs

            for j in range(len(self.T)):

                gas.TPX = self.T[j], p0, reactants
                f = ct.FreeFlame(gas, width=width)

                self.mu[i][j] = gas.viscosity
                self.kappa[i][j] = gas.thermal_conductivity

                # divide by gas constant according to NASA JANAF definitions
                self.cp[i][j] = gas.cp_mole/(ct.gas_constant)
                self.h[i][j]  = gas.enthalpy_mole/(ct.gas_constant*self.T[j])
                self.s[i][j]  = gas.entropy_mole/(ct.gas_constant)

                # both Cv definitions are required in Sutherland formulation
                self.cv_mass[i][j] = gas.cv_mass
                self.cv_mole[i][j] = gas.cv_mole
    

    def evaluate_mixture_properties(self,mixture_name, X):
        '''
        Evaluate according to the temperature sampling, defined in class declaration.
        X is molar concentration in cantera style formatting
        '''
        # Initialise a free flame object
        gas = ct.Solution(self.mechanismFile)
        p0 = ct.one_atm     # - not utilised
        width = 0.03        # - not utilised
        gas.TPX = 300.0, p0, X 
        f = ct.FreeFlame(gas, width=width)

        mu = np.zeros(len(self.T))
        kappa = np.zeros(len(self.T))
        cp = np.zeros(len(self.T))
        cv_mass = np.zeros(len(self.T))
        cv_mole = np.zeros(len(self.T))
        h = np.zeros(len(self.T))
        s = np.zeros(len(self.T))

        print("Evaluating thermophysical properties for mixture " + mixture_name + ": " + repr(X))
        
        for i in range(len(self.T)):

            gas.TPX = self.T[i], p0, X
            f = ct.FreeFlame(gas, width=width)

            mu[i] = gas.viscosity
            kappa[i] = gas.thermal_conductivity

            # divide by gas constant according to NASA JANAF definitions
            cp[i] = gas.cp_mole/(ct.gas_constant)
            h[i]  = gas.enthalpy_mole/(ct.gas_constant*self.T[i])
            s[i]  = gas.entropy_mole/(ct.gas_constant)

            # both Cv definitions are required in Sutherland formulation
            cv_mass[i] = gas.cv_mass
            cv_mole[i] = gas.cv_mole


        mix_tran_data = {"mu": mu, "kappa": kappa}
        mix_thermo_data = {"cp": cp, "h": h, "s": s, "cv_mass": cv_mass, "cv_mole": cv_mole}
        
        return mix_tran_data, mix_thermo_data


    def check_temperature_limits(self, sp_i):

        minT_ref = self.gas.species(self.gas.species_index(sp_i)).thermo.min_temp
        minT_user = np.min(self.T)
        maxT_ref = self.gas.species(self.gas.species_index(sp_i)).thermo.max_temp
        maxT_user = np.max(self.T)

        if( minT_ref > minT_user ):
            print("\tWARNING: Given min(T) for evaluation is out of original bounds:")
            print("\t\t" + sp_i + " : " + repr(minT_ref) + ' K (orig.) vs. ' + repr(minT_user) + ' K' )
        if( maxT_ref < maxT_user ):
            print("\tWARNING: Given max(T) for evaluation is out of original bounds:")
            print("\t\t" + sp_i +  " : " + repr(maxT_ref) + ' K (orig.) vs. ' + repr(maxT_user) + ' K' )

