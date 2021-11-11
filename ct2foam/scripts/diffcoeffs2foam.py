import numpy as np
import cantera as ct
from matplotlib import pyplot as plt

from ct2foam.thermo_transport import ct_properties
from ct2foam.thermo_transport import foam_writer





def diffcoeffs_to_foam():

    mech = 'gri30.yaml'
    fuel = 'CH4'
    oxidant = 'O2:1, N2:3.76'
    phi = 1.0
    output = "Dm_coeffs.foam"

    T_array = np.linspace(280, 3000, 32)
    p_array = np.linspace(0.1e5, 10e5, 32)

    T0 = 300
    p_unity = 1
    reactants = fuel + ':1,' + oxidant 

    gas = ct.Solution(mech)
    gas.TPX = T0, p_unity, reactants
    gas.set_equivalence_ratio(phi, fuel, oxidant)

    _, Dm = ct_properties.evaluate_diff_coeffs(gas, Trange=T_array)
    Dm_matrix = ct_properties.Dm_TP_matrix(p_array, Dm)

    mixture_model = True
    foam_writer.write_Dm_dict(output, gas.species_names, T_array, p_array, Dm_matrix, mixture_model)

    # - plot important species for checking T dependence
    ref_p = 1e5
    plt.plot(T_array, Dm[gas.species_index("H2"), :] / ref_p, label='H2')
    plt.plot(T_array, Dm[gas.species_index("H"), :] / ref_p, label='H')
    plt.plot(T_array, Dm[gas.species_index("O"), :] / ref_p, label='O')
    plt.legend()
    plt.xlabel("T [K]")
    plt.ylabel("Dm")
    plt.show()

if(__name__ == '__main__'):
    diffcoeffs_to_foam()



