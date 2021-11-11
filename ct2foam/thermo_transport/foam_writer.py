import sys
import numpy as np


def write_species_list(file_name, species_names):
    """
    Writes a species.foam file with an OpenFOAM formatted list of species.
    """
    with open(file_name, 'a') as output:
        output.write('species\n(\n')
        for sp_i in species_names:
            output.write('\t')
            output.write(sp_i)
            output.write('\n')
        output.write(');\n\n')


def write_reactions(file_name):
    """
    Writes an empty reactions.foam file with an OpenFOAM formatted list of 0 reactions.
    """
    with open(file_name, 'a') as output:
        output.write('// - This dictionary intentionally left blank.\n')
        output.write('reactions\n{\n')
        output.write('}\n\n')


def write_thermo_transport(file_name, name, MW, As, Ts, poly_mu, poly_kappa, logpoly_mu, logpoly_kappa, nasa7_Tmid, nasa7_Tlo, nasa7_Thi, nasa7_lo, nasa7_hi, elements=None):
    """
    Writes thermophysicalProperties file required dictionary entries according to given data.
    file_name: output file name.
    name: species name / mixture name
    MW: molecular weight
    As, Ts: Sutherland entries
    poly*: transport polynomial entries
    nasa7*: NASA7 polynomial coefficient entry data
    elements: a dictionary with the following syntax:  elements = {"C": 1, "H":1}
    """
    # Numpy-based polynomial fit has a reversed order to OpenFoam dictionary definition
    poly_mu_rev = np.copy(poly_mu)
    poly_mu_rev = np.flip(poly_mu_rev)
    logpoly_mu_rev = np.copy(logpoly_mu)
    logpoly_mu_rev = np.flip(logpoly_mu_rev)

    poly_kappa_rev = np.copy(poly_kappa)
    poly_kappa_rev = np.flip(poly_kappa_rev)
    logpoly_kappa_rev = np.copy(logpoly_kappa)
    logpoly_kappa_rev = np.flip(logpoly_kappa_rev)

    # Write the thermo output in openfoam format
    with open(file_name, 'a') as output:

        output.write(name + '\n{\n')
        output.write('\t')

        output.write('specie\n\t{\n')
        output.write('\t\tnMoles \t 1;\n')
        output.write('\t\tmolWeight \t'+repr(MW) + ';')
        output.write('\n\t}\n\n')

        output.write('\tthermodynamics')
        output.write('\n\t{\n')
        output.write('\t\tTlow\t\t' + repr(nasa7_Tlo) + ';\n')
        output.write('\t\tThigh\t\t' + repr(nasa7_Thi) + ';\n')
        output.write('\t\tTcommon\t\t' + repr(nasa7_Tmid) + ';\n')
        output.write('\t\tlowCpCoeffs\t(\t')
        # - NASA pol has 7 coeffs
        for wi in range(0, 7):
            output.write(repr(nasa7_lo[wi]))
            output.write(' ')
        output.write(' );\n')
        output.write('\t\thighCpCoeffs\t(\t')
        # - NASA pol has 7 coeffs
        for wi in range(0, 7):
            output.write(repr(nasa7_hi[wi]))
            output.write(' ')
        output.write(' );\n')
        output.write('\t}\n\n')

        output.write('\ttransport \n\t{\n')

        output.write('\t\tAs\t' + repr(As) + ';\n')
        output.write('\t\tTs\t' + repr(Ts) + ';\n')

        output.write('\t\tmuLogCoeffs<8>\t(\t')
        for wi in range(8):
            if(wi < len(logpoly_mu_rev)):
                output.write(repr(logpoly_mu_rev[wi]))
            else:
                output.write("0")
            output.write(' ')
        output.write(' );\n')

        output.write('\t\tmuCoeffs<8>\t(\t')
        for wi in range(8):
            if(wi < len(poly_mu_rev)):
                output.write(repr(poly_mu_rev[wi]))
            else:
                output.write("0") 
            output.write(' ')
        output.write(' );\n')

        output.write('\t\tkappaLogCoeffs<8>\t(\t')
        for wi in range(8):
            if(wi < len(logpoly_kappa_rev)):
                output.write(repr(logpoly_kappa_rev[wi]))
            else:
                output.write("0")
            output.write(' ')
        output.write(' );\n')
        output.write('\t\tkappaCoeffs<8>\t(\t')
        for wi in range(8):
            if(wi < len(poly_kappa_rev)):
                output.write(repr(poly_kappa_rev[wi]))
            else:
                output.write("0")
            output.write(' ')
        output.write(' );\n')
        output.write('\t}\n\n')

        if(elements is not None):
            output.write('\telements')
            output.write('\n\t{\n')
            # Variables for the elemental composition entry:
            for elem_i in elements.keys():
                output.write("\t\t" + elem_i + "\t" + repr(int(elements[elem_i])) + ';\n')
            output.write('\t}\n')
        output.write('}\n\n')

def create_uniformtable2_dict(name, low, high, values, indent="\t"):
    """
    Returns OpenFOAM dictionary input for UniformTable2 entry.
    Implementation follows the corresponding OpenFOAM implementation.
    /src/OpenFOAM/primitives/functions/Function2/UniformTable2/UniformTable2.H
    Input:
        low: 2-sized tuple or array corresponding to x and y min values.
        high: 2-sized tuple or array corresponding to x and y max values.
        values: [m,n] array where m and n refer to x and y array sizes.
        indent: indentation prior to the dictionary entry.
    Return: Dictionary string.
    """

    m, n = values.shape

    indent1 = indent + "\t"
    indent2 = indent1 + "\t"

    of_entry = indent + name + "\n"
    of_entry += indent + "{\n"
    of_entry += indent1 + "type \tuniformTable;\n"
    of_entry += indent1 + "low \t(" + repr(low[0]) + " " + repr(low[1]) + ");\n"
    of_entry += indent1 + "high \t(" + repr(high[0]) + " " + repr(high[1]) + ");\n"
    of_entry += indent1 + "values\n"
    of_entry += indent1 + repr(m) + " " + repr(n) + "\n"
    of_entry += indent1 + "(\n"
    for mi in range(m):
        of_entry += indent2 + "("
        for ni in range(n):
            of_entry += " " + repr(values[mi, ni])
        of_entry += " )\n"
    of_entry += indent1 + ");\n"
    of_entry += indent + "}\n"

    return of_entry

def write_Dm_dict(file, species_names, T, p, Dm, mixture_model=True):
    """
    Writes the "Dm" sub-dictionaries for multi-component transport
    models defined in thermophysicalTransport dictionary in OpenFOAM.
    See an example input below:
    laminar
    {
        model 	FickianFourier;
        mixtureDiffusionCoefficients yes;
        Dm // [m^2/s]
        {
            CH4
            {
                type 	uniformTable;
                low 	(10000.0 280.0);
                high 	(300000.0 3000.0);
                values
                64 128
                (
                    (D00 D01 D03 ...)
                    (D10 D11 D13 ...)
                    (D20 D21 D23 ...)
                )
            }
            .
            .
            .
        }
    }

    Input:
        file: output file path
        species_names: k-sized array of species names.
        T: n-sized array of temperature values.
        p: m-sized array of temperature values.
        Dm: [m,n] array where m and n refer to p and T array sizes.
        mixture_model: whether to use mixture-averaged formulation or not.
    """

    n = len(T)
    m = len(p)

    if((m != Dm.shape[1]) or (n != Dm.shape[2])):
        sys.exit("Error in write_Dm_dict(): input array shape mismatch.")

    with open(file, 'w', encoding='utf8') as output:

        if(mixture_model):
            output.write("\tmixtureDiffusionCoefficients yes;\n")
        else:
            output.write("\tmixtureDiffusionCoefficients no;\n")

        output.write("\tDm // [m^2/s]\n\t{\n")
        for i, spi in enumerate(species_names):
            Dmi = Dm[i]
            Dm_func = create_uniformtable2_dict(
                spi, (p[0], T[0]), (p[-1], T[-1]), Dmi, indent="\t\t")
            output.write(Dm_func)

        output.write('\t}\n')
        output.write('}')
