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

    #Write the thermo output in openfoam format
    with open(file_name,'a') as output:
        
        output.write(name+'\n{\n')
        output.write('\t')
        ############################################################################# 
        output.write('specie\n\t{\n')
        output.write('\t\tnMoles \t 1;\n')
        output.write('\t\tmolWeight \t'+str(MW)+';')
        output.write('\n\t}\n\n')
        #############################################################################

        ############################################################################# 
        output.write('\tthermodynamics')
        #############################################################################         
        output.write('\n\t{\n')
        output.write('\t\tTlow\t\t'+str(nasa7_Tlo)+';\n')
        output.write('\t\tThigh\t\t'+str(nasa7_Thi)+';\n')       
        output.write('\t\tTcommon\t\t'+str(nasa7_Tmid)+';\n')
        output.write('\t\tlowCpCoeffs\t(\t' )
        for wi in range(0,7): #NASA pol has 7 coeffs
            output.write(str(nasa7_lo[wi])) 
            output.write(' ')
        output.write(' );\n')
        output.write('\t\thighCpCoeffs\t(\t' )
        for wi in range(0,7):#NASA pol has 7 coeffs
            output.write(str(nasa7_hi[wi])) 
            output.write(' ')
        output.write(' );\n')
        output.write('\t}\n\n')
        
        ############################################################################# 
        output.write('\ttransport \n\t{\n')
        #############################################################################

        output.write('\t\tAs\t'+str(As)+';\n' )
        output.write('\t\tTs\t'+str(Ts)+';\n' )

        output.write('\t\tmuLogCoeffs<8>\t(\t' )
        for wi in range(8):
            if(wi < len(logpoly_mu_rev)):
                output.write(str(logpoly_mu_rev[wi])) 
            else:
                output.write("0") 
            output.write(' ')
        output.write(' );\n')

        output.write('\t\tmuCoeffs<8>\t(\t' )
        for wi in range(8):
            if(wi < len(poly_mu_rev)):
                output.write(str(poly_mu_rev[wi])) 
            else:
                output.write("0") 
            output.write(' ')
        output.write(' );\n')
        #############################################################################
        output.write('\t\tkappaLogCoeffs<8>\t(\t' )
        for wi in range(8):
            if(wi < len(logpoly_kappa_rev)):
                output.write(str(logpoly_kappa_rev[wi])) 
            else:
                output.write("0") 
            output.write(' ')
        output.write(' );\n')
        output.write('\t\tkappaCoeffs<8>\t(\t' )
        for wi in range(8):
            if(wi < len(poly_kappa_rev)):
                output.write(str(poly_kappa_rev[wi])) 
            else:
                output.write("0") 
            output.write(' ')
        output.write(' );\n')
        output.write('\t}\n\n')
        ############################################################################# 

        ############################################################################# 
        if(elements is not None):
            output.write('\telements')
            #############################################################################      
            output.write('\n\t{\n')   
            # Variables for the elemental composition entry:
            for elem_i in elements.keys():
                output.write("\t\t" + elem_i+"\t"+str(int(elements[elem_i]))+';\n')         
            output.write('\t}\n')
            #############################################################################
        output.write('}\n\n')    

