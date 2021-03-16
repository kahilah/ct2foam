import argparse
import numpy as np
import os

from thermo_transport import  ctBasedProperties as ctProps

if(__name__ == '__main__'):
    parser = argparse.ArgumentParser(description='Solve adiabatic free flame over multiple initial conditions.')
    parser.add_argument('-n','--ncores',  type=int, help='Number available cores for parallel execution.', default=1, required=False)
    parser.add_argument('-id','--rankID', type=int, help='Individual rank/task id.', default=0, required=False)
    args = parser.parse_args()


    


    thermo = ctProps.ctThermo("uli")
    thermo.sanity_checks()