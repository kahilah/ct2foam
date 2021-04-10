/*---------------------------------------------------------------------------*\
Application
    Test-ThermoMixture

Description
    This program relies on source code and libraries of OpenFOAM-8, 
    which is licensed under the GNU General Public License, 
    see <http://www.gnu.org/licenses/>.

    - Source code is inspired by the OpenFOAM Test-thermoMixture test application.
\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "IFstream.H"
#include "specie.H"
#include "perfectGas.H"
#include "hConstThermo.H"
#include "janafThermo.H"
#include "sensibleEnthalpy.H"
#include "thermo.H"
#include "constTransport.H"
#include "sutherlandTransport.H"
#include "polynomialTransport.H"

#include <iomanip>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    typedef constTransport
    <
        species::thermo
        <
            hConstThermo<perfectGas<specie>>,
            sensibleEnthalpy
        >
    > constThermo;

    typedef sutherlandTransport
    <
        species::thermo
        <
            hConstThermo<perfectGas<specie>>,
            sensibleEnthalpy
        >
    > sutherlandThermo;

    typedef sutherlandTransport
    <
        species::thermo
        <
            janafThermo<perfectGas<specie>>,
            sensibleEnthalpy
        >
    > sutherlandJanafTransport;

    typedef polynomialTransport
    <
        species::thermo
        <
            hConstThermo<perfectGas<specie>>,
            sensibleEnthalpy
        >
    > polynomialTransport;


    dictionary dict(IFstream("thermoDict")());

    constThermo tc(dict.subDict("specie1"));
    Info << "Universal gas constant:" << endl;
    std::cout << std::setprecision(15) << tc.R()*tc.W() << '\n' << endl;

    Info << "\nConstant:" << endl;
    Info << "specie1" << endl;
    Info<< "R = " << tc.R() << endl;
    Info<< "cp = " << tc.cp(1, 1) << endl;
    Info<< "Cp = " << tc.Cp(1, 1) << endl;
    Info<< "cv = " << tc.cv(1, 1) << endl;
    Info<< "Cv = " << tc.Cv(1, 1) << endl;
    Info<< "mu = " << tc.mu(1, 1) << endl;
    Info<< "kappa = " << tc.kappa(1, 1) << endl;

    Info << "\nspecie2" << endl;
    constThermo tc2(dict.subDict("specie2"));
    Info<< "R = " << tc2.R() << endl;
    Info<< "cp = " << tc2.cp(1, 1) << endl;
    Info<< "Cp = " << tc2.Cp(1, 1) << endl;
    Info<< "cv = " << tc2.cv(1, 1) << endl;
    Info<< "Cv = " << tc2.Cv(1, 1) << endl;
    Info<< "mu = " << tc2.mu(1, 1) << endl;
    Info<< "kappa = " << tc2.kappa(1, 1) << endl;


    Info << "\nSutherland:" << endl;
    Info << "specie3" << endl;    
    sutherlandThermo ts(dict.subDict("specie3"));
    Info<< "R = " << ts.R() << endl;
    Info<< "cp = " << ts.cp(1, 1) << endl;
    Info<< "Cp = " << ts.Cp(1, 1) << endl;
    Info<< "cv = " << ts.cv(1, 1) << endl;
    Info<< "Cv = " << ts.Cv(1, 1) << endl;
    Info<< "mu = " << ts.mu(1, 1) << endl;
    Info<< "kappa = " << ts.kappa(1, 1) << endl;

    Info << "\nPolynomial:" << endl;
    Info << "specie4" << endl;    
    polynomialTransport tp(dict.subDict("specie4"));
    scalar T = 400;
    Info<< "R = " << tp.R() << endl;
    Info<< "cp = " << tp.cp(1, T) << endl;
    Info<< "Cp = " << tp.Cp(1, T) << endl;
    Info<< "cv = " << tp.cv(1, T) << endl;
    Info<< "Cv = " << tp.Cv(1, T) << endl;
    Info<< "mu = " << tp.mu(1, T) << endl;
    Info<< "kappa = " << tp.kappa(1, T) << endl;


    Info << "\nH2O" << endl;    
    sutherlandJanafTransport ts2(dict.subDict("H2O"));
    T = 400;
    scalar p = 1e5; //=Pstd to have fair comparison with references.
    Info<< "R = " << ts2.R() << endl;
    Info<< "cp = " << ts2.cp(p, T) << endl;
    Info<< "Cp = " << ts2.Cp(p, T) << endl;
    Info<< "cv = " << ts2.cv(p, T) << endl;
    Info<< "Cv = " << ts2.Cv(p, T) << endl;
    Info<< "mu = " << ts2.mu(p, T) << endl;
    Info<< "kappa = " << ts2.kappa(p, T) << endl;
    Info<< "h = " << ts2.ha(p, T) << endl;
    Info<< "s = " << ts2.s(p, T) << endl;
    Info<< "\nEnd\n" << endl;

    return 0;
}