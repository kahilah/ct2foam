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
#include "logPolynomialTransport.H"

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
            janafThermo<perfectGas<specie>>,
            sensibleEnthalpy
        >
    > polynomialTransport;

    typedef logPolynomialTransport
    <
        species::thermo
        <
            janafThermo<perfectGas<specie>>,
            sensibleEnthalpy
        >
    > logPolynomialTransport;


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
    logPolynomialTransport tlp(dict.subDict("specie4"));
    Info<< "logmu = " << tlp.mu(1, T) << endl;
    Info<< "logkappa = " << tlp.kappa(1, T) << endl;


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


    Info << "\nH2" << endl;    
    sutherlandJanafTransport ts3(dict.subDict("H2"));
    T = 400;
    p = 1e5; //=Pstd to have fair comparison with references.
    Info<< "R = " << ts3.R() << endl;
    Info<< "cp = " << ts3.cp(p, T) << endl;
    Info<< "Cp = " << ts3.Cp(p, T) << endl;
    Info<< "cv = " << ts3.cv(p, T) << endl;
    Info<< "Cv = " << ts3.Cv(p, T) << endl;
    Info<< "mu = " << ts3.mu(p, T) << endl;
    Info<< "kappa = " << ts3.kappa(p, T) << endl;
    Info<< "h = " << ts3.ha(p, T) << endl;
    Info<< "s = " << ts3.s(p, T) << endl;


    // for testing python module writing functions
    dictionary dict_py(IFstream("thermoDict_H2")());
    Info << "\nH2 from python" << endl;    
    sutherlandJanafTransport ts_py(dict_py.subDict("H2"));
    Info<< "R = " << ts_py.R() << endl;
    Info<< "cp = " << ts_py.cp(p, T) << endl;
    Info<< "Cp = " << ts_py.Cp(p, T) << endl;
    Info<< "cv = " << ts_py.cv(p, T) << endl;
    Info<< "Cv = " << ts_py.Cv(p, T) << endl;
    Info<< "mu = " << ts_py.mu(p, T) << endl;
    Info<< "kappa = " << ts_py.kappa(p, T) << endl;
    Info<< "h = " << ts_py.ha(p, T) << endl;
    Info<< "s = " << ts_py.s(p, T) << endl;

    Info << "\nH2 polynomial transport from python" << endl;    
    polynomialTransport ts_py_poly(dict_py.subDict("H2"));
    Info<< "mu = " << ts_py_poly.mu(p, T) << endl;
    Info<< "kappa = " << ts_py_poly.kappa(p, T) << endl;

    Info << "\nH2 log-polynomial transport from python" << endl;    
    logPolynomialTransport ts_py_logpoly(dict_py.subDict("H2"));
    Info<< "mu = " << ts_py_logpoly.mu(p, T) << endl;
    Info<< "kappa = " << ts_py_logpoly.kappa(p, T) << endl;


    Info << "\nReference data against data written by the python module for H2" << endl;    

    Info<< "dR = " << mag(ts_py.R()-ts3.R())/ts3.R() << endl;
    Info<< "dmu = " << mag(ts_py.mu(p, T)-ts3.mu(p, T))/ts3.mu(p, T) << endl;
    Info<< "dkappa = " << mag(ts_py.kappa(p, T)-ts3.kappa(p, T))/ts3.kappa(p, T)  << endl;
    // Note, here sutherland vs poly due to lack of other reference. See NIST for actual absolute values
    Info<< "dmu_poly = " << mag(ts_py.mu(p, T)-ts_py_poly.mu(p, T))/ts_py_poly.mu(p, T) << endl;
    Info<< "dkappa_poly = " << mag(ts_py.kappa(p, T)-ts_py_poly.kappa(p, T))/ts_py_poly.kappa(p, T)  << endl;
    // here log-poly against normal poly for consistensy
    Info<< "dmu_logpoly = " << mag(ts_py_logpoly.mu(p, T)-ts_py_poly.mu(p, T))/ts_py_poly.mu(p, T) << endl;
    Info<< "dkappa_logpoly = " << mag(ts_py_logpoly.kappa(p, T)-ts_py_poly.kappa(p, T))/ts_py_poly.kappa(p, T)  << endl;

    Info<< "dcp = " << mag(ts_py.cp(p, T)-ts3.cp(p, T))/ts3.cp(p, T) << endl;
    Info<< "dCp = " << mag(ts_py.Cp(p, T)-ts3.Cp(p, T))/ts3.Cp(p, T) << endl;
    
    Info<< "\nEnd\n" << endl;


    return 0;
}