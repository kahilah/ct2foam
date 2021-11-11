import argparse
from pathlib import Path
from ct2foam.thermo_transport import ct_properties
from ct2foam.thermo_transport import ct2foam_utils as utils
from ct2foam.thermo_transport import warning_msg


def init_dirs(output_dir=None):
    if(output_dir is None):
        output_dir = Path.cwd()
        print("Using current directory as an output directory.")
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    return output_dir


def main():
    parser = argparse.ArgumentParser(
        description='Convert/refit cantera .cti -based transport ' +
        'and thermodynamic data into OpenFOAM format.')
    parser.add_argument('-i', '--input', type=str,
                        help='Mechanism path.', default=None, required=True)
    parser.add_argument(
        '-o', '--output', type=str, help='Output directory path.',
        default=None, required=False)
    parser.add_argument(
        '-n', '--name', type=str, help='Mixture name, e.g. "air".',
        default=None, required=True)
    parser.add_argument(
        '-m', '--mixture', type=str,
        help='Molecular mixture ratio in cantera style: "O2:1, N2:3.76" ',
        default=None, required=True)
    parser.add_argument(
        '-T', '--Tmid', type=float,
        help='Common temperature for NASA-7 thermodynamical fits.',
        default=1000.0, required=False)
    parser.add_argument(
        '-p', '--plot', action='store_true',
        help='Generate plots when available.', required=False)
    args = parser.parse_args()

    mech = args.input
    mix_name = args.name
    ct_mixture = args.mixture
    nasa7_Tmid = args.Tmid

    output_dir = init_dirs(args.output)
    warning_msg.init_mixture_log(output_dir, mech, mix_name, ct_mixture)

    data = ct_properties.ctThermoTransport(mech, Tmid=nasa7_Tmid)
    data.evaluate_mixture_properties(mix_name, ct_mixture)

    print("\nFitting transport properties:")
    transport_fits = utils.fit_ct_transport(data, poly_order=3)
    success = utils.transport_fit_quality(data, transport_fits, output_dir, plot=args.plot)
    if(not success):
        print("\nMixture transport fit has failed. See Figures/ for visualisation.\n")
    else:
        print("Success.")

    print("\nFitting thermodynamic properties:")
    thermo_fits = utils.fit_mixture_thermo(data)
    success = utils.nasa7_fit_quality(data, thermo_fits, output_dir, plot=args.plot)
    if(not success):
        print("NASA7 polynomial fit has failed. See Figures/ for visualisation.\n")
    else:
        print("Success.")

    print("\nWriting mixture into OpenFOAM format")
    species_file = Path(output_dir, "species.foam")
    thermo_file = Path(output_dir, "thermo.foam")
    reactions_file = Path(output_dir, "reactions.foam")
    utils.ctmix2foam_thermo_writer(
        species_file, thermo_file, reactions_file, data, transport_fits,
        thermo_fits)

    print("\nDone")


if(__name__ == '__main__'):
    main()
