"""
ct2foam CLI entry point.

Converts Cantera thermodynamic and transport data into OpenFOAM-compatible
thermophysicalProperties dictionaries. Re-fits NASA7 polynomials and
Sutherland/polynomial transport functions where needed, and evaluates
fit quality against Cantera reference data. A single common mid-point
temperature (Tmid) is enforced across all species, as OpenFOAM requires
all species to share the same NASA7 transition temperature.
Run via ``ct2foam --help`` or ``python -m ct2foam``.
"""

import argparse
from datetime import date
from pathlib import Path


from ct2foam.species import SpeciesList
from ct2foam.mixture import Mixture


def main():
    """Parse CLI arguments and run the ct2foam conversion pipeline."""
    parser = argparse.ArgumentParser(
            "Convert cantera-based transport and thermodynamic"
            " data into OpenFOAM format."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        help="Cantera mechanism (.yaml/.xml) file path.",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="(Optional) Output directory path. (default is current directory)",
        default=Path.cwd(),
        required=False,
    )
    parser.add_argument(
        "-T",
        "--Tmid",
        type=float,
        help="(Optional) Common temperature for NASA-7 polynomials.",
        default=1000.0,
        required=False,
    )
    parser.add_argument(
        "-Tl",
        "--Tlow",
        type=float,
        help="(Optional) Temperature low-limit for NASA-7 polynomials.",
        default=280.0,
        required=False,
    )
    parser.add_argument(
        "-Th",
        "--Thigh",
        type=float,
        help="(Optional) Temperature high-limit for NASA-7 polynomials.",
        default=3000.0,
        required=False,
    )
    parser.add_argument(
        "-n",
        "--mixture-name",
        type=str,
        help='(Optional) Mixture name, e.g. "air".',
        required=False,
    )
    parser.add_argument(
        "-m",
        "--mixture",
        type=str,
        help='(Optional) Molecular mixture ratio in cantera style: "O2:1, N2:3.76" ',
        required=False,
    )
    parser.add_argument(
        "-p",
        "--plot",
        action="store_true",
        help="(Optional) Generate plots for fitted variables.",
        required=False,
    )
    parser.add_argument(
        "-t",
        "--tol-nasa7",
        type=float,
        help="(Optional) Tolerance for NASA-7 polynomial fit error.",
        default=1e-2,
        required=False,
    )
    parser.add_argument(
        "-tc0",
        "--tol-nasa7-c0",
        type=float,
        help="(Optional) Tolerance for NASA-7 polynomial C0 continuity.",
        default=1e-6,
        required=False,
    )
    parser.add_argument(
        "-tt",
        "--tol-transport",
        type=float,
        help="(Optional) Tolerance for all transport function fits.",
        default=1e-1,
        required=False,
    )

    args = parser.parse_args()

    mechanism = args.input

    output_dir = Path(args.output)
    output_dir.mkdir(exist_ok=True)
    fig_dir = Path(output_dir, "Figures")
    fig_dir.mkdir(exist_ok=True)

    # For mixtures, we need both arguments to be defined
    if args.mixture and not args.mixture_name:
        parser.error("--mixture_name must be specified")

    today = date.today()
    print("Date: " + str(today.strftime("%B %d, %Y")))
    print("Using Mechanism: " + mechanism)

    # Generate output for a mixture
    if args.mixture_name:
        print("Creating mixture - " + str(args.mixture_name) + ": " + args.mixture)
        mixture = Mixture.from_ct(
            mechanism_file=mechanism,
            mixture_name=args.mixture_name,
            mixture=args.mixture,
            Tmin=args.Tlow,
            Tmax=args.Thigh,
            Tmid=args.Tmid,
            n=256,
            plot=True,
            fig_dir=fig_dir,
            tol_nasa7=args.tol_nasa7,
            tol_nasa7_c0=args.tol_nasa7_c0,
            tol_transport=args.tol_transport,
        )

        mixture.write_foam(output_dir)
        print("\nDone")
        return

    # Generate output for individual species
    species_list = SpeciesList.from_ct_mech(
        mechanism,
        Tmin=args.Tlow,
        Tmax=args.Thigh,
        Tmid=args.Tmid,
        n=256,
        plot=args.plot,
        fig_dir=fig_dir,
        tol_nasa7=args.tol_nasa7,
        tol_nasa7_c0=args.tol_nasa7_c0,
        tol_transport=args.tol_transport,
    )
    species_list.write_foam(output_dir)
    print("\nDone")


if __name__ == "__main__":
    main()
