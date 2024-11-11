from pathlib import Path
from datetime import date

# private variables
_log_file = "log.txt"


def init_log(dir: Path, mechanism: str):
    from datetime import date
    today = date.today()    
    log_file = Path(dir, _log_file)
    with open(log_file, 'w') as output:
        output.write("Data is generated: " + str(today.strftime("%B %d, %Y")))
        output.write("\nMechanism used: " + mechanism)
        output.write("\n")


def init_mixture_log(dir: Path, mechanism: str, mixture_name: str, mixture: str):
    today = date.today()
    log_file = Path(dir, _log_file)
    with open(log_file, 'w') as output:
        output.write("Data is generated: " + str(today.strftime("%B %d, %Y")))
        output.write("\nMechanism used: " + mechanism)
        output.write("\nMixture - " + str(mixture_name) + ": " + mixture)
        output.write("\n")


def warning_sutherland(output_dir: Path, specie_name: str, err_mu: float, std: float, err_kappa: float):
    warning = (
        "\n\tWarning : " + str(specie_name) +
        " fit error is larger than tolerance.\n" + "\tSutherland mu error: " +
        str(err_mu) + "\n" + "\tSutherland STD: " + str(std) + "\n" +
        "\tEuken kappa error: " + str(err_kappa))
    print(warning)
    with open(Path(output_dir, _log_file), 'a') as f:
        print(warning, file=f)


def warning_polynomial(output_dir: Path, specie_name: str, err_mu: float, err_kappa: float, log_type: bool):
    if(log_type):
        polynomial_type = "Polynomial"
    else:
        polynomial_type = "Log-polynomial"

    warning = (
        "\n\tWarning : " + str(specie_name) +
        " fit error is larger than tolerance.\n" + "\t" + polynomial_type +
        " mu error: " + str(err_mu) + "\n" + "\t" + polynomial_type +
        " kappa error: " + str(err_kappa))
    print(warning)
    with open(Path(output_dir, _log_file), 'a') as f:
        print(warning, file=f)


def warning_not_nasa7(thermo, specie_name, output_dir):
    type_name = thermo.get_thermo_fit_type(specie_name)
    warning = ("\n\tWarning, " + specie_name + " : " + type_name +
               " thermo type will be fitted to NASA 7-coefficient format." +
               "\tUser is recommended to visually check the fit quality.")
    print(warning)
    with open(Path(output_dir, _log_file), 'a') as f:
        print(warning, file=f)


def error_consistency(sp_i: str, output_dir: Path):
    error = (
        "\n\tError, " + sp_i +
        " : there is a mismatch between Cantera thermo data and NASA polynomials."
        + "\tExiting...")
    print(error)
    with open(Path(output_dir, _log_file), 'a') as f:
        print(error, file=f)


def warning_nasa7_err(name: str, rel_tol: float, err_cp: float, err_h: float, err_s: float, output_dir: Path):
    warning = (
        "\n\ttWarning : " + str(name) +
        " NASA7 fit error is larger than tolerance " + str(rel_tol) +
        "\tcp error: " + str(err_cp) + "\n" + "\th error: " + str(err_h) + "\n" +
        "\ts error: " + str(err_s))
    print(warning)
    with open(Path(output_dir, _log_file), 'a') as f:
        print(warning, file=f)


def warning_nasa7_continuity(name: str, output_dir: Path):
    warning = ("\n\tWarning : " + str(name) + ": the new NASA polynomial fit is not continuous.")
    print(warning)
    with open(Path(output_dir, _log_file), 'a') as f:
        print(warning, file=f)        
