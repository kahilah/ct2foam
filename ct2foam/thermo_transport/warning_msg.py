import os

#private variables
_log_file = "log.txt"


def init_log(dir, mechanism):
    from datetime import date
    today = date.today()    
    log_file = os.path.join(dir, _log_file)
    with open(log_file,'w') as output:
        output.write("Data is generated: " + repr(today.strftime("%B %d, %Y")))
        output.write("\nMechanism used: " + mechanism)
        output.write("\n")


def init_mixture_log(dir, mechanism, mixture_name, mixture):
    from datetime import date
    today = date.today()    
    log_file = os.path.join(dir, _log_file)
    with open(log_file,'w') as output:
        output.write("Data is generated: " + repr(today.strftime("%B %d, %Y")))
        output.write("\nMechanism used: " + mechanism)
        output.write("\nMixture - " + repr(mixture_name) + ": " + mixture)
        output.write("\n")


def warning_sutherland(output_dir, specie_name, err_mu, std, err_kappa):
    warning = (
            "\n\tWarning : " + repr(specie_name) + " fit error is larger than tolerance.\n"
            "\tSutherland mu error: " + repr(err_mu) + "\n"
            "\tSutherland STD: " + repr(std) + "\n"
            "\tEuken kappa error: " + repr(err_kappa)
    )
    print(warning)
    with open(os.path.join(output_dir, _log_file),'a') as f:
        print(warning, file=f)


def warning_polynomial(output_dir, specie_name, err_mu, err_kappa, log_type):
    if(log_type):
        polynomial_type = "Polynomial"
    else:
        polynomial_type = "Log-polynomial"

    warning = (
            "\n\tWarning : " + repr(specie_name) + " fit error is larger than tolerance.\n"
            "\t" + polynomial_type + " mu error: " + repr(err_mu) + "\n"
            "\t" + polynomial_type + " kappa error: " + repr(err_kappa)
    )
    print(warning)
    with open(os.path.join(output_dir, _log_file),'a') as f:
        print(warning, file=f)


def warning_not_nasa7(thermo, specie_name, output_dir):
        type_name = thermo.get_thermo_fit_type(specie_name)
        warning = (
            "\n\tWarning, " + specie_name + " : "+ type_name + " thermo type will be fitted to NASA 7-coefficient format."
            "\tUser is recommended to visually check the fit quality."
            )
        print(warning)
        with open(os.path.join(output_dir, _log_file),'a') as f:
            print(warning, file=f)


def error_consistency(sp_i, output_dir):
    error = (
        "\n\tError, " + sp_i + " : there is a mismatch between Cantera thermo data and NASA polynomials."
        "\tExiting..."
        )
    print(error)
    with open(os.path.join(output_dir, _log_file),'a') as f:
        print(error, file=f)


def warning_nasa7_err(name, rel_tol, err_cp, err_h, err_s, output_dir):
    warning = (
            "\n\ttWarning : " + repr(name) + " NASA7 fit error is larger than tolerance " + repr(rel_tol) +
            "\tcp error: " + repr(err_cp) + "\n"
            "\th error: " + repr(err_h) + "\n"
            "\ts error: " + repr(err_s)
    )
    print(warning)
    with open(os.path.join(output_dir, _log_file),'a') as f:
        print(warning, file=f)


def warning_nasa7_continuity(name, output_dir):
    warning = ("\n\tWarning : " + repr(name) + ": the new NASA polynomial fit is not continuous.")
    print(warning)
    with open(os.path.join(output_dir, _log_file),'a') as f:
        print(warning, file=f)        