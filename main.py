"""
Main TDTR pipeline script.

Features:
1. Specify an experimental data file path (--data-file / -f) or analyze synthetic data.
2. Lock-in amplifier phase shift error correction (--phase-shift-deg / --autocorrect-phase).
3. Optional calculation flags:
   - Sensitivity analysis (--calc-sens)
   - Uncertainty / Errorbar analysis (--calc-errors)
   - 2D Parameter correlation scan (--calc-corr)
"""

import argparse
import os
import time
import numpy as np
import tdtr


def run_pipeline(
    data_file: str = r"C:\Users\d-cahill\OneDrive - University of Illinois - Urbana\Documents\Data\psec\2026\jul2506\9880_cen",
    time_min_ps: float = 100.0,
    time_max_ps: float = 3600.0,
    fit_params_list: list = ('k3','k4'),  #layers are number 1,2,3,4...
    Xguess_list: list = (0.1,2000),
    calc_sens: bool = False,
    calc_errors: bool = False,
    calc_corr: bool = False,
    phase_shift_deg: float = 0.0,
    autocorrect_phase: bool = True,
    save_plots: bool = True,
    show_plots: bool = True,
):
    start_time = time.time()
    print("==================================================")
    print(" TDTR Analysis Pipeline")
    print("==================================================")

    # 1. Sample Setup
    # typical Example stack: Al abosprtion layer (1 nm) / Al (60 nm) / Interface (1 nm) / Substrate (1 mm)
    lambda_val = [1500.0, 150.0, 0.1, 1000.0]  # W/m-K
    C_val = [24.2e6, 2.42e6, 0.1e6, 1.8e6]  # J/m^3-K
    t_val = [1.0e-9, 60e-9, 1.0e-9, 1e-3]  # m
    eta_val = [1.0] * len(lambda_val)

    sample = tdtr.Sample.from_arrays(
        lambda_array=lambda_val,
        C_array=C_val,
        t_array=t_val,
        eta_array=eta_val,
        names=["Al_top", "Al", "Interface", "Substrate"],
    )

    # 2. Experimental Setup
    r_spot = 9.7e-6  # spot radius m
    exp_params = tdtr.ExperimentParams(
        r_pump=r_spot,
        r_probe=r_spot,
        f=9.3e6,  # 9.3 MHz
        tau_rep=1.0 / 74.88e6,  # 74.88 MHz rep rate
        A_pump=12e-3,  # 12 mW
        TCR=1e-4,
        nnodes=35,
    )

    # 3. Steady-State Heating
    dT_ss = tdtr.ss_heating(
        lambda_vec=sample.lambda_vec,
        C_vec=sample.C_vec,
        t_vec=sample.t_vec,
        eta_vec=sample.eta_vec,
        r=r_spot,
        absorbance=0.12 * 0.9,
        A_tot_powermeter=(8.0 + 2 * 4.0) * 1e-3,
    )
    print(f"\nSteady-state temperature rise: dT_SS = {dT_ss:.2f} K")

    # 4. Load Data (Real or Synthetic)
    if data_file and os.path.exists(data_file):
        print(f"\nLoading experimental data from: {data_file}")
        exp_data = tdtr.read_exp_data(
            filepath=data_file,
            time_unit='ps',
            autocorrect_phase=autocorrect_phase,
            phase_shift_deg=phase_shift_deg,
        )
        print(f"Applied lock-in phase correction angle: {np.degrees(exp_data.phase_angle):.3f} deg")

        # Extract fitting range
        tdelay_fit, ratio_fit_data = tdtr.extract_interior(
            exp_data.time_exp, exp_data.ratio, time_min_ps * 1e-12, time_max_ps * 1e-12
        )
        _, Vin_fit = tdtr.extract_interior(exp_data.time_exp, exp_data.Vin, time_min_ps * 1e-12, time_max_ps * 1e-12)
        _, Vout_fit = tdtr.extract_interior(exp_data.time_exp, exp_data.Vout, time_min_ps * 1e-12, time_max_ps * 1e-12)
        print(f"Extracted {len(tdelay_fit)} data points in [{time_min_ps}, {time_max_ps}] ps range.")
    else:
        if data_file:
            print(f"\nWarning: Specified data file '{data_file}' not found. Falling back to synthetic demonstration data.")
        else:
            print("\nNo experimental data file specified. Generating synthetic demonstration data.")

        tdelay_fit = np.logspace(np.log10(time_min_ps * 1e-12), np.log10(time_max_ps * 1e-12), 60)
        _, ratio_model_synthetic = tdtr.tdtr_refl(
            tdelay_fit,
            exp_params.TCR,
            exp_params.tau_rep,
            exp_params.f,
            sample.lambda_vec,
            sample.C_vec,
            sample.t_vec,
            sample.eta_vec,
            exp_params.r_pump,
            exp_params.r_probe,
            exp_params.A_pump,
            exp_params.nnodes,
        )
        np.random.seed(42)
        noise = np.random.normal(0, 0.005, size=len(ratio_model_synthetic))
        ratio_fit_data = ratio_model_synthetic * (1.0 + noise)
        Vin_fit = None
        Vout_fit = None

    # 5. Model Fitting
    if fit_params_list is None:
        fit_params_list = ["k3", "k4"]
    if Xguess_list is None:
        Xguess_list = [0.1, 140.0]

    print(f"\nRunning parameter fit for {fit_params_list}...")
    t_fit_start = time.time()
    fit_result = tdtr.fit_tdtr(
        ratio_data=ratio_fit_data,
        tdelay_data=tdelay_fit,
        sample=sample,
        exp_params=exp_params,
        fit_params=fit_params_list,
        Xguess=Xguess_list,
    )
    print(f"Fit completed in {time.time() - t_fit_start:.3f} s.")
    print("Fit Results:")
    for p_name, val in fit_result.params.items():
        print(f"  {p_name} = {val:.4f}")
    print(f"  Residual sum Z = {fit_result.Z_min:.6e}")

    if save_plots or show_plots:
        tdtr.plot_fit_result(
            tdelay_data=tdelay_fit,
            ratio_data=ratio_fit_data,
            ratio_model=fit_result.ratio_model,
            sample=sample,
            title="TDTR Parameter Fit Result",
            show=show_plots,
            save_path="tdtr_fit_result.png" if save_plots else None,
        )
        if save_plots:
            print("Saved 'tdtr_fit_result.png'")

    # 6. Optional: Sensitivity Calculation
    if calc_sens:
        print("\nCalculating sensitivities (OPTIONAL)...")
        t_sens_start = time.time()
        _, sens_dict = tdtr.calculate_sensitivities(
            tdelay=tdelay_fit,
            TCR=exp_params.TCR,
            tau_rep=exp_params.tau_rep,
            f=exp_params.f,
            lambda_vec=sample.lambda_vec,
            C_vec=sample.C_vec,
            t_vec=sample.t_vec,
            eta_vec=sample.eta_vec,
            r_pump=exp_params.r_pump,
            r_probe=exp_params.r_probe,
            A_pump=exp_params.A_pump,
            nnodes=exp_params.nnodes,
        )
        print(f"Sensitivities calculated in {time.time() - t_sens_start:.3f} s.")
        if save_plots or show_plots:
            tdtr.plot_sensitivities(
                tdelay=tdelay_fit,
                sensitivities=sens_dict,
                params_to_plot=["k2", "k3", "k4", "t2", "r_pump"],
                title="TDTR Logarithmic Sensitivities",
                show=show_plots,
                save_path="tdtr_sensitivities.png" if save_plots else None,
            )
            if save_plots:
                print("Saved 'tdtr_sensitivities.png'")

    # 7. Optional: Uncertainty Estimation
    if calc_errors:
        print("\nCalculating error bars / uncertainties (OPTIONAL)...")
        uncertainties = tdtr.compute_uncertainties(
            ratio_data=ratio_fit_data,
            tdelay_data=tdelay_fit,
            sample=sample,
            exp_params=exp_params,
            fit_params=fit_params_list,
            Vin_data=Vin_fit,
            Vout_data=Vout_fit,
        )
        for p_name, p_val, err_pct, err_abs in zip(
            uncertainties.param_names,
            uncertainties.Xsol,
            uncertainties.total_percent_err,
            uncertainties.total_abs_err,
        ):
            print(f"  {p_name} = {p_val:.4f} +/- {err_abs:.4f} ({err_pct * 100:.2f}%)")

    # 8. Optional: 2D Parameter Correlation Scan
    if calc_corr and len(fit_params_list) == 2:
        print("\nRunning 2D Parameter Correlation Scan (OPTIONAL)...")
        scan_results = tdtr.parameter_correlation_scan(
            ratio_data=ratio_fit_data,
            tdelay_data=tdelay_fit,
            sample=sample,
            exp_params=exp_params,
            fit_params=fit_params_list,
            scan_param=fit_params_list[0],
            refit_param=fit_params_list[1],
            scan_fraction=0.30,
            Nscan=15,
        )
        if save_plots or show_plots:
            tdtr.plot_correlation_scan(
                scan_results=scan_results,
                show=show_plots,
                save_path="tdtr_correlation_scan.png" if save_plots else None,
            )
            if save_plots:
                print("Saved 'tdtr_correlation_scan.png'")

    total_time = time.time() - start_time
    # print("==================================================")
    # print(f" Pipeline finished successfully in {total_time:.3f} s!")
    # print("==================================================")


def main():
    parser = argparse.ArgumentParser(description="TDTR Thermal Modeling and Data Analysis")
    parser.add_argument("-f", "--data-file", type=str, default=None, help="Path to experimental TDTR data file")
    parser.add_argument("--tmin", type=float, default=None, help="Min delay time for fitting in ps")
    parser.add_argument("--tmax", type=float, default=None, help="Max delay time for fitting in ps")
    parser.add_argument("--fit-params", nargs="+", default=None, help="Parameters to fit (e.g., k2 k4)")
    parser.add_argument("--Xguess", nargs="+", type=float, default=None, help="Initial guesses for fit parameters")
    parser.add_argument("--calc-sens", action="store_true", help="Enable ratio sensitivity calculation & plot")
    parser.add_argument("--calc-errors", action="store_true", help="Enable errorbar / uncertainty estimation")
    parser.add_argument("--calc-corr", action="store_true", help="Enable 2D parameter correlation scan & plot")
    parser.add_argument("--phase-shift-deg", type=float, default=None, help="Manual lock-in phase shift error angle (degrees)")
    parser.add_argument("--no-autocorrect-phase", action="store_false", dest="autocorrect_phase", help="Disable automatic phase jump correction")
    parser.add_argument("--no-show-plots", action="store_false", dest="show_plots", help="Disable displaying interactive plot windows")

    args = parser.parse_args()

    kwargs = {}
    if args.data_file is not None:
        kwargs["data_file"] = args.data_file
    if args.tmin is not None:
        kwargs["time_min_ps"] = args.tmin
    if args.tmax is not None:
        kwargs["time_max_ps"] = args.tmax
    if args.fit_params is not None:
        kwargs["fit_params_list"] = args.fit_params
    if args.Xguess is not None:
        kwargs["Xguess_list"] = args.Xguess
    if args.phase_shift_deg is not None:
        kwargs["phase_shift_deg"] = args.phase_shift_deg
    if not args.show_plots:
        kwargs["show_plots"] = False

    run_pipeline(
        calc_sens=args.calc_sens,
        calc_errors=args.calc_errors,
        calc_corr=args.calc_corr,
        autocorrect_phase=args.autocorrect_phase,
        **kwargs,
    )


if __name__ == "__main__":
    main()
