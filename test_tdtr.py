"""
Verification and accuracy test suite for TDTR Python implementation.
"""

import time
import numpy as np
import tdtr


def test_quadrature():
    print("Testing Gauss-Legendre quadrature...")
    nodes, weights = tdtr.lgwt(35, 0.0, 1.0)
    assert len(nodes) == 35
    assert len(weights) == 35
    integral_val = np.sum(weights * (nodes**3))
    assert np.isclose(integral_val, 0.25, atol=1e-10)
    print("  Gauss-Legendre quadrature test passed.")


def test_thermal_model():
    print("Testing thermal transfer matrix model and reflectance signals...")
    lambda_val = np.array([150.0, 1.0, 140.0])
    C_val = np.array([2.42e6, 0.1e6, 1.6e6])
    t_val = np.array([80e-9, 1e-9, 500e-6])
    eta_val = np.array([1.0, 1.0, 1.0])

    r_spot = 10e-6
    f = 9.3e6
    tau_rep = 1.0 / 80e6
    tdelay = np.logspace(-10, -9, 20)

    deltaR, ratio = tdtr.tdtr_refl(
        tdelay=tdelay,
        TCR=1e-4,
        tau_rep=tau_rep,
        f=f,
        lambda_vec=lambda_val,
        C_vec=C_val,
        t_vec=t_val,
        eta_vec=eta_val,
        r_pump=r_spot,
        r_probe=r_spot,
        A_pump=12e-3,
        nnodes=35,
    )

    assert len(ratio) == len(tdelay)
    assert np.all(np.isfinite(ratio))
    assert np.all(ratio > 0)
    print("  Thermal model test passed.")


def test_steady_state_heating():
    print("Testing steady-state temperature rise calculation...")
    lambda_val = np.array([150.0, 140.0])
    C_val = np.array([2.42e6, 1.6e6])
    t_val = np.array([80e-9, 500e-6])
    eta_val = np.array([1.0, 1.0])

    r_spot = 10e-6
    dT_ss = tdtr.ss_heating(
        lambda_vec=lambda_val,
        C_vec=C_val,
        t_vec=t_val,
        eta_vec=eta_val,
        r=r_spot,
        absorbance=0.1,
        A_tot_powermeter=12e-3,
    )
    assert np.isfinite(dT_ss)
    assert dT_ss > 0
    print(f"  Steady-state heating dT_SS = {dT_ss:.4f} K passed.")


def test_data_processing():
    print("Testing signal processing (phase correction and time-zero shift)...")
    time_raw = np.linspace(-30e-12, 100e-12, 100)
    Vin_raw = np.exp(-((time_raw - 5e-12) ** 2) / (2 * (10e-12) ** 2))
    Vin_raw[time_raw < 0] = 0.01
    Vout_raw = 0.2 * Vin_raw

    theta_true = 0.1
    Vin_mixed = Vin_raw * np.cos(theta_true) + Vout_raw * np.sin(theta_true)
    Vout_mixed = -Vin_raw * np.sin(theta_true) + Vout_raw * np.cos(theta_true)

    time_exp_shifted, t_half = tdtr.correct_time_zero(time_raw, Vin_mixed)
    assert np.isfinite(t_half)

    Vin_corr, Vout_corr, theta_found = tdtr.auto_correct_phase(time_exp_shifted, Vin_mixed, Vout_mixed)
    assert np.isfinite(theta_found)
    assert len(Vin_corr) == len(time_raw)
    print("  Data processing test passed.")


def test_acoustic_peaks():
    print("Testing acoustic echo baseline subtraction and peak detection...")
    t_ps = np.linspace(0, 100, 200)
    baseline_true = 0.5 - 0.002 * t_ps
    echo_true = 0.05 * np.exp(-((t_ps - 45.0) ** 2) / (2 * (2.0**2)))
    signal = baseline_true + echo_true

    echo_res = tdtr.find_acoustic_peaks(t_ps, signal, fit_window_ps=(35.0, 55.0), degree=1)
    assert np.isclose(echo_res["peak_time_ps"], 45.0, atol=1.0)
    print("  Acoustic echo peak detection test passed.")


def test_fitting():
    print("Testing parameter optimizer fitting...")
    lambda_val = np.array([150.0, 1.0, 140.0])
    C_val = np.array([2.42e6, 0.1e6, 1.6e6])
    t_val = np.array([80e-9, 1e-9, 500e-6])
    eta_val = np.array([1.0, 1.0, 1.0])

    sample = tdtr.Sample.from_arrays(lambda_val, C_val, t_val, eta_val)
    exp_params = tdtr.ExperimentParams(
        r_pump=10e-6, r_probe=10e-6, f=9.3e6, tau_rep=1.0 / 80e6, A_pump=12e-3
    )

    tdelay = np.logspace(-10, -9, 15)
    _, ratio_true = tdtr.tdtr_refl(
        tdelay,
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
    )

    fit_res = tdtr.fit_tdtr(
        ratio_data=ratio_true,
        tdelay_data=tdelay,
        sample=sample,
        exp_params=exp_params,
        fit_params=["k1"],
        Xguess=[120.0],
    )

    fitted_k1 = fit_res.params["k1"]
    assert np.isclose(fitted_k1, 150.0, rtol=1e-2)
    print(f"  Optimizer fit target 150.0 W/m-K -> recovered {fitted_k1:.2f} W/m-K. Test passed.")


if __name__ == "__main__":
    t0 = time.time()
    print("==================================================")
    print(" Running TDTR Test Suite")
    print("==================================================")
    test_quadrature()
    test_thermal_model()
    test_steady_state_heating()
    test_data_processing()
    test_acoustic_peaks()
    test_fitting()
    print(f"\nAll tests passed in {time.time() - t0:.3f} seconds!")
    print("==================================================")
